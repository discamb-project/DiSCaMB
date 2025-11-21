#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/IO/atom_type_io.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/MathUtilities/real_spherical_harmonics.h"
#include "discamb/HC_Model/hc_model_utilities.h"
#include "discamb/BasicUtilities/discamb_version.h"
#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/MathUtilities/statistics.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/AtomTyping/atom_typing_utilities.h"
#include "discamb/IO/cif_io.h"
#include "discamb/IO/shelx_io.h"
#include "discamb/IO/mol2_io.h"
#include "discamb/IO/xd_io.h"
#include "discamb/IO/xyz_io.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <deque>
#include <chrono>
#include <ctime>

#include <omp.h>

using namespace discamb;
using namespace std;


void getXdFilePaths(
    vector<pair<string, string> > &xdFilePath,
    vector<string> &ids,
    const string &path_part_1,
    const string &path_part_2)
{
    filesystem::path refCodeDirFolder = (path_part_1.empty() ? filesystem::current_path() : filesystem::path(path_part_1));

    ids.clear();

    for (auto &p : filesystem::directory_iterator(refCodeDirFolder))
    {

        if (filesystem::is_directory(p.status()))
        {
            filesystem::path xdMasPath = p.path();
            filesystem::path xdResPath = p.path();
            xdMasPath /= path_part_2;
            xdMasPath /= "xd.mas";
            xdResPath /= path_part_2;
            xdResPath /= "/xd.res";

            if (filesystem::exists(xdMasPath) && filesystem::exists(xdResPath))
                if (filesystem::is_regular_file(xdMasPath) && filesystem::is_regular_file(xdResPath))
                {

                    string pathStr = p.path().string();
                    ids.push_back(pathStr.substr(pathStr.rfind(filesystem::path::preferred_separator) + 1));
                    xdFilePath.push_back({ xdMasPath.string(), xdResPath.string() });
                }
        }
    }

}


string atomicNumbers2Formula(
    const vector<size_t> &atomicNumbers)
{
    map<string, size_t> formula;
    string symbol, formulaStr;
    for (auto z : atomicNumbers)
    {
        symbol = periodic_table::symbol(z);
        if (formula.find(symbol) != formula.end())
            formula[symbol]++;
        else
            formula[symbol] = 1;
    }
    for (auto x : formula)
    {
        formulaStr += x.first;
        if (x.second != 1)
            formulaStr += to_string(x.second);
    }
    return formulaStr;
}

struct Settings {
    DescriptorsSettings descriptorsSettings;
    double min_plm = 0.002;
    double nSigma = 1.0;
    size_t min_n_instaces = 3;
    string files_folder = string("");
    bool printAssignmentInfo = true;
    bool printMoreDetailedStats = true;
};



void getFilePaths(
    const string &cifs_folder,
    vector<string> &filePaths,
    vector<string> &ids,
    const string &suffix) // e.g. .cif or .res
{
    filesystem::path dirFolder = (cifs_folder.empty() ? filesystem::current_path() : filesystem::path(cifs_folder));

    ids.clear();
    filePaths.clear();

    for (auto &p : filesystem::directory_iterator(dirFolder))

        if (filesystem::is_regular_file(p.status()))
            if (p.path().extension() == suffix)
            {
                ids.push_back(p.path().stem().string());
                filePaths.push_back(p.path().string());
            }
}


struct Atoms1stNeighbours
{
    string symbol;
    string formula;
    bool operator > (const Atoms1stNeighbours &x) {
        pair<string, string> p1(symbol, formula), p2(x.symbol, x.formula);
        return p1 > p2;
    }

};

bool operator < (const Atoms1stNeighbours &x1, const Atoms1stNeighbours &x2) {
    pair<string, string> p1(x1.symbol, x1.formula), p2(x2.symbol, x2.formula);
    return p1 < p2;
}


bool less1(const Atoms1stNeighbours &v1, const Atoms1stNeighbours &v2)
{
    pair<string, string> p1(v1.symbol, v1.formula), p2(v2.symbol, v2.formula);
    return p1 < p2;
}

struct InstanceAtom1stNeighboursNeighbours
{
    size_t atomIdx, structureIdx;
    vector<Atoms1stNeighbours> neighbors;
};





bool less2(const InstanceAtom1stNeighboursNeighbours &v1, const InstanceAtom1stNeighboursNeighbours &v2)
{
    return v1.neighbors < v2.neighbors;
}



void printDetailedStats(
    const string &fName,
    const string &header,
    const vector<pair<string, vector<size_t> > > &multitypes,
    const vector< vector<vector<pair<size_t, size_t> > > > &typeInstances,
    const vector<string> &structureNames,
    const vector<StructureWithDescriptors> structuresDescriptors,
    size_t nAssignedAtoms,
    size_t nNotAssignedAtoms,
	size_t nConsideredStructures,
	size_t nCompletelyRecognizedStructures)
{
        ofstream out(fName);
        size_t nTypes, typeIdx, i, nInstances, nSubtypes, subtypeIdx, atomIdx, structureIdx, nNeighbours, neighbourIdx;
    
    
        if (!out.good())
            on_error::throwException(string("cannot print to file '") + fName + string("'"), __FILE__, __LINE__);
    
        string hash80(80, '#');
        out << hash80 << "\n" << header << hash80 << "\n\n"
            << "   content:\n"
            << "      (1) fraction and number of atoms with assigned atom type\n"
			<< "      (2) fraction and number of structures with all atoms with atom type assigned\n"
            << "      (3) list of types which instances were found\n"
            << "      (4) list of types which instances were not found\n"
            << "      (5) number of instances per type\n"
            << "      (6) list of instances for each type \n"
            << "\n\n" << hash80 << "\n\n\n";
    

        //------------------------------------
        
        vector<size_t> nInstancesPerMultitype;
        nTypes = multitypes.size();
        nInstancesPerMultitype.assign(nTypes, 0);
        for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {
            nSubtypes = typeInstances[typeIdx].size();
            for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
                nInstancesPerMultitype[typeIdx] += typeInstances[typeIdx][subtypeIdx].size();
        }
        
        //------------------------------------

        out << "\n(1) fraction and number of atoms with assigned atom type\n";
        double nAllAtoms = double(nAssignedAtoms + nNotAssignedAtoms);
        out << "    atom type assigned to " << setprecision(2) << fixed << nAssignedAtoms / nAllAtoms * 100.0 << "% of atoms (" << nAssignedAtoms << ")\n";
        cout << "\n\n    atom type assigned to " << setprecision(2) << fixed << nAssignedAtoms / nAllAtoms * 100.0 << "% of atoms (" << nAssignedAtoms << ")\n";

        cout << "\n\n    atom type not assigned to " << nNotAssignedAtoms << " atoms\n";

		out << "\n(2) fraction and number of structures with all atoms with atom type assigned\n";
		out << "    " << nCompletelyRecognizedStructures << " of " << nConsideredStructures << " with completely recognized atom types "
			<< " - " << setprecision(2) << fixed << double(nCompletelyRecognizedStructures) / double(nConsideredStructures) * 100.0 << "%\n";

        cout << "\n    fraction and number of structures with all atoms with atom type assigned\n";
        cout << "    " << nCompletelyRecognizedStructures << " of " << nConsideredStructures << " with completely recognized atom types "
            << " - " << setprecision(2) << fixed << double(nCompletelyRecognizedStructures) / double(nConsideredStructures) * 100.0 << "%\n";
    
        out << "(3) list of types which instances were found\n\n";
        i = 0;
        for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {

            if (nInstancesPerMultitype[typeIdx] > 0)
            {
                out << multitypes[typeIdx].first << " ";
                i++;
                if ((i + 1) % 10 == 0)
                    out << "\n";
            }
        }
        out << "\n";
    
    
    
        out << "(4) list of types which instances were not found\n\n";
        
        i = 0;
        for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {

            if (nInstancesPerMultitype[typeIdx] == 0)
            {
                out << multitypes[typeIdx].first << " ";
                i++;
                if ((i + 1) % 10 == 0)
                    out << "\n";
            }
        }
        out << "\n";
    
            
    


        out << "\n(5) number of instances per type\n\n";
        
    
        for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
            out << multitypes[typeIdx].first << " " << setw(8) << nInstancesPerMultitype[typeIdx] << "\n";
    
        out << "(6) list of instances for each type\n\n";
    
        vector<InstanceAtom1stNeighboursNeighbours> instanceAtomsNeighbours;
        vector<vector<size_t> > shells;
        vector<size_t> neighborNeighbours;

        // {first neighbour symbols, it neighbours formula}
        pair<string,vector<string> > firstNeighborsFormula;
        set<size_t> secondNeighbors;
        string formula;
        vector<size_t> atomicNumbers;
        
        
    for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        out << "\n" << "TYPE " << multitypes[typeIdx].first << "\n";
        nInstances = nInstancesPerMultitype[typeIdx];
        out << " number of instances " << nInstances << "\n";

        //formula_of_2nd_neighbours_and_idx.resize(nInstances);
        instanceAtomsNeighbours.clear();
        instanceAtomsNeighbours.resize(nInstances);
        nSubtypes = typeInstances[typeIdx].size();

        size_t instancesCounter = 0;
        for (subtypeIdx=0 ; subtypeIdx < nSubtypes ; subtypeIdx++)
        {
            nInstances = typeInstances[typeIdx][subtypeIdx].size();

            for (i = 0; i < nInstances; i++)
            {
                structureIdx = typeInstances[typeIdx][subtypeIdx][i].first;
                atomIdx = typeInstances[typeIdx][subtypeIdx][i].second;
                nNeighbours = structuresDescriptors[structureIdx].connectivity[atomIdx].size();

                instanceAtomsNeighbours[instancesCounter].atomIdx = atomIdx;
                instanceAtomsNeighbours[instancesCounter].structureIdx = structureIdx;
                instanceAtomsNeighbours[instancesCounter].neighbors.resize(nNeighbours);
                firstNeighborsFormula.second.clear();
                for (size_t j = 0; j < nNeighbours; j++)  // central atom 1-st neighbors
                {
                    neighbourIdx = structuresDescriptors[structureIdx].connectivity[atomIdx][j];
                    instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol = firstNeighborsFormula.first = periodic_table::symbol(structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].atomicNumber);
                    const vector<int> &neighbours = structuresDescriptors[structureIdx].connectivity[neighbourIdx];
                    atomicNumbers.clear();
                    for (auto &neighbor : neighbours) // 1-st neighbor neighbors
                        atomicNumbers.push_back(structuresDescriptors[structureIdx].atomDescriptors[neighbor].atomicNumber);

                    instanceAtomsNeighbours[instancesCounter].neighbors[j].formula = atomicNumbers2Formula(atomicNumbers);
                }
                sort(instanceAtomsNeighbours[instancesCounter].neighbors.begin(), instanceAtomsNeighbours[instancesCounter].neighbors.end(), less1);
                instancesCounter++;
            }

        }
            

        std::sort(instanceAtomsNeighbours.begin(), instanceAtomsNeighbours.end(),less2);

        nInstances = instanceAtomsNeighbours.size();
        //formula, instance
        map<string, vector<pair<size_t, size_t> > > formulaGrouppedInstances;
        for (i = 0; i < nInstances; i++)
        {

            formula.clear();
            // ----- 1-st neighbors neighbors formula as string
            for (size_t neighIdx = 0; neighIdx < instanceAtomsNeighbours[i].neighbors.size(); neighIdx++)
            {
                if (neighIdx != 0)
                    formula += ";";
                formula += instanceAtomsNeighbours[i].neighbors[neighIdx].symbol +
                    string("(") + instanceAtomsNeighbours[i].neighbors[neighIdx].formula + string(")");

            }

            // -------
            structureIdx = instanceAtomsNeighbours[i].structureIdx;
            atomIdx = instanceAtomsNeighbours[i].atomIdx;
            formulaGrouppedInstances[formula].push_back({ structureIdx, atomIdx });
        }
        vector<pair<size_t, string> > neighbourGroups;
        for (auto &formula : formulaGrouppedInstances)
            neighbourGroups.push_back({ formula.second.size(), formula.first });
        sort(neighbourGroups.begin(), neighbourGroups.end(), greater< pair<size_t, string> >());

        for (auto const &group : neighbourGroups)
        {
            const vector<pair<size_t, size_t> > &cases = formulaGrouppedInstances[group.second];
            out << "  " << group.first << " instances with neighbours formula: " << group.second << "\n";

            for (auto &c: cases)
            {
                string label = structureNames[c.first] + string(",") + structuresDescriptors[c.first].atomDescriptors[c.second].label; //stats[typeIdx].atomLabels[index].first + string(",") + stats[typeIdx].atomLabels[index].second;
                out << setw(label.size()+6) << label << "\n";
            }
        }
    }
    
        out.close();

}

void printMoreDetailedDetailedStats(
    const string &fName,
    const string &header,
    const vector<pair<string, vector<size_t> > > &multitypes,
    const vector< vector<vector<pair<size_t, size_t> > > > &typeInstances,
    const vector<string> &structureNames,
    const vector<StructureWithDescriptors> structuresDescriptors)
{
    ofstream out(fName);
    size_t nTypes, typeIdx, i, nInstances, nSubtypes, subtypeIdx, atomIdx, structureIdx, nNeighbours, neighbourIdx;


    if (!out.good())
        on_error::throwException(string("cannot print to file '") + fName + string("'"), __FILE__, __LINE__);

    string hash80(80, '#');
    out << hash80 << "\n" << header << hash80 << "\n\n"
        << "   content: list of instances for each type\n\n";


    //------------------------------------

    vector<size_t> nInstancesPerMultitype;
    nTypes = multitypes.size();
    nInstancesPerMultitype.assign(nTypes, 0);
    for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        nSubtypes = typeInstances[typeIdx].size();
        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
            nInstancesPerMultitype[typeIdx] += typeInstances[typeIdx][subtypeIdx].size();
    }

    //------------------------------------



    vector<InstanceAtom1stNeighboursNeighbours> instanceAtomsNeighbours;
    vector<vector<size_t> > shells;
    vector<size_t> neighborNeighbours;

    // {first neighbour symbols, it neighbours formula}
    pair<string, vector<string> > firstNeighborsFormula;
    set<size_t> secondNeighbors;
    string formula;
    vector<size_t> atomicNumbers;


    for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        out << "\n" << "TYPE " << multitypes[typeIdx].first << "\n";
        nInstances = nInstancesPerMultitype[typeIdx];
        out << " number of instances " << nInstances << "\n";

        //formula_of_2nd_neighbours_and_idx.resize(nInstances);
        instanceAtomsNeighbours.clear();
        instanceAtomsNeighbours.resize(nInstances);
        nSubtypes = typeInstances[typeIdx].size();
        size_t instancesCounter = 0;
        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
        {
            nInstances = typeInstances[typeIdx][subtypeIdx].size();

            for (i = 0; i < nInstances; i++)
            {
                structureIdx = typeInstances[typeIdx][subtypeIdx][i].first;
                atomIdx = typeInstances[typeIdx][subtypeIdx][i].second;
                nNeighbours = structuresDescriptors[structureIdx].connectivity[atomIdx].size();
                //instanceAtomsNeighbours[i].instanceIdx = i;
                instanceAtomsNeighbours[instancesCounter].atomIdx = atomIdx;
                instanceAtomsNeighbours[instancesCounter].structureIdx = structureIdx;
                instanceAtomsNeighbours[instancesCounter].neighbors.resize(nNeighbours);
                firstNeighborsFormula.second.clear();
                for (size_t j = 0; j < nNeighbours; j++)  // central atom 1-st neighbors
                {
                    neighbourIdx = structuresDescriptors[structureIdx].connectivity[atomIdx][j];
                    instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol = firstNeighborsFormula.first = periodic_table::symbol(structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].atomicNumber);
                    if (structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].planar == discamb::Tribool::True)
                        instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += ",planar";
                    if (structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].planar == discamb::Tribool::False)
                        instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += ",not planar";

                    vector<int> const & ringSizes = structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].planarRingsSizes;

                    if (!ringSizes.empty())
                    {
                        instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += ",rings(";
                        for (size_t k = 0; k < ringSizes.size(); k++)
                        {
                            if (k > 0)
                                instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += ",";
                            instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += to_string(ringSizes[k]);
                        }
                        instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += ")";
                    }

                    size_t n3Rings = structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].n3memberRings;
                    size_t n4Rings = structuresDescriptors[structureIdx].atomDescriptors[neighbourIdx].n4memberRings;

                    if (n3Rings!=0)
                        instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += string(",3 rings: ") + to_string(n3Rings);
                    if (n4Rings != 0)
                        instanceAtomsNeighbours[instancesCounter].neighbors[j].symbol += string(",4 rings: ") + to_string(n4Rings);

                    const vector<int> &neighbours = structuresDescriptors[structureIdx].connectivity[neighbourIdx];
                    atomicNumbers.clear();
                    for (auto &neighbor : neighbours) // 1-st neighbor neighbors
                        atomicNumbers.push_back(structuresDescriptors[structureIdx].atomDescriptors[neighbor].atomicNumber);



                    instanceAtomsNeighbours[instancesCounter].neighbors[j].formula = atomicNumbers2Formula(atomicNumbers);
                }
                sort(instanceAtomsNeighbours[instancesCounter].neighbors.begin(), instanceAtomsNeighbours[instancesCounter].neighbors.end(), less1);
                instancesCounter++;
            }

        }


        std::sort(instanceAtomsNeighbours.begin(), instanceAtomsNeighbours.end(), less2);

        nInstances = instanceAtomsNeighbours.size();
        //formula, instance
        map<string, vector<pair<size_t, size_t> > > formulaGrouppedInstances;
        for (i = 0; i < nInstances; i++)
        {

            formula.clear();
            // ----- 1-st neighbors neighbors formula as string
            for (size_t neighIdx = 0; neighIdx < instanceAtomsNeighbours[i].neighbors.size(); neighIdx++)
            {
                if (neighIdx != 0)
                    formula += ";";
                formula += instanceAtomsNeighbours[i].neighbors[neighIdx].symbol +
                    string(", (") + instanceAtomsNeighbours[i].neighbors[neighIdx].formula + string(")");

            }

            // -------
            structureIdx = instanceAtomsNeighbours[i].structureIdx;
            atomIdx = instanceAtomsNeighbours[i].atomIdx;
            formulaGrouppedInstances[formula].push_back({ structureIdx, atomIdx });
        }
        vector<pair<size_t, string> > neighbourGroups;
        for (auto &formula : formulaGrouppedInstances)
            neighbourGroups.push_back({ formula.second.size(), formula.first });
        sort(neighbourGroups.begin(), neighbourGroups.end(), greater< pair<size_t, string> >());

        for (auto const &group : neighbourGroups)
        {
            const vector<pair<size_t, size_t> > &cases = formulaGrouppedInstances[group.second];
            out << "  " << group.first << " instances with neighbours formula: " << group.second << "\n";

            for (auto &c : cases)
            {
                string label = structureNames[c.first] + string(",") + structuresDescriptors[c.first].atomDescriptors[c.second].label; //stats[typeIdx].atomLabels[index].first + string(",") + stats[typeIdx].atomLabels[index].second;
                out << setw(label.size() + 6) << label << "\n";
            }
        }
    }

    out.close();

}

string secondNeighboursLabel(
    const StructureWithDescriptors& structureDescriptors, 
    size_t &atomIdx)
{
    vector < vector<int> > neighbours;
    graph_algorithms::breadth_first_search(structureDescriptors.connectivity, atomIdx, neighbours, 2);
    if (neighbours.size() < 3)
        return string("0");
    size_t n = neighbours[2].size();
    set<string> symbols;
    for (size_t i = 0; i < n; i++)
        symbols.insert(periodic_table::symbol(structureDescriptors.atomDescriptors[neighbours[2][i]].atomicNumber));
    string result;
    for (auto& symbol : symbols)
        result += symbol;
    result += to_string(n);

    return result;

}


typedef tuple<size_t, string, string, string, size_t, size_t, string > UnassignedAtomDescriptors;

string formulaAsString(const multiset<int>& f)
{
    string formula;
    map<int, int> elementCount;
    for (int z : f)
    {
        if (elementCount.count(z) == 0)
            elementCount[z] = 1;
        else
            elementCount[z]++;
    }
    for (auto& item : elementCount)
    {
        formula += periodic_table::symbol(item.first);
        if(item.second>1)
            formula += to_string(item.second);
    }
    return formula;
}

void printStructureAssignedAndNot(
    const string& fName,
    const vector<string>& structureIds,
    const vector<vector<int> >& typeIndices,
    const vector<StructureWithDescriptors> &structureDescriptors)
{
    list<string> assigned, not_assigned, one_not_assigned;
    vector<string> one_not_assigned_atom_label;
    int i, n = structureIds.size();
    
    for (i = 0; i < n; i++)
    {
        int nUnrecognized = count(typeIndices[i].begin(), typeIndices[i].end(), -1);
        if (nUnrecognized == 0)
            assigned.push_back(structureIds[i]);
        else
            not_assigned.push_back(structureIds[i]);
        
        if (nUnrecognized == 1)
        {
            one_not_assigned.push_back(structureIds[i]);
            auto it = find(typeIndices[i].begin(), typeIndices[i].end(), -1);
            int atomIdx = distance(typeIndices[i].begin(), it);
            one_not_assigned_atom_label.push_back(structureDescriptors[i].atomDescriptors[atomIdx].label);
        }
    }

    ofstream out(fName);
    if (!out.good())
        on_error::throwException("cannot write to output file " + fName, __FILE__, __LINE__);

    out << " n assigned " << assigned.size() << endl;
    out << " n not assigned " << not_assigned.size() << endl;
    out << " n one not assigned " << assigned.size() << endl;
    
    out << "\n ASSIGNED \n\n";
    for (string& id : assigned)
        out << id << "\n";

    out << "\n\n NOT ASSIGNED \n\n";

    for (string& id : not_assigned)
        out << id << "\n";

    out << "\n\n ONE NOT ASSIGNED \n\n";
    i = 0;
    for (string& id : one_not_assigned)
        out << id << "   " << one_not_assigned_atom_label[i++] << "\n";


    out.close();
}

void printMoveAssigned(
    const string& fName,
    const vector<vector<string> >& filePaths,
    const vector<vector<int> >& typeIndices)
{
    ofstream out(fName);
    if (!out.good())
        on_error::throwException(string("cannot write file '") + fName + string("'"), __FILE__, __LINE__);

    int structureIdx, nStructures = typeIndices.size();

    for (structureIdx = 0; structureIdx < nStructures; structureIdx++)
        if (find(typeIndices[structureIdx].begin(), typeIndices[structureIdx].end(), -1) == typeIndices[structureIdx].end())
            for (const string& path : filePaths[structureIdx])
                out << "MOVE " << path << " recognized\n";

    out.close();
}


void printAssignmentInfo(
    const string &fName,
    const vector<string> &structureIds,
    const std::string &header,
    const DescriptorsSettings &descriptorsSettings,
    const vector<StructureWithDescriptors> &structureDescriptors,
    const vector<vector<int> > &typeIndices,
    const vector<AtomType> &types,
    //const vector<vector<UBDB_LocalCoordinateSystem<size_t> > > &lcs,
    const vector < vector <string> > &lcs,
    map< UnassignedAtomDescriptors, vector<string> > &sortedUnassignedAtoms,
    const vector<size_t> &nAtomsInAsymmetricUnit)
{
    sortedUnassignedAtoms.clear();
    UnassignedAtomDescriptors unassignedAtomDescriptors;
    ofstream out(fName);
    char xyz[] = { 'X', 'Y', 'Z' };
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fName + string("'"), __FILE__, __LINE__);

    size_t nStructures, nAtoms, atomIdx, structureIdx;
    string ringsSizesStr;
    vector<size_t> ringsSizes;
    string hash80(80, '#');
    nStructures = structureDescriptors.size();

    out << hash80 << "\n" << header << hash80 << "\n\n"
        << "       Structural Descriptors for each structure and atom and\n"
        << "       atom type assignemnt and local cordinate system information\n\n" << hash80 << "\n\n"
        << "       Number of structures " << nStructures << "\n\n" << hash80 << "\n\n";

    vector<string> labels;

    for (structureIdx = 0; structureIdx < nStructures; structureIdx++)
    {
        labels.clear();
        for (auto const &atomData : structureDescriptors[structureIdx].atomDescriptors)
            labels.push_back(atomData.label);

        out << "\n Structure " << structureIdx + 1 << " " << structureIds[structureIdx] << endl;
        out << "\n number of planar rings: " << structureDescriptors[structureIdx].planarRings.size() << endl;
        for (size_t ringIdx = 0; ringIdx < structureDescriptors[structureIdx].planarRings.size(); ringIdx++)
        {
            out << "\n  ring " << ringIdx + 1 << "\n    planarity esd " << structureDescriptors[structureIdx].planarRingsPlanarityEsd[ringIdx];
            out << "\n    atoms: ";
            for (size_t i = 0; i < structureDescriptors[structureIdx].planarRings[ringIdx].size(); i++)
                out << " " << labels[structureDescriptors[structureIdx].planarRings[ringIdx][i]];
            out << "\n";
        }

        //------------------------------------------------------------------------------------


        nAtoms = nAtomsInAsymmetricUnit[structureIdx];
        out << "\n Atomic Descriptors\n";
        out << "\n number of atoms " << nAtoms << endl;

        for (size_t atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {

            string symbol;
            string formulaStr = formulaAsString(structureDescriptors[structureIdx].atomDescriptors[atomIdx].neighborsFormula);
            string planarityStr;

            auto planarity = structureDescriptors[structureIdx].atomDescriptors[atomIdx].planar;
            if (planarity == discamb::Tribool::True)
                planarityStr = "planar";
            if (planarity == discamb::Tribool::False)
                planarityStr = "not planar";
            if (planarity == discamb::Tribool::True)
                planarityStr = "planar";
            if (planarity == discamb::Tribool::Undefined)
                planarityStr = "maybe planar";


            out << setw(8) << labels[atomIdx]
                << " atomic number " << setw(4) << structureDescriptors[structureIdx].atomDescriptors[atomIdx].atomicNumber
                << " , neighbours formula " << setw(12) << formulaStr
                << " , is planar " << setw(15) << planarityStr
                << " , planarity esd " << setw(10) << setprecision(6) << fixed << structureDescriptors[structureIdx].atomDescriptors[atomIdx].planarityDisplacementEsd
                << " , belongs to planar rings: ";
            
            ringsSizesStr.clear();

            if (structureDescriptors[structureIdx].atomDescriptors[atomIdx].planarRingsIndices.empty())
            {
                out << "  no ring size: ";
                ringsSizesStr = " - ";
            }
            else
            {
                out << " yes ring size: "; 
                ringsSizes.clear();
                const vector<int> &ringsIndices = structureDescriptors[structureIdx].atomDescriptors[atomIdx].planarRingsIndices;
                for (auto ringIdx : structureDescriptors[structureIdx].atomDescriptors[atomIdx].planarRingsIndices)
                    ringsSizes.push_back(structureDescriptors[structureIdx].planarRings[ringIdx].size());
                sort(ringsSizes.begin(), ringsSizes.end());
                for(auto ringSize: ringsSizes)
                {
                    if (!ringsSizesStr.empty())
                        ringsSizesStr += ",";
                    ringsSizesStr += to_string(ringSize);
                }
            }

            out << setw(11) << ringsSizesStr;

            // 3 and 4 member rings

            out << " , in " << structureDescriptors[structureIdx].atomDescriptors[atomIdx].n3memberRings << " 3-member rings "
                << " , in " << structureDescriptors[structureIdx].atomDescriptors[atomIdx].n4memberRings << " 4-member rings ";

            // collect information on unassigned atoms

            if (typeIndices[structureIdx][atomIdx] < 0)
            {
                // atomic number, neighbours, planarity, rings
                std::get<0>(unassignedAtomDescriptors) = structureDescriptors[structureIdx].atomDescriptors[atomIdx].atomicNumber;
                std::get<1>(unassignedAtomDescriptors) = formulaStr;
                std::get<2>(unassignedAtomDescriptors) = planarityStr;
                std::get<3>(unassignedAtomDescriptors) = ringsSizesStr;
                std::get<4>(unassignedAtomDescriptors) = structureDescriptors[structureIdx].atomDescriptors[atomIdx].n3memberRings;
                std::get<5>(unassignedAtomDescriptors) = structureDescriptors[structureIdx].atomDescriptors[atomIdx].n4memberRings;
                std::get<6>(unassignedAtomDescriptors) = secondNeighboursLabel(structureDescriptors[structureIdx], atomIdx);

                string atomLabel = structureIds[structureIdx] + string(",") + labels[atomIdx];
                sortedUnassignedAtoms[unassignedAtomDescriptors].push_back(atomLabel);
            }

            // ---------------------------------------

            out << " , neighbours ";
            const vector<int> &neighbourIndices = structureDescriptors[structureIdx].connectivity[atomIdx];
            vector<string> neighborsLabels;

            for (size_t nIdx = 0; nIdx < neighbourIndices.size(); nIdx++)
                neighborsLabels.push_back(structureDescriptors[structureIdx].atomDescriptors[neighbourIndices[nIdx]].label);
            sort(neighborsLabels.begin(), neighborsLabels.end());

            for (size_t nIdx = 0; nIdx < neighborsLabels.size(); nIdx++)
            {
                if (nIdx > 0)
                    out << ",";
                out << neighborsLabels[nIdx];
            }

            vector<string> neighborsTypes;
            out << " , neighbours types ";

            for (size_t nIdx = 0; nIdx < neighbourIndices.size(); nIdx++)
            {
                size_t neighbourIdx = neighbourIndices[nIdx];
                string typeName;
                if (neighbourIdx >= nAtomsInAsymmetricUnit[structureIdx])
                    typeName = "?";
                else
                    typeName = (typeIndices[structureIdx][neighbourIdx] >= 0 ? types[typeIndices[structureIdx][neighbourIdx]].id : string("-"));
                neighborsTypes.push_back( typeName );
            }

            sort(neighborsTypes.begin(), neighborsTypes.end());

            for (size_t nIdx = 0; nIdx < neighbourIndices.size(); nIdx++)
            {
                if (nIdx > 0)
                    out << ",";
                out << neighborsTypes[nIdx];
            }

            out << endl;
        } // for(size_t atomIdx = 0; atomIdx < structreDescriptors[moleculeIdx].atomDescriptors.size(); atomIdx++)

        out << "\n\n Atom Type Assignemnt\n\n";
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            string typeName = (typeIndices[structureIdx][atomIdx] >= 0 ? types[typeIndices[structureIdx][atomIdx]].id : string("-"));
            out << setw(8) << " " << labels[atomIdx] << " " << setw(8) << typeName << " ";


         // print lcs
            if (typeIndices[structureIdx][atomIdx] >= 0)
                out << lcs[structureIdx][atomIdx];// ubdbLcsAsString(lcs[structureIdx][atomIdx], labels);
            out << endl;

        }

    }//for (moleculeIdx = 0; moleculeIdx < nMolecules; moleculeIdx++)

}




void printSettingsInfo(
    const string &fName,
    const string &header,
    const Settings &settings)
{
    ofstream out(fName);

    if (!out.good())
        on_error::throwException(string("cannot print to file '") + fName + string("'"), __FILE__, __LINE__);

    string hash80(80, '#');
    hash80 += string("\n");
    
    out << hash80 << header << hash80 << "\n";

    
    out << "the following settings were used:"
        << "\n      ring considered not planar if it contains atoms with more than " << settings.descriptorsSettings.atomInRingMaxNeighbourCount << " neighbours"
        << "\n      ring considered not planar if it contains atoms with planarity above " << settings.descriptorsSettings.atomInRingPlanarityThreshold
        << "\n      ring considered not planar if its planarity is above " << settings.descriptorsSettings.ringPlanarityThreshold
        << "\n      atom considered planar if the planarity is below " << settings.descriptorsSettings.atomPlanarityThreshold
        << "\n      atoms are bonded if the interatomic distance is below sum of their covalent bonds plus " << settings.descriptorsSettings.covalentBondThreshold
        << "\n      minimal absolute value of P_lm to be included in the bank is " << settings.min_plm
        << "\n      P_lm is included in the bank if its absolute value is at least " << settings.nSigma << " times higher than the P_lm standard deviation"
        << "\n      minimal number of atom type instances for parameters update " << settings.min_n_instaces
        << "\n\n" << hash80 << "\n";
    out.close();
}


void printUnassignedAtomsInfo(
    const string fName,
    string header,
    map< UnassignedAtomDescriptors, vector<string> > &sortedUnassignedAtoms)
{
    ofstream out(fName);
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fName + string("'"), __FILE__, __LINE__);

    string hash80(80, '#');

    out << hash80 << "\n" << header << hash80 << "\n\n"
        << "       atoms with unassigned types sorted by descriptors\n\n" << hash80 << "\n\n";
    size_t i, n;
    for (auto &atoms : sortedUnassignedAtoms)
    {
        out << "atomic number " << setw(4) << get<0>(atoms.first)
            << " neighbours " << setw(12) << get<1>(atoms.first)
            << " planarity " << setw(15) << get<2>(atoms.first)
            << " planar rings with planar atoms " << setw(15) << get<3>(atoms.first) 
            << " ,3 member rings " << setw(3) << get<4>(atoms.first) 
            << " ,4 member rings " << setw(3) << get<5>(atoms.first) 
            << " , 2-nd neighbours " << get<6>(atoms.first) << "\n\n";
        
        n = atoms.second.size();
		out << " number of atoms : " << n << "\n";
        for (i = 0; i < n; i++)
        {
            out << atoms.second[i] << " ";
            if ((i + 1) % 8 == 0)
                out << "\n";
        }
        if (i % 8 != 0)
            out << "\n";
        out << "\n---------------------------------------------------------------------------\n\n";
    }
    
    out.close();
}

void printUnassignedAtomsInfoSortedByN(
	const string fName,
	string header,
	map< UnassignedAtomDescriptors, vector<string> > &sortedUnassignedAtoms)
{
	ofstream out(fName);
	if (!out.good())
		on_error::throwException(string("cannot write to log file '") + fName + string("'"), __FILE__, __LINE__);

	string hash80(80, '#');

	out << hash80 << "\n" << header << hash80 << "\n\n"
		<< "       atoms with unassigned types sorted by number of atoms in group\n\n" << hash80 << "\n\n";

	vector < pair< UnassignedAtomDescriptors, vector<string> > > unassignedAtomGroups;
	vector<pair<size_t, size_t> > groupMultiplicityAndIndex;

	size_t index = 0;
	for (auto &atoms : sortedUnassignedAtoms)
	{
		unassignedAtomGroups.push_back(atoms);
		groupMultiplicityAndIndex.push_back({atoms.second.size(), index});
		index++;
	}

	sort(groupMultiplicityAndIndex.begin(), groupMultiplicityAndIndex.end(), greater< pair<size_t, size_t> >());

	for (size_t i = 0; i < groupMultiplicityAndIndex.size(); i++)
	{
		index = groupMultiplicityAndIndex[i].second;

		out << "atomic number " << setw(4) << get<0>(unassignedAtomGroups[index].first) 
			<< " neighbours " << setw(12) << get<1>(unassignedAtomGroups[index].first)
			<< " planarity " << setw(15) << get<2>(unassignedAtomGroups[index].first)
			<< " planar rings with planar atoms " << setw(15) << get<3>(unassignedAtomGroups[index].first)
			<< " ,3 member rings " << setw(3) << get<4>(unassignedAtomGroups[index].first)
            << " ,4 member rings " << setw(3) << get<5>(unassignedAtomGroups[index].first)
            << " , 2-nd neighbours " << get<6>(unassignedAtomGroups[index].first) << "\n\n";

		size_t n = unassignedAtomGroups[index].second.size();
		out << " number of atoms : " << n << "\n";
		size_t j;
		for (j = 0; j < n; j++)
		{
			out << setw(20) << unassignedAtomGroups[index].second[j] << " ";
			if ((j + 1) % 8 == 0)
				out << "\n";
		}
		if ( j % 8 != 0)
			out << "\n";
		out << "\n---------------------------------------------------------------------------\n\n";
	}


	out.close();
}


void printUnassignedAtomsInfoSortedByNStructures(
    const string fName,
    string header,
    map< UnassignedAtomDescriptors, vector<string> >& sortedUnassignedAtoms)
{
    ofstream out(fName);
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fName + string("'"), __FILE__, __LINE__);

    string hash80(80, '#');

    out << hash80 << "\n" << header << hash80 << "\n\n"
        << "       atoms with unassigned types sorted by number of atoms in group\n\n" << hash80 << "\n\n";

    vector < pair< UnassignedAtomDescriptors, vector<string> > > unassignedAtomGroups;
    vector<pair<size_t, size_t> > groupMultiplicityAndIndex;

    size_t index = 0;
    for (auto& atoms : sortedUnassignedAtoms)
    {
        unassignedAtomGroups.push_back(atoms);
        set<string> uniqueStructures;
        for (auto& label : atoms.second)
        {
            vector<string> words;
            string_utilities::split(label, words, ',');
            uniqueStructures.insert(words[0]);
        }
        groupMultiplicityAndIndex.push_back({ uniqueStructures.size(), index });
        index++;
    }

    sort(groupMultiplicityAndIndex.begin(), groupMultiplicityAndIndex.end(), greater< pair<size_t, size_t> >());

    for (size_t i = 0; i < groupMultiplicityAndIndex.size(); i++)
    {
        index = groupMultiplicityAndIndex[i].second;

        out << "atomic number " << setw(4) << get<0>(unassignedAtomGroups[index].first)
            << " neighbours " << setw(12) << get<1>(unassignedAtomGroups[index].first)
            << " planarity " << setw(15) << get<2>(unassignedAtomGroups[index].first)
            << " planar rings with planar atoms " << setw(15) << get<3>(unassignedAtomGroups[index].first)
            << " ,3 member rings " << setw(3) << get<4>(unassignedAtomGroups[index].first)
            << " ,4 member rings " << setw(3) << get<5>(unassignedAtomGroups[index].first)
            << " , 2-nd neighbours " << get<6>(unassignedAtomGroups[index].first) << "\n\n";

        size_t n = unassignedAtomGroups[index].second.size();
        out << " number of structures : " << groupMultiplicityAndIndex[i].first << "\n";
        out << " number of atoms : " << n << "\n";
        size_t j;
        for (j = 0; j < n; j++)
        {
            out << setw(20) << unassignedAtomGroups[index].second[j] << " ";
            if ((j + 1) % 8 == 0)
                out << "\n";
        }
        if (j % 8 != 0)
            out << "\n";
        out << "\n---------------------------------------------------------------------------\n\n";
    }


    out.close();
}



void printUnassignedAtomsInfoSortedByNStructuresAndFormula(
    const string fName,
    string header,
    map< UnassignedAtomDescriptors, vector<string> >& _sortedUnassignedAtoms,
    bool hydrogen_only = false,
    bool hydrogen_bonded_only = false)
{
    if(hydrogen_only && hydrogen_bonded_only)
        on_error::throwException("printUnassignedAtomsInfoSortedByNStructuresAndFormula: "
            "only one of hydrogen_only and hydrogen_bonded_only can be set true", __FILE__, __LINE__);

    if (hydrogen_only || hydrogen_bonded_only)
    {
        map< UnassignedAtomDescriptors, vector<string> > filteredSortedUnassignedAtoms;
        if (hydrogen_only)
        {
            for (auto& item : _sortedUnassignedAtoms)
                if (get<0>(item.first) == 1) // hydrogen
                    filteredSortedUnassignedAtoms.insert(item);
        }
        if (hydrogen_bonded_only)
        {
            for (auto& item : _sortedUnassignedAtoms)
            {
                string neighbours = get<1>(item.first);
                if (neighbours.find("H")!=string::npos) // bonded to hydrogen
                    filteredSortedUnassignedAtoms.insert(item);
            }
        }
        printUnassignedAtomsInfoSortedByNStructuresAndFormula(fName, header, filteredSortedUnassignedAtoms);
        return;
    }

    ofstream out(fName);
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fName + string("'"), __FILE__, __LINE__);

    map<string, vector<string> > sortedUnassignedAtoms;

    for (auto& item : _sortedUnassignedAtoms)
    {
        string label = periodic_table::symbol(get<0>(item.first)) + " " + get<1>(item.first);
        for (auto& str : item.second)
            sortedUnassignedAtoms[label].push_back(str);
    }

    string hash80(80, '#');

    out << hash80 << "\n" << header << hash80 << "\n\n"
        << "       atoms with unassigned types sorted by number of atoms in group\n\n" << hash80 << "\n\n";

    vector < pair< string, vector<string> > > unassignedAtomGroups;
    vector<pair<int, int> > groupMultiplicityAndIndex;

    size_t index = 0;
    for (auto& atoms : sortedUnassignedAtoms)
    {
        unassignedAtomGroups.push_back(atoms);
        set<string> uniqueStructures;
        for (auto& label : atoms.second)
        {
            vector<string> words;
            string_utilities::split(label, words, ',');
            uniqueStructures.insert(words[0]);
        }
        groupMultiplicityAndIndex.push_back({ uniqueStructures.size(), index });
        index++;
    }

    sort(groupMultiplicityAndIndex.begin(), groupMultiplicityAndIndex.end(), greater< pair<size_t, size_t> >());

    for (size_t i = 0; i < groupMultiplicityAndIndex.size(); i++)
    {
        index = groupMultiplicityAndIndex[i].second;

        out << unassignedAtomGroups[index].first << "\n\n";

        size_t n = unassignedAtomGroups[index].second.size();
        out << " number of structures : " << groupMultiplicityAndIndex[i].first << "\n";
        out << " number of atoms : " << n << "\n";
        size_t j;
        sort(unassignedAtomGroups[index].second.begin(), unassignedAtomGroups[index].second.end());
        for (j = 0; j < n; j++)
        {
            out << setw(20) << unassignedAtomGroups[index].second[j] << " ";
            if ((j + 1) % 8 == 0)
                out << "\n";
        }
        if (j % 8 != 0)
            out << "\n";
        out << "\n---------------------------------------------------------------------------\n\n";
    }


    out.close();
}

void printUnassignedFormulaStats(
    const string fName,
    string header,
    map< UnassignedAtomDescriptors, vector<string> >& _sortedUnassignedAtoms)
{
    ofstream out(fName);
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fName + string("'"), __FILE__, __LINE__);

    map<string, vector<string> > sortedUnassignedAtoms;

    for (auto& item : _sortedUnassignedAtoms)
    {
        string label = periodic_table::symbol(get<0>(item.first)) + " " + get<1>(item.first);
        for (auto& str : item.second)
            sortedUnassignedAtoms[label].push_back(str);
    }

    string hash80(80, '#');

    out << hash80 << "\n" << header << hash80 << "\n\n"
        << "       atoms with unassigned types sorted by number of atoms in group\n\n" << hash80 << "\n\n";

    vector < pair< string, vector<string> > > unassignedAtomGroups;
    vector<pair<int, int> > groupMultiplicityAndIndex;

    size_t index = 0;
    for (auto& atoms : sortedUnassignedAtoms)
    {
        unassignedAtomGroups.push_back(atoms);
        set<string> uniqueStructures;
        for (auto& label : atoms.second)
        {
            vector<string> words;
            string_utilities::split(label, words, ',');
            uniqueStructures.insert(words[0]);
        }
        groupMultiplicityAndIndex.push_back({ uniqueStructures.size(), index });
        index++;
    }

    sort(groupMultiplicityAndIndex.begin(), groupMultiplicityAndIndex.end(), greater< pair<size_t, size_t> >());

    for (size_t i = 0; i < groupMultiplicityAndIndex.size(); i++)
    {
        index = groupMultiplicityAndIndex[i].second;

        out << unassignedAtomGroups[index].first          
            << " n str. " << groupMultiplicityAndIndex[i].first 
            << " n atoms " << unassignedAtomGroups[index].second.size() << "\n";
    }


    out.close();
}


void generateVersionInfo(string &header)
{
    string versionString = discamb_version::version();
    time_t time_now = std::chrono::system_clock::to_time_t(chrono::system_clock::now());

    header = string("\n\n file generated with typeQuest, version ") + versionString + string("\n") +
        string("      compiled at ") + string(__DATE__) + string(", ") + string(__TIME__) + string(".\n") +
        string("      executed at ") + string(ctime(&time_now)) + string("\n\n\n");

    cout << "\n  typeQuest, version " << versionString << "\n"
        << "  compiled at " << __DATE__ << ", " << __TIME__ << ".\n\n";

}
void read_cif(
    const string &fileName,
    size_t &nAtomsInAsymmetricUnit,
    StructureWithDescriptors &structure)
{
    vector<cif_io::DataSet> dataSets;
    cif_io::readCif(fileName, dataSets);
    Crystal crystal;
    cif_io::cifDataToCrystal(dataSets[0], crystal);
    nAtomsInAsymmetricUnit = crystal.atoms.size();
    vector<int> atomicNumbers;
    vector<string> labels;
    vector<Vector3d> positions;
    structural_properties::assymetricUnitWithNeighbours(crystal, atomicNumbers, positions, labels, 8, structure.settings.covalentBondThreshold);
    structure.set(atomicNumbers, positions, labels);
}

void read_shelx(
    const string &fileName,
    size_t &nAtomsInAsymmetricUnit,
    StructureWithDescriptors &structure)
{
    Crystal crystal;
    shelx_io::read(fileName, crystal);
    nAtomsInAsymmetricUnit = crystal.atoms.size();
    vector<int> atomicNumbers;
    vector<string> labels;
    vector<Vector3d> positions;
    structural_properties::assymetricUnitWithNeighbours(crystal, atomicNumbers, positions, labels, 8, structure.settings.covalentBondThreshold);
    structure.set(atomicNumbers, positions, labels);
}

bool readStructure(
    const vector<string> &files,
    const std::string &format,
    const CrystalAtomTypeAssigner &crystalAssigner,
    size_t atomTypeRange,
    size_t namedNeighboursRange,
    size_t &nAtomsInAsymmetricUnit,
    StructureWithDescriptors &structure,
    string &error_message)
{
    try {

        if (files.empty())
        {
            error_message = "no input file";
            return false;
        }

        vector<int> atomicNumbers;
        vector<string> labels;
        vector<Vector3d> positions;


        if (format == string("xyz"))
        {
            xyz_io::readXyz(files[0], atomicNumbers, positions);
            for (size_t i = 0; i < atomicNumbers.size(); i++)
                labels.push_back(periodic_table::symbol(atomicNumbers[i]) + to_string(i + 1));
            nAtomsInAsymmetricUnit = atomicNumbers.size();
            structure.set(atomicNumbers, positions, labels);
        }
        else
        {
            Crystal crystal;

            if (format == string("shelx"))
                shelx_io::read(files[0], crystal);
            if (format == string("cif"))
            {
                vector<cif_io::DataSet> dataSets;
                cif_io::readCif(files[0], dataSets);
                cif_io::cifDataToCrystal(dataSets[0], crystal);
            }
            if (format == string("xd"))
            {
                HC_ModelParameters hcParameters;
                vector<XdLocalCoordinateSystem> xdLcs;
                xd_io::read(files[0], files[1], hcParameters, crystal, xdLcs, false);
            }

            nAtomsInAsymmetricUnit = crystal.atoms.size();
            vector<int> typeIds;
            vector<LocalCoordinateSystem<AtomInCrystalID> > lcs;
            crystalAssigner.assign(crystal, typeIds, lcs, structure);
        }
    }
    catch (exception &e)
    {
        error_message = e.what();
        return false;
    }
    return true;
}

//bool readStructureAndAssign(
//    const vector<string> &files,
//    const std::string &format,
//    const CrystalAtomTypeAssigner &crystalAssigner,
//    const MolecularAtomTypeAssigner &molecularAssigner,
//    size_t atomTypeRange,
//    size_t namedNeighboursRange,
//    size_t &nAtomsInAsymmetricUnit,
//    StructureWithDescriptors &structure,
//    vector<int> &typeIds,
//    vector<string> &lcs,
//    string &error_message)
//{
//    try {
//
//        if (files.empty())
//        {
//            error_message = "no input file";
//            return false;
//        }
//
//        lcs.clear();
//
//        vector<int> atomicNumbers;
//        vector<string> labels;
//        vector<Vector3d> positions;
//        bool processingCrystalStructure = true;
//
//        if (format == string("xyz") || format == string("mol2"))
//        {
//            if (format == string("xyz"))
//            {
//                xyz_io::readXyz(files[0], atomicNumbers, positions);
//                for (size_t i = 0; i < atomicNumbers.size(); i++)
//                    labels.push_back(periodic_table::symbol(atomicNumbers[i]) + to_string(i + 1));
//            }
//            else
//            {
//                mol2_io::Mol2Data mol2data;
//                mol2_io::read(files[0], mol2data);
//                mol2_io::atomicNumbers(mol2data, atomicNumbers);
//                positions = mol2data.atomPosition;
//                labels = mol2data.atomName;
//            }
//
//            nAtomsInAsymmetricUnit = atomicNumbers.size();
//            structure.set(atomicNumbers, positions, labels);
//            vector<LocalCoordinateSystem<int> > lcsMolecule;
//            molecularAssigner.assign(structure, typeIds, lcsMolecule);
//            
//            for (auto &coordinateSystem : lcsMolecule)
//                lcs.push_back(ubdbLcsAsString(coordinateSystem, labels));
//        }
//        else
//        {
//            Crystal crystal;
//
//            if (format == string("shelx"))
//                shelx_io::read(files[0], crystal);
//            if (format == string("cif"))
//            {
//                vector<cif_io::DataSet> dataSets;
//                cif_io::readCif(files[0], dataSets);
//                cif_io::cifDataToCrystal(dataSets[0], crystal);
//            }
//            if (format == string("xd"))
//            {
//                HC_ModelParameters hcParameters;
//                vector<XdLocalCoordinateSystem> xdLcs;
//                xd_io::read(files[0], files[1], hcParameters, crystal, xdLcs, false);
//            }
//
//            nAtomsInAsymmetricUnit = crystal.atoms.size();
//            
//            vector<LocalCoordinateSystem<AtomInCrystalID> > lcsCrystal;
//            crystalAssigner.assign(crystal, typeIds, lcsCrystal, structure);
//
//            for (auto &atom : crystal.atoms)
//                labels.push_back(atom.label);
//
//            for (auto &coordinateSystem : lcsCrystal)
//                lcs.push_back(ubdbLcsAsString(coordinateSystem, labels));
//        }
//    }
//    catch (exception &e)
//    {
//        error_message = e.what();
//        return false;
//    }
//    return true;
//}

bool readStructureAndAssign(
    const vector<string>& files,
    const std::string& format,
    const CrystalAtomTypeAssigner& crystalAssigner,
    const MolecularAtomTypeAssigner& molecularAssigner,
    size_t atomTypeRange,
    size_t namedNeighboursRange,
    size_t& nAtomsInAsymmetricUnit,
    StructureWithDescriptors& structure,
    vector<int>& typeIds,
    vector<vector<int> >& multiTypeIds,
    vector<string>& lcs,
    string& error_message)
{
    try {

        if (files.empty())
        {
            error_message = "no input file";
            return false;
        }

        lcs.clear();

        vector<int> atomicNumbers;
        vector<string> labels;
        vector<Vector3d> positions;
        bool processingCrystalStructure = true;

        if (format == string("xyz") || format == string("mol2"))
        {
            if (format == string("xyz"))
            {
                xyz_io::readXyz(files[0], atomicNumbers, positions);
                for (size_t i = 0; i < atomicNumbers.size(); i++)
                    labels.push_back(periodic_table::symbol(atomicNumbers[i]) + to_string(i + 1));
            }
            else
            {
                mol2_io::Mol2Data mol2data;
                mol2_io::read(files[0], mol2data);
                mol2_io::atomicNumbers(mol2data, atomicNumbers);
                positions = mol2data.atomPosition;
                labels = mol2data.atomName;
            }

            nAtomsInAsymmetricUnit = atomicNumbers.size();
            structure.set(atomicNumbers, positions, labels);
            vector<LocalCoordinateSystem<int> > lcsMolecule;
            molecularAssigner.assign(structure, typeIds, lcsMolecule);
            molecularAssigner.assign_all_possible(structure, multiTypeIds);

            for (auto& coordinateSystem : lcsMolecule)
                lcs.push_back(ubdbLcsAsString(coordinateSystem, labels));
        }
        else
        {
            Crystal crystal;

            if (format == string("shelx"))
                shelx_io::read(files[0], crystal);
            if (format == string("cif"))
            {
                vector<cif_io::DataSet> dataSets;
                cif_io::readCif(files[0], dataSets);
                cif_io::cifDataToCrystal(dataSets[0], crystal);
            }
            if (format == string("xd"))
            {
                HC_ModelParameters hcParameters;
                vector<XdLocalCoordinateSystem> xdLcs;
                xd_io::read(files[0], files[1], hcParameters, crystal, xdLcs, false);
            }

            nAtomsInAsymmetricUnit = crystal.atoms.size();

            vector<LocalCoordinateSystem<AtomInCrystalID> > lcsCrystal;
            crystalAssigner.assign(crystal, typeIds, lcsCrystal, structure);
            crystalAssigner.assign_all_possible(crystal, multiTypeIds);

            for (auto& atom : crystal.atoms)
                labels.push_back(atom.label);

            for (auto& coordinateSystem : lcsCrystal)
                lcs.push_back(ubdbLcsAsString(coordinateSystem, labels));
        }
    }
    catch (exception& e)
    {
        error_message = e.what();
        return false;
    }
    return true;
}

void findMultitypes(
    const vector<AtomType> &types,
    vector<pair<string, vector<size_t> > > &multitypes)
{
    multitypes.clear();
    map<string, vector<size_t> > mtypes;

    for (size_t i = 0; i < types.size(); i++)
        mtypes[types[i].id].push_back(i);

    for (auto mtype : mtypes)
        multitypes.push_back({ mtype.first, mtype.second });
}

void findTypeInstances(
    const vector<pair<string, vector<size_t> > > &multitypes,
    const vector<vector<int> > &typeIndices,
    vector< vector<vector<pair<size_t, size_t> > > > &typeInstances,
    size_t total_n_subtypes)
{
    size_t multitypeIdx, subtypeIdx, nSubtypes, nMultiTypes = multitypes.size();
    typeInstances.clear();
    typeInstances.resize(nMultiTypes);
    for (multitypeIdx = 0; multitypeIdx < nMultiTypes; multitypeIdx++)
        typeInstances[multitypeIdx].resize(multitypes[multitypeIdx].second.size());

    vector<pair<size_t, size_t> > type2multitypeSubtype(total_n_subtypes);

    for (multitypeIdx = 0; multitypeIdx < nMultiTypes; multitypeIdx++)
    {
        nSubtypes = multitypes[multitypeIdx].second.size();

        for (subtypeIdx = 0; subtypeIdx < nSubtypes; subtypeIdx++)
            type2multitypeSubtype[multitypes[multitypeIdx].second[subtypeIdx]] = { multitypeIdx, subtypeIdx };
    }
    //
    size_t structureIdx, nStructures, atomIdx, nAtoms;

    nStructures = typeIndices.size();

    for (structureIdx = 0; structureIdx < nStructures; structureIdx++)
    {
        nAtoms = typeIndices[structureIdx].size();
        for(atomIdx=0;atomIdx<nAtoms; atomIdx++)
            if (typeIndices[structureIdx][atomIdx] >= 0)
            {
                multitypeIdx = type2multitypeSubtype[typeIndices[structureIdx][atomIdx]].first;
                subtypeIdx = type2multitypeSubtype[typeIndices[structureIdx][atomIdx]].second;
                typeInstances[multitypeIdx][subtypeIdx].push_back({ structureIdx, atomIdx });
            }
    }
    
}

void getFiles(
    const std::string &fileFormat,
    vector<vector<string> > &files,
    vector<string> &structureIds,
    const string &structuresFolder,
    const string &path_part_2)
{   
    vector<string> formats{ "cif","shelx","xyz","mol2","xd" };

    if (find(formats.begin(), formats.end(), fileFormat) == formats.end())
        on_error::throwException(string("invalid file format: '") + fileFormat + string(" expected one of the following: cif, shelx, xyz, xd"), __FILE__, __LINE__);

    files.clear();

    if (fileFormat == string("xd"))
    {
        vector<pair<string, string> > xdFilePath;
        getXdFilePaths(xdFilePath, structureIds, structuresFolder, path_part_2);
        for (auto &f : xdFilePath)
            files.push_back({ f.first,f.second });
    }
    else
    {
        vector<string> structureFiles;
        map<string, string> suffix{ {"cif",".cif"},{"xyz",".xyz"},{"mol2",".mol2"},{"shelx",".res"} };
        getFilePaths(structuresFolder, structureFiles, structureIds, suffix[fileFormat]);
        for (auto &f : structureFiles)
            files.push_back({ f });
    }

}

void printMultiTypeInfo(
    const string& fileName,
    const vector< vector<vector<int> > >& typeIndices,
    const vector<AtomType>& types)
{
    ofstream out(fileName);
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fileName + string("'"), __FILE__, __LINE__);
    size_t nTypes = types.size();
    size_t nStructures = typeIndices.size();
    vector<size_t> typeCounts(nTypes, 0);

    for (size_t structureIdx = 0; structureIdx < nStructures; structureIdx++)
    {
        size_t nAtoms = typeIndices[structureIdx].size();
        for (size_t atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            for (int i = 0; i < typeIndices[structureIdx][atomIdx].size(); i++)
            {
                int typeIdx = typeIndices[structureIdx][atomIdx][i];
                if (typeIdx >= 0)
                    typeCounts[typeIdx]++;
            }
        }
    }
    vector<pair<int, string> > countsAndTypeId;
    int typeIdMaxLength = 0;
    for (int i = 0; i < nTypes; i++)
    {
        countsAndTypeId.push_back({ (int)typeCounts[i], types[i].id });
        if (typeIdMaxLength < types[i].id.size())
            typeIdMaxLength = types[i].id.size();
    }
    sort(countsAndTypeId.begin(), countsAndTypeId.end(), greater< pair<int, string> >());

    out << "Type ID               Number of assigned atoms\n";
    out << "----------------------------------------------\n";
    for (size_t typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        out << setw(typeIdMaxLength + 2) << countsAndTypeId[typeIdx].second << " " << setw(10) << countsAndTypeId[typeIdx].first << "\n";
        //out << setw(20) << types[typeIdx].id << " " << setw(10) << typeCounts[typeIdx] << "\n";
    }
    out.close();

}

void printTypeInfo(
    const string& fileName,
    const vector<vector<int> >& typeIndices,
    const vector<AtomType>& types)
{
    ofstream out(fileName);
    if (!out.good())
        on_error::throwException(string("cannot write to log file '") + fileName + string("'"), __FILE__, __LINE__);
    size_t nTypes = types.size();
    size_t nStructures = typeIndices.size();
    vector<size_t> typeCounts(nTypes, 0);

    for (size_t structureIdx = 0; structureIdx < nStructures; structureIdx++)
    {
        size_t nAtoms = typeIndices[structureIdx].size();
        for (size_t atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int typeIdx = typeIndices[structureIdx][atomIdx];
            if (typeIdx >= 0)
                typeCounts[typeIdx]++;
        }
    }
    vector<pair<int, string> > countsAndTypeId;
    int typeIdMaxLength = 0;
    for (int i = 0; i < nTypes; i++)
    {
        countsAndTypeId.push_back({ (int)typeCounts[i], types[i].id });
        if(typeIdMaxLength< types[i].id.size())
            typeIdMaxLength = types[i].id.size();
    }
    sort(countsAndTypeId.begin(), countsAndTypeId.end(), greater< pair<int, string> >());

    out << "Type ID               Number of assigned atoms\n";
    out << "----------------------------------------------\n";
    for (size_t typeIdx = 0; typeIdx < nTypes; typeIdx++)
    {
        out << setw(typeIdMaxLength + 2) << countsAndTypeId[typeIdx].second << " " << setw(10) << countsAndTypeId[typeIdx].first << "\n";
        //out << setw(20) << types[typeIdx].id << " " << setw(10) << typeCounts[typeIdx] << "\n";
    }
    out.close();

}

int main(int argc, char *argv[])
{

    try {
        string versionString = discamb_version::version();
        time_t time_now = std::chrono::system_clock::to_time_t(chrono::system_clock::now());

        string header = string("\n\n file generated with typeStat, version ") + versionString + string("\n") +
                        string("      compiled at ") + string(__DATE__) + string(", ") + string(__TIME__) + string(".\n") +
                        string("      executed at ") + string(ctime(&time_now)) + string("\n\n\n");
    
        std::cout<< "\n  typeStat, version " << versionString << "\n"
            << "  compiled at " << __DATE__ << ", " << __TIME__ << ".\n\n";

        //Settings settings;
        
        string fileFormat = "shelx";
        string fileFolder;
        string bank_file;
        vector<string> arguments, options;
        map<string, string> optionsWithValues;
        parse_cmd::get_args_and_options(argc, argv, arguments, options, optionsWithValues);

        int nCores = 1;
        if (optionsWithValues.count("-n")>0)
            nCores = stoi(optionsWithValues["-n"]);

        bool do_not_write_assignment_info = parse_cmd::hasOption(options, "-no_ai");
        bool do_not_write_type_instances_detailed_info = parse_cmd::hasOption(options, "-no_tidi");
        bool do_not_write_type_instances_info = parse_cmd::hasOption(options, "-no_tii");
        if( parse_cmd::hasOption(options, "-no_bulky") )
        {
            do_not_write_assignment_info = true;
            do_not_write_type_instances_detailed_info = true;
            do_not_write_type_instances_info = true;
        }

        if (arguments.size() == 0)
        {
            string error_message = 
                "missing arguments, expected name of databank file and optionally\n" 
                "file format and name of folder with structure files, \n"
                "file formats can be cif, shelx, xd, mol2 or xyz - shelx is default\n"
                "possible options:\n"
                "-no_ai - assignment_info.txt is not written\n"
                "-no_tidi - type_instances_detailed_info.txt is not written\n"
                "-no_tii - type_instances_info.txt is not written\n"
                "-no_bulky - the three file mentioned above are not written\n";
            on_error::throwException(error_message, __FILE__, __LINE__);
        }
        else {
            bank_file = arguments[0];
            if (arguments.size() > 1)
                fileFormat = arguments[1];
            if (argc > 3)
                fileFolder = arguments[2];
        }

        vector<Crystal> crystals;
        vector<string> structureIds;
        vector<AtomType> types;


        vector<MolecularAtomTypeAssigner> molecularAssigners(nCores);
        vector<CrystalAtomTypeAssigner> crystalAssigners(nCores);

        BankSettings bankSettings;
        size_t nConsideredStructures, nCompletelyRecognizedStructures;

        //vector<vector<Crystal> > crystals_in_thread(nCores);
        //vector< vector<string> > structureIds_in_thread(nCores);
        vector<size_t> nConsideredStructures_in_thread(nCores), nCompletelyRecognizedStructures_in_thread(nCores);
        vector<vector<string> > validStructureIds_in_thread(nCores);
        vector< vector<StructureWithDescriptors> > structureDescriptors_in_thread(nCores);
        vector< vector < vector<string> > > lcs_in_thread(nCores);
        vector<vector<vector<int> > > typeIndices_in_thread(nCores);
        vector< vector<vector<vector<int> > > > multiTypeIndices_in_thread(nCores);
        vector< vector<size_t> > natomsInAsymmetricUnit_in_thread(nCores);
        vector< vector<string> > fileIoErrorMessages_in_thread(nCores);
        vector< vector<string> > invalidStructureReads_in_thread(nCores);

        /*
                
                        vector<vector<int> > typeIndices;    
        vector<size_t> natomsInAsymmetricUnit;
        vector < vector<string> > lcs;
        vector<vector<string> > filePaths;

                
                
                typeIndices.push_back(oneStructureTypeIndices);
                nConsideredStructures++;        
        */


        DescriptorsSettings descriptorsSettings;
        atom_type_io::readAtomTypes(bank_file, types, descriptorsSettings);

        vector<vector<int> > typesGeneralized;
        vector<vector<int> > hierarchyLevel;
        vector<pair<int, int> > equivalentTypes;
        vector<pair<int, vector<int> > > generalizedByMultipleTypesAtSameLevel;

        atom_typing_utilities::typeGeneralizationDiagnostics(
            typesGeneralized,
            hierarchyLevel,
            equivalentTypes,
            generalizedByMultipleTypesAtSameLevel);

        if (!equivalentTypes.empty())
        {
            cout << "Error in type bank: there are equivalent types, info written to equivalent_types.txt\n";
            ofstream equivalent_types_file("equivalent_types.txt");
            equivalent_types_file << "Equivalent atom types:\n";
            for (auto& type_pair : equivalentTypes)
                equivalent_types_file << types[type_pair.first].id << "   " << types[type_pair.second].id << "\n";
            equivalent_types_file.close();
            return 0;
        }
        
        //atom_typing_utilities::sortTypesByGenarality_LevelsAbove(types);
        atom_typing_utilities::sortTypesByGenarality_LevelsBelow(types);

        atom_type_io::writeAtomTypesTxt("atom_types.log", types);
        //bankReader.read(argv[1], types, typeHcParameters, bankSettings);
    
        cout << types.size() << " types\n";

        for (int i = 0; i < nCores; i++)
        {
            molecularAssigners[i].setDescriptorsSettings(descriptorsSettings);
            molecularAssigners[i].setAtomTypes(types);

            crystalAssigners[i].setDescriptorsSettings(bankSettings.descriptorsSettings);
            crystalAssigners[i].setAtomTypes(types);

        }

        //molecularAssigner.setDescriptorsSettings(descriptorsSettings);
        //molecularAssigner.setAtomTypes(types);
        //
        //crystalAssigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
        //crystalAssigner.setAtomTypes(types);

        int namedAtomsRange, planarRingRange, ring34Range, atomTypeRange;
        atomTypeRange = atom_typing_utilities::atomTypesRange(types, bankSettings.descriptorsSettings.maxPlanarRing, 
                                                       namedAtomsRange, planarRingRange, ring34Range);
    
        vector<vector<int> > typeIndices;    
        vector<vector<vector<int> > > multiTypeIndices;
        vector<size_t> natomsInAsymmetricUnit;
        vector < vector<string> > lcs;
        vector<vector<string> > filePaths;

        getFiles(fileFormat, filePaths, structureIds, fileFolder, string());
        

        int nStructures = filePaths.size();
        

        LocalCoordinateSystemCalculator lcs_calculator;
        vector<StructureWithDescriptors> structureDescriptors;


        //-------
        
        //Crystal crystal;
        //----------


        vector<size_t> atomToAssign;
        vector<string> fileIoErrorMessages;
        vector<string> invalidStructureReads;
        vector<string> validStructureIds;
        vector<vector<int> > structureIdxPerThread(nCores);
		nConsideredStructures = nCompletelyRecognizedStructures = 0;
        cout << "n threads = " << nCores << "\n";
        cout << "    processing:\n";
        omp_set_num_threads(nCores);


#pragma omp parallel for num_threads(nCores)        
        for (int structureIdx = 0; structureIdx < nStructures; structureIdx++)
        {
            string error_message;
            int thread_id = omp_get_thread_num();
            structureIdxPerThread[thread_id].push_back(structureIdx);
            if(thread_id==0)
                cout << "\r" << setw(15) << structureIds[structureIdx];

            StructureWithDescriptors oneStructureDescriptors;
            size_t nAtomsInOneStructureAsymmetricUnit;
            vector<int> oneStructureTypeIndices;
            vector<string> oneStructureLcs;
            vector<vector<int> > oneStructureMultiTypeIndices;

            oneStructureDescriptors.settings = bankSettings.descriptorsSettings;

            if(readStructureAndAssign(filePaths[structureIdx], fileFormat, crystalAssigners[thread_id],
                molecularAssigners[thread_id], atomTypeRange, namedAtomsRange,
                nAtomsInOneStructureAsymmetricUnit, oneStructureDescriptors, oneStructureTypeIndices,
                oneStructureMultiTypeIndices,oneStructureLcs, error_message))
            {
                /*
        
        
        
        vector<size_t> nConsideredStructures_in_thread(nCores), nCompletelyRecognizedStructures_in_thread(nCores);
        vector<vector<string> > validStructureIds_in_thread(nCores);
        vector< vector<StructureWithDescriptors> > structureDescriptors_in_thread(nCores);
        vector< vector < vector<string> > > lcs_in_thread(nCores);
        vector<vector<vector<int> > > typeIndices_in_thread(nCores);

                */
                //natomsInAsymmetricUnit.push_back(nAtomsInOneStructureAsymmetricUnit);
                multiTypeIndices_in_thread[thread_id].push_back(oneStructureMultiTypeIndices);
                natomsInAsymmetricUnit_in_thread[thread_id].push_back(nAtomsInOneStructureAsymmetricUnit);
                //validStructureIds.push_back(structureIds[structureIdx]);
                validStructureIds_in_thread[thread_id].push_back(structureIds[structureIdx]);
                structureDescriptors_in_thread[thread_id].push_back(oneStructureDescriptors);
                //structureDescriptors.push_back(oneStructureDescriptors);
                //lcs.push_back(oneStructureLcs);
                lcs_in_thread[thread_id].push_back(oneStructureLcs);
                typeIndices_in_thread[thread_id].push_back(oneStructureTypeIndices);
                nConsideredStructures_in_thread[thread_id]++;
				//nConsideredStructures++;
				
				if (find(oneStructureTypeIndices.begin(), oneStructureTypeIndices.end(), -1) == oneStructureTypeIndices.end())
                    nCompletelyRecognizedStructures_in_thread[thread_id]++;
                    //nCompletelyRecognizedStructures++;
            }
            else 
            {
                fileIoErrorMessages_in_thread[thread_id].push_back(error_message);
                invalidStructureReads_in_thread[thread_id].push_back(structureIds[structureIdx]);
            }

        }
        cout << "\n";
        for(int threadIdx=0; threadIdx<nCores; threadIdx++)
        {
            cout << "thread " << threadIdx << " processed " << structureIdxPerThread[threadIdx].size() << " structures\n";
            nConsideredStructures += nConsideredStructures_in_thread[threadIdx];
            nCompletelyRecognizedStructures += nCompletelyRecognizedStructures_in_thread[threadIdx];
            for(auto &id: validStructureIds_in_thread[threadIdx])
                validStructureIds.push_back(id);
            for(auto &desc: structureDescriptors_in_thread[threadIdx])
                structureDescriptors.push_back(desc);
            for(auto &lcs_one_thread: lcs_in_thread[threadIdx])
                lcs.push_back(lcs_one_thread);
            for(auto &ti: typeIndices_in_thread[threadIdx])
                typeIndices.push_back(ti);
            for(auto &natu: natomsInAsymmetricUnit_in_thread[threadIdx])
                natomsInAsymmetricUnit.push_back(natu);
            for(auto &msg: fileIoErrorMessages_in_thread[threadIdx])
                fileIoErrorMessages.push_back(msg);
            for(auto &invRead: invalidStructureReads_in_thread[threadIdx])
                invalidStructureReads.push_back(invRead);
            for(auto &mti: multiTypeIndices_in_thread[threadIdx])
                multiTypeIndices.push_back(mti);
        }

        validStructureIds.swap(structureIds);

        size_t nAssignedAtoms, nNotAssignedAtoms;
        nAssignedAtoms = nNotAssignedAtoms = 0;
        for (auto const &atomTypes : typeIndices)
            for (int typeIdx : atomTypes)
                if (typeIdx < 0)
                    nNotAssignedAtoms++;
                else
                    nAssignedAtoms++;


        vector< vector<vector<pair<size_t, size_t> > > > typeInstances;
        vector<pair<string, vector<size_t> > > multitypes; //[i].first - label [i].second - list of indices of types

        findMultitypes(types, multitypes);
        findTypeInstances(multitypes, typeIndices, typeInstances, types.size());

        if (!do_not_write_type_instances_info)
            printDetailedStats("type_instances_info.txt", header, multitypes, typeInstances, structureIds, structureDescriptors, nAssignedAtoms, nNotAssignedAtoms, nConsideredStructures, nCompletelyRecognizedStructures);

        if (!do_not_write_type_instances_detailed_info)
            printMoreDetailedDetailedStats("type_instances_detailed_info.txt", header, multitypes, typeInstances, structureIds, structureDescriptors);
        

        map< UnassignedAtomDescriptors, vector<string> > sortedUnassignedAtoms;

        if(!do_not_write_assignment_info)
            printAssignmentInfo("assignment_info.txt", structureIds, header, bankSettings.descriptorsSettings, structureDescriptors, typeIndices, types,
                        lcs, sortedUnassignedAtoms, natomsInAsymmetricUnit);

        printUnassignedAtomsInfo("unassigned_atoms_info.txt", header, sortedUnassignedAtoms);        
		printUnassignedAtomsInfoSortedByN("unassigned_atoms_sorted_by_n.txt", header, sortedUnassignedAtoms);
        printUnassignedAtomsInfoSortedByNStructures("unassigned_atoms_sorted_by_n_str.txt", header, sortedUnassignedAtoms);
        printUnassignedAtomsInfoSortedByNStructuresAndFormula("unassigned_atoms_sorted_by_n_str_and_formula.txt", header, sortedUnassignedAtoms);
        printUnassignedAtomsInfoSortedByNStructuresAndFormula("h_only_unassigned_atoms_sorted_by_n_str_and_formula.txt", header, sortedUnassignedAtoms, true);
        printUnassignedAtomsInfoSortedByNStructuresAndFormula("h_bonded_only_unassigned_atoms_sorted_by_n_str_and_formula.txt", header, sortedUnassignedAtoms, false, true);
        printUnassignedFormulaStats("unassigned_formula_stats.txt", header, sortedUnassignedAtoms);
        printStructureAssignedAndNot("structures_complete_assignemnet.txt", structureIds, typeIndices, structureDescriptors);
        printTypeInfo("types_stats.txt", typeIndices, types);
        printMultiTypeInfo("multitypes_stats.txt", multiTypeIndices, types);
        printMoveAssigned("move_assigned.bat", filePaths, typeIndices);

        cout << "\r";

        if (!fileIoErrorMessages.empty())
        {
            cout << "\n\n\n" << fileIoErrorMessages.size() << " structures have not been read sucessfully:\n\n";
            for (size_t i = 0; i < fileIoErrorMessages.size(); i++)
            {
                cout << invalidStructureReads[i] << "\n" << fileIoErrorMessages[i] << "\n\n";
                
            }
        }

        
        cout << "\n" << " DONE ";
        if (!fileIoErrorMessages.empty())
            cout << "\n\nNOTICE: some structures have not been read sucessfully, see above\n";
        cout << "\n";


    }
    catch (exception &e)
    {
        cout << e.what() << endl;
    }
}
