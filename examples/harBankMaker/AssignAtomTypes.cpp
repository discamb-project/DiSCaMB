#include "AssignAtomTypes.h"

#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/atom_typing_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/IO/atom_type_io.h"
#include "discamb/IO/mol2_io.h"
#include "discamb/IO/molecule_io.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "json.hpp"

#include <fstream>

using namespace discamb;
using namespace std;

AssignAtomTypes::AssignAtomTypes()
{

}

AssignAtomTypes::~AssignAtomTypes()
{

}

void AssignAtomTypes::readBank(
    const std::string& bankFile, 
    const std::string& selectedTypesFile, 
    std::vector< discamb::AtomType>& atomTypes,
    discamb::DescriptorsSettings& descriptorsSettings)
{
    if (bankFile.empty())
        on_error::throwException("missing \"bank file\" field in setting.json", __FILE__, __LINE__);

    vector<AtomType> types;
    atom_type_io::readAtomTypes(bankFile, types, descriptorsSettings);
    atom_typing_utilities::sortTypesByGenarality_LevelsBelow(types);


    // select atom types 
    


    if (selectedTypesFile.empty())
        atomTypes.swap(types);
    else
    {
        ifstream input(selectedTypesFile);
        string line;
        vector<string> words;
        vector<string> selectedTypeIds;
        map<string, int> type2Idx;
        for (int i = 0; i < types.size(); i++)
            type2Idx[types[i].id] = i;

        while (getline(input, line))
        {
            string_utilities::split(line, words);
            for (string& word : words)
            {
                auto it = type2Idx.find(word);
                if (it == type2Idx.end())
                    on_error::throwException("selected non-existing type, with id: '" + word + "'", __FILE__, __LINE__);
                else
                    atomTypes.push_back(types[it->second]);
            }

        }
        input.close();
    }

}

void AssignAtomTypes::set()
{
    nlohmann::json data;
    readSettings(data);

    // reads atom types
    string bankFile = data.value("bank file", string(""));
    if (bankFile.empty())
        on_error::throwException("expected \"bank file\" field in setting.json", __FILE__, __LINE__);

    string selectedAtomIdsFile = data.value("selected atom type ids file", string(""));

    readBank(bankFile, selectedAtomIdsFile, mAtomTypes, mDescriptorsSettings);

    //vector<AtomType> atomTypes;
    //atom_type_io::readAtomTypes(bankFile, atomTypes, mDescriptorsSettings);
    //
    //// select atom types 
    //
    //

    //if (selectedAtomIdsFile.empty())
    //    mAtomTypes.swap(atomTypes);
    //else
    //{
    //    input.open(selectedAtomIdsFile);
    //    string line;
    //    vector<string> words;
    //    vector<string> selectedTypeIds;
    //    map<string, int> type2Idx;
    //    for (int i = 0; i < atomTypes.size(); i++)
    //        type2Idx[atomTypes[i].id] = i;
    //    
    //    while (getline(input, line))
    //    {
    //        string_utilities::split(line, words);
    //        for (string& word : words)
    //        {
    //            auto it = type2Idx.find(word);
    //            if (it == type2Idx.end())
    //                on_error::throwException("selected non-existing type, with id: '" + word + "'", __FILE__, __LINE__);
    //            else
    //                mAtomTypes.push_back(atomTypes[it->second]);
    //        }

    //    }

    //}
    
    // sets mChosenResFolder
    
    string folder;
    if (mChosenResFolder.empty())
        folder = filesystem::current_path().string();
    else
        folder = mChosenResFolder.string();

    mChosenResFolder = data.value("chosen res folder", folder);

    // sets mMolFolder

    
    if (mMolFolder.empty())
        mMolFolder = filesystem::current_path().string();
    else
        mMolFolder = mMolFolder.string();

    mMolFolder = data.value("mol folder", folder);

    // sets mMinNumberOfInstances

    //mMinNumberOfInstances = data.value("min n type instances", mMinNumberOfInstances);

    // sets output

    mOutputFileName = data.value("type assignment log file", mOutputFileName);


}

void AssignAtomTypes::findFormulas(
    const mol2_io::Mol2Data& mol2Data,
    const std::vector<int>& atomicNumbers,
    std::vector<std::map<int, int> >& formulas)
{
    formulas.clear();

    vector<vector<int> > molecules;
    int nMolecules = *max_element(mol2Data.substructureIdx.begin(), mol2Data.substructureIdx.end());
    molecules.resize(nMolecules);
    formulas.resize(nMolecules);

    for (int i = 0; i < mol2Data.substructureIdx.size(); i++)
        molecules[mol2Data.substructureIdx[i] - 1].push_back(i);

    for (int molIdx = 0; molIdx < nMolecules; molIdx++)
    {
        vector<int> z;
        for (int idx : molecules[molIdx])
            z.push_back(atomicNumbers[idx]);
        basic_chemistry_utilities::getFormula(z, formulas[molIdx]);
    }
}

int AssignAtomTypes::findStructuralFormula(
    const StructuralFormula& sf,
    const std::vector<StructuralFormula>& sfs)
{
    for (int i = 0; i < sfs.size(); i++)
        if (structural_properties::molecularGraphIsomorphism(sf.z, sf.connectivity, sfs[i].z, sfs[i].connectivity))
            return i;
    return -1;
}

void AssignAtomTypes::findSubstructures(
    const discamb::mol2_io::Mol2Data& mol2Data,
    const std::vector<int>& atomicNumbers,
    std::vector<std::vector<int> >& atomIndices,
    std::vector<StructuralFormula>& structuralFormulas)
{
    int nSubstructures = *max_element(mol2Data.substructureIdx.begin(), mol2Data.substructureIdx.end());
    atomIndices.resize(nSubstructures);
    structuralFormulas.resize(nSubstructures);
    vector<int> structure2substructureIdx(atomicNumbers.size());
    
    for (int i = 0; i < mol2Data.atomId.size(); i++)
    {
        structure2substructureIdx[i] = atomIndices[mol2Data.substructureIdx[i]-1].size();
        atomIndices[mol2Data.substructureIdx[i]-1].push_back(i);
        structuralFormulas[mol2Data.substructureIdx[i]-1].z.push_back(atomicNumbers[i]);
    }

    for (int i = 0; i < structuralFormulas.size(); i++)
    {
        structuralFormulas[i].connectivity.clear();
        structuralFormulas[i].connectivity.resize(structuralFormulas[i].z.size());
    }

    for (const auto& bond : mol2Data.bonds)
    {
        int substructureIdx1 = mol2Data.substructureIdx[bond.first-1] - 1;
        int substructureIdx2 = mol2Data.substructureIdx[bond.second-1] - 1;

        if (substructureIdx1 != substructureIdx2)
            continue;
        int idx1 = structure2substructureIdx[bond.first-1];
        int idx2 = structure2substructureIdx[bond.second-1];
        structuralFormulas[substructureIdx1].connectivity[idx1].push_back(idx2);
        structuralFormulas[substructureIdx1].connectivity[idx2].push_back(idx1);
    }
    
}


void AssignAtomTypes::run()
{



    mAssigner.setAtomTypes(mAtomTypes);
    mAssigner.setDescriptorsSettings(mDescriptorsSettings);

    vector<string> resFiles;
    file_system_utilities::find_files("res", mChosenResFolder.string(), resFiles, false);

    ofstream out(mOutputFileName);
    vector<int> nOccurences(mAtomTypes.size(), 0);
    vector<int> nContainingMols(mAtomTypes.size(), 0);
    vector<std::set<int> > containingFormulaIndices(mAtomTypes.size());
    //vector<std::set<int> > 
    //[type idx]
    //vector<std::set<map<int, int> > > includingMolFormulas(mAtomTypes.size());
    vector<std::set<int> > inludingMolFormulaIdx(mAtomTypes.size());
    vector<StructuralFormula> structuralFormulas;

    int longestTypeId = 0;
    for (auto& type : mAtomTypes)
        if (type.id.size() > longestTypeId)
            longestTypeId = type.id.size();

    out << "n structures " << resFiles.size() << "\n";

    vector<StructuralFormula> uniqueStructuralFormulas;

    for (string& resFile : resFiles)
    {
        string fileNameCore = resFile.substr(0, resFile.size() - 4);
        string molFile = fileNameCore + ".mol2";
        filesystem::path molFilePath = mMolFolder / molFile;

        if (!filesystem::exists(molFilePath))
            on_error::throwException("missing file: '" + molFilePath.string() + "'", __FILE__, __LINE__);

        vector<bool> contains(mAtomTypes.size(), false);

        mol2_io::Mol2Data mol2Data;
        mol2_io::read(molFilePath.string(), mol2Data);

        StructureWithDescriptors structure;
        vector<int> atomicNumbers;
        mol2_io::atomicNumbers(mol2Data, atomicNumbers);

        structure.set(atomicNumbers, mol2Data.atomPosition, mol2Data.atomName);
        vector<int> typeId;
        vector<LocalCoordinateSystem<int> > lcs;
        mAssigner.assign(structure, typeId, lcs);

        vector<map<int, int> > formulas;
        findFormulas(mol2Data, atomicNumbers, formulas);

        out << "structure " << fileNameCore << "\n number of atoms " << typeId.size() << "\n";

        vector<vector<int> > atomIndices;
        vector<StructuralFormula> structuralFormulas;
        findSubstructures(mol2Data, atomicNumbers, atomIndices, structuralFormulas);

        out << " number of substructures " << atomIndices.size() << "\n";
        out << "   idx      formula      idx of uniqe str. formula\n";
        vector<int> formulaIdx(structuralFormulas.size());
        for (int i = 0; i < structuralFormulas.size(); i++)
        {
            int idx = -1;
            for (int j = 0; j < uniqueStructuralFormulas.size(); j++)
            {
                bool isomorphic =
                    structural_properties::molecularGraphIsomorphism(
                        structuralFormulas[i].z,
                        structuralFormulas[i].connectivity,
                        uniqueStructuralFormulas[j].z,
                        uniqueStructuralFormulas[j].connectivity);

                if (isomorphic)
                {
                    idx = j;
                    break;
                }
            }

            if (idx < 0)
            {
                idx = uniqueStructuralFormulas.size();
                uniqueStructuralFormulas.push_back(structuralFormulas[i]);
            }
            string formulaString = basic_chemistry_utilities::formulaAsString(formulas[i]);

            out << " " << setw(5) << i + 1 << "  " << setw(std::max<size_t>(13, formulaString.size())) << formulaString << "    " << idx << "\n";
            formulaIdx[i] = idx;
        }

        out << "  atom label    z    type id   substructure idx\n";
        int nAtoms = mol2Data.atomName.size();
        for (int i = 0; i < nAtoms; i++)
        {

            out << setw(12) << mol2Data.atomName[i] << "  " << setw(3) << atomicNumbers[i] << "  ";
            if (typeId[i] >= 0)
            {
                out << setw(longestTypeId + 2) << mAtomTypes[typeId[i]].id;
                nOccurences[typeId[i]]++;
                contains[typeId[i]] = true;
                //includingMolFormulas[typeId[i]].insert(formulas[mol2Data.substructureIdx[i] - 1]);
                inludingMolFormulaIdx[typeId[i]].insert(formulaIdx[mol2Data.substructureIdx[i] - 1]);

            }
            else
                out << setw(longestTypeId + 2) << "----";
            out << setw(10) << mol2Data.substructureIdx[i] << "\n";
        }



        for (int i = 0; i < contains.size(); i++)
            if (contains[i])
                nContainingMols[i]++;

    }
    out << "\n\n TYPE STATISTICS\n"
        << "  type id   n occurences    n containing   n with different\n"
        << "                             structures     formulas\n";


    vector<tuple<int, int, int, int> > typeStats;
    for (int i = 0; i < nOccurences.size(); i++)
    {
        out << setw(longestTypeId + 2) << mAtomTypes[i].id << "  " << setw(12) << nOccurences[i] << " "
            << setw(12) << nContainingMols[i] << setw(12) << inludingMolFormulaIdx[i].size() << "\n";
        typeStats.push_back({ nOccurences[i], nContainingMols[i], (int)inludingMolFormulaIdx[i].size(), i });
    }
    //

    out << "\n\n TYPE STATISTICS SORTED BY OCCURENCES\n"
        << "  type id   n occurences    n containing   n with different\n"
        << "                             structures     formulas\n";
    sort(typeStats.begin(), typeStats.end());
    reverse(typeStats.begin(), typeStats.end());
    for (int i = 0; i < nOccurences.size(); i++)
    {
        out << setw(longestTypeId + 2) << mAtomTypes[get<3>(typeStats[i])].id << "  " << setw(12) << get<0>(typeStats[i]) << " "
            << setw(12) << get<1>(typeStats[i]) << setw(12) << get<2>(typeStats[i]) << "\n";
    }

    // 
    typeStats.clear();
    out << "\n\n TYPE STATISTICS SORTED BY OCCURENCES including subtype contribution\n"
        << "  type id   n occurences    n containing   n with different\n"
        << "                             structures     formulas\n";
    vector<vector<int> > typesGeneralized;
    vector<vector<int> > hierarchyLevel;
    atom_typing_utilities::typeGeneralization(mAtomTypes, typesGeneralized, hierarchyLevel);

    for (int i = 0; i < nOccurences.size(); i++)
    {
        
        int nOcc = nOccurences[i];
        int nContaining = nContainingMols[i];
        int nFormulas = (int)inludingMolFormulaIdx[i].size();
        for(int generalizeIdx: typesGeneralized[i])
        {
            nOcc += nOccurences[generalizeIdx];
            nContaining += nContainingMols[generalizeIdx];
            nFormulas += (int)inludingMolFormulaIdx[generalizeIdx].size();
        }
        typeStats.push_back({ nOcc, nContaining, nFormulas, i });
    }

    sort(typeStats.begin(), typeStats.end());
    reverse(typeStats.begin(), typeStats.end());

    for (int i = 0; i < nOccurences.size(); i++)
    {
        out << setw(longestTypeId + 2) << mAtomTypes[get<3>(typeStats[i])].id << "  " << setw(12) << get<0>(typeStats[i]) << " "
            << setw(12) << get<1>(typeStats[i]) << setw(12) << get<2>(typeStats[i]) << "\n";
    }


    //

    auto it = min_element(nOccurences.begin(), nOccurences.end());

    out << "\n\nMin n occurences " << *it << " " << mAtomTypes[distance(nOccurences.begin(), it)].id << "\n";
    it = min_element(nContainingMols.begin(), nContainingMols.end());
    out << "Min n containing structures " << *it << " " << mAtomTypes[distance(nContainingMols.begin(), it)].id << "\n";

    out << "Min n containing structures with different structural formulas ";
    vector<pair<int, int> > nIndepFormulas_idx;
    for (int i = 0; i < inludingMolFormulaIdx.size(); i++)
        nIndepFormulas_idx.push_back({(int)inludingMolFormulaIdx[i].size(), i});
    auto it2 = min_element(nIndepFormulas_idx.begin(), nIndepFormulas_idx.end());
    out << it2->first << " " << mAtomTypes[it2->second].id << "\n";

    out.close();
}


