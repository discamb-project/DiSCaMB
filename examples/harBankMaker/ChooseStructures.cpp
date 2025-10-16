#include "ChooseStructures.h"

#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/IO/mol2_io.h"

#include "json.hpp"


#include <fstream>

using namespace discamb;
using namespace std;

namespace {

    int formulaDotProduct(
        const std::map<int, int>& formula1,
        const std::map<int, int>& formula2)
    {
        set<int> types;
        for (auto item : formula1)
            types.insert(item.first);
        for (auto item : formula2)
            types.insert(item.first);

        int result = 0;

        for (int type : types)
        {
            int n1 = 0;
            int n2 = 0;
            if (formula1.find(type) != formula1.end())
                n1 = formula1.find(type)->second;

            if (formula2.find(type) != formula2.end())
                n2 = formula2.find(type)->second;

            result += n1 * n2;
        }

        return result;
    }


}


ChooseStructures::ChooseStructures()
{

}

ChooseStructures::~ChooseStructures()
{

}

void ChooseStructures::set()
{
    nlohmann::json data;
    readSettings(data);

    string folder;
    if (mMolFolder.empty())
        folder = (filesystem::current_path()/string("mol")).string();
    else
        folder = mMolFolder.string();

    mMolFolder = data.value("mol folder", folder);

    mAllowLessThanTarget = data.value("allow less molecules than target", mAllowLessThanTarget);

    
    mAtomTypesLog = data.value("type assignment log file", mAtomTypesLog);
    mTargetNumberOfContainingMolecules = data.value("target n containing molecules", mTargetNumberOfContainingMolecules);

    mStructuresToExclude.clear();
    string s = data.value("exclude structures", string());
    string_utilities::split(s, mStructuresToExclude,',');



    //mStructuresToExclude
}

void ChooseStructures::readAtomTypeAssignement()
{
    ifstream in(mAtomTypesLog);
    
    if (!in.good())
        on_error::throwException("cannot read atom types assignment log file: '" + mAtomTypesLog + "'", __FILE__, __LINE__);

    string s, line, atomLabel, typeId;
    vector<string> words;
    int nStructures, nAtoms, structureIdx, atomIdx, z, substructureIdx, nSubstructures;

    in >> s >> s >> nStructures;

    
    mParentStructureNames.resize(nStructures);
    mSubstructures.clear();


    for (structureIdx = 0; structureIdx < nStructures; structureIdx++)
    {
        /*
structure ABEQUW.search1
 number of atoms 32
 number of substructures 1
   idx      formula      idx of uniqe str. formula
     1     H14C13N2O3    6
  atom label    z    type id   substructure idx
          O1    8     O201         1
        */
        int nSubstructuresBefore = mSubstructures.size();
        in  >> s >> mParentStructureNames[structureIdx]
            >> s >> s >> s >> nAtoms
            >> s >> s >> s >> nSubstructures;
        getline(in, line);
        getline(in, line);
        vector<string> substructureFormulaString(nSubstructures);
        vector<int> substructureStructuralFormulaIdx(nSubstructures);
        for (substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
            in >> s >> substructureFormulaString[substructureIdx] >> substructureStructuralFormulaIdx[substructureIdx];

        getline(in, line);
        getline(in, line);

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            
            in >> atomLabel >> z >> typeId >> substructureIdx;
            int totalSubstructureIdx = substructureIdx + nSubstructuresBefore-1;
            if (totalSubstructureIdx >= mSubstructures.size())
            {
                mSubstructures.resize(totalSubstructureIdx + 1);
                mSubstructures[totalSubstructureIdx].idxInMol2 = substructureIdx;
                mSubstructures[totalSubstructureIdx].parentStructure = structureIdx;
            }
            mSubstructures[totalSubstructureIdx].atomLabels.push_back(atomLabel);
            mSubstructures[totalSubstructureIdx].z.push_back(z);
            mSubstructures[totalSubstructureIdx].types.push_back(typeId);
            mSubstructures[totalSubstructureIdx].formulaIdx = substructureStructuralFormulaIdx[substructureIdx - 1];
            mSubstructures[totalSubstructureIdx].formulaString = substructureFormulaString[substructureIdx - 1];
        }
    }

    for (auto& structure : mSubstructures)
        basic_chemistry_utilities::getFormula(structure.z, structure.formula);


    getline(in, line);
    getline(in, line);
    getline(in, line);

    /*
 TYPE STATISTICS
  type id   n occurences    n containing structures
   O001             0            0

    */
    getline(in, line);
    getline(in, line);
    getline(in, line);
    mTypeStatistics.clear();
    vector<string> underrepresentedTypes;
    while (in.good())
    {
        getline(in, line);
        string_utilities::split(line, words);
        if (words.size() == 4)
        {
            mTypeStatistics.resize(mTypeStatistics.size() + 1);
            mTypeStatistics.back().nContainingStructures = stoi(words[2]);
            mTypeStatistics.back().nOccurences = stoi(words[1]);
            mTypeStatistics.back().typeLabel = words[0];
            int nDifferentMol = stoi(words[3]);
            mTypeStatistics.back().nOccurencesInMoleculesWithDiffrentFormula = nDifferentMol;
            if (nDifferentMol < mTargetNumberOfContainingMolecules)
                underrepresentedTypes.push_back(words[0]);
        }
    }

    in.close();

    if(!underrepresentedTypes.empty() && !mAllowLessThanTarget)
    {
        cout << "ERROR: some atom types occures in less chemical units with different structural formulas than expected:\n";
        for (auto s : underrepresentedTypes)
            cout << " " << s << "\n";
        exit(0);
    }
     

}

bool ChooseStructures::search_further(
    const vector<int> & nContainingMolecules)
    const
{
    for (int i = 0; i < nContainingMolecules.size(); i++)
        if (nContainingMolecules[i] < mTargetNumberOfContainingMolecules && nContainingMolecules[i] < mTypeStatistics[i].nOccurencesInMoleculesWithDiffrentFormula)
            return true;
    return false;

}

double ChooseStructures::formulaSimilarity(
    const std::map<int, int>& formula1,
    const std::map<int, int>& formula2)
{
    double f1 = formulaDotProduct(formula1, formula2) / sqrt(double(formulaDotProduct(formula1, formula1) * formulaDotProduct(formula2, formula2)));
    int n1 = 0;
    int n2 = 0;
    for (auto item : formula1)
        n1 += item.second;
    for (auto item : formula2)
        n2 += item.second;
    double f2 = abs(n1 - n2) / double(n1 + n2);
    return f1 * f2;
}

double ChooseStructures::calcScore(
    int structureIdx,
    const std::vector<std::vector<int> >& containingMolecules,
    const std::set<int>& typesInStructure,
    const vector<bool> &typeInstancesFound,
    const vector<bool> &formulaUsed)
    const
{
    double score = 0;

    

    for (int typeIdx: typesInStructure)
        //if(containingMolecules[typeIdx].size() < mTargetNumberOfContainingMolecules && containingMolecules[typeIdx].size() < mTypeStatistics[typeIdx].nOccurencesInMoleculesWithDiffrentFormula)
        if(!typeInstancesFound[typeIdx])
        {
            double averageMolSimilarity = 0;
            if (containingMolecules[typeIdx].size() == 0)
                averageMolSimilarity = 1.0;
            else
            {
                for (int molIdx : containingMolecules[typeIdx])
                    averageMolSimilarity += formulaSimilarity(mSubstructures[structureIdx].formula, mSubstructures[structureIdx].formula);
                averageMolSimilarity /= containingMolecules[typeIdx].size();
            }
            score += (mTargetNumberOfContainingMolecules - containingMolecules[typeIdx].size()) * (averageMolSimilarity + 0.1);
        }
    score /= mSubstructures[structureIdx].atomLabels.size();
    int formulaIdx = mSubstructures[structureIdx].formulaIdx;
    if (formulaUsed[formulaIdx])
        score /= 100.0;
    return score;
}

void ChooseStructures::readChosenMol(
    std::vector<std::string>& structureName,
    std::vector<int>& mol_idx)
{
    structureName.clear();
    mol_idx.clear();

    if (!filesystem::exists(filesystem::path("chosen_mol")))
        on_error::throwException("missing file 'chosen_mol'", __FILE__, __LINE__);

    ifstream input("chosen_mol");
    vector<string> words;
    string line;
    do
    {
        getline(input, line);
        string_utilities::split(line, words);
        if (words.size() > 2)
        {
            structureName.push_back(words[0]);
            mol_idx.push_back(stoi(words[1]));
        }
    } while (!words.empty());
    input.close();

}


void ChooseStructures::run()
{
    vector<int> chosenSubstructureIdx;

    readAtomTypeAssignement();

    int substructureIdx, nSubstructures = mSubstructures.size();
    int nTypes = mTypeStatistics.size();
    vector<double> score(nSubstructures);
    vector<int> nContainingMolecules(mTypeStatistics.size(), 0);
    vector< vector<int> > containingMolecules(mTypeStatistics.size());

    // target number of instances of given type already found
    vector<bool> typeInstancesFound(nTypes, false);

    int nStructuralFormulas = 0;
    for (substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
        if (mSubstructures[substructureIdx].formulaIdx + 1 > nStructuralFormulas)
            nStructuralFormulas = mSubstructures[substructureIdx].formulaIdx + 1;
    // formulaUsed[i] i-th structural formula already used
    vector<bool> formulaUsed(nStructuralFormulas, false);


    //

    vector<bool> searchable(nSubstructures, true);

    // exclude structures specified in settings file to be excluded
    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        bool exclude = false;
        string structureName = mParentStructureNames[mSubstructures[substructureIdx].parentStructure];

        for (string& structureToExclude : mStructuresToExclude)           
            if (structureName.find(structureToExclude) != string::npos)
                exclude = true;
        if (exclude)
            searchable[substructureIdx] = false;
    }
    

    // substructureTypes[i] list of unique types (as indices) in i-th structure
    vector<std::set<int> > substructureTypes(nSubstructures);
    map<string, int> typeName2Idx;
    for (int i = 0; i < mTypeStatistics.size(); i++)
        typeName2Idx[mTypeStatistics[i].typeLabel] = i;

    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
        for (string& typeLabel : mSubstructures[substructureIdx].types)
            if (typeLabel != "----")
                substructureTypes[substructureIdx].insert(typeName2Idx[typeLabel]);

    //

    //while (search_further(nContainingMolecules))
    while(find(typeInstancesFound.begin(), typeInstancesFound.end(), false)!= typeInstancesFound.end())
    {
        if (find(searchable.begin(), searchable.end(), true) == searchable.end())
            on_error::throwException("bug", __FILE__, __LINE__);

        for (substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
        {
            if (searchable[substructureIdx])
            {
                // check if has only types for which multiplicity has been satisfied
                bool areThereAtomsToUse = false;
                for (int t : substructureTypes[substructureIdx])
                    //if (nContainingMolecules[t] < mTargetNumberOfContainingMolecules &&
                    //    nContainingMolecules[t] < mTypeStatistics[t].nOccurencesInMoleculesWithDiffrentFormula)
                    if(!typeInstancesFound[t])
                        areThereAtomsToUse = true;
                if(!areThereAtomsToUse)
                {
                    searchable[substructureIdx] = false;
                    score[substructureIdx] = -1;
                }

                if (searchable[substructureIdx])
                {
                    score[substructureIdx] = calcScore(substructureIdx, containingMolecules, substructureTypes[substructureIdx], typeInstancesFound, formulaUsed);
                    if (score[substructureIdx] == 0.0)
                    {
                        searchable[substructureIdx] = false;
                        score[substructureIdx] = -1;
                    }
                }
            }
        }

        auto it = std::max_element(score.begin(), score.end());
        int idx = distance(score.begin(), it);
        cout << mParentStructureNames[mSubstructures[idx].parentStructure] + " " + to_string(mSubstructures[idx].idxInMol2) << "\n";
        searchable[idx] = false;
        score[idx] = -1;
        chosenSubstructureIdx.push_back(idx);
        formulaUsed[mSubstructures[idx].formulaIdx] = true;
        for (int typeIdx : substructureTypes[idx])
        {
            nContainingMolecules[typeIdx]++;
            containingMolecules[typeIdx].push_back(idx);
            if (nContainingMolecules[typeIdx] == mTargetNumberOfContainingMolecules)
                typeInstancesFound[typeIdx] = true;
        }
    }


    ofstream out("chosen_mol");
    for (int idx : chosenSubstructureIdx)
    {
        string formula;
        for (auto& item : mSubstructures[idx].formula)
        {
            formula += " " + periodic_table::symbol(item.first);
            if (item.second > 1)
                formula += to_string(item.second);
        }
        out << mParentStructureNames[mSubstructures[idx].parentStructure] << "  "
            << mSubstructures[idx].idxInMol2 << " "
            << formula << "\n";
    }

    out << "\nType containing molecules\n\n";

    for (int i=0;i<containingMolecules.size();i++)
    {
        out << "\n" << mTypeStatistics[i].typeLabel << " n mol " << containingMolecules[i].size() << "\n\n";
        for (int substructureIdx : containingMolecules[i])
            out << mParentStructureNames[mSubstructures[substructureIdx].parentStructure] << " " << mSubstructures[substructureIdx].idxInMol2 << "\n";
    }

    out.close();

    out.open("copy_chosen_mol.bat");
    out << "mkdir chosen_mol_dir\n";
    for (int i : chosenSubstructureIdx)
    {
        string structureName = mParentStructureNames[mSubstructures[i].parentStructure];
        out << "copy " << (mMolFolder / (structureName + ".mol2")).string() << " chosen_mol_dir\n";
    }
    out.close();
}


