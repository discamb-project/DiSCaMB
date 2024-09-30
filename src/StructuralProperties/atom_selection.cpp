#include "discamb/StructuralProperties/atom_selection.h"
#include "discamb/BasicUtilities/OnError.h"

#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/BasicUtilities/StringUtilities.h"

using namespace std;

namespace discamb{

    
    void AtomSubsetSelectionData::set(
        nlohmann::json const& data)
    {
        *this = AtomSubsetSelectionData();
        vector<string> words, words2;
        string s;

        

        if (data.is_string())
        {            
            crystal_structure_utilities::stringToAtomList(data.get<string>(), atomList, ';');

        }
        else
        {
            if (data.is_object())
            {
                if (data.find("atoms") == data.end())
                    on_error::throwException(string("invalid entry in aspher.json, should contain entry atoms, entry:\n") +
                        data.dump(), __FILE__, __LINE__);
                crystal_structure_utilities::stringToAtomList(data.find("atoms").value().get<string>(), atomList, ';');

                // "periodic directions": "(1,0,1),(0,1,0)"
                if (data.find("periodic directions") != data.end())
                {
                    string periodicDirectionsString = data.find("periodic directions").value().get<string>();
                    string_utilities::split(periodicDirectionsString, words, ')');

                    for (auto& word : words)
                    {
                        replace(word.begin(), word.end(), '(', ',');
                        string_utilities::split(word, words2, ',');
                        if (words2.size() != 3)
                            on_error::throwException(string("invalid string in \"periodic directions\" section in aspher.json: \"") +
                                periodicDirectionsString + string("\""), __FILE__, __LINE__);
                        periodicDirections.push_back(
                            Vector3i(stoi(words2[0]), stoi(words2[1]), stoi(words2[2])));
                    }
                }

                if (data.find("symmetry equivalent") != data.end())
                    allSymmetryEquivalent = data.find("symmetry equivalent").value().get<bool>();
            }
            else
                on_error::throwException(string("invalid entry in aspher.json:\n") + data.dump(), __FILE__, __LINE__);

        }


    }

namespace atom_selection {

    void select_subset(
        const UnitCellContent& ucContent,
        const vector< pair<string, string> >& set,
        const AtomSubsetSelectionData& subsetSelectionData,
        vector< pair<string, string> >& subset)
    {
        vector<UnitCellContent::AtomID> atomSet, atomSubset;
        UnitCellContent::AtomID atomId;
        for (auto const & atom : set)
        {
            ucContent.findAtom(atom.first, atom.second, atomId);
            atomSet.push_back(atomId);
        }
        select_subset(ucContent, atomSet, subsetSelectionData, atomSubset);
        pair<string, string> atomLabelAndSymmOp;
        for (auto& atom : atomSubset)
        {
            ucContent.interpreteAtomID(atom, atomLabelAndSymmOp.first, atomLabelAndSymmOp.second);
            subset.push_back(atomLabelAndSymmOp);
        }
    }

    void select_subset(
        const UnitCellContent& ucContent,
        const std::vector<UnitCellContent::AtomID>& atomSet,
        const AtomSubsetSelectionData& subsetSelectionData,
        std::vector<UnitCellContent::AtomID>& atomSubset)
    {
        int atomIdx, nAtoms = int(atomSet.size());
        vector<bool> selected(nAtoms, false);

        if (subsetSelectionData.allSymmetryEquivalent && !subsetSelectionData.periodicDirections.empty())
            on_error::throwException(
                "invalid specification of multipole cluster,"
                "defining \"symmetry equivalent\": true and "
                "\"periodic directions\" not allowed", __FILE__, __LINE__);

        if (subsetSelectionData.allSymmetryEquivalent)
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                string label, symmetryOperationString;
                ucContent.interpreteAtomID(atomSet[atomIdx], label, symmetryOperationString);
                for (auto const& atom : subsetSelectionData.atomList)
                    if (atom.first == label)
                        selected[atomIdx] = true;
            }
        else
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                UnitCellContent::AtomID selectorAtomId;

                for (auto const& atom : subsetSelectionData.atomList)
                {
                    ucContent.findAtom(atom.first, atom.second, selectorAtomId);

                    if (selectorAtomId.atomIndex == atomSet[atomIdx].atomIndex)
                        if (crystal_structure_utilities::isSuperLatticeNode(atomSet[atomIdx].unitCellPosition, selectorAtomId.unitCellPosition, subsetSelectionData.periodicDirections))
                            selected[atomIdx] = true;
                }
            }

        for (int idx = 0; idx < selected.size(); idx++)
            if (selected[idx])
                atomSubset.push_back(atomSet[idx]);

    }

    void remove_subset(
        const UnitCellContent& ucContent,
        std::vector<std::pair<std::string, std::string> >& set,
        const AtomSubsetSelectionData& subsetSelectionData)
    {

        vector<UnitCellContent::AtomID> atomSet;//, atomSubset;
        UnitCellContent::AtomID atomId;
        for (auto const& atom : set)
        {
            ucContent.findAtom(atom.first, atom.second, atomId);
            atomSet.push_back(atomId);
        }

        remove_subset(ucContent, atomSet, subsetSelectionData);

        pair<string, string> atomLabelAndSymmOp;
        vector< pair<string, string> > newSet;
        for (auto& atom : atomSet)
        {
            ucContent.interpreteAtomID(atom, atomLabelAndSymmOp.first, atomLabelAndSymmOp.second);
            newSet.push_back(atomLabelAndSymmOp);
        }
        set.swap(newSet);
    }

    void remove_subset(
        const UnitCellContent& ucContent,
        std::vector<UnitCellContent::AtomID>& set,
        const AtomSubsetSelectionData& subsetSelectionData)
    {
        vector<UnitCellContent::AtomID> atomsToRemove, newSet;
        select_subset(ucContent, set, subsetSelectionData, atomsToRemove);
        for (auto& atom : set)
            if (find(atomsToRemove.begin(), atomsToRemove.end(), atom) == atomsToRemove.end())
                newSet.push_back(atom);
        set.swap(newSet);
    }


    void merge(
        const UnitCellContent& ucContent,
        const std::vector<std::pair<std::string, std::string> >& set1,
        const std::vector<std::pair<std::string, std::string> >& set2,
        std::vector<std::pair<std::string, std::string> >& result)
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }

    void merge(
        const std::vector<UnitCellContent::AtomID>& set1,
        const std::vector<UnitCellContent::AtomID>& set2,
        const std::vector<UnitCellContent::AtomID>& result)
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }



} //namespace atom_selection

} //namespace discamb

