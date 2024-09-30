#include "discamb/QuantumChemistry/fragmentation.h"

#include "discamb/BasicChemistry/PeriodicTable.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/MathUtilities/GraphAlgorithms.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/BasicUtilities/utilities.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/Scattering/gar_utilities.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include <fstream>
#include <algorithm>

using namespace std;

namespace {


    struct AtomInSequence {
        int chain;
        int residue;
        int idxInResidue;
    };


}

namespace discamb{


    void FragmentConstructionData::set(
        const std::string& line)
    {
        vector<string> words;
        string_utilities::split(line, words);
        set(words);
    }

    /*
        e.g.
        {
            {"add": "C1 C2"},
            {
                "connect": "C3",
                "remove": "C4 C9",
                "remove and cap": "C5 C6"
            },
            {"connect": "C12 C17"}
        }

        AtomList include;
        AtomList connect;
        AtomList remove;
        AtomList removeAndReplaceWithHydrogen;


    */

    void FragmentConstructionData::clear()
    {
        include.clear();
        connect.clear();
        remove.clear();
        removeAndReplaceWithHydrogen.clear();
    }

    void FragmentConstructionData::set(
        const nlohmann::json& data)
    {
        clear();
        vector<string> words;

        AtomList atomList;

        if (data.find("add") != data.end())
        {    
            string_utilities::split(data.find("add")->get<string>(), words);
            crystal_structure_utilities::splitIntoAtomAndSymmOp(words, atomList, true);
            include.swap(atomList);
        }
        if (data.find("connect") != data.end())
        {
            string_utilities::split(data.find("connect")->get<string>(), words);
            crystal_structure_utilities::splitIntoAtomAndSymmOp(words, atomList, true);
            connect.swap(atomList);
        }
        if (data.find("remove") != data.end())
        {
            string_utilities::split(data.find("remove")->get<string>(), words);
            crystal_structure_utilities::splitIntoAtomAndSymmOp(words, atomList, true);
            remove.swap(atomList);
        }
        if (data.find("remove and cap") != data.end())
        {
            string_utilities::split(data.find("remove and cap")->get<string>(), words);
            crystal_structure_utilities::splitIntoAtomAndSymmOp(words, atomList, true);
            removeAndReplaceWithHydrogen.swap(atomList);
        }

    }

    void FragmentConstructionData::getAtomList(
        const UnitCellContent& ucContent,
        const std::vector<std::vector<UnitCellContent::AtomID> >& connectivity,
        std::vector<UnitCellContent::AtomID>& atomList,
        std::vector<std::pair<UnitCellContent::AtomID, UnitCellContent::AtomID> >& cappingHydrogens)
        const
    {
        atomList.clear();
        cappingHydrogens.clear();
        
        std::set <UnitCellContent::AtomID> uniqueSubsystemAtoms;
        vector<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > disconnectedBonds;
        vector<UnitCellContent::AtomID> fragment;
        vector<vector<UnitCellContent::AtomID> > modifiedConnectivity = connectivity;
        UnitCellContent::AtomID atomId;

        if (!connect.empty())
        {
            vector<UnitCellContent::AtomID> atomsToRemoveAndNotCap, atomsToRemoveAndCap, atomsToRemove;
            
            for (auto& atom : remove)
            {
                ucContent.findAtom(atom.first, atom.second, atomId);
                atomsToRemoveAndNotCap.push_back(atomId);
            }
            for (auto& atom : removeAndReplaceWithHydrogen)
            {
                ucContent.findAtom(atom.first, atom.second, atomId);
                atomsToRemoveAndCap.push_back(atomId);
            }
            atomsToRemove = atomsToRemoveAndNotCap;
            atomsToRemove.insert(atomsToRemove.end(), atomsToRemoveAndCap.begin(), atomsToRemoveAndCap.end());
            
            for (auto atom : connect)
            {
                ucContent.findAtom(atom.first, atom.second, atomId);
                structural_properties::findIncludingMolecule(atomId, atomsToRemove, connectivity,
                    fragment, disconnectedBonds, ucContent);
                uniqueSubsystemAtoms.insert(fragment.begin(), fragment.end());
            }
            

            for (auto bond : disconnectedBonds)
            {
                if (find(atomsToRemoveAndCap.begin(), atomsToRemoveAndCap.end(), bond.second) != atomsToRemoveAndCap.end())
                    cappingHydrogens.push_back(bond);
            }
        }
        
        for (auto& atom : include)
        {
            ucContent.findAtom(atom.first, atom.second, atomId);
            uniqueSubsystemAtoms.insert(atomId);
        }
        
        atomList.insert(atomList.end(), uniqueSubsystemAtoms.begin(), uniqueSubsystemAtoms.end());

    }


    void FragmentConstructionData::set(
        const std::vector<std::string>& words)
    {
        clear();
        *this = FragmentConstructionData();
        if (words.empty())
            return;
        std::set<string> validSectionNames { "$connect", "$remove", "$remove_and_cap" };

        if (words[0][0] == '$')
        {
            map<string, vector<string> > sections;
            string currentSection;
            for (auto const& word : words)
                if (word[0] == '$')
                {
                    currentSection = word;
                    if (validSectionNames.count(word) == 0)
                        on_error::throwException("invalid keyword in fragmentation definition: '" + word + "'", __FILE__, __LINE__);
                }
                else
                    sections[currentSection].push_back(word);
            for (auto& section : sections)
            {
                AtomList atomList;
                crystal_structure_utilities::splitIntoAtomAndSymmOp(section.second, atomList, true);

                if (section.first == "$connect")
                    connect = atomList;
                if (section.first == "$remove")
                    remove = atomList;
                if (section.first == "$remove_and_cap")
                    this->removeAndReplaceWithHydrogen = atomList;
            }
        }
        else
        {
            include.clear();
            if(words.size()==1)
                include.push_back({ words[0], "X,Y,Z" });
            if (words.size() == 2)
                include.push_back({ words[0], words[1] });
        }


    }


namespace fragmentation{

    void from_json(
        const nlohmann::json& data,
        const Crystal& crystal,
        vector < pair<string, string> >& atomList)
    {

    }

    void from_fragment_construction_data(
        const vector<vector< FragmentConstructionData> >& data,
        const Crystal& crystal,
        vector < pair<string, string> >& atomList)
    {

    }

    void from_fragment_data(
        const std::vector<FragmentData>& data,
        const Crystal& crystal,
        vector < pair<string, string> >& atomList)
    {
        int nFragments = data.size();
        vector<vector< FragmentConstructionData> > fragmentConstructionData(nFragments);

        for (int fragIdx = 0; fragIdx < nFragments; fragIdx++)
            fragmentConstructionData[fragIdx] = data[fragIdx].fragmentConstructionData;

        from_fragment_construction_data(fragmentConstructionData, crystal, atomList);
    }

    Vector3d capping_h_position(
        const Crystal& crystal, 
        const CappingHydrogen& cappingHydrogen)
    {
        SpaceGroupOperation bondedAtomSymmOp(cappingHydrogen.bondedAtomSymmOp);
        SpaceGroupOperation directingAtomSymmOp(cappingHydrogen.directingAtomSymmOp);
        return capping_h_position(crystal, cappingHydrogen.bondedAtom, bondedAtomSymmOp, cappingHydrogen.directingAtom, directingAtomSymmOp);
    }

    Vector3d capping_h_position(const Crystal& crystal,
        const std::string& bondedAtom, const SpaceGroupOperation& bondedAtomSymmOp,
        const std::string& directingAtom, const SpaceGroupOperation& directingAtomSymmOp)
    {
        Vector3d bondedAtomPosition = crystal_structure_utilities::atomPosition(bondedAtom, bondedAtomSymmOp, crystal);
        Vector3d directingAtomPosition = crystal_structure_utilities::atomPosition(directingAtom, directingAtomSymmOp, crystal);
        int directingAtomIdx = 0;
        for (int i = 0; i < crystal.atoms.size(); i++)
            if (crystal.atoms[i].label == directingAtom)
                directingAtomIdx = i;
        return structural_properties::capping_atom_position(bondedAtomPosition, directingAtomPosition, periodic_table::atomicNumber(crystal.atoms[directingAtomIdx].type));
    }

    
    void intermolecular(
        const Crystal& crystal,
        std::vector<std::vector<std::pair<std::string, std::string> > >& clusterAtoms,
        std::vector<std::string>& clustersLabels,
        std::vector<int>& clustersCharges,
        std::vector<int>& clustersSpinMultiplicity,
        int spinMultiplicityHint)
    {
        gar_utilities::findDeaultSubsystems(crystal, clusterAtoms, clustersLabels, clustersCharges, clustersSpinMultiplicity);
    }


    bool addedToFragmentsSet(
        vector<vector<int> > &fragments,
        const vector<int> &candidate,
        vector<vector<int> > &representedAtoms,
        int representedAtom,
        vector<int> &eliminatedSubgraphs,
        // defined only if some fragment is subgraph of candidate
        vector<int> &inclusionMap, // inclusionMap[i] index of fragment replacing i-th fragment, needed for handling representatives
        int &includingFragmentIdx)
    {
        eliminatedSubgraphs.clear();
        

        for (int fragIdx = 0; fragIdx < fragments.size(); fragIdx++)
            if (includes(fragments[fragIdx].begin(), fragments[fragIdx].end(), candidate.begin(), candidate.end()))
            {
                representedAtoms[fragIdx].push_back(representedAtom);
                includingFragmentIdx = fragIdx;
                return false;
            }
            else
                if (includes(candidate.begin(), candidate.end(),fragments[fragIdx].begin(), fragments[fragIdx].end()))
                    eliminatedSubgraphs.push_back(fragIdx);

        inclusionMap.resize(fragments.size());
        vector<vector<int> > newFragments, newRepresentedAtoms;
        
        int nRemovedSoFar=0;
        for (int fragIdx = 0; fragIdx < fragments.size(); fragIdx++)
        {
            // if is not eliminated
            if (find(eliminatedSubgraphs.begin(), eliminatedSubgraphs.end(), fragIdx) == eliminatedSubgraphs.end())
            {
                newFragments.push_back(fragments[fragIdx]);
                newRepresentedAtoms.push_back(representedAtoms[fragIdx]);
                inclusionMap[fragIdx] = fragIdx - nRemovedSoFar;
            }
            else
                nRemovedSoFar++;
        }

        newRepresentedAtoms.resize(newRepresentedAtoms.size() + 1);
        newRepresentedAtoms.back().push_back(representedAtom);


        for (int i : eliminatedSubgraphs)
        {
            inclusionMap[i] = newFragments.size();
            newRepresentedAtoms.back().insert(newRepresentedAtoms.back().end(), representedAtoms[i].begin(), representedAtoms[i].end());
        }
        newRepresentedAtoms.swap(representedAtoms);
        fragments.swap(newFragments);
        fragments.push_back(candidate);

        return true;
    }

    void findCappingHydrogens(
        const vector<vector<int> >& connectivity,
        const vector<int>& fragment,
        vector<pair<int, int> > & cappingHydrogens)
    {
        cappingHydrogens.clear();

        for (int atomIdx : fragment)
            for (int neighbour : connectivity[atomIdx])
                if (find(fragment.begin(), fragment.end(), neighbour) == fragment.end())
                    cappingHydrogens.push_back({ atomIdx, neighbour });
            
    }

    void findFragments(
        const vector<vector<int> >& connectivity, 
        const vector<vector<bool> >& isBondBreakable,
        const vector<int> &atomicNumbers,
        vector< vector < int> > &fragments,
        vector<vector<pair<int, int> > > &cappingHydrogens,
        vector<int>& representingFragment)
    {
        int range = 1;

        fragments.clear();
        representingFragment.clear();

        int atomIdx, nAtoms = connectivity.size();
        vector<vector<int> > subgraph;
        vector<int> fragment;
        vector<pair<int, int> > fragmentCappingHydrogens;
        int includingFragmentIdx;
        vector<int> eliminatedSubgraphs;
        vector<vector<int> > representedAtoms;
            


        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            if (atomicNumbers[atomIdx] != 1)
            {
                graph_algorithms::breadth_first_search(connectivity, isBondBreakable, atomIdx, subgraph, range);
                fragment.clear();
                for (auto v : subgraph)
                    fragment.insert(fragment.end(), v.begin(), v.end());
                sort(fragment.begin(), fragment.end());
                vector<int> inclusionMap;
                addedToFragmentsSet(fragments, fragment, representedAtoms, atomIdx, eliminatedSubgraphs, inclusionMap, includingFragmentIdx);
            }
            //if (addedToFragmentsSet(fragments, fragment, representedAtoms, atomIdx, eliminatedSubgraphs, inclusionMap, includingFragmentIdx))
            //{
            //    // others were subset of the fragment
            //    if (!eliminatedSubgraphs.empty())
            //        for (int& rf : representingFragment)
            //            rf = inclusionMap[rf];
            //    representingFragment.push_back(fragments.size());
            //}
            //else
            //    representingFragment.push_back(includingFragmentIdx);
        }

        representingFragment.resize(nAtoms);
        for (int fragIdx = 0; fragIdx < representedAtoms.size(); fragIdx++)
            for (int atIdx : representedAtoms[fragIdx])
                representingFragment[atIdx] = fragIdx;

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (atomicNumbers[atomIdx] == 1)
                representingFragment[atomIdx] = representingFragment[connectivity[atomIdx][0]];


        for (auto& fragment : fragments)
        {
            findCappingHydrogens(connectivity, fragment, fragmentCappingHydrogens);
            cappingHydrogens.push_back(fragmentCappingHydrogens);
        }



    }

}

}