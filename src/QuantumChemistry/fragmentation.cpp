#include "discamb/QuantumChemistry/fragmentation.h"

#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
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


    void FragmentPartConstructionData::set(
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

    void FragmentPartConstructionData::clear()
    {
        include.clear();
        connect.clear();
        remove.clear();
        removeAndReplaceWithHydrogen.clear();
    }

    void FragmentPartConstructionData::set(
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

    void FragmentPartConstructionData::getAtomList(
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


    void FragmentPartConstructionData::set(
        const std::vector<std::string>& words)
    {
        clear();
        *this = FragmentPartConstructionData();
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

    bool FragmentConstructionData::set(
        const nlohmann::json& data)
    {
        if (!data.is_object())
            return false;

        string errorMessage = "invalid format of JSON subsystems definition for Hirshfeld Atom Model calculations";

        *this = FragmentConstructionData();
        this->label = data.value("name", "");
        this->charge = data.value("charge", 0);
        this->spin_multiplicity = data.value("spin multiplicity", 1);

        if (data.find("atoms") == data.end())
            on_error::throwException(errorMessage, __FILE__, __LINE__);

        FragmentPartConstructionData fragmentPartConstructionData;

        if (data["atoms"].is_string())  //list of atoms
        {
            vector<string> words;
            string_utilities::split(data["atoms"].get<string>(), words);
            crystal_structure_utilities::splitIntoAtomAndSymmOp(words, fragmentPartConstructionData.include, true);
            this->fragmentPartConstructionData.push_back(fragmentPartConstructionData);
        }
        else
        {
            if (data["atoms"].is_object())
            {
                fragmentPartConstructionData.set(data["atoms"]);
                this->fragmentPartConstructionData.push_back(fragmentPartConstructionData);
            }
            else
            {
                if (data["atoms"].is_array())
                {
                    for (auto& group : data["atoms"])
                    {
                        fragmentPartConstructionData.set(group);
                        this->fragmentPartConstructionData.push_back(fragmentPartConstructionData);
                    }
                }
                else
                    on_error::throwException(errorMessage, __FILE__, __LINE__);
            }
        }


        return true;
    }

    void QmFragmentInCrystal::toXyzMol(
        const Crystal& crystal,
        std::vector<ChemicalElement>& elements,
        std::vector<Vector3d>& positions)
        const
    {
        positions.clear();
        elements.clear();

        map<string, int> label2idx;
        for (int atomIdx = 0; atomIdx < int(crystal.atoms.size()); atomIdx++)
            label2idx[crystal.atoms[atomIdx].label] = atomIdx;

        vector<int> atomicNumbersAsymmUnit;
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbersAsymmUnit);

        for (auto const& atom : atoms.atomList)
        {
            if (label2idx.find(atom.first) != label2idx.end())
            {
                positions.push_back(crystal_structure_utilities::atomPosition(atom.first, atom.second, crystal));
                elements.push_back(ChemicalElement(int(atomicNumbersAsymmUnit[label2idx[atom.first]])));
            }
            else
                on_error::throwException(string("invalid atom label '") + atom.first + string("'"), __FILE__, __LINE__);
        }
        for (auto const& capH : atoms.cappingHydrogens)
        {
            positions.push_back(fragmentation::capping_h_position(crystal, capH.bondedAtom, capH.bondedAtomSymmOp, capH.directingAtom, capH.directingAtomSymmOp));
            elements.push_back(ChemicalElement(1));
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
        const vector<vector< FragmentPartConstructionData> >& data,
        const Crystal& crystal,
        vector < pair<string, string> >& atomList)
    {

    }

    void make_qm_fragments(
        const Crystal& crystal,
        const std::vector<FragmentConstructionData>& fragmentConstructionData,
        std::vector<QmFragmentInCrystal>& qmFragments)
    {
        qmFragments.clear();
        
        // fragmentPartConstructionData to atom list

        UnitCellContent unitCellContent;
        unitCellContent.set(crystal);
        set<UnitCellContent::AtomID> fragmentAtoms;
        set<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > cappingHydrogens;
        vector<UnitCellContent::AtomID> atomList;
        vector<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > cappingHydrogenList;
        vector<vector<UnitCellContent::AtomID> > connectivity;

        structural_properties::calcUnitCellConnectivity(unitCellContent, connectivity, 0.4);
        string label, symmetryOperationStr;

        for (auto& fragmentData : fragmentConstructionData)
        {
            fragmentAtoms.clear();
            atomList.clear();
            QmFragmentInCrystal subsystem;
            subsystem.charge = fragmentData.charge;
            subsystem.label = fragmentData.label;
            subsystem.spin_multiplicity = fragmentData.spin_multiplicity;

            for (auto& fragmentPartConstructionData : fragmentData.fragmentPartConstructionData)
            {
                fragmentPartConstructionData.getAtomList(unitCellContent, connectivity, atomList, cappingHydrogenList);
                fragmentAtoms.insert(atomList.begin(), atomList.end());
                cappingHydrogens.insert(cappingHydrogenList.begin(), cappingHydrogenList.end());
            }
            for (auto& atom : fragmentAtoms)
            {
                unitCellContent.interpreteAtomID(atom, label, symmetryOperationStr);
                subsystem.atoms.atomList.push_back({ label, symmetryOperationStr });
            }
            for (auto& cappingH : cappingHydrogens)
            {
                subsystem.atoms.cappingHydrogens.resize(subsystem.atoms.cappingHydrogens.size() + 1);
                auto& capH = subsystem.atoms.cappingHydrogens.back();
                unitCellContent.interpreteAtomID(cappingH.first, label, symmetryOperationStr);
                capH.bondedAtom = label;
                capH.bondedAtomSymmOp = symmetryOperationStr;
                unitCellContent.interpreteAtomID(cappingH.second, label, symmetryOperationStr);
                capH.directingAtom = label;
                capH.directingAtomSymmOp = symmetryOperationStr;
            }
            qmFragments.push_back(subsystem);
        }

    }


    //void from_fragment_data(
    //    const std::vector<FragmentData>& data,
    //    const Crystal& crystal,
    //    vector < pair<string, string> >& atomList)
    //{
    //    int nFragments = data.size();
    //    vector<vector< FragmentPartConstructionData> > fragmentPartConstructionData(nFragments);

    //    for (int fragIdx = 0; fragIdx < nFragments; fragIdx++)
    //        fragmentPartConstructionData[fragIdx] = data[fragIdx].fragmentPartConstructionData;

    //    from_fragment_construction_data(fragmentPartConstructionData, crystal, atomList);
    //}

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

//    void find_fragments_atoms(
//        const Crystal& crystal,
//        std::vector<FragmentConstructionData>& fragmentConstructionData,
//        std::vector<FragmentAtoms>& fragmentsAtoms)
//    {
//        int fragmentIdx, nFragments = fragmentConstructionData.size();
//        fragmentsAtoms.clear();
//        fragmentsAtoms.resize(nFragments);
//        UnitCellContent unitCellContent(crystal);
//        vector<vector<UnitCellContent::AtomID> > connectivity;
//        structural_properties::calcUnitCellConnectivity(unitCellContent, connectivity, 0.4);
//        
//        for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
//        {
//            vector<UnitCellContent::AtomID> atomList;
//            vector<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > cappingHydrogens;
//            fragmentConstructionData[fragmentIdx].getAtomList(
//                unitCellContent,
//                connectivity,
//                atomList,
//                cappingHydrogens);
//            UnitCellContent::AtomID id;
//            string label, symmOpStr;
//            int nAtoms = atomList.size();
//            fragmentsAtoms[fragmentIdx].atomList.resize(nAtoms);
//
//            int nCappingHydrogens = cappingHydrogens.size();
//            fragmentsAtoms[fragmentIdx].cappingHydrogens.resize(nCappingHydrogens);
//            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
//            {
//                unitCellContent.interpreteAtomID(atomList[atomIdx], label, symmOpStr);
//                fragmentsAtoms[fragmentIdx].atomList[atomIdx].first = label;
//                fragmentsAtoms[fragmentIdx].atomList[atomIdx].second = symmOpStr;
//            }
//            for (int capHydrogenIdx = 0; capHydrogenIdx < nCappingHydrogens; capHydrogenIdx++)
//            {
//                unitCellContent.interpreteAtomID(cappingHydrogens[capHydrogenIdx].first, label, symmOpStr);
//                fragmentsAtoms[fragmentIdx].cappingHydrogens[capHydrogenIdx].bondedAtom = label;
//                fragmentsAtoms[fragmentIdx].cappingHydrogens[capHydrogenIdx].bondedAtomSymmOp = symmOpStr;
//                unitCellContent.interpreteAtomID(cappingHydrogens[capHydrogenIdx].second, label, symmOpStr);
//                fragmentsAtoms[fragmentIdx].cappingHydrogens[capHydrogenIdx].directingAtom = label;
//                fragmentsAtoms[fragmentIdx].cappingHydrogens[capHydrogenIdx].directingAtomSymmOp = symmOpStr;
//            }
//        }
//
//    }
//
//    void find_fragments_atoms(
//        const Crystal& crystal,
//        std::vector<FragmentData>& fragmentData,
//        std::vector<FragmentAtoms>& fragmentsAtoms)
//    {
//        vector<FragmentConstructionData> fragmentConstructionData;
//        for (auto& const data : fragmentData)
//            fragmentConstructionData.push_back(data.fragmentConstructionData);
//    }

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
