#include "discamb/Scattering/disordered_structure_fragments.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include<fstream>
#include<set>

using namespace std;

namespace discamb{
namespace disordered_structure_fragments{

    namespace {

        string convert_internal_altloc_label(
            const string& label)
        {
            string newLabel;
            size_t pos = label.find('.');
            if (pos == string::npos)
                newLabel = label;
            else
            {
                string dotChar = label.substr(pos, 2);
                newLabel = label.substr(0, pos) + label.substr(pos + 2) + dotChar;
            }
            return newLabel;
        }

        void correct_zero_weights_for_ordered(
            int nAtoms,
            const vector<int> &orderedAtoms,
            vector<vector<pair<int, double> > > & ordered_parts)
        {
            // if there is zero weights for some ordered atom set it to 1.0
            
            vector<bool> isOrdered(nAtoms, false);
            vector<double> totalWeight(nAtoms, 0.0);
            for (int i = 0; i < orderedAtoms.size(); i++)
                isOrdered[orderedAtoms[i]] = true;

            for (auto& part : ordered_parts)
                for (auto& atom : part)
                    totalWeight[atom.first] += atom.second;

            for (auto& part : ordered_parts)
                for (auto& atom : part)
                    if (isOrdered[atom.first] && totalWeight[atom.first] == 0.0)
                        atom.second = 1.0;
        }
    }

    void split_structure_with_altlocs(
        const Crystal& crystal,
        std::vector < char >& altLocs,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts)
    {

    }


//    void split_with_labels_0(
//        const Crystal& crystal,
//        std::vector < std::vector<std::pair<std::string, double> > >& _ordered_parts)
//    {
//        _ordered_parts.clear();
//        vector < vector<pair<int, double> > > ordered_parts;
//
//        //map<string, int> label2idx;
//        //for (int i = 0; i < crystal.atoms.size(); i++)
//        //    label2idx[crystal.atoms[i].label] = i;
//
//
//        map<char, int> disorderGroupIdxMap;
//        for (int i = 65; i <= 90; i++)
//            disorderGroupIdxMap[char(i)] = i-65;
//
//        set<string> unique_disorder_part_labels;
//        vector<string> words;
//        vector<int> atom_disorder_part;
//        //aletrnativeAtomGroups[i].first - disorder group 
//        //aletrnativeAtomGroups[i].second - atom index
//        map<string, vector<pair<int, int> > > aletrnativeAtomGroups;
//        vector<int> orderedAtoms;
//
//        for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
//        {
//            auto const& atom = crystal.atoms[atomIdx];
//
//            string_utilities::split(atom.label, words, '.');
//            int disorder_part = -1;
//            if (words.size() == 2)
//            {
//                if (words[1].size() > 1)
//                {
//                    string message = "only single character notation is supported for alternative "
//                        "conformation coding in atom labels, the label " + atom.label +
//                        " cannot be used for automatic fragments generation";
//                    on_error::throwException(message, __FILE__, __LINE__);
//                }
//                else
//                {
//                    char c = toupper(words[1][0]);
//                    if (disorderGroupIdxMap.count(c)==0)
//                    {
//                        string message = "only A-Z,a-z characters are supported for alternative "
//                            "conformation coding in atom labels, the label " + atom.label +
//                            " cannot be used for automatic fragments generation";
//                        on_error::throwException(message, __FILE__, __LINE__);
//                    }
//                    disorder_part = disorderGroupIdxMap[c];
//                }
//                aletrnativeAtomGroups[words[0]].push_back({ disorder_part , atomIdx });
//            }
//            else
//                orderedAtoms.push_back(atomIdx);
//
//            atom_disorder_part.push_back(disorder_part);
//        }
//
//        for (auto& group : aletrnativeAtomGroups)
//        {
//            sort(group.second.begin(), group.second.end());
//
//            //validate
////            for (int i = 0; i < group.second.size(); i++)
////                if (group.second[i].first != i)
////                {
////                    string message = "inconsistent labeling of alternative configurations for atoms " + group.first +
////                        " e.g. A, B, D (C is missing)";
////                    on_error::throwException(message, __FILE__, __LINE__);
////                }
//        }
//
//        int max_disorder_part = *max_element(atom_disorder_part.begin(), atom_disorder_part.end());
//        if (max_disorder_part <= 0)
//            return;
//
//        int nParts = max_disorder_part + 1;
//
//        ordered_parts.resize(nParts);
//        for (int partIdx = 0; partIdx < nParts; partIdx++)
//        {
//
//            Crystal c = crystal;
//            c.atoms.clear();
//
//            // add ordered part
//            for (int atomIdx : orderedAtoms)
//            {
//                c.atoms.push_back(crystal.atoms[atomIdx]);
//                ordered_parts[partIdx].push_back({ atomIdx, 0.0 });
//            }
//            // add disordered part
//            for (auto& group : aletrnativeAtomGroups)
//            {
//                int group_size = group.second.size();
//                int idx_in_group = partIdx;
//                if (group_size < partIdx + 1)
//                    idx_in_group = group_size - 1;
//                int atomIdx = group.second[idx_in_group].second;
//                c.atoms.push_back(crystal.atoms[atomIdx]);
//                ordered_parts[partIdx].push_back({ atomIdx, 1.0 });
//            }
//            // find connectivity and use in weight calculations
//            vector<vector<pair<int, string> > > connectivity;
//            structural_properties::asymmetricUnitConnectivity(c, connectivity, 0.4);
//            int nOrdered = orderedAtoms.size();
//            for (int i = 0; i < nOrdered; i++)
//            {
//                double weight = 1.0;
//                for (int j = 0; j < connectivity[i].size(); j++)
//                    if (connectivity[i][j].first >= nOrdered)
//                    {
//                        int neighbourIdxInOriginalCrystal = ordered_parts[partIdx][connectivity[i][j].first].first;
//                        weight *= crystal.atoms[neighbourIdxInOriginalCrystal].occupancy;
//                    }
//                ordered_parts[partIdx][i].second = weight;
//            }
//
//        }
//        
//        correct_zero_weights_for_ordered(crystal.atoms.size(), orderedAtoms, ordered_parts);
//
//        int nOrderedParts = ordered_parts.size();
//        _ordered_parts.resize(nOrderedParts);
//        for (int partIdx = 0; partIdx < nOrderedParts; partIdx++)
//        {
//            _ordered_parts[partIdx].resize(ordered_parts[partIdx].size());
//            for (int atomIdx = 0; atomIdx < ordered_parts[partIdx].size(); atomIdx++)
//            {
//                int idxInCrystal = ordered_parts[partIdx][atomIdx].first;
//                double weight = ordered_parts[partIdx][atomIdx].second;
//                _ordered_parts[partIdx][atomIdx] = { crystal.atoms[idxInCrystal].label, weight };
//            }
//        }
//        
//    }
    
    void split_with_labels_internal_altloc(
        const Crystal& _crystal,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts)
    {
        Crystal crystal = _crystal;
        map<string, string> new_2_old_label;
        for (auto& atom : crystal.atoms)
        {   
            string new_label = convert_internal_altloc_label(atom.label);
            new_2_old_label[new_label] = atom.label;
            atom.label = new_label;
        }
        split_with_labels(crystal, ordered_parts);
        for (auto& ordered_part : ordered_parts)
            for (auto& atom : ordered_part)
                atom.first = new_2_old_label[atom.first];

    }


    void split_with_labels(
        const Crystal& crystal,
        std::vector < std::vector<std::pair<std::string, double> > >& _ordered_parts)
    {
        split_with_labels_new_impl(crystal, _ordered_parts);
        return;

        _ordered_parts.clear();
        vector < vector<pair<int, double> > > ordered_parts;

        map<char, int> disorderGroupIdxMap;
        for (int i = 65; i <= 90; i++)
            disorderGroupIdxMap[char(i)] = i - 65;

        set<string> unique_disorder_part_labels;
        vector<string> words;
        vector<int> atom_disorder_part;
        //aletrnativeAtomGroups[label][i].first - disorder group 
        //aletrnativeAtomGroups[label][i].second - atom index
        map<string, vector<pair<int, int> > > aletrnativeAtomGroups;
        vector<int> orderedAtoms;
        int nAtoms = crystal.atoms.size();
        // set only for disordered atoms,
        // needed for weights calculatino
        vector<int> nContatiningConfigurations(nAtoms, 0);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            auto const& atom = crystal.atoms[atomIdx];

            string_utilities::split(atom.label, words, '.');
            int disorder_part = -1;
            if (words.size() == 2)
            {
                if (words[1].size() > 1)
                {
                    string message = "only single character notation is supported for alternative "
                        "conformation coding in atom labels, the label " + atom.label +
                        " cannot be used for automatic fragments generation";
                    on_error::throwException(message, __FILE__, __LINE__);
                }
                else
                {
                    char c = toupper(words[1][0]);
                    if (disorderGroupIdxMap.count(c) == 0)
                    {
                        string message = "only A-Z,a-z characters are supported for alternative "
                            "conformation coding in atom labels, the label " + atom.label +
                            " cannot be used for automatic fragments generation";
                        on_error::throwException(message, __FILE__, __LINE__);
                    }
                    disorder_part = disorderGroupIdxMap[c];
                }
                aletrnativeAtomGroups[words[0]].push_back({ disorder_part , atomIdx });
            }
            else
                orderedAtoms.push_back(atomIdx);

            atom_disorder_part.push_back(disorder_part);
        }

        for (auto& group : aletrnativeAtomGroups)
            sort(group.second.begin(), group.second.end());


        int max_disorder_part = *max_element(atom_disorder_part.begin(), atom_disorder_part.end());
        if (max_disorder_part <= 0)
            return;

        int nParts = max_disorder_part + 1;

        ordered_parts.resize(nParts);
        for (int partIdx = 0; partIdx < nParts; partIdx++)
        {

            Crystal c = crystal;
            c.atoms.clear();

            // add ordered part
            for (int atomIdx : orderedAtoms)
            {
                c.atoms.push_back(crystal.atoms[atomIdx]);
                ordered_parts[partIdx].push_back({ atomIdx, 0.0 });
            }
            // add disordered part
            for (auto& group : aletrnativeAtomGroups)
            {
                int group_size = group.second.size();
                int idx_in_group = 0;
                for (int groupElementIdx = 0; groupElementIdx < group_size; groupElementIdx++)
                    if (group.second[groupElementIdx].first == partIdx)
                        idx_in_group = groupElementIdx;
                int atomIdx = group.second[idx_in_group].second;
                nContatiningConfigurations[atomIdx]++;
                c.atoms.push_back(crystal.atoms[atomIdx]);
                ordered_parts[partIdx].push_back({ atomIdx, 1.0 });
            }
            // find connectivity and use in weight calculations
            vector<vector<pair<int, string> > > connectivity;
            structural_properties::asymmetricUnitConnectivity(c, connectivity, 0.4);
            int nOrdered = orderedAtoms.size();
            for (int i = 0; i < nOrdered; i++)
            {
                double weight = 1.0;
                //bool sameDisorderedPart = false;
                //bool neighbourDisorderedPart = false;
                //set<int> 
                for (int j = 0; j < connectivity[i].size(); j++)
                    if (connectivity[i][j].first >= nOrdered)
                    {
                        int neighbourIdxInOriginalCrystal = ordered_parts[partIdx][connectivity[i][j].first].first;
                        weight *= crystal.atoms[neighbourIdxInOriginalCrystal].occupancy;
                        if (nContatiningConfigurations[neighbourIdxInOriginalCrystal] > 0)
                            weight /= nContatiningConfigurations[neighbourIdxInOriginalCrystal];
                    }
                ordered_parts[partIdx][i].second = weight;
            }

        }

        correct_zero_weights_for_ordered(nAtoms, orderedAtoms, ordered_parts);

        int nOrderedParts = ordered_parts.size();
        _ordered_parts.resize(nOrderedParts);
        for (int partIdx = 0; partIdx < nOrderedParts; partIdx++)
        {
            _ordered_parts[partIdx].resize(ordered_parts[partIdx].size());
            for (int atomIdx = 0; atomIdx < ordered_parts[partIdx].size(); atomIdx++)
            {
                int idxInCrystal = ordered_parts[partIdx][atomIdx].first;
                double weight = ordered_parts[partIdx][atomIdx].second;
                _ordered_parts[partIdx][atomIdx] = { crystal.atoms[idxInCrystal].label, weight };
            }
        }

    }

    void split_with_altlocs(
        const Crystal& crystal,
        const std::vector<char>& _altlocs,
        const std::vector< std::vector<std::pair<int, std::string> > >& connectivity,
        std::vector< std::vector<std::pair<std::string, double> > >& _ordered_parts)
    {
        //_ordered_parts.clear();
        //std::vector<char> altlocs;
        //for(auto c: _altlocs)
        //    altlocs.push_back(toupper(c));

        //vector < vector<pair<int, double> > > ordered_parts;

        //map<char, int> disorderGroupIdxMap;
        //for (int i = 65; i <= 90; i++)
        //    disorderGroupIdxMap[char(i)] = i - 65;

        //set<string> unique_disorder_part_labels;
        //vector<int> atom_disorder_part;
        ////aletrnativeAtomGroups[label][i].first - disorder group 
        ////aletrnativeAtomGroups[label][i].second - atom index
        ////map<string, vector<pair<int, int> > > aletrnativeAtomGroups;
        //map<char, vector<int> > disorderedParts;
        //vector<int> orderedAtoms;
        //int nAtoms = crystal.atoms.size();
        //// set only for disordered atoms,
        //// needed for weights calculatino
        //vector<int> nContatiningConfigurations(nAtoms, 0);
        //for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        //{
        //    auto const& atom = crystal.atoms[atomIdx];

        //    int disorder_part = -1;
        //    if (altlocs[atomIdx] != ' ')
        //    {
        //        char altloc = altlocs[atomIdx];
        //        if (disorderGroupIdxMap.count(altloc) == 0)
        //        {
        //            string message = "only A-Z,a-z characters are supported for alternative "
        //                "conformation coding in atom labels, the label " + atom.label +
        //                " cannot be used for automatic fragments generation";
        //            on_error::throwException(message, __FILE__, __LINE__);
        //        }
        //        disorder_part = disorderGroupIdxMap[altloc];

        //        disorderedParts[altloc].push_back(atomIdx);
        //    }
        //    else
        //        orderedAtoms.push_back(atomIdx);

        //    atom_disorder_part.push_back(disorder_part);
        //}


        //int max_disorder_part = *max_element(atom_disorder_part.begin(), atom_disorder_part.end());
        //if (max_disorder_part <= 0)
        //{
        //    for(auto const &atom: crystal.atoms)
        //        _ordered_parts.push_back({ { atom.label, 1.0 } });
        //    return;
        //}



        ////int nParts = max_disorder_part + 1;

        ////ordered_parts.resize(nParts);
        ////[ordered atom idx][part idx][nerighbour idx]
        //int nOrdered = orderedAtoms.size();
        ////vector<vector<vector<int> > > orderedPartAtomDisorderedNeighbours(nOrdered, vector<vector<int> >(nParts));

        ////--new implementation
        //for (auto disordered_part: disorderedParts)
        //{

        //    Crystal c = crystal;
        //    c.atoms.clear();

        //    // add ordered part
        //    for (int atomIdx : orderedAtoms)
        //        c.atoms.push_back(crystal.atoms[atomIdx]);

        //    // add disordered part
        //    for (auto &atomIdx: disordered_part.second)
        //        c.atoms.push_back(crystal.atoms[atomIdx]);

 
        //    // find connectivity and use in weight calculations
        //    vector<vector<pair<int, string> > > connectivity_ordered;
        //    if(connectivity.empty())
        //        structural_properties::asymmetricUnitConnectivity(c, connectivity_ordered, 0.4);
        //    else
        //    {
        //        vector<int> idxInDisordered = orderedAtoms;
        //        idxInDisordered.insert(idxInDisordered.end(), disordered_part.second.begin(), disordered_part.second.end());
        //        map<int, int> disordered2orderedStructIdxMap;
        //        for(int i=0;i<idxInDisordered.size();i++)
        //            disordered2orderedStructIdxMap[idxInDisordered[i]] = i;
        //        connectivity_ordered.resize(idxInDisordered.size());
        //        for(int i=0;i<idxInDisordered.size();i++)
        //            for(auto const &neighbour: connectivity[idxInDisordered[i]])
        //                if (altlocs[neighbour.first] == ' ' || altlocs[neighbour.first] == disordered_part.first)
        //                {
        //                    int neighbour_idx = disordered2orderedStructIdxMap[neighbour.first];
        //                    string neighbour_symmetry = neighbour.second;
        //                    connectivity_ordered[i].push_back({ neighbour_idx, neighbour_symmetry });
        //                }
        //    }


        //    for (int i = 0; i < nOrdered; i++)
        //        for (int j = 0; j < connectivity_ordered[i].size(); j++)
        //            if (connectivity_ordered[i][j].first >= nOrdered)
        //            {
        //                int neighbourIdxInOriginalCrystal = ordered_parts[partIdx][connectivity[i][j].first].first;
        //                orderedPartAtomDisorderedNeighbours[i][partIdx].push_back(neighbourIdxInOriginalCrystal);
        //            }
        //}
        ////-- eof new implementation


        ////-- eof old implementation

        //// when ordered atom disordered neighbours in each ordered subsystem have the same altloc (e.g. only A, in another only B)
        //// in such a case weight for ordered atom is an average crystal.atoms[neighbourIdxInOriginalCrystal].occupancy
        //// otherwise it is a product of such terms

        //vector<bool> allDiorderedNeighboursFromTheSameGroup(nOrdered);

        //for (int i = 0; i < nOrdered; i++)
        //{
        //    bool sameGroup = true;
        //    for (int partIdx = 0; partIdx < nParts; partIdx++)
        //        if (!orderedPartAtomDisorderedNeighbours[i][partIdx].empty())
        //        {
        //            int disorderGroup = atom_disorder_part[orderedPartAtomDisorderedNeighbours[i][partIdx][0]];
        //            for (int j = 1; j < orderedPartAtomDisorderedNeighbours[i][partIdx].size(); j++)
        //                if (disorderGroup != atom_disorder_part[orderedPartAtomDisorderedNeighbours[i][partIdx][j]])
        //                    sameGroup = false;
        //        }
        //    allDiorderedNeighboursFromTheSameGroup[i] = sameGroup;
        //}

        //for (int partIdx = 0; partIdx < nParts; partIdx++)
        //    for (int i = 0; i < nOrdered; i++)
        //    {
        //        double weight = 1.0;
        //        if (allDiorderedNeighboursFromTheSameGroup[i])
        //        {
        //            weight = 0.0;
        //            for (int j = 0; j < orderedPartAtomDisorderedNeighbours[i][partIdx].size(); j++)
        //            {
        //                int neighbourIdxInOriginalCrystal = orderedPartAtomDisorderedNeighbours[i][partIdx][j];
        //                weight += crystal.atoms[neighbourIdxInOriginalCrystal].occupancy / nContatiningConfigurations[neighbourIdxInOriginalCrystal];
        //            }
        //            if (weight == 0.0)
        //                weight = 1.0;
        //        }
        //        else
        //            for (int j = 0; j < orderedPartAtomDisorderedNeighbours[i][partIdx].size(); j++)
        //            {
        //                int neighbourIdxInOriginalCrystal = orderedPartAtomDisorderedNeighbours[i][partIdx][j];
        //                weight *= crystal.atoms[neighbourIdxInOriginalCrystal].occupancy;
        //                if (nContatiningConfigurations[neighbourIdxInOriginalCrystal] > 0)
        //                    weight /= nContatiningConfigurations[neighbourIdxInOriginalCrystal];
        //            }

        //        ordered_parts[partIdx][i].second = weight;
        //    }

        //correct_zero_weights_for_ordered(nAtoms, orderedAtoms, ordered_parts);

        //int nOrderedParts = ordered_parts.size();
        //_ordered_parts.resize(nOrderedParts);
        //for (int partIdx = 0; partIdx < nOrderedParts; partIdx++)
        //{
        //    _ordered_parts[partIdx].resize(ordered_parts[partIdx].size());
        //    for (int atomIdx = 0; atomIdx < ordered_parts[partIdx].size(); atomIdx++)
        //    {
        //        int idxInCrystal = ordered_parts[partIdx][atomIdx].first;
        //        double weight = ordered_parts[partIdx][atomIdx].second;
        //        _ordered_parts[partIdx][atomIdx] = { crystal.atoms[idxInCrystal].label, weight };
        //    }
        //}


    }


    void split_with_labels_new_impl(
        const Crystal& crystal,
        std::vector< std::vector<std::pair<std::string, double> > >& _ordered_parts)
    {
        _ordered_parts.clear();
        vector < vector<pair<int, double> > > ordered_parts;

        map<char, int> disorderGroupIdxMap;
        for (int i = 65; i <= 90; i++)
            disorderGroupIdxMap[char(i)] = i - 65;

        set<string> unique_disorder_part_labels;
        vector<string> words;
        vector<int> atom_disorder_part;
        //aletrnativeAtomGroups[label][i].first - disorder group 
        //aletrnativeAtomGroups[label][i].second - atom index
        map<string, vector<pair<int, int> > > aletrnativeAtomGroups;
        vector<int> orderedAtoms;
        int nAtoms = crystal.atoms.size();
        // set only for disordered atoms,
        // needed for weights calculatino
        vector<int> nContatiningConfigurations(nAtoms, 0);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            auto const& atom = crystal.atoms[atomIdx];

            string_utilities::split(atom.label, words, '.');
            int disorder_part = -1;
            if (words.size() == 2)
            {
                if (words[1].size() > 1)
                {
                    string message = "only single character notation is supported for alternative "
                        "conformation coding in atom labels, the label " + atom.label +
                        " cannot be used for automatic fragments generation";
                    on_error::throwException(message, __FILE__, __LINE__);
                }
                else
                {
                    char c = toupper(words[1][0]);
                    if (disorderGroupIdxMap.count(c) == 0)
                    {
                        string message = "only A-Z,a-z characters are supported for alternative "
                            "conformation coding in atom labels, the label " + atom.label +
                            " cannot be used for automatic fragments generation";
                        on_error::throwException(message, __FILE__, __LINE__);
                    }
                    disorder_part = disorderGroupIdxMap[c];
                }
                aletrnativeAtomGroups[words[0]].push_back({ disorder_part , atomIdx });
            }
            else
                orderedAtoms.push_back(atomIdx);

            atom_disorder_part.push_back(disorder_part);
        }

        for (auto& group : aletrnativeAtomGroups)
            sort(group.second.begin(), group.second.end());


        int max_disorder_part = *max_element(atom_disorder_part.begin(), atom_disorder_part.end());
        if (max_disorder_part <= 0)
            return;

        int nParts = max_disorder_part + 1;

        ordered_parts.resize(nParts);
        //[ordered atom idx][part idx][nerighbour idx]
        int nOrdered = orderedAtoms.size();
        vector<vector<vector<int> > > orderedPartAtomDisorderedNeighbours(nOrdered, vector<vector<int> >(nParts));
        
        for (int partIdx = 0; partIdx < nParts; partIdx++)
        {
            
            Crystal c = crystal;
            c.atoms.clear();

            // add ordered part
            for (int atomIdx : orderedAtoms)
            {
                c.atoms.push_back(crystal.atoms[atomIdx]);
                ordered_parts[partIdx].push_back({ atomIdx, 0.0 });
            }
            // add disordered part
            for (auto& group : aletrnativeAtomGroups)
            {
                int group_size = group.second.size();
                int idx_in_group = 0;
                for (int groupElementIdx = 0; groupElementIdx < group_size; groupElementIdx++)
                    if (group.second[groupElementIdx].first == partIdx)
                        idx_in_group = groupElementIdx;
                int atomIdx = group.second[idx_in_group].second;
                nContatiningConfigurations[atomIdx]++;
                c.atoms.push_back(crystal.atoms[atomIdx]);
                ordered_parts[partIdx].push_back({ atomIdx, 1.0 });
            }
            // find connectivity and use in weight calculations
            vector<vector<pair<int, string> > > connectivity;
            structural_properties::asymmetricUnitConnectivity(c, connectivity, 0.4);
           

            for (int i = 0; i < nOrdered; i++)
                for (int j = 0; j < connectivity[i].size(); j++)
                    if (connectivity[i][j].first >= nOrdered)
                    {
                        int neighbourIdxInOriginalCrystal = ordered_parts[partIdx][connectivity[i][j].first].first;
                        orderedPartAtomDisorderedNeighbours[i][partIdx].push_back(neighbourIdxInOriginalCrystal);
                    }
        }
        
        // when ordered atom disordered neighbours in each ordered subsystem have the same altloc (e.g. only A, in another only B)
        // in such a case weight for ordered atom is an average crystal.atoms[neighbourIdxInOriginalCrystal].occupancy
        // otherwise it is a product of such terms

        vector<bool> allDiorderedNeighboursFromTheSameGroup(nOrdered);
        
        for (int i = 0; i < nOrdered; i++)
        {
            bool sameGroup = true;
            for (int partIdx = 0; partIdx < nParts; partIdx++)
                if (!orderedPartAtomDisorderedNeighbours[i][partIdx].empty())
                {
                    int disorderGroup = atom_disorder_part[orderedPartAtomDisorderedNeighbours[i][partIdx][0]];
                    for (int j = 1; j < orderedPartAtomDisorderedNeighbours[i][partIdx].size(); j++)
                        if (disorderGroup != atom_disorder_part[orderedPartAtomDisorderedNeighbours[i][partIdx][j]])
                            sameGroup = false;
                }
            allDiorderedNeighboursFromTheSameGroup[i] = sameGroup;
        }

        for (int partIdx = 0; partIdx < nParts; partIdx++)
            for (int i = 0; i < nOrdered; i++)
            {
                double weight = 1.0;
                if (allDiorderedNeighboursFromTheSameGroup[i])
                {
                    weight = 0.0;
                    for (int j = 0; j < orderedPartAtomDisorderedNeighbours[i][partIdx].size(); j++)
                    {
                        int neighbourIdxInOriginalCrystal = orderedPartAtomDisorderedNeighbours[i][partIdx][j];
                        weight += crystal.atoms[neighbourIdxInOriginalCrystal].occupancy / nContatiningConfigurations[neighbourIdxInOriginalCrystal];
                    }
                    if (weight == 0.0)
                        weight = 1.0;
                }
                else
                    for (int j = 0; j < orderedPartAtomDisorderedNeighbours[i][partIdx].size(); j++)
                    {
                        int neighbourIdxInOriginalCrystal = orderedPartAtomDisorderedNeighbours[i][partIdx][j]; 
                        weight *= crystal.atoms[neighbourIdxInOriginalCrystal].occupancy;
                        if (nContatiningConfigurations[neighbourIdxInOriginalCrystal] > 0)
                            weight /= nContatiningConfigurations[neighbourIdxInOriginalCrystal];
                    }

                ordered_parts[partIdx][i].second = weight;
            }

        correct_zero_weights_for_ordered(nAtoms, orderedAtoms, ordered_parts);

        int nOrderedParts = ordered_parts.size();
        _ordered_parts.resize(nOrderedParts);
        for (int partIdx = 0; partIdx < nOrderedParts; partIdx++)
        {
            _ordered_parts[partIdx].resize(ordered_parts[partIdx].size());
            for (int atomIdx = 0; atomIdx < ordered_parts[partIdx].size(); atomIdx++)
            {
                int idxInCrystal = ordered_parts[partIdx][atomIdx].first;
                double weight = ordered_parts[partIdx][atomIdx].second;
                _ordered_parts[partIdx][atomIdx] = { crystal.atoms[idxInCrystal].label, weight };
            }
        }

    }

    void fragments_from_json(
        const Crystal& crystal,
        const nlohmann::json& fragments_json,
        std::vector<Fragment>& fragments)
    {
        fragments.clear();

        if (!fragments_json.is_array())
            on_error::throwException("expected array when processing information on ordered fragmnets from json data",
                __FILE__, __LINE__);

        map<string, int> label2idx;
        for (int i = 0; i < crystal.atoms.size(); i++)
            label2idx[crystal.atoms[i].label] = i;

        for (auto const& fragmentJson : fragments_json)
        {
            Fragment fragment;
            for (auto const& weight_atom_symmOp : fragmentJson)
            {
                string atomLabel = weight_atom_symmOp[1];
                double weight = weight_atom_symmOp[0];
                string symmOp = "x,y,z";
                if(weight_atom_symmOp.size()==3)
                    symmOp = weight_atom_symmOp[2]; 
                
                int atomIdx = 0;
                if (label2idx.count(atomLabel) == 0)
                    on_error::throwException("unknown atom label '" + atomLabel + "' when processing information on ordered fragmnets from json data", __FILE__, __LINE__);
                fragment.atomList.push_back({ label2idx[atomLabel], symmOp });
                fragment.atomRelativeWeights.push_back(weight);
            }
            fragments.push_back(fragment);
        }

    }


    void substructures_from_json(
        //const Crystal& crystal,
        const nlohmann::json& disorder_parts,
        std::vector< std::vector<std::pair<std::string, double> > >& substructures)
    {
        substructures.clear();

        if (!disorder_parts.is_array())
            on_error::throwException("expected array when processing information on  ordered substructures from json data", 
                __FILE__, __LINE__);
        
        for (auto const& substructureJson : disorder_parts)
        {
            vector < pair<string, double> > substructure;
            for (auto const& atom_weight : substructureJson)
            {
                string atomLabel = atom_weight[0];
                double weight = atom_weight[1];
                substructure.push_back({ atomLabel, weight });
            }
            substructures.push_back(substructure);
        }
    }

    void split_structure_json_str(
        const Crystal& crystal,
        const nlohmann::json& disorder_groups_json,
        std::vector< std::vector<std::pair<std::string, double> > >& substructures)
    {
        substructures.clear();
        vector< vector< vector<int> > > disorder_groups_int;
        if (!disorder_groups_json.is_array())
            on_error::throwException("expected array when processing information on  ordered substructures from json data",
                __FILE__, __LINE__);

        map<string, int> label2int;
        for (int i = 0; i < crystal.atoms.size(); i++)
            label2int[crystal.atoms[i].label] = i;

        for (auto const& disorder_assembly : disorder_groups_json)
        {
            vector< vector<int> > assembly;
            for (auto const& disorder_group : disorder_assembly)
            {
                vector<int> part;
                for (auto const& atomIdx : disorder_group.items())
                    part.push_back(label2int[atomIdx.value().get<string>()]);
                assembly.push_back(part);
            }
            disorder_groups_int.push_back(assembly);
        }
        split_structure(crystal, disorder_groups_int, substructures);
    }


void split_structure(
    const Crystal& crystal,
    const vector< vector< vector<int> > > & disorder_groups,
    vector< vector<pair<string, double> > >& _substructures)
    {
        _substructures.clear();        
        int nParts = 0;
        for (auto const& group : disorder_groups)
            if (group.size() > nParts)
                nParts = group.size();
        

        set<int> orderedAtomsSet;
        int nAtoms = crystal.atoms.size();
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            orderedAtomsSet.insert(atomIdx);
        for (auto const& disorder_group : disorder_groups)
            for(auto const& same_altloc_group : disorder_group)
                for (int atom : same_altloc_group)
                    orderedAtomsSet.erase(atom);
        
        vector<int> orderedAtoms(orderedAtomsSet.begin(), orderedAtomsSet.end());

        vector< vector<pair<int, double> > > substructures;
        substructures.resize(nParts);
        for (int partIdx = 0; partIdx < nParts; partIdx++)
        {

            Crystal c = crystal;
            c.atoms.clear();

            // add ordered part
            for (int atomIdx : orderedAtoms)
            {
                c.atoms.push_back(crystal.atoms[atomIdx]);
                substructures[partIdx].push_back({ atomIdx, 1.0 });
            }
            // add disordered part
            for (auto const& alternative_disorder_groups : disorder_groups)
            {
                //number of configurations at given site
                //int max_group_idx = alternative_disorder_groups.size();
                //int groupIdx = min(partIdx, max_group_idx);
      
                int groupIdx = partIdx;
                if (groupIdx >= alternative_disorder_groups.size())
                    groupIdx = 0;

                for (int atomIdx : alternative_disorder_groups[groupIdx])
                {
                    c.atoms.push_back(crystal.atoms[atomIdx]);
                    substructures[partIdx].push_back({ atomIdx, 1.0 });
                }
            }
            // find connectivity and use in weight calculations
            vector<vector<pair<int, string> > > connectivity;
            structural_properties::asymmetricUnitConnectivity(c, connectivity, 0.4);
            int nOrdered = orderedAtoms.size();
            for (int i = 0; i < nOrdered; i++)
            {
                double weight = 1.0;
                for (int j = 0; j < connectivity[i].size(); j++)
                    if (connectivity[i][j].first >= nOrdered)
                    {
                        int neighbourIdxInSubstructure = connectivity[i][j].first;
                        int neighbourIdxInOriginalCrystal = substructures[partIdx][neighbourIdxInSubstructure].first;
                        weight *= crystal.atoms[neighbourIdxInOriginalCrystal].occupancy;
                    }
                substructures[partIdx][i].second = weight;
            }

        }

        correct_zero_weights_for_ordered(nAtoms, orderedAtoms, substructures);

        convert_ordered_parts_list(crystal, substructures, _substructures);
    }


    void convert_ordered_parts_list(
        const Crystal& crystal,
        const std::vector< std::vector<std::pair<int, double> > >& ordered_parts_int,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts_str)
    {
        int nOrderedParts = ordered_parts_int.size();
        ordered_parts_str.resize(nOrderedParts);
        for (int partIdx = 0; partIdx < nOrderedParts; partIdx++)
        {
            ordered_parts_str[partIdx].resize(ordered_parts_int[partIdx].size());
            for (int atomIdx = 0; atomIdx < ordered_parts_int[partIdx].size(); atomIdx++)
            {
                int idxInCrystal = ordered_parts_int[partIdx][atomIdx].first;
                double weight = ordered_parts_int[partIdx][atomIdx].second;
                ordered_parts_str[partIdx][atomIdx] = { crystal.atoms[idxInCrystal].label, weight };
            }
        }

    }


    void from_file(
        const string& fileName,
        vector< vector<pair<string, double> > >& ordered_parts)
    {
        ordered_parts.clear();

        ifstream in(fileName);
        if (!in.good())
            on_error::throwException("cannot open fragments file " + fileName, __FILE__, __LINE__);

        string line;

        while (getline(in, line))
        {
            if (line.empty())
                continue;
            if (line[0] == '#')
                continue; // skip comments

            vector < pair<string, double> > atomList;
            int nAtoms = stoi(line);
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                getline(in, line);
                vector<string> words;
                string_utilities::split(line, words, ' ');
                double weight = 1.0;
                if (words.size() == 2)
                    weight = stod(words[1]);
                atomList.push_back({ words[0], weight });
            }
            ordered_parts.push_back(atomList);
        }

    }

    void from_file(
        const Crystal& crystal,
        const std::string& fileName,
        std::vector<Fragment>& fragments)
    {
        fragments.clear();
        ifstream in(fileName);

        if (!in.good())
            on_error::throwException("cannot open disordered crystal fragments file " + fileName, __FILE__, __LINE__);

        map<string, int> label2idx;
        for (int i = 0; i < crystal.atoms.size(); i++)
            label2idx[crystal.atoms[i].label] = i;

        string line;
        int nFragments;
        in >> nFragments;
        fragments.resize(nFragments);
        for (int fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
        {
            int atomIdx, nAtoms;
            in >> nAtoms;
            getline(in, line);

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                double weight;
                string label;
                string symmOp = "x,y,z";
                vector<string> words;
                getline(in, line);
                string_utilities::split(line, words);
                if (words.size() < 2)
                    on_error::throwException("invalid format of disordered crystal fragments file", __FILE__, __LINE__);
                weight = stod(words[0]);
                label = words[1];
                if (words.size() == 3)
                    symmOp = words[2];
                fragments[fragmentIdx].atomList.push_back({ label2idx[label],symmOp });
                fragments[fragmentIdx].atomRelativeWeights.push_back(weight);
            }
        }
    }


}
}

