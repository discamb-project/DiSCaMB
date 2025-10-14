#pragma once

#include "discamb/CrystalStructure/Crystal.h"

#include "json.hpp"

namespace discamb{



namespace disordered_structure_fragments{

    struct Fragment {
        // list of atoms in fragment given by atom index in asymmetric unit + symmetry operation
        std::vector < std::pair<int, std::string> >  atomList;
        /*
        atomRelativeWeights[i] - relative weight for i-th atom i the fragment (defined by atomList[i])
        if relative weight is 0 i-th atom do not contribute to form factor of given atom

        if atom appears in multiple fragments then the contribution from
        each fragment is weighted by this factor (the weights are rescaled to sum up to 1)
        e.g the form factor for atom which appears in two fragments is given by:
        f_a = (w1 f_a_1 + w2 f_a_2)/(w1+w2)
        w1, w2 - relative weights,
        f_a_1, f_a_2 - atomic form factors for the atom in fragments 1 (f_a_1) and 2 (f_a_2),
        they may differ due to different local coordinate systems
        */
        std::vector<double> atomRelativeWeights;
    };

    //void fragments_from_labels(const Crystal& crystal, std::vector<Fragment>& fragments);
    void split_with_labels_0(
        const Crystal& crystal,
        std::vector< std::vector<std::pair<std::string, double> > > & ordered_parts);

    void split_with_labels(
        const Crystal& crystal,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts);
    // "H2.B   1    A    H105    Z N      1    A X CA     1    A"
    void split_with_labels_internal_altloc(
        const Crystal& crystal,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts);

    void split_with_altlocs(
        const Crystal& crystal,
        const std::vector<char>& altlocs,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts);


    void split_with_labels_new_impl(
        const Crystal& crystal,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts);

    /*
    populates ordered_parts which then can be used in TaamSfCalculator constructor
    TaamSfCalculator(const Crystal& crystal, const TaamSfCalculatorSettings &settings);
    by assigning settings.orderedSubcrystalAtoms = substructures
    ordered_parts is not populated if disorder_groups is empty
    */
    void split_structure(
        const Crystal& crystal,
        const std::vector< std::vector<std::vector<int> > >& disorder_groups,
        std::vector< std::vector<std::pair<std::string, double> > >& substructures);

    void split_structure_json_int(
        const Crystal& crystal,
        const nlohmann::json & disorder_groups,
        std::vector< std::vector<std::pair<std::string, double> > >& substructures);

    void split_structure_json_str(
        const Crystal& crystal,
        const nlohmann::json& disorder_groups,
        std::vector< std::vector<std::pair<std::string, double> > >& substructures);


    void substructures_from_json(
        const nlohmann::json& disorder_parts,
        std::vector< std::vector<std::pair<std::string, double> > >& substructures);

    void fragments_from_json(
        const Crystal& crystal,
        const nlohmann::json& fragments_json,
        std::vector<Fragment>& fragments);

    void convert_ordered_parts_list(
        const Crystal &crystal,
        const std::vector< std::vector<std::pair<int, double> > >& ordered_parts_int,
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts_str);

    /*
    void convert_list(
        const Crystal& crystal, 
        const std::vector< std::vector<std::pair<int, double> > >& ordered_parts_str, 
        std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts_int);
    */
    void from_file(const std::string& fileName, std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts);
    void from_file(const Crystal& crystal, const std::string& fileName, std::vector<Fragment>& fragments);
}
}

