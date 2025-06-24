#pragma once

#include "discamb/CrystalStructure/Crystal.h"

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
    void split_with_labels(const Crystal& crystal, std::vector< std::vector<std::pair<std::string, double> > > & ordered_parts);
    void from_file(const std::string& fileName, std::vector< std::vector<std::pair<std::string, double> > >& ordered_parts);
    void from_file(const Crystal& crystal, const std::string& fileName, std::vector<Fragment>& fragments);
}
}

