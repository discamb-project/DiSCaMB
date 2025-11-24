#pragma once

#include "json.hpp"

namespace discamb{
    struct MacromolecularStructuralInformation {
        std::vector<std::vector<std::pair<int, std::string> > > connectivity;
        std::vector<std::vector<std::pair<int, std::string> > > planes;
        std::vector<char> altlocs;
        std::vector<int> residueSequenceNumbers;
        std::vector<std::string> residueNames;
        /* like in pdb e.g. HB3, CA, OG1 */
        std::vector<std::string> atomNames; 
        void set(const nlohmann::json &data);
        void toOrderedSubcrystalAtoms(std::vector<std::vector<std::pair<std::string, double> > >& orderedSubcrystalAtoms) const;
    };
}

