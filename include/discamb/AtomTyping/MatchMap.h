#pragma once

#include <vector>
#include <map>
#include <string>

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    struct MatchMap
    {
        /**
        atomMatch[i] - number of atom in structure to which i-th atom in atom type (UBDB_AtomType) was mapped
        */
        std::vector<int> atomMatch;
        /**
        maps entries from UBDB_AtomType::ringLabels to index of ring in structure
        */
        std::map<std::string, int> labeledRingsMap;
        /**
        unnamedRingsMap[i][j] is an index of ring in structure corresponding to ring UBDB_AtomType::atoms[i].nonLabeledContainingRings[j]
        */
        std::vector<std::vector<int> > unnamedRingsMap;
    };
    /**@}*/
}