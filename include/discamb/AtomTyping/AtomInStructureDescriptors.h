#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicUtilities/Tribool.h"

#include <vector>
#include <cstddef>
#include <optional>
#include <set>


namespace discamb {

/**
* \addtogroup AtomTyping
* @{
*/


    struct AtomInStructureDescriptors {
        int atomicNumber;
        int nNeighbors;
        std::multiset<int> neighborsFormula;
        Tribool planar;
        double planarityDisplacementEsd;
        std::vector<int> planarRingsSizes;
        std::vector<int> planarRingsIndices;
        int n3memberRings, n4memberRings;
        std::vector<int> disorderGroups; // disorder groups which have to be present in the molecule containing the atom
        std::string label;
		Vector3d position;
    };

    /**@}*/

}
