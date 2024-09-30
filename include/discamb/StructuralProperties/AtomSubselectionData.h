#pragma once

#include <vector>
#include "discamb/MathUtilities/Vector3.h"

namespace discamb{

    struct AtomSubselectionData {
        std::vector<Vector3i> periodicDirections;
        bool allSymmetryEquivalent = false;
        /**
        atomList[i].first  - atomic label
        atomList[i].second - symmetry operation (in X,Y,Z format)
        */
        std::vector<std::pair<std::string, std::string> > atomList;

    };
}

