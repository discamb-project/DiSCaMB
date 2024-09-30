#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <string>

namespace discamb {

    struct MoleculeData {
        std::string comment;
        std::vector<std::string> atomLabels;
        std::vector <Vector3d> atomPositions;
        std::vector <int> atomicNumbers;
        std::vector<std::pair<int, int> > bonds;
    };

}

