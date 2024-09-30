#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <string>


namespace discamb {
    namespace rsp_io {
        void read(const std::string& fileName, std::vector<Vector3d>& q);
    }
}

