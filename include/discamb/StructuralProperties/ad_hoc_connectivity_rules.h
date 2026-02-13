#pragma once
#include "discamb/MathUtilities/Vector3.h"
#include <vector>
#include <string>

namespace discamb {

    namespace ad_hoc_connectivity_rules {

        enum class Preset{All, None};

        Preset preset_from_string(const std::string& s);
        std::string to_string(const Preset &preset);

        void apply_all(const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity);

        void apply(const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity,
            Preset preset);

        /**
        disconnects possible spurious contact between C1 and C3 
           C5
         /   \
        C1-C2-C3
         \   /
           C4
        */

        void disconnect_CC_in_trigonal_bypyramid(
            const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity);
    }

}
