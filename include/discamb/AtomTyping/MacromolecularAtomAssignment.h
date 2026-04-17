#pragma once

#include "discamb/AtomTyping/LocalCoordinateSystem.h"

#include <string>

namespace discamb{
    struct MacromolecularAtomAssignment{
        std::string type_name;
        LocalCoordinateSystem<std::pair<std::string, int> > lcs;
        void set(const std::string s);
    };
}

