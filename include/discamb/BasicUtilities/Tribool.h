#pragma once

#include <string>

namespace discamb{
    enum class Tribool {True, False, Undefined};

    inline std::string tribool_to_string(const Tribool& tribool)
    {
        if (tribool == Tribool::False)
            return "FALSE";
        if (tribool == Tribool::True)
            return "TRUE";
        return "UNDEFINED";
    }
}


