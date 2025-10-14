#pragma once

#include "json.hpp"

namespace discamb{
    struct MacromolecularStructuralInformation {
        std::vector<std::vector<std::pair<int, std::string> > > connectivity;
        std::vector<std::vector<std::pair<int, std::string> > > planes;
        std::vector<char> altlocs;
        void set(const nlohmann::json &data);
    };
}

