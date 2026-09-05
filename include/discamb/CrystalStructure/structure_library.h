#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include <set>

namespace discamb {

    namespace structure_library {
        
        bool getStructure(const std::string& name, Crystal& structure);
        bool hasStructure(const std::string& name);
        void structureList(std::set<std::string>& structureList);

    }
}
