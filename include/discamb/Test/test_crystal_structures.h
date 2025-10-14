#pragma once

#include "discamb/CrystalStructure/Crystal.h"

namespace discamb{
    namespace test_crystal_structures{
    
        std::vector<std::string> get_test_crystal_structure_names();
        void get_test_crystal_structure(const std::string &name, Crystal& crystal);

    }
}

