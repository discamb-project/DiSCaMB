#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


    namespace pdb_io {
        
        void write(const std::string& fileName, const std::vector<int>& atomicNumber, const std::vector<Vector3d>& position);
        void write(const std::string& fileName, const  std::vector<std::string>& symbol, const std::vector<Vector3d>& position);
    }
    /**@}*/
}

