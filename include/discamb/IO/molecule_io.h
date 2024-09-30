#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicChemistry/MoleculeData.h"

#include <vector>

namespace discamb {
    /**
    * \addtogroup IO IO
    * @{
    */

    namespace molecule_io {
        
        void read_molecular_structure(const std::string& fileName, std::vector<Vector3d>& positions, std::vector<int>& atomic_numbers);
        void read_molecular_structure(const std::string& fileName, MoleculeData &data);
    }
    /**@}*/
}
