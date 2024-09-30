#include <vector>
#include <string>

namespace discamb{

    namespace tham_io {

        void readAtomTypeElectronDensity(
            const std::string& fName,
            std::vector<double>& electronDensity,
            int& angularGridSize,
            int& radialGridSize,
            int& atomicNumber);

    }

}

