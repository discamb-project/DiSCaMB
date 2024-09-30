#include "discamb/IO/tham_io.h"
#include "discamb/BasicUtilities/OnError.h"

#include <fstream>

using namespace std;

namespace discamb {

    namespace tham_io {

        void readAtomTypeElectronDensity(
            const std::string& fName,
            std::vector<double>& electronDensity,
            int& angularGridSize,
            int& radialGridSize,
            int& atomicNumber)
        {
            ifstream in(fName);
            if (!in.good())
                on_error::throwException("cannot read file '" + fName + "'", __FILE__, __LINE__);
            string s;
            in >> s >> s >> atomicNumber >> s >> s >> angularGridSize >> s >> s >> radialGridSize;
            int gridSize = angularGridSize * radialGridSize;
            electronDensity.clear();
            electronDensity.resize(gridSize);

            for (int i = 0; i < gridSize; i++)
                in >> electronDensity[i];

            in.close();

        }

    }

}
