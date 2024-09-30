#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <string>
#include <vector>

namespace discamb{
    /**
    * \addtogroup IO IO
    * @{
    */

    namespace gaussian_io {

        struct Cube {
            int nx = 0;
            int ny = 0;
            int nz = 0;
            Vector3d origin, directions[3];
            std::vector<Vector3d> positions;
            std::vector<int> atomic_number;
            std::vector<double> atomic_charge;
            std::vector<std::vector<std::vector<double> > > data;
            std::string headerLine1, headerLine2;
        };

        void readCube(const std::string& fName, Cube& cube);
        //void readCubeHeader(const std::string &fName, Cube &cube);
        void writeCube(const std::string& fName, const Cube& cube);


    }
    /**@}*/
};
