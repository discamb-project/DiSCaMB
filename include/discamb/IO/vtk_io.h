#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/IO/gaussian_io.h"

#include <vector>
#include <string>

namespace discamb{

    /**
    * \addtogroup IO IO
    * @{
    */


	namespace vtk_io {

		void save_point_data(const std::string& fileName, const std::vector<Vector3d>& positions, const std::vector<double>& values, const std::string &data_name = std::string("data"));
        void save_regular_grid_data(
            const std::string& fileName, 
            const Vector3d& origin,
            const Vector3d& step1, //dx
            const Vector3d& step2, //dy
            const Vector3d& step3, //dx
            const std::vector< std::vector< std::vector<double> > >& values); //[x][y][z]
        void saveGaussianCubeData(const gaussian_io::Cube& cube, const std::string& vtkFileName, const std::string& molFileName = std::string());
	}
    /**@}*/
}
