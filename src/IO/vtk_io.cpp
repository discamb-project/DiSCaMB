#include "discamb/IO/vtk_io.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/Constants.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/IO/pdb_io.h"


#include <fstream>

using namespace std;

namespace discamb{
	namespace vtk_io {
		void save_point_data(
			const std::string& fileName,
			const std::vector<Vector3d>& positions,
			const std::vector<double>& values,
			const std::string& dataName)
		{
			ofstream out(fileName);

			if (!out.good())
				on_error::throwException(string("can not save to vtk file '") + fileName + string("'"), __FILE__, __LINE__);

			int nPoints = positions.size();
			if (positions.size() != values.size())
				on_error::throwException("number of points does not match number of values when attempting to write to vtk file", __FILE__, __LINE__);

			out << "# vtk DataFile Version 2.0\n"
				"atomic_density\n"
				"ASCII\n"
				"DATASET POLYDATA\n"
				"POINTS ";
			out << nPoints << " double\n";

			for (auto& point : positions)
				out << point[0] << " " << point[1] << " " << point[2] << "\n";

			out << "\nVERTICES " << nPoints << " " << 2 * nPoints << "\n";

			for (int pointIdx = 0; pointIdx < nPoints; pointIdx++)
				out << "1 " << pointIdx << "\n";

			out << "\nPOINT_DATA " << nPoints << "\n";
			out << "SCALARS " << dataName << " double 1\n"
				   "LOOKUP_TABLE default\n";

			for (auto x : values)
				out << x << "\n";
			out.close();

		}
        void save_regular_grid_data(
            const std::string& fileName,
            const Vector3d& origin,
            const Vector3d& step1, //dx
            const Vector3d& step2, //dy
            const Vector3d& step3, //dx
            const std::vector< std::vector< std::vector<double> > >& values) //[x][y][z]
        {
            ofstream out(fileName);

            if (!out.good())
                on_error::throwException(string("can not save to vtk file '") + fileName + string("'"), __FILE__, __LINE__);

            int nx, ny, nz;
            if (values.empty())
                on_error::throwException("invalid dimansions of interpolation grid (should be at least 1x1x1", __FILE__, __LINE__);
            if (values[0].empty())
                on_error::throwException("invalid dimansions of interpolation grid (should be at least 1x1x1", __FILE__, __LINE__);
            if (values[0][0].empty())
                on_error::throwException("invalid dimansions of interpolation grid (should be at least 1x1x1", __FILE__, __LINE__);
            nx = values.size();
            ny = values[0].size();
            nz = values[0][0].size();


            Vector3d r;
            int i, j, k;
            out << "# vtk DataFile Version 3.0\n"
                << "vtk output\n"
                << "ASCII\n"
                << "DATASET STRUCTURED_GRID\n"
                << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n"
                << "POINTS " << nx*ny*nz << " double\n";
            for (i = 0; i < nz; i++)
                for (j = 0; j < ny; j++)
                    for (k = 0; k < nx; k++)
                    {
                        r = origin + double(k) * step1 + double(j) * step2 + double(i) * step3;
                        out << r.x << " " << r.y << " " << r.z << "\n";
                    }

            out << "\n"
                << "POINT_DATA " << nx * ny * nz << "\n"
                << "SCALARS scalars double 1\n"
                << "LOOKUP_TABLE default\n";

            for (k = 0; k < nz; k++)
                for (j = 0; j < ny; j++)
                    for (i = 0; i < nx; i++)
                        out << values[i][j][k] << endl;

            out.close();

        }

        void saveGaussianCubeData(
            const gaussian_io::Cube& cube, 
            const std::string& vtkFileName,
            const std::string& molFileName)
        {
            vtk_io::save_regular_grid_data(vtkFileName, cube.origin / constants::Angstrom,
                cube.directions[0] / constants::Angstrom, cube.directions[1] / constants::Angstrom, cube.directions[2] / constants::Angstrom, cube.data);


            if (molFileName.empty())
                return;

            auto positions = cube.positions;
            for (auto& r : positions)
                r /= constants::Angstrom;

            if (molFileName.find(".pdb") == string::npos)
                xyz_io::writeXyz(molFileName, cube.atomic_number, positions);
            else
                pdb_io::write(molFileName, cube.atomic_number, positions);

        
        }

	}

}

