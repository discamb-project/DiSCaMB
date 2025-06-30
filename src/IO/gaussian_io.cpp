#include <string>
#include <vector>
#include <fstream>


#include "discamb/IO/gaussian_io.h"

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"


using namespace std;

namespace discamb{
    namespace gaussian_io {

        void readCube(const string& fName, Cube& cube)
        {
            ifstream in(fName);
            int i, j, k, nAtoms;
            string line;
            vector<string> words;
            Vector3d position;
            if (!in.good())
                on_error::throwException(string("cannot open file ") + fName, __FILE__, __LINE__);

            getline(in, cube.headerLine1);
            getline(in, cube.headerLine2);
            
            //in >> nAtoms >> cube.origin[0] >> cube.origin[1] >> cube.origin[2];
            getline(in, line);
            string_utilities::split(line, words);
            nAtoms = stoi(words[0]);
            cube.origin.set(stod(words[1]), stod(words[2]), stod(words[3]));

            in >> cube.nx >> cube.directions[0][0] >> cube.directions[0][1] >> cube.directions[0][2];
            in >> cube.ny >> cube.directions[1][0] >> cube.directions[1][1] >> cube.directions[1][2];
            in >> cube.nz >> cube.directions[2][0] >> cube.directions[2][1] >> cube.directions[2][2];

            cube.atomic_charge.resize(nAtoms);
            cube.atomic_number.resize(nAtoms);
            cube.positions.resize(nAtoms);

            for (i = 0; i < nAtoms; i++)
                in >> cube.atomic_number[i] >> cube.atomic_charge[i] >> cube.positions[i][0] >> cube.positions[i][1] >> cube.positions[i][2];

            cube.data.clear();
            cube.data.resize(cube.nx, vector<vector<double> >(cube.ny, vector<double>(cube.nz, 0)));

            for (i = 0; i < cube.nx; i++)
                for (j = 0; j < cube.ny; j++)
                    for (k = 0; k < cube.nz; k++)
                        in >> cube.data[i][j][k];

            in.close();
        }


        void writeCube(const string& fName, const Cube& cube)
        {
            ofstream out(fName);

            int i, j, k, n, nAtoms = cube.atomic_charge.size();

            out << cube.headerLine1 << endl;
            out << cube.headerLine2 << endl;

            out << nAtoms << "  " << cube.origin[0] << " " << cube.origin[1] << " " << cube.origin[2] << endl;

            out << cube.nx << " " << cube.directions[0][0] << " " << cube.directions[0][1] << " " << cube.directions[0][2] << endl;
            out << cube.ny << " " << cube.directions[1][0] << " " << cube.directions[1][1] << " " << cube.directions[1][2] << endl;
            out << cube.nz << " " << cube.directions[2][0] << " " << cube.directions[2][1] << " " << cube.directions[2][2] << endl;



            for (i = 0; i < nAtoms; i++)
                out << cube.atomic_number[i] << " " << cube.atomic_charge[i] << " " << cube.positions[i][0] << " " << cube.positions[i][1] << " " << cube.positions[i][2] << endl;

            n = 0;
            for (i = 0; i < cube.nx; i++)
                for (j = 0; j < cube.ny; j++)
                    for (k = 0; k < cube.nz; k++)
                    {
                        out << setw(18) << setprecision(9) << cube.data[i][j][k] << " ";
                        n++;
                        if (n % 6 == 0)
                            out << "\n";
                    }

            out.close();
        }
    }
} // namespace discamb
