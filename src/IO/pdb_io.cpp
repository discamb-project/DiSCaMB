#include "discamb/IO/pdb_io.h"
#include "discamb/BasicChemistry/periodic_table.h"


#include <fstream>
#include <iomanip>

using namespace std;

namespace discamb {

    namespace pdb_io {

        void write(
            const std::string& fileName,
            const  std::vector<std::string>& symbols,
            const std::vector<Vector3d>& positions)
        {
            int nAtoms = symbols.size();
            ofstream out(fileName);
            for (int i = 0; i < nAtoms; i++)
            {
                out << "HETATM" << setw(5) << i + 1 << setw(4) << symbols[i]
                    << "   LIG     1  ";
                for (int j = 0; j < 3; j++)
                    out << setw(8) << setprecision(3) << fixed << positions[i][j];

                out << "                       ";
                if (symbols[i].size() == 1)
                    out << " ";
                out << symbols[i] << "\n";
            }
            out.close();
        }

        void write(
            const std::string& fileName,
            const std::vector<int>& atomicNumber,
            const std::vector<Vector3d>& position)
        {

            vector<string> symbols;
            for (auto const& z : atomicNumber)
                symbols.push_back(periodic_table::symbol(z));
            write(fileName, symbols, position);
        }


    }
}

