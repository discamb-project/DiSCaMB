#include "discamb/IO/molecule_io.h"
#include "discamb/IO/mol2_io.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/BasicUtilities/on_error.h"

#include <filesystem>

using namespace std;

namespace discamb {

    namespace molecule_io {

        void read_molecular_structure(
            const std::string& fileName,
            MoleculeData& data)
        {
            if (!filesystem::exists(fileName))
                on_error::throwException("the file to open '" + fileName + "' does not exist", __FILE__, __LINE__);

            string extension = filesystem::path(fileName).extension().string();

            if (extension == ".xyz")
            {
                xyz_io::readXyz(fileName, data);
                return;
            }

            if (extension == ".mol2")
            {
                mol2_io::read(fileName, data);
                return;
            }

            on_error::throwException("not supported file extension '" + extension + "' encountered when attempting to read a file '" + fileName + "'", __FILE__, __LINE__);
        }

        void read_molecular_structure(
            const std::string& fileName,
            vector<Vector3d>& positions,
            vector<int>& atomic_numbers)
        {
            MoleculeData data;
            read_molecular_structure(fileName, data);
            positions = data.atomPositions;
            atomic_numbers = data.atomicNumbers;
        }

    }

}
