#include "discamb/IO/structure_io.h"
#include "discamb/IO/shelx_io.h"
#include "discamb/IO/cif_io.h"
#include "discamb/IO/NativeIAM_Reader.h"


#include <filesystem>

namespace fs = std::filesystem;

using namespace std;

namespace discamb {
    namespace structure_io {
        void read_structure(
            const std::string& fileName,
            Crystal& crystal)
        {
            string suffix = fs::path(fileName).extension().string();

            if (suffix == string(".cif"))
            {
                vector<cif_io::DataSet> data;
                cif_io::readCif(fileName, data);
                crystal.atoms.clear();
                int datasetIdx = 0;
                while (crystal.atoms.empty() && datasetIdx < data.size())
                {
                    cif_io::cifDataToCrystal(data[0], crystal);
                    datasetIdx++;
                }

                return;
            }

            if (suffix == string(".res") || suffix == string(".ins"))
                shelx_io::read(fileName, crystal);

            if (suffix == ".dis")
            {
                NativeIAM_Reader reader;
                ScatteringParameters params;
                reader.read(fileName, crystal, params);
            }
        }

        void write_structure(
            const std::string& fileName,
            const Crystal& crystal)
        {
            string suffix = fs::path(fileName).extension().string();

            if (suffix == string(".cif"))
            {
                cif_io::saveCif(fileName, crystal);
                return;
            }

            if (suffix == string(".res") || suffix == string(".ins"))
                shelx_io::save(fileName, crystal);

        }

    }
}

