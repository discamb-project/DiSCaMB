
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/SfCalculator.h"
#include <fstream>

#include "json.hpp"

using namespace std;
using namespace discamb;


void test_next_gen_taam(
    const string& structureFile,
    const string& hklFile)
{
}


int main(int argc, char* argv[])
{

    try {
        vector<string> arguments, options;
        parse_cmd::get_args_and_options(argc, argv, arguments, options);

        if (arguments.size() < 2)
            on_error::throwException("expected structure and hkl file\n", __FILE__, __LINE__);

        string structureFile = arguments[0];
        string hklFile = arguments[1];

        Crystal crystal;
        structure_io::read_structure(structureFile, crystal);

        nlohmann::json json_data;
        std::ifstream json_file("aspher.json");
        json_file >> json_data;
        vector<Vector3i> hkl;
        hkl_io::readHklIndices(hklFile, hkl);

        SfCalculator* sf_calc = SfCalculator::create(crystal, string("aspher.json"));
        vector<complex<double> > sf;
        sf_calc->calculateStructureFactors(crystal.atoms, hkl, sf);
        int i, n = hkl.size();
        ofstream out("sf");
        for (i = 0; i < n; i++)
            out << hkl[i] << " " << sf[i] << "\n";
        out.close();

        return 0;

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

