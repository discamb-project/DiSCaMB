#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/taam_utilities.h"
#include "discamb/Scattering/TscFileBasedSfCalculator.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include <fstream>

using namespace std;
using namespace discamb;


void calc_sf(
    const string& structure_file,
    const string& hkl_file)
{
    Crystal crystal;
    WallClockTimer timer;
    timer.start();
    structure_io::read_structure(structure_file, crystal);
    cout << "structure read in " << timer.stop() << " ms\n";


    vector<Vector3i> hkl;

    hkl_io::readHklIndices(hkl_file.c_str(), hkl);

    vector<complex<double> > sf;
    nlohmann::json json_data;
    ifstream jsonFileStream("aspher.json");

    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    timer.start();
    auto sf_calculator = SfCalculator::create(crystal, json_data);
    cout << "sf calculator created in " << timer.stop() << " ms\n";
    vector<bool> countAtom(crystal.atoms.size(), true);
    timer.start();
    cout << "calculating sf for " << crystal.atoms.size() << " atoms and " << hkl.size() << " reflections\n";
    sf_calculator->calculateStructureFactors(crystal.atoms, hkl, sf, countAtom);
    cout << "sf calculated in " << timer.stop() << " ms\n";
    ofstream out("cals_sf.txt");
    for (int i = 0; i < hkl.size(); i++)
        out << hkl[i] << " " << sf[i] << "\n";
    out.close();
    vector<TargetFunctionAtomicParamDerivatives> dTarget_dParam;
    vector<complex<double> > dTarget_dF;
    dTarget_dF.resize(hkl.size(),0.01);
    cout << "calculating sf and derivatives\n";
    timer.start();
    sf_calculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf, dTarget_dParam, dTarget_dF, countAtom);
    cout << "sf calculated in " << timer.stop() << " ms\n";
}

int main(int argc, char* argv[])
{
    try { 
        calc_sf(argv[1], argv[2]);
        return 0;

        tsc_differences(argv[1], argv[2]);
        return 0;
        
        string structFile = "";
        if (argc == 2)
            structFile = argv[1];
        test_connectivity(structFile);
        return 0;

        test_int_floor();
        return 0;

        check_args(argc, argv, 2, { "structure file", "bank file path" });
        testTypeAssignemntLog(argv[1], argv[2]);
        return 0;

        check_args(argc, argv, 2, { "structure file", "bank file path"});
        Crystal crystal;
        structure_io::read_structure(argv[1], crystal);
        TAAM_types(crystal, argv[2]);
        return 0;

        check_args(argc, argv, 1, { "structure file" });
        nGaussian(argv[1]);
        return 0;
        compare_sf_taam_qm(argc, argv);
        return 0;
        test_taam_disorder(argc, argv);

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

