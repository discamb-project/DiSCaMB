

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/Timer.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/agreement_factors.h"
#include "discamb/Scattering/SfCalculator.h"
#include <fstream>

using namespace std;
using namespace discamb;



int main(int argc, char* argv[])
{
    try {

        vector<string> arguments, options;
        parse_cmd::get_args_and_options(argc, argv, arguments, options);
        
        string message = "invalid input, expected arguments:\n"
            "  (1) structure file (cif or shelx ins/res)\n"
            "  (2) shelx hkl indices file\n"
            "  (3) json config file\n"
            "and optionally:\n"
            "  -s - for version with symmetry only\n"
            "  -d - do not no calculations for derivatives\n";

        if (arguments.size() != 3)
            on_error::throwException(message, __FILE__, __LINE__);
        
        bool onlyWithSymmetry = parse_cmd::hasOption(options, "-s");
        bool sf_and_sfAndDf = !parse_cmd::hasOption(options, "-d");

        //taam_parallel(argv[1], argv[2], onlyNew, sf_and_sfAndDf);
        string structureFile = arguments[0];
        string hklFile = arguments[1];
        string jsonFileName = arguments[2];


        Crystal crystal;
        structure_io::read_structure(structureFile, crystal);
        vector<vector<complex<double> > > formFactors;
        vector<Vector3i> hkl_all, hkl;
        hkl_io::readHklIndices(hklFile, hkl);

        int nHkl = hkl.size();
        int nAtoms = crystal.atoms.size();

        nlohmann::json json_data;
        ifstream jsonFileStream(jsonFileName);
        if (jsonFileStream.good())
            jsonFileStream >> json_data;
        jsonFileStream.close();


        json_data["def-val symmetry"] = false;
        shared_ptr<SfCalculator> sfCalculator = shared_ptr<SfCalculator>(SfCalculator::create(crystal, json_data));
        json_data["def-val symmetry"] = true;
        shared_ptr<SfCalculator> sfCalculatorSymm = shared_ptr<SfCalculator>(SfCalculator::create(crystal, json_data));



        vector<bool> countAtom(nAtoms, true);
        vector<complex<double> > sf_no_symm, sf_symm;
        vector<TargetFunctionAtomicParamDerivatives> dT_dp_no_symm, dT_dp_symm;
        vector<complex<double> > dT_dF(nHkl, complex<double>(1.0, 0.3));
        // x y z U occ
        SfDerivativesAtHkl derivatives_symm, derivatives_no_symm;
        WallClockTimer timer;

        // structure factors call 'calculateStructureFactors'

        if (sf_and_sfAndDf)
        {
            if (!onlyWithSymmetry)
            {
                timer.start();
                sfCalculator->calculateStructureFactors(crystal.atoms, hkl, sf_no_symm, countAtom);
                cout << "time = " << timer.stop() << " ms\n";
            }
            timer.start();
            sfCalculatorSymm->calculateStructureFactors(crystal.atoms, hkl, sf_symm, countAtom);
            cout << "time = " << timer.stop() << " ms\n";
            ofstream out("out");
            for (int i = 0; i < sf_symm.size(); i++)
                out << sf_symm[i] << " " << sf_no_symm[i] << "\n";
            cout << "agreement factor " << agreement_factors::value(sf_symm, sf_no_symm) << "\n";
        }

        // structure factors and derivatives:  'calculateStructureFactorsAndDerivatives'
        if (!onlyWithSymmetry)
        {
            timer.start();
            sfCalculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_no_symm, dT_dp_no_symm, dT_dF, countAtom);
            cout << "time = " << timer.stop() << " ms\n";
        }

        timer.start();

        sfCalculatorSymm->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_symm, dT_dp_symm, dT_dF, countAtom);
        cout << "time = " << timer.stop() << " ms\n";

        if (!onlyWithSymmetry)
            cout << "agreement factor " << agreement_factors::value(sf_no_symm, sf_symm) << "\n";




        bool size_match;
        double d_xyz_agreement_factor, d_adp_agreement_factor, d_occ_agreement_factor, d_fpfdp_agreement_factor;
        agreement_factors::for_derivatives(
            dT_dp_no_symm,
            dT_dp_symm,
            size_match,
            d_xyz_agreement_factor, d_adp_agreement_factor, d_occ_agreement_factor, d_fpfdp_agreement_factor);

        cout << "xyz : " << d_xyz_agreement_factor << "\n"
            << "adp : " << d_adp_agreement_factor << "\n"
            << "occ : " << d_occ_agreement_factor << "\n";


    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

