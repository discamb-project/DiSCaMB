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
#include "AnyScattererStructureFactorCalculator2.h"
#include <fstream>
#include <memory>

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
    shared_ptr<SfCalculator> sf_calculator = shared_ptr<SfCalculator>(SfCalculator::create(crystal, json_data));
    AnyScattererStructureFactorCalculator2 sf_calculator2(crystal);
    sf_calculator2.setSfCalculator(sf_calculator);

    cout << "sf calculator created in " << timer.stop() << " ms\n";
    
    timer.start();
    cout << "calculating sf for " << crystal.atoms.size() << " atoms and " << hkl.size() << " reflections\n";
    sf_calculator2.calculateStructureFactors(hkl, sf);
    cout << "sf calculated in " << timer.stop() << " ms\n";
    ofstream out("cals_sf.txt");
    for (int i = 0; i < hkl.size(); i++)
        out << hkl[i] << " " << sf[i] << "\n";
    out.close();
    /*
    vector<TargetFunctionAtomicParamDerivatives> dTarget_dParam;
    vector<complex<double> > dTarget_dF;
    dTarget_dF.resize(hkl.size(),0.01);
    cout << "calculating sf and derivatives\n";
    timer.start();
    sf_calculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf, dTarget_dParam, dTarget_dF, countAtom);
    cout << "sf calculated in " << timer.stop() << " ms\n";
    */
}

int main(int argc, char* argv[])
{
    try { 
        calc_sf(argv[1], argv[2]);
        return 0;


    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

