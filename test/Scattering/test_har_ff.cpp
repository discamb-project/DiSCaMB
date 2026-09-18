#include "discamb/Scattering/StockholderAtomSfCalculator.h"
#include "discamb/CrystalStructure/structure_library.h"
#include "discamb/QuantumChemistry/OrcaRunner.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include <string>
#include <iostream>

using namespace std;
using namespace discamb;

const string aspher_json_urea = R"({
    "model": "HAR",
    "electron scattering": false,
    "qm program": "orca",
    "qm method": "PBE",
    "basis set": "def2-SVP",
    "n cores": 2,
    "memory": "2GB",
    "partition settings": {
        "atoms file": "atomic_densities.txt"
    },
    "qm structure": [
        {
            "label": "subsystem_1",
            "charge": 0,
            "spin multiplicity": 1,
            "atoms": " C O N N,-X,-Y+1,Z H1 H2 H1,-X,-Y+1,Z H2,-X,-Y+1,Z"
        }
    ]
}
)";

//"atoms file": "C:\\Program Files (x86)\\dGui\\data\\atomic_densities.txt"

const vector<Vector3i> hkl{ {-9, -10, -3}, {-5, 5, 1}, {0, 0, 0}, {1, -1, 2}, {0, 0, 2} };

const vector < vector<complex<double> > > ff_reference{
    {{0.850048, 0.000566}, { 1.161227, 0.002038}, {1.037299, 0.000985}, { 0.002543, 0.001335}, { 0.003893, -0.000245}},
    {{1.513708,-0.000661}, {1.849211,0.002440}, {1.588482,-0.000555}, {0.039837,-0.000820}, {0.041385,0.001283} },
    {{5.849983, 0.000000}, { 8.359612, 0.000000}, { 7.147759, 0.000000}, { 0.865058, 0.000000}, { 0.881359, 0.000000}},
    {{2.985411,-0.041379}, {4.985679,0.091183}, {4.024863,0.033339}, {0.374709,-0.044979}, {0.376942,0.070183} },
    {{3.326942,-0.072419}, {5.597934,0.137092}, {4.597560,0.035254}, {0.455637,-0.064374}, {0.452510,0.100616} }
};

int main(int argc, char *argv[])
{
    file_system_utilities::NewFilesRemover newFilesRemover;
    bool passed = true;

    try {

        Crystal crystal;
        structure_library::getStructure("urea", crystal);
        nlohmann::json data = nlohmann::json::parse(aspher_json_urea);
        string orca_folder;
        if (!discamb::OrcaRunner::findOrcaFolder(orca_folder))
            return 1;
        data["qm folder"] = orca_folder;

        //settings.crystalFragments[0].atoms = 
        StockholderAtomSfCalculator calculator(crystal, data);
        vector<vector<complex<double> > > formFactors;
        int nAtoms = crystal.atoms.size();
        calculator.calculateFormFactors(hkl, formFactors, vector<bool>(nAtoms, true));

        for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
        {
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                if (abs(formFactors[hklIdx][atomIdx].imag() - ff_reference[hklIdx][atomIdx].imag()) > 1e-06)
                    passed = false;
                if (abs(formFactors[hklIdx][atomIdx].real() - ff_reference[hklIdx][atomIdx].real()) > 1e-06)
                    passed = false;
            }

        }
    }
    catch (...)
    {
        cout << "exeption thrown when executing test_har_ff " << endl;
        newFilesRemover.cleanNew();
        return 1;
    }
    
    newFilesRemover.cleanNew();

    cout << "test_har_ff ";
    if (!passed)
        cout << " not ";
    cout << "passed\n";

    if (passed)
        return 0;

	return 1;
}

