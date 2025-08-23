#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/SfCalculator.h"

#include "json.hpp"
#include <iostream>

using namespace discamb;
using namespace std;

int main(int argc, char *argv[])
{
  
    if(argc!=2)
    {
        cout<< "expected structure file (cif/ins/res) as an argument\n";
        exit(0);
    }

    Crystal crystal;
    structure_io::read_structure(argv[1], crystal);
    vector<Vector3i> hkls{ {1,0,0}, { 1, 2, 0}, {-2, 3, 5}, {-1, -1, -1}, {3, -6, 2} };

    nlohmann::json setting = nlohmann::json::parse(R"(
    {
        "model": "taam",
        "electron scattering": true,
        "assignmnet info": "type_assignment.log",
        "bank path": "MATTS2021databank.txt"
    }
    )");

    auto form_factor_calculator = SfCalculator::create_shared_ptr(crystal, setting);

    vector<vector<complex<double> > > atomic_form_factors;
    int nAtoms = crystal.atoms.size();
    vector<bool> include_atom(nAtoms, true);
    
    form_factor_calculator->calculateFormFactors(hkls, atomic_form_factors, include_atom);

    for (int hklIdx = 0; hklIdx < hkls.size(); hklIdx++)
    {
        cout << hkls[hklIdx];
        for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
            cout << " " << atomic_form_factors[hklIdx][atomIdx];
        cout << endl;
    }

}


