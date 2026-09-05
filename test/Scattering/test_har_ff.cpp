#include "discamb/Scattering/StockholderAtomSfCalculator.h"
#include "discamb/CrystalStructure/structure_library.h"
#include <string>
#include <iostream>

using namespace std;
using namespace discamb;

int main(int argc, char *argv[])
{
    Crystal crystal;
    structure_library::getStructure("urea", crystal);

    HirshfeldAtomModelSettings settings;
    settings.crystalFragments.resize(1);
    //settings.crystalFragments[0].atoms = 
    //StockholderAtomSfCalculator calculator(;
    
	return 0;
}
