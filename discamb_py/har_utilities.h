#include "discamb/Scattering/HirshfeldAtomModelSettings.h"
#include "CrystalStructure.h"

namespace har_utilities {

    struct Representative {

        int substructureIdx;
        int atomIdx;
        double weight;
    };

    std::vector<Representative> find_default_representatives(
        const CrystalStructure & crystal_structure, 
        const std::vector< std::vector<std::vector<int> > >& fragments);

}