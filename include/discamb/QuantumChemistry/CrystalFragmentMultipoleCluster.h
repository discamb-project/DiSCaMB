#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/StructuralProperties/atom_selection.h"

namespace discamb{
    /**
    * \addtogroup QuantumChemistry
    * @{
    */

    struct Settings {
        double multipoleClusterThreshold = 8.0;
        std::vector<std::pair<std::string, std::string> > predefinedMultipoleClusterAtoms;
        std::vector<AtomSubsetSelectionData> atomsToOmit;
        std::vector<std::pair<AtomSubsetSelectionData, double> > atomsWithCustomWeights;
    };

    class CrystalFragmentDistributedMultipoleCluster {
    public:
        CrystalFragmentDistributedMultipoleCluster();
        ~CrystalFragmentDistributedMultipoleCluster();

    };
    /**@}*/
}

