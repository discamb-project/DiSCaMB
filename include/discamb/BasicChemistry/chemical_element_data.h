#ifndef _DISCAMB_BASICCHEMISTRY_CHEMICAL_ELEMENT_DATA_H_
#define _DISCAMB_BASICCHEMISTRY_CHEMICAL_ELEMENT_DATA_H_

#include "discamb/BasicChemistry/PeriodicTable.h"

#include <vector>

namespace discamb {

    namespace chemical_element_data {
        // in Angstroms
        double covalentRadius(int atomic_number);
        // in Angstroms
        double vdwRadius(int atomicNumber);
        double standardAtomicMass(int atomicNumber); 
        int groundStateSpinMultiplicity(int atomicNumber);
        
    }

}

#endif