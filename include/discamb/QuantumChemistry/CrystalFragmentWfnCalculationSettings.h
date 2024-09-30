#pragma once

#include "discamb/QuantumChemistry/WaveFunctionCalculationData.h"
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"
#include "discamb/QuantumChemistry/fragmentation.h"

namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */

    struct CrystalWfnCalcWithFragmentsSettings {
        // qm program
        WaveFunctionDataGeneratorRunner::Program qmProgram;
        nlohmann::json qmProgramSpecificData;
        //hardware
        int nCpu;
        int memoryMB;
        // method
        std::string qmMethod;
        std::string relativisticMethod;
        // basis set
        std::string basisSet;
        std::map<int, std::string> atomicNumber2BasisSetMap;
        // input template
        std::string inputTemplate;

    };

    struct QmCrystalFagment {
        int charge;
        int spin_multilicity;
        std::vector<std::pair<std::string, std::string> > atoms;
        std::vector<CappingHydrogen> cappingHydrogens;
    };

    struct CrystalFragmentWfnCalculationSettings{
        // qm fragment
        QmCrystalFagment fragment;
        // names
        std::string jobName;
        std::string wfnFileName;
        // qm program
        WaveFunctionDataGeneratorRunner::Program qmProgram;
        nlohmann::json qmProgramSpecificData;
        //hardware
        int nCpu;
        int memoryMB;
        // method
        std::string qmMethod;
        std::string relativisticMethod;
        // basis set
        std::string basisSet;
        std::map<int, std::string> atomicNumber2BasisSetMap;
        // charge & spin 
        // point charges
        std::vector<double> pointChargeValue;
        std::vector<Vector3d> pointChargePosition;
        
        std::string inputTemplate;

    };

    /**@}*/
}

