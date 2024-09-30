#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicUtilities/utilities.h"
#include "discamb/QuantumChemistry/WaveFunctionFileFormat.h"
#include "json.hpp"

#include <vector>
#include <string>


namespace discamb{

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    struct QmSystem {
        std::vector<Vector3d> positions;
        std::vector<int> atomicNumbers;
        int charge;
        int spin_multilicity;
        std::vector<double> pointChargeValue;
        std::vector<Vector3d> pointChargePosition;
    };

    struct QmSettings {
        void set(const nlohmann::json& data);
        std::string qmMethod;
        std::string relativisticMethod;
        std::string basisSet;
        std::map<int, std::string> atomicNumber2BasisSetMap;
        std::string inputTemplate;
        bool tryToReadGuess;
    };

    struct WaveFunctionCalculationData {
        std::string jobName;
        std::string wfnFileName;
        HardwareResources hardware;
        QmSettings qmSettings;
        QmSystem qmSystem;
        std::map<int, std::string> atomIdx2BasisSetMap;
        nlohmann::json qmProgramSpecificData;
        WaveFunctionFileFormat wfnFormat = WaveFunctionFileFormat::unspecified;
    };
    /**@}*/
}

