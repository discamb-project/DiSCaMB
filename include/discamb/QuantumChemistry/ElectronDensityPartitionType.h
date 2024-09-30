#pragma once

#include <string>

namespace discamb{

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


enum class ElectronDensityPartitionType { Hirshfeld, IterativeHirshfeld, IterativeStockholder, MBIS, none};
ElectronDensityPartitionType stockholderPartitionTypeFromString(const std::string &type);
std::string  stockholderPartitionTypeAsString(ElectronDensityPartitionType type);

/**@}*/

}
