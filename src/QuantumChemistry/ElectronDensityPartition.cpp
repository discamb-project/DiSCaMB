#include "discamb/QuantumChemistry/ElectronDensityPartition.h"
#include "discamb/QuantumChemistry/HirshfeldPartition.h"
#include "discamb/QuantumChemistry/MbisPartition.h"
#include "discamb/BasicUtilities/StringUtilities.h"

using namespace std;

namespace discamb {

    ElectronDensityPartition::~ElectronDensityPartition()
    {
    }

    ElectronDensityPartition* ElectronDensityPartition::create(
        const std::string& partitionName)
    {
        string _partitionName = string_utilities::toLower(partitionName);
        ElectronDensityPartition* partition;
        if (_partitionName == "hirshfeld")
            return new HirshfeldPartition();
        if (_partitionName == "mbis")
            return new MbisPartition();
        return nullptr;
    }
}
