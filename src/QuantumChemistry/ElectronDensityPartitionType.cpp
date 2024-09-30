#include "discamb/QuantumChemistry/ElectronDensityPartitionType.h"
#include "discamb/BasicUtilities/OnError.h"
using namespace std;

/*
namespace discamb{
enum class StockholderPartitionType { Hirshfeld, IterativeHirshfeld, IterativeStockholder};
StockholderPartitionType stockholderPartitionTypeFromString(const std::string &type);
}
*/

namespace discamb{

//enum class StockholderPartitionType { Hirshfeld, IterativeHirshfeld, IterativeStockholder};

    ElectronDensityPartitionType stockholderPartitionTypeFromString(
        const std::string &_type)
    {
        string type;
        for (char c : _type)
            type += tolower(c);

        if (type == string("ih") || type == string("iterativehirshfeld") || type == string("iterative hirshfeld") || type == string("iterative_hirshfeld"))
            return ElectronDensityPartitionType::IterativeHirshfeld;
        if (type == string("h") || type == string("hirshfeld"))
            return ElectronDensityPartitionType::Hirshfeld;
        if (type == string("is") || type == string("iterativestockholder") || type == string("iterative stockholder") || type == string("iterative_stockholder"))
            return ElectronDensityPartitionType::IterativeStockholder;
        if (type == string("mbis") || type == string("minimal basis iterative stockholder") || type == string("minimal_basis_iterative_stockholder"))
            return ElectronDensityPartitionType::MBIS;
        on_error::throwException(string("unknown type of electron density partition:'") + _type + string("'"), __FILE__, __LINE__);
        return ElectronDensityPartitionType::Hirshfeld;
    }

    std::string  stockholderPartitionTypeAsString(
        ElectronDensityPartitionType  type)
    {
        if (type == ElectronDensityPartitionType::IterativeHirshfeld)
            return string("iterative Hirshfeld");
        if (type == ElectronDensityPartitionType::Hirshfeld)
            return string("Hirshfeld");
        if (type == ElectronDensityPartitionType::IterativeStockholder)
            return string("iterative stockholder");
        if (type == ElectronDensityPartitionType::IterativeStockholder)
            return string("minimal basis iterative stockholder");

        return string("unknown");
    }
}
