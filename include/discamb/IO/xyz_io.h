#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicChemistry/ChemicalElement.h"
#include "discamb/BasicChemistry/MoleculeData.h"

#include <vector>
#include <string>

namespace discamb{

    /**
    * \addtogroup IO IO
    * @{
    */


namespace xyz_io{
    void readXyz(const std::string &fileName, std::vector<std::string> &symbol, std::vector<Vector3d> &position);
    void readXyz(const std::string &fileName, std::vector<int> &atomicNumber, std::vector<Vector3d> &position);
    void readXyz(const std::string& fileName, std::vector<ChemicalElement>& element, std::vector<Vector3d>& position);
    void readXyz(const std::string& fileName, MoleculeData& data);
    void writeXyz(const std::string &fileName,const  std::vector<std::string> &symbol, const std::vector<Vector3d> &position);
    void writeXyz(const std::string &fileName, const std::vector<int> &atomicNumber, const std::vector<Vector3d> &position);
    void writeXyz(const std::string& fileName, const std::vector<ChemicalElement>& element, const std::vector<Vector3d>& position);
    void writeXyz(const std::string& fileName, const MoleculeData& data);
}
/**@}*/
}