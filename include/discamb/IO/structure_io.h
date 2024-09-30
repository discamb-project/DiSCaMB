#include "discamb/CrystalStructure/Crystal.h"

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


    namespace structure_io {
        void read_structure(const std::string& fileName, Crystal& crystal);
        void write_structure(const std::string& fileName, const Crystal& crystal);
    }
    /**@}*/
}
