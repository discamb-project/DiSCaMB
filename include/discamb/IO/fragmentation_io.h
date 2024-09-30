#include "discamb/QuantumChemistry/fragmentation.h"

#include <vector>
#include <string>

namespace discamb {
    /**
    * \addtogroup IO IO
    * @{
    */

    namespace fragmentation_io {

        //void read(const std::string& fileName, std::vector<FragmentData>& fragments);
        void readPlainText(const std::string& fileName, std::vector<FragmentData>& fragments);
        void readJson(const std::string& fileName, std::vector<FragmentData>& fragments);
        void readHirsFrag(const std::string& fileName, std::vector<FragmentData>& fragments);

        void parseCappingAtomInfo(
            const std::string word1,
            const std::string word2,
            std::string& bondedAtom,
            std::string& bondedAtomSymmOp,
            std::string& directingAtom,
            std::string& directingAtomSymmOp);

    }
    /**@}*/
}

