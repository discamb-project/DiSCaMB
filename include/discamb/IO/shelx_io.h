#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include <string>
#include <istream>

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


    namespace shelx_io {
        void save(const std::string &fName, const Crystal &crystal);
        void read(const std::string &fName, Crystal &crystal);
        void read(const std::string& fName, Crystal& crystal, std::map<std::string, std::string> &data);
        void read(std::istream& input, Crystal& crystal);
        void read(std::istream& input, Crystal& crystal, std::map<std::string, std::string>& data);

    }
    /**@}*/
}
