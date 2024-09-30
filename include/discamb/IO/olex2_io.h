#pragma once

#include "discamb/CrystalStructure/CrystalVarianceCovarianceMatrix.h"

#include <string>
#include <vector>

namespace discamb
{
    /**
    * \addtogroup IO IO
    * @{
    */

    namespace olex2_io {

        void read_vcov(const std::string& fName, std::vector<std::vector<double> >& data, std::vector<std::string>& dataLabels);
        void read_vcov_npy(const std::string& fName, std::vector<std::vector<double> >& data, std::vector<std::string>& dataLabels);
        void read_vcov(const std::string& fName, CrystalVarianceCovarianceMatrix& vcov);
        void read_vcov_npy(const std::string& fName, CrystalVarianceCovarianceMatrix& vcov);
    }
    /**@}*/
}

