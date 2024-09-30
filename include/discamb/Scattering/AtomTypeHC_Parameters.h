#pragma once

#include <string>
#include <utility>
#include <vector>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    struct AtomTypeHC_Parameters
    {
        double p_val = 1.0;
        double p_val_sigma = 0.0;
        double kappa = 1.0;
        double kappa_sigma = 0.0;
        double kappa_prime = 1.0;
        double kappa_prime_sigma = 0.0;
        std::string symmetry;
        std::vector<double> p_lms;
        std::vector<double> p_lms_sigma;
        std::vector<std::pair<int, int> > p_lm_indices;
    };


    /**@}*/

}
