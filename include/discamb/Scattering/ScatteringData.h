
#ifndef _DISCAMB_SCATTERING_SCATTERINGDATA_H_
#define _DISCAMB_SCATTERING_SCATTERINGDATA_H_

#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    struct ReflectionsData 
    {
        std::vector<Vector3i> hkl;
        std::vector<double> intensity_calculated;
        std::vector<double> intensity_measured;
        std::vector<double> intensity_sigma;
        std::vector<std::string> scale_group_code;
        std::vector<char> observed_status;

        void clear()
        {
            hkl.clear();
            intensity_calculated.clear();
            intensity_measured.clear();
            intensity_sigma.clear();
            scale_group_code.clear();
            observed_status.clear();
        }
    };
    /** @}*/
}

#endif