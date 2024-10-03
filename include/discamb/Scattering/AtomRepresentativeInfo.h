#pragma once

#include <string>

namespace discamb{

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    struct AtomRepresentativeInfo{
        int fragmentIdx = 0;
        int idxInSubsystem = 0;
        std::string atomLabel;
        std::string symmetryCode;
        std::string transformType;// identity, invert - invert of symmetry operation, symmetry - symmetry operation
        std::string transformSetting;
        bool isWeightFixed = true;
        double fixedWeightValue = 1.0;
        std::string weightRepresentativeLabel;
        int weightRepresentativeIndex = 0;
    };

    /** @}*/

}

