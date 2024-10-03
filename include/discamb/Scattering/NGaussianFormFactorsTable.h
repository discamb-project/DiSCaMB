#ifndef _DISCAMB_SCATTERING_NGAUSSIANFORMFACTORSTABLE_H_
#define _DISCAMB_SCATTERING_NGAUSSIANFORMFACTORSTABLE_H_
#include "NGaussianFormFactor.h"

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


namespace n_gaussian_form_factors_table
{
    /** Possible values for table: "Waasmeier-Kirfel", "IT92", "electron-IT", "electron-cctbx"*/
    NGaussianFormFactor getFormFactor(const std::string &label,const std::string &table=std::string("Waasmeier-Kirfel"));
    bool hasFormFactor(const std::string &label, const std::string &table = std::string("Waasmeier-Kirfel"));
}

/** @}*/
}


#endif /*_DISCAMB_SCATTERING_NGAUSSIANFORMFACTORSTABLE_H_*/
