#ifndef _DISCAMB_SCATTERING_SF_CALCDATATYPES_H_
#define _DISCAMB_SCATTERING_SF_CALCDATATYPES_H_

#include "Real.h"
#include "discamb/MathUtilities/Vector3.h"
#include <vector>
#include <complex>
#include <string>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


/**
  \brief Structure for storing derivatives of a target function with respect to atomic parameters.

  Stores data for parameter of one atom.
  */

struct TargetFunctionAtomicParamDerivatives
{
    Vector3<REAL> atomic_position_derivatives;
    std::vector<REAL> adp_derivatives;
    REAL occupancy_derivatives = 0;
};

/** 
\brief Structure for storing derivatives of a structure factor with respect to atomic parameters.

Stores data for parameters of muliple atoms. E.g. atomicPostionDerivatives[i] stores derivatives 
with respect to coordinates of (i+1)-th atom.
*/

struct SfDerivativesAtHkl
{
    std::vector<Vector3<std::complex<REAL> > > atomicPostionDerivatives;
    std::vector<std::vector<std::complex<REAL> > > adpDerivatives;
    std::vector<std::complex<REAL> > occupancyDerivatives;
};

struct ScatteringParameters
{
    std::string parameterisation; // e.g. IAM 
    std::vector<std::string> atomTypeSymbols;
    std::vector<std::complex<double> > atomTypeAnomalousDispersion; 
    std::vector<int> atomToTypeMap;

    static void groupAtomTypes(const std::vector<std::string> &atomType, 
                               std::vector<std::string> &types, 
                               std::vector<int> &atomToTypeMap);
};
/** @}*/
} // namespace discamb


#endif /*_DISCAMB_SCATTERING_SF_CALCDATATYPES_H_*/
