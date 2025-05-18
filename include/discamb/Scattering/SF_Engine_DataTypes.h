#ifndef _DISCAMB_SCATTERING_SF_ENGINE_DATATYPES_H_
#define _DISCAMB_SCATTERING_SF_ENGINE_DATATYPES_H_

#include "discamb/MathUtilities/Matrix3.h"
#include "Real.h"
#include <vector>
#include <complex>
#include <string>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    /** 
    \brief structures used in structure factor calculations intended mainly for internal use
    */

    

    namespace sf_engine_data_types {




        /**
          \brief structure intended mainly for internal use

          Structure related to description of the electron density based on wave-function for single atom.
          Characteristic to chemical element. 
          */


        struct HC_WfnParam
        {
            std::string label;

            /**
              parameters of slater functions describing core electron density
              electron_density(r) = sum_{k} core_coeff[k] * r^core_pow[k] * exp( - core_exp[k]*r )
              */
            std::vector<REAL> core_coeff;
            std::vector<REAL> core_exp;
            std::vector<int> core_pow;

            /**
              parameters of slater functions describing spherical valence electron density
              (in order to describe fully this term also parameters characteristic to atom type
              are necessary - stored in structure HC_TypeParam : HC_TypeParam::kappa_spherical
              and HC_TypeParam::p_val)
              */

            std::vector<REAL> valence_coeff;
            std::vector<REAL> valence_exp;
            std::vector<int> valence_pow;

            /**
              one of factor in exponent of a exponential function in the deformation valence density term
              (the other one is HC_TypeParam::kappa_def_valence)
              */

            REAL def_valence_exp = 1.0;

            /**
              powers or r in radial part of deformation valence density term indexed with l
              */

            std::vector<int> def_valence_pow;

            /** anomalous_scattering */

            std::complex<REAL> anomalous_scattering;
        };

        /**

          \brief structure intended mainly for internal use

          Structure related to data describing atom type
          (in case of protein structure probably there can be let say 100 different atom types)
          */

        struct HC_TypeParam
        {
            /** coefficients in the multipolar expansion part of the atomic density (deformation valence density) */
            std::vector<std::vector<REAL> > p_lm;

            /** prefactor for spherical valence part of the electron density */
            REAL p_val;

            /** scale factor for deformation valence part of the electron density */
            REAL kappa_def_valence;

            /** scale factor for spherical valence part of the electron density */
            REAL kappa_spherical;
        };


        /**

          \brief structure intended mainly for internal use

          See classes discamb::SpaceGroupOperation and discamb::SpaceGroup for general structures for handling symmetry in DiSCaMB.

          Crystallographic symmetry operation - it transforms vector r in the following way:
          r' = rotation*r + translation
          rotation is an orthogonal matrix (not necessarily corresponding to rotation)
          */

        struct SymmetryOperation
        {
            Matrix3<REAL>  rotation;
            Vector3<REAL> translation;
        };
    } // namespace sf_engine_data_types

/** @}*/

} // namespace discamb

#endif /*_DISCAMB_SCATTERING_SF_ENGINE_DATATYPES_H_*/
