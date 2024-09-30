#ifndef _DISCAMB_SCATTERING_MULTIPOLE_SCATTERING_H_
#define _DISCAMB_SCATTERING_MULTIPOLE_SCATTERING_H_

#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <complex>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    /** \brief Tools for calculation of scattering from potential described with multipole expansion.

    Includes calculation of 'angular part' of scattering factor and calculation of scattering factor for 
    deformation-valence term in Hansen-Coppens model. 
    
    Note: these functions are not used in DiSCaMB for speed optimized calculations for Hansen-Coppens multipole model.
    
    */

    namespace multipole_scattering {

    
    /**
    * \addtogroup Scattering Scattering
    * @{
    */
    
    
        /** \brief Calculate scattering factor for deformation valence term in Hansen-Coppens model.
        
        Version with \f$ \kappa ' \f$ as one of parameters.

        \param p_lms coefficients \f$ P_{l,m} \f$ of density normalized spherical harmonics in
               deformation-valence term (see \ref hcd1 "eq. hc.d.1" ). p_lm [i][j] corresponds to \f$ P_{i, -i+j} \f$
               e.g. for l=2 (p_lm[2]) the parameters are stored in the following order:
               \f$ P_{2,-2} \f$ , \f$ P_{2,-1} \f$ , \f$ P_{2,0} \f$ , \f$ P_{2,1} \f$ , \f$ P_{2,2} \f$ .
        \param kappa parameter \f$ \kappa ' \f$ of atomic electron density (see \ref hcd1 "eq. hc.d.1" ). 
        \param zeta exponent (\f$ \zeta \f$) in the deformation valence radial function  (see eq. \ref hcd12 "hc.d.12")
        \param powerR exponents (\f$ n_l \f$) of r in the deformation valence radial function  (see eq. \ref hcd12 "hc.d.12")
        \param h vector in resiprocal space in Cartesian coordinates (in \f$ {\si{\angstrom}}^{-1}\f$)
        */
        std::complex<double> hcDeformationValence(
            const std::vector<std::vector<double> > &p_lms, 
            double kappa,
            double zeta,
            const std::vector<int> &powerR,
            const Vector3d &h);

        /** \brief Calculate scattering factor for deformation valence term in Hansen-Coppens model.

        Version without \f$ \kappa ' \f$.

        \param p_lms coefficients \f$ P_{l,m} \f$ of density normalized spherical harmonics in
        deformation-valence term (see \ref hcd1 "eq. hc.d.1" ). p_lm [i][j] corresponds to \f$ P_{i, -i+j} \f$
        e.g. for l=2 (p_lm[2]) the parameters are stored in the following order:
        \f$ P_{2,-2} \f$ , \f$ P_{2,-1} \f$ , \f$ P_{2,0} \f$ , \f$ P_{2,1} \f$ , \f$ P_{2,2} \f$ .
        \param zeta exponent (\f$ \zeta \f$) in the deformation valence radial function  (see eq. \ref hcd12 "hc.d.12")
        \param powerR exponents (\f$ n_l \f$) of r in the deformation valence radial function  (see eq. \ref hcd12 "hc.d.12")
        \param h vector in resiprocal space in Cartesian coordinates (in \f$ {\si{\angstrom}}^{-1}\f$)
        */

        std::complex<double> hcDeformationValence(
            const std::vector<std::vector<double> > &p_lms,
            double zeta,
            const std::vector<int> &powerR,
            const Vector3d &h);



        std::complex<double> multipoleScattering(
            int l, 
            int m,
            const Vector3d &h,
            bool densityNormalized);

        std::complex<double> multipoleScattering(
            int l, 
            const std::vector<double> &multipliers,
            const Vector3d &h,
            bool densityNormalized);

        std::complex<double> multipoleScattering(
            const std::vector<std::vector<double> > &multipliers,
            const Vector3d &h,
            bool densityNormalized);

    }
    /** @}*/
}

#endif /* _DISCAMB_SCATTERING_MULTIPOLE_SCATTERING_H_ */

