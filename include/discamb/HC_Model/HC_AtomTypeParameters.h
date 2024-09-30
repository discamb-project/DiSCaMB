#ifndef _DISCAMB_HC_MODEL_HC_ATOMTYPEPARAMETERS_H_
#define _DISCAMB_HC_MODEL_HC_ATOMTYPEPARAMETERS_H_

#include "HC_WfnData.h"


#include <vector>

namespace discamb {

/** \ingroup HC_Model */

/** 
    \brief Container for the parameters of atomic electron density which are not related to the atomic wave function 
    and atom/ion electron configuration.

*/

struct HC_AtomTypeParameters
{
    /** 
     \brief Coefficients \f$ P_{l,m} \f$ of density normalized spherical harmonics in 
     deformation-valence term (see \ref hcd1 "eq. hc.d.1" ).

     \p p_lm [i][j] corresponds to \f$ P_{i, -i+j} \f$
     e.g. for l=2 (p_lm[2]) the parameters are stored in the following order:
     \f$ P_{2,-2} \f$ , \f$ P_{2,-1} \f$ , \f$ P_{2,0} \f$ , \f$ P_{2,1} \f$ , \f$ P_{2,2} \f$ .
     
    */

    std::vector<std::vector<double> > p_lm;

    /** \brief Parameter \f$ P_{val} \f$ of atomic electron density (see \ref hcd1 "eq. hc.d.1" ). */
    double p_val;

    /** \brief Parameter \f$ \kappa \f$ of atomic electron density (see \ref hcd1 "eq. hc.d.1" ). */
    double kappa_spherical_valence;

    /** \brief Parameter \f$ \kappa ' \f$ of atomic electron density (see \ref hcd1 "eq. hc.d.1" ). */
    double kappa_deformation_valence;

    /** \brief Constructs object with \f$ P_{val}=0 \f$ , \f$ \kappa=1 \f$ and \f$ \kappa ' = 1 \f$ .*/
    HC_AtomTypeParameters()
    {
        p_val = 0.0;
        kappa_spherical_valence = 1.0;
        kappa_deformation_valence = 1.0;
    }

};

}//namespace discamb 

#endif /*_DISCAMB_HC_MODEL_HC_ATOMTYPEPARAMETERS_H_*/


