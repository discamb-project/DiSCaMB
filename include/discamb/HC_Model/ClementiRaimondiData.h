#ifndef _DISCAMB_HC_MODEL_CLEMENTIRAIMONDIDATA_H_
#define _DISCAMB_HC_MODEL_CLEMENTIRAIMONDIDATA_H_

#include "HC_WfnData.h"
#include "HC_AtomTypeParameters.h"

namespace discamb {

    /** \defgroup HC_Model HC_Model
    \brief Hansen-Coppens model parameterization.
    */


    /** 
        \ingroup HC_Model
        \brief Exponents of single zeta slater type atomic orbitals for isolated atoms (from E. Clementi and D. L. Raimondi publication \cite Clementi_Raimondi_1963).

        
        */
    namespace clementi_raimondi {
        /** 
        \ingroup HC_Model 

        \brief Provides exponents of single zeta atomic orbitals.

        Fills \p exponents with single zeta exponents (in \f$ {\si{\angstrom}}^{-1}\f$) of 
        atomic orbitals for atoms defined by parameter \p atomicNumber. 
        Exponent corresponding to orbital with principal number \a n and azimuthal number \a l 
        is saved to \p exponents <tt>\[\a n - 1\]\[\a l\] </tt>. */
        void getExponents(int atomicNumber, std::vector<std::vector<double> > &exponents);
    }

} //namespace discamb

#endif /*_DISCAMB_HC_MODEL_CLEMENTIRAIMONDIDATA_H_*/
