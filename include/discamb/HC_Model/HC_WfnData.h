
#ifndef _DISCAMB_HC_MODEL_HC_WFNDATA_H_
#define _DISCAMB_HC_MODEL_HC_WFNDATA_H_

//#include "Real.h"

#include <vector>
#include <string>
#include <complex>

namespace discamb {

struct HC_WfnBankEntry;

/** \ingroup HC_Model*/
namespace sto_atomic_wfn
{
    /**
    * \addtogroup HC_Model
    * @{
    */


    /**
    \ref hcd7 "normalization factor (Eq. (hc.d.7))"
    for radial part (\f$ \chi(r) \f$) of slater type orbital:
    
    \f[
    \chi(r) = \mathcal{N} (\alpha,k) r^{k} e^{-\alpha r} \\
    1 = \int \chi^2 r^2 dr
    \f]
    */
    inline double stoNormalizationFactor(int k, double alpha)
    {
        static const double factorial[] = { 1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600 };
        double result = pow(2 * alpha, int(k + 1))*sqrt(2 * alpha / factorial[2 * (k + 1)]);
        return result;
    }


    /**       
    Calculates \ref hcd12 "normalization factor (Eq. (hc.d.12))"
    for density normalized slater type radial function \f$ N r^n exp(-\zeta r) \f$ normalization factor for "density" normalization:
    \param zeta exponent in the radial function
    \param n power of r in radial function (max \p n = 9)
    */

    inline double stoDensityNormalizationFactor(int n, double zeta)
    {
        double invFactorial[] = { 1.0 / 1.0, 1.0 / 1.0, 1.0 / 2.0, 1.0 / 6.0, 1.0 / 24.0, 1.0 / 120.0, 1.0 / 720.0, 1.0 / 5040.0, 1.0 / 40320.0,
                                1.0 / 362880.0, 1.0 / 3628800.0, 1.0 / 39916800 };

        return pow(zeta, double(n + 3)) * invFactorial[ n + 2 ];
    }

    /*
    static const int TOTAL_DENSITY = 0;
    static const int CORE_DENSITY = 1;
    static const int VALENCE_DENSITY = 2;
    */
    /**
    Enumeration type used in function atomicStoWfnToSphericalDensity 
    */
    enum ElectronDensityType {
        TOTAL_DENSITY = 0,
        CORE_DENSITY = 1,
        VALENCE_DENSITY = 2
    };

    /** 
    Finds expression for spherically averaged atomic electron density (total, core or valence)
    \f[
    \rho(r) = \sum_v d_v r^{n_v} e^{-\beta_v r} 
    \f]
    (see also \ref hcd8 "eq. hc.d.8-9")
    \param wfn atomic wave function related data (including occupancy and assignment of orbitals to core/valence type)
    \param coefficients  the coefficients \f$ d_v \f$ of the linear combination
    \param exponents exponents \f$ \beta_v \f$ (in \f$ {\si{\angstrom}}^{-1}\f$)
    \param powerR powers \f$ n_v \f$ of \f$ r \f$ 
    \param densityType speciies type of density - total, core or valence
    */                                      
    
    void atomicStoWfnToSphericalDensity(
            const HC_WfnBankEntry &wfn,
            std::vector<double> &coefficients,
            std::vector<double> &exponents,
            std::vector<int> &powerR, 
            ElectronDensityType densityType);

    /**@}*/
}

/**
* \addtogroup HC_Model
* @{
*/

/**
Represents radial function of Slater type atomic orbital (see also \ref hcd5 "eqs. hc.d.5-7" )
being a combination of slater type functions:
\f[
\mathcal{R}(r) = \sum_p d_p r^{k_p} e^{-\alpha_p r}
\f]
*/

struct SlaterTypeAtomicOrbitalRdf
{
    /** Principal number of the orbital (e.g. 2 for 2S or 2P orbital)*/
    int principal_number = 1; 
    /** Azimuthal number of the orbital (e.g. 0 for s type orbitals, 1 for p type orbitals)*/
    int azimuthal_number = 0;
    /** Powers of r in the Slater type function (\f$ k_p \f$ in the above equation)*/
    std::vector<int> power_r;
    /** Exponents in the Slater type function (\f$ \alpha_p \f$ in the above equation )*/
    std::vector<double> exponents;
    /** Coefficients of the linear combination (\f$ d_p \f$ in the above equation ). The coefficients 
    already include normalization factors of Slater type function 
    (they correspond to products of \f$ c_p \f$ and \f$ \mathcal{N} (\alpha_p,k_p)\f$ 
    from eqs. \ref hcd5 "hc.d.5-6")*/
    std::vector<double> coefficients; 
};

/** Represents Hartree-Fock atomic wave function - orbitals and 
    occupancy including assignment to core and valence orbitals. */

struct HC_WfnBankEntry
{
    /** Label for the atom/ion */
    std::string wfnAtomTypeName;
    /** Atomic number */
    int atomic_number;
    /** Charge */
    int charge;
    /** Slater type orbitals (in fact their radial functions) */
    std::vector<SlaterTypeAtomicOrbitalRdf> orbitals;
    /** Orbital_occupancy[i] is an occupancy for orbital orbitals[i] */
    std::vector<int> orbital_occupancy;
    /** Indices of core orbitals - i.e. orbital[core_orbitals_indices[i]] is core orbital. */
    std::vector<int> core_orbitals_indices;
    /** Indices of valence orbitals - i.e. orbital[valence_orbitals_indices[i]] is valence orbital. */
    std::vector<int> valence_orbitals_indices;
    /** Constructs default empty wave function with label X, atomic number 0 and charge 0*/
    HC_WfnBankEntry(): wfnAtomTypeName("X"), atomic_number(0), charge(0){}
};

/** Represents atomic data in Hansen-Coppens model usually specific to atom/ion 
i.e. information on atomic wave function including assignment to core and valence orbitals,
exponents (\f$ \zeta \f$) and powers (\f$ n_l \f$) of r for radial functions
of deformation valence term (see eqs. \ref hcd11 "hc.d.11-12") and anomalous dispersion term.
*/
struct  HC_WfnType: HC_WfnBankEntry
{
    /** Assignment operator.*/
    HC_WfnType &operator=(const HC_WfnBankEntry &entry){HC_WfnBankEntry::operator=(entry);return *this;}
    /** Exponent (\f$ \zeta \f$) in the deformation valence radial function  (see eq. \ref hcd12 "hc.d.12")*/
    double deformation_valence_exponent;
    /** Exponent (\f$ n_l \f$) of r in the deformation valence radial function  (see eq. \ref hcd12 "hc.d.12")*/
    std::vector<int> deformation_valence_power;
    /** Anomalous dispersion term. */
    std::complex<double> anomalous_scattering;
    /** Default contructor.*/
    HC_WfnType(): deformation_valence_exponent(1.0){}
};

/**@}*/


}// namespace discamb

#endif /*_DISCAMB_HC_MODEL_HC_WFNDATA_H_*/


