#ifndef _DISCAMB_HC_MODEL_DEFORMATIONVALENCEPARAMETERS_H_
#define _DISCAMB_HC_MODEL_DEFORMATIONVALENCEPARAMETERS_H_

#include <string>
#include <vector>
#include <map>
#include <utility>

namespace discamb {

/** 
\ingroup HC_Model 

\brief A container for deformation-valence parameters. 

Involves the parameters (\f$ n_l \f$ and \f$ \zeta \f$) of deformation-valence term
(\f$ \rho_{def-val}(\mathbf{r}) \f$) in Hansen-Coppens model:
\f[
\rho_{def-val}(\mathbf{r}) = \sum_l r^{n_l} e^{ - \zeta \kappa^{'} r} \sum_m P_{lm} Y_{lm} (\theta , \phi)
\f]
*/

class DeformationValenceParameters
{
public:
    
    /** \brief Indicates one of two built in parameterizations of the deformation valence parameters.*/
    enum class ParametersType {
        /**
        Default parameterization based on algorithms described in \ref algorithms_for_parameterization "Algorithms for parameters assignment section". See also 
        \ref xd "comparison with XD2006" parameters.
        */
        STANDARD, 
        /**
        Parameters used in UBDB (University Bufallo Data Bank) - the same as for STANDARD except for
        P (powers 4, 6, 6, 6, 6) and S (powers 1, 2, 4, 6, 8).
        */
        UBDB,
        /**
        Setting for Su, Coppens & Macchi wavefunction parameters for atoms and ions up to Cs+.
        */
        SCM};
    DeformationValenceParameters();
    /**
    \brief Constructs one of two buil in parameter sets.
    
    The parameter \p type defines which set is used.
    \sa ParametersType
    */
    DeformationValenceParameters(ParametersType type);
    ~DeformationValenceParameters();

    /**
    \brief Set deformation valence parameters

    \param types - labels for atom types
    \param zeta - exponents (\f$\zeta\f$) for the corresponding types (units \f$ 1/ {\si{\angstrom}} \f$) 
    \param powers_r - powers (\f$n_l\f$) of r in the radial functions
    */
    void setParameters(const std::vector<std::string> &types,const std::vector<double> &zeta,
                       const std::vector<std::vector<int> > &powers_r);
    /** \brief Sets parameters to of two buil in parameter sets 

    \sa ParametersType
    */
    void set(ParametersType type);

    /**
    \brief Set deformation valence parameters for specified type.

    \param type - labels for the atom/ion type
    \param zeta - exponent (\f$\zeta\f$) for the corresponding type (units \f$ 1/ {\si{\angstrom}} \f$)
    \param powers_r - powers (\f$n_l\f$) of r in the radial function
    */
    
    void setParameter(const std::string &type, double zeta,const std::vector<int> &powers_r);

    /** 
    \brief Get deformation valence parameters for specified type.
    
    Involves possibility for specification of 
    non default valence orbitals, electron configuration and single zeta exponents:
    \param type - label for the atom/ion type
    \param zeta - exponent (\f$\zeta\f$), units \f$ 1/ {\si{\angstrom}} \f$
    \param powers_r - powers (\f$n_l\f$) of r in the radial function
    \param valenceOrbitals - valenceOrbitals[i].first is principal number of the orbital while valenceOrbitals[i].second is its azimuthal number
    \param configuration - electron configuration configuration[i][j] corresponds to occupancy if 
                           (j+1)-th subshell of (i+1)-th shell, if not given then default is used 
    \param exponents - exponents[i-1][j-1] corresponds to exponent of single zeta atomic orbital for j-th subshell of i-th shell 
                       (expected in units \f$ 1/ {\si{\angstrom}} \f$), if not given then default is used 
    */
    bool getParameters(const std::string &type, double &zeta, std::vector<int> &powers_r,
                       const std::vector<std::pair<int, int> > &valenceOrbitals,
                       const std::vector<std::vector<int> > &configuration = std::vector<std::vector<int> >(),
                       const std::vector<std::vector<double> > &exponents = std::vector<std::vector<double> >()) const;

    /**
    \brief Get deformation valence parameters for specified type.

    \param type - label for the atom/ion type
    \param zeta - exponent (\f$\zeta\f$), units \f$ 1/ {\si{\angstrom}} \f$
    \param powers_r - powers (\f$n_l\f$) of r in the radial function
    */

    bool getParameters(const std::string &type, double &zeta, std::vector<int> &powers_r) const;

    /**
    \brief Get deformation valence parameters for all types in the container. 

    \param types - labels for atom/ion types
    \param zeta - exponents (\f$\zeta\f$), units \f$ 1/ {\si{\angstrom}} \f$
    \param powers_r - powers (\f$n_l\f$) of r in the radial function
    */

    void getParameters(std::vector<std::string> &types, std::vector<double> &zeta,
                       std::vector<std::vector<int> > &powers_r) const;
private:
    // mParameters[label].first - zeta, mParameters[label].second - powers of r 
    std::map<std::string, std::pair<double, std::vector<int> > > mParameters;
    void setUbdbParameterization();
    void setStandardParameterization();
    void setSCMParameterization();
};

} // namespace discamb

#endif /*_DISCAMB_HC_MODEL_DEFORMATIONVALENCEPARAMETERS_H_*/

