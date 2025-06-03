#ifndef _DISCAMB_SCATTERING_HansenCoppens_SF_Engine4_H_
#define _DISCAMB_SCATTERING_HansenCoppens_SF_Engine4_H_


#include "Real.h"
#include "SF_Engine_DataTypes.h"
#include "SF_CalcDataTypes.h"
#include "discamb/MathUtilities/real_spherical_harmonics.h"
#include "NGaussianFormFactor.h"

#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include "NGaussianFormFactorsTable.h"

#include <cassert>
#include <cmath>
#include <fstream>
#include <complex>

namespace discamb{

    /**
    * \addtogroup Scattering Scattering
    * @{
    */

/** \brief Engine for structure factor calculations with CPU intended mainly to be used internally in DiSCaMB
           (while HansenCoppensStructureFactorCalculator is intended to be more convenient interface for 
           structure factor calculations). 

*/

class HansenCoppens_SF_Engine4
{
public:
    HansenCoppens_SF_Engine4();
    ~HansenCoppens_SF_Engine4();


/**
  \brief Calculates Hansen-Coppens multipole model structure factors and derivatives 
  of a target function with respect to atomic parameters.

  \param wfn_parameters sets of wave-function atomic parameters 
  \param type_parameters sets of parameters related to atom types
  \param atom_to_wfn_map atom_to_wfn_map[i] is an index of wfn_parameters set corresponding to i-th atom
  \param atom_to_type_map atom_to_type_map[i] is an index of type_parameters set corresponding to i-th atom
  \param atomicPositions positions of atoms in asymmetric unit (in Cartesian coordinate system, Angstroms)
  \param atomic_displacement_parameters parameters describing atomic displacement (in Cartesian coordinate system)
  \param atomic_occupancy atomic occupancy factor
  \param atomic_multiplicity_factor 1/(number of symmetry operations which transfer atom into itself)
  \param local_coordinate_systems local coordinate systems, one per atom, changes with each change of atomic positions

  \param symmetry_operations - list of crystal symmetry operations (not incuding these which can be generated with
                               inversion center and pure translations)
  \param centrosymmetric - specifies if the space group is centrosymmetric
  \param inversionTranslation - specifies trnaslation vector of inversion center i.e. vector t in inversion operation {-1|t}
  \param h_vectors - list of points where structure factor has to be calculated
  \param f f[i] is a value of the structure factor at h_vectors[i] (has to be calculated)
  \param dTarget_dparam - derivatives of target functions with respect to atomic parameters
  \param dTarget_df derivatives of target function \f$ T \f$ with respect to structure factors \f$ F_k = A_k + iB_k \f$,
  \f$A_k, B_k \in \mathbb{R}\f$ given as:
  \f$ \frac{\partial T}{\partial A_k} + i \frac{\partial T}{\partial B_k}\f$
  (necessary to calculate dTarget_dparam without knowing target function (via chain rule)) 
  \param include_atom_contribution - spcifies which atoms contribution should be included
  \param nCores - number of threads to use
*/

void calculateSF(
    const UnitCell &unitCell,
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_factor,
    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
    const std::vector<Matrix3i>& symmetry_operations_rotation_matrix,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &h_vectors,
    const std::vector<Vector3i >& hkl_indices,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &dTarget_df,
    const std::vector<bool> &include_atom_contribution,
	int nThreads,
    const DerivativesSelector& derivativesSwitch,
    bool electron = false,
    const std::vector<int>& atomic_number = std::vector<int>());


void calculateSphericalTermsInFormFactors(
	const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
	const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
	const std::vector <double> h,
	std::vector< std::vector<REAL> >& f_core,
	std::vector< std::vector<REAL> >& f_sph_valence,
	const std::vector<int>& type_2_wfn_type,
	const std::vector<std::vector<REAL> >& def_val_slater_normalization,
	const std::vector<int>& typeMaxL);





/**
\param f_spherical - for each type spherical valence + core
*/

void calculateFormFactors(
	const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
	const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
	const std::vector<double>& f_spherical, // for each type spherical valence + core
	const std::vector<int>& atom_to_wfn_map,
	const std::vector<int>& atom_to_type_map,
	//const std::vector<Vector3<REAL> >& atomicPositions,
	//const std::vector<std::vector<REAL> >& atomic_displacement_parameters,
	//const std::vector<REAL>& atomic_occupancy,
	//const std::vector<REAL>& atomic_multiplicity_factor,
	const std::vector<Matrix3<REAL> >& local_coordinate_systems,
	const Vector3<REAL>& h_vectors,
	std::vector<std::complex<REAL> >& formFactors,
	const std::vector<bool>& includeAtom,
	const std::vector<int>& type_2_wfn_type,
	const std::vector<std::vector<REAL> >& def_val_slater_normalization,
	const std::vector<int>& typeMaxL);




/*
		   virtual void calculateFrac(
				const Vector3i &hkl,
				std::vector<std::complex<double> > &formFactors,
				const std::vector<bool> &includeAtom) const;

			virtual void calculateCart(
				const Vector3d &hkl,
				std::vector<std::complex<double> > &formFactors,
				const std::vector<bool> &includeAtom) const;
*/

/*
\brief Calculates structure factors and derivatives 
of a target function with respect to atomic parameters
within spherical atom approximation (so called independent atom model (IAM)).

It uses the same code as multipole model version. 
The spherical core electron density constribution to scattering is substituted 
with the independent atom form factor and other terms parameterized to give zero contribution.
It was intended for testing purposes.



\param atomicType scatterer types
\param atomTypeAnomalousScattering anomalous scattering for given atom type
\param atom_to_type_map atom_to_type_map[i] is an index of \p atomicType corresponding to i - th atom
\param atomicPositions positions of atoms in asymmetric unit(in Cartesian coordinate system, Angstroms)
\param atomic_displacement_parameters parameters describing atomic displacement(in Cartesian coordinate system)
\param atomic_occupancy atomic occupancy factor
\param atomic_multiplicity_factor 1 / (number of symmetry operations which transfer atom into itself)
\param symmetry_operations - list of crystal symmetry operations(not incuding these which can be generated with
    inversion center and pure translations)
\param centrosymmetric - specifies if the space group is centrosymmetric
\param inversionTranslation - specifies trnaslation vector of inversion center i.e.vector t in inversion operation{ -1 | t }
\param h_vectors - list of points where structure factor has to be calculated
\param f f[i] is a value of the structure factor at h_vectors[i](has to be calculated)
\param dTarget_dparam - derivatives of target functions with respect to atomic parameters
\param dTarget_df derivatives of target function \f$ T \f$ with respect to structure factors \f$ F_k = A_k + iB_k \f$,
                 \f$A_k, B_k \in \mathbb{ R }\f$ given as :
                 \f$ \frac{ \partial T }{\partial A_k} +i \frac{ \partial T }{\partial B_k}\f$
(necessary to calculate dTarget_dparam without knowing target function(via chain rule))
\param include_atom_contribution - spcifies which atoms contribution should be included
\param nCores - number of threads to use
*/

void calculateSF_IAM(
    const UnitCell &unitCell,
    const std::vector<std::string> &atomicType,
    const std::vector<std::complex<REAL> > &atomTypeAnomalousScattering,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_weight,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &h_vectors,
    const std::vector<Vector3i >& hkl_indices,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &dTarget_df,
    const std::vector<bool> &include_atom_contribution,
    int nCores = 1);

/**
\brief Calculates structure factors and their derivatives with respect to structural parameters.

For each hkl vector an user provided object \p singleHklResultsProcessor function is called
which collects the results for given h vector.

\param wfn_parameters sets of wave-function atomic parameters
\param type_parameters sets of parameters related to atom types
\param atom_to_wfn_map atom_to_wfn_map[i] is an index of wfn_parameters set corresponding to i-th atom
\param atom_to_type_map atom_to_type_map[i] is an index of type_parameters set corresponding to i-th atom
\param atomicPositions positions of atoms in asymmetric unit (in Cartesian coordinate system, Angstroms)
\param atomic_displacement_parameters parameters describing atomic displacement (in Cartesian coordinate system)
\param atomic_occupancy atomic occupancy factor
\param atomic_multiplicity_weight 1/(number of symmetry operations which transfer atom into itself)
\param local_coordinate_systems local coordinate systems, one per atom, changes with each change of atomic positions
\param symmetry_operations - list of crystal symmetry operations (not incuding these which can be generated with
inversion center and pure translations)
\param centrosymmetric - specifies if the space group is centrosymmetric
\param inversionTranslation - specifies trnaslation vector of inversion center i.e. vector t in inversion operation {-1|t}
\param h_vectors - list of points where structure factor has to be calculated
\param centeringMultiplier multipliers related to centering (e.g. 2 for all reflections in C centered unit cell for 
                           all non-extinct reflections)
\param include_atom_contribution flag for counting atom contribution
\param singleHklResultsProcessor object of user defined type intended for collecting results,
    for calculation for each Miller index its member function:
    \code{.cpp}
    void onPerHklCalculation(
    int hkl_idx,
    std::complex<double> &structureFactor,
    discamb::SfDerivativesAtHkl &derivatives);
    \endcode
    is called (hkl_idx starts from 0). See e.g. Tests/correctness/test_1/GradientDataCollector.h for an example implementation.
*/

template<typename T>
void calculateSF(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_weight,
    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &h_vectors,
    const std::vector<double> &centeringMultiplier,
    const std::vector<bool> &include_atom_contribution,
    T &singleHklResultsProcessor);


/**

\brief Calculates structure factors and their derivatives with respect to structural parameters
       within spherical atom approximation (so called independent atom model (IAM)).

For each hkl vector an user provided object \p singleHklResultsProcessor function is called
which collects the results for given h vector.

\param atomicType scatterer types
\param atomTypeAnomalousScattering anomalous scattering for given atom type
\param atom_to_type_map atom_to_type_map[i] is an index of \p atomicType corresponding to i - th atom
\param atomicPositions positions of atoms in asymmetric unit (in Cartesian coordinate system, Angstroms)
\param atomic_displacement_parameters parameters describing atomic displacement (in Cartesian coordinate system)
\param atomic_occupancy atomic occupancy factor
\param atomic_multiplicity_weight 1/(number of symmetry operations which transfer atom into itself)
\param symmetry_operations - list of crystal symmetry operations (not incuding these which can be generated with
inversion center and pure translations)
\param centrosymmetric - specifies if the space group is centrosymmetric
\param inversionTranslation - specifies trnaslation vector of inversion center i.e. vector t in inversion operation {-1|t}
\param h_vectors - list of points where structure factor has to be calculated
\param centeringMultiplier multipliers related to centering (e.g. 2 for all reflections in C centered unit cell for
all non-extinct reflections)
\param include_atom_contribution flag for counting atom contribution
\param singleHklResultsProcessor object of user defined type intended for collecting results,
for calculation for each Miller index its member function:
\code{.cpp}
void onPerHklCalculation(
int hkl_idx,
std::complex<double> &structureFactor,
discamb::SfDerivativesAtHkl &derivatives);
\endcode
is called (hkl_idx starts from 0). See e.g. Tests/correctness/test_1/GradientDataCollector.h for an example implementation.
*/


template<typename T>
void calculateSF_IAM(
    const std::vector<std::string> &atomicType,
    const std::vector<std::complex<REAL> > &atomTypeAnomalousScattering,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_weight,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &h_vectors,
    const std::vector<double> &centeringMultiplier,
    const std::vector<bool> &include_atom_contribution,
    T &singleHklResultsProcessor);

/**

*/

void pre_hkl_loop_sf_calc(
	const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
	const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
	const std::vector<int>& atom_to_wfn_map,
	const std::vector<int>& atom_to_type_map,
	std::vector<int>& type_2_wfn_type,
	std::vector<std::vector<REAL> >& def_val_slater_normalization,
	std::vector<int>& typeMaxL
);


static std::complex<double> calculateDeformationValence(
    const std::vector<std::vector<REAL> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
    const std::vector<REAL>& g_functions_and_slater_normalization,
    int maxL,
    const std::vector<std::vector<double> >& sphericalHarmonics);


private:
    /*
    3 - IAM
    */
    //int mSphericalTermCalculationMode;

    bool mUseIAM;
    std::vector<std::string> mIamAtomType;
    std::vector<int> mAtomToIamTypeMap;
    std::vector<NGaussianFormFactor> mIamFormFactors;
    

    // [thread][l][l+m]
    std::vector<std::vector<std::vector<double> > > mSphericalHarmonicsData;

    static void calculateGlobalCoordinatesPlm(
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
        const std::vector<int>& atom_to_type_map,
        const std::vector<Matrix3<REAL> >& local_coordinate_systems,
        std::vector< std::vector<std::vector<double> > >& atomPlms);

    // atom_selection[i] if i-th atom only nonzero Plms are P10 and P20
    static void select_P10P20_atoms(
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters, 
        const std::vector<int>& atom_to_type_map, 
        std::vector<bool> &atom_selection);

    void calculateSF_SerialAcentric(
        const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
        const std::vector<int> &atom_to_wfn_map,
        const std::vector<int> &atom_to_type_map,
        const std::vector<Vector3<REAL> > &atomicPositions,
        const std::vector< std::vector<Vector3d> >& r_atom_symm,
        const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
        const std::vector<REAL> &atomic_occupancy,
        const std::vector<REAL> &atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> > &local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
        //const std::vector<Vector3<REAL> > &h_vectors,
        const std::vector<std::vector<Vector3<REAL> > > & h_vector_lines,
        const std::vector< std::vector< Vector3i > > & h_vector_lines_int,
        const std::vector < std::vector <int> > &line_to_orginal_hkl_list_idx,
        int line_direction,
        const Vector3d &lineStepCart,
        std::vector<std::complex<REAL> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        //const std::vector<std::complex<REAL> > &dTarget_df,
        const std::vector< std::vector<std::complex<REAL> > >& dTarget_df_lines,
        const std::vector<bool> &include_atom_contribution,
        const std::vector<int> &type_2_wfn_type,
        const std::vector<std::vector<REAL> > &def_val_slater_normalization,
        std::vector<std::vector<double> > &sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int> &typeMaxL,
        const std::vector< std::vector<std::vector<double> > > &atomPlms,
        int &executionTime,
        const DerivativesSelector &derivativesSwitch,
        bool electron,
        const std::vector<int>& atomic_number);

    void electronScatteringAt000(
        const std::vector<int>& atomic_numbers,
        std::vector<double>& f);

    void calculateSF_SerialAcentric0(
        const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
        const std::vector<int> &atom_to_wfn_map,
        const std::vector<int> &atom_to_type_map,
        const std::vector<Vector3<REAL> > &atomicPositions,
        const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
        const std::vector<REAL> &atomic_occupancy,
        const std::vector<REAL> &atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> > &local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
        const std::vector<Vector3<REAL> > &h_vectors,
        std::vector<std::complex<REAL> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        const std::vector<std::complex<REAL> > &dTarget_df,
        const std::vector<bool> &include_atom_contribution,
        const std::vector<int> &type_2_wfn_type,
        const std::vector<std::vector<REAL> > &def_val_slater_normalization,
        std::vector<std::vector<double> > &sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int> &typeMaxL);


    void calculateSF_SerialCentrosymmetric(
        const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
        const std::vector<int> &atom_to_wfn_map,
        const std::vector<int> &atom_to_type_map,
        const std::vector<Vector3<REAL> > &atomicPositions,
        const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
        const std::vector<REAL> &atomic_occupancy,
        const std::vector<REAL> &atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> > &local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
        const std::vector<Vector3<REAL> > &h_vectors,
        std::vector<std::complex<REAL> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        const std::vector<std::complex<REAL> > &dTarget_df,
        const std::vector<bool> &include_atom_contribution,
        const std::vector<int> &type_2_wfn_type,
        const std::vector<std::vector<REAL> > &def_val_slater_normalization,
        std::vector<std::vector<double> > &sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int> &typeMaxL);

    void calculateSF_SerialCentrosymmetric_ordered_hkl(
        const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
        const std::vector<int>& atom_to_wfn_map,
        const std::vector<int>& atom_to_type_map,
        const std::vector<Vector3<REAL> >& atomicPositions,
        const std::vector< std::vector<Vector3d> >& r_atom_symm,
        const std::vector<std::vector<REAL> >& atomic_displacement_parameters,
        const std::vector<REAL>& atomic_occupancy,
        const std::vector<REAL>& atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> >& local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation>& symmetry_operations,

        const std::vector<std::vector<Vector3<REAL> > >& h_vector_lines,
        const std::vector< std::vector< Vector3i > >& h_vector_lines_int,
        const std::vector < std::vector <int> >& line_to_orginal_hkl_list_idx,
        int line_direction,
        const Vector3d& lineStepCart,

        std::vector<std::complex<REAL> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        //const std::vector<std::complex<REAL> >& dTarget_df,
        const std::vector< std::vector<std::complex<REAL> > >& dTarget_df_lines,
        const std::vector<bool>& include_atom_contribution,
        const std::vector<int>& type_2_wfn_type,
        const std::vector<std::vector<REAL> >& def_val_slater_normalization,
        std::vector<std::vector<double> >& sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int>& typeMaxL,
        const std::vector< std::vector<std::vector<double> > >& atomPlms,
        const DerivativesSelector& derivativesSwitch);

    void calculateSF_SerialCentrosymmetric_ordered_hkl2(
        const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
        const std::vector<int>& atom_to_wfn_map,
        const std::vector<int>& atom_to_type_map,
        const std::vector<Vector3<REAL> >& atomicPositions,
        const std::vector< std::vector<Vector3d> >& r_atom_symm,
        const std::vector<std::vector<REAL> >& atomic_displacement_parameters,
        const std::vector<REAL>& atomic_occupancy,
        const std::vector<REAL>& atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> >& local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation>& symmetry_operations,

        const std::vector<std::vector<Vector3<REAL> > >& h_vector_lines,
        const std::vector< std::vector< Vector3i > >& h_vector_lines_int,
        const std::vector < std::vector <int> >& line_to_orginal_hkl_list_idx,
        int line_direction,
        const Vector3d& lineStepCart,

        std::vector<std::complex<REAL> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        //const std::vector<std::complex<REAL> >& dTarget_df,
        const std::vector< std::vector<std::complex<REAL> > >& dTarget_df_lines,
        const std::vector<bool>& include_atom_contribution,
        const std::vector<int>& type_2_wfn_type,
        const std::vector<std::vector<REAL> >& def_val_slater_normalization,
        std::vector<std::vector<double> >& sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int>& typeMaxL,
        const std::vector< std::vector<std::vector<double> > >& atomPlms,
        const DerivativesSelector& derivativesSwitch,
        bool electron,
        const std::vector<int>& atomic_number);

    void calculateSF_SerialSymmetryCenterNotAtOrigin2(
        const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
        const std::vector<int>& atom_to_wfn_map,
        const std::vector<int>& atom_to_type_map,
        const std::vector<Vector3<REAL> >& atomicPositions,
        const std::vector< std::vector<Vector3d> >& r_atom_symm,
        const std::vector<std::vector<REAL> >& atomic_displacement_parameters,
        const std::vector<REAL>& atomic_occupancy,
        const std::vector<REAL>& atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> >& local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation>& symmetry_operations,

        const std::vector<std::vector<Vector3<REAL> > >& h_vector_lines,
        const std::vector< std::vector< Vector3i > >& h_vector_lines_int,
        const std::vector < std::vector <int> >& line_to_orginal_hkl_list_idx,
        int line_direction,
        const Vector3d& lineStepCart,

        std::vector<std::complex<REAL> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        //const std::vector<std::complex<REAL> >& dTarget_df,
        const std::vector< std::vector<std::complex<REAL> > >& dTarget_df_lines,
        const std::vector<bool>& include_atom_contribution,
        const std::vector<int>& type_2_wfn_type,
        const std::vector<std::vector<REAL> >& def_val_slater_normalization,
        std::vector<std::vector<double> >& sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int>& typeMaxL,
        const std::vector< std::vector<std::vector<double> > >& atomPlms,
        const DerivativesSelector& derivativesSwitch);

    void calculateSF_SerialSymmetryCenterNotAtOrigin(
        const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
        const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
        const std::vector<int> &atom_to_wfn_map,
        const std::vector<int> &atom_to_type_map,
        const std::vector<Vector3<REAL> > &atomicPositions,
        const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
        const std::vector<REAL> &atomic_occupancy,
        const std::vector<REAL> &atomic_multiplicity_weight,
        const std::vector<Matrix3<REAL> > &local_coordinate_systems,
        const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
        const Vector3<REAL> &inversionTranslation,
        const std::vector<Vector3<REAL> > &h_vectors,
        std::vector<std::complex<REAL> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        const std::vector<std::complex<REAL> > &dTarget_df,
        const std::vector<bool> &include_atom_contribution,
        const std::vector<int> &type_2_wfn_type,
        const std::vector<std::vector<REAL> > &def_val_slater_normalization,
        std::vector<std::vector<double> > &sphericalHarmonicsData,// memory for spherical harmonics data
        const std::vector<int> &typeMaxL,
        const DerivativesSelector& derivativesSwitch,
        bool electron,
        const std::vector<int>& atomic_number);



    template<int l_max>
    static std::complex<REAL> combine_multipolar_terms(
        const std::vector<std::vector<REAL> > &p_lm,
        const std::vector<REAL> &radial,
        const std::vector<std::vector<REAL> > & sphericalHarmonicBuffer);

    std::complex<REAL> calculateDeformationValence(
        const std::vector<std::vector<REAL> > &p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
        const std::vector<REAL> &g_functions_and_slater_normalization,
        const Matrix3<REAL>  &local_coordinates_system,
        const Vector3<REAL> &normalized_h_vector,
        int maxL,
        std::vector<std::vector<double> > & sphericalHarmonicBuffer);


    /* 
    point_group_idx - 0 for undefined point group, then numeration follows order in CrystalStructure/crystallographic_point_groups.cpp
    */
    void calculateDeformationValence(
        int point_group_idx,
        const std::vector<std::vector<REAL> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
        const std::vector<REAL>& g_functions_and_slater_normalization,
        int maxL,
        const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
        const std::vector<int> &symmetryOperationsOrder,
        std::vector<std::complex<double> > &defVal);

    void calculateDefVal_general(
        const std::vector<std::vector<REAL> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
        const std::vector<REAL>& g_functions_and_slater_normalization,
        int maxL,
        const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
        const std::vector<int>& symmetryOperationsOrder,
        std::vector<std::complex<double> >& defVal);

    void calculateDefVal_pg_1(
        const std::vector<std::vector<REAL> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
        const std::vector<REAL>& g_functions_and_slater_normalization,
        int maxL,
        const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
        const std::vector<int>& symmetryOperationsOrder,
        std::vector<std::complex<double> >& defVal);

    void calculateDefVal_pg_222(
        const std::vector<std::vector<REAL> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
        const std::vector<REAL>& g_functions_and_slater_normalization,
        int maxL,
        const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
        const std::vector<int>& symmetryOperationsOrder,
        std::vector<std::complex<double> >& defVal);


    void pre_atom_loop_sf_calc(
        //in:
        const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,
        const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams,
        const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
        const std::vector<int> &type_2_wfn,
        const std::vector<std::vector<REAL> > &def_val_slater_normalization,
        const Vector3<REAL> &hVector,
        REAL hVectorLength,
        //out:
        std::vector<REAL> &wfn_spherical_core_sf,
        std::vector<REAL> &wfn_spherical_valence_sf,
        std::vector<std::vector<REAL> > &g_functions_and_slater_norm,
        std::vector<Vector3<REAL> > &rotated_h,
        std::vector<Vector3<REAL> > &rotated_normalized_h,
        std::vector<REAL> &translation_factor,
        std::vector<std::vector<REAL> > &adpMultipliers
    );

    
    inline REAL calc_temperature_factor(
        const Vector3<REAL> &h,
        REAL hVectorLength,
        const std::vector<REAL> &adp)
    {
        REAL temperature_factor = 1;
        if (adp.size() == 1)
            temperature_factor = exp(-hVectorLength * hVectorLength * adp[0]);
        else
        {
            const REAL temperature_factor_exponent =
                h[0] * (adp[0] * h[0] + adp[3] * h[1] + adp[4] * h[2]) +
                h[1] * (adp[3] * h[0] + adp[1] * h[1] + adp[5] * h[2]) +
                h[2] * (adp[4] * h[0] + adp[5] * h[1] + adp[2] * h[2]);

            temperature_factor = exp(-temperature_factor_exponent);
        }

        return temperature_factor;
    }



    void add_contribution_to_position_derivatives(
        Vector3<REAL> &position_derivatives,
        const std::complex<REAL> dTarget_dF,
        const std::complex<REAL> &atomic_f,
        const Vector3<REAL> &h);


    void add_contribution_to_adp_derivatives(
        std::vector<std::complex<REAL> > &adp_derivatives,
        const std::complex<REAL> &dTarget_dF,
        const std::complex<REAL> &atomic_f,
        const Vector3<REAL> &h);


    void process_adp_derivatives(
        std::complex<REAL> *pre_derivatives,
        const std::complex<REAL> &atomic_f,
        const Vector3<REAL> &h,
        REAL h_length,
        int n_adp_components);


    void add_contribution_to_occupancy_derivative(
        REAL &occupancy_derivative,
        const std::complex<REAL> &dTarget_dF,
        const std::complex<REAL> &atomic_f_divided_by_occupancy);


};



//#############################################################################
//
//                    TEMPLATE IMPLEMENTATION
//
//#############################################################################

template<typename T>
void HansenCoppens_SF_Engine4::calculateSF(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &typeParameters,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &adps,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_weight,
    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation, 
    const std::vector<Vector3<REAL> > &hVectors,
    const std::vector<double> &centeringMultipliers,
    const std::vector<bool> &include_atom_contribution,
    T &singleHklResultsProcessor)
{
    //cout << "calls HansenCoppens_SF_Engine4_6::calculateSF" << endl;

    std::complex<REAL> structureFactor, inversionTranslationPhaseFactor;
    SfDerivativesAtHkl derivatives;
    double centeringMultiplier;
    bool inversionNotAtOrigin = false;
    int nSymmOps = symOps.size();
    std::vector<std::vector<REAL> > adpMultipliers(nSymmOps, std::vector<REAL>(6));
    if(centrosymmetric)
        if( inversionTranslation != Vector3d(0.0,0.0,0.0) )
            inversionNotAtOrigin = true;
    

    // allocate memory buffers for spherical harmonics calculations

    mSphericalHarmonicsData.resize(1);
    int maxNL, nL, typeIdx, nTypes = typeParameters.size();
    maxNL = 0;
    for (typeIdx = 0; typeIdx<nTypes; typeIdx++)
    {
        nL = typeParameters[typeIdx].p_lm.size();
        if (nL>maxNL)
            maxNL = nL;
    }

    mSphericalHarmonicsData[0].resize(maxNL);
    for (int l = 0; l<maxNL; l++)
        mSphericalHarmonicsData[0][l].resize(2 * l + 1);
    
    std::vector<std::vector<REAL> > def_val_slater_normalization;
    std::vector<int> type_2_wfn_type(typeParameters.size(), 100);
    std::vector<int> type_max_L;
    pre_hkl_loop_sf_calc(wfnParameters, typeParameters, atom_to_wfn_map, atom_to_type_map,
                         type_2_wfn_type, def_val_slater_normalization, type_max_L);

    //

    int nHklVectors = hVectors.size();
    int nAtoms = atom_to_wfn_map.size();


    // allocate derivatives container

    derivatives.adpDerivatives.resize(nAtoms);
    derivatives.atomicPostionDerivatives.resize(nAtoms);
    derivatives.occupancyDerivatives.resize(nAtoms);

    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        derivatives.adpDerivatives[atomIdx].resize(adps[atomIdx].size());



    const std::complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
    const std::complex<REAL> sqrtMinusOne(0,1);
    std::complex<REAL> fAtomSphericalAndAnomalous, term1, term2;
    std::complex<REAL> transformedAtomF, xyzDerivativesMultiplier;
    const REAL two_pi = 2 * REAL(M_PI);
    const REAL two_pi_sqare = 2 * REAL(M_PI*M_PI);
    std::vector<Vector3<REAL> > rotated_h(symOps.size());
    std::vector<Vector3<REAL> > rotated_normalized_h(symOps.size());
    std::vector<REAL> translation_factor(symOps.size());

    std::vector<REAL> wfn_spherical_core_sf(wfnParameters.size());
    std::vector<REAL> wfn_spherical_valence_sf(typeParameters.size());
    std::vector<std::vector<REAL> > g_functions_and_slater_norm(typeParameters.size(), std::vector<REAL>(5));

    

    //


    for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
    {
        REAL hVectorLength = sqrt(hVectors[hklIndex] * hVectors[hklIndex]);
        structureFactor = 0.0;
        centeringMultiplier = centeringMultipliers[hklIndex];

        pre_atom_loop_sf_calc(
            //in:
            wfnParameters, typeParameters, symOps, type_2_wfn_type, def_val_slater_normalization, hVectors[hklIndex], hVectorLength,
            //out:
            wfn_spherical_core_sf, wfn_spherical_valence_sf, g_functions_and_slater_norm,
            rotated_h, rotated_normalized_h, translation_factor, adpMultipliers);

        if(inversionNotAtOrigin)
            inversionTranslationPhaseFactor = exp(two_pi_i*(hVectors[hklIndex] * inversionTranslation));

        for (int atomIdx = 0; atomIdx < atomicPositions.size(); atomIdx++)
        {
            std::complex<REAL> adp_derivatives[6];
            std::complex<REAL> xyz_derivatives[3];
            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;

            if (!include_atom_contribution[atomIdx])
                continue;
            std::complex<REAL> atomicFContrib = 0;

            int atomWfnIdx = atom_to_wfn_map[atomIdx];
            int atomTypeIdx = atom_to_type_map[atomIdx];




            REAL atomWeight = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx] * centeringMultiplier;

            // calculates h vector direction independent part of atomic scattering factor
            //                (identical for all symmetry dependent atoms)

            REAL atom_f_core = wfn_spherical_core_sf[atomWfnIdx];
            REAL atom_f_sph_val = wfn_spherical_valence_sf[atomTypeIdx];
            atom_f_sph_val *= typeParameters[atomTypeIdx].p_val;


            std::complex<REAL> anomalousScattering = wfnParameters[atomWfnIdx].anomalous_scattering;
            fAtomSphericalAndAnomalous = atom_f_core + atom_f_sph_val + anomalousScattering;

            // end of h direction independent part of atomic scattering factor calculation


            // prepare derivatives data

            int n_adp_components = adps[atomIdx].size();
            for (int i = 0; i<n_adp_components; i++)
                adp_derivatives[i] = 0.0;
            xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

            derivatives.occupancyDerivatives[atomIdx] = 0;

            for (int symOpIdx = 0; symOpIdx < nSymmOps ; symOpIdx++)
            {

                const Vector3<REAL> & rotated_h_ref = rotated_h[symOpIdx];

                atomic_phase_factor_phase = two_pi*(rotated_h_ref*atomicPositions[atomIdx] +
                                                    translation_factor[symOpIdx]);

                atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                std::complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);


                std::complex<REAL> atom_f_def_val = calculateDeformationValence(typeParameters[atomTypeIdx].p_lm,
                                                                           g_functions_and_slater_norm[atomTypeIdx],
                                                                           local_coordinate_systems[atomIdx],
                                                                           rotated_normalized_h[symOpIdx],
                                                                           type_max_L[atomTypeIdx], mSphericalHarmonicsData[0]);
                // temperature factor
                REAL temperature_factor;
                adps[atomIdx].empty() ?
                    temperature_factor = 1.0 :
                    temperature_factor = calc_temperature_factor(rotated_h_ref, hVectorLength, adps[atomIdx]);

                //-----------

                if(!centrosymmetric)
                {
                    transformedAtomF = (fAtomSphericalAndAnomalous + atom_f_def_val)
                                        * atomic_position_phase_factor * temperature_factor;

                    xyz_derivatives[0] += rotated_h_ref[0] * transformedAtomF;
                    xyz_derivatives[1] += rotated_h_ref[1] * transformedAtomF;
                    xyz_derivatives[2] += rotated_h_ref[2] * transformedAtomF;

                }
                else
                {
                    if (inversionNotAtOrigin)
                    {
                    
                        term1 = (fAtomSphericalAndAnomalous + atom_f_def_val) * atomic_position_phase_factor;
                        term2 = (fAtomSphericalAndAnomalous + conj(atom_f_def_val)) * conj(atomic_position_phase_factor) *
                                inversionTranslationPhaseFactor;

                        transformedAtomF = temperature_factor * ( term1 + term2 );

                        xyzDerivativesMultiplier = temperature_factor*(term1 - term2);

                    }
                    else
                    {
                        transformedAtomF = 2 * temperature_factor * ( fAtomSphericalAndAnomalous * atomic_phase_factor_real +
                                                                      atom_f_def_val.real() * atomic_phase_factor_real -
                                                                      atom_f_def_val.imag() * atomic_phase_factor_im );

                        xyzDerivativesMultiplier = ((fAtomSphericalAndAnomalous + atom_f_def_val)*atomic_position_phase_factor -
                            (fAtomSphericalAndAnomalous + conj(atom_f_def_val))*conj(atomic_position_phase_factor))*temperature_factor;

                        xyzDerivativesMultiplier = 2.0 * (fAtomSphericalAndAnomalous*sqrtMinusOne*atomic_position_phase_factor.imag()+
                                                          sqrtMinusOne*(atom_f_def_val*atomic_position_phase_factor).imag())*temperature_factor;
                    }

                    xyz_derivatives[0] += rotated_h_ref[0] * xyzDerivativesMultiplier;
                    xyz_derivatives[1] += rotated_h_ref[1] * xyzDerivativesMultiplier;
                    xyz_derivatives[2] += rotated_h_ref[2] * xyzDerivativesMultiplier;


                }

                atomicFContrib += atomWeight * transformedAtomF;

                if (n_adp_components == 1)
                    adp_derivatives[0] -= hVectorLength * hVectorLength * transformedAtomF;
                else if (n_adp_components == 6)
                {
                    adp_derivatives[0] -= rotated_h_ref[0] * rotated_h_ref[0] * transformedAtomF;
                    adp_derivatives[1] -= rotated_h_ref[1] * rotated_h_ref[1] * transformedAtomF;
                    adp_derivatives[2] -= rotated_h_ref[2] * rotated_h_ref[2] * transformedAtomF;
                    adp_derivatives[3] -= 2 * rotated_h_ref[0] * rotated_h_ref[1] * transformedAtomF;
                    adp_derivatives[4] -= 2 * rotated_h_ref[0] * rotated_h_ref[2] * transformedAtomF;
                    adp_derivatives[5] -= 2 * rotated_h_ref[1] * rotated_h_ref[2] * transformedAtomF;
                }

                derivatives.occupancyDerivatives[atomIdx] += transformedAtomF * atomic_multiplicity_weight[atomIdx] * centeringMultiplier;
            } // symmetry operations

            if (include_atom_contribution[atomIdx])
            {
                structureFactor += atomicFContrib;
                //realFContrib += atomicFContrib.real();
                //imagFContrib += atomicFContrib.imag();


                //complex<REAL> aux(dTarget_df[hklIndex] * atomWeight);
                // adp
                if (n_adp_components > 0)
                {
                    if (n_adp_components == 1)
                        derivatives.adpDerivatives[atomIdx][0] = two_pi_sqare * atomWeight * adp_derivatives[0];
                       
                    else
                        for (int i = 0; i<6; ++i)
                            derivatives.adpDerivatives[atomIdx][i] = two_pi_sqare * atomWeight * adp_derivatives[i];
                            
                }
                // xyz
                std::complex<REAL> aux( atomWeight * two_pi_i );
                
                for (int i = 0; i<3; ++i)
                    derivatives.atomicPostionDerivatives[atomIdx][i] = aux*xyz_derivatives[i];
                    //dTarget_dparam[atomIdx].atomic_position_derivatives[i] += (aux*xyz_derivatives[i]).real();
            }

            

        } // symetrically independent atoms

        singleHklResultsProcessor.onPerHklCalculation(hklIndex, structureFactor, derivatives);

    } // h vectors



}

template<typename T>
void HansenCoppens_SF_Engine4::calculateSF_IAM(
    const std::vector<std::string> &atomicType,
    const std::vector<std::complex<REAL> > &atomTypeAnomalousScattering,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_factor,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &h_vectors,
    const std::vector<double> &centeringMultiplier,
    const std::vector<bool> &include_atom_contribution,
    T &singleHklResultsProcessor)
{
    mUseIAM = true;
    mIamAtomType = atomicType;
    mAtomToIamTypeMap = atom_to_type_map;


    std::vector<sf_engine_data_types::HC_WfnParam> wfnParams(atomicType.size());
    std::vector<sf_engine_data_types::HC_TypeParam> typeParams(1);
    std::vector<int> atomToWfnMap = atom_to_type_map;
    std::vector<int> atomToTypeMap(atomicPositions.size(), 0);
    Matrix3d identity;
    identity.setToIdentity();
    std::vector<Matrix3d> localCoordinateSystems(atomicPositions.size(), identity);

    int iamTypeIdx, nIamTypes = atomicType.size();

    mIamFormFactors.resize(nIamTypes);

    for (iamTypeIdx = 0; iamTypeIdx < nIamTypes; iamTypeIdx++)
    {
        wfnParams[iamTypeIdx].anomalous_scattering = atomTypeAnomalousScattering[iamTypeIdx];
        if (n_gaussian_form_factors_table::hasFormFactor(mIamAtomType[iamTypeIdx]))
            mIamFormFactors[iamTypeIdx] = n_gaussian_form_factors_table::getFormFactor(mIamAtomType[iamTypeIdx]);
        else
            on_error::throwException(
                std::string("request for Gaussian type atomic form factor parameter for unknown atom type: ")
                + mIamAtomType[iamTypeIdx], __FILE__, __LINE__);
    }


    calculateSF(wfnParams, typeParams, atomToWfnMap, atomToTypeMap, atomicPositions,
                atomic_displacement_parameters, atomic_occupancy, atomic_multiplicity_factor,
                localCoordinateSystems, symmetry_operations, centrosymmetric, inversionTranslation,
                h_vectors, centeringMultiplier, include_atom_contribution, singleHklResultsProcessor);
}


template<int l_max>
std::complex<REAL> HansenCoppens_SF_Engine4::combine_multipolar_terms(
    const std::vector<std::vector<REAL> > &p_lm,
    const std::vector<REAL> &radial,
    const std::vector<std::vector<REAL> > &sphericalHarmonics)
{
    REAL resultReal = 0;
    REAL resultImag = 0;
    REAL angular;


    resultReal += radial[0] * p_lm[0][0] * sphericalHarmonics[0][0];


    if (l_max > 0)
    {

        angular = p_lm[1][0] * sphericalHarmonics[1][0] +
            p_lm[1][1] * sphericalHarmonics[1][1] +
            p_lm[1][2] * sphericalHarmonics[1][2];

        resultImag += radial[1] * angular;


        if (l_max > 1)
        {


            //angular = 0;
            //for (int i = 0; i<5; i++)
            //    angular += p_lm[2][i] * sphericalHarmonics[2][i];

            angular = p_lm[2][0] * sphericalHarmonics[2][0] +
                p_lm[2][1] * sphericalHarmonics[2][1] +
                p_lm[2][2] * sphericalHarmonics[2][2] +
                p_lm[2][3] * sphericalHarmonics[2][3] +
                p_lm[2][4] * sphericalHarmonics[2][4];


            resultReal -= radial[2] * angular;

            if (l_max > 2)
            {

                //angular = 0;
                //for (int i = 0; i<7; i++)
                //    angular += p_lm[3][i] * sphericalHarmonics[3][i];

                angular = p_lm[3][0] * sphericalHarmonics[3][0] +
                    p_lm[3][1] * sphericalHarmonics[3][1] +
                    p_lm[3][2] * sphericalHarmonics[3][2] +
                    p_lm[3][3] * sphericalHarmonics[3][3] +
                    p_lm[3][4] * sphericalHarmonics[3][4] +
                    p_lm[3][5] * sphericalHarmonics[3][5] +
                    p_lm[3][6] * sphericalHarmonics[3][6];


                resultImag -= radial[3] * angular;


                if (l_max > 3)
                {


                    //angular = 0;
                    //for (int i = 0; i<9; i++)
                    //    angular += p_lm[4][i] * sphericalHarmonics[4][i];

                    angular = p_lm[4][0] * sphericalHarmonics[4][0] +
                        p_lm[4][1] * sphericalHarmonics[4][1] +
                        p_lm[4][2] * sphericalHarmonics[4][2] +
                        p_lm[4][3] * sphericalHarmonics[4][3] +
                        p_lm[4][4] * sphericalHarmonics[4][4] +
                        p_lm[4][5] * sphericalHarmonics[4][5] +
                        p_lm[4][6] * sphericalHarmonics[4][6] +
                        p_lm[4][7] * sphericalHarmonics[4][7] +
                        p_lm[4][8] * sphericalHarmonics[4][8];


                    resultReal += radial[4] * angular;

                } // l_max>3
            } // l_max>2
        } // l_max>1
    } // l_max>0


    return 4 * M_PI* std::complex<REAL>(resultReal, resultImag);//d*REAL(4)*REAL(M_PI)*complex<REAL>(resultReal, resultImag);
}
/** @}*/
} // namespace discamb

#endif /*_DISCAMB_SCATTERING_HansenCoppens_SF_Engine4_H_*/
