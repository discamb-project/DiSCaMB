#include "discamb/Scattering/HansenCoppens_SF_Engine.h"
#include "discamb/Scattering/SlaterTypeOrbitalScattering.h"
#include "discamb/HC_Model/HC_WfnData.h"
#include "discamb/Scattering/NGaussianFormFactorsTable.h"
#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/Timer.h"

#include <cmath>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <iomanip>

#if defined(_OPENMP)
#include <omp.h>
#endif


using namespace std;

namespace discamb {

HansenCoppens_SF_Engine::HansenCoppens_SF_Engine()
{
    mUseIAM = false;
}

HansenCoppens_SF_Engine::~HansenCoppens_SF_Engine()
{
}



inline void HansenCoppens_SF_Engine::add_contribution_to_occupancy_derivative(
    REAL &occupancy_derivative,
    const complex<REAL> &dTarget_dF,
    const complex<REAL> &atomic_f_divided_by_occupancy)
{
    occupancy_derivative += (dTarget_dF * atomic_f_divided_by_occupancy).real();
}

inline void HansenCoppens_SF_Engine::add_contribution_to_position_derivatives(
    Vector3<REAL> &position_derivatives,
    const complex<REAL> dTarget_dF,
    const complex<REAL> &atomic_f,
    const Vector3<REAL> &h)
{
    static const complex<REAL> two_pi_i = REAL(2*REAL(M_PI))*complex<REAL>(0,1);
    complex<REAL> df_dparam;

    for(int k=0;k<3;k++) {
        df_dparam = two_pi_i * h[k] * atomic_f;
        position_derivatives[k] += (dTarget_dF * df_dparam).real();
    }
}


inline void HansenCoppens_SF_Engine::add_contribution_to_adp_derivatives(
    std::vector<std::complex<REAL> > &adp_derivatives,
    const std::complex<REAL> &dTarget_dF,
    const std::complex<REAL> &atomic_f,
    const Vector3<REAL> &h)
{
    complex<REAL> df_dparam;
    REAL hVectorLength = sqrt(h*h);

    if(adp_derivatives.size() == 1) {
        df_dparam = -hVectorLength * hVectorLength * atomic_f;
        adp_derivatives[0] += dTarget_dF * df_dparam;

    }
    else {
        for (int k = 0; k < 3; k++) {
            df_dparam = -h[k] * h[k] * atomic_f;
            adp_derivatives[k] += dTarget_dF * df_dparam;
        }

        // U_12
        df_dparam = -2*h[0] * h[1] * atomic_f;
        adp_derivatives[3] += dTarget_dF * df_dparam;

        // U_13
        df_dparam = -2*h[0] * h[2] * atomic_f;
        adp_derivatives[4] += dTarget_dF * df_dparam;

        // U_23
        df_dparam = -2*h[1] * h[2] * atomic_f;
        adp_derivatives[5] += dTarget_dF * df_dparam;
    }
}

inline void HansenCoppens_SF_Engine::process_adp_derivatives( std::complex<REAL> *pre_derivatives,
                                     const std::complex<REAL> &atomic_f,
                                     const Vector3<REAL> &h,
                                     REAL h_length,
                                     int n_adp_components)
{

    if(n_adp_components==1)
    {
        pre_derivatives[0] -= h_length*h_length*atomic_f;
        return;
    }
    
     
    pre_derivatives[0] -= h[0]*h[0]*atomic_f;
    pre_derivatives[1] -= h[1]*h[1]*atomic_f;
    pre_derivatives[2] -= h[2]*h[2]*atomic_f;
    pre_derivatives[3] -= 2*h[0]*h[1]*atomic_f;
    pre_derivatives[4] -= 2*h[0]*h[2]*atomic_f;
    pre_derivatives[5] -= 2*h[1]*h[2]*atomic_f;
}


std::complex<REAL> HansenCoppens_SF_Engine::calculateDeformationValence(
    const std::vector<std::vector<REAL> > &p_lm,
    const std::vector<REAL> &g_functions_and_slater_normalization,
    const Matrix3<REAL>  &local_coordinates_system,
    const Vector3<REAL> &normalized_h_vector,
    int max_l,
    std::vector<std::vector<double> > &sphericalHarmonicBuffer)
{
    if(max_l<0)
        return 0;

    const Matrix3<REAL>  &lcs = local_coordinates_system;
    const REAL x = (lcs(0, 0)*normalized_h_vector(0) + lcs(1, 0)*normalized_h_vector(1) + lcs(2, 0)*normalized_h_vector(2));//hRotated[0];
    const REAL y = (lcs(0, 1)*normalized_h_vector(0) + lcs(1, 1)*normalized_h_vector(1) + lcs(2, 1)*normalized_h_vector(2));//hRotated[1];
    const REAL z = (lcs(0, 2)*normalized_h_vector(0) + lcs(1, 2)*normalized_h_vector(1) + lcs(2, 2)*normalized_h_vector(2));//hRotated[2];

    Vector3d h(x,y,z);

    switch(max_l)
    {
    case 0:
        real_spherical_harmonics::getDensityNormalized<0>(h, sphericalHarmonicBuffer);
        return combine_multipolar_terms<0>(p_lm,g_functions_and_slater_normalization, sphericalHarmonicBuffer);
    case 1:
        real_spherical_harmonics::getDensityNormalized<1>(h, sphericalHarmonicBuffer);
        return combine_multipolar_terms<1>(p_lm,g_functions_and_slater_normalization, sphericalHarmonicBuffer);
    case 2:
        real_spherical_harmonics::getDensityNormalized<2>(h, sphericalHarmonicBuffer);
        return combine_multipolar_terms<2>(p_lm,g_functions_and_slater_normalization, sphericalHarmonicBuffer);
    case 3:
        real_spherical_harmonics::getDensityNormalized<3>(h, sphericalHarmonicBuffer);
        return combine_multipolar_terms<3>(p_lm,g_functions_and_slater_normalization, sphericalHarmonicBuffer);
    case 4:
        real_spherical_harmonics::getDensityNormalized<4>(h, sphericalHarmonicBuffer);
        return combine_multipolar_terms<4>(p_lm,g_functions_and_slater_normalization, sphericalHarmonicBuffer);
    default:
        return 0;
    }
}



void HansenCoppens_SF_Engine::pre_hkl_loop_sf_calc(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    std::vector<int> &type_2_wfn_type,
    std::vector<std::vector<REAL> > &def_val_slater_normalization,
    std::vector<int> &typeMaxL)
{
    int nWfnTypes = wfn_parameters.size();
    int nTypes = type_parameters.size();
    int nAtoms = atom_to_wfn_map.size();
    int i,j,nL;

    

    if (mUseIAM)
    {
        typeMaxL.assign(nTypes,-1);
        return;
    }

    
	type_2_wfn_type.resize(nTypes);
    for( int atomIdx = 0 ; atomIdx < nAtoms ; atomIdx++ )
    {
            int atomWfnIdx = atom_to_wfn_map[atomIdx];
            int atomTypeIdx = atom_to_type_map[atomIdx];
            type_2_wfn_type[atomTypeIdx] = atomWfnIdx;
    }
        
    def_val_slater_normalization.resize(nWfnTypes);

    
    
    for(i=0;i<nWfnTypes;i++)
    {
        nL = wfn_parameters[i].def_valence_pow.size();
        def_val_slater_normalization[i].resize(nL);
        for(j=0;j<nL;j++)
            def_val_slater_normalization[i][j] = 
            sto_atomic_wfn::stoDensityNormalizationFactor(wfn_parameters[i].def_valence_pow[j], wfn_parameters[i].def_valence_exp);
    }

    typeMaxL.resize(nTypes);
    int maxL_FromPlm;
    
    for(i=0;i<nTypes;i++)
    {
        maxL_FromPlm = -1;
        for(int l=0;l<type_parameters[i].p_lm.size();l++)
        {
    
            for(j=0;j<2*l+1;j++)
                if(type_parameters[i].p_lm[l][j]!=0.0)
                    maxL_FromPlm = int(l);
        }

        typeMaxL[i] = std::min(4, int(wfn_parameters[type_2_wfn_type[i]].def_valence_pow.size()) - 1);
        typeMaxL[i] = std::min(typeMaxL[i],maxL_FromPlm);
    }
}


void HansenCoppens_SF_Engine::calculateSF_IAM(
    const std::vector<std::string> &atomicType,
    const std::vector<std::complex<REAL> > &atomTypeAnomalousScattering,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_factor,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetryOperations,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &hVectors,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &dTarget_df,
    const std::vector<bool> &include_atom_contribution,
    int nThreads)
{
    mUseIAM = true;
    mIamAtomType = atomicType;
    mAtomToIamTypeMap = atom_to_type_map;

    
    vector<sf_engine_data_types::HC_WfnParam> wfnParams(atomicType.size());
    vector<sf_engine_data_types::HC_TypeParam> typeParams(1);
    vector<int> atomToWfnMap = atom_to_type_map;
    vector<int> atomToTypeMap(atomicPositions.size(),0);
    Matrix3d idenity;
    idenity.setToIdentity();
    vector<Matrix3d> localCoordinateSystems(atomicPositions.size(), idenity);

    int iamTypeIdx,nIamTypes = atomicType.size();

    mIamFormFactors.resize(nIamTypes);
    
    for( iamTypeIdx = 0 ; iamTypeIdx < nIamTypes ; iamTypeIdx++ )
    {
        wfnParams[iamTypeIdx].anomalous_scattering = atomTypeAnomalousScattering[iamTypeIdx];
        if(n_gaussian_form_factors_table::hasFormFactor(mIamAtomType[iamTypeIdx]))
            mIamFormFactors[iamTypeIdx] = n_gaussian_form_factors_table::getFormFactor(mIamAtomType[iamTypeIdx]);
        else
            on_error::throwException(
                string("request for Gaussian type atomic form factor parameter for unknown atom type: ")
                + mIamAtomType[iamTypeIdx], __FILE__, __LINE__);
    }
    
    

    calculateSF(wfnParams,typeParams, atomToWfnMap, atomToTypeMap,atomicPositions,
                atomic_displacement_parameters, atomic_occupancy,atomic_multiplicity_factor,
                localCoordinateSystems, symmetryOperations, centrosymmetric, inversionTranslation, 
                hVectors, f, dTarget_dparam, dTarget_df, include_atom_contribution, nThreads);
    
    
}


void HansenCoppens_SF_Engine::pre_atom_loop_sf_calc(
        //in:
            const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,
            const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams,
            const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
            const std::vector<int> &type_2_wfn,
            const std::vector<std::vector<REAL> > &def_val_slater_normalization,
            const Vector3<REAL> &hVector,
            REAL hVectorLength,
        //out:
            vector<REAL> &wfn_spherical_core_sf,
            vector<REAL> &wfn_spherical_valence_sf,
            vector<vector<REAL> > &g_functions_and_slater_norm,
            vector<Vector3<REAL> > &rotated_h,
            vector<Vector3<REAL> > &rotated_normalized_h,
            std::vector<REAL> &translation_factor,
            std::vector<std::vector<REAL> > &adp_multipliers)
{


    for( int symmOpIdx = 0 ; symmOpIdx< symOps.size() ; symmOpIdx++ )
    {
        translation_factor[symmOpIdx] = hVector*symOps[symmOpIdx].translation;
        rotated_h[symmOpIdx] = hVector*symOps[symmOpIdx].rotation;
        rotated_normalized_h[symmOpIdx] = rotated_h[symmOpIdx]/hVectorLength;

        // sets mAdpMultipliers
        Vector3<REAL> &h = rotated_h[symmOpIdx];
        REAL *adpMultipliers = &adp_multipliers[symmOpIdx][0];
        
        adpMultipliers[0] = h.x*h.x;
        adpMultipliers[1] = h.y*h.y;
        adpMultipliers[2] = h.z*h.z;
        adpMultipliers[3] = 2.0*h.x*h.y;
        adpMultipliers[4] = 2.0*h.x*h.z;
        adpMultipliers[5] = 2.0*h.y*h.z;
    }

    if(mUseIAM)
    {
        
        for (int i = 0, n = wfn_spherical_core_sf.size(); i<n; i++)
            wfn_spherical_core_sf[i] = mIamFormFactors[i].calculate_h(hVectorLength);
        return;
    }
    else
        for( int wfnTypeIdx = 0 ; wfnTypeIdx < wfnParams.size() ; wfnTypeIdx++ )
            wfn_spherical_core_sf[ wfnTypeIdx] = 
                sto_scattering::scatteringSphericalDensity( wfnParams[wfnTypeIdx].core_coeff,
                                                        wfnParams[wfnTypeIdx].core_exp,
                                                        wfnParams[wfnTypeIdx].core_pow,
                                                        hVectorLength);


    int nTypes = typeParams.size();

    for( int typeIdx = 0 ; typeIdx < nTypes ; typeIdx++)
    {
        int wfnTypeIdx = type_2_wfn[ typeIdx ];
        
        wfn_spherical_valence_sf[ typeIdx ] = 
            sto_scattering::scatteringSphericalDensity( wfnParams[wfnTypeIdx].valence_coeff,
                                                        wfnParams[wfnTypeIdx].valence_exp,
                                                        wfnParams[wfnTypeIdx].valence_pow,
                                                        hVectorLength / typeParams[typeIdx].kappa_spherical);
        
        int nL = wfnParams[wfnTypeIdx].def_valence_pow.size();

        const vector<int>& def_valence_pow = wfnParams[wfnTypeIdx].def_valence_pow;
        
        if(nL>0)
            g_functions_and_slater_norm[typeIdx][0] = def_val_slater_normalization[wfnTypeIdx][0] * 
                                                      sto_scattering::gFunction<0>(int(def_valence_pow[0])+2,
                                                                   hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                                   wfnParams[wfnTypeIdx].def_valence_exp);
        if(nL>1)
            g_functions_and_slater_norm[typeIdx][1] = def_val_slater_normalization[wfnTypeIdx][1] * 
                                                      sto_scattering::gFunction<1>(int(def_valence_pow[1])+2,
                                                                   hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                                   wfnParams[wfnTypeIdx].def_valence_exp);

        if(nL>2)
            g_functions_and_slater_norm[typeIdx][2] = def_val_slater_normalization[wfnTypeIdx][2] * 
                                                      sto_scattering::gFunction<2>(int(def_valence_pow[2])+2,
                                                                   hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                                   wfnParams[wfnTypeIdx].def_valence_exp);
        if(nL>3)
            g_functions_and_slater_norm[typeIdx][3] = def_val_slater_normalization[wfnTypeIdx][3] * 
                                                      sto_scattering::gFunction<3>(int(def_valence_pow[3])+2,
                                                                   hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                                   wfnParams[wfnTypeIdx].def_valence_exp);
        if(nL>4)
            g_functions_and_slater_norm[typeIdx][4] = def_val_slater_normalization[wfnTypeIdx][4] * 
                                                      sto_scattering::gFunction<4>(int(def_valence_pow[4])+2,
                                                                   hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                                   wfnParams[wfnTypeIdx].def_valence_exp);

    }

}


void HansenCoppens_SF_Engine::calculateFormFactors(
    const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
    const std::vector<double>& f_spherical, // for each type spherical valence + core
    const std::vector<int>& atom_to_wfn_map,
    const std::vector<int>& atom_to_type_map,
    const std::vector<Matrix3<REAL> >& local_coordinate_systems,
    const Vector3<REAL>& h_vector,
    std::vector<std::complex<REAL> >& form_factors,
    const std::vector<bool>& include_atom,
	const std::vector<int> &type_2_wfn_type,
	const std::vector<std::vector<REAL> > &def_val_slater_normalization,
	const std::vector<int> &typeMaxL)
{
	//--------

	mSphericalHarmonicsData.resize(1);
	mSphericalHarmonicsData[0].resize(5);
	for (int i = 0; i < 5; i++)
		mSphericalHarmonicsData[0][i].resize(2 * i + 1);

	//--------


	REAL hVectorLength;
	
	//hVectorLength2 = h_vector * h_vector;
	hVectorLength = sqrt(h_vector * h_vector);
	Vector3<REAL> normalized_h = h_vector/ hVectorLength;
	

	int atomWfnIdx, atomTypeIdx;
	complex<REAL> atom_f_def_val, aux;

	//--

	vector<vector<REAL> > g_functions_and_slater_norm(type_parameters.size(), vector<REAL>(5));
	
	int nAtoms;
	nAtoms = atom_to_type_map.size();

	form_factors.resize(nAtoms);

	//

	bool hkl000  = (hVectorLength < 1e-10);

	int nTypes = type_parameters.size();

	for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
	{
		int wfnTypeIdx = type_2_wfn_type[typeIdx];

		int nL = wfn_parameters[wfnTypeIdx].def_valence_pow.size();

        const vector<int>& def_valence_pow = wfn_parameters[wfnTypeIdx].def_valence_pow;

		if (nL > 0)
			g_functions_and_slater_norm[typeIdx][0] = def_val_slater_normalization[wfnTypeIdx][0] *
			sto_scattering::gFunction<0>(int(def_valence_pow[0]) + 2,
				hVectorLength / type_parameters[typeIdx].kappa_def_valence,
				wfn_parameters[wfnTypeIdx].def_valence_exp);
		if (nL > 1)
			g_functions_and_slater_norm[typeIdx][1] = def_val_slater_normalization[wfnTypeIdx][1] *
			sto_scattering::gFunction<1>(int(def_valence_pow[1]) + 2,
				hVectorLength / type_parameters[typeIdx].kappa_def_valence,
				wfn_parameters[wfnTypeIdx].def_valence_exp);

		if (nL > 2)
			g_functions_and_slater_norm[typeIdx][2] = def_val_slater_normalization[wfnTypeIdx][2] *
			sto_scattering::gFunction<2>(int(def_valence_pow[2]) + 2,
				hVectorLength / type_parameters[typeIdx].kappa_def_valence,
				wfn_parameters[wfnTypeIdx].def_valence_exp);
		if (nL > 3)
			g_functions_and_slater_norm[typeIdx][3] = def_val_slater_normalization[wfnTypeIdx][3] *
			sto_scattering::gFunction<3>(int(def_valence_pow[3]) + 2,
				hVectorLength / type_parameters[typeIdx].kappa_def_valence,
				wfn_parameters[wfnTypeIdx].def_valence_exp);
		if (nL > 4)
			g_functions_and_slater_norm[typeIdx][4] = def_val_slater_normalization[wfnTypeIdx][4] *
			sto_scattering::gFunction<4>(int(def_valence_pow[4]) + 2,
				hVectorLength / type_parameters[typeIdx].kappa_def_valence,
				wfn_parameters[wfnTypeIdx].def_valence_exp);

	}

	//------------- end of pre_atom_loop_sf_calc

	for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
	{

		if (!include_atom[atomIdx])
		{
			form_factors[atomIdx] = 0;
			continue;
		}

		atomWfnIdx = atom_to_wfn_map[atomIdx];
		atomTypeIdx = atom_to_type_map[atomIdx];

        if (hkl000)
            atom_f_def_val = 0;
        else
            atom_f_def_val = calculateDeformationValence(type_parameters[atomTypeIdx].p_lm,
				g_functions_and_slater_norm[atomTypeIdx],
				local_coordinate_systems[atomIdx],
				normalized_h,
				typeMaxL[atomTypeIdx], mSphericalHarmonicsData[0]);

		form_factors[atomIdx] = atom_f_def_val + f_spherical[atomTypeIdx];

		} 

}

void HansenCoppens_SF_Engine::calculateSphericalTermsInFormFactors(
	const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
	const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
	const std::vector <double> h,
	std::vector< std::vector<REAL> >& f_core,
	std::vector< std::vector<REAL> >& f_sph_valence,
	const std::vector<int>& type_2_wfn_type,
	const std::vector<std::vector<REAL> >& def_val_slater_normalization,
	const std::vector<int>& typeMaxL)
{

	//--
	int nTypes, nWfnTypes;
	nTypes = type_parameters.size();
	nWfnTypes = wfn_parameters.size();
	vector<REAL> wfn_spherical_core_sf(nWfnTypes);
	vector<REAL> wfn_spherical_valence_sf(nTypes);
	vector<vector<REAL> > g_functions_and_slater_norm(nTypes, vector<REAL>(5));

	int nH = h.size();
	

	f_core.resize(nWfnTypes,vector<double>(nH));
	f_sph_valence.resize(nTypes,vector<double>(nH));

	for (int hIndex = 0; hIndex < nH; hIndex++)
	{

		for (int wfnTypeIdx = 0; wfnTypeIdx < nWfnTypes; wfnTypeIdx++)
			//wfn_spherical_core_sf[wfnTypeIdx] =
			f_core[wfnTypeIdx][hIndex]=
			sto_scattering::scatteringSphericalDensity(wfn_parameters[wfnTypeIdx].core_coeff,
				wfn_parameters[wfnTypeIdx].core_exp,
				wfn_parameters[wfnTypeIdx].core_pow,
				h[hIndex]);

		for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
		{
			int wfnTypeIdx = type_2_wfn_type[typeIdx];

			//wfn_spherical_valence_sf[typeIdx] =
			f_sph_valence[typeIdx][hIndex] = 
				sto_scattering::scatteringSphericalDensity(wfn_parameters[wfnTypeIdx].valence_coeff,
					wfn_parameters[wfnTypeIdx].valence_exp,
					wfn_parameters[wfnTypeIdx].valence_pow,
					h[hIndex] / type_parameters[typeIdx].kappa_spherical);

			//wfn_spherical_valence_sf[typeIdx] *= type_parameters[typeIdx].p_val;
			f_sph_valence[typeIdx][hIndex] *= type_parameters[typeIdx].p_val;
		}

	}


}


//void HansenCoppens_SF_Engine::calculateFormFactors(
//	const std::vector<sf_engine_data_types::HC_WfnParam>& wfn_parameters,
//	const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
//	const std::vector<int>& atom_to_wfn_map,
//	const std::vector<int>& atom_to_type_map,
//	const std::vector<Matrix3<REAL> >& local_coordinate_systems,
//	const std::vector<Vector3<REAL> >& h_vectors,
//	std::vector< std::vector<std::complex<REAL> > >& f)
//{
//	// allocate memory buffers for spherical harmonics calculations
//
//	mSphericalHarmonicsData.resize(1);
//	int maxNL, nL, typeIdx, nTypes = type_parameters.size();
//	maxNL = 0;
//	for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
//	{
//		nL = type_parameters[typeIdx].p_lm.size();
//		if (nL > maxNL)
//			maxNL = nL;
//	}
//
//	mSphericalHarmonicsData[0].resize(maxNL);
//	for (int l = 0; l < maxNL; l++)
//		mSphericalHarmonicsData[0][l].resize(2 * l + 1);
//	
//	//
//
//	vector<vector<REAL> > def_val_slater_normalization;
//	vector<int> type_2_wfn(type_parameters.size(), 100);
//	vector<int> type_max_L;
//	pre_hkl_loop_sf_calc(wfnParams, typeParams, atom_to_wfn_map, atom_to_type_map,
//		type_2_wfn, def_val_slater_normalization, type_max_L);
//
//
//}


void HansenCoppens_SF_Engine::calculateSF(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,
    const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomicPositions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_factor,
    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
    bool centrosymmetric,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &hVectors,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &dTarget_df,
    const std::vector<bool> &include_atom_contribution,
    int nThreads)
{
    int threadIdx;
    WallClockTimer timer;
    
    if( nThreads<1 )
        on_error::throwException(string("wrong number of cores/thread specified: ")+string_utilities::convertToString(nThreads),
                                 __FILE__,__LINE__);
    
#ifndef _OPENMP   
    nThreads = 1;
#endif
    // allocate memory buffers for spherical harmonics calculations

    mSphericalHarmonicsData.resize(nThreads);
    int maxNL,nL,typeIdx,nTypes = typeParams.size();
    maxNL = 0;
    for(typeIdx=0;typeIdx<nTypes;typeIdx++)
    {
        nL = typeParams[typeIdx].p_lm.size();
        if(nL>maxNL)
            maxNL = nL;
    }

    for(threadIdx=0;threadIdx<nThreads;threadIdx++)
    {
        mSphericalHarmonicsData[threadIdx].resize(maxNL);
        for(int l=0;l<maxNL;l++)
            mSphericalHarmonicsData[threadIdx][l].resize(2*l+1);
    }
    
    //

    vector<vector<REAL> > def_val_slater_normalization;
    vector<int> type_2_wfn(typeParams.size(), 100);
    vector<int> type_max_L;
    pre_hkl_loop_sf_calc(wfnParams, typeParams, atom_to_wfn_map, atom_to_type_map,
                         type_2_wfn, def_val_slater_normalization, type_max_L);

    //

    int nHklVectors = hVectors.size();
    int atomIdx, nAtoms = atom_to_wfn_map.size();

    // set structure factors and derivatives to zero

    f.assign(nHklVectors,0.0);

    dTarget_dparam.resize(nAtoms);

    for (atomIdx = 0; atomIdx < nAtoms; ++atomIdx)
    {
        dTarget_dparam[atomIdx].adp_derivatives.assign(atomic_displacement_parameters[atomIdx].size(), 0.0);
        dTarget_dparam[atomIdx].atomic_position_derivatives = Vector3d(0, 0, 0);
        dTarget_dparam[atomIdx].occupancy_derivatives = 0.0;
    }


    // declare arrays for storing partial results for each thread

    vector< vector< complex< REAL > > > perThreadSF;
    vector< vector< Vector3< REAL > > > perThreadHklVectors;
    vector< vector< complex< REAL > > > perThreadTarget_dF;
    vector<vector<TargetFunctionAtomicParamDerivatives> > perThread_dTarget_dParam;


    // assign HKL vectors to threads

    perThreadSF.resize(nThreads,f);
    perThreadHklVectors.resize(nThreads);
    perThread_dTarget_dParam.resize(nThreads, dTarget_dparam);
    perThreadTarget_dF.resize(nThreads);


    
    scattering_utilities::divideSet(hVectors, nThreads, perThreadHklVectors);
    scattering_utilities::divideSet(dTarget_df, nThreads, perThreadTarget_dF);
    
    /*
    if(nHklVectors % nThreads == 0)
        maxN_HklVectorsPerThread = nHklVectors / nThreads;
    else
        maxN_HklVectorsPerThread = nHklVectors / nThreads + 1;
    
    for( threadIdx=0 ; threadIdx<nThreads; threadIdx++ )
        if(threadIdx != nThreads - 1 )
        {
            perThreadHklVectors[threadIdx].insert(perThreadHklVectors[threadIdx].end(),
                                                  hVectors.begin() + threadIdx * maxN_HklVectorsPerThread,
                                                  hVectors.begin() + (threadIdx + 1) * maxN_HklVectorsPerThread);
            perThreadTarget_dF[threadIdx].insert( perThreadTarget_dF[threadIdx].end(),
                                                  dTarget_df.begin() + threadIdx * maxN_HklVectorsPerThread,
                                                  dTarget_df.begin() + (threadIdx + 1) * maxN_HklVectorsPerThread);
        }
        else
        {
            perThreadHklVectors[threadIdx].insert(perThreadHklVectors[threadIdx].end(),
                                                  hVectors.begin() + threadIdx * maxN_HklVectorsPerThread,
                                                  hVectors.end());

            perThreadTarget_dF[threadIdx].insert( perThreadTarget_dF[threadIdx].end(),
                                                  dTarget_df.begin() + threadIdx * maxN_HklVectorsPerThread,
                                                  dTarget_df.end());
        }
     */
    // set number of threads 
#if defined(_OPENMP)
    omp_set_dynamic(0);
    omp_set_num_threads(int(nThreads));
#endif

    // run parallel calculations 
    timer.start();
#pragma omp parallel
    {
        int threadId;
#if defined(_OPENMP)
        threadId = omp_get_thread_num();
#else
        threadId = 0;
#endif

        if(centrosymmetric)
        {
            if(inversionTranslation == Vector3<REAL>(0, 0, 0))
                calculateSF_SerialCentrosymmetric(
                    wfnParams, typeParams, atom_to_wfn_map, atom_to_type_map,
                    atomicPositions, atomic_displacement_parameters, atomic_occupancy,
                    atomic_multiplicity_factor, local_coordinate_systems, symOps,
                    perThreadHklVectors[threadId], perThreadSF[threadId],
                    perThread_dTarget_dParam[threadId], perThreadTarget_dF[threadId], include_atom_contribution,
                    type_2_wfn, def_val_slater_normalization, mSphericalHarmonicsData[threadId], type_max_L);
            else
                calculateSF_SerialSymmetryCenterNotAtOrigin(
                    wfnParams, typeParams, atom_to_wfn_map, atom_to_type_map,
                    atomicPositions, atomic_displacement_parameters, atomic_occupancy,
                    atomic_multiplicity_factor, local_coordinate_systems, symOps, inversionTranslation,
                    perThreadHklVectors[threadId], perThreadSF[threadId],
                    perThread_dTarget_dParam[threadId], perThreadTarget_dF[threadId], include_atom_contribution,
                    type_2_wfn, def_val_slater_normalization, mSphericalHarmonicsData[threadId], type_max_L);

        }
        else
            calculateSF_SerialAcentric(
                wfnParams, typeParams, atom_to_wfn_map, atom_to_type_map,
                atomicPositions, atomic_displacement_parameters, atomic_occupancy,
                atomic_multiplicity_factor, local_coordinate_systems, symOps,
                perThreadHklVectors[threadId], perThreadSF[threadId],
                perThread_dTarget_dParam[threadId], perThreadTarget_dF[threadId], include_atom_contribution,
                type_2_wfn, def_val_slater_normalization, mSphericalHarmonicsData[threadId], type_max_L);

    }

    double t = timer.stop();

    
    scattering_utilities::combineScatteringFactorsSets(perThreadSF, f);
    scattering_utilities::merge_dTarget_dParameterSets(perThread_dTarget_dParam, dTarget_dparam);
    
   

} //calculateSF_parallel_2




//---------------

void HansenCoppens_SF_Engine::calculateSF_SerialAcentric(
    const std::vector<sf_engine_data_types::HC_WfnParam> &_wfnParams,
    const std::vector<sf_engine_data_types::HC_TypeParam> &_typeParams,
    const std::vector<int> &_atom_to_wfn_map,
    const std::vector<int> &_atom_to_type_map,
    const std::vector<Vector3<REAL> > &_atomic_positions,
    const std::vector<std::vector<REAL> > &_atomic_displacement_parameters,
    const std::vector<REAL> &_atomic_occupancy,
    const std::vector<REAL> &_atomic_multiplicity_weight,
    const std::vector<Matrix3<REAL> > &_local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &_symOps,
    const std::vector<Vector3<REAL> > &_hVectors,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &_dTarget_df,
    const std::vector<bool> &_include_atom_contribution,
    const std::vector<int> &_type_2_wfn_type,
    const std::vector<std::vector<REAL> > &_def_val_slater_normalization,
    std::vector<std::vector<double> > &_sphericalHarmonicsData,
    const std::vector<int> &_typeMaxL)
{
    // local copies

    vector<sf_engine_data_types::HC_WfnParam> wfnParams = _wfnParams;
    vector<sf_engine_data_types::HC_TypeParam> typeParams = _typeParams;
    vector<int> atom_to_wfn_map = _atom_to_wfn_map;
    vector<int> atom_to_type_map = _atom_to_type_map;
    vector<Vector3<REAL> > atomic_positions = _atomic_positions;
    vector<std::vector<REAL> > atomic_displacement_parameters = _atomic_displacement_parameters;
    vector<REAL> atomic_occupancy = _atomic_occupancy;
    vector<REAL> atomic_multiplicity_weight = _atomic_multiplicity_weight;
    vector<Matrix3<REAL> > local_coordinate_systems = _local_coordinate_systems;
    vector<sf_engine_data_types::SymmetryOperation> symOps = _symOps;
    vector<Vector3<REAL> > hVectors = _hVectors;
    vector<std::complex<REAL> > dTarget_df = _dTarget_df;
    vector<bool> include_atom_contribution = _include_atom_contribution;
    vector<int> type_2_wfn_type = _type_2_wfn_type;
    vector<std::vector<REAL> > def_val_slater_normalization = _def_val_slater_normalization;
    vector<std::vector<double> > sphericalHarmonicsData = _sphericalHarmonicsData;
    vector<int> typeMaxL = _typeMaxL;

    //--


    const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
    complex<REAL> fAtomSphericalAndAnomalous;
    const REAL two_pi = 2 * REAL(M_PI);
    const REAL two_pi_sqare = 2 * REAL(M_PI*M_PI);
    REAL temperature_factor, hVectorLength, hVectorLength2, multiplier;
    REAL const *adps;
    vector<Vector3<REAL> > rotated_h(symOps.size());
    vector<Vector3<REAL> > rotated_normalized_h(symOps.size());
    vector<REAL> translation_factor(symOps.size());
    complex<REAL> unweightedTransformedAtomFF, unweightedTransformedAtomFF_Sum, dTargetDf;

    complex<REAL> adp_derivatives[6];
    complex<REAL> xyz_derivatives[3];
    REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;
    REAL realFContrib, imagFContrib;

    int atomWfnIdx, atomTypeIdx;
    REAL atomWeight; // = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];
    REAL atom_f_core, atom_f_sph_val; 
    complex<REAL> anomalousScattering, atom_f_def_val, aux;

    //--

    vector<REAL> wfn_spherical_core_sf(wfnParams.size());
    vector<REAL> wfn_spherical_valence_sf(typeParams.size());
    vector<vector<REAL> > g_functions_and_slater_norm(typeParams.size(), vector<REAL>(5));
    int nSymmOps = symOps.size();
    int n_adp_components, nAtoms, nHklVectors = hVectors.size();
    vector<vector<REAL> > adpMultipliers(nSymmOps,vector<double>(6));
    nAtoms = atomic_positions.size();
    

    // set dTarget_dparam to zero ..

    f.assign(nHklVectors, 0.0);

    dTarget_dparam.resize(atomic_positions.size());

    for (int atom_index = 0; atom_index < atomic_positions.size(); atom_index++)
    {
        dTarget_dparam[atom_index].adp_derivatives.assign(atomic_displacement_parameters[atom_index].size(), 0.0);
        dTarget_dparam[atom_index].atomic_position_derivatives = Vector3d(0, 0, 0);
        dTarget_dparam[atom_index].occupancy_derivatives = 0.0;
    }

    //

    bool hkl000;

    for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++) 
    {
        hVectorLength2 = hVectors[hklIndex] * hVectors[hklIndex];
        hVectorLength = sqrt(hVectorLength2);

        hkl000 = (hVectorLength < 1e-10);

        dTargetDf = dTarget_df[hklIndex];

        pre_atom_loop_sf_calc(
            //in:
            wfnParams, typeParams, symOps, type_2_wfn_type, def_val_slater_normalization, hVectors[hklIndex], hVectorLength,
            //out:
            wfn_spherical_core_sf, wfn_spherical_valence_sf, g_functions_and_slater_norm,
            rotated_h, rotated_normalized_h, translation_factor, adpMultipliers);

        realFContrib = 0;
        imagFContrib = 0;


        for (int atomIdx = 0; atomIdx <nAtoms; atomIdx++)
        {

            if (!include_atom_contribution[atomIdx])
                continue;

            atomWfnIdx = atom_to_wfn_map[atomIdx];
            atomTypeIdx = atom_to_type_map[atomIdx];

            atomWeight = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];
            atom_f_core = wfn_spherical_core_sf[atomWfnIdx];
            atom_f_sph_val = wfn_spherical_valence_sf[atomTypeIdx];
            atom_f_sph_val *= typeParams[atomTypeIdx].p_val;
            anomalousScattering = wfnParams[atomWfnIdx].anomalous_scattering;
            fAtomSphericalAndAnomalous = atom_f_core + atom_f_sph_val + anomalousScattering;

            n_adp_components = atomic_displacement_parameters[atomIdx].size();

            if(!atomic_displacement_parameters[atomIdx].empty())
                adps = &atomic_displacement_parameters[atomIdx][0];

            if( n_adp_components == 6 )
                for (int i = 0; i< 6 ; i++)
                    adp_derivatives[i] = 0.0;

            xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

            unweightedTransformedAtomFF_Sum = 0.0;

            for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
            {

                const Vector3<REAL> & rotated_h_ref = rotated_h[symOpIdx];

                atomic_phase_factor_phase = two_pi*(rotated_h_ref*atomic_positions[atomIdx] +
                                                    translation_factor[symOpIdx]);

                atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);

                if (hkl000)
                    atom_f_def_val = 0;
                else
                    atom_f_def_val =  calculateDeformationValence(typeParams[atomTypeIdx].p_lm,
                                                                 g_functions_and_slater_norm[atomTypeIdx],
                                                                 local_coordinate_systems[atomIdx],
                                                                 rotated_normalized_h[symOpIdx],
                                                                 typeMaxL[atomTypeIdx],sphericalHarmonicsData);

                if(n_adp_components==6)
                {
                    double *multipliers = &adpMultipliers[symOpIdx][0];
                    temperature_factor = exp(- multipliers[0] * adps[0] - multipliers[1] * adps[1]
                                             - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                                             - multipliers[4] * adps[4] - multipliers[5] * adps[5]);

                    unweightedTransformedAtomFF = (fAtomSphericalAndAnomalous + atom_f_def_val)
                                                   * atomic_position_phase_factor * temperature_factor;

                    adp_derivatives[0] -= multipliers[0] * unweightedTransformedAtomFF;
                    adp_derivatives[1] -= multipliers[1] * unweightedTransformedAtomFF;
                    adp_derivatives[2] -= multipliers[2] * unweightedTransformedAtomFF;
                    adp_derivatives[3] -= multipliers[3] * unweightedTransformedAtomFF;
                    adp_derivatives[4] -= multipliers[4] * unweightedTransformedAtomFF;
                    adp_derivatives[5] -= multipliers[5] * unweightedTransformedAtomFF;

                }
                else
                    unweightedTransformedAtomFF = (fAtomSphericalAndAnomalous + atom_f_def_val)
                                                  * atomic_position_phase_factor;
                
                unweightedTransformedAtomFF_Sum += unweightedTransformedAtomFF;
                
                xyz_derivatives[0] += rotated_h_ref[0] * unweightedTransformedAtomFF;
                xyz_derivatives[1] += rotated_h_ref[1] * unweightedTransformedAtomFF;
                xyz_derivatives[2] += rotated_h_ref[2] * unweightedTransformedAtomFF;


            } // symmetry operations
            
            
            if(n_adp_components == 1)
            {
                temperature_factor = exp(-hVectorLength * hVectorLength * (*adps));
                unweightedTransformedAtomFF_Sum *= temperature_factor;
                dTarget_dparam[atomIdx].adp_derivatives[0] -= hVectorLength2 * (dTargetDf.real() * unweightedTransformedAtomFF_Sum.real() -
                                                                                dTargetDf.imag() * unweightedTransformedAtomFF_Sum.imag());
            }
            else 
                if(n_adp_components == 6)
                    for (int i = 0; i<6; ++i)
                        dTarget_dparam[atomIdx].adp_derivatives[i] += dTargetDf.real()*adp_derivatives[i].real()-
                                                                      dTargetDf.imag()*adp_derivatives[i].imag();

            
            dTarget_dparam[atomIdx].occupancy_derivatives += unweightedTransformedAtomFF_Sum.real()*dTargetDf.real() -
                                                             unweightedTransformedAtomFF_Sum.imag()*dTargetDf.imag();

            realFContrib += unweightedTransformedAtomFF_Sum.real() * atomWeight;
            imagFContrib += unweightedTransformedAtomFF_Sum.imag() * atomWeight;

            n_adp_components == 1 ? aux = temperature_factor * dTargetDf : aux = dTargetDf;

            for (int i = 0; i<3; ++i)
                dTarget_dparam[atomIdx].atomic_position_derivatives[i] -= aux.real()*xyz_derivatives[i].imag()+
                                                                          aux.imag()*xyz_derivatives[i].real();
 
       } // symmetrically independent atoms
        
        f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);

    } // h vectors

    for(int atomIdx=0; atomIdx<nAtoms ; atomIdx++)
    {
        atomWeight = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];
        dTarget_dparam[atomIdx].occupancy_derivatives *= atomic_multiplicity_weight[atomIdx];
        multiplier = two_pi_sqare * atomWeight;
        for(int i=0;i<dTarget_dparam[atomIdx].adp_derivatives.size();i++)
            dTarget_dparam[atomIdx].adp_derivatives[i] *= multiplier;
        multiplier = two_pi * atomWeight;
        dTarget_dparam[atomIdx].atomic_position_derivatives.x *= multiplier;
        dTarget_dparam[atomIdx].atomic_position_derivatives.y *= multiplier;
        dTarget_dparam[atomIdx].atomic_position_derivatives.z *= multiplier;
    }

} // calculateSF_SerialAcentric




void HansenCoppens_SF_Engine::calculateSF_SerialCentrosymmetric(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,
    const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomic_positions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_weight,
    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
    const std::vector<Vector3<REAL> > &hVectors,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &dTarget_df,
    const std::vector<bool> &include_atom_contribution,
    const std::vector<int> &type_2_wfn_type,
    const std::vector<std::vector<REAL> > &def_val_slater_normalization,
    std::vector<std::vector<double> > &sphericalHarmonicsData,
    const std::vector<int> &typeMaxL)
{
    const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
    complex<REAL> fAtomSphericalAndAnomalous,f0,f_dispersion;
    const REAL two_pi = 2 * REAL(M_PI);
    const REAL two_pi_sqare = 2 * REAL(M_PI*M_PI);
    vector<Vector3<REAL> > rotated_h(symOps.size());
    vector<Vector3<REAL> > rotated_normalized_h(symOps.size());
    vector<REAL> translation_factor(symOps.size());

    vector<REAL> wfn_spherical_core_sf(wfnParams.size());
    vector<REAL> wfn_spherical_valence_sf(typeParams.size());
    vector<vector<REAL> > g_functions_and_slater_norm(typeParams.size(), vector<REAL>(5));
    int nSymmOps = symOps.size();
    int nHklVectors = hVectors.size();
    vector<vector<REAL> > adpMultipliers(nSymmOps,vector<REAL>(6));
    REAL sumTeperatureRealPartPhase;
    complex<REAL> transformedAtomContribWithoutWeight, xyzDerivativesMultiplier;
    // set dTarget_dparam to zero ..

    f.assign(nHklVectors, 0.0);

    dTarget_dparam.resize(atomic_positions.size());

    for (int atom_index = 0; atom_index < atomic_positions.size(); atom_index++)
    {
        dTarget_dparam[atom_index].adp_derivatives.assign(atomic_displacement_parameters[atom_index].size(), 0.0);
        dTarget_dparam[atom_index].atomic_position_derivatives = Vector3d(0, 0, 0);
        dTarget_dparam[atom_index].occupancy_derivatives = 0.0;
    }

    //
    bool hkl000;

    for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
    {
        REAL hVectorLength = sqrt(hVectors[hklIndex] * hVectors[hklIndex]);
        hkl000 = (hVectorLength < 1e-10);

        pre_atom_loop_sf_calc(
            //in:
            wfnParams, typeParams, symOps, type_2_wfn_type, def_val_slater_normalization, hVectors[hklIndex], hVectorLength,
            //out:
            wfn_spherical_core_sf, wfn_spherical_valence_sf, g_functions_and_slater_norm,
            rotated_h, rotated_normalized_h, translation_factor, adpMultipliers);



        REAL realFContrib = 0;
        REAL imagFContrib = 0;

        

        for (int atomIdx = 0; atomIdx < atomic_positions.size(); atomIdx++)
        {
            complex<REAL> adp_derivatives[6];
            complex<REAL> xyz_derivatives[3];
            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;

            if (!include_atom_contribution[atomIdx])
                continue;
            complex<REAL> atomicFContribWithoutWeight = 0;

            int atomWfnIdx = atom_to_wfn_map[atomIdx];
            int atomTypeIdx = atom_to_type_map[atomIdx];




            REAL atomWeight = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];


            REAL atom_f_core = wfn_spherical_core_sf[atomWfnIdx];
            REAL atom_f_sph_val = wfn_spherical_valence_sf[atomTypeIdx];
            atom_f_sph_val *= typeParams[atomTypeIdx].p_val;

            complex<REAL> anomalousScattering = wfnParams[atomWfnIdx].anomalous_scattering;
            fAtomSphericalAndAnomalous = atom_f_core + atom_f_sph_val + anomalousScattering;

            // end of h direction independent part of atomic scattering factor calculation

            // prepare derivatives data

            int n_adp_components = atomic_displacement_parameters[atomIdx].size();
            for (int i = 0; i<n_adp_components; i++)
                adp_derivatives[i] = 0.0;
            xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

            // loop over symetry operations

            sumTeperatureRealPartPhase = 0;

            for (int symOpIdx = 0; symOpIdx < symOps.size(); symOpIdx++)
            {

                const Vector3<REAL> & rotated_h_ref = rotated_h[symOpIdx];

                atomic_phase_factor_phase = two_pi*(rotated_h_ref*atomic_positions[atomIdx] +
                                                    translation_factor[symOpIdx]);

                atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);

                complex<REAL> atom_f_def_val;



                if (hkl000)
                    atom_f_def_val = 0;
                else
                    atom_f_def_val = calculateDeformationValence(typeParams[atomTypeIdx].p_lm,
                        g_functions_and_slater_norm[atomTypeIdx],
                        local_coordinate_systems[atomIdx],
                        rotated_normalized_h[symOpIdx],
                        typeMaxL[atomTypeIdx], sphericalHarmonicsData);


                // temperature factor
                REAL temperature_factor;
                atomic_displacement_parameters[atomIdx].empty() ?
                    temperature_factor = 1.0 :
                    temperature_factor = calc_temperature_factor(rotated_h_ref, hVectorLength, atomic_displacement_parameters[atomIdx]);


                complex<REAL> transformedAtomF = 2 * temperature_factor * (fAtomSphericalAndAnomalous*atomic_phase_factor_real +
                                          atom_f_def_val.real() * atomic_phase_factor_real -
                                          atom_f_def_val.imag() * atomic_phase_factor_im);


                atomicFContribWithoutWeight += transformedAtomF;

                // add contribution of dF[h]/dAtomicParameter to dTargetFunction/dParameter
         

                xyzDerivativesMultiplier = ((fAtomSphericalAndAnomalous + atom_f_def_val)*temperature_factor*atomic_position_phase_factor - 
                                           (fAtomSphericalAndAnomalous + conj(atom_f_def_val))*temperature_factor*conj(atomic_position_phase_factor))*
                                           two_pi_i;

                xyz_derivatives[0] += rotated_h_ref[0] * xyzDerivativesMultiplier;
                xyz_derivatives[1] += rotated_h_ref[1] * xyzDerivativesMultiplier;
                xyz_derivatives[2] += rotated_h_ref[2] * xyzDerivativesMultiplier;

                //process_adp_derivatives(adp_derivatives, transformedAtomF, rotated_h_ref, hVectorLength, n_adp_components);

                if(n_adp_components==6)
                {
                    adp_derivatives[0] -= rotated_h_ref[0] * rotated_h_ref[0] * transformedAtomF;
                    adp_derivatives[1] -= rotated_h_ref[1] * rotated_h_ref[1] * transformedAtomF;
                    adp_derivatives[2] -= rotated_h_ref[2] * rotated_h_ref[2] * transformedAtomF;
                    adp_derivatives[3] -= 2 * rotated_h_ref[0] * rotated_h_ref[1] * transformedAtomF;
                    adp_derivatives[4] -= 2 * rotated_h_ref[0] * rotated_h_ref[2] * transformedAtomF;
                    adp_derivatives[5] -= 2 * rotated_h_ref[1] * rotated_h_ref[2] * transformedAtomF;
                }
                


            } // symmetry operations

            if (include_atom_contribution[atomIdx])
            {
                dTarget_dparam[ atomIdx ].occupancy_derivatives += 
                    ( atomicFContribWithoutWeight * dTarget_df[ hklIndex ] ).real() * 
                    atomic_multiplicity_weight[ atomIdx ];

                realFContrib += atomicFContribWithoutWeight.real()*atomWeight;
                imagFContrib += atomicFContribWithoutWeight.imag()*atomWeight;

                complex<REAL> aux(dTarget_df[hklIndex] * atomWeight);
                // adp
                if (n_adp_components > 0)
                {
                    if (n_adp_components == 1)
                        dTarget_dparam[atomIdx].adp_derivatives[0] += -hVectorLength * hVectorLength * two_pi_sqare *
                                                                        (atomicFContribWithoutWeight * aux ).real();
                    else
                        for (int i = 0; i<6; ++i)
                            dTarget_dparam[atomIdx].adp_derivatives[i] += two_pi_sqare*(aux*adp_derivatives[i]).real();
                }
                // xyz

                for (int i = 0; i<3; ++i)
                    dTarget_dparam[atomIdx].atomic_position_derivatives[i] += (aux*xyz_derivatives[i]).real();
            }

        } // symetrically independent atoms
        f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);

    } // h vectors


} // calculateSF_SerialCentrosymmetric




void HansenCoppens_SF_Engine::calculateSF_SerialSymmetryCenterNotAtOrigin(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,
    const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL> > &atomic_positions,
    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<REAL> &atomic_multiplicity_weight,
    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
    const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL> > &hVectors,
    std::vector<std::complex<REAL> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL> > &dTarget_df,
    const std::vector<bool> &include_atom_contribution,
    const std::vector<int> &type_2_wfn_type,
    const std::vector<std::vector<REAL> > &def_val_slater_normalization,
    std::vector<std::vector<double> > &sphericalHarmonicsData,
    const std::vector<int> &typeMaxL)
{
    const complex<REAL> two_pi_i(0, 2 * REAL(M_PI));
    complex<REAL> fAtomSphericalAndAnomalous, f0, f_dispersion, term1, term2;
    const REAL two_pi = 2 * REAL(M_PI);
    const REAL two_pi_sqare = 2 * REAL(M_PI*M_PI);
    vector<Vector3<REAL> > rotated_h(symOps.size());
    vector<Vector3<REAL> > rotated_normalized_h(symOps.size());
    vector<REAL> translation_factor(symOps.size());

    vector<REAL> wfn_spherical_core_sf(wfnParams.size());
    vector<REAL> wfn_spherical_valence_sf(typeParams.size());
    vector<vector<REAL> > g_functions_and_slater_norm(typeParams.size(), vector<REAL>(5));
    int nSymmOps = symOps.size();
    vector<vector<REAL> > adpMultipliers(nSymmOps,vector<REAL>(6));
    int nHklVectors = hVectors.size();
    REAL sumTeperatureRealPartPhase;
    complex<REAL> transformedAtomContribWithoutWeight, xyzDerivativesMultiplier, inversionTranslationPhaseFactor;
    // set dTarget_dparam to zero ..

    //DEBUG

    vector<vector<complex<double> > > sf(atomic_positions.size(),vector<complex<double> >(nSymmOps));

    //END DEBUG

    f.assign(nHklVectors, 0.0);

    dTarget_dparam.resize(atomic_positions.size());

    for (int atom_index = 0; atom_index < atomic_positions.size(); atom_index++)
    {
        dTarget_dparam[atom_index].adp_derivatives.assign(atomic_displacement_parameters[atom_index].size(), 0.0);
        dTarget_dparam[atom_index].atomic_position_derivatives = Vector3d(0, 0, 0);
        dTarget_dparam[atom_index].occupancy_derivatives = 0.0;
    }

    //

    bool hkl000;
    for (int hklIndex = 0; hklIndex < nHklVectors; hklIndex++)
    {
        REAL hVectorLength = sqrt(hVectors[hklIndex] * hVectors[hklIndex]);
        hkl000 = (hVectorLength < 1e-10);

        pre_atom_loop_sf_calc(
            //in:
            wfnParams, typeParams, symOps, type_2_wfn_type, def_val_slater_normalization, hVectors[hklIndex], hVectorLength,
            //out:
            wfn_spherical_core_sf, wfn_spherical_valence_sf, g_functions_and_slater_norm,
            rotated_h, rotated_normalized_h, translation_factor, adpMultipliers);

        inversionTranslationPhaseFactor = exp(two_pi_i*(hVectors[hklIndex]*inversionTranslation));

        REAL realFContrib = 0;
        REAL imagFContrib = 0;


        for (int atomIdx = 0; atomIdx < atomic_positions.size(); atomIdx++)
        {
            complex<REAL> adp_derivatives[6];
            complex<REAL> xyz_derivatives[3];
            REAL atomic_phase_factor_real, atomic_phase_factor_im, atomic_phase_factor_phase;

            if (!include_atom_contribution[atomIdx])
                continue;
            complex<REAL> atomicFContribWithoutWeight = 0;

            int atomWfnIdx = atom_to_wfn_map[atomIdx];
            int atomTypeIdx = atom_to_type_map[atomIdx];




            REAL atomWeight = atomic_occupancy[atomIdx] * atomic_multiplicity_weight[atomIdx];


            REAL atom_f_core = wfn_spherical_core_sf[atomWfnIdx];
            REAL atom_f_sph_val = wfn_spherical_valence_sf[atomTypeIdx];
            atom_f_sph_val *= typeParams[atomTypeIdx].p_val;

            complex<REAL> anomalousScattering = wfnParams[atomWfnIdx].anomalous_scattering;
            fAtomSphericalAndAnomalous = atom_f_core + atom_f_sph_val + anomalousScattering;

            // end of h direction independent part of atomic scattering factor calculation

            // prepare derivatives data

            int n_adp_components = atomic_displacement_parameters[atomIdx].size();
            for (int i = 0; i<n_adp_components; i++)
                adp_derivatives[i] = 0.0;
            xyz_derivatives[0] = xyz_derivatives[1] = xyz_derivatives[2] = 0;

            // loop over symmetry operations

            sumTeperatureRealPartPhase = 0;

            for (int symOpIdx = 0; symOpIdx < symOps.size(); symOpIdx++)
            {

                const Vector3<REAL> & rotated_h_ref = rotated_h[symOpIdx];

                atomic_phase_factor_phase = two_pi*(rotated_h_ref*atomic_positions[atomIdx] +
                                                    translation_factor[symOpIdx]);

                atomic_phase_factor_real = cos(atomic_phase_factor_phase);
                atomic_phase_factor_im = sin(atomic_phase_factor_phase);

                complex<REAL> atomic_position_phase_factor(atomic_phase_factor_real, atomic_phase_factor_im);


                complex<REAL> atom_f_def_val = 
                    hkl000 ? 0.0 : 
                             calculateDeformationValence(typeParams[atomTypeIdx].p_lm,
                                                         g_functions_and_slater_norm[atomTypeIdx],
                                                         local_coordinate_systems[atomIdx],
                                                         rotated_normalized_h[symOpIdx],
                                                         typeMaxL[atomTypeIdx], sphericalHarmonicsData);

                // temperature factor
                REAL temperature_factor;
                atomic_displacement_parameters[atomIdx].empty() ?
                    temperature_factor = 1.0 :
                    temperature_factor = calc_temperature_factor(rotated_h_ref, hVectorLength, atomic_displacement_parameters[atomIdx]);

                term1 = (fAtomSphericalAndAnomalous + atom_f_def_val) * atomic_position_phase_factor;
                term2 = (fAtomSphericalAndAnomalous + conj(atom_f_def_val)) * conj(atomic_position_phase_factor) *
                        inversionTranslationPhaseFactor;

                complex<REAL> transformedAtomF = temperature_factor * ( term1 + term2 );

                atomicFContribWithoutWeight += transformedAtomF;

                if(hklIndex==5)
                    sf[atomIdx][symOpIdx] = transformedAtomF;

                xyzDerivativesMultiplier = temperature_factor*(term1 - term2)*two_pi_i;

                xyz_derivatives[0] += rotated_h_ref[0] * xyzDerivativesMultiplier;
                xyz_derivatives[1] += rotated_h_ref[1] * xyzDerivativesMultiplier;
                xyz_derivatives[2] += rotated_h_ref[2] * xyzDerivativesMultiplier;

                if (n_adp_components == 6)
                {
                    adp_derivatives[0] -= rotated_h_ref[0] * rotated_h_ref[0] * transformedAtomF;
                    adp_derivatives[1] -= rotated_h_ref[1] * rotated_h_ref[1] * transformedAtomF;
                    adp_derivatives[2] -= rotated_h_ref[2] * rotated_h_ref[2] * transformedAtomF;
                    adp_derivatives[3] -= 2 * rotated_h_ref[0] * rotated_h_ref[1] * transformedAtomF;
                    adp_derivatives[4] -= 2 * rotated_h_ref[0] * rotated_h_ref[2] * transformedAtomF;
                    adp_derivatives[5] -= 2 * rotated_h_ref[1] * rotated_h_ref[2] * transformedAtomF;
                }

            } // symmetry operations

            if (include_atom_contribution[atomIdx])
            {
                dTarget_dparam[atomIdx].occupancy_derivatives +=
                    (atomicFContribWithoutWeight * dTarget_df[hklIndex]).real() *
                    atomic_multiplicity_weight[atomIdx];

                realFContrib += atomicFContribWithoutWeight.real()*atomWeight;
                imagFContrib += atomicFContribWithoutWeight.imag()*atomWeight;

                complex<REAL> aux(dTarget_df[hklIndex] * atomWeight);
                // adp
                if (n_adp_components > 0)
                {
                    if (n_adp_components == 1)
                        //                        df_dparam = -hVectorLength * hVectorLength * atomic_f;
                        //                  adp_derivatives[0] += dTarget_dF * df_dparam;
                        //dTarget_dparam[atomIdx].adp_derivatives[0] += -hVectorLength * hVectorLength * two_pi_sqare *
                          //                                             (atomicFContribWithoutWeight * aux * dTarget_df[hklIndex]).real();
                        dTarget_dparam[atomIdx].adp_derivatives[0] += -hVectorLength * hVectorLength * two_pi_sqare *
                                                                     (atomicFContribWithoutWeight * aux ).real();
                    else
                        for (int i = 0; i<6; ++i)
                            dTarget_dparam[atomIdx].adp_derivatives[i] += two_pi_sqare*(aux*adp_derivatives[i]).real();
                }
                // xyz
                //aux *= 2.0 * two_pi_i;
                for (int i = 0; i<3; ++i)
                    dTarget_dparam[atomIdx].atomic_position_derivatives[i] += (aux*xyz_derivatives[i]).real();
            }

        } // symetrically independent atoms
          //exit(0);
        f[hklIndex] = complex<REAL>(realFContrib, imagFContrib);

    } // h vectors


} // calculateSF_SerialSymmetryCenterNotAtOrigin

} // namespace discamb
