#include "discamb/Scattering/HansenCoppens_SF_Engine4.h"

#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/Timer.h"

#include "discamb/HC_Model/HC_WfnData.h"

#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/MathUtilities/SphConverter.h"

#include "discamb/Scattering/NGaussianFormFactorsTable.h"
#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/Scattering/SlaterTypeOrbitalScattering.h"


#include <cmath>
#include <algorithm>

#include <cassert>
#include <iomanip>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <array>
#include <iostream>
#include <ctime>


using namespace std;

namespace discamb {

    HansenCoppens_SF_Engine4::HansenCoppens_SF_Engine4()
    {
        mUseIAM = false;
    }

    HansenCoppens_SF_Engine4::~HansenCoppens_SF_Engine4()
    {
    }



    inline void HansenCoppens_SF_Engine4::add_contribution_to_occupancy_derivative(
        REAL &occupancy_derivative,
        const complex<REAL> &dTarget_dF,
        const complex<REAL> &atomic_f_divided_by_occupancy)
    {
        occupancy_derivative += (dTarget_dF * atomic_f_divided_by_occupancy).real();
    }

    inline void HansenCoppens_SF_Engine4::add_contribution_to_position_derivatives(
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


    inline void HansenCoppens_SF_Engine4::add_contribution_to_adp_derivatives(
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

    inline void HansenCoppens_SF_Engine4::process_adp_derivatives( std::complex<REAL> *pre_derivatives,
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

    std::complex<double> HansenCoppens_SF_Engine4::calculateDeformationValence(
        const std::vector<std::vector<REAL> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
    const std::vector<REAL>& g_functions_and_slater_normalization,
    //const Matrix3<REAL>& local_coordinates_system,
    int max_l,
    std::vector<std::vector<double> >& sphericalHarmonics)
    {
        if (max_l < 0)
            return 0;

        switch (max_l)
        {
            case 0:
                return combine_multipolar_terms<0>(p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
            case 1:
                return combine_multipolar_terms<1>(p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
            case 2:
                return combine_multipolar_terms<2>(p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
            case 3:
                return combine_multipolar_terms<3>(p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
            case 4:
                return combine_multipolar_terms<4>(p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
            default:
                return 0;
        }

    }


    std::complex<REAL> HansenCoppens_SF_Engine4::calculateDeformationValence(
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



    void HansenCoppens_SF_Engine4::pre_hkl_loop_sf_calc(
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


    void HansenCoppens_SF_Engine4::calculateSF_IAM(
        const UnitCell& unitCell,
        const std::vector<std::string> &atomicType,
        const std::vector<std::complex<REAL> > &atomTypeAnomalousScattering,
        const std::vector<int> &atom_to_type_map,
        const std::vector<Vector3<REAL> > &atomicPositions,
        const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
        const std::vector<REAL> &atomic_occupancy,
        const std::vector<std::complex<REAL> >& anomalous_dispersion,
        const std::vector<REAL> &atomic_multiplicity_factor,
        const std::vector<sf_engine_data_types::SymmetryOperation> &symmetryOperations,
        bool centrosymmetric,
        const Vector3<REAL> &inversionTranslation,
        const std::vector<Vector3<REAL> > &hVectors,
        const std::vector<Vector3i >& hkl_indices,
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


        DerivativesSelector derivativesSwitch;
        calculateSF(unitCell, wfnParams, typeParams, atomToWfnMap, atomToTypeMap, atomicPositions,
                    atomic_displacement_parameters, atomic_occupancy, anomalous_dispersion, atomic_multiplicity_factor,
                    localCoordinateSystems, symmetryOperations, centrosymmetric, inversionTranslation,
                    hVectors, hkl_indices, f, dTarget_dparam, dTarget_df, include_atom_contribution, nThreads, derivativesSwitch);


    }


    void HansenCoppens_SF_Engine4::pre_atom_loop_sf_calc(
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

        if (mUseIAM)
        {

            for (int i = 0, n = wfn_spherical_core_sf.size(); i < n; i++)
                wfn_spherical_core_sf[i] = mIamFormFactors[i].calculate_h(hVectorLength);
            return;
        }
        else
            for (int wfnTypeIdx = 0; wfnTypeIdx < wfnParams.size(); wfnTypeIdx++)
                wfn_spherical_core_sf[wfnTypeIdx] =
                sto_scattering::scatteringSphericalDensity( wfnParams[wfnTypeIdx].core_coeff,
                                                            wfnParams[wfnTypeIdx].core_exp,
                                                            wfnParams[wfnTypeIdx].core_pow,
                                                            hVectorLength);


                int nTypes = typeParams.size();

            for( int typeIdx = 0 ; typeIdx < nTypes ; typeIdx++)
            {
                int wfnTypeIdx = type_2_wfn[ typeIdx ];

                wfn_spherical_valence_sf[typeIdx] =
                sto_scattering::scatteringSphericalDensity( wfnParams[wfnTypeIdx].valence_coeff,
                                                            wfnParams[wfnTypeIdx].valence_exp,
                                                            wfnParams[wfnTypeIdx].valence_pow,
                                                            hVectorLength / typeParams[typeIdx].kappa_spherical);

                int nL = wfnParams[wfnTypeIdx].def_valence_pow.size();

                const vector<int>& def_valence_pow = wfnParams[wfnTypeIdx].def_valence_pow;

                if (nL > 0)
                    g_functions_and_slater_norm[typeIdx][0] =  def_val_slater_normalization[wfnTypeIdx][0] *
                    sto_scattering::gFunction<0>(int(def_valence_pow[0])+2,
                                                 hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                 wfnParams[wfnTypeIdx].def_valence_exp);
                    if (nL > 1)
                        g_functions_and_slater_norm[typeIdx][1] =  def_val_slater_normalization[wfnTypeIdx][1] *
                        sto_scattering::gFunction<1>(int(def_valence_pow[1])+2,
                                                     hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                     wfnParams[wfnTypeIdx].def_valence_exp);

                        if (nL > 2)
                            g_functions_and_slater_norm[typeIdx][2] = def_val_slater_normalization[wfnTypeIdx][2] *
                            sto_scattering::gFunction<2>(int(def_valence_pow[2])+2,
                                                         hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                         wfnParams[wfnTypeIdx].def_valence_exp);
                            if (nL > 3)
                                g_functions_and_slater_norm[typeIdx][3] = def_val_slater_normalization[wfnTypeIdx][3] *
                                sto_scattering::gFunction<3>(int(def_valence_pow[3])+2,
                                                             hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                             wfnParams[wfnTypeIdx].def_valence_exp);
                                if (nL > 4)
                                    g_functions_and_slater_norm[typeIdx][4] = def_val_slater_normalization[wfnTypeIdx][4] *
                                    sto_scattering::gFunction<4>(int(def_valence_pow[4])+2,
                                                                 hVectorLength / typeParams[typeIdx].kappa_def_valence,
                                                                 wfnParams[wfnTypeIdx].def_valence_exp);

            }

    }


    void HansenCoppens_SF_Engine4::calculateFormFactors(
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

    void HansenCoppens_SF_Engine4::calculateSphericalTermsInFormFactors(
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



    void HansenCoppens_SF_Engine4::calculateGlobalCoordinatesPlm(
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
        const std::vector<int>& atom_to_type_map,
        const std::vector<Matrix3<REAL> >& local_coordinate_systems,// rows are vectors
        std::vector< std::vector<std::vector<double> > > & atomPlms)
    {
        int maxL = 4;
        vector<vector<double> > den2wfn;
        real_spherical_harmonics::getDensityToWfnMultipliers(maxL, den2wfn);

        int nTypes = type_parameters.size();
        vector<vector<vector<double> > > typePlmWfn(nTypes);
        for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {
            int typeMaxL = type_parameters[typeIdx].p_lm.size() - 1;
            typePlmWfn[typeIdx] = type_parameters[typeIdx].p_lm;
            for (int l = 0; l <= typeMaxL; l++)
                for (int i = 0; i < 2 * l + 1; i++)
                {
                    int abs_m = abs(l - i);
                    typePlmWfn[typeIdx][l][i] *= den2wfn[l][abs_m];
                }

        }

        SphConverter sphConverter;
        vector<vector<vector<double> > > conversionMatrices;
        sphConverter.setMaxL(maxL);

        vector<vector<double> > localCoordinates(3, vector<double>(3));
        vector<vector<double> > cartesianCoordinates{ {1.0,0.0,0.0}, {0.0,1.0,0.0},{0.0,0.0,1.0} };


        int nAtoms = atom_to_type_map.size();
        atomPlms.resize(nAtoms);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    localCoordinates[i][j] = local_coordinate_systems[atomIdx](i, j);

            sphConverter.convert(localCoordinates, cartesianCoordinates, conversionMatrices);
            int atomType = atom_to_type_map[atomIdx];
            int typeMaxL = type_parameters[atomType].p_lm.size() - 1;
            atomPlms[atomIdx].resize(typeMaxL + 1);
            for (int l = 0; l <= typeMaxL; l++)
            {
                atomPlms[atomIdx][l].resize(2 * l + 1);
                for (int i = 0; i < 2 * l + 1; i++)
                {
                    atomPlms[atomIdx][l][i] = 0.0;
                    for (int j = 0; j < 2 * l + 1; j++)
                        atomPlms[atomIdx][l][i] += conversionMatrices[l][i][j] * typePlmWfn[atomType][l][j];

                    int abs_m = abs(l - i);
                    atomPlms[atomIdx][l][i] /= den2wfn[l][abs_m];
                }
            }
        }
    }
    /*
     * void calculateSF(
     *    const UnitCell &unitCell,
     *    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
     *    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
     *    const std::vector<int> &atom_to_wfn_map,
     *    const std::vector<int> &atom_to_type_map,
     *    const std::vector<Vector3<REAL> > &atomicPositions,
     *    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
     *    const std::vector<REAL> &atomic_occupancy,
     *    const std::vector<REAL> &atomic_multiplicity_factor,
     *    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
     *    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
     *    bool centrosymmetric,
     *    const Vector3<REAL> &inversionTranslation,
     *    const std::vector<Vector3<REAL> > &h_vectors,
     *    const std::vector<Vector3i >& hkl_indices,
     *    std::vector<std::complex<REAL> > &f,
     *    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
     *    const std::vector<std::complex<REAL> > &dTarget_df,
     *    const std::vector<bool> &include_atom_contribution,
     *    int nThreads);
     */

    inline void printStep(std::string name){
        std::cout << "calculateSF: after - " << name << " - time = " << std::clock() << std::endl;
    }

    inline void printInLoop(std::string name){
        //std::cout << "calculateSF-main: after - " << name << " - time = " << std::clock() << std::endl;
    }

    template<typename T>
    inline void infrequentValueLog(std::string name, T value){
        //std::cout << name << " = " << value << std::endl;
    }

    template<typename T>
    inline void valueLog(std::string name, T value){
        //std::cout << name << " = " << value << std::endl;
    }

    constexpr std::array<int, 3> binSize = {8, 8, 8};
    //constexpr bool virtLinePhaseFlag = true; // maybe add later (the ability to change to false)
    //constexpr bool virtLineTemperatureFlag = true; // maybe add later (the ability to change to false)
    constexpr bool virtHklPhaseFlag = false;
    constexpr bool virtHklTemperatureFlag = false;
    constexpr bool virtHklTMulFlag = true;
    //constexpr bool symOpExtractionFlag = true; // maybe add later (the ability to change to false)
    //constexpr bool fSymDeduplicationFlag = false; // cctbx already dedupes // maybe add later (the ability to change to true)
    //constexpr bool versorDeduplicationFlag = true; // maybe add later (the ability to change to true) since it should only speed it up ~1.2 times
    //constexpr bool lengthDeduplicationFlag = true; // maybe add later (the ability to change to true)
    //constexpr bool symOpOffsetDeduplication = false; // maybe add later (the ability to change to true)

    constexpr double two_pi = 2.0*M_PI;
    constexpr double two_pi_squared = two_pi*M_PI;
    constexpr double four_pi_squared = two_pi*two_pi;
    constexpr double four_pi = 4.0*M_PI;

    inline REAL square(REAL x){
        return x*x;
    }

    inline bool closeToZero(REAL x){
        return -1e-12<x and x<1e-12;
    }

    inline bool closeToZero(Vector3<REAL> vec){
        return (
            closeToZero(vec[0]) and
            closeToZero(vec[1]) and
            closeToZero(vec[2]));
    }

    inline bool closeToZero(Matrix3<REAL> mat){
        return (
            closeToZero({mat(0,0), mat(0,1), mat(0,2)}) and
            closeToZero({mat(1,0), mat(1,1), mat(1,2)}) and
            closeToZero({mat(2,0), mat(2,1), mat(2,2)}));
    }

    inline Matrix3<REAL> U(const std::vector<REAL> &adps){
        assert(adps.size()==6);
        return Matrix3<REAL>(
            adps[0], adps[3], adps[4],
            adps[3], adps[1], adps[5],
            adps[4], adps[5], adps[2]);
    }

    inline int gcd(int a, int b){
        while (b != 0){
            int t = b;
            b = a % b;
            a = t;
        }
        return a;
    }

    inline REAL sqrt(REAL x){
        return std::sqrt(x); // TODO make it use a more precise approach (some new c++ versions may not have support for double here)
    }

    inline REAL pow(REAL x, int n){
        return std::pow(x, n); // TODO make it use a more precise approach (some new c++ versions may not have support for double here)
    }
    /* this breaks the results by ~10^-5 (tyrosine) TODO make a more robust replacement (for now use the compiler default)
     * REAL cos(REAL x){
     *    return std::cosf(std::fmodf(x, two_pi));
}

REAL sin(REAL x){
return std::sinf(std::fmodf(x, two_pi));
}
*/
    void HansenCoppens_SF_Engine4::calculateSF(
        const UnitCell &unitCell,
        const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams, // per wfn
        const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams, // per type
        const std::vector<int> &atom_to_wfn_map,
        const std::vector<int> &atom_to_type_map,
        const std::vector<Vector3<REAL> > &atomicPositions, // per atom
        const std::vector<std::vector<REAL> > &atomic_displacement_parameters, // per atom and already premultiplied by two_pi_squared
        const std::vector<REAL> &atomic_occupancy,
        const std::vector<std::complex<REAL> >& anomalous_dispersion, // per atom // TODO use
        const std::vector<REAL> &atomic_multiplicity_factor,
        const std::vector<Matrix3<REAL> > &local_coordinate_systems, // per atom
        const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
        bool centrosymmetric, // TODO use
        const Vector3<REAL> &inversionTranslation, // TODO use
        const std::vector<Vector3<REAL> > &hVectors,
        const std::vector<Vector3i >& hkl_indices, // unused
        std::vector<std::complex<REAL> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam, // per atom // TODO generate
        const std::vector<std::complex<REAL> > &dTarget_df, // per hkl // TODO use
        const std::vector<bool> &include_atom_contribution, // per atom
        int nThreads,
        const DerivativesSelector& derivativesSwitch, // TODO use
        bool electron, // TODO use
        const std::vector<int>& atomic_numbers) // TODO use
    {
        printStep("calculateSF start");
        const int trueNAtoms = atom_to_wfn_map.size();
        std::vector<int> usedAtomIndices;
        usedAtomIndices.clear();
        for (int i=0; i<trueNAtoms; i++){
            if (include_atom_contribution[i])
                usedAtomIndices.emplace_back(i);
        }
        const int hklCount = hVectors.size();
        const int nAtoms = usedAtomIndices.size();

        infrequentValueLog("nAtoms", nAtoms);
        for (int atom=0; atom<nAtoms; atom++){
            valueLog("usedAtomIndices[atom]", usedAtomIndices[atom]);
            for (int i=0; i<atomic_displacement_parameters[usedAtomIndices[atom]].size(); i++){
                valueLog("atomic_displacement_parameters[usedAtomIndices[atom]]",atomic_displacement_parameters[usedAtomIndices[atom]][i]);
            }
        }

        const int nSymOps = symOps.size();

        std::vector<Matrix3<REAL>> symOpMults;
        std::vector<int> symOpToMult;
        symOpToMult.resize(nSymOps);
        symOpMults.emplace_back(symOps[0].rotation);
        valueLog("symOpIdx", 0);
        valueLog("current(0,0)", symOps[0].rotation(0,0));
        valueLog("current(0,1)", symOps[0].rotation(0,1));
        valueLog("current(0,2)", symOps[0].rotation(0,2));
        valueLog("current(1,0)", symOps[0].rotation(1,0));
        valueLog("current(1,1)", symOps[0].rotation(1,1));
        valueLog("current(1,2)", symOps[0].rotation(1,2));
        valueLog("current(2,0)", symOps[0].rotation(2,0));
        valueLog("current(2,1)", symOps[0].rotation(2,1));
        valueLog("current(2,2)", symOps[0].rotation(2,2));
        for (int symOpIdx = 1; symOpIdx<nSymOps; symOpIdx++){
            auto &current = symOps[symOpIdx].rotation;
            valueLog("symOpIdx", symOpIdx);
            valueLog("current(0,0)", current(0,0));
            valueLog("current(0,1)", current(0,1));
            valueLog("current(0,2)", current(0,2));
            valueLog("current(1,0)", current(1,0));
            valueLog("current(1,1)", current(1,1));
            valueLog("current(1,2)", current(1,2));
            valueLog("current(2,0)", current(2,0));
            valueLog("current(2,1)", current(2,1));
            valueLog("current(2,2)", current(2,2));
            bool isInMults = false;
            for (int multIdx = 0; multIdx<symOpMults.size(); multIdx++){
                if (closeToZero(current - symOpMults[multIdx])){
                    isInMults = true;
                    symOpToMult[symOpIdx] = multIdx;
                    break;
                }
            }
            if (not isInMults){
                symOpToMult[symOpIdx] = symOpMults.size();
                symOpMults.emplace_back(current);
            }
        }
        const int symOpMultCount = symOpMults.size();
        infrequentValueLog("symOpMultCount", symOpMultCount);

        std::vector<Vector3<REAL>> symOpOffsets;
        std::vector<int> symOpToOffset;
        symOpToOffset.resize(nSymOps);
        symOpOffsets.emplace_back(symOps[0].translation);
        valueLog("symOpIdx", 0);
        valueLog("current[0]", symOps[0].translation[0]);
        valueLog("current[1]", symOps[0].translation[1]);
        valueLog("current[2]", symOps[0].translation[2]);
        for (int symOpIdx = 1; symOpIdx<nSymOps; symOpIdx++){
            auto &current = symOps[symOpIdx].translation;
            valueLog("symOpIdx", symOpIdx);
            valueLog("current[0]", current[0]);
            valueLog("current[1]", current[1]);
            valueLog("current[2]", current[2]);
            bool isInOffsets = false;
            for (int offsetIdx = 0; offsetIdx<symOpOffsets.size(); offsetIdx++){
                if (closeToZero(current - symOpOffsets[offsetIdx])){
                    isInOffsets = true;
                    symOpToOffset[symOpIdx] = offsetIdx;
                    break;
                }
            }
            if (not isInOffsets){
                symOpToOffset[symOpIdx] = symOpOffsets.size();
                symOpOffsets.emplace_back(current);
            }
        }
        const int symOpOffsetCount = symOpOffsets.size();
        infrequentValueLog("symOpOffsetCount", symOpOffsetCount);

        printStep("symOp deduplication");

        Vector3<int> minHkl = hkl_indices[0];
        Vector3<int> maxHkl = hkl_indices[0];
        Vector3<int> step = {0, 0, 0};
        //Vector3<REAL> minH = hVectors[0];
        //Vector3<REAL> maxH = hVectors[0];
        for (int hklIdx = 1; hklIdx<hklCount; hklIdx++){
            const Vector3<int> &currentHkl = hkl_indices[hklIdx];
            const Vector3<REAL> &currentH = hVectors[hklIdx];
            for (int i = 0; i<3; i++){
                if (currentHkl[i]<minHkl[i]) {
                    minHkl[i] = currentHkl[i];
                    //minH[i] = currentH[i];
                }
                if (currentHkl[i]>maxHkl[i]) {
                    maxHkl[i] = currentHkl[i];
                    //maxH[i] = currentH[i];
                }
                const int diff = std::abs(currentHkl[i]-hkl_indices[0][i]);
                if (step[i]==0)
                    step[i]=diff;
                else if ((diff % step[i]) != 0) {
                    step[i] = gcd(diff, step[i]);
                };
            }
        }

        for (int i=0; i<3; i++){
            if (step[i]==0)
                step[i]=1;
        }

        infrequentValueLog("step[0]", step[0]);
        infrequentValueLog("step[1]", step[1]);
        infrequentValueLog("step[2]", step[2]);

        /*
         *    Vector3<REAL> scale;
         *    for (int i =0; i<3; i++){
         *        if (((maxHkl[i] - minHkl[i]) == 0) or closeToZero(maxH[i] - minH[i])){
         *          if ((maxHkl[i] == 0) or closeToZero(maxH[i])){
         *            scale[i]=1.0;
    } else {
        scale[i]=maxH[i]/maxHkl[i];
    }
    }
    else{
        scale[i]=(maxH[i] - minH[i]) / (maxHkl[i] - minHkl[i]);
    }
    }

    infrequentValueLog("scale[0]", scale[0]);
    infrequentValueLog("scale[1]", scale[1]);
    infrequentValueLog("scale[2]", scale[2]);
    */

        Vector3<int> nHkls;
        for (int i = 0; i<3; i++){
            const int diff = maxHkl[i]-minHkl[i];
            assert(diff % step[i] == 0);
            nHkls[i]=diff/step[i] + 1;
        }

        Vector3<int> nBins;
        for (int i = 0; i<3; i++){
            if ((nHkls[i]%binSize[i]) == 0)
                nBins[i] = nHkls[i] / binSize[i];
            else
                nBins[i] = (nHkls[i] / binSize[i]) + 1;
        }

        const int totalNBins = nBins[0]*nBins[1]*nBins[2];


        printStep("hkl binning");

        // [h*nk*nl + k*nl + l]
        std::vector<bool> isBinUsed;
        isBinUsed.resize(totalNBins, false);
        for (int hklIdx = 0; hklIdx<hklCount; hklIdx++){
            const Vector3<int> &currentHkl = hkl_indices[hklIdx];
            Vector3<int> binIdx;
            for (int i=0; i<3; i++){
                //assert((currentHkl[i] % step[i]) == 0)
                binIdx[i] = (currentHkl[i] - minHkl[i]) / (step[i] * binSize[i]);
            }
            isBinUsed[(binIdx[0]*nBins[1] + binIdx[1])*nBins[2] + binIdx[2]] = true;
        }
        std::vector<int> allBinsMap;
        allBinsMap.resize(nBins[0]*nBins[1]*nBins[2]);
        std::vector<Vector3<int>> usedBins;
        usedBins.clear();
        for (int h = 0; h<nBins[0]; h++){
            for (int k = 0; k<nBins[1]; k++){
                for (int l = 0; l<nBins[2]; l++){
                    if (isBinUsed[(h*nBins[1] + k)*nBins[2] + l]){
                        allBinsMap[(h*nBins[1] + k)*nBins[2] + l] = usedBins.size();
                        usedBins.emplace_back(h, k, l);
                    }
                }
            }
        }

        const int binCount = usedBins.size();

        infrequentValueLog("binCount", binCount);

        printStep("bin deduplication");

        vector<vector<Vector3d> > r_atom_rot(nAtoms,vector<Vector3<REAL>>(symOpMultCount));
        #pragma omp parallel for num_threads(nThreads) collapse(2)
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++){
            for (int symOpMultIdx = 0; symOpMultIdx < symOpMultCount; symOpMultIdx++){
                r_atom_rot[atomIdx][symOpMultIdx] = symOpMults[symOpMultIdx]*atomicPositions[usedAtomIndices[atomIdx]];
            }
        }

        ReciprocalLatticeUnitCell recUnitCell(unitCell);
        std::array<Vector3<REAL>, 3> stepCartesian;
        recUnitCell.fractionalToCartesian({step[0], 0.0, 0.0}, stepCartesian[0]);
        recUnitCell.fractionalToCartesian({0.0, step[1], 0.0}, stepCartesian[1]);
        recUnitCell.fractionalToCartesian({0.0, 0.0, step[2]}, stepCartesian[2]);

        std::vector<Vector3<REAL>> bin000Cartesian;
        bin000Cartesian.resize(binCount);
        #pragma omp parallel for num_threads(nThreads)
        for (int binIdx = 0; binIdx<binCount; binIdx++){
            Vector3<int> bin000;
            for (int i=0; i<3; i++){
                bin000[i] = usedBins[binIdx][i] * step[i] * binSize[i] + minHkl[i];
                infrequentValueLog("bin000[i]", bin000[i]);
            }
            recUnitCell.fractionalToCartesian(bin000, bin000Cartesian[binIdx]);
            infrequentValueLog("bin000Cartesian[binIdx][0]", bin000Cartesian[binIdx][0]);
            infrequentValueLog("bin000Cartesian[binIdx][1]", bin000Cartesian[binIdx][1]);
            infrequentValueLog("bin000Cartesian[binIdx][2]", bin000Cartesian[binIdx][2]);
        }

        printStep("fractional to cartesian conversion");

        // [atom][symOpMult | 0]
        std::vector<std::vector<Matrix3<REAL>>> rotatedUs;
        rotatedUs.resize(nAtoms);
        for (int atom=0; atom<nAtoms; atom++){
            rotatedUs[atom].resize((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 6) ? symOpMultCount : 0);
        }
        #pragma omp parallel for num_threads(nThreads)
        for (int atom=0; atom<nAtoms; atom++){
            if (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 6) {
                /*auto ftcMatrixT = recUnitCell.getFractionalToCartesianMatrix();
                 *            ftcMatrixT.transpose();
                 *            const auto orig = recUnitCell.getFractionalToCartesianMatrix() * U(atomic_displacement_parameters[usedAtomIndices[atom]]) * ftcMatrixT;
                 */
                valueLog("atom", atom);
                for (int i=0; i<6; i++){
                    valueLog("atomic_displacement_parameters[usedAtomIndices[atom]][i]", atomic_displacement_parameters[usedAtomIndices[atom]][i]);
                }
                const auto orig = U(atomic_displacement_parameters[usedAtomIndices[atom]]);
                for (int symOpIdx = 0; symOpIdx<symOpMultCount; symOpIdx++){
                    Matrix3<REAL> symOpT;
                    for (int i=0; i<3; i++){
                        for(int j=0; j<3; j++){
                            symOpT(i, j) = symOpMults[symOpIdx](j, i);
                        }
                    }
                    rotatedUs[atom][symOpIdx] = symOpT * orig * symOpMults[symOpIdx];
                    /*
                     *                // trying to do M = symOp * orig * symOp^T
                     *                // because orig^T = orig and M^T = M
                     *                // it is enough to compute only 6/9 of acc and 6/9 of M = P * symOp^T
                     *                  std::array<Matrix3<REAL>, 3> acc;
                     *                  for (int x=0; x<3; x++){
                     *                      for (int y=0; y<3; y++){
                     *                          for (int i=y; i<3; i++){
                     *                              acc[i](x, y) = symOp(i, y) * orig(x, i);
                }
                }
                }
                Matrix3<REAL> P;
                for (int x=0; x<3; x++){
                    for (int y=0; y<3; y++){
                        for (int i=0; i<y; i++){
                            P(x, y) += acc[y](x, i);
                }
                for (int i=y; i<3; i++){
                    P(x, y) += acc[i](x, y);
                }
                }
                }
                Matrix3<REAL> M;
                for (int x=0; x<3; x++){
                    for (int y=0; y<=x; y++){
                        for (int i=0; i<3; i++){
                            M(x, y) += P(i, y) * symOp(i, x);
                }
                }
                }
                rotatedUs[atom][symOpIdx] = Matrix3<REAL>(
                    M(0, 0), M(1, 0), M(2, 0),
                    M(1, 0), M(1, 1), M(2, 1),
                    M(2, 0), M(2, 1), M(2, 2));
                    */
                }
            }
        }
        printStep("rotating Us");


        // [bin][atom][symOpMult | 1]
        std::vector<std::vector<std::vector<REAL>>> temperatureFactorRoots;
        temperatureFactorRoots.resize(binCount);
        for (int binIdx = 0; binIdx<binCount; binIdx++){
            temperatureFactorRoots[binIdx].resize(nAtoms);
            for (int atom = 0; atom<nAtoms; atom++){
                temperatureFactorRoots[binIdx][atom].resize((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount);
            }
        }
        #pragma omp parallel for num_threads(nThreads) collapse(2)
        for (int binIdx = 0; binIdx<binCount; binIdx++){
            for (int atom = 0; atom<nAtoms; atom++){
                bool iso = (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1);
                if (iso) {
                    REAL T1iso = std::exp(-atomic_displacement_parameters[usedAtomIndices[atom]][0] * (
                        square(bin000Cartesian[binIdx][0]) +
                        square(bin000Cartesian[binIdx][1]) +
                        square(bin000Cartesian[binIdx][2])) );
                    temperatureFactorRoots[binIdx][atom][0] = T1iso;
                } else {
                    for (int symOpIdx = 0; symOpIdx<symOpMultCount; symOpIdx++){
                        Matrix3<REAL> currentU = rotatedUs[atom][symOpIdx];
                        REAL T1 = std::exp(-(
                            square(bin000Cartesian[binIdx][0])*currentU(0,0) +
                            square(bin000Cartesian[binIdx][1])*currentU(1,1) +
                            square(bin000Cartesian[binIdx][2])*currentU(2,2) + 2.0*(
                                bin000Cartesian[binIdx][0]*bin000Cartesian[binIdx][1]*currentU(0,1) +
                                bin000Cartesian[binIdx][0]*bin000Cartesian[binIdx][2]*currentU(0,2) +
                                bin000Cartesian[binIdx][1]*bin000Cartesian[binIdx][2]*currentU(1,2))));
                        temperatureFactorRoots[binIdx][atom][symOpIdx] = T1;
                    }
                }
            }
        }

        printStep("temperature factor roots");

        // [bin][atom][symOp]
        std::vector<std::vector<std::vector<std::complex<REAL>>>> phaseFactorRoots;
        phaseFactorRoots.resize(binCount);
        for (int binIdx = 0; binIdx<binCount; binIdx++){
            phaseFactorRoots[binIdx].resize(nAtoms);
            for (int atom = 0; atom<nAtoms; atom++){
                phaseFactorRoots[binIdx][atom].resize(nSymOps);
            }
        }
        #pragma omp parallel for num_threads(nThreads) collapse(3)
        for (int binIdx = 0; binIdx<binCount; binIdx++){
            for (int atom = 0; atom<nAtoms; atom++){
                for (int symOpIdx = 0; symOpIdx<nSymOps; symOpIdx++){
                    const REAL phase_angle_root = two_pi * (r_atom_rot[atom][symOpToMult[symOpIdx]] + symOpOffsets[symOpToOffset[symOpIdx]]) * bin000Cartesian[binIdx];
                    const std::complex<REAL> result = { cos(phase_angle_root), sin(phase_angle_root) };
                    phaseFactorRoots[binIdx][atom][symOpIdx] = result;
                }
            }
        }

        printStep("phase factor roots");

        // [dir][atom][symOp][n]
        std::array<std::vector<std::vector<std::vector<std::complex<REAL>>>>, 3> phaseFactorMults;
        for (int i = 0; i<3; i++){
            phaseFactorMults[i].resize(nAtoms);
            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                phaseFactorMults[i][atomIdx].resize(nSymOps);
                for (int symOpIdx = 0; symOpIdx<nSymOps; symOpIdx++){
                    phaseFactorMults[i][atomIdx][symOpIdx].resize(binSize[i]);
                }
            }
        }
        #pragma omp parallel for num_threads(nThreads) collapse(3)
        for (int i = 0; i < 3; i++){
            for (int atom = 0; atom<nAtoms; atom++){
                for (int symOpIdx = 0; symOpIdx<nSymOps; symOpIdx++){
                    const REAL phase_angle_mult = two_pi * (r_atom_rot[atom][symOpToMult[symOpIdx]] + symOpOffsets[symOpToOffset[symOpIdx]]) * stepCartesian[i];
                    const std::complex<REAL> single = { cos(phase_angle_mult), sin(phase_angle_mult) };
                    std::complex<REAL> acc = 1.0;
                    for (int j = 0; j<binSize[i]; j++){
                        phaseFactorMults[i][atom][symOpIdx][j] = acc;
                        acc*=single;
                    }
                }
            }
        }

        printStep("phase factor multipliers");

        // [h*nk*nl + k*nl + l][atom][symOp]
        std::vector<std::vector<std::vector<std::complex<REAL>>>> virtHklPhase;
        if constexpr (virtHklPhaseFlag) {
            virtHklPhase.resize(binSize[0]*binSize[1]*binSize[2]);
            for (int hklIdx = 0; hklIdx<(binSize[0]*binSize[1]*binSize[2]); hklIdx++){
                virtHklPhase[hklIdx].resize(nAtoms);
                for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                    virtHklPhase[hklIdx][atomIdx].resize(nSymOps);
                }
            }
            #pragma omp parallel for num_threads(nThreads) collapse(4)
            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                for (int h = 0; h<binSize[0]; h++){
                    for (int k = 0; k<binSize[1]; k++){
                        for (int l = 0; l<binSize[2]; l++){
                            for (int symOpIdx = 0; symOpIdx<nSymOps; symOpIdx++){
                                virtHklPhase[(h*binSize[1] + k)*binSize[2] + l][atomIdx][symOpIdx] =
                                phaseFactorMults[0][atomIdx][symOpIdx][h] *
                                phaseFactorMults[1][atomIdx][symOpIdx][k] *
                                phaseFactorMults[2][atomIdx][symOpIdx][l];
                            }
                        }
                    }
                }
            }
        }

        printStep("virtual hkl phase factors");

        // [dir][atom][symOpMult | 1][n]
        std::array<std::vector<std::vector<std::vector<REAL>>>, 3> temperatureFactorMults;
        for (int i=0; i<3; i++){
            temperatureFactorMults[i].resize(nAtoms);
            for (int atom=0; atom<nAtoms; atom++){
                temperatureFactorMults[i][atom].resize((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount);
                for (int symOp=0; symOp<((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount); symOp++){
                    temperatureFactorMults[i][atom][symOp].resize(binSize[i]);
                }
            }
        }
        #pragma omp parallel for num_threads(nThreads) collapse(2)
        for (int i=0; i<3; i++){
            for (int atom=0; atom<nAtoms; atom++){
                bool iso = (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1);
                if (iso) {
                    REAL c = std::exp(
                        -atomic_displacement_parameters[usedAtomIndices[atom]][0] * (
                            square(stepCartesian[i][0]) +
                            square(stepCartesian[i][1]) +
                            square(stepCartesian[i][2])));
                    REAL mult = 1.0;
                    REAL acc = c;
                    REAL c_pow = square(c);
                    temperatureFactorMults[i][atom][0][0] = 1.0;
                    for (int j=1; j<binSize[i]; j++){
                        //mult*=pow(c, j*2 - 1);
                        mult*=acc;
                        temperatureFactorMults[i][atom][0][j] = mult;
                        acc*=c_pow;
                    }
                } else {
                    for (int symOp=0; symOp<symOpMultCount; symOp++){
                        Matrix3<REAL> currentU = rotatedUs[atom][symOp];
                        REAL c = std::exp(-(
                            square(stepCartesian[i][0])*currentU(0, 0) +
                            square(stepCartesian[i][1])*currentU(1, 1) +
                            square(stepCartesian[i][2])*currentU(2, 2) + 2.0*(
                                stepCartesian[i][0]*stepCartesian[i][1]*currentU(0, 1) +
                                stepCartesian[i][0]*stepCartesian[i][2]*currentU(0, 2) +
                                stepCartesian[i][1]*stepCartesian[i][2]*currentU(1, 2))));
                        REAL mult = 1.0;
                        REAL acc = c;
                        REAL c_pow = square(c);
                        temperatureFactorMults[i][atom][symOp][0] = 1.0;
                        for (int j=1; j<binSize[i]; j++){
                            //mult*=pow(c, j*2 - 1);
                            mult*=acc;
                            temperatureFactorMults[i][atom][symOp][j] = mult;
                            acc*=c_pow;
                        }
                    }
                }
            }
        }

        printStep("temperature factor multipliers");

        // [dir][atom][symOpMult | 1][n][m]
        std::array<std::vector<std::vector<std::vector<std::vector<REAL>>>>, 3> temperatureFactorMultsSquare;
        for (int i=0; i<3; i++){
            temperatureFactorMultsSquare[i].resize(nAtoms);
            for (int atom=0; atom<nAtoms; atom++){
                temperatureFactorMultsSquare[i][atom].resize((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount);
                for (int symOp=0; symOp<((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount); symOp++){
                    temperatureFactorMultsSquare[i][atom][symOp].resize(binSize[(i+1)%3]);
                    for (int n=0; n<binSize[(i+1)%3];n++){
                        temperatureFactorMultsSquare[i][atom][symOp][n].resize(binSize[(i+2)%3]);
                    }
                }
            }
        }
        #pragma omp parallel for num_threads(nThreads) collapse(2)
        for (int i=0; i<3; i++){
            for (int atom=0; atom<nAtoms; atom++){
                bool iso = (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1);
                int ns = (i+1)%3;
                int ms = (i+2)%3;
                if (iso) {
                    REAL c = std::exp(
                        -2.0*atomic_displacement_parameters[usedAtomIndices[atom]][0] * (
                            stepCartesian[ns][0]*stepCartesian[ms][0] +
                            stepCartesian[ns][1]*stepCartesian[ms][1] +
                            stepCartesian[ns][2]*stepCartesian[ms][2]));
                    REAL mult = 1.0;
                    for (int n=0; n<binSize[ns]; n++){
                        REAL acc = 1.0;
                        for (int m=0; m<binSize[ms]; m++){
                            temperatureFactorMultsSquare[i][atom][0][n][m]=acc;
                            acc *= mult;
                        }
                        mult *= c;
                    }
                } else {
                    for (int symOp=0; symOp<symOpMultCount; symOp++){
                        Matrix3<REAL> currentU = rotatedUs[atom][symOp];
                        REAL c = std::exp(-2.0*(
                            stepCartesian[ns][0]*stepCartesian[ms][0]*currentU(0,0) +
                            stepCartesian[ns][1]*stepCartesian[ms][1]*currentU(1,1) +
                            stepCartesian[ns][2]*stepCartesian[ms][2]*currentU(2,2) +
                            (stepCartesian[ns][0]*stepCartesian[ms][1] + stepCartesian[ns][1]*stepCartesian[ms][0])*currentU(0,1) +
                            (stepCartesian[ns][0]*stepCartesian[ms][2] + stepCartesian[ns][2]*stepCartesian[ms][0])*currentU(0,2) +
                            (stepCartesian[ns][1]*stepCartesian[ms][2] + stepCartesian[ns][2]*stepCartesian[ms][1])*currentU(1,2)));
                        REAL mult = 1.0;
                        for (int n=0; n<binSize[ns]; n++){
                            REAL acc = 1.0;
                            for (int m=0; m<binSize[ms]; m++){
                                temperatureFactorMultsSquare[i][atom][symOp][n][m]=acc;
                                acc *= mult;
                            }
                            mult *= c;
                        }
                    }
                }
            }
        }

        printStep("temperature factor square multipliers");

        // [h*nk*nl + k*nl + l][atom][symOpMult | 1]
        std::vector<std::vector<std::vector<REAL>>> virtHklTemperature;
        if constexpr (virtHklTemperatureFlag) {
            virtHklTemperature.resize(binSize[0]*binSize[1]*binSize[2]);
            for (int hklIdx = 0; hklIdx<(binSize[0]*binSize[1]*binSize[2]); hklIdx++){
                virtHklTemperature[hklIdx].resize(nAtoms);
                for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                    virtHklTemperature[hklIdx][atomIdx].resize((atomic_displacement_parameters[usedAtomIndices[atomIdx]].size() == 1) ? 1 : symOpMultCount);
                }
            }
            #pragma omp parallel for num_threads(nThreads) collapse(4)
            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                for (int h = 0; h<binSize[0]; h++){
                    for (int k = 0; k<binSize[1]; k++){
                        for (int l = 0; l<binSize[2]; l++){
                            for (int symOp = 0; symOp<((atomic_displacement_parameters[usedAtomIndices[atomIdx]].size() == 1) ? 1 : symOpMultCount); symOp++){
                                virtHklTemperature[(h*binSize[1] + k)*binSize[2] + l][atomIdx][symOp] =
                                temperatureFactorMults[0][atomIdx][symOp][h] *
                                temperatureFactorMults[1][atomIdx][symOp][k] *
                                temperatureFactorMults[2][atomIdx][symOp][l] *
                                temperatureFactorMultsSquare[0][atomIdx][symOp][k][l] *
                                temperatureFactorMultsSquare[1][atomIdx][symOp][l][h] *
                                temperatureFactorMultsSquare[2][atomIdx][symOp][h][k];
                            }
                        }
                    }
                }
            }
        }

        printStep("virtual hkl temperature factors");

        // [bin][atom][symOpMult | 1]
        std::vector<std::vector<std::vector<Vector3<REAL>>>> perBinTemperatureFactorMult;
        perBinTemperatureFactorMult.resize(binCount);
        for (int binIdx=0; binIdx<binCount; binIdx++){
            perBinTemperatureFactorMult[binIdx].resize(nAtoms);
            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                perBinTemperatureFactorMult[binIdx][atomIdx].resize((atomic_displacement_parameters[usedAtomIndices[atomIdx]].size() == 1)?1:symOpMultCount);
            }
        }
        #pragma omp parallel for num_threads(nThreads) collapse(3)
        for (int binIdx=0; binIdx<binCount; binIdx++){
            for (int atom=0; atom<nAtoms; atom++){
                for (int i=0; i<3; i++){
                    if (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1){
                        perBinTemperatureFactorMult[binIdx][atom][0][i] = std::exp(
                            -2.0 * atomic_displacement_parameters[usedAtomIndices[atom]][0] *
                            (bin000Cartesian[binIdx] * stepCartesian[i]));
                    } else {
                        for (int symOp = 0; symOp<symOpMultCount; symOp++){
                            Matrix3<REAL> currentU = rotatedUs[atom][symOp];
                            perBinTemperatureFactorMult[binIdx][atom][symOp][i] = std::exp(-2.0 * (
                                bin000Cartesian[binIdx][0] * stepCartesian[i][0] * currentU(0, 0) +
                                bin000Cartesian[binIdx][1] * stepCartesian[i][1] * currentU(1, 1) +
                                bin000Cartesian[binIdx][2] * stepCartesian[i][2] * currentU(2, 2) +
                                (bin000Cartesian[binIdx][0] * stepCartesian[i][1] + bin000Cartesian[binIdx][1] * stepCartesian[i][0]) * currentU(0, 1) +
                                (bin000Cartesian[binIdx][0] * stepCartesian[i][2] + bin000Cartesian[binIdx][2] * stepCartesian[i][0]) * currentU(0, 2) +
                                (bin000Cartesian[binIdx][1] * stepCartesian[i][2] + bin000Cartesian[binIdx][2] * stepCartesian[i][1]) * currentU(1, 2)));
                        }
                    }
                }
            }
        }

        printStep("per bin temperature factor multipliers");

        // [bin][dir][n][atom][symOpMult | 1]
        std::vector<std::array<std::vector<std::vector<std::vector<REAL>>>, 3>> virtHklTMul;
        if constexpr (virtHklTMulFlag) {
            virtHklTMul.resize(binCount);
            for (int binIdx=0; binIdx<binCount; binIdx++){
                for (int i=0; i<3; i++){
                    virtHklTMul[binIdx][i].resize(binSize[i]);
                    for (int n=0; n<binSize[i];n++){
                        virtHklTMul[binIdx][i][n].resize(nAtoms);
                        for (int atom=0; atom<nAtoms; atom++){
                            virtHklTMul[binIdx][i][n][atom].resize((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount);
                        }
                    }
                }
            }
            #pragma omp parallel for num_threads(nThreads) collapse(3)
            for (int binIdx=0; binIdx<binCount; binIdx++){
                for (int i=0; i<3; i++){
                    for (int atom=0; atom<nAtoms; atom++){
                        for (int symOp = 0; symOp<((atomic_displacement_parameters[usedAtomIndices[atom]].size() == 1) ? 1 : symOpMultCount); symOp++){
                            REAL mult = 1.0;
                            for (int n=0; n<binSize[i]; n++){
                                virtHklTMul[binIdx][i][n][atom][symOp] = mult;
                                mult *= perBinTemperatureFactorMult[binIdx][atom][symOp][i];
                            }
                        }
                    }
                }
            }
        }

        printStep("virtual hkl temperature multipliers");

        std::vector<int> usedWfns;
        usedWfns.clear();
        std::vector<int> usedTypes;
        usedTypes.clear();
        std::vector<std::array<int, 2>> usedWfnTypeCombo;
        usedWfnTypeCombo.clear();
        std::vector<int> atomToUsedWfnTypeCombo;
        atomToUsedWfnTypeCombo.resize(nAtoms);
        for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
            const int trueAtomIdx = usedAtomIndices[atomIdx];
            const int wfnIdx = atom_to_wfn_map[trueAtomIdx];
            const int typeIdx = atom_to_type_map[trueAtomIdx];

            int currentUsedWfnIdx = -1;
            for (int usedWfnIdx = 0; usedWfnIdx<usedWfns.size(); usedWfnIdx++){
                if (usedWfns[usedWfnIdx] == wfnIdx){
                    currentUsedWfnIdx = usedWfnIdx;
                    break;
                }
            }
            if (currentUsedWfnIdx<0){
                currentUsedWfnIdx = usedWfns.size();
                usedWfns.emplace_back(wfnIdx);
            }

            int currentUsedTypeIdx = -1;
            for (int usedTypeIdx = 0; usedTypeIdx<usedTypes.size(); usedTypeIdx++){
                if (usedTypes[usedTypeIdx] == typeIdx){
                    currentUsedTypeIdx = usedTypeIdx;
                    break;
                }
            }
            if (currentUsedTypeIdx<0){
                currentUsedTypeIdx = usedTypes.size();
                usedTypes.emplace_back(typeIdx);
            }

            int currentUsedComboIdx = -1;
            for (int usedComboIdx = 0; usedComboIdx<usedWfnTypeCombo.size(); usedComboIdx++){
                if ((usedWfnTypeCombo[usedComboIdx][0] == currentUsedWfnIdx) and (usedWfnTypeCombo[usedComboIdx][1] == currentUsedTypeIdx)){
                    currentUsedComboIdx = usedComboIdx;
                    break;
                }
            }
            if (currentUsedComboIdx<0){
                currentUsedComboIdx = usedWfnTypeCombo.size();
                std::array<int, 2> pair;
                pair[0] = currentUsedWfnIdx;
                pair[1] = currentUsedTypeIdx;
                usedWfnTypeCombo.emplace_back(pair);
            }

            atomToUsedWfnTypeCombo[atomIdx] = currentUsedComboIdx;
        }

        const int wfnCount = usedWfns.size();
        const int typeCount = usedTypes.size();
        const int comboCount = usedWfnTypeCombo.size();

        printStep("wfn and type deduplication");

        // [wfn][l]
        std::vector<std::vector<std::complex<REAL>>> N;
        N.resize(wfnParams.size());
        for (int wfnIdx=0; wfnIdx<wfnCount; wfnIdx++){
            N[usedWfns[wfnIdx]].resize(wfnParams[usedWfns[wfnIdx]].def_valence_pow.size());
        }
        #pragma omp parallel for num_threads(nThreads)
        for (int wfnIdx=0; wfnIdx<wfnCount; wfnIdx++){
            const auto &wfn = wfnParams[usedWfns[wfnIdx]];
            const int nl = wfn.def_valence_pow.size();
            for (int l=0; l<nl; l++){
                std::complex<REAL> perL;
                const int rest = l%4;
                if (rest==0){
                    perL = {four_pi, 0.0};
                } else if (rest==1){
                    perL = {0.0, four_pi};
                } else if (rest==2){
                    perL = {-four_pi, 0.0};
                } else {
                    perL = {0.0, -four_pi};
                }
                N[usedWfns[wfnIdx]][l] = perL *
                sto_atomic_wfn::stoDensityNormalizationFactor(
                    wfn.def_valence_pow[l],
                    wfn.def_valence_exp);
            }
        }

        printStep("density normalization factor");

        f.resize(hklCount);
        #pragma omp parallel for num_threads(nThreads) schedule(guided)
        for (int hklIdx = 0; hklIdx<hklCount; hklIdx++){
            printInLoop("begin");
            valueLog("hklIdx", hklIdx);
            valueLog("hkl_indices[hklIdx][0]", hkl_indices[hklIdx][0]);
            valueLog("hkl_indices[hklIdx][1]", hkl_indices[hklIdx][1]);
            valueLog("hkl_indices[hklIdx][2]", hkl_indices[hklIdx][2]);
            Vector3<int> origBin;
            for (int i=0; i<3; i++){
                origBin[i] = (hkl_indices[hklIdx][i] - minHkl[i]) / (step[i] * binSize[i]);
            }
            int binIdx = allBinsMap[(origBin[0]*nBins[1] + origBin[1])*nBins[2] + origBin[2]];

            printInLoop("bin index");
            valueLog("binIdx", binIdx);

            std::complex<REAL> f_acc = 0.0;

            Vector3<int> offset;
            for (int i=0; i<3; i++){
                offset[i] = ((hkl_indices[hklIdx][i] - minHkl[i]) / step[i]) - (usedBins[binIdx][i] * binSize[i]);
                valueLog("offset[i]", offset[i]);
            }

            printInLoop("offset");

            Vector3<REAL> cartesianH;
            recUnitCell.fractionalToCartesian({
                hkl_indices[hklIdx][0],
                hkl_indices[hklIdx][1],
                hkl_indices[hklIdx][2]}, cartesianH);

            valueLog("cartesianH[0]", cartesianH[0]);
            valueLog("cartesianH[1]", cartesianH[1]);
            valueLog("cartesianH[2]", cartesianH[2]);

            printInLoop("fractional to cartesian conversion");

            REAL hLength = 0.0;
            for (int i=0; i<3; i++){
                hLength += square(cartesianH[i]);
            }
            hLength = sqrt(hLength);

            printInLoop("length of h");

            std::vector<REAL> f_core;
            f_core.resize(wfnCount);
            #pragma omp simd
            for (int wfnIdx=0; wfnIdx<wfnCount; wfnIdx++){
                const auto &wfn = wfnParams[usedWfns[wfnIdx]];
                const int kMax = wfn.core_coeff.size();
                for (int k=0; k<kMax; k++)
                    f_core[wfnIdx] += wfn.core_coeff[k] * sto_scattering::gFunction(0, wfn.core_pow[k] +2, hLength, wfn.core_exp[k]); // TODO find out why pow+2 in all gFunction pow
                    f_core[wfnIdx] *= four_pi;
            }

            printInLoop("core factor");

            std::vector<REAL> val;
            val.resize(comboCount);
            #pragma omp simd
            for (int comboIdx=0; comboIdx<comboCount; comboIdx++){
                const auto &combo = usedWfnTypeCombo[comboIdx];
                const auto &wfn = wfnParams[usedWfns[combo[0]]];
                const auto &type = typeParams[usedTypes[combo[1]]];
                const int kMax = wfn.valence_coeff.size();
                const auto h = hLength/type.kappa_spherical;
                for (int k=0; k<kMax; k++)
                    val[comboIdx] += wfn.valence_coeff[k] * sto_scattering::gFunction(0, wfn.valence_pow[k] +2, h, wfn.valence_exp[k]);
                val[comboIdx] *= type.p_val * four_pi;
            }

            printInLoop("valence component");

            Vector3<REAL> hVersor;
            for (int i=0; i<3; i++){
                hVersor[i] = closeToZero(hLength) ? 0.0 : cartesianH[i]/hLength;
            }

            printInLoop("h versor");

            std::vector<std::vector<std::complex<REAL>>> virtHklPhaseCurrent;
            if constexpr (virtHklPhaseFlag)
                virtHklPhaseCurrent = virtHklPhase[(offset[0]*binSize[1] + offset[1])*binSize[2] + offset[2]];

            printInLoop("virtual hkl phase factors retreival");

            std::vector<std::vector<REAL>> virtHklTemperatureCurrent;
            if constexpr (virtHklTemperatureFlag)
                virtHklTemperatureCurrent = virtHklTemperature[(offset[0]*binSize[1] + offset[1])*binSize[2] + offset[2]];

            printInLoop("virtual hkl temperature factors retreival");

            std::array<std::vector<std::vector<REAL>>, 3> virtHklTMulCurrent;
            if constexpr (virtHklTMulFlag){
                for (int i=0; i<3; i++)
                    virtHklTMulCurrent[i] = virtHklTMul[binIdx][i][offset[i]];
            }

            printInLoop("virtual hkl temperature multipliers retreival");

            #pragma omp simd
            for (int atomIdx = 0; atomIdx<nAtoms; atomIdx++){
                std::complex<REAL> perAtomF = 0.0;
                std::vector<std::complex<REAL>> symOpFMult;
                symOpFMult.resize(symOpMultCount);
                for (int symOpIdx = 0; symOpIdx<nSymOps; symOpIdx++){
                    std::complex<REAL> localF = phaseFactorRoots[binIdx][atomIdx][symOpIdx];
                    valueLog("0 - localF", localF);
                    if constexpr (virtHklPhaseFlag){
                        localF *= virtHklPhaseCurrent[atomIdx][symOpIdx];
                    } else {
                        localF *=
                        phaseFactorMults[0][atomIdx][symOpIdx][offset[0]] *
                        phaseFactorMults[1][atomIdx][symOpIdx][offset[1]] *
                        phaseFactorMults[2][atomIdx][symOpIdx][offset[2]];
                    }
                    valueLog("1 - localF", localF);
                    symOpFMult[symOpToMult[symOpIdx]] += localF;
                }

                //const auto hVersorLocal = local_coordinate_systems[usedAtomIndices[atomIdx]] * hVersor;

                for (int symOpIdx=0; symOpIdx<symOpMultCount; symOpIdx++){

                    bool iso = atomic_displacement_parameters[usedAtomIndices[atomIdx]].size() == 1;

                    REAL localF = temperatureFactorRoots[binIdx][atomIdx][iso?0:symOpIdx];
                    valueLog("2 - localF", localF);
                    if constexpr (virtHklTemperatureFlag) {
                        localF *= virtHklTemperatureCurrent[atomIdx][iso?0:symOpIdx];
                    } else {
                        localF *=
                        temperatureFactorMults[0][atomIdx][iso?0:symOpIdx][offset[0]] *
                        temperatureFactorMults[1][atomIdx][iso?0:symOpIdx][offset[1]] *
                        temperatureFactorMults[2][atomIdx][iso?0:symOpIdx][offset[2]] *
                        temperatureFactorMultsSquare[0][atomIdx][iso?0:symOpIdx][offset[1]][offset[2]] *
                        temperatureFactorMultsSquare[1][atomIdx][iso?0:symOpIdx][offset[2]][offset[0]] *
                        temperatureFactorMultsSquare[2][atomIdx][iso?0:symOpIdx][offset[0]][offset[1]];
                    }
                    valueLog("3 - localF", localF);

                    for (int i=0; i<3; i++){
                        if constexpr (virtHklTMulFlag)
                            localF *= virtHklTMulCurrent[i][atomIdx][iso?0:symOpIdx];
                        else
                            localF *= pow(perBinTemperatureFactorMult[binIdx][atomIdx][iso?0:symOpIdx][i], offset[i]);
                    }

                    valueLog("4 - localF", localF);

                    const int wfnIdx = atom_to_wfn_map[usedAtomIndices[atomIdx]];
                    const auto &wfn = wfnParams[wfnIdx];
                    const auto &type = typeParams[atom_to_type_map[usedAtomIndices[atomIdx]]];

                    const auto h = (hVersor * symOpMults[symOpIdx]) * local_coordinate_systems[usedAtomIndices[atomIdx]];
                    valueLog("(square(h[0]) + square(h[1]) + square(h[2]))", (square(h[0]) + square(h[1]) + square(h[2])));

                    const int nl = std::min(wfn.def_valence_pow.size(), type.p_lm.size());
                    std::complex<REAL> dval;
                    for (int l=0; l<nl; l++){
                        // may be ordered differently than in publication because the publication doesn't seem to have a consistent ordering of arguments passed to g
                        const REAL multPerL = sto_scattering::gFunction(l, wfn.def_valence_pow[l] +2, hLength / type.kappa_def_valence, wfn.def_valence_exp);

                        REAL sumPerM = 0.0;
                        for (int m=-l; m<=l; m++){
                            sumPerM += type.p_lm[l][m+l] * real_spherical_harmonics::densityNormalized(h, l, m);
                        }
                        dval += N[wfnIdx][l] * (multPerL * sumPerM);
                    }

                    valueLog("symOpFMult[symOpIdx]", symOpFMult[symOpIdx]);
                    valueLog("dval", dval);
                    valueLog("f_core[usedWfnTypeCombo[atomToUsedWfnTypeCombo[atomIdx]][0]]", f_core[usedWfnTypeCombo[atomToUsedWfnTypeCombo[atomIdx]][0]]);
                    valueLog("val[atomToUsedWfnTypeCombo[atomIdx]]", val[atomToUsedWfnTypeCombo[atomIdx]]);

                    perAtomF += symOpFMult[symOpIdx] * localF * (dval + 1.0 * f_core[usedWfnTypeCombo[atomToUsedWfnTypeCombo[atomIdx]][0]] + val[atomToUsedWfnTypeCombo[atomIdx]]);
                }
                f_acc += perAtomF * atomic_occupancy[usedAtomIndices[atomIdx]] * atomic_multiplicity_factor[usedAtomIndices[atomIdx]];
            }
            printInLoop("atom loop");
            f[hklIdx] = f_acc;
            printInLoop("returning results");
        }

        printStep("main loop");

        // TODO add centrosymmetry support
        // TODO derivatives
        // TODO find out how to implement implement wfn.anomalous_scattering, electron
        /*
         *    bool mUseIAM;
         *    std::vector<std::string> mIamAtomType;
         *    std::vector<int> mAtomToIamTypeMap;
         *    std::vector<NGaussianFormFactor> mIamFormFactors;
         */ // TODO use this
    } //calculateSF_parallel_2

    void HansenCoppens_SF_Engine4::electronScatteringAt000(
        const std::vector<int>& atomic_numbers,
        std::vector<double>& f)
    {
        map<int, int> z_2_ff_idx;
        set<int> unique_z(atomic_numbers.begin(), atomic_numbers.end());
        vector<int> unique_z_vec(unique_z.begin(), unique_z.end());
        for (int i = 0; i < unique_z_vec.size(); i++)
            z_2_ff_idx[unique_z_vec[i]] = i;

        vector<double> ff_type;
        for (int i = 0; i < unique_z_vec.size(); i++)
            ff_type.push_back(n_gaussian_form_factors_table::getFormFactor(periodic_table::symbol(unique_z_vec[i]), "electron-IT").calculate_h(0.0));

        int atomIdx, nAtoms = atomic_numbers.size();
        f.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            f[atomIdx] = ff_type[z_2_ff_idx[atomic_numbers[atomIdx]]];
    }




    //---------------
    /*
     * void calculateSF(
     *    const UnitCell &unitCell,
     *    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
     *    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
     *    const std::vector<int> &atom_to_wfn_map,
     *    const std::vector<int> &atom_to_type_map,
     *    const std::vector<Vector3<REAL> > &atomicPositions,
     *    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
     *    const std::vector<REAL> &atomic_occupancy,
     *    const std::vector<REAL> &atomic_multiplicity_factor,
     *    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
     *    const std::vector<sf_engine_data_types::SymmetryOperation> &symmetry_operations,
     *    bool centrosymmetric,
     *    const Vector3<REAL> &inversionTranslation,
     *    const std::vector<Vector3<REAL> > &h_vectors,
     *    const std::vector<Vector3i >& hkl_indices,
     *    std::vector<std::complex<REAL> > &f,
     *    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
     *    const std::vector<std::complex<REAL> > &dTarget_df,
     *    const std::vector<bool> &include_atom_contribution,
     *    int nThreads);
     *
     */

    void HansenCoppens_SF_Engine4::select_P10P20_atoms(
        const std::vector<sf_engine_data_types::HC_TypeParam>& type_parameters,
        const std::vector<int>& atom_to_type_map,
        std::vector<bool>& atom_selection)
    {
        int nTypes = type_parameters.size();
        int nAtoms = atom_to_type_map.size();
        vector<bool> pz_dz_type(nTypes, false);
        for (int typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {
            auto const& plms = type_parameters[typeIdx].p_lm;
            if (plms.size() == 3)
                if (plms[0][0] == 0.0)
                {
                    pz_dz_type[typeIdx] = true;
                    for (int l = 1; l <= 2; l++)
                        for (int i = 0; i < 2 * l + 1; i++)
                            if (plms[l][i] != 0.0)
                                if (l != i)
                                    pz_dz_type[typeIdx] = false;
                }
        }
        atom_selection.resize(nAtoms);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            atom_selection[atomIdx] = pz_dz_type[atom_to_type_map[atomIdx]];

    }

} // namespace discamb
