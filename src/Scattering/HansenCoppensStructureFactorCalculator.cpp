#ifndef M_PI
#define 	M_PI   3.14159265358979323846
#endif


#include "discamb/Scattering/HansenCoppensStructureFactorCalculator.h"

#include "discamb/config.h"

#ifdef BUILD_FOR_GPU
    #include "cuda/HC_SF_GPU.h" 
#endif

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/HC_Model/HC_WfnData.h"
#include "discamb/MathUtilities/real_spherical_harmonics.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/Scattering/scattering_utilities.h"

#include "discamb/Scattering/HansenCoppens_SF_Engine2.h"
#include "discamb/Scattering/HansenCoppens_SF_Engine3.h"
#include <fstream>
#include <iostream>
#include <exception>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif


using namespace std;


namespace discamb {

namespace spc = structural_parameters_convention;

HansenCoppensStructureFactorCalculator::HansenCoppensStructureFactorCalculator()
{
	mWithGradients = true;
	mEngineType = CPU;
	mN_Threads = 1;
	mPlmWfnNormalized = false;
	mWallClockTime = 0;
	mCpuTime = 0;
	mSingleH_DataInitialized = false;
    mD_AdpConvention = structural_parameters_convention::AdpConvention::U_cif;
    mD_XyzConvention = structural_parameters_convention::XyzCoordinateSystem::fractional;

}

HansenCoppensStructureFactorCalculator::HansenCoppensStructureFactorCalculator(
    const Crystal &crystal, 
    const HC_ModelParameters &parameters)
{
    mWithGradients = true;
    mEngineType = CPU;
    mN_Threads = 1;
    mPlmWfnNormalized = false;
    mWallClockTime = 0;
    mCpuTime = 0;
	mSingleH_DataInitialized = false;

    mD_AdpConvention = crystal.adpConvention;
    mD_XyzConvention = crystal.xyzCoordinateSystem;

    setModel(crystal, parameters);
}


HansenCoppensStructureFactorCalculator::~HansenCoppensStructureFactorCalculator()
{
}


void HansenCoppensStructureFactorCalculator::getCalculationTime(
    double &cpuTime,
    double &wallClockTime)
const
{
    cpuTime = mCpuTime;
    wallClockTime = mWallClockTime;
}

void HansenCoppensStructureFactorCalculator::setN_Threads(
    int n)
{
    mN_Threads = n;
    #ifdef _OPENMP
        omp_set_num_threads(n);
    #endif

}

void HansenCoppensStructureFactorCalculator::setEngineType(
    HC_SF_EngineId engineId)
{
    mEngineType = engineId;
    setEngineSpecificSpaceGroupRepresentation();
}

void HansenCoppensStructureFactorCalculator::setEngineSpecificSpaceGroupRepresentation()
{
    Matrix3d frac2cart = mUnitCell.getFractionalToCartesianMatrix();
    Matrix3d cart2frac = mUnitCell.getCartesianToFractionalMatrix();
    Matrix3d rotationFractional;
    Vector3d translationFractional;
    int symmOpIndex, nSymmOp;
    bool engineHandleInversionCenter; 
    
    engineHandleInversionCenter = (mEngineType == CPU || mEngineType == CPU_IAM);

    nSymmOp = mSpaceGroup.nSymmetryOperationsInSubset();

    if(engineHandleInversionCenter)
        mSymmetryOperations.resize( nSymmOp );
    else
    {
        if(mIsCentrosymmetric)
            mSymmetryOperations.resize( 2 * nSymmOp );
        else
            mSymmetryOperations.resize( nSymmOp );
    }

    for (symmOpIndex = 0; symmOpIndex < nSymmOp; symmOpIndex++)
    {
        mSpaceGroup.getSpaceGroupOperation(0, 0, symmOpIndex).get(rotationFractional, translationFractional);
        mSymmetryOperations[ symmOpIndex ].rotation = frac2cart * rotationFractional * cart2frac;
        mSymmetryOperations[ symmOpIndex ].translation = frac2cart * translationFractional;
    }

    if(!engineHandleInversionCenter && mIsCentrosymmetric)
        for (symmOpIndex = 0; symmOpIndex < nSymmOp; symmOpIndex++)
        {
            mSpaceGroup.getSpaceGroupOperation(0, 1, symmOpIndex).get(rotationFractional, translationFractional);
            mSymmetryOperations[ symmOpIndex + nSymmOp ].rotation = frac2cart * rotationFractional * cart2frac;
            mSymmetryOperations[ symmOpIndex + nSymmOp ].translation = frac2cart * translationFractional;
        }



}

void HansenCoppensStructureFactorCalculator::setModel(
    const Crystal &crystal,
    const HC_ModelParameters &parameters)
{
	mSingleH_DataInitialized = false;

    mUnitCell = crystal.unitCell;
    mConverter.set(mUnitCell);
    mAdpInputType = crystal.adpConvention;
    mXyzInputCoordinateSystem = crystal.xyzCoordinateSystem;


    crystal_structure_utilities::atomicNumbers(crystal, mAtomicNumbers);

    // space group related things

    mSpaceGroup = crystal.spaceGroup;
    mIsCentrosymmetric = mSpaceGroup.isCentrosymmetric();
    Vector3<CrystallographicRational> inversionCenterTranslationFract;
    if (mIsCentrosymmetric)
    {
        mSpaceGroup.inversionCenterTranslation(inversionCenterTranslationFract);
        mUnitCell.fractionalToCartesian(inversionCenterTranslationFract, mInversionCenterTranslation);
    }
    else
        mInversionCenterTranslation.set(0, 0, 0);
    
    setEngineSpecificSpaceGroupRepresentation();

    // -----------

    int wfnTypeIndex, nWfnTypes = parameters.wfn_parameters.size();
    vector<double> orb_coeff, orb_exp;
    vector<int> orb_pow;
    int orbitalIndex, nOrbitals, term, nTerms;
    double totalValenceOrbitalOccupancy;

    mWfnParameters.clear();
    mWfnParameters.resize(nWfnTypes);
    mWfnTypeLabels.resize(nWfnTypes);

    for (wfnTypeIndex = 0; wfnTypeIndex < nWfnTypes; wfnTypeIndex++)
    {
        HC_WfnType const & wfn_type = parameters.wfn_parameters[wfnTypeIndex];

        mWfnTypeLabels[wfnTypeIndex] = wfn_type.wfnAtomTypeName;


        sto_atomic_wfn::atomicStoWfnToSphericalDensity(wfn_type,
                                                       mWfnParameters[wfnTypeIndex].core_coeff,
                                                       mWfnParameters[wfnTypeIndex].core_exp,
                                                       mWfnParameters[wfnTypeIndex].core_pow,
                                                       sto_atomic_wfn::CORE_DENSITY);

        sto_atomic_wfn::atomicStoWfnToSphericalDensity(wfn_type,
                                                       mWfnParameters[wfnTypeIndex].valence_coeff,
                                                       mWfnParameters[wfnTypeIndex].valence_exp,
                                                       mWfnParameters[wfnTypeIndex].valence_pow,
                                                       sto_atomic_wfn::VALENCE_DENSITY);


        nOrbitals = wfn_type.valence_orbitals_indices.size();
        totalValenceOrbitalOccupancy = 0;

        for (int i = 0; i < nOrbitals; i++)
        {
            orbitalIndex = wfn_type.valence_orbitals_indices[i];
            totalValenceOrbitalOccupancy += wfn_type.orbital_occupancy[orbitalIndex];
        }


        nTerms = mWfnParameters[wfnTypeIndex].valence_coeff.size();
        for (term = 0; term<nTerms; term++)
            mWfnParameters[wfnTypeIndex].valence_coeff[term] /= totalValenceOrbitalOccupancy;

        mWfnParameters[wfnTypeIndex].anomalous_scattering = wfn_type.anomalous_scattering;
        mWfnParameters[wfnTypeIndex].def_valence_exp = wfn_type.deformation_valence_exponent;
        mWfnParameters[wfnTypeIndex].label = wfn_type.wfnAtomTypeName;
        mWfnParameters[wfnTypeIndex].def_valence_pow = wfn_type.deformation_valence_power;
    }

    int typeIndex, nTypes;
    int m, l, nL;

    nTypes = parameters.type_parameters.size();

    mTypeParameters.resize(nTypes);

    for (typeIndex = 0; typeIndex<nTypes; typeIndex++)
    {
        nL = int(parameters.type_parameters[typeIndex].p_lm.size());
        mTypeParameters[typeIndex].p_lm.resize(nL);
        for (l = 0; l<nL; l++)
        {
            mTypeParameters[typeIndex].p_lm[l].resize(2 * l + 1);
            for (m = -l; m <= l; m++)
                mTypeParameters[typeIndex].p_lm[l][l + m] = parameters.type_parameters[typeIndex].p_lm[l][l + m];
        }
        mTypeParameters[typeIndex].kappa_def_valence = parameters.type_parameters[typeIndex].kappa_deformation_valence;
        mTypeParameters[typeIndex].kappa_spherical = parameters.type_parameters[typeIndex].kappa_spherical_valence;
        mTypeParameters[typeIndex].p_val = parameters.type_parameters[typeIndex].p_val;
    }

    copyVector(parameters.atom_to_wfn_map, mAtomToWfnTypeMap);
    copyVector(parameters.atom_to_type_map, mAtomToAtomTypeMap);



    int atomIndex, nAtoms;
    nAtoms = crystal.atoms.size();
    mAtomicPositions.resize(crystal.atoms.size());

    mAtomic_displacement_parameters.resize(nAtoms);
    for (atomIndex = 0; atomIndex<nAtoms; atomIndex++)
        mAtomic_displacement_parameters[atomIndex].resize(crystal.atoms[atomIndex].adp.size());

    mAtomicOccupancy.resize(nAtoms);
    mAtomicMultiplicityFactor.resize(nAtoms);
    mAtomicMultiplicityWeights.resize(nAtoms);
    double nSymmOpsAsDouble = double(crystal.spaceGroup.nSymmetryOperations());

    for(int i=0;i<nAtoms;i++)
    {
        mAtomicMultiplicityFactor[i] = crystal.atoms[i].multiplicity;
        mAtomicMultiplicityWeights[i] = crystal.atoms[i].multiplicity / nSymmOpsAsDouble;
    }
	mReciprocalUnitCell.set(crystal.unitCell);

    mReciprocalUnitCell.fractionalToCartesian(Vector3d(1, 0, 0), mA_Star);
    mReciprocalUnitCell.fractionalToCartesian(Vector3d(0, 1, 0), mB_Star);
    mReciprocalUnitCell.fractionalToCartesian(Vector3d(0, 0, 1), mC_Star);


}



void HansenCoppensStructureFactorCalculator::calculateStructureFactors(
const std::vector<AtomInCrystal> &atoms,
const std::vector<Matrix3d> &localCoordinateSystems,
const std::vector<Vector3i> &hkl,
std::vector<std::complex<double> > &f) 

{
    vector<bool> count_atoms(atoms.size(),true);
    calculateStructureFactors(atoms, localCoordinateSystems, hkl,f,count_atoms);
}

void HansenCoppensStructureFactorCalculator::calculateStructureFactors(
//const Crystal &crystal,
const std::vector<AtomInCrystal> &atoms,
const std::vector<Matrix3d> &localCoordinateSystems,
const std::vector<Vector3i> &hkl,
std::vector<std::complex<double> > &f,
const std::vector<bool> &countAtomContribution)

{
    mWithGradients = false;
    vector<complex<double> > fake_dTarget_df(hkl.size(),1.0);

    DerivativesSelector derivativesSelector;
    derivativesSelector.d_adp = false;
    derivativesSelector.d_anom = false;
    derivativesSelector.d_occ = false;
    derivativesSelector.d_xyz = false;

    vector<TargetFunctionAtomicParamDerivatives> dTarget_dparam;
    calculateStructureFactorsAndDerivatives(atoms, localCoordinateSystems,hkl,f,dTarget_dparam,fake_dTarget_df,countAtomContribution, derivativesSelector);
    mWithGradients = true;
}

void HansenCoppensStructureFactorCalculator::useCPU()
{
    setEngineType(CPU);
}

void HansenCoppensStructureFactorCalculator::useGPU()
{
    setEngineType(GPU);
}

void HansenCoppensStructureFactorCalculator::useCPU_IAM()
{
    setEngineType(CPU_IAM); 
}


// CPU, GPU or CPU_IAM
void HansenCoppensStructureFactorCalculator::setCalculationMode(
    const std::string &engineType)
{
    if(engineType == string("CPU"))
        useCPU();
    else
    {
        if(engineType == string("GPU"))
            useGPU();
        else
        {
            if(engineType == string("CPU_IAM"))
                useCPU_IAM();
            else
                on_error::throwException("unknown calculation mode for multipolar structure factor engine", __FILE__, __LINE__);
        }
    }
}

void HansenCoppensStructureFactorCalculator::calculateStructureFactorsAndDerivatives(
    const std::vector<AtomInCrystal>& atoms,
    const std::vector<Matrix3d>& localCoordinateSystems,
    const std::vector<Vector3i>& hkl,
    std::vector<std::complex<double> >& f,
    std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
    const std::vector<std::complex<double> >& _dTarget_df,
    const std::vector<bool>& countAtomContribution,
    const DerivativesSelector& derivativesSelector)
{
    vector<complex<double> > dTarget_df = _dTarget_df;
    vector<double> sfMultipliers;
    char latticeCentering;
    bool obverse;

    latticeCentering = mSpaceGroup.centering(obverse);
    scattering_utilities::centeringSfMultipliers(latticeCentering, hkl, sfMultipliers, obverse);

    int hklIdx, nHkl = hkl.size();

    if (dTarget_df.size() != hkl.size())
        on_error::throwException("Inconsistent size of hkl vectors set and (d target/ d F) table when calculating multipolar structure factors.", __FILE__, __LINE__);

    for (hklIdx = 0; hklIdx < nHkl; ++hklIdx)
        dTarget_df[hklIdx] = conj(dTarget_df[hklIdx]) * sfMultipliers[hklIdx];


    prepareDataForEngine(atoms, hkl, f, dTarget_dparam);


    vector<Matrix3i> rotations;
    vector<Vector3d> translations;
    Vector3d a, b, c, a_star, b_star, c_star;
    Vector3<CrystallographicRational> translation;
    ReciprocalLatticeUnitCell reciprocalLattice;
    Matrix3d rotation;

    switch (mEngineType)
    {
    case GPU:
    {
        Crystal crystal;
        crystal.atoms = atoms;
        crystal.spaceGroup = mSpaceGroup;
        crystal.unitCell = mUnitCell;

        mCpuTimer.start();
        mWallClockTimer.start();
        calculate_using_GPU_batched(crystal, localCoordinateSystems, hkl, f, dTarget_dparam, dTarget_df);


        mCpuTime = mCpuTimer.stop();
        mWallClockTime = mWallClockTimer.stop();

        break;
    }
    case CPU:
    {
        HansenCoppens_SF_Engine3 engine;

        engine.calculateSF(mUnitCell, mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mAtomicPositions,
            mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
            localCoordinateSystems, mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
            mHKL_Cartesian, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution, mN_Threads, derivativesSelector, 
            mElectronScattering, mAtomicNumbers);

        //HansenCoppens_SF_Engine2 engine;

        //engine.calculateSF(mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mAtomicPositions,
        //    mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
        //    localCoordinateSystems, mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
        //    mHKL_Cartesian, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution, mN_Threads);

        break;
    }
    case CPU_IAM:
    {
        HansenCoppens_SF_Engine engine;
        int wfnTypeIdx, nWfnTypes = mWfnParameters.size();
        vector<complex<double> > anomalousScattering(nWfnTypes);

        for (wfnTypeIdx = 0; wfnTypeIdx < nWfnTypes; wfnTypeIdx++)
            anomalousScattering[wfnTypeIdx] = mWfnParameters[wfnTypeIdx].anomalous_scattering;
        mCpuTimer.start();
        mWallClockTimer.start();

        engine.calculateSF_IAM(mWfnTypeLabels, anomalousScattering, mAtomToWfnTypeMap, mAtomicPositions,
            mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
            mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
            mHKL_Cartesian, f, dTarget_dparam, dTarget_df, countAtomContribution,
            mN_Threads);

        mCpuTime = mCpuTimer.stop();
        mWallClockTime = mWallClockTimer.stop();

        break;
    }
    default:
        on_error::throwException("wrong engine type", __FILE__, __LINE__);
    }

    for (hklIdx = 0; hklIdx < nHkl; ++hklIdx)
        f[hklIdx] *= sfMultipliers[hklIdx];

    convertDerivatives(dTarget_dparam);
}

void HansenCoppensStructureFactorCalculator::setElectronScattering(
    bool electronScattering)
{
    mElectronScattering = electronScattering;
}

void HansenCoppensStructureFactorCalculator::calculateStructureFactorsAndDerivatives( 
const std::vector<AtomInCrystal> &atoms,
const std::vector<Matrix3d> &localCoordinateSystems,
const std::vector<Vector3i> &hkl,
std::vector<std::complex<double> > &f,
std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
const std::vector<std::complex<double> > &_dTarget_df,
const std::vector<bool> &countAtomContribution)

{
    DerivativesSelector selector;
    calculateStructureFactorsAndDerivatives(
        atoms,
        localCoordinateSystems,
        hkl,
        f,
        dTarget_dparam,
        _dTarget_df,
        countAtomContribution,
        selector);

    //vector<complex<double> > dTarget_df = _dTarget_df;
    //vector<double> sfMultipliers;
    //char latticeCentering;
    //bool obverse;
    //
    //latticeCentering = mSpaceGroup.centering(obverse);
    //scattering_utilities::centeringSfMultipliers(latticeCentering, hkl, sfMultipliers, obverse);
    //
    //int hklIdx, nHkl =  hkl.size();

    //if(dTarget_df.size() != hkl.size())
    //    on_error::throwException("Inconsistent size of hkl vectors set and (d target/ d F) table when calculating multipolar structure factors.",__FILE__,__LINE__);
    //
    //for( hklIdx = 0 ; hklIdx<nHkl ; ++hklIdx)
    //    dTarget_df[hklIdx] = conj(dTarget_df[hklIdx]) * sfMultipliers[hklIdx];

    //
    //prepareDataForEngine(atoms,hkl,f,dTarget_dparam);
    //

    //vector<Matrix3i> rotations;
    //vector<Vector3d> translations;
    //Vector3d a, b, c, a_star, b_star, c_star;
    //Vector3<CrystallographicRational> translation;
    //ReciprocalLatticeUnitCell reciprocalLattice;
    //Matrix3d rotation;
    //DerivativesSelector selector;
    //switch(mEngineType)
    //{
    //    case GPU:
    //    {
    //        Crystal crystal;
    //        crystal.atoms = atoms;
    //        crystal.spaceGroup = mSpaceGroup;
    //        crystal.unitCell = mUnitCell;

    //        mCpuTimer.start();
    //        mWallClockTimer.start();
    //        calculate_using_GPU_batched(crystal, localCoordinateSystems, hkl, f, dTarget_dparam, dTarget_df);
    //        

    //        mCpuTime = mCpuTimer.stop();
    //        mWallClockTime = mWallClockTimer.stop();

    //        break;
    //    }
    //    case CPU:
    //    {
    //        HansenCoppens_SF_Engine3 engine;
    //        
    //        engine.calculateSF(mUnitCell, mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mAtomicPositions,
    //                           mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
    //                           localCoordinateSystems, mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
    //                           mHKL_Cartesian, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution, mN_Threads);

    //        //HansenCoppens_SF_Engine2 engine;

    //        //engine.calculateSF(mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mAtomicPositions,
    //        //    mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
    //        //    localCoordinateSystems, mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
    //        //    mHKL_Cartesian, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution, mN_Threads);

    //        break;
    //    }   
    //    case CPU_IAM:
    //    {
    //        HansenCoppens_SF_Engine engine;
    //        int wfnTypeIdx, nWfnTypes = mWfnParameters.size();
    //        vector<complex<double> > anomalousScattering( nWfnTypes );

    //        for(wfnTypeIdx=0 ; wfnTypeIdx<nWfnTypes ; wfnTypeIdx++ )
    //            anomalousScattering[wfnTypeIdx] = mWfnParameters[wfnTypeIdx].anomalous_scattering;
    //        mCpuTimer.start();
    //        mWallClockTimer.start();

    //        engine.calculateSF_IAM(mWfnTypeLabels, anomalousScattering, mAtomToWfnTypeMap,mAtomicPositions,
    //                           mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
    //                           mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
    //                               mHKL_Cartesian, f, dTarget_dparam, dTarget_df, countAtomContribution,
    //                           mN_Threads);

    //        mCpuTime = mCpuTimer.stop();
    //        mWallClockTime = mWallClockTimer.stop();

    //        break;
    //    }
    //    default:
    //        on_error::throwException("wrong engine type",__FILE__,__LINE__);
    //} 
    //
    //for (hklIdx = 0; hklIdx<nHkl; ++hklIdx)
    //    f[hklIdx] *= sfMultipliers[hklIdx];
    //
    //convertDerivatives(dTarget_dparam);

}

void HansenCoppensStructureFactorCalculator::convertDerivatives(
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam) 
const
{
    int atomIdx, nAtoms = dTarget_dparam.size();
    vector<complex<double> > adpIn(6), adpOut(6);
    Vector3<complex<double> > xyzIn, xyzOut;
    int i;

    if(mD_AdpConvention != spc::AdpConvention::U_cart)
    {
        for(atomIdx=0; atomIdx<nAtoms; atomIdx++)
            if(dTarget_dparam[atomIdx].adp_derivatives.size()==6)
            {
                for( i=0 ; i<6 ; i++)
                    adpIn[i] = dTarget_dparam[atomIdx].adp_derivatives[i];
                mConverter.convertDerivativesADP(adpIn, adpOut, spc::AdpConvention::U_cart, mD_AdpConvention);
                for ( i = 0; i<6; i++)
                    dTarget_dparam[atomIdx].adp_derivatives[i] = adpOut[i].real();
            }
    }
    if(mD_XyzConvention != spc::XyzCoordinateSystem::cartesian)
        for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
        {
            for (i = 0; i<3; i++)
                xyzIn[i] = dTarget_dparam[atomIdx].atomic_position_derivatives[i];
            mConverter.convertDerivativesXyz(xyzIn, xyzOut, spc::XyzCoordinateSystem::cartesian, mD_XyzConvention);
            for (i = 0; i<3; i++)
                dTarget_dparam[atomIdx].atomic_position_derivatives[i] = xyzOut[i].real();
        }
}

void HansenCoppensStructureFactorCalculator::test(
    const Crystal &crystal, 
    const std::vector<Matrix3d> &localCoordinateSystems,
    const std::vector<Vector3i> &hkl, 
    std::vector<std::complex<double> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<double> > &dTarget_df)
{
    HansenCoppens_SF_Engine engine;

    mCpuTimer.start();
    mWallClockTimer.start();
    vector<bool> countAtomContribution(crystal.atoms.size());

    int i, nHkl = hkl.size();
    vector<Vector3d> hkls(nHkl);
    ReciprocalLatticeUnitCell rl(crystal.unitCell);
    for(i=0;i<nHkl;i++)
        rl.fractionalToCartesian(hkl[i],hkls[i]);

    engine.calculateSF(mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mAtomicPositions,
                       mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
                       localCoordinateSystems, mSymmetryOperations, mIsCentrosymmetric, mInversionCenterTranslation,
                       hkls, f, dTarget_dparam, dTarget_df, countAtomContribution, mN_Threads);

    mCpuTime = mCpuTimer.stop();
    mWallClockTime = mWallClockTimer.stop();
}

void HansenCoppensStructureFactorCalculator::calculate_using_GPU_batched(
    const Crystal &crystal, 
    const std::vector<Matrix3d> &localCoordinateSystems,
    const std::vector<Vector3i> &hkl,
    std::vector<std::complex<double> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dParam,
    const std::vector<std::complex<double> > &dTarget_df)
{
    // only calls to GPU part are commented out for no GPU build
    // so if some incompatibility in API appears we can see it
    // even without GPU build

#ifndef BUILD_FOR_GPU

    string errorMessage = string("requested multipolar structure factors calculation with GPU, ") +
        string("but DiSCaMB library was compiled without support for GPU, ");

    on_error::throwException(errorMessage, __FILE__, __LINE__);

#endif


    bool hasGpu = false;
    int gpuGlobalMemory, gpuMemoryToAlloc;
    string gpuName;

    // to avoid compiler complaints on unitialized
    // variables used in the case of non GPU build 
    gpuGlobalMemory = gpuMemoryToAlloc = 1;

#ifdef BUILD_FOR_GPU
    hc_sf_gpu::gpuInfo(hasGpu, gpuGlobalMemory, gpuName);
#endif

    if(!hasGpu)
        on_error::throwException("GPU computations requested but no GPU has been detected", __FILE__, __LINE__);

#ifdef BUILD_FOR_GPU
    gpuMemoryToAlloc = hc_sf_gpu::gpuOnlyMemory(hkl.size(), crystal.atoms.size());
#endif 

    int nBatches = int(1.2 * gpuMemoryToAlloc) / gpuGlobalMemory + 1;
        

    //------------------------

    vector< vector< complex< REAL > > > batchedSF;
    vector<vector<Vector3i> > batchedHkl;
    vector<vector<complex<double> > > batched_dTarget_dF;
    vector<vector<TargetFunctionAtomicParamDerivatives> > batched_dTarget_dParam;


    // assign HKL vectors to threads

    batchedSF.resize(nBatches);
    batchedHkl.resize(nBatches);
    batched_dTarget_dParam.resize(nBatches, dTarget_dParam);
    batched_dTarget_dF.resize(nBatches);


    scattering_utilities::divideSet(hkl, nBatches, batchedHkl);
    scattering_utilities::divideSet(dTarget_df, nBatches, batched_dTarget_dF);

    for(int i=0;i<nBatches; i++)
    {
        batchedSF[i].resize(batchedHkl[i].size());
        calculate_using_GPU(crystal, localCoordinateSystems, batchedHkl[i], batchedSF[i], 
                            batched_dTarget_dParam[i], batched_dTarget_dF[i]);
    }

    scattering_utilities::combineScatteringFactorsSets(batchedSF, f);

    scattering_utilities::merge_dTarget_dParameterSets(batched_dTarget_dParam, dTarget_dParam);

}



void HansenCoppensStructureFactorCalculator::calculate_using_GPU(
    const Crystal &crystal,
    const std::vector<Matrix3d> &localCoordinateSystems,
    const std::vector<Vector3i> &hkl, 
    std::vector<std::complex<double> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<double> > &dTarget_df)
{
    #ifndef BUILD_FOR_GPU

        string errorMessage = string("requested multipolar structure factors calculation with GPU, ") +
        string("but DiSCaMB library was compiled without support for GPU, ");

        on_error::throwException(errorMessage, __FILE__, __LINE__);

    #endif


    stringstream buffer;

    saveRawData(crystal, localCoordinateSystems, hkl, buffer);

    //--------------------------------------------------------
    int atomIdx, nAtoms = crystal.atoms.size();
    int nHkl = hkl.size();
    vector<double> f_real(nHkl), f_imag(nHkl), dTarget_df_real(nHkl), dTarget_df_imag(nHkl);
    vector<double> occupancyDerivatives(nAtoms), xyzDerivatives(3*nAtoms), adpDerivatives(6*nAtoms);


    for(int i=0;i<nHkl;i++)
    {
        dTarget_df_real[i] = dTarget_df[i].real();
        dTarget_df_imag[i] = dTarget_df[i].imag();
    }
#ifdef BUILD_FOR_GPU    
    hc_sf_gpu::calculate_multipolar_structure_factors_gpu(
        buffer, dTarget_df_real, dTarget_df_imag, f_real, f_imag, xyzDerivatives, adpDerivatives, occupancyDerivatives);
#endif
    
    //------------------------------------------------------
    for(int i=0;i<nHkl;i++)
        f[i] = complex<double>(f_real[i],f_imag[i]);
    
    double two_pi_2 = 2 * M_PI*M_PI;

    for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
    {

        if (dTarget_dparam[atomIdx].adp_derivatives.size() == 1)
        {
            dTarget_dparam[atomIdx].adp_derivatives[0] = adpDerivatives[6 * atomIdx] + adpDerivatives[6 * atomIdx + 1] + adpDerivatives[6 * atomIdx + 2];
        }
        else
        {
            for (int i = 0; i<dTarget_dparam[atomIdx].adp_derivatives.size(); i++)
                dTarget_dparam[atomIdx].adp_derivatives[i] = adpDerivatives[6 * atomIdx + i];
        }
            
        if (dTarget_dparam[atomIdx].adp_derivatives.size() == 6)
            for (int i = 0; i<6; i++)
                dTarget_dparam[atomIdx].adp_derivatives[i] *= two_pi_2;

        if (dTarget_dparam[atomIdx].adp_derivatives.size() == 1)
            dTarget_dparam[atomIdx].adp_derivatives[0] *= two_pi_2;

        for (int i = 0; i<3; i++)
            dTarget_dparam[atomIdx].atomic_position_derivatives[i] = xyzDerivatives[3*atomIdx+i];

        dTarget_dparam[atomIdx].occupancy_derivatives = occupancyDerivatives[atomIdx];
    }
    

}


double HansenCoppensStructureFactorCalculator::slaterNormalizationFactor(
int powerR,
double exponent)
{
    static const double factorial[] = {1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800,479001600};
    double result = pow(2*exponent,int(powerR+1))*sqrt(2*exponent/factorial[2*(powerR+1)]);
    return result;
}


void HansenCoppensStructureFactorCalculator::addOrbitalDensity(
const SlaterTypeAtomicOrbitalRdf &orbital,
int occupancy,
std::vector<double> &density_coeff,
std::vector<int> &density_powers,
std::vector<double> &density_exponents)
{
    int i,j,nPrimitives = orbital.coefficients.size();
    double multiplier;

    for(i=0;i<nPrimitives;i++)
        for(j=0;j<=i;j++)
        {
            i == j ? multiplier = 1 : multiplier = 2;
            multiplier *= occupancy * orbital.coefficients[i] * orbital.coefficients[j];
            density_coeff.push_back(multiplier);
            density_powers.push_back(orbital.power_r[i]+orbital.power_r[j]);
            density_exponents.push_back(orbital.exponents[i]+orbital.exponents[j]);
        }
}

void HansenCoppensStructureFactorCalculator::setDerivativesConvention(
    spc::XyzCoordinateSystem dXyzConvention,
    spc::AdpConvention dAdpConvention)
{
    mD_XyzConvention = dXyzConvention;
    mD_AdpConvention = dAdpConvention;
}

void HansenCoppensStructureFactorCalculator::getDerivativesConvention(
    spc::XyzCoordinateSystem &dXyzConvention,
    spc::AdpConvention &dAdpConvention)
const
{
    dXyzConvention = mD_XyzConvention;
    dAdpConvention = mD_AdpConvention;

}

void HansenCoppensStructureFactorCalculator::prepareDataForEngine_part(
    const std::vector<AtomInCrystal> &atoms,
    const std::vector<Vector3i> &hkl)
const
{
    int atom_index, nAtoms, nAdpComponents;
    double two_pi2 = 2.0*M_PI*M_PI;
    Vector3d x, y, z;
    ReciprocalLatticeUnitCell reciprocalLattice(mUnitCell);
    nAtoms = mAtomicPositions.size();
    spc::XyzCoordinateSystem xyzCoordinatesType;
    spc::AdpConvention adpType;
    

    if(atoms.empty())
        return;


    xyzCoordinatesType = mXyzInputCoordinateSystem; 
    adpType = mAdpInputType;


    // convert ADP
    for (atom_index = 0; atom_index<nAtoms; atom_index++)
    {
        nAdpComponents = atoms[atom_index].adp.size();

        if (nAdpComponents == 6)
        {
            if (adpType == spc::AdpConvention::U_cart)
                mAtomic_displacement_parameters[atom_index] = atoms[atom_index].adp;
            else
                mConverter.convertADP(atoms[atom_index].adp, mAtomic_displacement_parameters[atom_index], adpType, spc::AdpConvention::U_cart);
            
            for (int i = 0; i<6; i++)
                mAtomic_displacement_parameters[atom_index][i] *= two_pi2;
        }
        else
        {
            if (nAdpComponents == 1)
            {
                mAtomic_displacement_parameters[atom_index] = atoms[atom_index].adp;
                mAtomic_displacement_parameters[atom_index][0] *= two_pi2;
            }
        }
    }

    // convert atomic positions
    
    if (xyzCoordinatesType == spc::XyzCoordinateSystem::cartesian)
        for (atom_index = 0; atom_index<nAtoms; atom_index++)
            mAtomicPositions[atom_index] = atoms[atom_index].coordinates;
    else
        for (atom_index = 0; atom_index<nAtoms; atom_index++)
            mConverter.xyzFractionalToCartesian(atoms[atom_index].coordinates, mAtomicPositions[atom_index]);

    // set occupancy

    for (atom_index = 0; atom_index<nAtoms; atom_index++)
        mAtomicOccupancy[atom_index] = atoms[atom_index].occupancy;

    //


    int i, n = hkl.size();
    Vector3d a_star, b_star, c_star;
    mHKL_Cartesian.resize(n);

    reciprocalLattice.fractionalToCartesian(Vector3d(1, 0, 0), a_star);
    reciprocalLattice.fractionalToCartesian(Vector3d(0, 1, 0), b_star);
    reciprocalLattice.fractionalToCartesian(Vector3d(0, 0, 1), c_star);

    for (i = 0; i<n; i++)
        mHKL_Cartesian[i] = a_star*double(hkl[i][0]) + b_star*double(hkl[i][1]) + c_star*double(hkl[i][2]);

    if (mEngineType != CPU && !mPlmWfnNormalized)
    {
        vector<vector<double> > densityToWfnNormalization;
        real_spherical_harmonics::getDensityToWfnMultipliers(4, densityToWfnNormalization);

        int l, m, nL;
        int typeIndex, nTypes = mTypeParameters.size();
        for (typeIndex = 0; typeIndex<nTypes; typeIndex++)
        {
            nL = mTypeParameters[typeIndex].p_lm.size();

            for (l = 0; l<nL; l++)
                for (m = -l; m <= l; m++)
                    mTypeParameters[typeIndex].p_lm[l][l + m] *= densityToWfnNormalization[l][abs(m)];

        }
        mPlmWfnNormalized = true;
    }

    if (mEngineType == CPU && mPlmWfnNormalized)
    {
        vector<vector<double> > densityToWfnMultipliers;//wfnToDensityNormalization;
        real_spherical_harmonics::getDensityToWfnMultipliers(4, densityToWfnMultipliers);
        //real_spherical_harmonics::get

        int l, m, nL;
        int typeIndex, nTypes = mTypeParameters.size();
        for (typeIndex = 0; typeIndex<nTypes; typeIndex++)
        {
            nL = mTypeParameters[typeIndex].p_lm.size();

            for (l = 0; l<nL; l++)
                for (m = -l; m <= l; m++)
                    mTypeParameters[typeIndex].p_lm[l][l + m] /= densityToWfnMultipliers[l][abs(m)];

        }
        mPlmWfnNormalized = false;
    }

}


void HansenCoppensStructureFactorCalculator::prepareDataForEngine(
    const std::vector<AtomInCrystal> &atoms,
    const std::vector<Vector3i> &hkl,
    std::vector<std::complex<double> > &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam)
    const
{
    prepareDataForEngine_part(atoms, hkl);
    int atom_index, nAtoms, nAdpComponents;
    nAtoms = mAtomicPositions.size();

    dTarget_dparam.resize(nAtoms);
    f.resize(hkl.size());

    for (atom_index = 0; atom_index<nAtoms; atom_index++)
    {
        nAdpComponents = atoms[atom_index].adp.size();

        if (nAdpComponents == 6)
            dTarget_dparam[atom_index].adp_derivatives.resize(6);
        else if (nAdpComponents == 1)
            dTarget_dparam[atom_index].adp_derivatives.resize(1);
    }

}


void HansenCoppensStructureFactorCalculator::copyVector(
const std::vector<int> &from,
std::vector<int> &to)
{
    int i,n=from.size();
    to.resize(n);
    for(i=0;i<n;i++)
        to[i]=from[i];
}

void HansenCoppensStructureFactorCalculator::saveRawData(
    const Crystal &crystal,
    const std::vector<Matrix3d> &localCoordinateSystems, 
    const std::vector<Vector3i> &hkl,
    const std::vector<std::complex<double> > &dTarget_df,
    const std::vector<bool> &countAtomContribution,
    std::ostream &out) 
{

    int i, n;

    saveRawData(crystal, localCoordinateSystems, hkl, out);

    out << setprecision(14);

    out << endl << "dTarget / df" << endl;
    out << endl;

    n = dTarget_df.size();

    for (i = 0; i<n; i++)
        out << dTarget_df[i].real() << " " << dTarget_df[i].imag() << endl;

    out << endl << "atom contribution flag" << endl << endl;

    n = countAtomContribution.size();

    for (i = 0; i<n; i++)
        if (countAtomContribution[i])
            out << "1" << endl;
        else
            out << "0" << endl;

}



void HansenCoppensStructureFactorCalculator::saveRawData(
    const Crystal &crystal,
    const std::vector<Matrix3d> &localCoordinateSystems,
    const std::vector<Vector3i> &hkl,
    std::ostream &out)
{
    HC_SF_EngineId engineId = mEngineType;
    setEngineType(GPU);
    prepareDataForEngine_part(crystal.atoms,hkl);
    mEngineType = engineId;

    int atomIdx, nAtoms = mAtomicPositions.size();
    
    for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
    {
        if (mAtomic_displacement_parameters[atomIdx].size() == 6)
        {
            mAtomic_displacement_parameters[atomIdx][3] *= 2.0;
            mAtomic_displacement_parameters[atomIdx][4] *= 2.0;
            mAtomic_displacement_parameters[atomIdx][5] *= 2.0;
        }
    }


    //------------------------------------------------------
    
    int i, j, k, n, p, maxL;
    vector<TargetFunctionAtomicParamDerivatives> x;
    vector<complex<double> > f;
    prepareDataForEngine(crystal.atoms, hkl, f, x);

    out << setprecision(14);

    out << "n_symmetry_operations " << mSymmetryOperations.size() << endl;
    out << "n_wfn_parameter_sets " << mWfnParameters.size() << endl;
    out << "n_atom_types " << mTypeParameters.size() << endl;
    out << "n_atoms " << mAtomicOccupancy.size() << endl;
    out << "n_h_vectors " << mHKL_Cartesian.size() << endl;


    out << endl << "symmetry operations" << endl << endl;

    n = mSymmetryOperations.size();

    for (i = 0; i<n; i++)
    {
        for (j = 0; j<3; j++)
        {
            for (k = 0; k<3; k++)
                out << " " << mSymmetryOperations[i].rotation(j, k);
            out << endl;
        }

        for (j = 0; j<3; j++)
            out << " " << mSymmetryOperations[i].translation[j];
        out << endl << endl;
    }

    out << endl << "wfn parameter sets " << endl << endl;

    n = mWfnParameters.size();
    for (i = 0; i<n; i++)
    {
        // core - coeff,pow,exp
        p = mWfnParameters[i].core_coeff.size();
        out << p << endl;
        for (j = 0; j<p; j++)
            out << 4 * M_PI * mWfnParameters[i].core_coeff[j] << " " << mWfnParameters[i].core_pow[j]
            << " " << mWfnParameters[i].core_exp[j] << endl;

        // spherical valence

        p = mWfnParameters[i].valence_coeff.size();
        out << p << endl;
        for (j = 0; j<p; j++)
            out << 4 * M_PI * mWfnParameters[i].valence_coeff[j] << " " << mWfnParameters[i].valence_pow[j]
            << " " << mWfnParameters[i].valence_exp[j] << endl;

        // deformation valence
        out << mWfnParameters[i].def_valence_exp << endl;
        p = mWfnParameters[i].def_valence_pow.size();
        out << p;
        for (j = 0; j<p; j++)
            out << " " << mWfnParameters[i].def_valence_pow[j];
        out << endl;
    }

    out << "atom type parameter sets " << endl << endl;

    n = mTypeParameters.size();
    for (i = 0; i<n; i++)
    {
        maxL = mTypeParameters[i].p_lm.size() - 1;
        out << mTypeParameters[i].p_val << " " << mTypeParameters[i].kappa_spherical
            << " " << mTypeParameters[i].kappa_def_valence << " " << maxL << endl;

        for (j = 0; j <= maxL; j++)
        {
            for (k = 0; k<2 * j + 1; k++)
                out << mTypeParameters[i].p_lm[j][k] << " ";
            out << endl;
        }

        out << endl;
    }

    out << "atomic data" << endl;

    n = mAtomicPositions.size();

    for (i = 0; i<n; i++)
    {
        // wfn_index type_index
        // position occupancy mutliplicity_f
        // adps
        // lcs

        out << endl;
        out << mAtomToWfnTypeMap[i] << " " << mAtomToAtomTypeMap[i] << endl;



        for (j = 0; j<3; j++)
            out << mAtomicPositions[i][j] << " ";

        out << mAtomicOccupancy[i] << " " << mAtomicMultiplicityWeights[i] << endl;

        p = mAtomic_displacement_parameters[i].size();

        out << p;
        for (j = 0; j<p; j++)
            out << " " << mAtomic_displacement_parameters[i][j];
        out << endl;
        for (j = 0; j<3; j++)
        {
            for (k = 0; k<3; k++)
                out << " " << localCoordinateSystems[i](j, k);
            out << endl;
        }

    }

    out << endl << "h vectors" << endl << endl;
    
    n = hkl.size();
    ReciprocalLatticeUnitCell rl(crystal.unitCell);
    Vector3d hklCartesian,hklFractional;
    for (i = 0; i<n; i++)
    {
        hklFractional = hkl[i];
        rl.fractionalToCartesian(hklFractional, hklCartesian);
        out << hklCartesian[0] << " " << hklCartesian[1] << " " << hklCartesian[2] << endl;
    }

}

void HansenCoppensStructureFactorCalculator::setDataForCpuEngineCallsForSingleHkl()
{
	if (mEngineType != CPU)
		on_error::throwException("using functionality available only for CPU calcluations of multipole model form factors", __FILE__, __LINE__);
	
	int i, n;
	HansenCoppens_SF_Engine engine;
	double step = 0.001;
	int nInterpolationPoints = 10000;

	engine.pre_hkl_loop_sf_calc(mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mType2WfnType, mDefValSlaterNormalization, mTypeMaxL);

	vector<double> h(nInterpolationPoints);
	for (i = 0; i < nInterpolationPoints; i++)
		h[i] = i * step;

	vector<vector<double> > f_core, f_sph_val;
	engine.calculateSphericalTermsInFormFactors(mWfnParameters, mTypeParameters, h, f_core, f_sph_val, mType2WfnType, mDefValSlaterNormalization, mTypeMaxL);

	n = f_core.size();
	mF_CoreInterpolators.resize(n);
	for (i = 0; i < n; i++)
		mF_CoreInterpolators[i].set(f_core[i], step, 0.0);

	n = f_sph_val.size();
	mF_SphericalValenceInterpolators.resize(n);
	for (i = 0; i < n; i++)
		mF_SphericalValenceInterpolators[i].set(f_sph_val[i], step, 0.0);

	mSingleH_DataInitialized = true;
}


void HansenCoppensStructureFactorCalculator::calculateFormFactors(
	const std::vector<Matrix3d>& localCoordinateSystems,
	const Vector3i& hkl,
	std::vector<std::complex<double> >& f,
	const std::vector<bool>& countAtom)
{
	//
	if(!mSingleH_DataInitialized)
		setDataForCpuEngineCallsForSingleHkl();
	//
	HansenCoppens_SF_Engine engine;
	Vector3d h_cart;
	mReciprocalUnitCell.fractionalToCartesian(hkl, h_cart);
	int nAtoms = localCoordinateSystems.size();
	f.resize(nAtoms);
	int i, nTypes = mType2WfnType.size();
	int nWfnTypes = mWfnParameters.size();
	vector<double> f_sph(nTypes), f_core(nWfnTypes);
	double h = sqrt(h_cart * h_cart);
	for (i = 0; i < nWfnTypes; i++)
		f_core[i] = mF_CoreInterpolators[i](h);
	for (i = 0; i < nTypes; i++)
		f_sph[i] = f_core[mType2WfnType[i]] + mF_SphericalValenceInterpolators[i](h);

	engine.calculateFormFactors(mWfnParameters, mTypeParameters, f_sph, mAtomToWfnTypeMap, mAtomToAtomTypeMap,
		                        localCoordinateSystems, h_cart, f, countAtom, mType2WfnType, mDefValSlaterNormalization, mTypeMaxL);
	//vector<double> 
}

void HansenCoppensStructureFactorCalculator::calculateFormFactorsCart(
    const std::vector<Matrix3d>& localCoordinateSystems,
    const Vector3d& h_cart,
    std::vector<std::complex<double> >& f,
    const std::vector<bool>& countAtom)
{
    //
    if (!mSingleH_DataInitialized)
        setDataForCpuEngineCallsForSingleHkl();
    //
    HansenCoppens_SF_Engine engine;
    //Vector3d h_cart;
    //mReciprocalUnitCell.fractionalToCartesian(hkl, h_cart);
    int nAtoms = localCoordinateSystems.size();
    f.resize(nAtoms);
    int i, nTypes = mType2WfnType.size();
    int nWfnTypes = mWfnParameters.size();
    vector<double> f_sph(nTypes), f_core(nWfnTypes);
    double h = sqrt(h_cart * h_cart);
    for (i = 0; i < nWfnTypes; i++)
        f_core[i] = mF_CoreInterpolators[i](h);
    for (i = 0; i < nTypes; i++)
        f_sph[i] = f_core[mType2WfnType[i]] + mF_SphericalValenceInterpolators[i](h);

    engine.calculateFormFactors(mWfnParameters, mTypeParameters, f_sph, mAtomToWfnTypeMap, mAtomToAtomTypeMap,
        localCoordinateSystems, h_cart, f, countAtom, mType2WfnType, mDefValSlaterNormalization, mTypeMaxL);
    //vector<double> 
}


 

} // namespace discamb
