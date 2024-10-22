#ifndef _DISCAMB_SCATTERING_HansenCoppensStructureFactorCalculatorDev_H_
#define _DISCAMB_SCATTERING_HansenCoppensStructureFactorCalculatorDev_H_

#include "SF_Engine_DataTypes.h"
#include "SF_CalcDataTypes.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
#include "discamb/BasicUtilities/Timer.h"
#include "HansenCoppens_SF_Engine.h"
#include "scattering_utilities.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/CrystalStructure/StructuralParametersConverter.h"
#include "discamb/MathUtilities/NaturalCubicSpline.h"



#include <vector>
#include <map>
#include <complex>
#include <sstream>



namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */

    /** 
    \brief Hansen Coppens model structure factors calculation. 
    */

class HansenCoppensStructureFactorCalculatorDev
{
    HansenCoppensStructureFactorCalculatorDev();

public:
    /** Calculate structure factors and their derivatives*/
    HansenCoppensStructureFactorCalculatorDev(const Crystal &crystal, const HC_ModelParametersDev &parameters);
    
    ~HansenCoppensStructureFactorCalculatorDev();

    /** Sets number of threads to be used in the case of CPU calculations,
        does not effect GPU calculation or CPU when compiled with compilers with no openMP support
    */
    void setN_Threads(int n);
    
    /** Sets model of the electron density (note that the model is also set in constructor).*/
    void setModel(const Crystal &crystal, const HC_ModelParametersDev &parameters);

    

    /**
    Sets convention for derivatives calculation:
    \param dXyzConvention specifies if structure factor derivatives with respect to coordinates should be calculated for
                          fractional or cartesian coordinates.
    \param dAdpConvention specifies if structure factor derivatives with respect to aniosotropic ADPs should be 
                          calculated for \f$U_{cart}\f$, \f$U_{cif}\f$ or \f$U^{*}\f$

    */
    
    void setDerivativesConvention(structural_parameters_convention::XyzCoordinateSystem dXyzConvention, 
                                  structural_parameters_convention::AdpConvention dAdpConvention);

    /**
    Gets convention for derivatives calculation.
    \sa setDerivativesConvention
    */
    void getDerivativesConvention(structural_parameters_convention::XyzCoordinateSystem &dXyzConvention,
                                  structural_parameters_convention::AdpConvention &dAdpConvention) const;


    /** Main part of calculations will be run on CPU. */
    void useCPU();
    /** Main part of calculations will be run on GPU. */
    void useGPU();
    /** \brief Main part of calculations will be run on CPU with use of spherical atomic form factors.
       
       Parameterization published by D. Waasmaier & A. Kirfel \cite Waasmaier will be applied 
       and SDS fit from CCTBX \cite cctbx for hydrogen.
    */
    void useCPU_IAM();
    // CPU, GPU or CPU_IAM
    void setCalculationMode(const std::string &engineType);
    
    /** \brief Calculates structure factors.

    \param atoms carry atomic structural parameters
    \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to   
                                  normalized x, y and z directions in the new coordinate system
    \param hkl Miller indices 
    \param f the resulting structure factors
    \param countAtomContribution flags for counting atom contribution to structure factor calculations, size should match 
                                 number of atoms
    */
    void calculateStructureFactors( const std::vector<AtomInCrystal> &atoms,
                                    const std::vector<Matrix3d> &localCoordinateSystems,
                                    const std::vector<Vector3i> &hkl,
                                    std::vector<std::complex<double> > &f,
                                    const std::vector<bool> &countAtomContribution);

    /** \brief Calculates structure factors.

    \param atoms carry atomic structural parameters
    \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
    normalized x, y and z directions in the new coordinate system
    \param hkl Miller indices
    \param f the resulting structure factors
    */

    void calculateStructureFactors(
             const std::vector<AtomInCrystal> &atoms,
             const std::vector<Matrix3d> &localCoordinateSystems,
             const std::vector<Vector3i> &hkl,
             std::vector<std::complex<double> > &f) ;

    /** \brief Calculates structure factors and derivatives of target function with respet to structural parameters.

    \param atoms carry atomic structural parameters (positions, ADPs, occupancy)
    \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
    normalized x, y and z directions in the new coordinate system
    \param hkl Miller indices
    \param f the resulting structure factors
    \param dTarget_dparam derivatives of target function with respet to structural parameters 
    \param dTarget_df derivatives of target function \f$ T \f$ with respect to structure factors \f$ F_k = A_k + iB_k \f$, 
                      \f$A_k, B_k \in \mathbb{R}\f$ given as: 
                      \f$ \frac{\partial T}{\partial A_k} + i \frac{\partial T}{\partial B_k}\f$
                      
    \param countAtomContribution flags for counting atom contribution to structure factor calculations, size should match
    number of atoms
    */
    void calculateStructureFactorsAndDerivatives( 
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Matrix3d> &localCoordinateSystems,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution);

	/**

	*/

	void calculateFormFactors( 
		const std::vector<Matrix3d>& localCoordinateSystems,
		const Vector3i& hkl,
		std::vector<std::complex<double> >& f,
		const std::vector<bool>& countAtom);

    void calculateFormFactorsCart(
        const std::vector<Matrix3d>& localCoordinateSystems,
        const Vector3d& hkl,
        std::vector<std::complex<double> >& f,
        const std::vector<bool>& countAtom);


    /** 
        \brief Calculates structure factors and their derivatives with respect to structural parameters. 
               
        For each hkl vector an user provided object \p singleHklResultsProcessor function is called 
        which collects the results for given h vector.
               
        \param atoms carry atomic structural parameters (positions, ADPs, occupancy)
        \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
        normalized x, y and z directions in the new coordinate system
        \param hkl Miller indices
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
    void calculateStructureFactorsAndDerivatives(
                 const std::vector<AtomInCrystal> &atoms,
                 const std::vector<Matrix3d> &localCoordinateSystems, 
                 const std::vector<Vector3i> &hkl,
                 T &singleHklResultsProcessor) const;

    /**
    \brief Calculates structure factors and their derivatives with respect to structural parameters.


    \param atoms carry atomic structural parameters (positions, ADPs, occupancy)
    \param localCoordinateSystems atomic local systems of coordinates, columns of the matrices should corresspond to
    normalized x, y and z directions in the new coordinate system
    \param hkl Miller indices
    \param singleHklResultsProcessor object of user defined type intended for collecting results,
    for calculation for each Miller index its member function:
    \code{.cpp}
    void onPerHklCalculation(
    int hkl_idx,
    std::complex<double> &structureFactor,
    discamb::SfDerivativesAtHkl &derivatives);
    \endcode
    is called (hkl_idx starts from 0). See e.g. Tests/correctness/test_1/GradientDataCollector.h for an example implementation.
    \param countAtomContribution (if present) flags for counting atom contribution to structure factor calculations, size should match
    number of atoms

    */


    template<typename T>
    void calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Matrix3d> &localCoordinateSystems,
        const std::vector<Vector3i> &hkl, 
        T &singleHklResultsProcessor,
        const std::vector<bool> &countAtomContribution) const;
    /** \brief Information on execution time of last structure factor calculation.*/
    void getCalculationTime(double &cpuTime,double &wallClockTime) const;

    /** \brief Defined for debug purposes. */

    void saveRawData(const Crystal &crystal, const std::vector<Matrix3d> &localCoordinateSystems, const std::vector<Vector3i> &hkl,
                     const std::vector<std::complex<double> > &dTarget_df, const std::vector<bool> &countAtomContribution,
                     std::ostream &output);

private:

    enum HC_SF_EngineId
    {            
        CPU = 0, 
        GPU = 1, 
        CPU_IAM = 2 
    };

	/* some data used by HansenCoppens_SF_Engine, needed for an efficient atomic form factor calculation for one h*/
	void setDataForCpuEngineCallsForSingleHkl();

	bool mSingleH_DataInitialized;
	std::vector< NaturalCubicSpline > mF_CoreInterpolators, mF_SphericalValenceInterpolators;
	std::vector<int> mType2WfnType;
	std::vector<std::vector<REAL> > mDefValSlaterNormalization;
	std::vector<int> mTypeMaxL;

	/* end of some data used by HansenCoppens_SF_Engine.. */

    void setEngineType(HC_SF_EngineId engineId);
    std::vector<std::string> mWfnTypeLabels;
    mutable CpuTimer mCpuTimer;
    mutable WallClockTimer mWallClockTimer;
    mutable double mWallClockTime, mCpuTime;
    int mN_Threads;
    mutable bool mPlmWfnNormalized;
    UnitCell mUnitCell;
	ReciprocalLatticeUnitCell mReciprocalUnitCell;
    SpaceGroup mSpaceGroup;
    void setEngineSpecificSpaceGroupRepresentation();
    structural_parameters_convention::AdpConvention mD_AdpConvention,mAdpInputType;
    structural_parameters_convention::XyzCoordinateSystem mD_XyzConvention, mXyzInputCoordinateSystem;
    StructuralParametersConverter mConverter;

    void test(const Crystal &crystal, const std::vector<Matrix3d> &localCoordinateSystems,
              const std::vector<Vector3i> &hkl, std::vector<std::complex<double> > &f,
              std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
              const std::vector<std::complex<double> > &dTarget_df);

    void calculate_using_GPU(const Crystal &crystal, const std::vector<Matrix3d> &localCoordinateSystems,
                             const std::vector<Vector3i> &hkl, std::vector<std::complex<double> > &f,
                             std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
                             const std::vector<std::complex<double> > &dTarget_df);

    void calculate_using_GPU_batched(
                             const Crystal &crystal, const std::vector<Matrix3d> &localCoordinateSystems,
                             const std::vector<Vector3i> &hkl, std::vector<std::complex<double> > &f,
                             std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
                             const std::vector<std::complex<double> > &dTarget_df);

    // unfolded atoms - atoms with multiple atom types are represented as multiple atoms

    std::vector<int> mUnfoldedAtom2Atom;
    std::vector<std::vector<int> > mAtomUnfoldedAtoms;

    mutable HC_SF_EngineId mEngineType;
    std::vector<sf_engine_data_types::HC_WfnParam> mWfnParameters;
    mutable std::vector<sf_engine_data_types::HC_TypeParam> mTypeParameters;


    std::vector<int> mUnfoldedAtomToWfnTypeMap;
    std::vector<int> mUnfoldedAtomToAtomTypeMap;
    std::vector < std::vector <double> > mAtomConfigurationWeight;
    
    mutable std::vector<Vector3d> mUnfoldedAtomicPositions;
    mutable bool mWithGradients;
    mutable std::vector<std::vector<double> > mUnfoldedAtomic_displacement_parameters;
    mutable std::vector<double> mUnfoldedAtomicOccupancy;
    mutable std::vector<double> mUnfoldedAtomicMultiplicityFactor;
    mutable std::vector<double> mUnfoldedAtomicMultiplicityWeights;
    std::vector<sf_engine_data_types::SymmetryOperation> mSymmetryOperations;
    bool mIsCentrosymmetric;
    Vector3d mInversionCenterTranslation;
    mutable std::vector<Vector3d> mHKL_Cartesian;
    Vector3d mA_Star,mB_Star,mC_Star;

    static double slaterNormalizationFactor(int powerR,double exponent);
    static void addOrbitalDensity(const SlaterTypeAtomicOrbitalRdf &orbital,int occupancy,std::vector<double> &density_coeff,
                                  std::vector<int> &density_powers,std::vector<double> &density_exponents);
    //static void copyVector(const std::vector<int> &from,std::vector<int> &to);
    

    void prepareDataForEngine(const std::vector<AtomInCrystal> &atoms, const std::vector<Vector3i> &hkl,
                              std::vector<std::complex<double> > &f,
                              std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam) const;
    
    void prepareDataForEngine_part(const std::vector<AtomInCrystal> &atoms, const std::vector<Vector3i> &hkl) const;

    void convertDerivatives(std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam) const;
        
    template<typename T>
    class WrapperForSingleHklResultProcessor
    {
    private:
        T &mSingleHklResultProcessor;
        const StructuralParametersConverter &mConverter;
        SfDerivativesAtHkl mConvertedDerivatives;
        std::vector<bool> mU_Aniso;
        structural_parameters_convention::XyzCoordinateSystem mXyzIn, mXyzOut;
        structural_parameters_convention::AdpConvention mAdpIn, mAdpOut;
        void convert(const SfDerivativesAtHkl &derivatives);

    public:
        WrapperForSingleHklResultProcessor(T &t, const StructuralParametersConverter &c, const std::vector<int> &nAdpComponents);

        ~WrapperForSingleHklResultProcessor();

        void setConventions(
             structural_parameters_convention::XyzCoordinateSystem xyzIn,
             structural_parameters_convention::XyzCoordinateSystem xyzOut,
             structural_parameters_convention::AdpConvention adpIn,
             structural_parameters_convention::AdpConvention adpOut);

        void onPerHklCalculation(int hkl_idx, std::complex<double> &structureFactor, SfDerivativesAtHkl &derivatives);
    };

    void saveRawData(const Crystal &crystal, const std::vector<Matrix3d> &localCoordinateSystems,
                     const std::vector<Vector3i> &hkl, std::ostream &output);


};

//-------------------------------
//
//    TEMPLATE IMPLEMENTATION
//
//-------------------------------


template<typename T>
void HansenCoppensStructureFactorCalculatorDev::calculateStructureFactorsAndDerivatives(
    const std::vector<AtomInCrystal> &atoms,
    const std::vector<Matrix3d> &localCoordinateSystems,
    const std::vector<Vector3i> &hkl, 
    T &singleHklResultsProcessor) 
const
{
    // set wrapper for singleHklResultsProcessor enabling required
    // coordinate system choice and type of ADP parameterization

    int atomIdx, nAtoms = atoms.size();
    std::vector<int> nAdpComponents(nAtoms);
    for(atomIdx=0; atomIdx<nAtoms; atomIdx++)
        nAdpComponents[atomIdx] = atoms[atomIdx].adp.size();
        
    WrapperForSingleHklResultProcessor<T> wrapper(singleHklResultsProcessor, mConverter, nAdpComponents);
                                                                                         
    wrapper.setConventions( structural_parameters_convention::XyzCoordinateSystem::cartesian, 
                            mD_XyzConvention,
                            structural_parameters_convention::AdpConvention::U_cart,
                            mD_AdpConvention);
    
    // set calculations for all atoms
    std::vector<bool> countAtomContribution(atoms.size(),true);

    // calculate

    calculateStructureFactorsAndDerivatives(
        atoms,
        localCoordinateSystems,
        hkl,
        wrapper,
        countAtomContribution);

}

template<typename T>
void HansenCoppensStructureFactorCalculatorDev::calculateStructureFactorsAndDerivatives(
    const std::vector<AtomInCrystal> &atoms,
    const std::vector<Matrix3d> &localCoordinateSystems,
    const std::vector<Vector3i> &hkl, 
    T &singleHklResultsProcessor,
    const std::vector<bool> &countAtomContribution) 
const
{
    HansenCoppens_SF_Engine engine;
    prepareDataForEngine_part(atoms, hkl);
    std::vector<double> centeringMultipliers;
    char latticeCentering;
    bool obverse;
    latticeCentering = mSpaceGroup.centering(obverse);
    scattering_utilities::centeringSfMultipliers(latticeCentering, hkl, centeringMultipliers, obverse);
    Vector3<CrystallographicRational> inversionCenterTranslationFractional;
    Vector3d inversionCenterTranslation;
    mSpaceGroup.inversionCenterTranslation(inversionCenterTranslationFractional);
    mUnitCell.fractionalToCartesian(inversionCenterTranslationFractional, inversionCenterTranslation); 
    if(mEngineType!=CPU_IAM)
    {
        mCpuTimer.start();
        mWallClockTimer.start();

        engine.calculateSF(mWfnParameters, mTypeParameters, mAtomToWfnTypeMap, mAtomToAtomTypeMap, mAtomicPositions,
                           mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
                           localCoordinateSystems, mSymmetryOperations, mSpaceGroup.isCentrosymmetric(),
                           inversionCenterTranslation, mHKL_Cartesian,
                           centeringMultipliers, countAtomContribution, singleHklResultsProcessor);
        mCpuTime = mCpuTimer.stop();
        mWallClockTime = mWallClockTimer.stop();

    }
    else
    {
        int wfnTypeIdx, nWfnTypes = mWfnParameters.size();
        std::vector<std::complex<double> > anomalousScattering(nWfnTypes);

        for (wfnTypeIdx = 0; wfnTypeIdx<nWfnTypes; wfnTypeIdx++)
            anomalousScattering[wfnTypeIdx] = mWfnParameters[wfnTypeIdx].anomalous_scattering;

        mCpuTimer.start();
        mWallClockTimer.start();

        engine.calculateSF_IAM(mWfnTypeLabels, anomalousScattering, mAtomToWfnTypeMap, mAtomicPositions,
                               mAtomic_displacement_parameters, mAtomicOccupancy, mAtomicMultiplicityWeights,
                               mSymmetryOperations, mSpaceGroup.isCentrosymmetric(), 
                               inversionCenterTranslation, mHKL_Cartesian, centeringMultipliers,
                               countAtomContribution, singleHklResultsProcessor);
        mCpuTime = mCpuTimer.stop();
        mWallClockTime = mWallClockTimer.stop();
    }

    
}

template<typename T>
HansenCoppensStructureFactorCalculatorDev::WrapperForSingleHklResultProcessor<T>::WrapperForSingleHklResultProcessor(
    T &t, 
    const StructuralParametersConverter &c, 
    const std::vector<int> &nAdpComponents) : mSingleHklResultProcessor(t), mConverter(c)
{
    int atomIdx, nAtoms = nAdpComponents.size();
    mConvertedDerivatives.atomicPostionDerivatives.resize(nAtoms);
    mConvertedDerivatives.adpDerivatives.resize(nAtoms);
    mConvertedDerivatives.occupancyDerivatives.resize(nAtoms);
    mU_Aniso.resize(nAtoms);
    for(atomIdx=0; atomIdx<nAtoms; atomIdx++)
    {
        mConvertedDerivatives.adpDerivatives[atomIdx].resize(nAdpComponents[atomIdx]);
        mU_Aniso[atomIdx] = ( 6 == nAdpComponents[atomIdx] );
    }
    mXyzIn = structural_parameters_convention::XyzCoordinateSystem::cartesian;
    mXyzOut = structural_parameters_convention::XyzCoordinateSystem::fractional;
    mAdpIn = structural_parameters_convention::AdpConvention::U_cart;
    mAdpOut = structural_parameters_convention::AdpConvention::U_cif;
}

template<typename T>
HansenCoppensStructureFactorCalculatorDev::WrapperForSingleHklResultProcessor<T>::~WrapperForSingleHklResultProcessor()
{}

template<typename T>
void HansenCoppensStructureFactorCalculatorDev::WrapperForSingleHklResultProcessor<T>::setConventions(
structural_parameters_convention::XyzCoordinateSystem xyzIn,
structural_parameters_convention::XyzCoordinateSystem xyzOut,
structural_parameters_convention::AdpConvention adpIn,
structural_parameters_convention::AdpConvention adpOut)
{
    mXyzIn = xyzIn;
    mXyzOut = xyzOut;
    mAdpIn = adpIn;
    mAdpOut = adpOut;
}

template<typename T>
void HansenCoppensStructureFactorCalculatorDev::WrapperForSingleHklResultProcessor<T>::convert(
    const SfDerivativesAtHkl &derivatives)
{
    int atomIdx, nAtoms = mU_Aniso.size();
    for(atomIdx=0;atomIdx<nAtoms;atomIdx++)
    {
        mConverter.convertDerivativesXyz(derivatives.atomicPostionDerivatives[atomIdx],
                                         mConvertedDerivatives.atomicPostionDerivatives[atomIdx],
                                         mXyzIn, mXyzOut);
        if(mU_Aniso[atomIdx])
            mConverter.convertDerivativesADP(derivatives.adpDerivatives[atomIdx], 
                                             mConvertedDerivatives.adpDerivatives[atomIdx],
                                             mAdpIn, mAdpOut);
        else
            mConvertedDerivatives.adpDerivatives[atomIdx] = derivatives.adpDerivatives[atomIdx];
    }
    mConvertedDerivatives.occupancyDerivatives = derivatives.occupancyDerivatives;
}

template<typename T>
void HansenCoppensStructureFactorCalculatorDev::WrapperForSingleHklResultProcessor<T>::onPerHklCalculation(
    int hkl_idx, 
    std::complex<double> &structureFactor, 
    SfDerivativesAtHkl &derivatives)
{
    convert(derivatives);
    mSingleHklResultProcessor.onPerHklCalculation(hkl_idx, structureFactor, mConvertedDerivatives);
}
/** @}*/
} // namespace discamb

#endif /*_DISCAMB_SCATTERING_HansenCoppensStructureFactorCalculatorDev_H_*/
