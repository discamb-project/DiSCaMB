#include "discamb/Scattering/NumericalDensityIamSfCalculator.h"
#include "discamb/Scattering/NumericalDensityIamFfCalculator.h"

using namespace std;

namespace discamb {
    
    NumericalDensityIamSfCalculator::NumericalDensityIamSfCalculator() {};

    NumericalDensityIamSfCalculator::NumericalDensityIamSfCalculator(
        const Crystal& crystal,
        const nlohmann::json& data)
    {
        //std::shared_ptr<AnyScattererStructureFactorCalculator2> mSfCalculator;
        mSfCalculator = make_shared<AnyScattererStructureFactorCalculator2>(crystal);
        mFormFactorsCalculator = shared_ptr<NumericalDensityIamFfCalculator>(
            new NumericalDensityIamFfCalculator(crystal, data));
        mSfCalculator->setAtomicFormfactorManager(mFormFactorsCalculator);
    }

    NumericalDensityIamSfCalculator::~NumericalDensityIamSfCalculator() {};

    void NumericalDensityIamSfCalculator::getModelInformation(
        std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        
    }

    void NumericalDensityIamSfCalculator::setAnomalous(
        const std::vector<std::complex<double> >& anomalous)
    {
    }

    void NumericalDensityIamSfCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution)
    {
        mSfCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
    }

    void NumericalDensityIamSfCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        const std::vector<bool>& countAtomContribution)
    {
        mSfCalculator->calculateStructureFactors(hkl, f, countAtomContribution);
    }

    void NumericalDensityIamSfCalculator::update(
        const std::vector<AtomInCrystal>& atoms)
    {

    }

    void NumericalDensityIamSfCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i& hkl,
        std::complex<double>& scatteringFactor,
        discamb::SfDerivativesAtHkl& derivatives,
        const std::vector<bool>& countAtomContribution)
    {
        mSfCalculator->calculateStructureFactorsAndDerivatives(
            hkl, 
            scatteringFactor, 
            derivatives, 
            countAtomContribution);
    }

        // formFactors[hkl idx][atom idx]
    void NumericalDensityIamSfCalculator::calculateFormFactors(
        const Vector3i& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mSfCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }

        // formFactors[hkl idx][atom idx]
    void NumericalDensityIamSfCalculator::calculateFormFactorsCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom) const
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }

    


    

}

