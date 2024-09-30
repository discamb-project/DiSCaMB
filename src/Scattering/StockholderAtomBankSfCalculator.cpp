#include "discamb/Scattering/StockholderAtomBankSfCalculator.h"
#include "discamb/Scattering/StockholderAtomBankFormFactorCalculationManager.h"

using namespace std;

namespace discamb 
{

    StockholderAtomBankSfCalculator::StockholderAtomBankSfCalculator(
        const Crystal& crystal, 
        const nlohmann::json& data)
    {
        string bankFile;
        bankFile = data.value("bank file", bankFile);

        int nThreads = data.value("n threads", 1);
        nThreads = data.value("n cores", nThreads);

        bool useSphericalHarmonicsExpansion = true;
        useSphericalHarmonicsExpansion = data.value("spherical harmonics expansion", useSphericalHarmonicsExpansion);

        std::shared_ptr< AtomicFormFactorCalculationsManager> manager = std::shared_ptr< AtomicFormFactorCalculationsManager>
            (new StockholderAtomBankFormFactorCalculationManager(crystal, bankFile, nThreads, useSphericalHarmonicsExpansion));
        mCalculator = shared_ptr<AnyScattererStructureFactorCalculator>(new AnyScattererStructureFactorCalculator(crystal));
        mCalculator->setAtomicFormfactorManager(manager);
    }

    StockholderAtomBankSfCalculator::~StockholderAtomBankSfCalculator()
    {

    }


    void StockholderAtomBankSfCalculator::getModelInformation(
        std::vector<std::pair<std::string, std::string> >& modelInfo) 
        const
    {
        modelInfo.clear();
        modelInfo.push_back({ "SCATTERING MODEL", "THAM" });
    }

    void StockholderAtomBankSfCalculator::setAnomalous(const std::vector<std::complex<double> >& anomalous)
    {
        mCalculator->setAnoumalous(anomalous);
    }

    void StockholderAtomBankSfCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution)
    {
        mCalculator->update(atoms);
        mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
    }

    void StockholderAtomBankSfCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        const std::vector<bool>& countAtomContribution)
    {
        mCalculator->update(atoms);

    }

    void StockholderAtomBankSfCalculator::update(const std::vector<AtomInCrystal>& atoms)
    {
        mCalculator->update(atoms);
    }

    void StockholderAtomBankSfCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i& hkl,
        std::complex<double>& scatteringFactor,
        discamb::SfDerivativesAtHkl& derivatives,
        const std::vector<bool>& countAtomContribution)
    {
        mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives, countAtomContribution);
    }

    void StockholderAtomBankSfCalculator::calculateFormFactors(
        const Vector3i& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    void StockholderAtomBankSfCalculator::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    StockholderAtomBankSfCalculator::StockholderAtomBankSfCalculator()
    {

    }


}
