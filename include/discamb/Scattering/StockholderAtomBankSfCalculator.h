#pragma once

#include "SfCalculator.h"
#include "AnyScattererStructureFactorCalculator.h"
#include "AtomicFormFactorCalculationsManager.h"

#include <memory>

namespace discamb {


    class StockholderAtomBankSfCalculator : public SfCalculator
    {
    public:

        StockholderAtomBankSfCalculator(const Crystal& crystal, const nlohmann::json& data);
        ~StockholderAtomBankSfCalculator();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;

        virtual void setAnomalous(const std::vector<std::complex<double> >& anomalous);
        void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& dTarget_df,
            const std::vector<bool>& countAtomContribution);

        virtual void calculateStructureFactors(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            const std::vector<bool>& countAtomContribution);

        virtual void update(const std::vector<AtomInCrystal>& atoms);

        virtual void calculateStructureFactorsAndDerivatives(
            const Vector3i& hkl,
            std::complex<double>& scatteringFactor,
            discamb::SfDerivativesAtHkl& derivatives,
            const std::vector<bool>& countAtomContribution);

        virtual void calculateFormFactors(const Vector3i& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual void calculateFormFactors(const std::vector<Vector3i>& hkl, std::vector< std::vector<std::complex<double> > >& formFactors, const std::vector<bool>& includeAtom) const;
    private:
        StockholderAtomBankSfCalculator();
        std::shared_ptr<AnyScattererStructureFactorCalculator> mCalculator;
        //std::shared_ptr< AtomicFormFactorCalculationsManager> mManager;
    };

}

