#pragma once
#include "SfCalculator.h"

#include <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class CombinedStructureFactorCalculator : public SfCalculator
    {
    public:
        CombinedStructureFactorCalculator(
            std::vector<std::shared_ptr<SfCalculator> > &calculators,
            std::vector<std::vector<int> > &atoms);
        CombinedStructureFactorCalculator(const Crystal &crystal, const nlohmann::json &data);
        virtual ~CombinedStructureFactorCalculator();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const {};
        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous);
        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution);

        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& dTarget_df,
            const std::vector<bool>& countAtomContribution,
            const DerivativesSelector &derivativesSelector);

        virtual void calculateStructureFactors(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            const std::vector<bool> &countAtomContribution);

        virtual void update(const std::vector<AtomInCrystal> &atoms);

        virtual void calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &scatteringFactor,
            SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution);

        void set(std::vector<std::shared_ptr<SfCalculator> > &calculators,
            std::vector<std::vector<int> > &atoms);

		virtual void calculateFormFactors(const Vector3i& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual void calculateFormFactorsCart(const Vector3d& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;

    private:
        std::vector<std::shared_ptr<SfCalculator> > mCalculators;
        // [i][j] - is j-th atom included in i-th calculator calculations
        std::vector<std::vector<bool> > mCalculatorsAtomsInclusion;
        // [i][j] index of j-th atom which contribution os calculated with i-th
        std::vector<std::vector<int> > mCalculatorsAtoms;
        std::vector<std::complex<double> > mAnomalous;
        std::vector<std::complex<double> > mPartialScatteringFactorOneHkl;
        std::vector< std::vector<std::complex<double> > > mPartialScatteringFactorMultiHkl;
        std::vector< std::vector<TargetFunctionAtomicParamDerivatives> > mPartial_dTarget_dparam;
        std::vector<SfDerivativesAtHkl> mPartialDerivatives;
        

    };
    /** @}*/
}
