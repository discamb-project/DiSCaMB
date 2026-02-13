#include "discamb/Scattering/SfCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"

namespace discamb {

    class NumericalDensityIamSfCalculator: public SfCalculator
    {
        NumericalDensityIamSfCalculator();
        std::shared_ptr<AnyScattererStructureFactorCalculator2> mSfCalculator;
        std::shared_ptr<AtomicFormFactorCalculationsManager> mFormFactorsCalculator;
    public:
        NumericalDensityIamSfCalculator(const Crystal& crystal, const nlohmann::json& data);
        ~NumericalDensityIamSfCalculator();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;

        virtual void setAnomalous(const std::vector<std::complex<double> >& anomalous);

        virtual void calculateStructureFactorsAndDerivatives(
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

        // formFactors[hkl idx][atom idx]
        virtual void calculateFormFactors(const Vector3i& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;

        // formFactors[hkl idx][atom idx]
        virtual void calculateFormFactorsCart(const Vector3d& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;

        virtual std::string name() const { return std::string("numerical IAM"); }
                

    };

}
