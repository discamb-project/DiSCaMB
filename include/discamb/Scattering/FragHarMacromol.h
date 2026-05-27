#pragma once

#include "discamb/Scattering/StockholderAtomSfCalculator.h"
#include "discamb/StructuralProperties/MacromolecularStructuralInformation.h"

#include <memory>
#include <fstream>
#include <map>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */



    class FragHarMacromol : public SfCalculator
    {

    public:

        FragHarMacromol(
            const Crystal& crystal, 
            const HirshfeldAtomModelSettings& settings, 
            const MacromolecularStructuralInformation &macromolInfo,
            bool electronScattering = false, 
            const std::string& jobname = "ham");

        FragHarMacromol(const Crystal& crystal, const nlohmann::json& data);

        virtual ~FragHarMacromol();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;
        virtual std::string name() const { return "frag_har_macromol"; }
        virtual void setAnomalous(const std::vector<std::complex<double> >& anomalous);


        // ignores count atom contrib if iam atoms present
        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& dTarget_df,
            const std::vector<bool>& countAtomContribution);

        // ignores count atom contrib if iam atoms present
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
        std::unique_ptr<StockholderAtomSfCalculator> mSfCalculator;
        void set(const Crystal& crystal, const HirshfeldAtomModelSettings& settings, const MacromolecularStructuralInformation& macromolInfo, bool electronScattering = false, const std::string& jobname = "ham");
        void make_fragments(
            const Crystal &crystal, 
            const MacromolecularStructuralInformation& macromolInfo, 
            std::vector<QmFragmentInCrystal> &fragments,
            std::vector<std::vector<AtomRepresentativeInfo> > &representatives);
    };
    /** @}*/
}
