#pragma once

#include "SfCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"

#include <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class AnyIamCalculator : public SfCalculator
    {
        Crystal mCrystal;
        //AnyScattererStructureFactorCalculator2 *mCalculator2;
        //AnyScattererStructureFactorCalculator* mCalculator;
        std::shared_ptr < AnyScattererStructureFactorCalculator2> mCalculator2;
        std::shared_ptr < AnyScattererStructureFactorCalculator> mCalculator;
		std::shared_ptr<AtomicFormFactorCalculationsManager> mManager;
        void set(const Crystal &crystal, bool electronScattering, const std::string& table, const std::string &algorithm);
        std::vector<std::pair<std::string, std::string> > mModelInfo;
        bool mUseLineAlgorithm=false;
        
    public:
        AnyIamCalculator(const Crystal &crystal, bool electronScattering = false, const std::string & table = std::string());
        AnyIamCalculator(const Crystal &crystal, const nlohmann::json &data);
        virtual ~AnyIamCalculator();

        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous);

        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;

        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution);

        virtual void calculateStructureFactors(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            const std::vector<bool> &countAtomContribution);


        virtual void update(const std::vector<AtomInCrystal> &atoms);

        virtual void calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &scatteringFactor,
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution);

		virtual void calculateFormFactors(
			const Vector3i& hkl,
			std::vector<std::complex<double> >& formFactors,
			const std::vector<bool>& includeAtom) const;


    };

    /** @}*/
}
