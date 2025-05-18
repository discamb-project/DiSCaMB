#pragma once

#include "SfCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/Scattering/HansenCoppensStructureFactorCalculator.h"

#include <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class AnyHcCalculator : public SfCalculator
    {
        Crystal mCrystal;
        AnyScattererStructureFactorCalculator *mCalculator;
		std::shared_ptr<AtomicFormFactorCalculationsManager> mManager;
		// introduced for large molecules 
		HansenCoppensStructureFactorCalculator* mHcCalculator;
		bool mUseImplementationForLargeMolecules;
		std::vector<std::shared_ptr<LocalCoordinateSystemInCrystal> > mLcs;
		std::vector<Matrix3d> mLcsMatrices;
        
		//
        void set(const Crystal &crystal,
            const HC_ModelParameters &parameters,
            const std::vector<std::shared_ptr<LocalCoordinateSystemInCrystal> > &lcs,
            bool electronScattering,
            bool integerChargeSpherical,
			bool implementationForLargeMolecules,
			int nThreads,
            bool frozenLcs);
    public:
        AnyHcCalculator(const Crystal &crystal,
            const HC_ModelParameters &parameters,
            const std::vector<std::shared_ptr<LocalCoordinateSystemInCrystal> > &lcs,
            bool electronScattering = false,
            bool neutralSpherical = false,
			bool implementationForLargeMolecules = false,
			int nThreads = 1,
            bool frozenLcs = false);
        AnyHcCalculator(const Crystal &crystal, const nlohmann::json &data);
        virtual ~AnyHcCalculator();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;
        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous);
		virtual void setN_threads(int n);

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
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution);

		virtual void calculateFormFactors(
			const Vector3i& hkl, 
			std::vector<std::complex<double> >& formFactors, 
			const std::vector<bool>& includeAtom) const;
    };
    /** @}*/
}