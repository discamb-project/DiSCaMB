#pragma once


#include "discamb/Scattering/SfCalculator.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"


namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


	class TscFileBasedSfCalculator : public SfCalculator {

		AnyScattererStructureFactorCalculator* mCalculator;
		TscFileBasedSfCalculator();
	public:

		TscFileBasedSfCalculator(const Crystal& crystal, const nlohmann::json& data);

		virtual ~TscFileBasedSfCalculator();
        
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const {};
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




	};
    /** @}*/
}
