#ifndef _DISCAMB_EXAMPLES_REFINE_SFCALCULATOR_H_
#define _DISCAMB_EXAMPLES_REFINE_SFCALCULATOR_H_

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/config.h"
#include "discamb/Scattering/SF_CalcDataTypes.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "json.hpp"

#include <set>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class SfCalculator {
        
    public:
        
		SfCalculator();

        virtual ~SfCalculator() {};

		virtual void setN_threads(int n);

        static SfCalculator *create(const Crystal &crystal,
                                    const nlohmann::json &data = nlohmann::json());

        static SfCalculator* create(const Crystal& crystal, const std::string jsonFileName);


        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const = 0;

        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous) =0;

        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& dTarget_df,
            const std::vector<bool>& countAtomContribution) = 0;

        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution,
            const DerivativesSelector &selector);

        virtual void calculateStructureFactors(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            const std::vector<bool> &countAtomContribution)=0;

        virtual void update(const std::vector<AtomInCrystal> &atoms)=0;

        virtual void calculateStructureFactorsAndDerivatives(
                const Vector3i &hkl,
                std::complex<double> &scatteringFactor,
                discamb::SfDerivativesAtHkl &derivatives,
                const std::vector<bool> &countAtomContribution)=0;

	 	virtual void calculateFormFactors(const Vector3i& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const = 0;
        virtual void calculateFormFactors(const std::vector<Vector3i> &hkl, std::vector< std::vector<std::complex<double> > >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual void calculateFormFactorsCart(const Vector3d& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;

        virtual std::string name() const { return std::string("undefined"); }
        virtual void allowedNames(std::set<std::string>& names) const { names = { this->name() }; }

        static std::string uniqueModelName(const std::string &name);
        static void modelNames(std::map<std::string, std::set<std::string> >& names);
        static bool equivalentModelNames(const std::string& name1, const std::string& name2);
        static void jsonCurrentFormat(const nlohmann::json& inputData, std::string& modelName, nlohmann::json& modelDetails);
        
	protected:

		int mN_threads;
        template<typename T>
        static T getFromJSON(const nlohmann::json &data, const std::string &key)
        {
            if (data.find(key) == data.end())
            {
                std::string error_message = std::string("expected element '") + key +
                    std::string("' in JSON object defining scattering model");

                discamb::on_error::throwException(error_message, __FILE__, __LINE__);
            }
            return data.find(key).value().get<T>();
        }
        
        static void add_to_log(int n, const std::string &s = std::string());


    };

    /** @}*/

}


#endif
