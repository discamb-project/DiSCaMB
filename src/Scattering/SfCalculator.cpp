
#include "discamb/Scattering/SfCalculator.h"

#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/Scattering/AnyHcCalculator.h"
#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
#include "discamb/Scattering/TscFileBasedSfCalculator.h"
#include "discamb/Scattering/StockholderAtomSfCalculator.h"
#include "discamb/Scattering/StockholderAtomBankSfCalculator.h"


#include <map>

using namespace std;

namespace discamb {

	SfCalculator::SfCalculator()
	{
		mN_threads = 1;
	}
	
	void SfCalculator::setN_threads(
		int n)
	{
		mN_threads = n;
	}

    void SfCalculator::add_to_log(
        int n, 
        const string &s)
    {
        //ofstream out("discamb_log", std::ofstream::out | std::ofstream::app);
        //out << n << " " << s << endl;
        //out.close();
    }

    void SfCalculator::modelNames(
        std::map<std::string, std::set<std::string> >& names)
    {
        names = std::map<std::string, std::set<std::string> >(
            {
                {"iam",{"iam"}},
                //          {"legacy_har", {"legacy_har"}},
                {"gar", { "har", "multi_har", "gar"}},
                {"har bank",{"har bank"}},
                {"hc",{"hc"}},
                {"sab",{"sab"}},
                {"tsc",{"tsc"}},
                {"matts",{"matts", "taam", "ubdb"}},
                {"tham",{"tham"}}
             });
        
    }

    bool SfCalculator::equivalentModelNames(
        const std::string& name1,
        const std::string& name2)
    {
        std::map<std::string, std::set<std::string> > names;
        modelNames(names);


        for (auto& p : names)
            if (find(p.second.begin(), p.second.end(), name1) != p.second.end())
                if (find(p.second.begin(), p.second.end(), name2) != p.second.end())
                    return true;
        return false;
    }

    std::string SfCalculator::uniqueModelName(const std::string& _name)
    {
        string name;
        string_utilities::toLower(_name, name);

        std::map<std::string, std::set<std::string> > names;
        modelNames(names);

        for (auto& nameSet : names)
            if (nameSet.second.find(name) != nameSet.second.end())
                return nameSet.first;
        
        return string();
    }

    void SfCalculator::jsonCurrentFormat(
        const nlohmann::json& inputData,
        std::string& modelName,
        nlohmann::json& modelDetails)
    {
        if (inputData.find("form factor engine") != inputData.end())
        {
            modelName = getFromJSON<string>(inputData["form factor engine"], "type");
            modelDetails = getFromJSON<nlohmann::json>(inputData["form factor engine"], "data");
        }
        else
        {
            if (inputData.find("model") != inputData.end())
            {
                modelName = inputData["model"].get<string>();
                modelDetails = inputData;
            }
            else
            {
                if (inputData.find("type") != inputData.end())
                {
                    modelName = inputData["type"].get<string>();
                    modelDetails = inputData["data"];
                }
                else
                    on_error::throwException("incorrect json data passed encountered when trying to create form factors calculator", __FILE__, __LINE__);
            }
        }

    }

    void SfCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution,
        const DerivativesSelector& derivativesSwitch)
    {
        calculateStructureFactorsAndDerivatives(atoms, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
    }

    void SfCalculator::calculateFormFactorsCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }

    SfCalculator *SfCalculator::create(
        const Crystal &crystal,
        const nlohmann::json &data)
    {
        string type;
        nlohmann::json engineData;

        jsonCurrentFormat(data, type, engineData);

        add_to_log(__LINE__, string(__FILE__) + string(" : engine type - ") + type);              
        
        type = uniqueModelName(type);

        if (string("iam") == type)
            return new AnyIamCalculator(crystal, engineData);
        if (string("gar") == type)
            return new StockholderAtomSfCalculator(crystal, engineData);
        if (string("hc") == type)
            return new AnyHcCalculator(crystal, engineData);
		if (string("tsc") == type)
			return new TscFileBasedSfCalculator(crystal, engineData);
        if (string("matts") == type)
            return new HcAtomBankStructureFactorCalculator(crystal, engineData);
        if (string("tham") == type)
            return new StockholderAtomBankSfCalculator(crystal, engineData);
        string errorMessage =
            string("Invalid argument when creating SfCalculator - unknown type of the calculator: '") + type +
            string("', allowed values: iam, legacy_har, gar, har, multi_har, sab, tsc, matts, taam and ubdb");
        
        on_error::throwException(errorMessage, __FILE__, __LINE__);

        return NULL;
    }

    void SfCalculator::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        int hklIdx, nHkl = hkl.size();
        vector<complex<double> > oneH_FormFactors;
        formFactors.clear();
        formFactors.resize(nHkl);
        //cout << "calculate form factors" << endl;
        int procent = 0;

        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            calculateFormFactors(hkl[hklIdx], oneH_FormFactors, includeAtom);
            formFactors[hklIdx] = oneH_FormFactors;
            //if((hklIdx+1)%100 == 0)
              //  cout << hklIdx + 1 << " / " << nHkl << endl;
            //newProcent = (hklIdx + 1) * 100 / nHkl;
            //if (newProcent > procent)
            //{
            //    procent = newProcent;
            //    //cout << "\r" << procent;
            //    cout << procent << endl;
            //}
        }
        //cout << "\ndone" << endl;
    }

}

