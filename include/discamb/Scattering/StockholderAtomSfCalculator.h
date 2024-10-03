#include "discamb/CrystalStructure/UnitCellContent.h"

#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/SfCalculator.h"
#include "discamb/Scattering/StockholderAtomFormFactorCalcManager.h"                                 
//#include "discamb/Scattering/StockholderAtomSfCalculatorSettings.h"
//#include "discamb/Scattering/HirshfeldAtomModelSettings.h"

#include <memory>
#include <fstream>
#include <map>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


#ifndef SF_CALCULATOR_USE_CLUSTER
#define SF_CALCULATOR_USE_CLUSTER
#endif

class StockholderAtomSfCalculator : public SfCalculator
{

    public:

        StockholderAtomSfCalculator(
            const Crystal &crystal, const HirshfeldAtomModelSettings& settings, bool electronScattering=false, const std::string& jobname = "ham");

        StockholderAtomSfCalculator(const Crystal& crystal, const nlohmann::json& data);

        virtual ~StockholderAtomSfCalculator();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;
        virtual std::string name() const { return "gar"; }
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


        void printSubsystemsToXyz(const Crystal crystal,
            const std::vector<ham_settings::QmFragmentInCrystal>& subsystems);


        void printSubsystemsToMol2(const Crystal crystal,
            const std::vector<ham_settings::QmFragmentInCrystal>& subsystems);


        virtual void calculateFormFactors(const Vector3i& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual void calculateFormFactors(const std::vector<Vector3i>& hkl, std::vector< std::vector<std::complex<double> > >& formFactors, const std::vector<bool>& includeAtom) const;
    private:
        HirshfeldAtomModelSettings mHamSettings;

        std::vector<std::pair<std::string, std::string> > mModelInfo;



        //bool mPrerun = false;
        /**
        Intended to collects atom selection data from "remove" and "unit weight" entries
        in "cluster specific" sections of "multipoles" in aspher.json
        */
        //void collectAtomSelectionData(
        //    const std::string& fieldName, 
        //    const nlohmann::json& data, 
        //    std::vector<AtomSelecetionData>& selectionData);

        /**
        DistributedMultipoleSettings
        */
        //void collectDistributedMultipoleSettingsData(const nlohmann::json& data,
        //    DistributedMultipoleSettings& settings);


        /**
        converst string (in format like "H2;H6,1-X,Y,3/2-Z;C6,X,-Y,Z") to list of atoms
        */

        //void stringToAtomList(const std::string& s, std::vector<std::pair<std::string, std::string> >& atomList);

        std::shared_ptr<StockholderAtomFormFactorCalcManager> mFormFactorManager;

        std::shared_ptr<AtomicFormFactorCalculationsManager> mManager;
        std::shared_ptr<AnyIamCalculator> mIamCalculator;
        AnyScattererStructureFactorCalculator* mCalculator;
        //Crystal mCrystalHar, mCrystalIam;
        Crystal mCrystal;
        bool mSaveWfn;
        std::string mWfnFileName;
        std::string mJobName;


        void set(const Crystal& crystal, const HirshfeldAtomModelSettings& settings, bool electronScattering = false, const std::string& jobname = "ham");

    };
    /** @}*/
}
