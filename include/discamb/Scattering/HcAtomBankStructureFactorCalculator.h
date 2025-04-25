#include "discamb/Scattering/SfCalculator.h"
#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"
#include "discamb/Scattering/AtomTypeHC_Parameters.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/Scattering/CombinedStructureFactorCalculator.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"
#include "discamb/AtomTyping/LocalCoordinateSystem.h"

#include <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class HcAtomBankStructureFactorCalculator : public SfCalculator
    {
        // indices of atoms for which there were no data in Clementi-Roetti dataset
        std::vector<int> mIamAtoms;
        std::shared_ptr<SfCalculator> mHcCalculator, mIamCalculator;
        std::shared_ptr<CombinedStructureFactorCalculator> mCalculator;
        //
        std::vector<std::shared_ptr<CombinedStructureFactorCalculator> > mCalculators;
        std::vector < std::shared_ptr<SfCalculator> > mHcCalculators, mIamCalculators;
        int mN_Threads = 1;
        double mTotalCharge = 0.0;
        std::vector<std::pair<std::string, std::string> > mModelInfo;
        std::string mAlgorithm = "standard";
        //std::string mAssignmentInfo;
    public:
        /**
        if assignemntInfoFile is empty the assignment info is not printed
        */
        HcAtomBankStructureFactorCalculator(
            const Crystal &crystal,
            const std::vector<AtomType> &atomTypes,
            const std::vector<AtomTypeHC_Parameters> &parameters,
            bool electronScattering = false,
            const DescriptorsSettings &settings= DescriptorsSettings(),
            const std::string &assignemntInfoFile = std::string(),
            const std::string &assignmentCsvFile = std::string(),
            const std::string& parametersInfoFile = std::string(),
            const std::string& multipolarCif = std::string(),
            int nThreads = 1,
            double unitCellCharge = 0,
            bool scaleToMatchCharge = true,
            const std::string& iamTable = std::string(),
            bool iamElectronScattering = false,
            bool frozen_lcs = false,
            const std::string &algorithm = "standard");
        
        HcAtomBankStructureFactorCalculator(const Crystal &crystal, const nlohmann::json &data);

        HcAtomBankStructureFactorCalculator(
            const Crystal& crystal,
            const nlohmann::json& data,
            const std::string & bankString/*,
            std::string &assignemntInfo,
            bool generateAssignementInfo*/);

        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;

        void set(const Crystal& crystal,
            const nlohmann::json& data,
            const std::string& bankString/*,
            bool generateAssignementInfo=false*/);


        void set(
            const Crystal &crystal,
            const std::vector<AtomType> &atomTypes,
            const std::vector<AtomTypeHC_Parameters> &parameters,
            bool electronScattering = false,
            const DescriptorsSettings &settings = DescriptorsSettings(),
            const std::string &assignemntInfoFile = std::string(),
            const std::string &assignmentCsvFile = std::string(),
            const std::string& parametersInfoFile = std::string(),
            const std::string& multipolarCif = std::string(),
            int nThreads = 1,
            double unitCellCharge = 0,
            bool scaleToMatchCharge = true,
            const std::string& iamTable = std::string(),
            bool iamElectronScattering = false,
            bool frozen_lcs = false,
            const std::string& algorithm = "standard"/*,
            bool generateAssignmentInfo = false*/);
        

        virtual ~HcAtomBankStructureFactorCalculator();
        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous);
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

		virtual void calculateFormFactors(const Vector3i& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual void calculateFormFactors(const std::vector<Vector3i>& hkl, std::vector< std::vector<std::complex<double> > >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual void calculateFormFactorsCart(const Vector3d& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;
        //virtual void calculateFormFactors(const std::vector<Vector3d>& hkl, std::vector< std::vector<std::complex<double> > >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual std::string name() const { return std::string("MATTS"); }

    private:
        HC_ModelParameters mHcParameters;

        void printAssignedMultipolarParameters(
            const std::string& fName, 
            const Crystal& crystal,
            const std::vector < LocalCoordinateSystem<AtomInCrystalID> > &lcs,
            const std::vector<AtomType>& atomTypes, 
            const std::vector<AtomTypeHC_Parameters>& bankParameters, 
            const HC_ModelParameters &multipoleModelPalameters,
            const std::vector<int> &assignedTypeIdx);

        //void printHcParameters(std::ofstream &out, const AtomTypeHC_Parameters& bankParameters);
        void printHcParameters(std::ofstream& out, const HC_ModelParameters& multipoleModelPalameters, int atomIdx);
    };
    /** @}*/
}