#pragma once

#include "discamb/Scattering/SfCalculator.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"

#include <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class TaamSfCalculatorMultiOrderedImpl : public SfCalculator
    {
        std::vector<std::shared_ptr< HcAtomBankStructureFactorCalculator> > mCalculators;
        std::vector<std::vector<int> > mSubcrystal2CrystalAtomsMap;
        std::vector<std::vector<double> > mSubcrystalAtomsWeight;
        struct AtomReconstructionData {
            //[instance].first - subcrystal index
            //[iinsstance].second - atom in subcrystal
            std::vector<std::pair<int, int> > atoms;
            std::vector<double> weights;
        };
        
        std::vector<AtomReconstructionData> mReconstructionData;
        void setReconstructionData(
            const Crystal& crystal,
            const std::vector < std::vector <std::pair<std::string, double> > >& atomList);
        std::vector<Crystal> mSubcrystals;
        Crystal mCrystal;
        /* 
        the subcrystals which should be updated on update(atoms)
        used by methods that do no take current atomic parameters
        const std::vector<AtomInCrystal> &atoms
        as an argument
        */
        std::vector<Crystal> mSubcrystalsForUpdate;
        ReciprocalLatticeUnitCell mReciprocalLatticeUnitCell;
        void findSubstructureContributingAtoms(
            const std::vector<bool>& countAtomContribution,
            std::vector< std::vector<bool> >& substructureAtomContribution) const;
    public:
        /**
        if assignemntInfoFile is empty the assignment info is not printed
        */
        TaamSfCalculatorMultiOrderedImpl(
            const Crystal &crystal,
            const std::vector < std::vector <std::pair<std::string, double> > >& atomList,
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
        
        TaamSfCalculatorMultiOrderedImpl(const Crystal &crystal, const nlohmann::json &data);
        TaamSfCalculatorMultiOrderedImpl(const Crystal& crystal, const std::string &jsonFile);

        TaamSfCalculatorMultiOrderedImpl(
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
            const std::vector < std::vector <std::pair<std::string, double> > >& atomList,
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
        

        virtual ~TaamSfCalculatorMultiOrderedImpl();
        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous);

        using SfCalculator::calculateStructureFactorsAndDerivatives;

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
            const DerivativesSelector& selector);


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
        virtual std::string name() const { return std::string("TaamMultiOrder"); }

        //static void readSubstructuresFromFile(const std::string& fragmentsFile, std::vector<std::vector <std::pair<std::string, double> > >& fragmentAtoms);
    };
    /** @}*/
}