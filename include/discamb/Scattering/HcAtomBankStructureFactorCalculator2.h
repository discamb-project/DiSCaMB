#pragma once

#include "discamb/Scattering/SfCalculator.h"
#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"
#include "discamb/Scattering/AtomTypeHC_Parameters.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/HC_Model/SlaterOrbitalWfnData.h"
#include "discamb/Scattering/CombinedStructureFactorCalculator.h"
#include "discamb/Scattering/disordered_structure_fragments.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"
#include "discamb/AtomTyping/LocalCoordinateSystem.h"

#include <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */



    /*
    TAAM bank calclator with disorder handling
    */

    class HcAtomBankStructureFactorCalculator2 : public SfCalculator
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
        std::vector<disordered_structure_fragments::Fragment> mFragments;
        Crystal mExtendedCrystal;
        Crystal mCrystal;
        
        void setExtendedCrystal(const Crystal &crystal, const std::vector<disordered_structure_fragments::Fragment>& fragments);
        /**
        assigns atom types both for case with and without fragments definitions provided
        if taamFragments is not empty then
        output vector elements are in the same order as in mExtendedeCrystal
        i.e. loop over frgmants (loop over atoms in fragments (add if weight>0))
        else the output vector elements are in the same order as in crystal
        */
        void assignAtomTypes(
            const Crystal& crystal,
            const std::vector<disordered_structure_fragments::Fragment>& taamFragments,
            const std::vector<AtomType>& atomTypes,
            const DescriptorsSettings& settings,
            std::vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs,
            std::vector<int>& types,
            const std::string &assignemntInfoFile,
            const std::string& assignemntInfoCsvFile);

        std::vector<std::vector<int> > mNormal2extended;
        // mNormal2extendedByFrag[i][j] - idex in extended crystal the i-th atom of j-th fragment corresponds to
        // not defined if taamFragments is empty
        //std::vector<std::map<int, int> > mNormal2extendedByFrag;
        std::vector<int> mExtended2normal;
        std::vector<std::vector<int> > mFragmentToExtendedAtoms;
        AtomInCrystalID regularToExtendedCrystalAtomId(const AtomInCrystalID& atom) const;
        void updateExtendedCrystalAtoms(const std::vector<AtomInCrystal>& atoms);
        mutable std::vector<bool> mExtendedCrystalAtomsContribution;
        void updateExtendedCrystalAtomsContribution(const std::vector<bool>& countAtomContribution) const;
        mutable std::vector<TargetFunctionAtomicParamDerivatives> mExtended_dTarget_dparam;
        mutable SfDerivativesAtHkl mExtendedDerivatives1hkl;
        mutable std::vector < std::vector<std::complex<double> > > mExtendedFormFactors;
        mutable std::vector<std::complex<double> > mExtendedFormFactors1Hkl;
        void mergeExtendedDerivatives(
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<TargetFunctionAtomicParamDerivatives>& extended_dTarget_dparam) const;
        // [hkl index][atom index]
        void mergeExtendedFormFactors(
            const std::vector < std::vector<std::complex<double> > > &extended, 
            const std::vector<bool>& includeAtom,
            std::vector < std::vector<std::complex<double> > > &merged) const;
        // [atom index]
        void mergeExtendedFormFactors(
            const std::vector<std::complex<double> >& extended,
            const std::vector<bool>& includeAtom,
            std::vector < std::complex<double> >& merged) const;

        //order the same as in mNormal2extended
        // f[i] = sum( f_extended[mNormal2extended[i][j]] * mAtomContributionWeights[i][j]) for j=0..n
        std::vector<std::vector<double> > mAtomContributionWeights;
        //std::string mAssignmentInfo;
    public:
        /**
        if assignemntInfoFile is empty the assignment info is not printed
        */
        HcAtomBankStructureFactorCalculator2(
            const Crystal &crystal,
            const std::vector<AtomType> &atomTypes,
            const std::vector<AtomTypeHC_Parameters> &parameters,
            SlaterOrbitalWfnData::WfnDataBank slaterWavefunctionsDatabankId,
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
            const std::string &algorithm = "standard",
            const std::vector<disordered_structure_fragments::Fragment> &taamFragments = std::vector<disordered_structure_fragments::Fragment>(0));
        
        HcAtomBankStructureFactorCalculator2(const Crystal &crystal, const nlohmann::json &data);
        HcAtomBankStructureFactorCalculator2(const Crystal& crystal, const std::string &jsonFileName);

        HcAtomBankStructureFactorCalculator2(
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
            SlaterOrbitalWfnData::WfnDataBank slaterWavefunctionsDatabankId,
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
            const std::string& algorithm = "standard",
            const std::vector<disordered_structure_fragments::Fragment>& taamFragments = std::vector<disordered_structure_fragments::Fragment>(0));
        

        virtual ~HcAtomBankStructureFactorCalculator2();
        //static void readTaamFragments(const Crystal &crystal, const std::string& fileName, std::vector<disordered_structure_fragments::Fragment>& taamFragments);
        virtual void setAnomalous(const std::vector<std::complex<double> > & anomalous);
        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution);

        using SfCalculator::calculateStructureFactorsAndDerivatives;

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