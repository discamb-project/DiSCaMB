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

    struct TaamSfCalculatorSettings {
        void set(const Crystal &crystal, const nlohmann::json& data, const std::string& bankText=std::string());
        void set(const Crystal& crystal, const std::string &file, const std::string &bankText = std::string());

        std::vector<AtomType> atomTypes;
        std::vector<AtomTypeHC_Parameters> parameters;
        SlaterOrbitalWfnData::WfnDataBank slaterWavefunctionsDatabankId = SlaterOrbitalWfnData::WfnDataBank::CR;
        bool electronScattering = false;
        DescriptorsSettings descriptorsSettings;
        std::string assignmentInfoFile;
        std::string assignmentCsvFile;
        std::string parametersInfoFile;
        std::string multipolarCif;
        int nThreads = 1;
        double unitCellCharge = 0.0;
        bool scaleToMatchCharge = true;
        std::string iamTable;
        bool iamElectronScattering = false;
        bool frozen_lcs = false;
        bool def_val_symm = false;
        bool splitWithLabels = false;
        //works with labels like (altloc=B): "H2.B   1    A    H105    Z N      1    A X CA     1    A"
        bool splitWithInternalAltlocLabels = false;
        std::string algorithm = "standard";
        std::vector < std::vector <std::pair<std::string, double> > > orderedSubcrystalAtoms;
        std::vector<disordered_structure_fragments::Fragment> taamFragments;
    };



    /*
    TAAM bank calculator 
    */

    class TaamSfCalculator : public SfCalculator
    {
        std::shared_ptr<SfCalculator> mImplementation;

    public:
        /**
        if assignemntInfoFile is empty the assignment info is not printed
        */
        TaamSfCalculator(
            const Crystal& crystal,
            const TaamSfCalculatorSettings &settings);

        TaamSfCalculator(const Crystal& crystal, const nlohmann::json& data);
        TaamSfCalculator(const Crystal& crystal, const std::string& jsonFileName);


        virtual ~TaamSfCalculator();

        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const;

        void set(
            const Crystal& crystal,
            const TaamSfCalculatorSettings& settings);

        //static void readTaamFragments(const Crystal &crystal, const std::string& fileName, std::vector<disordered_structure_fragments::Fragment>& taamFragments);
        virtual void setAnomalous(const std::vector<std::complex<double> >& anomalous);
        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& dTarget_df,
            const std::vector<bool>& countAtomContribution);

        virtual void calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal>& atoms,
            const std::vector<Vector3i>& hkl,
            std::vector<std::complex<double> >& f,
            std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
            const std::vector<std::complex<double> >& dTarget_df,
            const std::vector<bool>& countAtomContribution,
            const DerivativesSelector& selector);


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
        virtual void calculateFormFactorsCart(const Vector3d& hkl, std::vector<std::complex<double> >& formFactors, const std::vector<bool>& includeAtom) const;
        //virtual void calculateFormFactors(const std::vector<Vector3d>& hkl, std::vector< std::vector<std::complex<double> > >& formFactors, const std::vector<bool>& includeAtom) const;
        virtual std::string name() const { return std::string("TAAM"); }
    };
}

