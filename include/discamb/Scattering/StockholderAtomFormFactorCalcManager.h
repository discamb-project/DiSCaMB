#pragma once

#include "AtomicFormFactorCalculationsManager.h"
#include "AtomRepresentativeInfo.h"
#include "discamb/Scattering/NGaussianFormFactor.h"
#include "AtomicDensityTransform.h"
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"
#include "discamb/QuantumChemistry/distributed_multipoles.h"
#include "discamb/QuantumChemistry/ProatomDB.h"
#include "discamb/QuantumChemistry/ElectronDensityPartitionType.h"
#include "discamb/Scattering/NumericalSphericalAtomFormFactor.h"
#include "discamb/Scattering/HirshfeldAtomModelSettings.h"
#include "discamb/QuantumChemistry/HirshfeldPartition.h"
#include "discamb/QuantumChemistry/distributed_multipoles.h"


#include "json.hpp"

#include <memory>
#include <set>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    /*
    
    File naming scheme:

    wfnFileName = 
        mJobName + "_" + mClusterLabels[i] + "_wfn_" + numberTo000string(indexAddedToFileNames) + "." + format;

    densityFileName = mJobName + string("_") + mClusterLabels[clusterIdx] +
                string("_ed_") + numberTo000string(indexAddedToFileNames) + ".h5";
    */
    class StockholderAtomFormFactorCalcManager : public AtomicFormFactorCalculationsManager {
    public:
        StockholderAtomFormFactorCalcManager(const Crystal &crystal, const HirshfeldAtomModelSettings &settings);
        virtual ~StockholderAtomFormFactorCalcManager();
        virtual void getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo) const {};
        // not implemented
        void set(const Crystal& crystal, const HirshfeldAtomModelSettings& settings);
        // not implemented
        void setModel(const HirshfeldAtomModelSettings& settings);
        
        virtual void update(const std::vector<AtomInCrystal>& atoms);

        virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i& hkl) const;
        virtual std::complex<double> calculateCart(int atomIdx, const Vector3d& hkl) const;

        virtual void calculateFrac(
            const std::vector<Vector3d>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;

        virtual void calculateFrac(
            const std::vector<Vector3i>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;


        virtual void calculateCart(
            const std::vector <Vector3d>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;

        enum class EE_SFC_CalcOption {
			from_scratch, 
			from_converged_wfn, 
            from_converged_multipoles,
			from_non_converged_wfn, 
			from_not_converged_multipoles};

        void calculateAtomicDensities();

    private:
        HirshfeldAtomModelSettings mSettings;
        
        void setOrbital(int idx);     
        
        void printDiagnosticInfo();
        /*
        clusterAtoms[i] ith atom of the cluster used in calculations, .first - atom label, .second - symmetry card
        representatives[i] - list of atoms in cluster to represent i-th atom in asymmetric unit
        */
        void setSubsystems(const std::vector<ham_settings::QmFragmentInCrystal>& subsystems,
            const std::vector<std::vector<AtomRepresentativeInfo> >& representatives);

        void setSubsystemWfnCalcData(
            int subsystem, 
            WaveFunctionCalculationData& data, 
            const std::vector<double> pointCharges,
            const std::vector<Vector3d> pointChargePosition,
            const std::string& outputWfnFileName);       



        void setSphericalAtomicDensities(const ProatomDB& db);
        
        int mOrbitalIdx = -1;
        //int mN_Threads;
        std::vector<int> mAtomicNumbers;

        void setMultipolesfromTaam();

        std::string mJobName;
		//##### Quantum Chemistry ######

		std::shared_ptr<WaveFunctionDataGeneratorRunner> mWfnGeneratorRunner;
        std::vector< WaveFunctionCalculationData> mWfnGeneratorData;

        //##### Spherical part #########################

        std::vector<NumericalSphericalAtomFormFactor> mSphericalFormFactors;
        std::vector<SphericalAtomicDensity> mSphericalDensities;
        std::vector<int> mAtomToSphericalTypeMap;
        std::vector<int> mSphericalTypeAtomicNumber;

        bool mUseSpehricalDensitySubstraction;               

        // [cluster idx][atom idx][value]
        std::vector < std::vector<std::vector<double> > > mAtomInClusterElectronDensity;

        // [z][point idx]
        std::vector < std::vector<Vector3d> > mElementGridPoints;
        // [z][point idx]
        std::vector < std::vector<double> > mElementIntegrationWeights;
        // to be executed on useHorton(), setCluster() && setGrid()
        void setElementalIntegrationGrids();

        // allocate mIntegrationGridPoints, mAtomicDensities, mIntegrationWeights
        // with no allocation of the innermost container (os sie=number of grid points)
        // e.g. mAtomicDensities.resize(nAtoms,vector<..>(mRepresentatives[atom].size())
        void allocateElectronDensityContainers();

		//#### MULTIPOLES ########
        //bool mStartWithTaamMultipoles = true;
        bool mTryToReadMultipolesFromFile = true;       

        void setDistributedMultipoleCenters();
		std::vector< std::vector<std::pair<int, SpaceGroupOperation> > > mMultipoleBearingAtoms;
        std::vector< std::vector< std::optional<double> > > mCustomWeightMultipoles;
        
        void calculateMultipoles();
        void calculateMultipoles(int atomIdx, int z, ElectricMultipoles&multipoles);

        WaveFunctionFileFormat getFormat() const;
            

		std::vector<ElectricMultipoles> mMultipoles, mMultipolesPrevious;
		bool multipoleCalculationConverged() const;

		static void setMultipolesTo0(std::vector<ElectricMultipoles>& multipoles);
            

		//-------------

		int mElectrostaticEmbeddingSCF_CycleIndex;

		static void saveMultipolesToFile(const std::string fileName, const std::vector<ElectricMultipoles>& multipoles);
		static void readMultipolesFromFile(const std::string fileName, std::vector<ElectricMultipoles>& multipoles);

		// ---------------------
        void calculateAsymmetricUnitPointChargeMultipoleRepresentation(
            std::vector< std::vector<double> >& chargeAsymm,
            std::vector< std::vector<Vector3d> >& positionAsymm);


        
        void calculatePointChargeMultipoleRepresentation(
            int clusterIdx,
            const std::vector< std::vector<double> >& chargeAsymm,
            const std::vector< std::vector<Vector3d> >& positionAsymm,
            const std::vector<Vector3d>& asymmUnitAtomCartesian,
            std::vector<double>& charge,
            std::vector<Vector3d>& position);


		void calculatePointChargeMultipoleRepresentation(std::vector<double>& charge, std::vector<Vector3d>& position);

        void calculatePointChargeMultipoleRepresentation(
            std::vector< std::vector<double> >& charge, 
            std::vector< std::vector<Vector3d> >& position);

        std::map<std::string, int> mAtomLabel2IdxInAsymmUnit;


        //#### CLUSTERS ########

        std::vector< std::set<int> > mSubsystemAtomsWhichAreRepresentatives;

        
        std::vector<std::vector<std::shared_ptr<AtomicDensityTransform> > > mTransforms;
        void setTransforms();

        void setAndPrintRepresentatives(const std::vector < std::vector<AtomRepresentativeInfo> > &representatives);

        std::vector<std::vector<int> > mSubsystemAtomicNumbers;
        std::vector<std::vector<Vector3d> > mSubsystemAtomicPositions; 

        // [subsystem idx][print out order idx] .first and .second
        // corresponds to entry in mRepresentatives[.first][.second]
        std::vector<std::vector<std::pair<int, int> > > mListOfRepresentativesInSubsystem;
        std::vector < std::vector<int> > mAtomsInSubsystemWhichAreRepresentatives;
        //######################

        ReciprocalLatticeUnitCell mReciprocalSpaceUnitCell;


        void electronDensityFromWfn(int indexAddedToFileNames);
        void electronDensityFromWfn(int indexAddedToFileNames, int clusterIdx);
        void electronDensityFromWfn(int indexAddedToFileNames, int clusterIdx, int molecularOrbitalIdx);

        void setPositions(const Crystal &crystal);
        Crystal mCrystal;
        void calculateWfns(int indexAddedToFileNames);

        // previousStep is 0 if there were no previous step
		void calculateElectronDensityInSelfConsistentMultipoleEmbedding(
            const std::vector<ElectricMultipoles>& multipoles = std::vector<ElectricMultipoles>(0),
            bool multipolesConverged = false,
            int previousStep = 0);

        
        bool tryToSetMultipolesFromFile(int& previousStep, std::string& fName);
        bool tryToFindAndReadMultipolesFile(
            int& previousStep, 
            std::string& fName,
            std::vector< ElectricMultipoles> &multipoles) const;

        void calculateElectronDensityWithPredefinedMultipoles(
            bool useTaamMultipolesIfNoOtherProvided = true,
            bool tryToReadMultipolesFromFile = true,
            const std::vector<ElectricMultipoles>& multipoles = std::vector<ElectricMultipoles>(0));

        void calculateFracWithSphericalHarmonics(
            const std::vector<Vector3d>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;



		// converts n to 3 character string with leading zeros 
		// e.g. 12 woulb be converted to 012, 8 to 008, 123 to 123 and 1234 to 1234
		std::string numberTo000string(int n);

        void setElectronDensityPartition(
            int clusterIdx, 
            const std::string& wfnFile, 
            const std::string& wfnFormat, 
            std::shared_ptr<ElectronDensityPartition>& partition);
    };
    /** @}*/
}



