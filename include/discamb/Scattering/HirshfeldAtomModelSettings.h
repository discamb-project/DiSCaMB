#pragma once


#include "discamb/BasicChemistry/ChemicalElement.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/QuantumChemistry/ElectronDensityPartitionType.h"
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"
#include "discamb/QuantumChemistry/distributed_multipoles.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/StructuralProperties/atom_selection.h"
#include "AtomRepresentativeInfo.h"

#include "json.hpp"

#include <map>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    //Hirshfeld Atom model settings

    namespace ham_settings {

        struct PartitionData {
            void set(const nlohmann::json& data);
            ElectronDensityPartitionType type = ElectronDensityPartitionType::Hirshfeld;
            nlohmann::json partitionSpecificData;
            double power = 1.0;
            std::optional<double> edCalcAtomIncludeRange;
        };

        struct WfnCalculation {
            
            void set(const nlohmann::json& data);
            WaveFunctionDataGeneratorRunner::Program qmProgram;
            nlohmann::json qmProgramSpecificData;
            QmSettings qmSettings;
            
            enum class Restart {
                from_scratch, from_converged_wfn, from_converged_density,
                from_converged_multipoles, from_non_converged_wfn, from_not_converged_density,
                from_not_converged_multipoles
            };
            
            Restart restart;
            int restartStep;
            bool useDistributedMultipoles = true;
            double multipoleClusterThreshold = 8.0;
            std::optional<HardwareResources> hardwareSettings;
        };

        struct CrystalFragmentWfnCalculation {
            void set(const nlohmann::json& data, double clusterThreshold);
            std::optional<WfnCalculation> wfnCalculation;
            std::map<int, std::string> atomicIdx2BasisSetMap;
            std::optional<DistributedMultipoleCentersSettings> distributedMultipoleCluster;
        };

        struct QmFragmentInCrystal {
            std::string label = "system";
            FragmentAtoms atoms;
            int charge = 0;
            int spin_multiplicity = 1;
            void toXyzMol(const Crystal& crystal, std::vector<ChemicalElement>& elements, std::vector<Vector3d>& positions) const;
        };


        struct MultipoleExpansionCalculation {
            void set(const nlohmann::json& data);
            bool calculate = false;
            int multipoleExpansionLevel = 1;
            double multipoleChargeDistance = 0.01;
            double multipoleConvergenceThreshold = 0.003;
            double clusterThreshold = 8.0;
            int maxN_Steps = 1;
            bool startWithTaam = true;
            bool lookForMultipoleFile = true;
            std::string previousMultipolesFile;
        };

        struct Diagnostics {
            void set(const nlohmann::json& data);
            bool prerun;
            bool printSubsystemsToXyz;
            bool printMultipoleClustersToXyz;
            bool printSubsystemsToMol2;
            bool printMultipoleClusterToMol2;
        };


        struct IntegrationGrid {
            int angularGridSize = 590;
            int radialGridSize = 75;
        };

        struct FormFactorsCalculation
        {
            void set(const nlohmann::json& data);
            IntegrationGrid integrationGrid;
            std::map<int, IntegrationGrid> elementSpecificIntegrationGrid;
            
        };
        /**
        subsystems are set with default method if no specified in data
        */
        void setCrystalFragments(
            const nlohmann::json& data, 
            const Crystal &crystal,
            std::vector<QmFragmentInCrystal>& crystalFragments);

        void setRepresentatives(const nlohmann::json& data,
            const Crystal& crystal,
            const std::vector<QmFragmentInCrystal>& subsystems,
            std::vector<std::vector<AtomRepresentativeInfo> > &representatives);

        void setSubsystemsWfnCalculation(
            const nlohmann::json& data,
            const std::vector<QmFragmentInCrystal>& subsystems,
            double clusterThreshold,
            std::vector< CrystalFragmentWfnCalculation > & subsystemWfnCalculation);

        //void subsystemWfnCalcData(
        //    const Crystal& crystal, 
        //    const Subsystem& subsystem, 
        //    const WfnCalculation &wfnCalc,
        //    const SubsystemWfnCalculation &subsystemWfnCalc,
        //    WaveFunctionCalculationData &wfnData);

        void wfnCalcData(
            const WfnCalculation& wfnCalc,
            WaveFunctionCalculationData& wfnData);

    }

    struct HirshfeldAtomModelSettings
    {
        void set(const nlohmann::json &data, const Crystal& crystal);
        std::vector<ham_settings::QmFragmentInCrystal> crystalFragments;
        ham_settings::WfnCalculation wfnCalculation;
        std::vector<ham_settings::CrystalFragmentWfnCalculation> fragmentWfnCalculation;
        ham_settings::PartitionData electronDensityPartition;
        ham_settings::MultipoleExpansionCalculation multipoleExpansion;
        ham_settings::FormFactorsCalculation formFactorsCalculation;
        HardwareResources hardware{ 1,1000 }; // 1 core 1000MB
        ham_settings::Diagnostics diagnostics;
        std::vector<std::vector<AtomRepresentativeInfo> > representatives;
        int sphericalHarmonicsExpansionLevel = -1; // if less than 0 then no expansion
    };
    /** @}*/
}