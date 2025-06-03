#pragma once

#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/Scattering/AtomTypeHC_Parameters.h"
#include "discamb/AtomTyping/AtomType.h"

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    namespace taam_utilities {
       
       void type_assignment_to_unscaled_HC_parameters(
            const std::vector<AtomTypeHC_Parameters> &bankMultipoleParameters,
            const std::vector<int> &atomTypeAssignment,
            const std::vector<int> &atomicNumbers,
            HC_ModelParameters &parameters,
            bool notAssignedRepresentedWithSlaters,
            std::vector<int> &nonMultipolarAtoms);

        void type_assignment_to_HC_parameters(
            const std::vector<AtomTypeHC_Parameters>& bankMultipoleParameters,
            const std::vector<int>& atomTypeAssignment,
            const std::vector<double> & multiplicityTimesOccupancy, 
            const std::vector<int>& atomicNumbers,
            double totalCharge,
            HC_ModelParameters& parameters,
            bool notAssignedRepresentedWithSlaters,
            std::vector<int>& nonMultipolarAtoms);

        struct TaamAtomicChargeInfo {
            std::vector<double> atomicChargesBeforeScaling;
            std::vector<double> atomicChargesAfterScaling;
        };

        void type_assignment_to_HC_parameters(
            const std::vector<AtomTypeHC_Parameters>& bankMultipoleParameters,
            const std::vector<int>& atomTypeAssignment,
            const std::vector<double>& multiplicityTimesOccupancy,
            const std::vector<int>& atomicNumbers,
            double totalCharge,
            HC_ModelParameters& parameters,
            bool notAssignedRepresentedWithSlaters,
            std::vector<int>& nonMultipolarAtoms,
            TaamAtomicChargeInfo &atomicChargesInfo);


        void void_atom_type(
            HC_WfnType &wfnType, 
            HC_AtomTypeParameters &type_parameters);

        // & multiplicityTimesOccupancy = 
        // how many atom sites in unit cell times x occupancy
        void electroneutrality_Faerman_Price(
            std::vector<double>& typePval,
            const std::vector<double>& typePvalSigma,
            const std::vector<int>& atomToTypeMap,
            const std::vector<double> & multiplicityTimesOccupancy,
            const std::vector <int> atomicNumbers,
            double targetCharge);

        void electroneutrality_Faerman_Price(
            std::vector<double>& typePval,
            const std::vector<double>& typePvalSigma,
            const std::vector<int>& atomToTypeMap,
            const std::vector<double>& multiplicityTimesOccupancy,
            const std::vector <int> atomicNumbers,
            double targetCharge,
            TaamAtomicChargeInfo& atomicChargesInfo);


        

        int atomTypeRange(
            const AtomType &type,
            int maxPlanarRing,
            int &namedNeighboursRange,
            int &planarRingsRange,
            int &ring34range);

        int atomTypesRange(
            const std::vector<AtomType> &type,
            int maxPlanarRing,
            int &namedNeighboursRange,
            int &planarRingsRange,
            int &ring34range);

    }
    /** @}*/
}

