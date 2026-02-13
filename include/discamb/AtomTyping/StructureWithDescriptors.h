#pragma once

#include <map>

#include "AtomInStructureDescriptors.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicChemistry/MoleculeData.h"
#include "discamb/StructuralProperties/ad_hoc_connectivity_rules.h"
//#include "Molecule/Molecule.h"

#include "json.hpp"

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    struct DescriptorsSettings {
		/** threshold used to detect bond between atoms (A and B), they are bonded if the interatomic R distance fulfills the following condition:
		R <= covalent radious (atom A) + covalent radious (atom B) + threshold */
		double covalentBondThreshold = 0.4;
        double atomPlanarityThreshold = 0.1;
        double ringPlanarityThreshold = 0.1;
        double atomInRingPlanarityThreshold = 0.1;
		int atomInRingMaxNeighbourCount = 3;
        int maxPlanarRing = 8; // not set by user
        ad_hoc_connectivity_rules::Preset addHocConnectivityRulesPreset = 
            ad_hoc_connectivity_rules::Preset::None;

        std::map<std::pair<int, int>, double > maxBondLengthAromaticRing;/* = {
            {{6,6}, 1.45},
            {{6,7}, 1.41}
        };*/

        void readFromJson(const nlohmann::json& data);

        //double maxCcDistanceAromaticRing = 1.45;
        //double maxCnDistanceAromaticRing = 1.44;
    };


    struct StructureWithDescriptors
    {

        void set(const MoleculeData &moleculeData);

        void set(const std::vector<Vector3d>& positions, 
                 const std::vector< std::vector<int> >& connectivity,
                 const std::vector< std::vector<int> >& planes);


        void set(const std::vector<int> &atomicNumbers,
                 const std::vector<Vector3d> &positions,
                 const std::vector<std::string> &labels = std::vector<std::string>(0));

        void set(const std::vector<int> &atomicNumbers,
                 const std::vector<Vector3d> &positions, 
                 const std::vector<std::vector<int> > &predefinedConnectivity, 
                 const  std::vector<std::string> &labels = std::vector<std::string>(0));

        struct DataForCrystalFragment{
            std::vector<std::vector<int> > shells;
            std::vector<int> mapToCoreShellAtoms;
            int namedNeighboursRange = 0;
        };


        void setFromCrystalFragment(const std::vector<int> &atomicNumbers,
            const std::vector<Vector3d> &positions,
            const DataForCrystalFragment &data,
            const std::vector<std::string> &labels = std::vector<std::string>(0));

        void setFromCrystalFragment(const std::vector<int> &atomicNumbers,
            const std::vector<Vector3d> &positions,
            const std::vector<std::vector<int> > &predefinedConnectivity,
            const DataForCrystalFragment &data,
            const  std::vector<std::string> &labels = std::vector<std::string>(0));



		DescriptorsSettings settings;
        std::vector<AtomInStructureDescriptors> atomDescriptors;
        std::vector<std::vector<int> > connectivity;
        std::vector<std::vector<double> > bondLengths; // for choosing reference atom in LCS definition
        std::vector<std::vector<int> > planarRings;
        //std::vector<std::vector<int> > threeMemberRings;
        //std::vector<std::vector<int> > fourMemberRings;
        std::vector<double> planarRingsPlanarityEsd;
    };
    /**@}*/
}