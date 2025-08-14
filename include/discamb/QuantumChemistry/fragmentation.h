#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/BasicChemistry/ChemicalElement.h"


#include "json.hpp"

#include <optional>

namespace discamb{

    /**
    * \addtogroup QuantumChemistry
    * @{
    */

    // fragment subset recipe
    struct FragmentPartConstructionData {
        typedef std::vector < std::pair<std::string, std::string> > AtomList;
        AtomList include;
        AtomList connect;
        AtomList remove;
        AtomList removeAndReplaceWithHydrogen;
        void set(const std::string& line);
        void set(const nlohmann::json &data);
        void set(const std::vector<std::string>& words);
        void clear();

        void getAtomList(
            const UnitCellContent &ucContent, 
            const std::vector<std::vector<UnitCellContent::AtomID> >& connectivity, 
            std::vector<UnitCellContent::AtomID> &atomList,
            std::vector<std::pair<UnitCellContent::AtomID, UnitCellContent::AtomID> >& cappingHydrogens) const;
    };

    struct FragmentConstructionData {
        std::vector<FragmentPartConstructionData> fragmentPartConstructionData;
        std::string label;
        int spin_multiplicity;
        int charge;
        bool set(const nlohmann::json& data);
    };

    struct CappingHydrogen {
        std::string bondedAtom;
        std::string bondedAtomSymmOp;
        std::string directingAtom;
        std::string directingAtomSymmOp;
    };

    //struct DistributedMultipoleCentersSettings {
    //    /**

    //    */
    //    double multipoleClusterThreshold = 8.0;
    //    std::vector<std::pair<std::string, std::string> > predefinedMultipoleClusterAtoms;
    //    std::vector<AtomSubsetSelectionData> atomsToOmit;
    //    std::vector<std::pair<AtomSubsetSelectionData, double> > atomsWithCustomWeights;
    //};

    struct FragmentAtoms {
        std::vector < std::pair<std::string, std::string> > atomList;
        std::vector <CappingHydrogen> cappingHydrogens;
    };

    struct QmFragmentInCrystal {
        std::string label = "system";
        FragmentAtoms atoms;
        // if defined then the definition of sites based on range/threshold is not applied
        std::vector < std::pair<std::string, std::string> > customMultipoleSites;
        int charge = 0;
        int spin_multiplicity = 1;
        void toXyzMol(const Crystal& crystal, std::vector<ChemicalElement>& elements, std::vector<Vector3d>& positions) const;
    };


namespace fragmentation{
    //void fragmentsData2fragments(const std::vector<FragmentData>& fragmentsData, );
    void make_qm_fragments(
        const Crystal& crystal,
        const std::vector<FragmentConstructionData>& fragmentConstructionData,
        std::vector<QmFragmentInCrystal>& qmFragments);

//    void find_fragments_atoms(
//        const Crystal& crystal,
//        std::vector<FragmentData>& fragmentData,
//        std::vector<FragmentAtoms>& fragmentsAtoms);


    Vector3d capping_h_position(const Crystal& crystal,
        const std::string& bondedAtom, const SpaceGroupOperation& bondedAtomSymmOp,
        const std::string& directingAtom, const SpaceGroupOperation& directingAtomSymmOp);

    Vector3d capping_h_position(const Crystal& crystal, const CappingHydrogen& cappingHydrogen);

    void from_input_file(
        const std::string& fileName,
        const Crystal& crystal,
        std::vector<int>& charge,
        std::vector<int>& spinMultiplicity,
        std::vector<std::string>& clusterLabels,
        std::vector< std::vector<std::pair<std::string, std::string> > >& clustersAtoms,
        std::vector<std::map<int, std::string> >& atomIdx2BasisSetMap);

    void with_hirshfrag(
        const Crystal& crystal,
        const std::string& hirshfragFile,
        const std::string& hirshfragFolder,
        std::vector<int>& charge,
        std::vector<int>& spinMultiplicity,
        std::vector<std::string>& clusterLabels,
        std::vector< std::vector<std::pair<std::string, std::string> > >& clustersAtoms);
    
    void intermolecular(
        const Crystal& crystal,
        std::vector<std::vector<std::pair<std::string, std::string> > >& clusterAtoms,
        std::vector<std::string>& clustersLabels,
        std::vector<int>& clustersCharges,
        std::vector<int>& clustersSpinMultiplicity,
        int spinMultiplicityHint = 0);



}
/**@}*/
}