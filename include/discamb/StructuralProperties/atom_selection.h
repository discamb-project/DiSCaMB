#pragma once

#include <vector>
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/CrystalStructure/UnitCellContent.h"

#include "json.hpp"

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    struct AtomSubsetSelectionData {
        std::vector<Vector3i> periodicDirections;
        bool allSymmetryEquivalent = false;
        /**
        atomList[i].first  - atomic label
        atomList[i].second - symmetry operation (in X,Y,Z format)
        */
        std::vector<std::pair<std::string, std::string> > atomList;
        void set(nlohmann::json const& jsonData);
    };

    namespace atom_selection {

        void select_subset(
            const UnitCellContent &ucContent,
            const std::vector<std::pair<std::string, std::string> >& set,
            const AtomSubsetSelectionData& subsetSelectionData,
            std::vector<std::pair<std::string, std::string> >& subset);

        void select_subset(
            const UnitCellContent& ucContent,
            const std::vector<UnitCellContent::AtomID>& set,
            const AtomSubsetSelectionData& subsetSelectionData,
            std::vector<UnitCellContent::AtomID>& subset);

        void remove_subset(
            const UnitCellContent& ucContent,
            std::vector<std::pair<std::string, std::string> >& set,
            const AtomSubsetSelectionData& subsetSelectionData);

        void remove_subset(
            const UnitCellContent& ucContent,
            std::vector<UnitCellContent::AtomID>& set,
            const AtomSubsetSelectionData& subsetSelectionData);


        void merge(
            const UnitCellContent& ucContent,
            const std::vector<std::pair<std::string, std::string> >& set1,
            const std::vector<std::pair<std::string, std::string> >& set2,
            std::vector<std::pair<std::string, std::string> >& result);

        void merge(
            const std::vector<UnitCellContent::AtomID>& set1,
            const std::vector<UnitCellContent::AtomID>& set2,
            const std::vector<UnitCellContent::AtomID>& result);



    } //namespace atom_selection
    /** @}*/
} //namespace discamb

