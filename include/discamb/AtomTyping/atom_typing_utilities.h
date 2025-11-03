#include "AtomType.h"

#include <utility>

namespace discamb{

/**
* \addtogroup AtomTyping
* @{
*/

namespace atom_typing_utilities{

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
        
        // less general first
        // the less levels above the more general type
        void sortTypesByGenarality_LevelsAbove(std::vector<AtomType>& type);

        // less general first
        // the more levels below the more general type
        void sortTypesByGenarality_LevelsBelow(std::vector<AtomType>& type);


        void typeGeneralization(
            const std::vector<AtomType>& type,
            std::vector<std::vector<int> >& typesGeneralized,
            std::vector<std::vector<int> >& hierarchyLevel);

        /*
        generalizedByMultipleTypesAtSameLevel[i].first - generalized type index
        generalizedByMultipleTypesAtSameLevel[i].second - list of types that generalize the type at the same hierarchy
                                                          level (one above the type level)
        */
        void typeGeneralizationDiagnostics(
            const std::vector<std::vector<int> >& typesGeneralized,
            const std::vector<std::vector<int> >& hierarchyLevel,
            std::vector<std::pair<int, int> > &equivalentTypes,
            std::vector<std::pair<int, std::vector<int> > > &generalizedByMultipleTypesAtSameLevel);
        
        void findGeneralizingTypes(
            const std::vector<std::vector<int> >& typesGeneralized,
            std::vector<std::vector<int> >& generalizedBy); 
        
        inline std::string planarity_as_string(const Tribool& planarity)
        {
            if (planarity == Tribool::False)
                return "NOT_PLANAR";
            else if (planarity == Tribool::True)
                return "PLANAR";

            return "ANY_PLANARITY";
        }

}
/**@}*/
}


