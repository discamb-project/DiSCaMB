#include "AtomType.h"

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
        void sortTypesByGenarality(std::vector<AtomType>& type);

        
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


