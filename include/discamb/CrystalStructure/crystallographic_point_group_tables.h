#pragma once

#include "discamb/MathUtilities/Matrix3.h"

#include <vector>
#include <iosfwd>


namespace discamb {

/**
 * \addtogroup CrystalStructure
 * @{
 */


/** \brief Crystallographic point group operation tables. */

    namespace crystallographic_point_group_tables {
    
        /**
         * \brief Returns the point group operation table for the given point group.
         *
         * The point group is specified by its Hermann–Mauguin symbol.
         *
         * \param pointGroupSymbol The Hermann–Mauguin symbol of the point group.
         * \return A vector of matrices representing the point group.
         */
        std::vector<Matrix3i> getPointGroupOperations(const std::string& pointGroupSymbol);
        /**
         * \brief Returns the point group operation table for the given point group.
         *
         * The point group is specified by its Hermann–Mauguin symbol.
         *
         * \param pointGroupSymbol The Hermann–Mauguin symbol of the point group.
         * \return A vector of strings representing the point group.
         */

        std::vector <std::string> getPointGroupOperationsStr(const std::string& pointGroupSymbol);
        /**
         * \brief Returns the list of all available point groups.
         *
         * \return A vector of strings containing the names of all available point groups.
         */
        std::vector<std::string> getAvailablePointGroups();


        /**
         * \brief Finds the point group based on the symmetry operations.
         *
         * \param symmOps A vector of symmetry operations represented as matrices.
         * \return The Hermann–Mauguin symbol of the point group, or an empty string if not found.
         */
        std::string findPointGroup(const std::vector <std::string> &symmOps);
        std::string findPointGroup(const std::vector<Matrix3i>& symmOps);
        /**
        canonicalOrder[i] - the index symmOps[i] have on the discamb internal list for this point grouop
        */
        std::string findPointGroup(const std::vector<Matrix3i>& symmOps, std::vector<int> &canonicalOrder);
        std::string findPointGroup(const std::vector<std::string>& symmOps, std::vector<int>& canonicalOrder);
        
    }


/**@}*/

}// namespace discamb





