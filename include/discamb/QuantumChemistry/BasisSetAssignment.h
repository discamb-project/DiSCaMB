#pragma once

#include "discamb_dev/QuantumChemistry/BasisSetId.h"

#include <string>
#include <map>
#include <utility>

namespace discamb {

    class BasisSetAssignment {
        public:
            BasisSetAssignment();
            BasisSetAssignment(const std::string &basisSetName);
            BasisSetAssignment(const std::string& basisSetId);
            ~BasisSetAssignment();
            
            void set(const std::string& name);
            void set(const BasisSetId& basisSetId);
            void get(const std::string& basisSetId) const;
            
            void setByElement(const std::map<std::string, std::string>  &symbolAndBasis);
            void setByElement(const std::map<std::string, BasisSetId> &symbolAndId);
            /*assigns values only for elements for which element specific basis set was set. */
            void getByElement(std::map<std::string, BasisSetId>& symbolAndId) const;

            void setByAtomIdx(const std::map < size_t, std::string >& atomIdxAndBasis);
            void setByAtomIdx(const std::map < size_t, BasisSetId >& atomIdxAndBasisId);
            /*assigns values only for atoms for which atom specific basis set was set. */
            void getByAtomIdx(std::map < size_t, BasisSetId >& atomIdxAndBasisId) const;
        private:
            BasisSetId mBasis;
            /** Defined only for those elements, for which it differs from general specs in mBasis. */
            std::map<std::string, BasisSetId> mBasisByElement;
            /** Defined only for those elements, for which it differs from general specs in mBasis
                or if entry in mBasisByElement for the corresponding element exists, from spces in mBasisByElement. */
            std::map < size_t, BasisSetId > mBasisByAtomIdx;
    };

}