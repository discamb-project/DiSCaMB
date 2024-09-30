#pragma once

#include <vector>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */

    class MixedRadixNumberIterator {
        std::vector<int> mRadices;
        std::vector<int> mValues;

    public:
        MixedRadixNumberIterator();
        MixedRadixNumberIterator(const std::vector<int> &radices);
        void get(std::vector<int>& value) const;
        bool increment();
    };
    /** @}*/
}

