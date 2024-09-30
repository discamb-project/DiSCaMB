#include "discamb_dev/MathUtilities/MixedRadixNumberIterator.h"

using namespace std;

namespace discamb {

    MixedRadixNumberIterator::MixedRadixNumberIterator()
    {

    }

    MixedRadixNumberIterator::MixedRadixNumberIterator(
        vector<size_t>& radices)
    {
        mRadices = radices;
        mValues.assign(radices.size(), 0);
    }


    void MixedRadixNumberIterator::get(
        vector<size_t>& value)
        const
    {
        value = mValues;
    }
 

    bool MixedRadixNumberIterator::increment()
    {
        bool incremented = false;

        return incremented;
    }


}