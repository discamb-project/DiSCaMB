#include "discamb/MathUtilities/MixedRadixNumberIterator.h"

using namespace std;

namespace discamb {

    MixedRadixNumberIterator::MixedRadixNumberIterator()
    {

    }

    MixedRadixNumberIterator::MixedRadixNumberIterator(
        const vector<int> & radices)
    {
        mRadices = radices;
        mValues.assign(radices.size(), 0);
    }


    void MixedRadixNumberIterator::get(
        vector<int>& value)
        const
    {
        value = mValues;
    }
 

    bool MixedRadixNumberIterator::increment()
    {
        bool incremented = false;
        int index=0, n = mRadices.size();
        while(!incremented && (index<n))
            if(mValues[index]+1<mRadices[index])
            { 
                incremented = true;
                mValues[index]++;
            }
            else
            {
                mValues[index] = 0;
                index++;
            }

        return incremented;
    }


}