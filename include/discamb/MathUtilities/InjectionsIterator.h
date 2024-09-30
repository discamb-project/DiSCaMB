#pragma once

#include <cstddef>
#include <vector>


namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    /**
    Represents all injections of n-element set onto m-element set.
    */

    class InjectionsIterator
    {


    public:
        InjectionsIterator();
        InjectionsIterator(int larger, int smaller);
        ~InjectionsIterator();
        void set(int larger, int smaller);
        int operator()(int i) const;
        int largerSetSize() const;
        int smallerSetSize() const;
        bool next();
        void setToStart();
    private:
        int mLarge, mSmall;
        std::vector<bool> mSubset;
        std::vector<int> mPermutation;
        void setFirstPermutation();
    };
    /** @}*/
}

