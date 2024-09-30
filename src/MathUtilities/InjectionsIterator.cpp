#include "discamb/MathUtilities/InjectionsIterator.h"

#include <algorithm>

using namespace std;

namespace discamb {

    void InjectionsIterator::setFirstPermutation()
    {
        int i, j = 0;
        for (i = 0; i < mLarge; i++)
            if (mSubset[i])
                mPermutation[j++] = i;
    }


    int InjectionsIterator::largerSetSize()
        const
    {
        return mLarge;
    }


    int InjectionsIterator::smallerSetSize()
        const
    {
        return mSmall;
    }


    InjectionsIterator::InjectionsIterator()
    {
        mLarge = 1;
        mSmall = 1;
        mPermutation.push_back(0);
    }

    InjectionsIterator::InjectionsIterator(
        int larger,
        int smaller)
    {
        set(larger, smaller);
    }

    InjectionsIterator::~InjectionsIterator()
    {
    }

    void InjectionsIterator::set(
        int larger,
        int smaller)
    {
        mLarge = larger;
        mSmall = smaller;
        mSubset.resize(mLarge);
        mPermutation.resize(mSmall);
        setToStart();
    }

    int InjectionsIterator::operator()(
        int i)
        const
    {
        return mPermutation[i];
    }

    bool InjectionsIterator::next()
    {
        if (!std::next_permutation(mPermutation.begin(), mPermutation.end()))
        {
            if (std::next_permutation(mSubset.begin(), mSubset.end()))
                setFirstPermutation();
            else
                return false;
        }
        return true;
    }

    void InjectionsIterator::setToStart()
    {
        std::fill(mSubset.end() - mSmall, mSubset.end(), true);
        std::fill(mSubset.begin(), mSubset.begin() + mLarge - mSmall, false);
        setFirstPermutation();
    }

}


