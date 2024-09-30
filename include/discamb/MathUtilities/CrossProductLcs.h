#pragma once
#include "Vector3.h"

#include <string>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    class CrossProductLcs
    {
    public:
        CrossProductLcs(const std::string &dir1, const std::string &dir2, bool isR = true);
        ~CrossProductLcs();
        void calculate(const Vector3d &v1, const Vector3d &v2, Vector3d &x, Vector3d &y, Vector3d &z);
    private:
      
        // 1 or -1
        double mAxesOrderRelatedHandedenessCorrection;
        // 1 for R, -1 for L
        double mHandedeness;
        // 1 for X, 2 for Y, 3 for Z
        int mCoordinate1, mCoordinate2, mCoordinate3;
        //can be 1 or -1
        double mCoordinate1VectorMultiplier, mDirection12Multuiplier;

    };
    /** @}*/
}
