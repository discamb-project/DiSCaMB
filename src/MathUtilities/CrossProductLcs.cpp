#include "discamb/MathUtilities/CrossProductLcs.h"

#include <map>
#include <cmath>

using namespace std;

namespace discamb {

    CrossProductLcs::CrossProductLcs(
        const string &dir1, 
        const string &dir2,
        bool isR)
    {
        map<string, int> directionIndexMap{ {"X",1 }, {"Y",2}, {"Z",3}, {"-X",-1}, {"-Y",-2}, {"-Z",-3} };

        
        mCoordinate1VectorMultiplier = directionIndexMap[dir1] / abs(directionIndexMap[dir1]);
        mCoordinate1 = abs(directionIndexMap[dir1]);
        mDirection12Multuiplier = directionIndexMap[dir2] / abs(directionIndexMap[dir2]);
        mCoordinate2 = abs(directionIndexMap[dir2]);
        mCoordinate3 = 6 - mCoordinate1 - mCoordinate2;

        // 1 for R, -1 for L
        isR ? mHandedeness = 1 : mHandedeness = -1;

        mAxesOrderRelatedHandedenessCorrection = 1 - 2 * ((2 + mCoordinate2 - mCoordinate1) % 3);

    }

    CrossProductLcs::~CrossProductLcs()
    {

    }

    void CrossProductLcs::calculate(
        const Vector3d &v1, 
        const Vector3d &v2, 
        Vector3d &x, 
        Vector3d &y, 
        Vector3d &z)
    {
        Vector3d point_1_position, point_2_position, point_3_position, reference_atom_position;
        Vector3d direction[3], direction_12;
        Vector3d *axes[] = { &x,&y,&z };

        direction[0] = mCoordinate1VectorMultiplier * v1;
        direction[0] /= sqrt(direction[0] * direction[0]);

        direction_12 = mDirection12Multuiplier * v2;

        direction[2] = mAxesOrderRelatedHandedenessCorrection * mHandedeness * cross_product(direction[0], direction_12);
        direction[2] /= sqrt(direction[2] * direction[2]);

        direction[1] = mAxesOrderRelatedHandedenessCorrection * mHandedeness * cross_product(direction[2], direction[0]);

        *axes[mCoordinate1 - 1] = direction[0];
        *axes[mCoordinate2 - 1] = direction[1];
        *axes[mCoordinate3 - 1] = direction[2];
    }



}