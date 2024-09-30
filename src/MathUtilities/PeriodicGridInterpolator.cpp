#include "discamb/MathUtilities/PeriodicGridInterpolator.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb
{
    
    namespace {
        double interpolateInCube(
            double x,
            double y,
            double z,
            double v000,
            double v001,
            double v010,
            double v011,
            double v100,
            double v101,
            double v110,
            double v111)
        {
            double v00, v01, v10, v11, v1, v0;
            v00 = x * v100 + (1 - x) * v000;
            v10 = x * v110 + (1 - x) * v010;
            v01 = x * v101 + (1 - x) * v001;
            v11 = x * v111 + (1 - x) * v011;
            v0 = y * v10 + (1 - y) * v00;
            v1 = y * v11 + (1 - y) * v01;
            return v1 * z + v0 * (1 - z);
        }

    }

    PeriodicGridInterpolator::PeriodicGridInterpolator()
    {
        mV.assign(2,vector<vector<double>>(2,vector<double>(2,0)));
        mDxFrac = 1.0; 
        mDyFrac = 1.0;
        mDzFrac = 1.0;
        mCartesian2Fractional.setToIdentity();
    }

    PeriodicGridInterpolator::~PeriodicGridInterpolator()
    {
    }

        
    void PeriodicGridInterpolator::set(
        const std::vector<std::vector<std::vector<double> > >& values,
        const Vector3d& aVec,
        const Vector3d& bVec,
        const Vector3d& cVec,
        bool periodicValueIncluded)
    {
        Matrix3d cartesian2fracrtional,
                 fractional2cartesian(
                   aVec[0], bVec[0], cVec[0],
                   aVec[1], bVec[1], cVec[1], 
                   aVec[2], bVec[2], cVec[2]);

        algebra3d::invert3d(fractional2cartesian, cartesian2fracrtional);

        set(cartesian2fracrtional, values, periodicValueIncluded);

    }

    void PeriodicGridInterpolator::set(
        const Matrix3d& cartesian2fractional,
        const std::vector<std::vector<std::vector<double> > >& values,
        bool periodicValueIncluded)
    {
        mV = values;

        mCartesian2Fractional = cartesian2fractional;

        if (values.empty())
            on_error::throwException("invalid dimansions of interpolation grid (should be at least 1x1x1", __FILE__, __LINE__);
        if (values[0].empty())
            on_error::throwException("invalid dimansions of interpolation grid (should be at least 1x1x1", __FILE__, __LINE__);
        if (values[0][0].empty())
            on_error::throwException("invalid dimansions of interpolation grid (should be at least 1x1x1", __FILE__, __LINE__);
        mNx = values.size();
        mNy = values[0].size();
        mNz = values[0][0].size();
        if (periodicValueIncluded)
        {
            mNx--;
            mNy--;
            mNz--;
        }

        mDxFrac = 1.0 / mNx;
        mDyFrac = 1.0 / mNy;
        mDzFrac = 1.0 / mNz;

    }


    void PeriodicGridInterpolator::set(
        const std::vector<std::vector<std::vector<double> > >& values,
        const UnitCell& unitCell,
        bool periodicValueIncluded)
    {
        set(unitCell.getCartesianToFractionalMatrix(), values, periodicValueIncluded);
    }

    double PeriodicGridInterpolator::value(
        const Vector3d& r)
    {
        Vector3d frac;

        double xCube, yCube, zCube;
        double x01, y01, z01;
        int x0, y0, z0;

        frac = mCartesian2Fractional* r;
        
        
        // move to 0-1
        x01 = frac[0] - floor(frac[0]);
        y01 = frac[1] - floor(frac[1]);
        z01 = frac[2] - floor(frac[2]);
        // origin of cube in grid
        x0 = int(x01 / mDxFrac);
        y0 = int(y01 / mDyFrac);
        z0 = int(z01 / mDzFrac);
        // fractional position in cube
        xCube = x01 / mDxFrac - x0;
        yCube = y01 / mDyFrac - y0;
        zCube = z01 / mDzFrac - z0;
        double result = interpolateInCube(
            xCube, yCube, zCube,
            mV [(x0    )%mNx] [(y0    )%mNy] [(z0    )%mNz],
            mV [(x0    )%mNx] [(y0    )%mNy] [(z0 + 1)%mNz],
            mV [(x0    )%mNx] [(y0 + 1)%mNy] [(z0    )%mNz],
            mV [(x0    )%mNx] [(y0 + 1)%mNy] [(z0 + 1)%mNz],
            mV [(x0 + 1)%mNx] [(y0    )%mNy] [(z0    )%mNz],
            mV [(x0 + 1)%mNx] [(y0    )%mNy] [(z0 + 1)%mNz],
            mV [(x0 + 1)%mNx] [(y0 + 1)%mNy] [(z0    )%mNz],
            mV [(x0 + 1)%mNx] [(y0 + 1)%mNy] [(z0 + 1)%mNz]);

        return result;
    }
    
    //    std::vector<std::vector<std::vector<double> > > mValues;
    //        UnitCell mUnitCell;
    
}