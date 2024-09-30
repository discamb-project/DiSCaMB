#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/CrystalStructure/UnitCell.h"
#include <vector>

namespace discamb
{

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    class PeriodicGridInterpolator {
    public:
        PeriodicGridInterpolator();
        ~PeriodicGridInterpolator();
        //periodicValueIncluded if values from both sides of bounary are included - e.g. f(0,0,0) and f (1,0,0)
        void set(
            const std::vector<std::vector<std::vector<double> > > &values, 
            const Vector3d &a,
            const Vector3d& b,
            const Vector3d& c,
            bool periodicValueIncluded);

        void set(
            const std::vector<std::vector<std::vector<double> > >& values,
            const UnitCell &unitCell,
            bool periodicValueIncluded);
        double value(const Vector3d& r);
    private:
        std::vector<std::vector<std::vector<double> > > mV;
        //UnitCell mUnitCell;
        Matrix3d mCartesian2Fractional;
        double mDxFrac, mDyFrac, mDzFrac;
        int mNx, mNy, mNz;
        void set(const Matrix3d& cartesian2fractional, 
            const std::vector<std::vector<std::vector<double> > >& values,
            bool periodicValueIncluded);
    };
    /** @}*/
}

