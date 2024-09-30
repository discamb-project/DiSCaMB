#include "discamb/Scattering/NumericalSphericalAtomFormFactor.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include "discamb/BasicUtilities/Constants.h"

#include <cmath>

using namespace std;

namespace discamb {

    NumericalSphericalAtomFormFactor::NumericalSphericalAtomFormFactor(
        const SphericalAtomicDensity& sphericalDensity)
    {
        set(sphericalDensity);
    }

    NumericalSphericalAtomFormFactor::~NumericalSphericalAtomFormFactor() 
    {
    }

    void NumericalSphericalAtomFormFactor::set(
        const SphericalAtomicDensity& sphericalDensity)
    {
        mStart = 0.0;
        mStep = 0.001;
        int i, n = 6001;
        mValues.resize(n);
        double h, f, a = constants::Angstrom;

        vector<double> densities(5000); // in e/A^3
        double two_pi = 2.0 * M_PI;

        for (int i = 0; i < 5000; i++)
            densities[i] = sphericalDensity.calculate(a * i * 0.001) * a * a * a;

        for (i = 0; i < n; i++)
        {
            h = mStart + i * mStep;
            f = 0;
            double r, ff = 0;

            if (h < 1e-10)
                for (int j = 1; j < 5000; j++)
                {
                    r = j * 0.001;
                    ff += r * r * densities[j];
                }
            else
                for (int j = 1; j < 5000; j++)
                {
                    r = j * 0.001;
                    ff += r * r * densities[j] * sin(two_pi * h * r) / (two_pi * h * r);
                }
            ff *= 2 * two_pi * 0.001;
                
            mValues[i] = ff;
        }
    }

    double NumericalSphericalAtomFormFactor::calculate(
        double h)
        const 
    {
        int n = static_cast<int>((h-mStart) / mStep);
        if (n < 0)
            return mValues[0];
        if (n >= 6000)
            return 0.0;
        double d_left, d_right, w_left, w_right;
        d_left = h - n * mStep;
        w_left = 1 - d_left / mStep;
        d_right = (int(n) + 1) * mStep - h;
        w_right = 1 - d_right / mStep;
        return w_left * mValues[n] + w_right * mValues[int(n) + 1];
    }


}
