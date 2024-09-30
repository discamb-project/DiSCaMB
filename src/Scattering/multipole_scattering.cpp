#include "discamb/Scattering/multipole_scattering.h"
#include "discamb/Scattering/SlaterTypeOrbitalScattering.h"
#include "discamb/MathUtilities/RealSphericalHarmonics.h"
#include "discamb/HC_Model/HC_WfnData.h"

#include <cmath>

using namespace std;

namespace discamb {

    namespace multipole_scattering {

        std::complex<double> hcDeformationValence(
            const std::vector<std::vector<double> > &p_lms,
            double kappa,
            double zeta,
            const std::vector<int> &powerR,
            const Vector3d &h)
        {
            return hcDeformationValence(p_lms, zeta, powerR, h/kappa);
        }

        std::complex<double> hcDeformationValence(
            const std::vector<std::vector<double> > &p_lms,
            double zeta,
            const std::vector<int> &powerR,
            const Vector3d &h)
        {
            if (p_lms.empty())
                return 0.0;

            double radial;
            complex<double> angular, result = 0.0;

            double hNorm = sqrt(h*h);
            Vector3d hNormalized = h / hNorm;

            int maxL = int(p_lms.size()) - 1;
            
            for (int l = 0; l <= maxL; l++)
            {
                radial = sto_scattering::gFunction(l, powerR[l]+2, hNorm, zeta) * sto_atomic_wfn::stoDensityNormalizationFactor(powerR[l], zeta);
                angular = multipoleScattering(l, p_lms[l], hNormalized, true);
                result += radial*angular;
            }

            return result;
        }




        std::complex<double> multipoleScattering(
            int l,
            int m,
            const Vector3d &h,
            bool densityNormalized)
        {
            double amplitude;
            
            densityNormalized ?
                amplitude = real_spherical_harmonics::densityNormalized(h, l, m) :
                amplitude = real_spherical_harmonics::wfnNormalized(h, l, m);

            amplitude *= 4.0 * M_PI;

            // accounts for i power l

            int imaginary = 1, negative = 2;
            if (l & negative)
                amplitude = -amplitude;
            if (l & imaginary)
                return complex<double>(0, amplitude);
            return amplitude;
        }


        std::complex<double> multipoleScattering(
            int l,
            const std::vector<double> &multipliers,
            const Vector3d &h,
            bool densityNormalized)
        {
            complex<double> result;

            for (int m = -l; m <= l; m++)
                result += multipliers[m+l] * multipoleScattering(l, m, h, densityNormalized);
            
            return result;
        }

        std::complex<double> multipoleScattering(
            const std::vector<std::vector<double> > &multipliers,
            const Vector3d &h,
            bool densityNormalized)
        {
            if (multipliers.empty())
                return 0.0;

            int maxL = int(multipliers.size() - 1);
            complex<double> result;

            for (int l = 0; l <= maxL; l++)
                result += multipoleScattering(l, multipliers[l], h, densityNormalized);
            
            return result;
        }

    }
}
