#include "discamb/Scattering/SlaterTypeOrbitalScattering.h"
#include "discamb/MathUtilities/math_utilities.h"

namespace discamb {

namespace sto_scattering {

REAL scatteringSphericalDensity(
    const std::vector<REAL> &coefficients,
    const std::vector<REAL> &exponents,
    const std::vector<int> &powers,
    const REAL &h_length)
{
    REAL result = 0;

    for (int k = 0; k < coefficients.size(); k++) {
        REAL g_factor = gFunction<0>(powers[k] + 2, h_length, exponents[k]);
        result += coefficients[k] * g_factor;
    }
    return result * 4 * M_PI;
}

REAL gFunction(int l, const int n, REAL const h, REAL const Z)
{
    if (l == 0)
        return gFunction<0>(n, h, Z);
    if (l == 1)
        return gFunction<1>(n, h, Z);
    if (l == 2)
        return gFunction<2>(n, h, Z);
    if (l == 3)
        return gFunction<3>(n, h, Z);
    if (l == 4)
        return gFunction<4>(n, h, Z);
    return 0.0;

}

}
} // namespace discamb