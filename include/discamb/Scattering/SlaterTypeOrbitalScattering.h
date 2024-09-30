#ifndef _DISCAMB_SCATTERING_SLATERTYPEORBITALSCATTERING_H_
#define _DISCAMB_SCATTERING_SLATERTYPEORBITALSCATTERING_H_

#include "Real.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include <vector>
#include <cassert>
#include <cstddef>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


namespace sto_scattering{

/**
calculate scattering factor from a spherical density given by linear combinatin of slater type functions
and centered at (0,0,0):
electron_density(r) = sum_{k} coeff[k] * r^pow[k] * exp( - exp[k]*r )
*/

REAL scatteringSphericalDensity(
    const std::vector<REAL> &coeff,
    const std::vector<REAL> &exp,
    const std::vector<int> &pow,
    const REAL &h_length);

/**
returns value of the following integral:
g(l,n,h,a) = Int_0^inf r^n exp(-a*r) j_l(2*pi*h*r) dr
where j_l is bessel function
*/

template<int l>
REAL gFunction( const int n, REAL const h, REAL const Z);



REAL gFunction(int l, const int n, REAL const h, REAL const Z);

}


//------------ IMPLEMENTATION -----------------


template<int l>
REAL sto_scattering::gFunction(
    const int n,
    REAL const h,
    REAL const Z)
{
    const REAL K = (REAL)(2.0*REAL(M_PI))*h; // K and Z symbols like in Coppens book

    const REAL K_pow2 = K*K;
    const REAL K_pow3 = K_pow2 * K;
    const REAL K_pow4 = K_pow2 * K_pow2;
    const REAL K_pow5 = K_pow2 * K_pow3;
    const REAL K_pow6 = K_pow4 * K_pow2;
    const REAL K_pow7 = K_pow6 * K;
    const REAL K_pow8 = K_pow4 * K_pow4;
    const REAL K_pow9 = K_pow8 * K;

    const REAL Z_pow2 = Z*Z;
    const REAL Z_pow3 = Z_pow2 * Z;
    const REAL Z_pow4 = Z_pow2 * Z_pow2;
    const REAL Z_pow5 = Z_pow2 * Z_pow3;
    const REAL Z_pow6 = Z_pow4 * Z_pow2;
    const REAL Z_pow7 = Z_pow6 * Z;
    const REAL Z_pow8 = Z_pow4 * Z_pow4;
    const REAL Z_pow9 = Z_pow8 * Z;

    const REAL d = K_pow2 + Z_pow2;

    const REAL d_inv = (REAL)1.0 / d;
    const REAL d_inv_pow2 = d_inv * d_inv;
    const REAL d_inv_pow4 = d_inv_pow2 * d_inv_pow2;
    const REAL d_inv_pow8 = d_inv_pow4 * d_inv_pow4;
    const REAL d_inv_pow3 = d_inv_pow2 * d_inv;
    const REAL d_inv_pow5 = d_inv_pow4 * d_inv;
    const REAL d_inv_pow6 = d_inv_pow4 * d_inv_pow2;
    const REAL d_inv_pow7 = d_inv_pow6 * d_inv;
    const REAL d_inv_pow9 = d_inv_pow8 * d_inv;
    const REAL d_inv_pow10 = d_inv_pow4 * d_inv_pow6;

    REAL value = (REAL)0;
    assert(l >= 0 && l <= 4);

    if (l == 0) {
        assert(n >= 2 && n <= 10);
        if (n == 2)
            value = 2 * Z*d_inv_pow2;
        else if (n == 3)
            value = 2 * (3 * Z_pow2 - K_pow2)*d_inv_pow3;
        else if (n == 4)
            value = 24 * Z*(Z_pow2 - K_pow2)*d_inv_pow4;
        else if (n == 5)
            value = 24 * (5 * Z_pow4 - 10 * K_pow2*Z_pow2 + K_pow4)*d_inv_pow5;//poprawionywspolczynnk!
        else if (n == 6)
            value = 240 * Z*(K_pow2 - 3 * Z_pow2)*(3 * K_pow2 - Z_pow2)*d_inv_pow6;
        else if (n == 7)
            value = 720 * (7 * Z_pow6 - 35 * K_pow2*Z_pow4 + 21 * K_pow4*Z_pow2 - K_pow6)*d_inv_pow7;
        else if (n == 8)
            value = 40320 * (Z_pow7 - 7 * K_pow2*Z_pow5 + 7 * K_pow4*Z_pow3 - K_pow6*Z)*d_inv_pow8;
        else if (n == 9)
            value = (362880 * Z_pow8 - 3386880 * K_pow2*Z_pow6 + 5080320 * K_pow4*Z_pow4 - 1451520 * K_pow6*Z_pow2 + 40320 * K_pow8)*d_inv_pow9;
        else if (n == 10)
            value = (3628800 * Z_pow9 - 43545600 * K_pow2*Z_pow7 + 91445760 * K_pow4*Z_pow5 - 43545600 * K_pow6*Z_pow3 + 3628800 * K_pow8*Z)*d_inv_pow10;
    }

    if (l == 1) {
        assert(n >= 3 && n <= 10);
        if (n == 3)
            value = 8 * K*Z*d_inv_pow3;
        else if (n == 4)
            value = 8 * K*(5 * Z_pow2 - K_pow2)*d_inv_pow4;
        else if (n == 5)
            value = 48 * K*Z*(5 * Z_pow2 - 3 * K_pow2)*d_inv_pow5;
        else if (n == 6)
            value = 48 * K*(35 * Z_pow4 - 42 * K_pow2*Z_pow2 + 3 * K_pow4)*d_inv_pow6;
        else if (n == 7)
            value = 1920 * K*Z*(7 * Z_pow4 - 14 * K_pow2*Z_pow2 + 3 * K_pow4)*d_inv_pow7;
        else if (n == 8)
            value = 5760 * K*(21 * Z_pow6 - 63 * K_pow2*Z_pow4 + 27 * K_pow4*Z_pow2 - K_pow6)*d_inv_pow8;
        else if (n == 9)
            value = (1209600 * K*Z_pow7 - 5080320 * K_pow3*Z_pow5 + 3628800 * K_pow5*Z_pow3 - 403200 * K_pow7*Z)*d_inv_pow9;
        else if (n == 10)
            value = (13305600 * K*Z_pow8 - 74511360 * K_pow3*Z_pow6 + 79833600 * K_pow5*Z_pow4 - 17740800 * K_pow7*Z_pow2 + 403200 * K_pow9)*d_inv_pow10;

    }

    if (l == 2) {
        assert(n >= 4 && n <= 10);
        if (n == 4)
            value = 48 * K_pow2*Z*d_inv_pow4;
        else if (n == 5)
            value = 48 * K_pow2*(7 * Z_pow2 - K_pow2)*d_inv_pow5;
        else if (n == 6)
            value = 384 * K_pow2*Z*(7 * Z_pow2 - 3 * K_pow2)*d_inv_pow6;
        else if (n == 7)
            value = 1152 * K_pow2*(21 * Z_pow4 - 18 * K_pow2*Z_pow2 + K_pow4)*d_inv_pow7;
        else if (n == 8)
            value = 11520 * K_pow2*Z*(21 * Z_pow4 - 30 * K_pow2*Z_pow2 + 5 * K_pow4)*d_inv_pow8;
        else if (n == 9)
            value = (2661120 * K_pow2*Z_pow6 - 5702400 * K_pow4*Z_pow4 + 1900800 * K_pow6*Z_pow2 - 57600 * K_pow8)*d_inv_pow9;
        else if (n == 10)
            value = (31933440 * K_pow2*Z_pow7 - 95800320 * K_pow4*Z_pow5 + 53222400 * K_pow6*Z_pow3 - 4838400 * K_pow8*Z)*d_inv_pow10;
    }
    if (l == 3) {
        assert(n >= 5 && n <= 10);
        if (n == 5)
            value = (384 * K_pow3*Z)*d_inv_pow5;
        else if (n == 6)
            value = (3456 * K_pow3*Z_pow2 - 384 * K_pow5)*d_inv_pow6;
        else if (n == 7)
            value = (34560 * K_pow3*Z_pow3 - 11520 * K_pow5*Z)*d_inv_pow7;
        else if (n == 8)
            value = (380160 * K_pow3*Z_pow4 - 253440 * K_pow5*Z_pow2 + 11520 * K_pow7)*d_inv_pow8;
        else if (n == 9)
            value = (4561920 * K_pow3*Z_pow5 - 5068800 * K_pow5*Z_pow3 + 691200 * K_pow7*Z)*d_inv_pow9;
        else if (n == 10)
            value = (59304960 * K_pow3*Z_pow6 - 98841600 * K_pow5*Z_pow4 + 26956800 * K_pow7*Z_pow2 - 691200 * K_pow9)*d_inv_pow10;
    }
    if (l == 4) {
        assert(n >= 6 && n <= 10);
        if (n == 6)
            value = (3840 * K_pow4*Z)*d_inv_pow6;
        else if (n == 7)
            value = (42240 * K_pow4*Z_pow2 - 3840 * K_pow6)*d_inv_pow7;
        else if (n == 8)
            value = (506880 * K_pow4*Z_pow3 - 138240 * K_pow6*Z)*d_inv_pow8;
        else if (n == 9)
            value = (6589440 * K_pow4*Z_pow4 - 3594240 * K_pow6*Z_pow2 + 138240 * K_pow8)*d_inv_pow9;
        else if (n == 10)
            value = (92252160 * K_pow4*Z_pow5 - 83865600 * K_pow6*Z_pow3 + 9676800 * K_pow8*Z)*d_inv_pow10;
    }

    return value;
}

/** @}*/

} // namespace discamb
#endif /*_DISCAMB_SCATTERING_SLATERTYPEORBITALSCATTERING_H_*/
