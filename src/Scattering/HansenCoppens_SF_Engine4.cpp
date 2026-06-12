#include "discamb/Scattering/HansenCoppens_SF_Engine4.h"

#include <algorithm>
#include <cassert>
#include <cereal/archives/binary.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <ostream>
#include <vector>

#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/Timer.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/HC_Model/HC_WfnData.h"
#include "discamb/MathUtilities/SphConverter.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/Scattering/NGaussianFormFactorsTable.h"
#include "discamb/Scattering/SlaterTypeOrbitalScattering.h"
#include "discamb/Scattering/scattering_utilities.h"

#if defined(_OPENMP)
#include <omp.h>
#endif

#ifdef SYCL
#include <sycl/sycl.hpp>

#include "sycl/access/access.hpp"
#include "sycl/buffer.hpp"
#include "sycl/detail/builtins/builtins.hpp"
#include "sycl/exception.hpp"

namespace vecnd_detail {

template <typename T>
struct is_list : std::false_type {};

template <typename T, typename Alloc>
struct is_list<std::vector<T, Alloc>> : std::true_type {};

template <typename T, std::size_t N>
struct is_list<std::array<T, N>> : std::true_type {};

template <typename T>
struct list_depth : std::integral_constant<std::size_t, 0> {};

template <typename T, typename Alloc>
struct list_depth<std::vector<T, Alloc>>
    : std::integral_constant<std::size_t, 1 + list_depth<T>::value> {};

template <typename T, std::size_t N>
struct list_depth<std::array<T, N>>
    : std::integral_constant<std::size_t, 1 + list_depth<T>::value> {};

template <typename T>
struct deepest_value_type {
    using type = T;
};

template <typename T, typename Alloc>
struct deepest_value_type<std::vector<T, Alloc>> {
    using type = typename deepest_value_type<T>::type;
};

template <typename T, std::size_t N>
struct deepest_value_type<std::array<T, N>> {
    using type = typename deepest_value_type<T>::type;
};

template <typename Vec, typename Out>
inline void flatten_impl(const Vec &list, Out &out) {
    if constexpr (is_list<typename Vec::value_type>::value) {
        for (const auto &sub : list) flatten_impl(sub, out);
    } else {
        for (const auto &el : list) out.push_back(el);
    }
}

template <typename Vec>
inline auto flatten_all(const Vec &list) {
    using value_type = typename deepest_value_type<Vec>::type;
    std::vector<value_type> out;
    flatten_impl(list, out);
    return out;
}

template <typename Vec, std::size_t N>
inline void offsets_impl(const Vec &list,
                         std::array<std::vector<int>, N> &offsets,
                         std::array<int, N> &running, std::size_t level = 0) {
    using Inner = typename Vec::value_type;
    if constexpr (is_list<Inner>::value) {
        for (const auto &sub : list) {
            offsets[level].push_back(running[level]);
            running[level] += static_cast<int>(sub.size());
            offsets_impl(sub, offsets, running, level + 1);
        }
    }
}

template <std::size_t Depth>
struct offsets_data {
    std::vector<int> flat;
    std::array<int, Depth - 1> level_offsets;
};

template <typename Vec>
inline auto build_offsets(const Vec &list) {
    constexpr std::size_t Depth = list_depth<Vec>::value;
    static_assert(Depth >= 2,
                  "vecnd_to_buffers requires nested std::vector depth >= 2");
    offsets_data<Depth> result;
    std::array<std::vector<int>, Depth - 1> per_level;
    std::array<int, Depth - 1> running{};
    offsets_impl(list, per_level, running, 0);
    int total = 0;
    for (std::size_t level = 0; level < Depth - 1; ++level) {
        result.level_offsets[level] = total;
        total += static_cast<int>(per_level[level].size());
    }
    result.flat.reserve(total);
    for (std::size_t level = 0; level < Depth - 1; ++level)
        result.flat.insert(result.flat.end(),
                           per_level[level].begin(),
                           per_level[level].end());
    return result;
}

template <typename OffsetAccessor, std::size_t N>
inline int index_flat(const OffsetAccessor &offsets_ax,
                      const std::array<int, N> &level_offsets,
                      const std::array<int, N + 1> &indices) {
    int idx = indices[0];
    for (std::size_t level = 0; level < N; ++level)
        idx = offsets_ax[level_offsets[level] + idx] + indices[level + 1];
    return idx;
}

}  // namespace vecnd_detail

#define vecnd_to_buffers(list, buf_name_prefix)                              \
    using buf_name_prefix##_list_type = std::decay_t<decltype(list)>;        \
    constexpr std::size_t buf_name_prefix##_depth =                          \
        vecnd_detail::list_depth<buf_name_prefix##_list_type>::value;        \
    auto buf_name_prefix##_value_vec = vecnd_detail::flatten_all(list);      \
    sycl::buffer<typename vecnd_detail::deepest_value_type<                  \
        buf_name_prefix##_list_type>::type>                                  \
        buf_name_prefix##_value_buf(buf_name_prefix##_value_vec);            \
    auto buf_name_prefix##_offsets_data = vecnd_detail::build_offsets(list); \
    auto buf_name_prefix##_level_offsets =                                   \
        buf_name_prefix##_offsets_data.level_offsets;                        \
    sycl::buffer<int> buf_name_prefix##_offsets_buf(                         \
        buf_name_prefix##_offsets_data.flat);

#define vecnd_buffer_accessors(buf_name_prefix)             \
    sycl::accessor buf_name_prefix##_value_ax(              \
        buf_name_prefix##_value_buf, cgh, sycl::read_only); \
    sycl::accessor buf_name_prefix##_offset_ax(             \
        buf_name_prefix##_offsets_buf, cgh, sycl::read_only);

#define index_vecnd_buffer(buf_name_prefix, ...)          \
    (buf_name_prefix##_value_ax[vecnd_detail::index_flat( \
        buf_name_prefix##_offset_ax,                      \
        buf_name_prefix##_level_offsets,                  \
        std::array<int, buf_name_prefix##_depth>{__VA_ARGS__})])
#endif

#include <array>
#include <ctime>
#include <iostream>

using namespace std;

namespace discamb {

#ifdef SYCL
template <int l>
inline REAL gFunction_sycl_impl(const int n, REAL const h, REAL const Z) {
    const REAL K =
        (REAL)(2.0 * REAL(M_PI)) * h;  // K and Z symbols like in Coppens book

    const REAL K_pow2 = K * K;
    const REAL K_pow3 = K_pow2 * K;
    const REAL K_pow4 = K_pow2 * K_pow2;
    const REAL K_pow5 = K_pow2 * K_pow3;
    const REAL K_pow6 = K_pow4 * K_pow2;
    const REAL K_pow7 = K_pow6 * K;
    const REAL K_pow8 = K_pow4 * K_pow4;
    const REAL K_pow9 = K_pow8 * K;

    const REAL Z_pow2 = Z * Z;
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
            value = 2 * Z * d_inv_pow2;
        else if (n == 3)
            value = 2 * (3 * Z_pow2 - K_pow2) * d_inv_pow3;
        else if (n == 4)
            value = 24 * Z * (Z_pow2 - K_pow2) * d_inv_pow4;
        else if (n == 5)
            value = 24 * (5 * Z_pow4 - 10 * K_pow2 * Z_pow2 + K_pow4) *
                    d_inv_pow5;  // poprawionywspolczynnk!
        else if (n == 6)
            value = 240 * Z * (K_pow2 - 3 * Z_pow2) * (3 * K_pow2 - Z_pow2) *
                    d_inv_pow6;
        else if (n == 7)
            value = 720 *
                    (7 * Z_pow6 - 35 * K_pow2 * Z_pow4 + 21 * K_pow4 * Z_pow2 -
                     K_pow6) *
                    d_inv_pow7;
        else if (n == 8)
            value = 40320 *
                    (Z_pow7 - 7 * K_pow2 * Z_pow5 + 7 * K_pow4 * Z_pow3 -
                     K_pow6 * Z) *
                    d_inv_pow8;
        else if (n == 9)
            value = (362880 * Z_pow8 - 3386880 * K_pow2 * Z_pow6 +
                     5080320 * K_pow4 * Z_pow4 - 1451520 * K_pow6 * Z_pow2 +
                     40320 * K_pow8) *
                    d_inv_pow9;
        else if (n == 10)
            value = (3628800 * Z_pow9 - 43545600 * K_pow2 * Z_pow7 +
                     91445760 * K_pow4 * Z_pow5 - 43545600 * K_pow6 * Z_pow3 +
                     3628800 * K_pow8 * Z) *
                    d_inv_pow10;
    }

    if (l == 1) {
        assert(n >= 3 && n <= 10);
        if (n == 3)
            value = 8 * K * Z * d_inv_pow3;
        else if (n == 4)
            value = 8 * K * (5 * Z_pow2 - K_pow2) * d_inv_pow4;
        else if (n == 5)
            value = 48 * K * Z * (5 * Z_pow2 - 3 * K_pow2) * d_inv_pow5;
        else if (n == 6)
            value = 48 * K * (35 * Z_pow4 - 42 * K_pow2 * Z_pow2 + 3 * K_pow4) *
                    d_inv_pow6;
        else if (n == 7)
            value = 1920 * K * Z *
                    (7 * Z_pow4 - 14 * K_pow2 * Z_pow2 + 3 * K_pow4) *
                    d_inv_pow7;
        else if (n == 8)
            value = 5760 * K *
                    (21 * Z_pow6 - 63 * K_pow2 * Z_pow4 + 27 * K_pow4 * Z_pow2 -
                     K_pow6) *
                    d_inv_pow8;
        else if (n == 9)
            value = (1209600 * K * Z_pow7 - 5080320 * K_pow3 * Z_pow5 +
                     3628800 * K_pow5 * Z_pow3 - 403200 * K_pow7 * Z) *
                    d_inv_pow9;
        else if (n == 10)
            value = (13305600 * K * Z_pow8 - 74511360 * K_pow3 * Z_pow6 +
                     79833600 * K_pow5 * Z_pow4 - 17740800 * K_pow7 * Z_pow2 +
                     403200 * K_pow9) *
                    d_inv_pow10;
    }

    if (l == 2) {
        assert(n >= 4 && n <= 10);
        if (n == 4)
            value = 48 * K_pow2 * Z * d_inv_pow4;
        else if (n == 5)
            value = 48 * K_pow2 * (7 * Z_pow2 - K_pow2) * d_inv_pow5;
        else if (n == 6)
            value = 384 * K_pow2 * Z * (7 * Z_pow2 - 3 * K_pow2) * d_inv_pow6;
        else if (n == 7)
            value = 1152 * K_pow2 *
                    (21 * Z_pow4 - 18 * K_pow2 * Z_pow2 + K_pow4) * d_inv_pow7;
        else if (n == 8)
            value = 11520 * K_pow2 * Z *
                    (21 * Z_pow4 - 30 * K_pow2 * Z_pow2 + 5 * K_pow4) *
                    d_inv_pow8;
        else if (n == 9)
            value = (2661120 * K_pow2 * Z_pow6 - 5702400 * K_pow4 * Z_pow4 +
                     1900800 * K_pow6 * Z_pow2 - 57600 * K_pow8) *
                    d_inv_pow9;
        else if (n == 10)
            value = (31933440 * K_pow2 * Z_pow7 - 95800320 * K_pow4 * Z_pow5 +
                     53222400 * K_pow6 * Z_pow3 - 4838400 * K_pow8 * Z) *
                    d_inv_pow10;
    }
    if (l == 3) {
        assert(n >= 5 && n <= 10);
        if (n == 5)
            value = (384 * K_pow3 * Z) * d_inv_pow5;
        else if (n == 6)
            value = (3456 * K_pow3 * Z_pow2 - 384 * K_pow5) * d_inv_pow6;
        else if (n == 7)
            value = (34560 * K_pow3 * Z_pow3 - 11520 * K_pow5 * Z) * d_inv_pow7;
        else if (n == 8)
            value = (380160 * K_pow3 * Z_pow4 - 253440 * K_pow5 * Z_pow2 +
                     11520 * K_pow7) *
                    d_inv_pow8;
        else if (n == 9)
            value = (4561920 * K_pow3 * Z_pow5 - 5068800 * K_pow5 * Z_pow3 +
                     691200 * K_pow7 * Z) *
                    d_inv_pow9;
        else if (n == 10)
            value = (59304960 * K_pow3 * Z_pow6 - 98841600 * K_pow5 * Z_pow4 +
                     26956800 * K_pow7 * Z_pow2 - 691200 * K_pow9) *
                    d_inv_pow10;
    }
    if (l == 4) {
        assert(n >= 6 && n <= 10);
        if (n == 6)
            value = (3840 * K_pow4 * Z) * d_inv_pow6;
        else if (n == 7)
            value = (42240 * K_pow4 * Z_pow2 - 3840 * K_pow6) * d_inv_pow7;
        else if (n == 8)
            value =
                (506880 * K_pow4 * Z_pow3 - 138240 * K_pow6 * Z) * d_inv_pow8;
        else if (n == 9)
            value = (6589440 * K_pow4 * Z_pow4 - 3594240 * K_pow6 * Z_pow2 +
                     138240 * K_pow8) *
                    d_inv_pow9;
        else if (n == 10)
            value = (92252160 * K_pow4 * Z_pow5 - 83865600 * K_pow6 * Z_pow3 +
                     9676800 * K_pow8 * Z) *
                    d_inv_pow10;
    }

    return value;
}

inline REAL gFunction_sycl(int l, const int n, REAL const h, REAL const Z) {
    if (l == 0) return gFunction_sycl_impl<0>(n, h, Z);
    if (l == 1) return gFunction_sycl_impl<1>(n, h, Z);
    if (l == 2) return gFunction_sycl_impl<2>(n, h, Z);
    if (l == 3) return gFunction_sycl_impl<3>(n, h, Z);
    if (l == 4) return gFunction_sycl_impl<4>(n, h, Z);
    return 0.0;
}

double polynomialSycl(const sycl::vec<REAL, 3> &v,  // noralized 3D vector,
                      int l, int m) {
    REAL x = v[0];
    REAL y = v[1];
    REAL z = v[2];
    switch (l) {
        case 0:
            return 1.0;
        case 1:
            if (m == -1)
                return y;
            else {
                if (m == 0)
                    return z;
                else
                    return x;
            }
            return 0;
        case 2:
            switch (m) {
                case -2:
                    return x * y;
                case -1:
                    return y * z;
                case 0:
                    return 3 * z * z - 1;
                case 1:
                    return x * z;
                case 2:
                    return x * x - y * y;
                default:
                    return 0;
            }
        case 3:
            switch (m) {
                case -3:
                    return (3 * x * x - y * y) * y;
                case -2:
                    return x * y * z;
                case -1:
                    return y * (5 * z * z - 1);
                case 0:
                    return z * (5 * z * z - 3);
                case 1:
                    return x * (5 * z * z - 1);
                case 2:
                    return (x * x - y * y) * z;
                case 3:
                    return (x * x - 3 * y * y) * x;
                default:
                    return 0;
            }
        case 4:
            switch (m) {
                case -4:
                    return x * y * (x * x - y * y);
                case -3:
                    return (3 * x * x - y * y) * y * z;
                case -2:
                    return x * y * (7 * z * z - 1);
                case -1:
                    return y * z * (7 * z * z - 3);
                case 0:
                    return z * z * (35 * z * z - 30) + 3;
                case 1:
                    return x * z * (7 * z * z - 3);
                case 2:
                    return (x * x - y * y) * (7 * z * z - 1);
                case 3:
                    return (x * x - 3 * y * y) * x * z;
                case 4:
                    return x * x * (x * x - 3 * y * y) -
                           y * y * (3 * x * x - y * y);
                default:
                    return 0;
            }
        default:
            return 0;
    }  // maxL>0
    return 0.0;
}

double densityNormalizedSycl(const sycl::vec<REAL, 3> &normalizedVector3D,
                             int l, int m) {
    return real_spherical_harmonics::densityNormalizationMultipliers[l][l + m] *
           polynomialSycl(normalizedVector3D, l, m);
}

#endif

HansenCoppens_SF_Engine4::HansenCoppens_SF_Engine4() { mUseIAM = false; }

HansenCoppens_SF_Engine4::~HansenCoppens_SF_Engine4() {}

inline void HansenCoppens_SF_Engine4::add_contribution_to_occupancy_derivative(
    REAL &occupancy_derivative, const complex<REAL> &dTarget_dF,
    const complex<REAL> &atomic_f_divided_by_occupancy) {
    occupancy_derivative += (dTarget_dF * atomic_f_divided_by_occupancy).real();
}

inline void HansenCoppens_SF_Engine4::add_contribution_to_position_derivatives(
    Vector3<REAL> &position_derivatives, const complex<REAL> dTarget_dF,
    const complex<REAL> &atomic_f, const Vector3<REAL> &h) {
    static const complex<REAL> two_pi_i =
        REAL(2 * REAL(M_PI)) * complex<REAL>(0, 1);
    complex<REAL> df_dparam;

    for (int k = 0; k < 3; k++) {
        df_dparam = two_pi_i * h[k] * atomic_f;
        position_derivatives[k] += (dTarget_dF * df_dparam).real();
    }
}

inline void HansenCoppens_SF_Engine4::add_contribution_to_adp_derivatives(
    std::vector<std::complex<REAL>> &adp_derivatives,
    const std::complex<REAL> &dTarget_dF, const std::complex<REAL> &atomic_f,
    const Vector3<REAL> &h) {
    complex<REAL> df_dparam;
    REAL hVectorLength = sqrt(h * h);

    if (adp_derivatives.size() == 1) {
        df_dparam = -hVectorLength * hVectorLength * atomic_f;
        adp_derivatives[0] += dTarget_dF * df_dparam;

    } else {
        for (int k = 0; k < 3; k++) {
            df_dparam = -h[k] * h[k] * atomic_f;
            adp_derivatives[k] += dTarget_dF * df_dparam;
        }

        // U_12
        df_dparam = -2 * h[0] * h[1] * atomic_f;
        adp_derivatives[3] += dTarget_dF * df_dparam;

        // U_13
        df_dparam = -2 * h[0] * h[2] * atomic_f;
        adp_derivatives[4] += dTarget_dF * df_dparam;

        // U_23
        df_dparam = -2 * h[1] * h[2] * atomic_f;
        adp_derivatives[5] += dTarget_dF * df_dparam;
    }
}

inline void HansenCoppens_SF_Engine4::process_adp_derivatives(
    std::complex<REAL> *pre_derivatives, const std::complex<REAL> &atomic_f,
    const Vector3<REAL> &h, REAL h_length, int n_adp_components) {
    if (n_adp_components == 1) {
        pre_derivatives[0] -= h_length * h_length * atomic_f;
        return;
    }

    pre_derivatives[0] -= h[0] * h[0] * atomic_f;
    pre_derivatives[1] -= h[1] * h[1] * atomic_f;
    pre_derivatives[2] -= h[2] * h[2] * atomic_f;
    pre_derivatives[3] -= 2 * h[0] * h[1] * atomic_f;
    pre_derivatives[4] -= 2 * h[0] * h[2] * atomic_f;
    pre_derivatives[5] -= 2 * h[1] * h[2] * atomic_f;
}

std::complex<double> HansenCoppens_SF_Engine4::calculateDeformationValence(
    const std::vector<std::vector<REAL>>
        &p_lm,  // coefficients for multipolar terms (with wavefunction
                // normalization of spherical harmonics)
    const std::vector<REAL> &g_functions_and_slater_normalization,
    // const Matrix3<REAL>& local_coordinates_system,
    int max_l, std::vector<std::vector<double>> &sphericalHarmonics) {
    if (max_l < 0) return 0;

    switch (max_l) {
        case 0:
            return combine_multipolar_terms<0>(
                p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
        case 1:
            return combine_multipolar_terms<1>(
                p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
        case 2:
            return combine_multipolar_terms<2>(
                p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
        case 3:
            return combine_multipolar_terms<3>(
                p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
        case 4:
            return combine_multipolar_terms<4>(
                p_lm, g_functions_and_slater_normalization, sphericalHarmonics);
        default:
            return 0;
    }
}

std::complex<REAL> HansenCoppens_SF_Engine4::calculateDeformationValence(
    const std::vector<std::vector<REAL>> &p_lm,
    const std::vector<REAL> &g_functions_and_slater_normalization,
    const Matrix3<REAL> &local_coordinates_system,
    const Vector3<REAL> &normalized_h_vector, int max_l,
    std::vector<std::vector<double>> &sphericalHarmonicBuffer) {
    if (max_l < 0) return 0;

    const Matrix3<REAL> &lcs = local_coordinates_system;
    const REAL x = (lcs(0, 0) * normalized_h_vector(0) +
                    lcs(1, 0) * normalized_h_vector(1) +
                    lcs(2, 0) * normalized_h_vector(2));  // hRotated[0];
    const REAL y = (lcs(0, 1) * normalized_h_vector(0) +
                    lcs(1, 1) * normalized_h_vector(1) +
                    lcs(2, 1) * normalized_h_vector(2));  // hRotated[1];
    const REAL z = (lcs(0, 2) * normalized_h_vector(0) +
                    lcs(1, 2) * normalized_h_vector(1) +
                    lcs(2, 2) * normalized_h_vector(2));  // hRotated[2];

    Vector3d h(x, y, z);

    switch (max_l) {
        case 0:
            real_spherical_harmonics::getDensityNormalized<0>(
                h, sphericalHarmonicBuffer);
            return combine_multipolar_terms<0>(
                p_lm,
                g_functions_and_slater_normalization,
                sphericalHarmonicBuffer);
        case 1:
            real_spherical_harmonics::getDensityNormalized<1>(
                h, sphericalHarmonicBuffer);
            return combine_multipolar_terms<1>(
                p_lm,
                g_functions_and_slater_normalization,
                sphericalHarmonicBuffer);
        case 2:
            real_spherical_harmonics::getDensityNormalized<2>(
                h, sphericalHarmonicBuffer);
            return combine_multipolar_terms<2>(
                p_lm,
                g_functions_and_slater_normalization,
                sphericalHarmonicBuffer);
        case 3:
            real_spherical_harmonics::getDensityNormalized<3>(
                h, sphericalHarmonicBuffer);
            return combine_multipolar_terms<3>(
                p_lm,
                g_functions_and_slater_normalization,
                sphericalHarmonicBuffer);
        case 4:
            real_spherical_harmonics::getDensityNormalized<4>(
                h, sphericalHarmonicBuffer);
            return combine_multipolar_terms<4>(
                p_lm,
                g_functions_and_slater_normalization,
                sphericalHarmonicBuffer);
        default:
            return 0;
    }
}

void HansenCoppens_SF_Engine4::pre_hkl_loop_sf_calc(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map, std::vector<int> &type_2_wfn_type,
    std::vector<std::vector<REAL>> &def_val_slater_normalization,
    std::vector<int> &typeMaxL) {
    int nWfnTypes = wfn_parameters.size();
    int nTypes = type_parameters.size();
    int nAtoms = atom_to_wfn_map.size();
    int i, j, nL;

    if (mUseIAM) {
        typeMaxL.assign(nTypes, -1);
        return;
    }

    type_2_wfn_type.resize(nTypes);
    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
        int atomWfnIdx = atom_to_wfn_map[atomIdx];
        int atomTypeIdx = atom_to_type_map[atomIdx];
        type_2_wfn_type[atomTypeIdx] = atomWfnIdx;
    }

    def_val_slater_normalization.resize(nWfnTypes);

    for (i = 0; i < nWfnTypes; i++) {
        nL = wfn_parameters[i].def_valence_pow.size();
        def_val_slater_normalization[i].resize(nL);
        for (j = 0; j < nL; j++)
            def_val_slater_normalization[i][j] =
                sto_atomic_wfn::stoDensityNormalizationFactor(
                    wfn_parameters[i].def_valence_pow[j],
                    wfn_parameters[i].def_valence_exp);
    }

    typeMaxL.resize(nTypes);
    int maxL_FromPlm;

    for (i = 0; i < nTypes; i++) {
        maxL_FromPlm = -1;
        for (int l = 0; l < type_parameters[i].p_lm.size(); l++) {
            for (j = 0; j < 2 * l + 1; j++)
                if (type_parameters[i].p_lm[l][j] != 0.0) maxL_FromPlm = int(l);
        }

        typeMaxL[i] = std::min(
            4,
            int(wfn_parameters[type_2_wfn_type[i]].def_valence_pow.size()) - 1);
        typeMaxL[i] = std::min(typeMaxL[i], maxL_FromPlm);
    }
}

void HansenCoppens_SF_Engine4::calculateSF_IAM(
    const UnitCell &unitCell, const std::vector<std::string> &atomicType,
    const std::vector<std::complex<REAL>> &atomTypeAnomalousScattering,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL>> &atomicPositions,
    const std::vector<std::vector<REAL>> &atomic_displacement_parameters,
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<std::complex<REAL>> &anomalous_dispersion,
    const std::vector<REAL> &atomic_multiplicity_factor,
    const std::vector<sf_engine_data_types::SymmetryOperation>
        &symmetryOperations,
    bool centrosymmetric, const Vector3<REAL> &inversionTranslation,
    const std::vector<Vector3<REAL>> &hVectors,
    const std::vector<Vector3i> &hkl_indices,
    std::vector<std::complex<REAL>> &f,
    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
    const std::vector<std::complex<REAL>> &dTarget_df,
    const std::vector<bool> &include_atom_contribution, int nThreads) {
    mUseIAM = true;
    mIamAtomType = atomicType;
    mAtomToIamTypeMap = atom_to_type_map;

    vector<sf_engine_data_types::HC_WfnParam> wfnParams(atomicType.size());
    vector<sf_engine_data_types::HC_TypeParam> typeParams(1);
    vector<int> atomToWfnMap = atom_to_type_map;
    vector<int> atomToTypeMap(atomicPositions.size(), 0);
    Matrix3d idenity;
    idenity.setToIdentity();
    vector<Matrix3d> localCoordinateSystems(atomicPositions.size(), idenity);

    int iamTypeIdx, nIamTypes = atomicType.size();

    mIamFormFactors.resize(nIamTypes);

    for (iamTypeIdx = 0; iamTypeIdx < nIamTypes; iamTypeIdx++) {
        wfnParams[iamTypeIdx].anomalous_scattering =
            atomTypeAnomalousScattering[iamTypeIdx];
        if (n_gaussian_form_factors_table::hasFormFactor(
                mIamAtomType[iamTypeIdx]))
            mIamFormFactors[iamTypeIdx] =
                n_gaussian_form_factors_table::getFormFactor(
                    mIamAtomType[iamTypeIdx]);
        else
            on_error::throwException(string("request for Gaussian type atomic "
                                            "form factor parameter for "
                                            "unknown atom type: ") +
                                         mIamAtomType[iamTypeIdx],
                                     __FILE__,
                                     __LINE__);
    }

    DerivativesSelector derivativesSwitch;
    calculateSF(unitCell,
                wfnParams,
                typeParams,
                atomToWfnMap,
                atomToTypeMap,
                atomicPositions,
                atomic_displacement_parameters,
                atomic_occupancy,
                anomalous_dispersion,
                atomic_multiplicity_factor,
                localCoordinateSystems,
                symmetryOperations,
                centrosymmetric,
                inversionTranslation,
                hVectors,
                hkl_indices,
                f,
                dTarget_dparam,
                dTarget_df,
                include_atom_contribution,
                nThreads,
                derivativesSwitch);
}

void HansenCoppens_SF_Engine4::pre_atom_loop_sf_calc(
    // in:
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,
    const std::vector<sf_engine_data_types::HC_TypeParam> &typeParams,
    const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
    const std::vector<int> &type_2_wfn,
    const std::vector<std::vector<REAL>> &def_val_slater_normalization,
    const Vector3<REAL> &hVector, REAL hVectorLength,
    // out:
    vector<REAL> &wfn_spherical_core_sf, vector<REAL> &wfn_spherical_valence_sf,
    vector<vector<REAL>> &g_functions_and_slater_norm,
    vector<Vector3<REAL>> &rotated_h,
    vector<Vector3<REAL>> &rotated_normalized_h,
    std::vector<REAL> &translation_factor,
    std::vector<std::vector<REAL>> &adp_multipliers) {
    for (int symmOpIdx = 0; symmOpIdx < symOps.size(); symmOpIdx++) {
        translation_factor[symmOpIdx] = hVector * symOps[symmOpIdx].translation;
        rotated_h[symmOpIdx] = hVector * symOps[symmOpIdx].rotation;
        rotated_normalized_h[symmOpIdx] = rotated_h[symmOpIdx] / hVectorLength;

        // sets mAdpMultipliers
        Vector3<REAL> &h = rotated_h[symmOpIdx];
        REAL *adpMultipliers = &adp_multipliers[symmOpIdx][0];

        adpMultipliers[0] = h.x * h.x;
        adpMultipliers[1] = h.y * h.y;
        adpMultipliers[2] = h.z * h.z;
        adpMultipliers[3] = 2.0 * h.x * h.y;
        adpMultipliers[4] = 2.0 * h.x * h.z;
        adpMultipliers[5] = 2.0 * h.y * h.z;
    }

    if (mUseIAM) {
        for (int i = 0, n = wfn_spherical_core_sf.size(); i < n; i++)
            wfn_spherical_core_sf[i] =
                mIamFormFactors[i].calculate_h(hVectorLength);
        return;
    } else
        for (int wfnTypeIdx = 0; wfnTypeIdx < wfnParams.size(); wfnTypeIdx++)
            wfn_spherical_core_sf[wfnTypeIdx] =
                sto_scattering::scatteringSphericalDensity(
                    wfnParams[wfnTypeIdx].core_coeff,
                    wfnParams[wfnTypeIdx].core_exp,
                    wfnParams[wfnTypeIdx].core_pow,
                    hVectorLength);

    int nTypes = typeParams.size();

    for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        int wfnTypeIdx = type_2_wfn[typeIdx];

        wfn_spherical_valence_sf[typeIdx] =
            sto_scattering::scatteringSphericalDensity(
                wfnParams[wfnTypeIdx].valence_coeff,
                wfnParams[wfnTypeIdx].valence_exp,
                wfnParams[wfnTypeIdx].valence_pow,
                hVectorLength / typeParams[typeIdx].kappa_spherical);

        int nL = wfnParams[wfnTypeIdx].def_valence_pow.size();

        const vector<int> &def_valence_pow =
            wfnParams[wfnTypeIdx].def_valence_pow;

        if (nL > 0)
            g_functions_and_slater_norm[typeIdx][0] =
                def_val_slater_normalization[wfnTypeIdx][0] *
                sto_scattering::gFunction<0>(
                    int(def_valence_pow[0]) + 2,
                    hVectorLength / typeParams[typeIdx].kappa_def_valence,
                    wfnParams[wfnTypeIdx].def_valence_exp);
        if (nL > 1)
            g_functions_and_slater_norm[typeIdx][1] =
                def_val_slater_normalization[wfnTypeIdx][1] *
                sto_scattering::gFunction<1>(
                    int(def_valence_pow[1]) + 2,
                    hVectorLength / typeParams[typeIdx].kappa_def_valence,
                    wfnParams[wfnTypeIdx].def_valence_exp);

        if (nL > 2)
            g_functions_and_slater_norm[typeIdx][2] =
                def_val_slater_normalization[wfnTypeIdx][2] *
                sto_scattering::gFunction<2>(
                    int(def_valence_pow[2]) + 2,
                    hVectorLength / typeParams[typeIdx].kappa_def_valence,
                    wfnParams[wfnTypeIdx].def_valence_exp);
        if (nL > 3)
            g_functions_and_slater_norm[typeIdx][3] =
                def_val_slater_normalization[wfnTypeIdx][3] *
                sto_scattering::gFunction<3>(
                    int(def_valence_pow[3]) + 2,
                    hVectorLength / typeParams[typeIdx].kappa_def_valence,
                    wfnParams[wfnTypeIdx].def_valence_exp);
        if (nL > 4)
            g_functions_and_slater_norm[typeIdx][4] =
                def_val_slater_normalization[wfnTypeIdx][4] *
                sto_scattering::gFunction<4>(
                    int(def_valence_pow[4]) + 2,
                    hVectorLength / typeParams[typeIdx].kappa_def_valence,
                    wfnParams[wfnTypeIdx].def_valence_exp);
    }
}

void HansenCoppens_SF_Engine4::calculateFormFactors(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<double>
        &f_spherical,  // for each type spherical valence + core
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Matrix3<REAL>> &local_coordinate_systems,
    const Vector3<REAL> &h_vector,
    std::vector<std::complex<REAL>> &form_factors,
    const std::vector<bool> &include_atom,
    const std::vector<int> &type_2_wfn_type,
    const std::vector<std::vector<REAL>> &def_val_slater_normalization,
    const std::vector<int> &typeMaxL) {
    //--------

    mSphericalHarmonicsData.resize(1);
    mSphericalHarmonicsData[0].resize(5);
    for (int i = 0; i < 5; i++) mSphericalHarmonicsData[0][i].resize(2 * i + 1);

    //--------

    REAL hVectorLength;

    // hVectorLength2 = h_vector * h_vector;
    hVectorLength = sqrt(h_vector * h_vector);
    Vector3<REAL> normalized_h = h_vector / hVectorLength;

    int atomWfnIdx, atomTypeIdx;
    complex<REAL> atom_f_def_val, aux;

    //--

    vector<vector<REAL>> g_functions_and_slater_norm(type_parameters.size(),
                                                     vector<REAL>(5));

    int nAtoms;
    nAtoms = atom_to_type_map.size();

    form_factors.resize(nAtoms);

    //

    bool hkl000 = (hVectorLength < 1e-10);

    int nTypes = type_parameters.size();

    for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        int wfnTypeIdx = type_2_wfn_type[typeIdx];

        int nL = wfn_parameters[wfnTypeIdx].def_valence_pow.size();

        const vector<int> &def_valence_pow =
            wfn_parameters[wfnTypeIdx].def_valence_pow;

        if (nL > 0)
            g_functions_and_slater_norm[typeIdx][0] =
                def_val_slater_normalization[wfnTypeIdx][0] *
                sto_scattering::gFunction<0>(
                    int(def_valence_pow[0]) + 2,
                    hVectorLength / type_parameters[typeIdx].kappa_def_valence,
                    wfn_parameters[wfnTypeIdx].def_valence_exp);
        if (nL > 1)
            g_functions_and_slater_norm[typeIdx][1] =
                def_val_slater_normalization[wfnTypeIdx][1] *
                sto_scattering::gFunction<1>(
                    int(def_valence_pow[1]) + 2,
                    hVectorLength / type_parameters[typeIdx].kappa_def_valence,
                    wfn_parameters[wfnTypeIdx].def_valence_exp);

        if (nL > 2)
            g_functions_and_slater_norm[typeIdx][2] =
                def_val_slater_normalization[wfnTypeIdx][2] *
                sto_scattering::gFunction<2>(
                    int(def_valence_pow[2]) + 2,
                    hVectorLength / type_parameters[typeIdx].kappa_def_valence,
                    wfn_parameters[wfnTypeIdx].def_valence_exp);
        if (nL > 3)
            g_functions_and_slater_norm[typeIdx][3] =
                def_val_slater_normalization[wfnTypeIdx][3] *
                sto_scattering::gFunction<3>(
                    int(def_valence_pow[3]) + 2,
                    hVectorLength / type_parameters[typeIdx].kappa_def_valence,
                    wfn_parameters[wfnTypeIdx].def_valence_exp);
        if (nL > 4)
            g_functions_and_slater_norm[typeIdx][4] =
                def_val_slater_normalization[wfnTypeIdx][4] *
                sto_scattering::gFunction<4>(
                    int(def_valence_pow[4]) + 2,
                    hVectorLength / type_parameters[typeIdx].kappa_def_valence,
                    wfn_parameters[wfnTypeIdx].def_valence_exp);
    }

    //------------- end of pre_atom_loop_sf_calc

    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
        if (!include_atom[atomIdx]) {
            form_factors[atomIdx] = 0;
            continue;
        }

        atomWfnIdx = atom_to_wfn_map[atomIdx];
        atomTypeIdx = atom_to_type_map[atomIdx];

        if (hkl000)
            atom_f_def_val = 0;
        else
            atom_f_def_val = calculateDeformationValence(
                type_parameters[atomTypeIdx].p_lm,
                g_functions_and_slater_norm[atomTypeIdx],
                local_coordinate_systems[atomIdx],
                normalized_h,
                typeMaxL[atomTypeIdx],
                mSphericalHarmonicsData[0]);

        form_factors[atomIdx] = atom_f_def_val + f_spherical[atomTypeIdx];
    }
}

void HansenCoppens_SF_Engine4::calculateSphericalTermsInFormFactors(
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<double> h, std::vector<std::vector<REAL>> &f_core,
    std::vector<std::vector<REAL>> &f_sph_valence,
    const std::vector<int> &type_2_wfn_type,
    const std::vector<std::vector<REAL>> &def_val_slater_normalization,
    const std::vector<int> &typeMaxL) {
    //--
    int nTypes, nWfnTypes;
    nTypes = type_parameters.size();
    nWfnTypes = wfn_parameters.size();
    vector<REAL> wfn_spherical_core_sf(nWfnTypes);
    vector<REAL> wfn_spherical_valence_sf(nTypes);
    vector<vector<REAL>> g_functions_and_slater_norm(nTypes, vector<REAL>(5));

    int nH = h.size();

    f_core.resize(nWfnTypes, vector<double>(nH));
    f_sph_valence.resize(nTypes, vector<double>(nH));

    for (int hIndex = 0; hIndex < nH; hIndex++) {
        for (int wfnTypeIdx = 0; wfnTypeIdx < nWfnTypes; wfnTypeIdx++)
            // wfn_spherical_core_sf[wfnTypeIdx] =
            f_core[wfnTypeIdx][hIndex] =
                sto_scattering::scatteringSphericalDensity(
                    wfn_parameters[wfnTypeIdx].core_coeff,
                    wfn_parameters[wfnTypeIdx].core_exp,
                    wfn_parameters[wfnTypeIdx].core_pow,
                    h[hIndex]);

        for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
            int wfnTypeIdx = type_2_wfn_type[typeIdx];

            // wfn_spherical_valence_sf[typeIdx] =
            f_sph_valence[typeIdx][hIndex] =
                sto_scattering::scatteringSphericalDensity(
                    wfn_parameters[wfnTypeIdx].valence_coeff,
                    wfn_parameters[wfnTypeIdx].valence_exp,
                    wfn_parameters[wfnTypeIdx].valence_pow,
                    h[hIndex] / type_parameters[typeIdx].kappa_spherical);

            // wfn_spherical_valence_sf[typeIdx] *=
            // type_parameters[typeIdx].p_val;
            f_sph_valence[typeIdx][hIndex] *= type_parameters[typeIdx].p_val;
        }
    }
}

void HansenCoppens_SF_Engine4::calculateGlobalCoordinatesPlm(
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Matrix3<REAL>>
        &local_coordinate_systems,  // rows are vectors
    std::vector<std::vector<std::vector<double>>> &atomPlms) {
    int maxL = 4;
    vector<vector<double>> den2wfn;
    real_spherical_harmonics::getDensityToWfnMultipliers(maxL, den2wfn);

    int nTypes = type_parameters.size();
    vector<vector<vector<double>>> typePlmWfn(nTypes);
    for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        int typeMaxL = type_parameters[typeIdx].p_lm.size() - 1;
        typePlmWfn[typeIdx] = type_parameters[typeIdx].p_lm;
        for (int l = 0; l <= typeMaxL; l++)
            for (int i = 0; i < 2 * l + 1; i++) {
                int abs_m = abs(l - i);
                typePlmWfn[typeIdx][l][i] *= den2wfn[l][abs_m];
            }
    }

    SphConverter sphConverter;
    vector<vector<vector<double>>> conversionMatrices;
    sphConverter.setMaxL(maxL);

    vector<vector<double>> localCoordinates(3, vector<double>(3));
    vector<vector<double>> cartesianCoordinates{
        {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

    int nAtoms = atom_to_type_map.size();
    atomPlms.resize(nAtoms);

    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                localCoordinates[i][j] =
                    local_coordinate_systems[atomIdx](i, j);

        sphConverter.convert(
            localCoordinates, cartesianCoordinates, conversionMatrices);
        int atomType = atom_to_type_map[atomIdx];
        int typeMaxL = type_parameters[atomType].p_lm.size() - 1;
        atomPlms[atomIdx].resize(typeMaxL + 1);
        for (int l = 0; l <= typeMaxL; l++) {
            atomPlms[atomIdx][l].resize(2 * l + 1);
            for (int i = 0; i < 2 * l + 1; i++) {
                atomPlms[atomIdx][l][i] = 0.0;
                for (int j = 0; j < 2 * l + 1; j++)
                    atomPlms[atomIdx][l][i] += conversionMatrices[l][i][j] *
                                               typePlmWfn[atomType][l][j];

                int abs_m = abs(l - i);
                atomPlms[atomIdx][l][i] /= den2wfn[l][abs_m];
            }
        }
    }
}
/*
 * void calculateSF(
 *    const UnitCell &unitCell,
 *    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
 *    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
 *    const std::vector<int> &atom_to_wfn_map,
 *    const std::vector<int> &atom_to_type_map,
 *    const std::vector<Vector3<REAL> > &atomicPositions,
 *    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
 *    const std::vector<REAL> &atomic_occupancy,
 *    const std::vector<REAL> &atomic_multiplicity_factor,
 *    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
 *    const std::vector<sf_engine_data_types::SymmetryOperation>
 * &symmetry_operations, bool centrosymmetric, const Vector3<REAL>
 * &inversionTranslation, const std::vector<Vector3<REAL> > &h_vectors, const
 * std::vector<Vector3i >& hkl_indices, std::vector<std::complex<REAL> > &f,
 *    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
 *    const std::vector<std::complex<REAL> > &dTarget_df,
 *    const std::vector<bool> &include_atom_contribution,
 *    int nThreads);
 */

inline void printStep(std::string name) {
    std::cout << "calculateSF: after - " << name << " - time = " << std::clock()
              << std::endl;
}

inline void printInLoop(std::string name) {
    // std::cout << "calculateSF-main: after - " << name << " - time = " <<
    // std::clock() << std::endl;
}

template <typename T>
inline void infrequentValueLog(std::string name, T value) {
    // std::cout << name << " = " << value << std::endl;
}

template <typename T>
inline void valueLog(std::string name, T value) {
    // std::cout << name << " = " << value << std::endl;
}

constexpr std::array<int, 3> binSize = {8, 8, 8};
// constexpr bool virtLinePhaseFlag = true; // maybe add later (the ability to
// change to false) constexpr bool virtLineTemperatureFlag = true; // maybe add
// later (the ability to change to false)
constexpr bool virtHklPhaseFlag = false;
constexpr bool virtHklTemperatureFlag = false;
constexpr bool virtHklTMulFlag = false;
// constexpr bool symOpExtractionFlag = true; // maybe add later (the ability to
// change to false) constexpr bool fSymDeduplicationFlag = false; // cctbx
// already dedupes // maybe add later (the ability to change to true) constexpr
// bool versorDeduplicationFlag = true; // maybe add later (the ability to
// change to true) since it should only speed it up ~1.2 times constexpr bool
// lengthDeduplicationFlag = true; // maybe add later (the ability to change to
// true) constexpr bool symOpOffsetDeduplication = false; // maybe add later
// (the ability to change to true)

constexpr double two_pi = 2.0 * M_PI;
constexpr double two_pi_squared = two_pi * M_PI;
constexpr double four_pi_squared = two_pi * two_pi;
constexpr double four_pi = 4.0 * M_PI;

inline REAL square(REAL x) { return x * x; }

inline bool closeToZero(REAL x) { return -1e-12 < x and x < 1e-12; }

inline bool closeToZero(Vector3<REAL> vec) {
    return (closeToZero(vec[0]) and closeToZero(vec[1]) and
            closeToZero(vec[2]));
}

inline bool closeToZero(Matrix3<REAL> mat) {
    return (closeToZero({mat(0, 0), mat(0, 1), mat(0, 2)}) and
            closeToZero({mat(1, 0), mat(1, 1), mat(1, 2)}) and
            closeToZero({mat(2, 0), mat(2, 1), mat(2, 2)}));
}

inline Matrix3<REAL> U(const std::vector<REAL> &adps) {
    assert(adps.size() == 6);
    return Matrix3<REAL>(adps[0],
                         adps[3],
                         adps[4],
                         adps[3],
                         adps[1],
                         adps[5],
                         adps[4],
                         adps[5],
                         adps[2]);
}

inline int gcd(int a, int b) {
    while (b != 0) {
        int t = b;
        b = a % b;
        a = t;
    }
    return a;
}

inline REAL sqrt(REAL x) {
    return std::sqrt(x);  // TODO make it use a more precise approach (some new
                          // c++ versions may not have support for double here)
}

inline REAL pow(REAL x, int n) {
    return std::pow(x,
                    n);  // TODO make it use a more precise approach (some new
                         // c++ versions may not have support for double here)
}
template <typename T>
inline void save_to_file(const std::string &filename, const T &value) {
    std::ofstream os(filename, std::ios::binary);
    if (!os) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    cereal::BinaryOutputArchive ar(os);
    ar(value);
}
/* this breaks the results by ~10^-5 (tyrosine) TODO make a more robust
replacement (for now use the compiler default)
 * REAL cos(REAL x){
 *    return std::cosf(std::fmodf(x, two_pi));
}

REAL sin(REAL x){
return std::sinf(std::fmodf(x, two_pi));
}
*/

#include <algorithm>
#include <complex>
#include <iomanip>
#include <iostream>

// -----------------------------------------------------------------------------
// Helper print utilities
// -----------------------------------------------------------------------------

template <typename T>
void printVector3(const Vector3<T> &v, const std::string &name) {
    std::cout << name << " = (" << v.x << ", " << v.y << ", " << v.z << ")\n";
}

template <typename T>
void printMatrix3(const Matrix3<T> &m, const std::string &name) {
    std::cout << name << " =\n";

    for (int r = 0; r < 3; ++r) {
        std::cout << "    [ ";
        for (int c = 0; c < 3; ++c) {
            std::cout << std::setw(12) << m(r, c) << " ";
        }
        std::cout << "]\n";
    }
}

template <typename T>
void printStdVector(const std::vector<T> &v, const std::string &name,
                    size_t start = 0, size_t count = 5) {
    std::cout << name << " (size = " << v.size() << ", range = [" << start
              << ":" << (start + count) << "]) : ";

    if (start >= v.size()) {
        std::cout << "<out of range>\n";
        return;
    }

    size_t end = std::min(start + count, v.size());

    for (size_t i = start; i < end; ++i) {
        std::cout << v[i];

        if (i + 1 < end) std::cout << ", ";
    }

    if (end < v.size()) std::cout << " ...";

    std::cout << "\n";
}
template <typename T>
void printVector3Array(const std::vector<Vector3<T>> &v,
                       const std::string &name, size_t start = 0,
                       size_t count = 5) {
    std::cout << name << " (size = " << v.size() << ")\n";

    if (start >= v.size()) {
        std::cout << "  <out of range>\n";
        return;
    }

    size_t end = std::min(start + count, v.size());

    for (size_t i = start; i < end; ++i) {
        std::cout << "  [" << i << "] = (" << v[i].x << ", " << v[i].y << ", "
                  << v[i].z << ")\n";
    }
}

template <typename T>
void printComplex(const std::complex<T> &c, const std::string &name) {
    std::cout << name << " = " << c.real() << " + " << c.imag() << "i\n";
}

void HansenCoppens_SF_Engine4::calculateSF(
    const UnitCell &unitCell,
    const std::vector<sf_engine_data_types::HC_WfnParam> &wfnParams,  // per wfn
    const std::vector<sf_engine_data_types::HC_TypeParam>
        &typeParams,  // per type
    const std::vector<int> &atom_to_wfn_map,
    const std::vector<int> &atom_to_type_map,
    const std::vector<Vector3<REAL>> &atomicPositions,  // per atom
    const std::vector<std::vector<REAL>>
        &atomic_displacement_parameters,  // per atom and already premultiplied
                                          // by two_pi_squared
    const std::vector<REAL> &atomic_occupancy,
    const std::vector<std::complex<REAL>>
        &anomalous_dispersion,  // per atom // TODO use
    const std::vector<REAL> &atomic_multiplicity_factor,
    const std::vector<Matrix3<REAL>> &local_coordinate_systems,  // per atom
    const std::vector<sf_engine_data_types::SymmetryOperation> &symOps,
    bool centrosymmetric,                       // TODO use
    const Vector3<REAL> &inversionTranslation,  // TODO use
    const std::vector<Vector3<REAL>> &hVectors,
    const std::vector<Vector3i> &hkl_indices,  // unused
    std::vector<std::complex<REAL>> &f,
    std::vector<TargetFunctionAtomicParamDerivatives>
        &dTarget_dparam,  // per atom // TODO generate
    const std::vector<std::complex<REAL>> &dTarget_df,   // per hkl // TODO use
    const std::vector<bool> &include_atom_contribution,  // per atom
    int nThreads,
    const DerivativesSelector &derivativesSwitch,  // TODO use
    bool electron,                                 // TODO use
    const std::vector<int> &atomic_numbers)        // TODO use
{
    std::cout << "\n";
    std::cout << "=====================================================\n";
    std::cout << "                 INPUT DEBUG DUMP\n";
    std::cout << "=====================================================\n\n";

    // -------------------------------------------------------------------------
    // Unit cell
    // -------------------------------------------------------------------------

    std::cout << "---------------- UNIT CELL ----------------\n";

    std::cout << "a      = " << unitCell.a() << "\n"
              << "b      = " << unitCell.b() << "\n"
              << "c      = " << unitCell.c() << "\n"
              << "alpha  = " << unitCell.alpha() << "\n"
              << "beta   = " << unitCell.beta() << "\n"
              << "gamma  = " << unitCell.gamma() << "\n\n";

    printMatrix3(unitCell.getFractionalToCartesianMatrix(),
                 "Fractional -> Cartesian");

    std::cout << "\n";

    printMatrix3(unitCell.getCartesianToFractionalMatrix(),
                 "Cartesian -> Fractional");

    // -------------------------------------------------------------------------
    // Sizes
    // -------------------------------------------------------------------------

    std::cout << "\n---------------- CONTAINER SIZES ----------------\n";

    std::cout << "wfnParams.size()                    = " << wfnParams.size()
              << "\n"
              << "typeParams.size()                   = " << typeParams.size()
              << "\n"
              << "atomicPositions.size()              = "
              << atomicPositions.size() << "\n"
              << "symOps.size()                       = " << symOps.size()
              << "\n"
              << "hVectors.size()                     = " << hVectors.size()
              << "\n"
              << "f.size()                            = " << f.size() << "\n"
              << "dTarget_dparam.size()               = "
              << dTarget_dparam.size() << "\n";

    // -------------------------------------------------------------------------
    // Wavefunction params
    // -------------------------------------------------------------------------

    if (!wfnParams.empty()) {
        std::cout << "\n---------------- FIRST WFN PARAM ----------------\n";

        const auto &w = wfnParams.front();

        std::cout << "label = " << w.label << "\n";

        // ---- multiple ranges ----

        printStdVector(w.core_coeff, "core_coeff [0:5]", 0, 5);
        printStdVector(w.core_coeff, "core_coeff [10:15]", 10, 5);

        printStdVector(w.core_exp, "core_exp [0:5]", 0, 5);
        printStdVector(w.core_exp, "core_exp [10:15]", 10, 5);

        printStdVector(w.core_pow, "core_pow [0:5]", 0, 5);

        printStdVector(w.valence_coeff, "valence_coeff [0:5]", 0, 5);
        printStdVector(w.valence_coeff, "valence_coeff [10:15]", 10, 5);

        printStdVector(w.valence_exp, "valence_exp [0:5]", 0, 5);
        printStdVector(w.valence_pow, "valence_pow [10:15]", 10, 5);

        std::cout << "def_valence_exp = " << w.def_valence_exp << "\n";

        printStdVector(w.def_valence_pow, "def_valence_pow [0:5]", 0, 5);

        printComplex(w.anomalous_scattering, "anomalous_scattering");
    }

    // -------------------------------------------------------------------------
    // Type params
    // -------------------------------------------------------------------------

    if (!typeParams.empty()) {
        std::cout << "\n---------------- FIRST TYPE PARAM ----------------\n";

        const auto &t = typeParams.front();

        std::cout << "p_val               = " << t.p_val << "\n"
                  << "kappa_def_valence   = " << t.kappa_def_valence << "\n"
                  << "kappa_spherical     = " << t.kappa_spherical << "\n";

        std::cout << "p_lm rows           = " << t.p_lm.size() << "\n";

        if (!t.p_lm.empty()) {
            printStdVector(t.p_lm[0], "p_lm[0] [0:5]", 0, 5);

            printStdVector(t.p_lm[0], "p_lm[0] [10:15]", 10, 5);
        }
    }

    // -------------------------------------------------------------------------
    // Atom data
    // -------------------------------------------------------------------------

    std::cout << "\n---------------- ATOM DATA ----------------\n";

    // print different atom ranges
    std::vector<size_t> atomStarts = {0, 10, 100};

    for (size_t start : atomStarts) {
        if (start >= atomicPositions.size()) continue;

        size_t end = std::min(start + 3, atomicPositions.size());

        std::cout << "\n========== ATOMS " << start << " -> " << (end - 1)
                  << " ==========\n";

        for (size_t i = start; i < end; ++i) {
            std::cout << "\nAtom #" << i << "\n";

            printVector3(atomicPositions[i], "position");

            if (i < atomic_numbers.size())
                std::cout << "atomic number = " << atomic_numbers[i] << "\n";

            if (i < atom_to_wfn_map.size())
                std::cout << "wfn map       = " << atom_to_wfn_map[i] << "\n";

            if (i < atom_to_type_map.size())
                std::cout << "type map      = " << atom_to_type_map[i] << "\n";

            if (i < atomic_occupancy.size())
                std::cout << "occupancy     = " << atomic_occupancy[i] << "\n";

            if (i < atomic_multiplicity_factor.size())
                std::cout << "multiplicity  = " << atomic_multiplicity_factor[i]
                          << "\n";

            if (i < anomalous_dispersion.size()) {
                printComplex(anomalous_dispersion[i], "anomalous_dispersion");
            }

            if (i < atomic_displacement_parameters.size()) {
                printStdVector(
                    atomic_displacement_parameters[i], "ADP [0:6]", 0, 6);

                printStdVector(
                    atomic_displacement_parameters[i], "ADP [6:12]", 6, 6);
            }

            if (i < local_coordinate_systems.size()) {
                printMatrix3(local_coordinate_systems[i],
                             "local_coordinate_system");
            }

            if (i < include_atom_contribution.size()) {
                std::cout << "include contribution = " << std::boolalpha
                          << include_atom_contribution[i] << "\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Symmetry operations
    // -------------------------------------------------------------------------

    if (!symOps.empty()) {
        std::cout << "\n---------------- SYMMETRY OPS ----------------\n";

        size_t n = std::min<size_t>(symOps.size(), 3);

        for (size_t i = 0; i < n; ++i) {
            std::cout << "\nSymmetry op #" << i << "\n";

            printMatrix3(symOps[i].rotation, "rotation");

            printVector3(symOps[i].translation, "translation");
        }
    }

    // -------------------------------------------------------------------------
    // h vectors
    // -------------------------------------------------------------------------

    if (!hVectors.empty()) {
        std::cout << "\n---------------- H VECTORS ----------------\n";

        printVector3Array(hVectors, "hVectors [0:5]", 0, 5);

        printVector3Array(hVectors, "hVectors [100:105]", 100, 5);

        printVector3Array(hVectors, "hVectors [1000:1005]", 1000, 5);
    }

    // -------------------------------------------------------------------------
    // HKL
    // -------------------------------------------------------------------------

    if (!hkl_indices.empty()) {
        std::cout << "\n---------------- HKL INDICES ----------------\n";

        std::vector<size_t> starts = {0, 100, 1000};

        for (size_t start : starts) {
            if (start >= hkl_indices.size()) continue;

            size_t end = std::min(start + 5, hkl_indices.size());

            std::cout << "\nHKL range [" << start << ":" << end << "]\n";

            for (size_t i = start; i < end; ++i) {
                std::cout << "HKL[" << i << "] = (" << hkl_indices[i].x << ", "
                          << hkl_indices[i].y << ", " << hkl_indices[i].z
                          << ")\n";
            }
        }
    }

    // -------------------------------------------------------------------------
    // Structure factors
    // -------------------------------------------------------------------------

    if (!f.empty()) {
        std::cout << "\n---------------- STRUCTURE FACTORS ----------------\n";

        std::vector<size_t> starts = {0, 100, 1000};

        for (size_t start : starts) {
            if (start >= f.size()) continue;

            size_t end = std::min(start + 5, f.size());

            std::cout << "\nf range [" << start << ":" << end << "]\n";

            for (size_t i = start; i < end; ++i) {
                printComplex(f[i], "f[" + std::to_string(i) + "]");
            }
        }
    }

    // -------------------------------------------------------------------------
    // Derivatives
    // -------------------------------------------------------------------------

    std::cout << "\n---------------- DERIVATIVES SETTINGS ----------------\n";

    std::cout << "d_xyz  = " << derivativesSwitch.d_xyz << "\n"
              << "d_adp  = " << derivativesSwitch.d_adp << "\n"
              << "d_occ  = " << derivativesSwitch.d_occ << "\n"
              << "d_anom = " << derivativesSwitch.d_anom << "\n";

    // -------------------------------------------------------------------------
    // Misc
    // -------------------------------------------------------------------------

    std::cout << "\n---------------- MISC ----------------\n";

    std::cout << "centrosymmetric = " << std::boolalpha << centrosymmetric
              << "\n";

    printVector3(inversionTranslation, "inversionTranslation");

    std::cout << "electron = " << electron << "\n";

    std::cout << "nThreads = " << nThreads << "\n";

    std::cout << "\n=====================================================\n";
    std::cout << "                  END DEBUG DUMP\n";
    std::cout << "=====================================================\n";

#ifndef TEST

    save_to_file("unit_cell.bin", unitCell);
    save_to_file("wfn_parameters.bin", wfnParams);
    save_to_file("type_parameters.bin", typeParams);
    save_to_file("atom_to_wfn_map.bin", atom_to_wfn_map);
    save_to_file("atom_to_type_map.bin", atom_to_type_map);
    save_to_file("atomic_positions.bin", atomicPositions);
    save_to_file("atomic_displacement_parameters.bin",
                 atomic_displacement_parameters);
    save_to_file("atomic_occupancy.bin", atomic_occupancy);
    save_to_file("anomalous_dispersion.bin", anomalous_dispersion);
    save_to_file("atomic_multiplicity_factor.bin", atomic_multiplicity_factor);
    save_to_file("local_coordinate_systems.bin", local_coordinate_systems);
    save_to_file("symmetry_operations.bin", symOps);
    save_to_file("centrosymmetric.bin", centrosymmetric);
    save_to_file("inversion_translation.bin", inversionTranslation);
    save_to_file("h_vectors.bin", hVectors);
    save_to_file("hkl_indices.bin", hkl_indices);
    save_to_file("f.bin", f);
    save_to_file("dtarget_dparam.bin", dTarget_dparam);
    save_to_file("dtarget_df.bin", dTarget_df);
    save_to_file("include_atom_contribution.bin", include_atom_contribution);
    save_to_file("n_threads.bin", nThreads);
    save_to_file("derivatives_switch.bin", derivativesSwitch);
    save_to_file("electron.bin", electron);
    save_to_file("atomic_number.bin", atomic_numbers);
#endif

    printStep("calculateSF start");
    const int trueNAtoms = atom_to_wfn_map.size();
    std::vector<int> usedAtomIndices;
    usedAtomIndices.clear();
    for (int i = 0; i < trueNAtoms; i++) {
        if (include_atom_contribution[i]) usedAtomIndices.emplace_back(i);
    }
    const int hklCount = hVectors.size();
    const int nAtoms = usedAtomIndices.size();

    infrequentValueLog("nAtoms", nAtoms);
    for (int atom = 0; atom < nAtoms; atom++) {
        valueLog("usedAtomIndices[atom]", usedAtomIndices[atom]);
        for (int i = 0;
             i < atomic_displacement_parameters[usedAtomIndices[atom]].size();
             i++) {
            valueLog("atomic_displacement_parameters[usedAtomIndices[atom]]",
                     atomic_displacement_parameters[usedAtomIndices[atom]][i]);
        }
    }

    const int nSymOps = symOps.size();

    std::vector<Matrix3<REAL>> symOpMults;
    std::vector<int> symOpToMult;
    symOpToMult.resize(nSymOps);
    symOpMults.emplace_back(symOps[0].rotation);
    valueLog("symOpIdx", 0);
    valueLog("current(0,0)", symOps[0].rotation(0, 0));
    valueLog("current(0,1)", symOps[0].rotation(0, 1));
    valueLog("current(0,2)", symOps[0].rotation(0, 2));
    valueLog("current(1,0)", symOps[0].rotation(1, 0));
    valueLog("current(1,1)", symOps[0].rotation(1, 1));
    valueLog("current(1,2)", symOps[0].rotation(1, 2));
    valueLog("current(2,0)", symOps[0].rotation(2, 0));
    valueLog("current(2,1)", symOps[0].rotation(2, 1));
    valueLog("current(2,2)", symOps[0].rotation(2, 2));
    for (int symOpIdx = 1; symOpIdx < nSymOps; symOpIdx++) {
        auto &current = symOps[symOpIdx].rotation;
        valueLog("symOpIdx", symOpIdx);
        valueLog("current(0,0)", current(0, 0));
        valueLog("current(0,1)", current(0, 1));
        valueLog("current(0,2)", current(0, 2));
        valueLog("current(1,0)", current(1, 0));
        valueLog("current(1,1)", current(1, 1));
        valueLog("current(1,2)", current(1, 2));
        valueLog("current(2,0)", current(2, 0));
        valueLog("current(2,1)", current(2, 1));
        valueLog("current(2,2)", current(2, 2));
        bool isInMults = false;
        for (int multIdx = 0; multIdx < symOpMults.size(); multIdx++) {
            if (closeToZero(current - symOpMults[multIdx])) {
                isInMults = true;
                symOpToMult[symOpIdx] = multIdx;
                break;
            }
        }
        if (not isInMults) {
            symOpToMult[symOpIdx] = symOpMults.size();
            symOpMults.emplace_back(current);
        }
    }
    const int symOpMultCount = symOpMults.size();
    infrequentValueLog("symOpMultCount", symOpMultCount);

    std::vector<Vector3<REAL>> symOpOffsets;
    std::vector<int> symOpToOffset;
    symOpToOffset.resize(nSymOps);
    symOpOffsets.emplace_back(symOps[0].translation);
    valueLog("symOpIdx", 0);
    valueLog("current[0]", symOps[0].translation[0]);
    valueLog("current[1]", symOps[0].translation[1]);
    valueLog("current[2]", symOps[0].translation[2]);
    for (int symOpIdx = 1; symOpIdx < nSymOps; symOpIdx++) {
        auto &current = symOps[symOpIdx].translation;
        valueLog("symOpIdx", symOpIdx);
        valueLog("current[0]", current[0]);
        valueLog("current[1]", current[1]);
        valueLog("current[2]", current[2]);
        bool isInOffsets = false;
        for (int offsetIdx = 0; offsetIdx < symOpOffsets.size(); offsetIdx++) {
            if (closeToZero(current - symOpOffsets[offsetIdx])) {
                isInOffsets = true;
                symOpToOffset[symOpIdx] = offsetIdx;
                break;
            }
        }
        if (not isInOffsets) {
            symOpToOffset[symOpIdx] = symOpOffsets.size();
            symOpOffsets.emplace_back(current);
        }
    }
    const int symOpOffsetCount = symOpOffsets.size();
    infrequentValueLog("symOpOffsetCount", symOpOffsetCount);

    printStep("symOp deduplication");

    Vector3<int> minHkl = hkl_indices[0];
    Vector3<int> maxHkl = hkl_indices[0];
    Vector3<int> step = {0, 0, 0};
    // Vector3<REAL> minH = hVectors[0];
    // Vector3<REAL> maxH = hVectors[0];
    for (int hklIdx = 1; hklIdx < hklCount; hklIdx++) {
        const Vector3<int> &currentHkl = hkl_indices[hklIdx];
        const Vector3<REAL> &currentH = hVectors[hklIdx];
        for (int i = 0; i < 3; i++) {
            if (currentHkl[i] < minHkl[i]) {
                minHkl[i] = currentHkl[i];
                // minH[i] = currentH[i];
            }
            if (currentHkl[i] > maxHkl[i]) {
                maxHkl[i] = currentHkl[i];
                // maxH[i] = currentH[i];
            }
            const int diff = std::abs(currentHkl[i] - hkl_indices[0][i]);
            if (step[i] == 0)
                step[i] = diff;
            else if ((diff % step[i]) != 0) {
                step[i] = gcd(diff, step[i]);
            };
        }
    }

    for (int i = 0; i < 3; i++) {
        if (step[i] == 0) step[i] = 1;
    }

    infrequentValueLog("step[0]", step[0]);
    infrequentValueLog("step[1]", step[1]);
    infrequentValueLog("step[2]", step[2]);

    /*
     *    Vector3<REAL> scale;
     *    for (int i =0; i<3; i++){
     *        if (((maxHkl[i] - minHkl[i]) == 0) or closeToZero(maxH[i] -
  minH[i])){
     *          if ((maxHkl[i] == 0) or closeToZero(maxH[i])){
     *            scale[i]=1.0;
  } else {
    scale[i]=maxH[i]/maxHkl[i];
  }
  }
  else{
    scale[i]=(maxH[i] - minH[i]) / (maxHkl[i] - minHkl[i]);
  }
  }

  infrequentValueLog("scale[0]", scale[0]);
  infrequentValueLog("scale[1]", scale[1]);
  infrequentValueLog("scale[2]", scale[2]);
  */

    Vector3<int> nHkls;
    for (int i = 0; i < 3; i++) {
        const int diff = maxHkl[i] - minHkl[i];
        assert(diff % step[i] == 0);
        nHkls[i] = diff / step[i] + 1;
    }

    Vector3<int> nBins;
    for (int i = 0; i < 3; i++) {
        if ((nHkls[i] % binSize[i]) == 0)
            nBins[i] = nHkls[i] / binSize[i];
        else
            nBins[i] = (nHkls[i] / binSize[i]) + 1;
    }

    const int totalNBins = nBins[0] * nBins[1] * nBins[2];

    printStep("hkl binning");

    // [h*nk*nl + k*nl + l]
    std::vector<bool> isBinUsed;
    isBinUsed.resize(totalNBins, false);
    for (int hklIdx = 0; hklIdx < hklCount; hklIdx++) {
        const Vector3<int> &currentHkl = hkl_indices[hklIdx];
        Vector3<int> binIdx;
        for (int i = 0; i < 3; i++) {
            // assert((currentHkl[i] % step[i]) == 0)
            binIdx[i] = (currentHkl[i] - minHkl[i]) / (step[i] * binSize[i]);
        }
        isBinUsed[(binIdx[0] * nBins[1] + binIdx[1]) * nBins[2] + binIdx[2]] =
            true;
    }
    std::vector<int> allBinsMap;
    allBinsMap.resize(nBins[0] * nBins[1] * nBins[2]);
    std::vector<Vector3<int>> usedBins;
    usedBins.clear();
    for (int h = 0; h < nBins[0]; h++) {
        for (int k = 0; k < nBins[1]; k++) {
            for (int l = 0; l < nBins[2]; l++) {
                if (isBinUsed[(h * nBins[1] + k) * nBins[2] + l]) {
                    allBinsMap[(h * nBins[1] + k) * nBins[2] + l] =
                        usedBins.size();
                    usedBins.emplace_back(h, k, l);
                }
            }
        }
    }

    const int binCount = usedBins.size();

    infrequentValueLog("binCount", binCount);

    printStep("bin deduplication");

    vector<vector<Vector3d>> r_atom_rot(nAtoms,
                                        vector<Vector3<REAL>>(symOpMultCount));
#pragma omp parallel for num_threads(nThreads) collapse(2)
    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
        for (int symOpMultIdx = 0; symOpMultIdx < symOpMultCount;
             symOpMultIdx++) {
            r_atom_rot[atomIdx][symOpMultIdx] =
                symOpMults[symOpMultIdx] *
                atomicPositions[usedAtomIndices[atomIdx]];
        }
    }

    ReciprocalLatticeUnitCell recUnitCell(unitCell);
    std::array<Vector3<REAL>, 3> stepCartesian;
    recUnitCell.fractionalToCartesian({(double)step[0], 0.0, 0.0},
                                      stepCartesian[0]);
    recUnitCell.fractionalToCartesian({0.0, (double)step[1], 0.0},
                                      stepCartesian[1]);
    recUnitCell.fractionalToCartesian({0.0, 0.0, (double)step[2]},
                                      stepCartesian[2]);

    std::vector<Vector3<REAL>> bin000Cartesian;
    bin000Cartesian.resize(binCount);
#pragma omp parallel for num_threads(nThreads)
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        Vector3<int> bin000;
        for (int i = 0; i < 3; i++) {
            bin000[i] = usedBins[binIdx][i] * step[i] * binSize[i] + minHkl[i];
            infrequentValueLog("bin000[i]", bin000[i]);
        }
        recUnitCell.fractionalToCartesian(bin000, bin000Cartesian[binIdx]);
        infrequentValueLog("bin000Cartesian[binIdx][0]",
                           bin000Cartesian[binIdx][0]);
        infrequentValueLog("bin000Cartesian[binIdx][1]",
                           bin000Cartesian[binIdx][1]);
        infrequentValueLog("bin000Cartesian[binIdx][2]",
                           bin000Cartesian[binIdx][2]);
    }

    printStep("fractional to cartesian conversion");

    // [atom][symOpMult | 0]
    std::vector<std::vector<Matrix3<REAL>>> rotatedUs;
    rotatedUs.resize(nAtoms);
    for (int atom = 0; atom < nAtoms; atom++) {
        rotatedUs[atom].resize(
            (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 6)
                ? symOpMultCount
                : 0);
    }
#pragma omp parallel for num_threads(nThreads)
    for (int atom = 0; atom < nAtoms; atom++) {
        if (atomic_displacement_parameters[usedAtomIndices[atom]].size() == 6) {
            /*auto ftcMatrixT = recUnitCell.getFractionalToCartesianMatrix();
             *            ftcMatrixT.transpose();
             *            const auto orig =
             * recUnitCell.getFractionalToCartesianMatrix() *
             * U(atomic_displacement_parameters[usedAtomIndices[atom]]) *
             * ftcMatrixT;
             */
            valueLog("atom", atom);
            for (int i = 0; i < 6; i++) {
                valueLog(
                    "atomic_displacement_parameters[usedAtomIndices[atom]][i]",
                    atomic_displacement_parameters[usedAtomIndices[atom]][i]);
            }
            const auto orig =
                U(atomic_displacement_parameters[usedAtomIndices[atom]]);
            for (int symOpIdx = 0; symOpIdx < symOpMultCount; symOpIdx++) {
                Matrix3<REAL> symOpT;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        symOpT(i, j) = symOpMults[symOpIdx](j, i);
                    }
                }
                rotatedUs[atom][symOpIdx] =
                    symOpT * orig * symOpMults[symOpIdx];
                /*
                 *                // trying to do M = symOp * orig * symOp^T
                 *                // because orig^T = orig and M^T = M
                 *                // it is enough to compute only 6/9 of acc and
            6/9 of M = P * symOp^T
                 *                  std::array<Matrix3<REAL>, 3> acc;
                 *                  for (int x=0; x<3; x++){
                 *                      for (int y=0; y<3; y++){
                 *                          for (int i=y; i<3; i++){
                 *                              acc[i](x, y) = symOp(i, y) *
            orig(x, i);
            }
            }
            }
            Matrix3<REAL> P;
            for (int x=0; x<3; x++){
                for (int y=0; y<3; y++){
                    for (int i=0; i<y; i++){
                        P(x, y) += acc[y](x, i);
            }
            for (int i=y; i<3; i++){
                P(x, y) += acc[i](x, y);
            }
            }
            }
            Matrix3<REAL> M;
            for (int x=0; x<3; x++){
                for (int y=0; y<=x; y++){
                    for (int i=0; i<3; i++){
                        M(x, y) += P(i, y) * symOp(i, x);
            }
            }
            }
            rotatedUs[atom][symOpIdx] = Matrix3<REAL>(
                M(0, 0), M(1, 0), M(2, 0),
                M(1, 0), M(1, 1), M(2, 1),
                M(2, 0), M(2, 1), M(2, 2));
                */
            }
        }
    }
    printStep("rotating Us");

    // [bin][atom][symOpMult | 1]
    std::vector<std::vector<std::vector<REAL>>> temperatureFactorRoots;
    temperatureFactorRoots.resize(binCount);
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        temperatureFactorRoots[binIdx].resize(nAtoms);
        for (int atom = 0; atom < nAtoms; atom++) {
            temperatureFactorRoots[binIdx][atom].resize(
                (atomic_displacement_parameters[usedAtomIndices[atom]].size() ==
                 1)
                    ? 1
                    : symOpMultCount);
        }
    }
#pragma omp parallel for num_threads(nThreads) collapse(2)
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        for (int atom = 0; atom < nAtoms; atom++) {
            bool iso =
                (atomic_displacement_parameters[usedAtomIndices[atom]].size() ==
                 1);
            if (iso) {
                REAL T1iso = std::exp(
                    -atomic_displacement_parameters[usedAtomIndices[atom]][0] *
                    (square(bin000Cartesian[binIdx][0]) +
                     square(bin000Cartesian[binIdx][1]) +
                     square(bin000Cartesian[binIdx][2])));
                temperatureFactorRoots[binIdx][atom][0] = T1iso;
            } else {
                for (int symOpIdx = 0; symOpIdx < symOpMultCount; symOpIdx++) {
                    Matrix3<REAL> currentU = rotatedUs[atom][symOpIdx];
                    REAL T1 = std::exp(-(
                        square(bin000Cartesian[binIdx][0]) * currentU(0, 0) +
                        square(bin000Cartesian[binIdx][1]) * currentU(1, 1) +
                        square(bin000Cartesian[binIdx][2]) * currentU(2, 2) +
                        2.0 *
                            (bin000Cartesian[binIdx][0] *
                                 bin000Cartesian[binIdx][1] * currentU(0, 1) +
                             bin000Cartesian[binIdx][0] *
                                 bin000Cartesian[binIdx][2] * currentU(0, 2) +
                             bin000Cartesian[binIdx][1] *
                                 bin000Cartesian[binIdx][2] * currentU(1, 2))));
                    temperatureFactorRoots[binIdx][atom][symOpIdx] = T1;
                }
            }
        }
    }

    printStep("temperature factor roots");

    // [bin][atom][symOp]
    std::vector<std::vector<std::vector<std::complex<REAL>>>> phaseFactorRoots;
    phaseFactorRoots.resize(binCount);
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        phaseFactorRoots[binIdx].resize(nAtoms);
        for (int atom = 0; atom < nAtoms; atom++) {
            phaseFactorRoots[binIdx][atom].resize(nSymOps);
        }
    }
#pragma omp parallel for num_threads(nThreads) collapse(3)
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        for (int atom = 0; atom < nAtoms; atom++) {
            for (int symOpIdx = 0; symOpIdx < nSymOps; symOpIdx++) {
                const REAL phase_angle_root =
                    two_pi *
                    (r_atom_rot[atom][symOpToMult[symOpIdx]] +
                     symOpOffsets[symOpToOffset[symOpIdx]]) *
                    bin000Cartesian[binIdx];
                const std::complex<REAL> result = {cos(phase_angle_root),
                                                   sin(phase_angle_root)};
                phaseFactorRoots[binIdx][atom][symOpIdx] = result;
            }
        }
    }

    printStep("phase factor roots");

    // [dir][atom][symOp][n]
    std::array<std::vector<std::vector<std::vector<std::complex<REAL>>>>, 3>
        phaseFactorMults;
    for (int i = 0; i < 3; i++) {
        phaseFactorMults[i].resize(nAtoms);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            phaseFactorMults[i][atomIdx].resize(nSymOps);
            for (int symOpIdx = 0; symOpIdx < nSymOps; symOpIdx++) {
                phaseFactorMults[i][atomIdx][symOpIdx].resize(binSize[i]);
            }
        }
    }
#pragma omp parallel for num_threads(nThreads) collapse(3)
    for (int i = 0; i < 3; i++) {
        for (int atom = 0; atom < nAtoms; atom++) {
            for (int symOpIdx = 0; symOpIdx < nSymOps; symOpIdx++) {
                const REAL phase_angle_mult =
                    two_pi *
                    (r_atom_rot[atom][symOpToMult[symOpIdx]] +
                     symOpOffsets[symOpToOffset[symOpIdx]]) *
                    stepCartesian[i];
                const std::complex<REAL> single = {cos(phase_angle_mult),
                                                   sin(phase_angle_mult)};
                std::complex<REAL> acc = 1.0;
                for (int j = 0; j < binSize[i]; j++) {
                    phaseFactorMults[i][atom][symOpIdx][j] = acc;
                    acc *= single;
                }
            }
        }
    }

    printStep("phase factor multipliers");

    // [h*nk*nl + k*nl + l][atom][symOp]
    std::vector<std::vector<std::vector<std::complex<REAL>>>> virtHklPhase;
    if constexpr (virtHklPhaseFlag) {
        virtHklPhase.resize(binSize[0] * binSize[1] * binSize[2]);
        for (int hklIdx = 0; hklIdx < (binSize[0] * binSize[1] * binSize[2]);
             hklIdx++) {
            virtHklPhase[hklIdx].resize(nAtoms);
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
                virtHklPhase[hklIdx][atomIdx].resize(nSymOps);
            }
        }
#pragma omp parallel for num_threads(nThreads) collapse(4)
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            for (int h = 0; h < binSize[0]; h++) {
                for (int k = 0; k < binSize[1]; k++) {
                    for (int l = 0; l < binSize[2]; l++) {
                        for (int symOpIdx = 0; symOpIdx < nSymOps; symOpIdx++) {
                            virtHklPhase[(h * binSize[1] + k) * binSize[2] +
                                         l][atomIdx][symOpIdx] =
                                phaseFactorMults[0][atomIdx][symOpIdx][h] *
                                phaseFactorMults[1][atomIdx][symOpIdx][k] *
                                phaseFactorMults[2][atomIdx][symOpIdx][l];
                        }
                    }
                }
            }
        }
    }

    printStep("virtual hkl phase factors");

    // [dir][atom][symOpMult | 1][n]
    std::array<std::vector<std::vector<std::vector<REAL>>>, 3>
        temperatureFactorMults;
    for (int i = 0; i < 3; i++) {
        temperatureFactorMults[i].resize(nAtoms);
        for (int atom = 0; atom < nAtoms; atom++) {
            temperatureFactorMults[i][atom].resize(
                (atomic_displacement_parameters[usedAtomIndices[atom]].size() ==
                 1)
                    ? 1
                    : symOpMultCount);
            for (int symOp = 0;
                 symOp < ((atomic_displacement_parameters[usedAtomIndices[atom]]
                               .size() == 1)
                              ? 1
                              : symOpMultCount);
                 symOp++) {
                temperatureFactorMults[i][atom][symOp].resize(binSize[i]);
            }
        }
    }
#pragma omp parallel for num_threads(nThreads) collapse(2)
    for (int i = 0; i < 3; i++) {
        for (int atom = 0; atom < nAtoms; atom++) {
            bool iso =
                (atomic_displacement_parameters[usedAtomIndices[atom]].size() ==
                 1);
            if (iso) {
                REAL c = std::exp(
                    -atomic_displacement_parameters[usedAtomIndices[atom]][0] *
                    (square(stepCartesian[i][0]) + square(stepCartesian[i][1]) +
                     square(stepCartesian[i][2])));
                REAL mult = 1.0;
                REAL acc = c;
                REAL c_pow = square(c);
                temperatureFactorMults[i][atom][0][0] = 1.0;
                for (int j = 1; j < binSize[i]; j++) {
                    // mult*=pow(c, j*2 - 1);
                    mult *= acc;
                    temperatureFactorMults[i][atom][0][j] = mult;
                    acc *= c_pow;
                }
            } else {
                for (int symOp = 0; symOp < symOpMultCount; symOp++) {
                    Matrix3<REAL> currentU = rotatedUs[atom][symOp];
                    REAL c = std::exp(
                        -(square(stepCartesian[i][0]) * currentU(0, 0) +
                          square(stepCartesian[i][1]) * currentU(1, 1) +
                          square(stepCartesian[i][2]) * currentU(2, 2) +
                          2.0 * (stepCartesian[i][0] * stepCartesian[i][1] *
                                     currentU(0, 1) +
                                 stepCartesian[i][0] * stepCartesian[i][2] *
                                     currentU(0, 2) +
                                 stepCartesian[i][1] * stepCartesian[i][2] *
                                     currentU(1, 2))));
                    REAL mult = 1.0;
                    REAL acc = c;
                    REAL c_pow = square(c);
                    temperatureFactorMults[i][atom][symOp][0] = 1.0;
                    for (int j = 1; j < binSize[i]; j++) {
                        // mult*=pow(c, j*2 - 1);
                        mult *= acc;
                        temperatureFactorMults[i][atom][symOp][j] = mult;
                        acc *= c_pow;
                    }
                }
            }
        }
    }

    printStep("temperature factor multipliers");

    // [dir][atom][symOpMult | 1][n][m]
    std::array<std::vector<std::vector<std::vector<std::vector<REAL>>>>, 3>
        temperatureFactorMultsSquare;
    for (int i = 0; i < 3; i++) {
        temperatureFactorMultsSquare[i].resize(nAtoms);
        for (int atom = 0; atom < nAtoms; atom++) {
            temperatureFactorMultsSquare[i][atom].resize(
                (atomic_displacement_parameters[usedAtomIndices[atom]].size() ==
                 1)
                    ? 1
                    : symOpMultCount);
            for (int symOp = 0;
                 symOp < ((atomic_displacement_parameters[usedAtomIndices[atom]]
                               .size() == 1)
                              ? 1
                              : symOpMultCount);
                 symOp++) {
                temperatureFactorMultsSquare[i][atom][symOp].resize(
                    binSize[(i + 1) % 3]);
                for (int n = 0; n < binSize[(i + 1) % 3]; n++) {
                    temperatureFactorMultsSquare[i][atom][symOp][n].resize(
                        binSize[(i + 2) % 3]);
                }
            }
        }
    }
#pragma omp parallel for num_threads(nThreads) collapse(2)
    for (int i = 0; i < 3; i++) {
        for (int atom = 0; atom < nAtoms; atom++) {
            bool iso =
                (atomic_displacement_parameters[usedAtomIndices[atom]].size() ==
                 1);
            int ns = (i + 1) % 3;
            int ms = (i + 2) % 3;
            if (iso) {
                REAL c = std::exp(
                    -2.0 *
                    atomic_displacement_parameters[usedAtomIndices[atom]][0] *
                    (stepCartesian[ns][0] * stepCartesian[ms][0] +
                     stepCartesian[ns][1] * stepCartesian[ms][1] +
                     stepCartesian[ns][2] * stepCartesian[ms][2]));
                REAL mult = 1.0;
                for (int n = 0; n < binSize[ns]; n++) {
                    REAL acc = 1.0;
                    for (int m = 0; m < binSize[ms]; m++) {
                        temperatureFactorMultsSquare[i][atom][0][n][m] = acc;
                        acc *= mult;
                    }
                    mult *= c;
                }
            } else {
                for (int symOp = 0; symOp < symOpMultCount; symOp++) {
                    Matrix3<REAL> currentU = rotatedUs[atom][symOp];
                    REAL c = std::exp(
                        -2.0 * (stepCartesian[ns][0] * stepCartesian[ms][0] *
                                    currentU(0, 0) +
                                stepCartesian[ns][1] * stepCartesian[ms][1] *
                                    currentU(1, 1) +
                                stepCartesian[ns][2] * stepCartesian[ms][2] *
                                    currentU(2, 2) +
                                (stepCartesian[ns][0] * stepCartesian[ms][1] +
                                 stepCartesian[ns][1] * stepCartesian[ms][0]) *
                                    currentU(0, 1) +
                                (stepCartesian[ns][0] * stepCartesian[ms][2] +
                                 stepCartesian[ns][2] * stepCartesian[ms][0]) *
                                    currentU(0, 2) +
                                (stepCartesian[ns][1] * stepCartesian[ms][2] +
                                 stepCartesian[ns][2] * stepCartesian[ms][1]) *
                                    currentU(1, 2)));
                    REAL mult = 1.0;
                    for (int n = 0; n < binSize[ns]; n++) {
                        REAL acc = 1.0;
                        for (int m = 0; m < binSize[ms]; m++) {
                            temperatureFactorMultsSquare[i][atom][symOp][n][m] =
                                acc;
                            acc *= mult;
                        }
                        mult *= c;
                    }
                }
            }
        }
    }

    printStep("temperature factor square multipliers");

    // [h*nk*nl + k*nl + l][atom][symOpMult | 1]
    std::vector<std::vector<std::vector<REAL>>> virtHklTemperature;
    if constexpr (virtHklTemperatureFlag) {
        virtHklTemperature.resize(binSize[0] * binSize[1] * binSize[2]);
        for (int hklIdx = 0; hklIdx < (binSize[0] * binSize[1] * binSize[2]);
             hklIdx++) {
            virtHklTemperature[hklIdx].resize(nAtoms);
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
                virtHklTemperature[hklIdx][atomIdx].resize(
                    (atomic_displacement_parameters[usedAtomIndices[atomIdx]]
                         .size() == 1)
                        ? 1
                        : symOpMultCount);
            }
        }
#pragma omp parallel for num_threads(nThreads) collapse(4)
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            for (int h = 0; h < binSize[0]; h++) {
                for (int k = 0; k < binSize[1]; k++) {
                    for (int l = 0; l < binSize[2]; l++) {
                        for (int symOp = 0;
                             symOp < ((atomic_displacement_parameters
                                           [usedAtomIndices[atomIdx]]
                                               .size() == 1)
                                          ? 1
                                          : symOpMultCount);
                             symOp++) {
                            virtHklTemperature[(h * binSize[1] + k) *
                                                   binSize[2] +
                                               l][atomIdx][symOp] =
                                temperatureFactorMults[0][atomIdx][symOp][h] *
                                temperatureFactorMults[1][atomIdx][symOp][k] *
                                temperatureFactorMults[2][atomIdx][symOp][l] *
                                temperatureFactorMultsSquare[0][atomIdx][symOp]
                                                            [k][l] *
                                temperatureFactorMultsSquare[1][atomIdx][symOp]
                                                            [l][h] *
                                temperatureFactorMultsSquare[2][atomIdx][symOp]
                                                            [h][k];
                        }
                    }
                }
            }
        }
    }

    printStep("virtual hkl temperature factors");

    // [bin][atom][symOpMult | 1]
    std::vector<std::vector<std::vector<Vector3<REAL>>>>
        perBinTemperatureFactorMult;
    perBinTemperatureFactorMult.resize(binCount);
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        perBinTemperatureFactorMult[binIdx].resize(nAtoms);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            perBinTemperatureFactorMult[binIdx][atomIdx].resize(
                (atomic_displacement_parameters[usedAtomIndices[atomIdx]]
                     .size() == 1)
                    ? 1
                    : symOpMultCount);
        }
    }
#pragma omp parallel for num_threads(nThreads) collapse(3)
    for (int binIdx = 0; binIdx < binCount; binIdx++) {
        for (int atom = 0; atom < nAtoms; atom++) {
            for (int i = 0; i < 3; i++) {
                if (atomic_displacement_parameters[usedAtomIndices[atom]]
                        .size() == 1) {
                    perBinTemperatureFactorMult[binIdx][atom][0][i] = std::exp(
                        -2.0 *
                        atomic_displacement_parameters[usedAtomIndices[atom]]
                                                      [0] *
                        (bin000Cartesian[binIdx] * stepCartesian[i]));
                } else {
                    for (int symOp = 0; symOp < symOpMultCount; symOp++) {
                        Matrix3<REAL> currentU = rotatedUs[atom][symOp];
                        perBinTemperatureFactorMult[binIdx][atom][symOp][i] =
                            std::exp(-2.0 *
                                     (bin000Cartesian[binIdx][0] *
                                          stepCartesian[i][0] * currentU(0, 0) +
                                      bin000Cartesian[binIdx][1] *
                                          stepCartesian[i][1] * currentU(1, 1) +
                                      bin000Cartesian[binIdx][2] *
                                          stepCartesian[i][2] * currentU(2, 2) +
                                      (bin000Cartesian[binIdx][0] *
                                           stepCartesian[i][1] +
                                       bin000Cartesian[binIdx][1] *
                                           stepCartesian[i][0]) *
                                          currentU(0, 1) +
                                      (bin000Cartesian[binIdx][0] *
                                           stepCartesian[i][2] +
                                       bin000Cartesian[binIdx][2] *
                                           stepCartesian[i][0]) *
                                          currentU(0, 2) +
                                      (bin000Cartesian[binIdx][1] *
                                           stepCartesian[i][2] +
                                       bin000Cartesian[binIdx][2] *
                                           stepCartesian[i][1]) *
                                          currentU(1, 2)));
                    }
                }
            }
        }
    }

    printStep("per bin temperature factor multipliers");

    // [bin][dir][n][atom][symOpMult | 1]
    std::vector<std::array<std::vector<std::vector<std::vector<REAL>>>, 3>>
        virtHklTMul;
    if constexpr (virtHklTMulFlag) {
        virtHklTMul.resize(binCount);
        for (int binIdx = 0; binIdx < binCount; binIdx++) {
            for (int i = 0; i < 3; i++) {
                virtHklTMul[binIdx][i].resize(binSize[i]);
                for (int n = 0; n < binSize[i]; n++) {
                    virtHklTMul[binIdx][i][n].resize(nAtoms);
                    for (int atom = 0; atom < nAtoms; atom++) {
                        virtHklTMul[binIdx][i][n][atom].resize(
                            (atomic_displacement_parameters
                                 [usedAtomIndices[atom]]
                                     .size() == 1)
                                ? 1
                                : symOpMultCount);
                    }
                }
            }
        }
#pragma omp parallel for num_threads(nThreads) collapse(3)
        for (int binIdx = 0; binIdx < binCount; binIdx++) {
            for (int i = 0; i < 3; i++) {
                for (int atom = 0; atom < nAtoms; atom++) {
                    for (int symOp = 0;
                         symOp <
                         ((atomic_displacement_parameters[usedAtomIndices[atom]]
                               .size() == 1)
                              ? 1
                              : symOpMultCount);
                         symOp++) {
                        REAL mult = 1.0;
                        for (int n = 0; n < binSize[i]; n++) {
                            virtHklTMul[binIdx][i][n][atom][symOp] = mult;
                            mult *= perBinTemperatureFactorMult[binIdx][atom]
                                                               [symOp][i];
                        }
                    }
                }
            }
        }
    }

    printStep("virtual hkl temperature multipliers");

    std::vector<int> usedWfns;
    usedWfns.clear();
    std::vector<int> usedTypes;
    usedTypes.clear();
    std::vector<std::array<int, 2>> usedWfnTypeCombo;
    usedWfnTypeCombo.clear();
    std::vector<int> atomToUsedWfnTypeCombo;
    atomToUsedWfnTypeCombo.resize(nAtoms);
    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
        const int trueAtomIdx = usedAtomIndices[atomIdx];
        const int wfnIdx = atom_to_wfn_map[trueAtomIdx];
        const int typeIdx = atom_to_type_map[trueAtomIdx];

        int currentUsedWfnIdx = -1;
        for (int usedWfnIdx = 0; usedWfnIdx < usedWfns.size(); usedWfnIdx++) {
            if (usedWfns[usedWfnIdx] == wfnIdx) {
                currentUsedWfnIdx = usedWfnIdx;
                break;
            }
        }
        if (currentUsedWfnIdx < 0) {
            currentUsedWfnIdx = usedWfns.size();
            usedWfns.emplace_back(wfnIdx);
        }

        int currentUsedTypeIdx = -1;
        for (int usedTypeIdx = 0; usedTypeIdx < usedTypes.size();
             usedTypeIdx++) {
            if (usedTypes[usedTypeIdx] == typeIdx) {
                currentUsedTypeIdx = usedTypeIdx;
                break;
            }
        }
        if (currentUsedTypeIdx < 0) {
            currentUsedTypeIdx = usedTypes.size();
            usedTypes.emplace_back(typeIdx);
        }

        int currentUsedComboIdx = -1;
        for (int usedComboIdx = 0; usedComboIdx < usedWfnTypeCombo.size();
             usedComboIdx++) {
            if ((usedWfnTypeCombo[usedComboIdx][0] == currentUsedWfnIdx) and
                (usedWfnTypeCombo[usedComboIdx][1] == currentUsedTypeIdx)) {
                currentUsedComboIdx = usedComboIdx;
                break;
            }
        }
        if (currentUsedComboIdx < 0) {
            currentUsedComboIdx = usedWfnTypeCombo.size();
            std::array<int, 2> pair;
            pair[0] = currentUsedWfnIdx;
            pair[1] = currentUsedTypeIdx;
            usedWfnTypeCombo.emplace_back(pair);
        }

        atomToUsedWfnTypeCombo[atomIdx] = currentUsedComboIdx;
    }

    const int wfnCount = usedWfns.size();
    const int typeCount = usedTypes.size();
    const int comboCount = usedWfnTypeCombo.size();

    printStep("wfn and type deduplication");

    // [wfn][l]
    std::vector<std::vector<std::complex<REAL>>> N;
    N.resize(wfnParams.size());
    for (int wfnIdx = 0; wfnIdx < wfnCount; wfnIdx++) {
        N[usedWfns[wfnIdx]].resize(
            wfnParams[usedWfns[wfnIdx]].def_valence_pow.size());
    }
#pragma omp parallel for num_threads(nThreads)
    for (int wfnIdx = 0; wfnIdx < wfnCount; wfnIdx++) {
        const auto &wfn = wfnParams[usedWfns[wfnIdx]];
        const int nl = wfn.def_valence_pow.size();
        for (int l = 0; l < nl; l++) {
            std::complex<REAL> perL;
            const int rest = l % 4;
            if (rest == 0) {
                perL = {four_pi, 0.0};
            } else if (rest == 1) {
                perL = {0.0, four_pi};
            } else if (rest == 2) {
                perL = {-four_pi, 0.0};
            } else {
                perL = {0.0, -four_pi};
            }
            N[usedWfns[wfnIdx]][l] =
                perL * sto_atomic_wfn::stoDensityNormalizationFactor(
                           wfn.def_valence_pow[l], wfn.def_valence_exp);
        }
    }

    printStep("density normalization factor");

    f.resize(hklCount);

    dTarget_dparam.resize(trueNAtoms);
    for (int atomIdx = 0; atomIdx < trueNAtoms; atomIdx++) {
        dTarget_dparam[atomIdx].adp_derivatives.assign(
            atomic_displacement_parameters[atomIdx].size(), 0.0);
        dTarget_dparam[atomIdx].atomic_position_derivatives =
            Vector3<REAL>(0, 0, 0);
        dTarget_dparam[atomIdx].occupancy_derivatives = 0.0;
    }

#ifdef SYCL
    sycl::queue queue;

#define Vector3_to_vec(type, vector) \
    ([vector]() {                    \
        sycl::vec<type, 3> tmp;      \
        tmp[0] = vector.x;           \
        tmp[1] = vector.y;           \
        tmp[2] = vector.z;           \
        return tmp;                  \
    }())

#define vector_of_Vector3_to_buffer(type, value, name_prefix)        \
    std::vector<sycl::vec<type, 3>> name_prefix##_vec(value.size()); \
    for (size_t i = 0; i < value.size(); ++i) {                      \
        Vector3<type> vec = value[i];                                \
        name_prefix##_vec[i] = Vector3_to_vec(type, vec);            \
    }                                                                \
    sycl::buffer<sycl::vec<type, 3>> name_prefix##_buf(name_prefix##_vec);

#define vector_of_objects_to_buffer_of_values_for_property( \
    type, list, prop, name_prefix)                          \
    std::vector<type> name_prefix##_vec(list.size());       \
    for (size_t i = 0; i < list.size(); ++i)                \
        name_prefix##_vec[i] = list[i].prop;                \
    sycl::buffer<type> name_prefix##_buf(name_prefix##_vec);

#define flattened(list)                                              \
    ([list]() {                                                      \
        std::remove_reference<decltype(list)>::type::value_type tmp; \
        for (auto &vec : list)                                       \
            for (auto &el : vec) tmp.push_back(el);                  \
        return tmp;                                                  \
    }())

#define offsets(list)                  \
    ([list]() {                        \
        std::vector<int> offsets;      \
        int offset = 0;                \
        for (auto &vec : list) {       \
            offsets.push_back(offset); \
            offset += vec.size();      \
        }                              \
        return offsets;                \
    }())

#define vecvec_to_buffers(list, buf_name_prefix)                             \
    auto buf_name_prefix##_tmp = flattened(list);                            \
    sycl::buffer<                                                            \
        std::remove_reference<decltype(list)>::type::value_type::value_type> \
        buf_name_prefix##_value_buf(buf_name_prefix##_tmp);                  \
    std::vector<int> buf_name_prefix##_offset_tmp = offsets(list);           \
    sycl::buffer<int> buf_name_prefix##_offset_buf(                          \
        buf_name_prefix##_offset_tmp);

#define vector_times_matrix(vector, col0, col1, col2)  \
    ([vector, col0, col1, col2]() {                    \
        decltype(vector) tmp(0);                       \
        decltype(vector) tmp2 = vector * col0;         \
        for (int i = 0; i < 3; ++i) tmp[0] += tmp2[i]; \
        tmp2 = vector * col1;                          \
        for (int i = 0; i < 3; ++i) tmp[1] += tmp2[i]; \
        tmp2 = vector * col2;                          \
        for (int i = 0; i < 3; ++i) tmp[2] += tmp2[i]; \
        return tmp;                                    \
    }())

#define vec_of_matrix_to_columns(list)                             \
    ([list]() {                                                    \
        std::array<std::vector<sycl::vec<REAL, 3>>, 3> res;        \
        for (auto matrix : list)                                   \
            for (int i = 0; i < 3; ++i) {                          \
                sycl::vec<REAL, 3> tmp;                            \
                for (int j = 0; j < 3; ++j) tmp[j] = matrix(j, i); \
                res[i].push_back(tmp);                             \
            }                                                      \
        return res;                                                \
    }())

#define vec_of_matrix_to_column_buffers(list, buf_name_prefix)      \
    auto buf_name_prefix##_tmp = vec_of_matrix_to_columns(list);    \
    sycl::buffer<sycl::vec<REAL, 3>> buf_name_prefix##_column0_buf( \
        buf_name_prefix##_tmp[0]);                                  \
    sycl::buffer<sycl::vec<REAL, 3>> buf_name_prefix##_column1_buf( \
        buf_name_prefix##_tmp[1]);                                  \
    sycl::buffer<sycl::vec<REAL, 3>> buf_name_prefix##_column2_buf( \
        buf_name_prefix##_tmp[2]);

    sycl::range<1> job_num(hklCount);

    vector_of_Vector3_to_buffer(int, hkl_indices, hkl);
    vector_of_Vector3_to_buffer(int, usedBins, used_bins);

    sycl::buffer<std::complex<REAL>> f_buf(hklCount);
    sycl::buffer<REAL, 2> f_core_buf(sycl::range(hklCount, wfnCount));
    sycl::buffer<REAL, 2> val_buf(sycl::range(hklCount, comboCount));
    sycl::buffer<std::complex<REAL>, 2> sym_op_f_mult_buf(
        sycl::range(hklCount, symOpMultCount));
    sycl::buffer<std::array<int, 2>> used_wfn_type_combo_buf(usedWfnTypeCombo);
    vector_of_objects_to_buffer_of_values_for_property(
        REAL, typeParams, kappa_spherical, type_kappa_spherical);
    vector_of_objects_to_buffer_of_values_for_property(
        REAL, typeParams, p_val, type_p_val);
    vector_of_objects_to_buffer_of_values_for_property(
        REAL, typeParams, kappa_def_valence, type_kappa_def_valence);
    vector_of_objects_to_buffer_of_values_for_property(
        REAL, wfnParams, def_valence_exp, wfn_def_valence_exp);
    sycl::buffer<int> atom_to_used_wfn_type_combo_buf(atomToUsedWfnTypeCombo);
    sycl::buffer<REAL> atomic_occupancy_buf(atomic_occupancy);
    sycl::buffer<REAL> atomic_multiplicity_factor_buf(
        atomic_multiplicity_factor);
    sycl::buffer<int> used_atom_indices_buf(usedAtomIndices);
    sycl::buffer<int> used_wfns_buf(usedWfns);
    sycl::buffer<int> used_types_buf(usedTypes);
    sycl::buffer<int> atom_to_wfn_map_buf(atom_to_wfn_map);
    sycl::buffer<int> atom_to_type_map_buf(atom_to_type_map);
    sycl::buffer<int> all_bins_map_buf(allBinsMap);
    sycl::buffer<int> sym_op_to_mult_buf(symOpToMult);

    std::vector<std::vector<std::vector<REAL>>> type_p_lm_vec;
    for (auto &type : typeParams) type_p_lm_vec.push_back(type.p_lm);
    vecnd_to_buffers(type_p_lm_vec, type_p_lm);

    vecnd_to_buffers(temperatureFactorMultsSquare[0],
                     temperature_factor_mults_square0);
    vecnd_to_buffers(temperatureFactorMultsSquare[1],
                     temperature_factor_mults_square1);
    vecnd_to_buffers(temperatureFactorMultsSquare[2],
                     temperature_factor_mults_square2);
    vecnd_to_buffers(virtHklTMul, virt_hkl_t_mul);

    vecvec_to_buffers(atomic_displacement_parameters,
                      atomic_displacement_parameters);
    size_t atomic_displacement_parameter_count =
        atomic_displacement_parameters_tmp.size();
    size_t atomic_displacement_parameters_size =
        atomic_displacement_parameters.size();

    vecvec_to_buffers(N, n);

    vec_of_matrix_to_column_buffers(symOpMults, sym_op_mults);
    vec_of_matrix_to_column_buffers(local_coordinate_systems,
                                    local_coordinate_systems);

    vecnd_to_buffers(temperatureFactorRoots, temperature_factor_roots);

    std::vector<std::vector<std::vector<sycl::vec<REAL, 3>>>>
        per_bin_temperature_factor_mult_sycl;
    for (auto &vec1 : perBinTemperatureFactorMult) {
        std::vector<std::vector<sycl::vec<REAL, 3>>> tmp1;
        for (auto &vec2 : vec1) {
            std::vector<sycl::vec<REAL, 3>> tmp2;
            for (auto el : vec2) tmp2.push_back(Vector3_to_vec(REAL, el));
            tmp1.push_back(tmp2);
        }
        per_bin_temperature_factor_mult_sycl.push_back(tmp1);
    }
    vecnd_to_buffers(per_bin_temperature_factor_mult_sycl,
                     per_bin_temperature_factor_mult);

    vecnd_to_buffers(temperatureFactorMults[0], temperature_factor_mults0);
    vecnd_to_buffers(temperatureFactorMults[1], temperature_factor_mults1);
    vecnd_to_buffers(temperatureFactorMults[2], temperature_factor_mults2);
    vecnd_to_buffers(phaseFactorMults[0], phase_factor_mults0);
    vecnd_to_buffers(phaseFactorMults[1], phase_factor_mults1);
    vecnd_to_buffers(phaseFactorMults[2], phase_factor_mults2);

    vecnd_to_buffers(phaseFactorRoots, phase_factor_roots);
    vecnd_to_buffers(virtHklTemperature, virt_hkl_temperature);
    vecnd_to_buffers(virtHklPhase, virt_hkl_phase);

    sycl::vec<int, 3> s_min_hkl = Vector3_to_vec(int, minHkl);
    sycl::vec<int, 3> s_step = Vector3_to_vec(int, step);
    sycl::vec<int, 3> s_n_bins = Vector3_to_vec(int, nBins);

    Matrix3d f2c = recUnitCell.getFractionalToCartesianMatrix();
    sycl::vec<REAL, 3> f2c_0;
    for (int i = 0; i < 3; ++i) f2c_0[i] = f2c(0, i);
    sycl::vec<REAL, 3> f2c_1;
    for (int i = 0; i < 3; ++i) f2c_1[i] = f2c(1, i);
    sycl::vec<REAL, 3> f2c_2;
    for (int i = 0; i < 3; ++i) f2c_2[i] = f2c(2, i);

    // Wfn params
    struct wfn_param_meta {
        uint32_t core_coeff_off;
        uint32_t core_coeff_sz;

        uint32_t core_exp_off;
        uint32_t core_exp_sz;

        uint32_t core_pow_off;
        uint32_t core_pow_sz;

        uint32_t valence_coeff_off;
        uint32_t valence_coeff_sz;

        uint32_t valence_exp_off;
        uint32_t valence_exp_sz;

        uint32_t valence_pow_off;
        uint32_t valence_pow_sz;

        uint32_t def_valence_pow_off;
        uint32_t def_valence_pow_sz;
    };

    std::vector<wfn_param_meta> wfn_param_metas;
    std::vector<REAL> core_coeffs;
    std::vector<REAL> core_exps;
    std::vector<int> core_pows;
    std::vector<REAL> valence_coeffs;
    std::vector<REAL> valence_exps;
    std::vector<int> valence_pows;
    std::vector<int> def_valence_pows;

    for (uint32_t i = 0; i < wfnCount; i++) {
        wfn_param_meta meta{};
        const auto &wfn = wfnParams[usedWfns[i]];
#define PROCESS_PARAM(name)                    \
    uint32_t name##_sz = wfn.name.size();      \
    meta.name##_sz = name##_sz;                \
    meta.name##_off = name##s.size();          \
    for (uint32_t j = 0; j < name##_sz; j++) { \
        name##s.push_back(wfn.name[j]);        \
    }

        PROCESS_PARAM(core_coeff);
        PROCESS_PARAM(core_exp);
        PROCESS_PARAM(core_pow);
        PROCESS_PARAM(valence_coeff);
        PROCESS_PARAM(valence_exp);
        PROCESS_PARAM(valence_pow);
        PROCESS_PARAM(def_valence_pow);

#undef PROCESS_PARAM
        wfn_param_metas.push_back(meta);
    }

    sycl::buffer<wfn_param_meta> wfn_param_metas_buff(wfn_param_metas);
    sycl::buffer<REAL> core_coeffs_buff(core_coeffs);
    sycl::buffer<REAL> core_exps_buff(core_exps);
    sycl::buffer<int> core_pows_buff(core_pows);
    sycl::buffer<REAL> valence_coeffs_buff(valence_coeffs);
    sycl::buffer<REAL> valence_exps_buff(valence_exps);
    sycl::buffer<int> valence_pows_buff(valence_pows);
    sycl::buffer<int> def_valence_pows_buff(def_valence_pows);

    std::vector<int> type_p_lm_sizes_vec(typeParams.size());
    for (int i = 0; i < typeParams.size(); ++i)
        type_p_lm_sizes_vec[i] = typeParams[i].p_lm.size();
    sycl::buffer<int> type_p_lm_sizes_buf(type_p_lm_sizes_vec);

    queue.submit([&](sycl::handler &cgh) {
        sycl::accessor f_ax(f_buf, cgh, sycl::write_only);
        sycl::accessor f_core_ax(f_core_buf, cgh, sycl::read_write);
        sycl::accessor val_ax(val_buf, cgh, sycl::read_write);
        sycl::accessor sym_op_f_mult_ax(
            sym_op_f_mult_buf, cgh, sycl::read_write);
        sycl::accessor hkl_ax(hkl_buf, cgh, sycl::read_only);
        sycl::accessor used_bins_ax(used_bins_buf, cgh, sycl::read_only);
        sycl::accessor used_wfn_type_combo_ax(
            used_wfn_type_combo_buf, cgh, sycl::read_only);
        sycl::accessor type_kappa_spherical_ax(
            type_kappa_spherical_buf, cgh, sycl::read_only);
        sycl::accessor type_p_val_ax(type_p_val_buf, cgh, sycl::read_only);
        sycl::accessor type_kappa_def_valence_ax(
            type_kappa_def_valence_buf, cgh, sycl::read_only);
        sycl::accessor wfn_def_valence_exp_ax(
            wfn_def_valence_exp_buf, cgh, sycl::read_only);
        sycl::accessor atom_to_used_wfn_type_combo_ax(
            atom_to_used_wfn_type_combo_buf, cgh, sycl::read_only);
        sycl::accessor atomic_occupancy_ax(
            atomic_occupancy_buf, cgh, sycl::read_only);
        sycl::accessor atomic_multiplicity_factor_ax(
            atomic_multiplicity_factor_buf, cgh, sycl::read_only);
        sycl::accessor used_atom_indices_ax(
            used_atom_indices_buf, cgh, sycl::read_only);
        sycl::accessor used_wfns_ax(used_wfns_buf, cgh, sycl::read_only);
        sycl::accessor used_types_ax(used_types_buf, cgh, sycl::read_only);
        sycl::accessor atomic_displacement_parameters_value_ax(
            atomic_displacement_parameters_value_buf, cgh, sycl::read_only);
        sycl::accessor atomic_displacement_parameters_offset_ax(
            atomic_displacement_parameters_offset_buf, cgh, sycl::read_only);
        sycl::accessor atom_to_wfn_map_ax(
            atom_to_wfn_map_buf, cgh, sycl::read_only);
        sycl::accessor atom_to_type_map_ax(
            atom_to_type_map_buf, cgh, sycl::read_only);
        sycl::accessor n_value_ax(n_value_buf, cgh, sycl::read_only);
        sycl::accessor n_offset_ax(n_offset_buf, cgh, sycl::read_only);
        sycl::accessor sym_op_mults_column0_ax(
            sym_op_mults_column0_buf, cgh, sycl::read_only);
        sycl::accessor sym_op_mults_column1_ax(
            sym_op_mults_column1_buf, cgh, sycl::read_only);
        sycl::accessor sym_op_mults_column2_ax(
            sym_op_mults_column2_buf, cgh, sycl::read_only);
        sycl::accessor local_coordinate_systems_column0_ax(
            local_coordinate_systems_column0_buf, cgh, sycl::read_only);
        sycl::accessor local_coordinate_systems_column1_ax(
            local_coordinate_systems_column1_buf, cgh, sycl::read_only);
        sycl::accessor local_coordinate_systems_column2_ax(
            local_coordinate_systems_column2_buf, cgh, sycl::read_only);
        sycl::accessor all_bins_map_ax(all_bins_map_buf, cgh, sycl::read_only);
        sycl::accessor sym_op_to_mult_ax(
            sym_op_to_mult_buf, cgh, sycl::read_only);
        sycl::accessor type_p_lm_sizes_ax(
            type_p_lm_sizes_buf, cgh, sycl::read_only);
        vecnd_buffer_accessors(temperature_factor_roots);
        vecnd_buffer_accessors(per_bin_temperature_factor_mult);
        vecnd_buffer_accessors(temperature_factor_mults0);
        vecnd_buffer_accessors(temperature_factor_mults1);
        vecnd_buffer_accessors(temperature_factor_mults2);
        vecnd_buffer_accessors(phase_factor_mults0);
        vecnd_buffer_accessors(phase_factor_mults1);
        vecnd_buffer_accessors(phase_factor_mults2);
        vecnd_buffer_accessors(phase_factor_roots);
        vecnd_buffer_accessors(virt_hkl_temperature);
        vecnd_buffer_accessors(virt_hkl_phase);
        vecnd_buffer_accessors(type_p_lm);
        vecnd_buffer_accessors(temperature_factor_mults_square0);
        vecnd_buffer_accessors(temperature_factor_mults_square1);
        vecnd_buffer_accessors(temperature_factor_mults_square2);
        vecnd_buffer_accessors(virt_hkl_t_mul);

        // Wfn params
        sycl::accessor wfn_param_metas_ax(
            wfn_param_metas_buff, cgh, sycl::read_only);
        sycl::accessor core_coeffs_ax(core_coeffs_buff, cgh, sycl::read_only);
        sycl::accessor core_exps_ax(core_exps_buff, cgh, sycl::read_only);
        sycl::accessor core_pows_ax(core_pows_buff, cgh, sycl::read_only);
        sycl::accessor valence_coeffs_ax(
            valence_coeffs_buff, cgh, sycl::read_only);
        sycl::accessor valence_exps_ax(valence_exps_buff, cgh, sycl::read_only);
        sycl::accessor valence_pows_ax(valence_pows_buff, cgh, sycl::read_only);
        sycl::accessor def_valence_pows_ax(
            def_valence_pows_buff, cgh, sycl::read_only);

        cgh.parallel_for<class calculate_sf>(job_num, [=](sycl::id<1> job_id) {
            size_t id = job_id.get(0);

            sycl::vec<int, 3> hkl = hkl_ax[id];

            sycl::vec<int, 3> orig_bin;
            for (int i = 0; i < 3; i++)
                orig_bin[i] =
                    (hkl[i] - s_min_hkl[i]) / (s_step[i] * binSize[i]);

            int binIdx =
                all_bins_map_ax[(orig_bin[0] * s_n_bins[1] + orig_bin[1]) *
                                    s_n_bins[2] +
                                orig_bin[2]];

            std::complex<REAL> f_acc = 0;

            sycl::vec<int, 3> offset;
            for (int i = 0; i < 3; i++)
                offset[i] = ((hkl[i] - s_min_hkl[i]) / s_step[i]) -
                            (used_bins_ax[binIdx][i] * binSize[i]);

            sycl::vec<REAL, 3> cartesian_h(0);
            sycl::vec<REAL, 3> hkl_real{hkl[0], hkl[1], hkl[2]};
            sycl::vec<REAL, 3> tmp = f2c_0 * hkl_real;
            cartesian_h[0] = tmp[0] + tmp[1] + tmp[2];
            tmp = f2c_1 * hkl_real;
            cartesian_h[1] = tmp[0] + tmp[1] + tmp[2];
            tmp = f2c_2 * hkl_real;
            cartesian_h[2] = tmp[0] + tmp[1] + tmp[2];

            REAL h_length = sycl::sqrt(cartesian_h[0] * cartesian_h[0] +
                                       cartesian_h[1] * cartesian_h[1] +
                                       cartesian_h[2] * cartesian_h[2]);

            auto f_core = f_core_ax[id];

            for (int i = 0; i < wfnCount; ++i) {
                wfn_param_meta wfn_meta = wfn_param_metas_ax[i];
                const int kMax = wfn_meta.core_coeff_sz;

                for (int k = 0; k < kMax; k++) {
                    REAL core_coeff_k =
                        core_coeffs_ax[wfn_meta.core_coeff_off + k];
                    REAL core_pow_k = core_pows_ax[wfn_meta.core_pow_off + k];
                    REAL core_exp_k = core_exps_ax[wfn_meta.core_exp_off + k];
                    f_core[i] +=
                        core_coeff_k *
                        gFunction_sycl(0,
                                       core_pow_k + 2,
                                       h_length,
                                       core_exp_k);  // TODO find out why pow+2
                                                     // in all gFunction pow
                }
                f_core[i] *= four_pi;
            }

            auto val = val_ax[id];
            for (int i = 0; i < comboCount; ++i) {
                const auto &combo =
                    used_wfn_type_combo_ax[i];  // Use extract combo
                const auto type_i = used_types_ax[combo[1]];
                wfn_param_meta wfn_meta = wfn_param_metas_ax[combo[0]];
                const int kMax = wfn_meta.valence_coeff_sz;
                const auto h = h_length / type_kappa_spherical_ax[type_i];
                for (int k = 0; k < kMax; k++) {
                    REAL valence_pow_k =
                        valence_pows_ax[wfn_meta.valence_pow_off + k];
                    REAL valence_exp_k =
                        valence_exps_ax[wfn_meta.valence_exp_off + k];
                    REAL valence_coeff_k =
                        valence_coeffs_ax[wfn_meta.valence_coeff_off + k];
                    val[i] +=
                        valence_coeff_k *
                        gFunction_sycl(0, valence_pow_k + 2, h, valence_exp_k);
                }
                val[i] *= type_p_val_ax[type_i] * four_pi;
            }

            sycl::vec<REAL, 3> h_versor = closeToZero(h_length)
                                              ? sycl::vec<REAL, 3>{0}
                                              : cartesian_h / h_length;

            int virtHklPhaseCurrent;
            if constexpr (virtHklPhaseFlag)
                virtHklPhaseCurrent =
                    (offset[0] * binSize[1] + offset[1]) * binSize[2] +
                    offset[2];

            int virtHklTemperatureCurrent;
            if constexpr (virtHklTemperatureFlag)
                virtHklTemperatureCurrent =
                    (offset[0] * binSize[1] + offset[1]) * binSize[2] +
                    offset[2];

            for (int atom_i = 0; atom_i < nAtoms; ++atom_i) {
                std::complex<REAL> perAtomF = 0.0;
                auto symOpFMult = sym_op_f_mult_ax[id];
                for (int i = 0; i < symOpMultCount; ++i) symOpFMult[i] = 0;
                for (int symOpIdx = 0; symOpIdx < nSymOps; symOpIdx++) {
                    std::complex<REAL> localF = index_vecnd_buffer(
                        phase_factor_roots, binIdx, atom_i, symOpIdx);
                    if constexpr (virtHklPhaseFlag) {
                        localF *= index_vecnd_buffer(virt_hkl_phase,
                                                     virtHklPhaseCurrent,
                                                     atom_i,
                                                     symOpIdx);
                    } else {
                        localF *= index_vecnd_buffer(phase_factor_mults0,
                                                     atom_i,
                                                     symOpIdx,
                                                     offset[0]) *
                                  index_vecnd_buffer(phase_factor_mults1,
                                                     atom_i,
                                                     symOpIdx,
                                                     offset[1]) *
                                  index_vecnd_buffer(phase_factor_mults2,
                                                     atom_i,
                                                     symOpIdx,
                                                     offset[2]);
                    }
                    symOpFMult[sym_op_to_mult_ax[symOpIdx]] += localF;
                }

                for (int symOpIdx = 0; symOpIdx < symOpMultCount; symOpIdx++) {
                    bool iso = (used_atom_indices_ax[atom_i] + 1 ==
                                        atomic_displacement_parameters_size
                                    ? atomic_displacement_parameter_count
                                    : atomic_displacement_parameters_offset_ax
                                          [used_atom_indices_ax[atom_i] + 1]) -
                                   atomic_displacement_parameters_offset_ax
                                       [used_atom_indices_ax[atom_i]] ==
                               1;

                    REAL localF = index_vecnd_buffer(temperature_factor_roots,
                                                     binIdx,
                                                     atom_i,
                                                     iso ? 0 : symOpIdx);
                    if constexpr (virtHklTemperatureFlag) {
                        localF *= index_vecnd_buffer(virt_hkl_temperature,
                                                     virtHklTemperatureCurrent,
                                                     atom_i,
                                                     iso ? 0 : symOpIdx);
                    } else {
                        localF *=
                            index_vecnd_buffer(temperature_factor_mults0,
                                               atom_i,
                                               iso ? 0 : symOpIdx,
                                               offset[0]) *
                            index_vecnd_buffer(temperature_factor_mults1,
                                               atom_i,
                                               iso ? 0 : symOpIdx,
                                               offset[1]) *
                            index_vecnd_buffer(temperature_factor_mults2,
                                               atom_i,
                                               iso ? 0 : symOpIdx,
                                               offset[2]) *
                            index_vecnd_buffer(temperature_factor_mults_square0,
                                               atom_i,
                                               iso ? 0 : symOpIdx,
                                               offset[1],
                                               offset[2]) *
                            index_vecnd_buffer(temperature_factor_mults_square1,
                                               atom_i,
                                               iso ? 0 : symOpIdx,
                                               offset[2],
                                               offset[0]) *
                            index_vecnd_buffer(temperature_factor_mults_square2,
                                               atom_i,
                                               iso ? 0 : symOpIdx,
                                               offset[0],
                                               offset[1]);
                    }

                    for (int i = 0; i < 3; i++) {
                        if constexpr (virtHklTMulFlag)
                            localF *= index_vecnd_buffer(virt_hkl_t_mul,
                                                         binIdx,
                                                         i,
                                                         offset[i],
                                                         atom_i,
                                                         iso ? 0 : symOpIdx);
                        else
                            localF *=
                                sycl::pow(index_vecnd_buffer(
                                              per_bin_temperature_factor_mult,
                                              binIdx,
                                              atom_i,
                                              iso ? 0 : symOpIdx)[i],
                                          offset[i]);
                    }

                    const int wfn_i =
                        atom_to_wfn_map_ax[used_atom_indices_ax[atom_i]];
                    const auto type_i =
                        atom_to_type_map_ax[used_atom_indices_ax[atom_i]];

                    auto sym_op_mult_col0 = sym_op_mults_column0_ax[symOpIdx];
                    auto sym_op_mult_col1 = sym_op_mults_column1_ax[symOpIdx];
                    auto sym_op_mult_col2 = sym_op_mults_column2_ax[symOpIdx];
                    auto local_coordinate_system_col0 =
                        local_coordinate_systems_column0_ax
                            [used_atom_indices_ax[atom_i]];
                    auto local_coordinate_system_col1 =
                        local_coordinate_systems_column1_ax
                            [used_atom_indices_ax[atom_i]];
                    auto local_coordinate_system_col2 =
                        local_coordinate_systems_column2_ax
                            [used_atom_indices_ax[atom_i]];
                    auto h_times_sym = vector_times_matrix(h_versor,
                                                           sym_op_mult_col0,
                                                           sym_op_mult_col1,
                                                           sym_op_mult_col2);
                    const auto h =
                        vector_times_matrix(h_times_sym,
                                            local_coordinate_system_col0,
                                            local_coordinate_system_col1,
                                            local_coordinate_system_col2);

                    const int nl = std::min(
                        (int)wfn_param_metas_ax[wfn_i].def_valence_pow_sz,
                        type_p_lm_sizes_ax[type_i]);
                    std::complex<REAL> dval;
                    const int used_wfn_i = used_wfn_type_combo_ax
                        [atom_to_used_wfn_type_combo_ax[atom_i]][0];
                    for (int l = 0; l < nl; l++) {
                        // may be ordered differently than in publication
                        // because the publication doesn't seem to have a
                        // consistent ordering of arguments passed to g
                        REAL def_valence_pow_l =
                            def_valence_pows_ax[wfn_param_metas_ax[used_wfn_i]
                                                    .def_valence_pow_off +
                                                l];
                        const REAL multPerL = gFunction_sycl(
                            l,
                            def_valence_pow_l + 2,
                            h_length / type_kappa_def_valence_ax[type_i],
                            wfn_def_valence_exp_ax[wfn_i]);

                        REAL sumPerM = 0.0;
                        for (int m = -l; m <= l; m++) {
                            sumPerM += index_vecnd_buffer(
                                           type_p_lm, type_i, l, m + l) *
                                       densityNormalizedSycl(h, l, m);
                        }
                        dval += n_value_ax[n_offset_ax[wfn_i] + l] *
                                (multPerL * sumPerM);
                    }

                    perAtomF +=
                        symOpFMult[symOpIdx] * localF *
                        (dval +
                         1.0 *
                             f_core[used_wfn_type_combo_ax
                                        [atom_to_used_wfn_type_combo_ax[atom_i]]
                                        [0]] +
                         val[atom_to_used_wfn_type_combo_ax[atom_i]]);
                }
                f_acc +=
                    perAtomF *
                    atomic_occupancy_ax[used_atom_indices_ax[atom_i]] *
                    atomic_multiplicity_factor_ax[used_atom_indices_ax[atom_i]];
            }
            f_ax[id] = f_acc;
        });
    });

    sycl::host_accessor f_ax(f_buf, sycl::read_only);

    for (int i = 0; i < hklCount; ++i) f[i] = f_ax[i];
#else
#pragma omp parallel for num_threads(nThreads) schedule(guided)
    for (int hklIdx = 0; hklIdx < hklCount; hklIdx++) {
        printInLoop("begin");
        valueLog("hklIdx", hklIdx);
        valueLog("hkl_indices[hklIdx][0]", hkl_indices[hklIdx][0]);
        valueLog("hkl_indices[hklIdx][1]", hkl_indices[hklIdx][1]);
        valueLog("hkl_indices[hklIdx][2]", hkl_indices[hklIdx][2]);
        Vector3<int> origBin;
        for (int i = 0; i < 3; i++) {
            origBin[i] =
                (hkl_indices[hklIdx][i] - minHkl[i]) / (step[i] * binSize[i]);
        }
        int binIdx =
            allBinsMap[(origBin[0] * nBins[1] + origBin[1]) * nBins[2] +
                       origBin[2]];

        printInLoop("bin index");
        valueLog("binIdx", binIdx);

        std::complex<REAL> f_acc = 0.0;

        Vector3<int> offset;
        for (int i = 0; i < 3; i++) {
            offset[i] = ((hkl_indices[hklIdx][i] - minHkl[i]) / step[i]) -
                        (usedBins[binIdx][i] * binSize[i]);
            valueLog("offset[i]", offset[i]);
        }

        printInLoop("offset");

        Vector3<REAL> cartesianH;
        recUnitCell.fractionalToCartesian({hkl_indices[hklIdx][0],
                                           hkl_indices[hklIdx][1],
                                           hkl_indices[hklIdx][2]},
                                          cartesianH);

        valueLog("cartesianH[0]", cartesianH[0]);
        valueLog("cartesianH[1]", cartesianH[1]);
        valueLog("cartesianH[2]", cartesianH[2]);

        printInLoop("fractional to cartesian conversion");

        REAL hLength = 0.0;
        for (int i = 0; i < 3; i++) {
            hLength += square(cartesianH[i]);
        }
        hLength = sqrt(hLength);

        printInLoop("length of h");

        std::vector<REAL> f_core;
        f_core.resize(wfnCount);
#pragma omp simd
        for (int wfnIdx = 0; wfnIdx < wfnCount; wfnIdx++) {
            const auto &wfn = wfnParams[usedWfns[wfnIdx]];
            const int kMax = wfn.core_coeff.size();
            for (int k = 0; k < kMax; k++)
                f_core[wfnIdx] +=
                    wfn.core_coeff[k] *
                    sto_scattering::gFunction(
                        0,
                        wfn.core_pow[k] + 2,
                        hLength,
                        wfn.core_exp[k]);  // TODO find out why pow+2
                                           // in all gFunction pow
            f_core[wfnIdx] *= four_pi;
        }

        printInLoop("core factor");

        std::vector<REAL> val;
        val.resize(comboCount);
#pragma omp simd
        for (int comboIdx = 0; comboIdx < comboCount; comboIdx++) {
            const auto &combo = usedWfnTypeCombo[comboIdx];
            const auto &wfn = wfnParams[usedWfns[combo[0]]];
            const auto &type = typeParams[usedTypes[combo[1]]];
            const int kMax = wfn.valence_coeff.size();
            const auto h = hLength / type.kappa_spherical;
            for (int k = 0; k < kMax; k++)
                val[comboIdx] +=
                    wfn.valence_coeff[k] *
                    sto_scattering::gFunction(
                        0, wfn.valence_pow[k] + 2, h, wfn.valence_exp[k]);
            val[comboIdx] *= type.p_val * four_pi;
        }

        printInLoop("valence component");

        Vector3<REAL> hVersor;
        for (int i = 0; i < 3; i++) {
            hVersor[i] = closeToZero(hLength) ? 0.0 : cartesianH[i] / hLength;
        }

        printInLoop("h versor");

        std::vector<std::vector<std::complex<REAL>>> virtHklPhaseCurrent;
        if constexpr (virtHklPhaseFlag)
            virtHklPhaseCurrent =
                virtHklPhase[(offset[0] * binSize[1] + offset[1]) * binSize[2] +
                             offset[2]];

        printInLoop("virtual hkl phase factors retreival");

        std::vector<std::vector<REAL>> virtHklTemperatureCurrent;
        if constexpr (virtHklTemperatureFlag)
            virtHklTemperatureCurrent =
                virtHklTemperature[(offset[0] * binSize[1] + offset[1]) *
                                       binSize[2] +
                                   offset[2]];

        printInLoop("virtual hkl temperature factors retreival");

        std::array<std::vector<std::vector<REAL>>, 3> virtHklTMulCurrent;
        if constexpr (virtHklTMulFlag) {
            for (int i = 0; i < 3; i++)
                virtHklTMulCurrent[i] = virtHklTMul[binIdx][i][offset[i]];
        }

        printInLoop("virtual hkl temperature multipliers retreival");

#pragma omp simd
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++) {
            std::complex<REAL> perAtomF = 0.0;
            std::vector<std::complex<REAL>> symOpFMult;
            symOpFMult.resize(symOpMultCount);
            for (int symOpIdx = 0; symOpIdx < nSymOps; symOpIdx++) {
                std::complex<REAL> localF =
                    phaseFactorRoots[binIdx][atomIdx][symOpIdx];
                valueLog("0 - localF", localF);
                if constexpr (virtHklPhaseFlag) {
                    localF *= virtHklPhaseCurrent[atomIdx][symOpIdx];
                } else {
                    localF *=
                        phaseFactorMults[0][atomIdx][symOpIdx][offset[0]] *
                        phaseFactorMults[1][atomIdx][symOpIdx][offset[1]] *
                        phaseFactorMults[2][atomIdx][symOpIdx][offset[2]];
                }
                valueLog("1 - localF", localF);
                symOpFMult[symOpToMult[symOpIdx]] += localF;
            }

            // const auto hVersorLocal =
            // local_coordinate_systems[usedAtomIndices[atomIdx]] * hVersor;

            bool iso = atomic_displacement_parameters[usedAtomIndices[atomIdx]]
                           .size() == 1;
            const int wfnIdx = atom_to_wfn_map[usedAtomIndices[atomIdx]];
            const auto &wfn = wfnParams[wfnIdx];
            const auto &type =
                typeParams[atom_to_type_map[usedAtomIndices[atomIdx]]];
            const auto anomalous =
                anomalous_dispersion.empty()
                    ? wfn.anomalous_scattering
                    : anomalous_dispersion[usedAtomIndices[atomIdx]];

            std::array<std::complex<REAL>, 6> d_adp_p;
            std::array<std::complex<REAL>, 3> d_xyz_p;

            for (int symOpIdx = 0; symOpIdx < symOpMultCount; symOpIdx++) {
                REAL localF =
                    temperatureFactorRoots[binIdx][atomIdx][iso ? 0 : symOpIdx];
                valueLog("2 - localF", localF);
                if constexpr (virtHklTemperatureFlag) {
                    localF *=
                        virtHklTemperatureCurrent[atomIdx][iso ? 0 : symOpIdx];
                } else {
                    localF *=
                        temperatureFactorMults[0][atomIdx][iso ? 0 : symOpIdx]
                                              [offset[0]] *
                        temperatureFactorMults[1][atomIdx][iso ? 0 : symOpIdx]
                                              [offset[1]] *
                        temperatureFactorMults[2][atomIdx][iso ? 0 : symOpIdx]
                                              [offset[2]] *
                        temperatureFactorMultsSquare[0][atomIdx]
                                                    [iso ? 0 : symOpIdx]
                                                    [offset[1]][offset[2]] *
                        temperatureFactorMultsSquare[1][atomIdx]
                                                    [iso ? 0 : symOpIdx]
                                                    [offset[2]][offset[0]] *
                        temperatureFactorMultsSquare[2][atomIdx]
                                                    [iso ? 0 : symOpIdx]
                                                    [offset[0]][offset[1]];
                }
                valueLog("3 - localF", localF);

                for (int i = 0; i < 3; i++) {
                    if constexpr (virtHklTMulFlag)
                        localF *=
                            virtHklTMulCurrent[i][atomIdx][iso ? 0 : symOpIdx];
                    else
                        localF *= pow(
                            perBinTemperatureFactorMult[binIdx][atomIdx]
                                                       [iso ? 0 : symOpIdx][i],
                            offset[i]);
                }

                valueLog("4 - localF", localF);

                const auto h =
                    (hVersor * symOpMults[symOpIdx]) *
                    local_coordinate_systems[usedAtomIndices[atomIdx]];
                valueLog("(square(h[0]) + square(h[1]) + square(h[2]))",
                         (square(h[0]) + square(h[1]) + square(h[2])));

                const int nl =
                    std::min(wfn.def_valence_pow.size(), type.p_lm.size());
                std::complex<REAL> dval;
                for (int l = 0; l < nl; l++) {
                    // may be ordered differently than in publication because
                    // the publication doesn't seem to have a consistent
                    // ordering of arguments passed to g
                    const REAL multPerL = sto_scattering::gFunction(
                        l,
                        wfn.def_valence_pow[l] + 2,
                        hLength / type.kappa_def_valence,
                        wfn.def_valence_exp);

                    REAL sumPerM = 0.0;
                    for (int m = -l; m <= l; m++) {
                        sumPerM += type.p_lm[l][m + l] *
                                   real_spherical_harmonics::densityNormalized(
                                       h, l, m);
                    }
                    dval += N[wfnIdx][l] * (multPerL * sumPerM);
                }

                valueLog("symOpFMult[symOpIdx]", symOpFMult[symOpIdx]);
                valueLog("dval", dval);
                valueLog(
                    "f_core[usedWfnTypeCombo[atomToUsedWfnTypeCombo[atomIdx]]["
                    "0]]",
                    f_core[usedWfnTypeCombo[atomToUsedWfnTypeCombo[atomIdx]]
                                           [0]]);
                valueLog("val[atomToUsedWfnTypeCombo[atomIdx]]",
                         val[atomToUsedWfnTypeCombo[atomIdx]]);

                perAtomF +=
                    symOpFMult[symOpIdx] * localF *
                    (dval +
                     f_core[usedWfnTypeCombo[atomToUsedWfnTypeCombo[atomIdx]]
                                            [0]] +
                     val[atomToUsedWfnTypeCombo[atomIdx]] + anomalous);

                const auto h_rot = (cartesianH * symOpMults[symOpIdx]);
                if (derivativesSwitch.d_adp and not iso) {
                    const std::array<REAL, 6> d_adp_p_part = {
                        h_rot[0] * h_rot[0],
                        h_rot[1] * h_rot[1],
                        h_rot[2] * h_rot[2],
                        h_rot[0] * h_rot[1] * 2,
                        h_rot[0] * h_rot[2] * 2,
                        h_rot[1] * h_rot[2] * 2};
                    for (int i = 0; i < 6; i++) {
                        d_adp_p[i] +=
                            d_adp_p_part[i] * symOpFMult[symOpIdx] * localF *
                            (dval +
                             f_core[usedWfnTypeCombo
                                        [atomToUsedWfnTypeCombo[atomIdx]][0]] +
                             val[atomToUsedWfnTypeCombo[atomIdx]] + anomalous);
                    }
                }
                if (derivativesSwitch.d_xyz) {
                    for (int i = 0; i < 3; i++) {
                        d_xyz_p[i] +=
                            h_rot[i] * symOpFMult[symOpIdx] * localF *
                            (dval +
                             f_core[usedWfnTypeCombo
                                        [atomToUsedWfnTypeCombo[atomIdx]][0]] +
                             val[atomToUsedWfnTypeCombo[atomIdx]] + anomalous);
                    }
                }
            }
            f_acc += perAtomF * atomic_occupancy[usedAtomIndices[atomIdx]] *
                     atomic_multiplicity_factor[usedAtomIndices[atomIdx]];

            if (derivativesSwitch.d_xyz) {
                for (int i = 0; i < 3; i++) {
                    dTarget_dparam[usedAtomIndices[atomIdx]]
                        .atomic_position_derivatives[i] -=
                        (dTarget_df[hklIdx].real() * d_xyz_p[i].imag() +
                         dTarget_df[hklIdx].imag() * d_xyz_p[i].real()) *
                        atomic_occupancy[usedAtomIndices[atomIdx]] *
                        atomic_multiplicity_factor[usedAtomIndices[atomIdx]] *
                        two_pi;
                }
            }
            if (derivativesSwitch.d_adp) {
                if (iso) {
                    const auto d_adp_part = perAtomF * square(hLength);
                    dTarget_dparam[usedAtomIndices[atomIdx]]
                        .adp_derivatives[0] +=
                        (d_adp_part.imag() * dTarget_df[hklIdx].imag() -
                         d_adp_part.real() * dTarget_df[hklIdx].real()) *
                        two_pi_squared *
                        atomic_occupancy[usedAtomIndices[atomIdx]] *
                        atomic_multiplicity_factor[usedAtomIndices[atomIdx]];
                } else {
                    for (int i = 0; i < 6; i++)
                        dTarget_dparam[usedAtomIndices[atomIdx]]
                            .adp_derivatives[i] +=
                            (d_adp_p[i].imag() * dTarget_df[hklIdx].imag() -
                             d_adp_p[i].real() * dTarget_df[hklIdx].real()) *
                            two_pi_squared *
                            atomic_occupancy[usedAtomIndices[atomIdx]] *
                            atomic_multiplicity_factor
                                [usedAtomIndices[atomIdx]];
                }
            }
            if (derivativesSwitch.d_occ) {
                const auto d_occ_part =
                    perAtomF *
                    atomic_multiplicity_factor[usedAtomIndices[atomIdx]];
                dTarget_dparam[usedAtomIndices[atomIdx]]
                    .occupancy_derivatives +=
                    d_occ_part.real() * dTarget_df[hklIdx].real() -
                    d_occ_part.imag() * dTarget_df[hklIdx].imag();
            }
        }
        printInLoop("atom loop");
        f[hklIdx] = f_acc;
        printInLoop("returning results");
    }
#endif

    printStep("main loop");

    // TODO add centrosymmetry support
    // TODO derivatives
    // TODO find out how to implement implement wfn.anomalous_scattering,
    // electron
    /*
     *    bool mUseIAM;
     *    std::vector<std::string> mIamAtomType;
     *    std::vector<int> mAtomToIamTypeMap;
     *    std::vector<NGaussianFormFactor> mIamFormFactors;
     */ // TODO use this
}  // calculateSF_parallel_2

void HansenCoppens_SF_Engine4::electronScatteringAt000(
    const std::vector<int> &atomic_numbers, std::vector<double> &f) {
    map<int, int> z_2_ff_idx;
    set<int> unique_z(atomic_numbers.begin(), atomic_numbers.end());
    vector<int> unique_z_vec(unique_z.begin(), unique_z.end());
    for (int i = 0; i < unique_z_vec.size(); i++)
        z_2_ff_idx[unique_z_vec[i]] = i;

    vector<double> ff_type;
    for (int i = 0; i < unique_z_vec.size(); i++)
        ff_type.push_back(
            n_gaussian_form_factors_table::getFormFactor(
                periodic_table::symbol(unique_z_vec[i]), "electron-IT")
                .calculate_h(0.0));

    int atomIdx, nAtoms = atomic_numbers.size();
    f.resize(nAtoms);
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        f[atomIdx] = ff_type[z_2_ff_idx[atomic_numbers[atomIdx]]];
}

//---------------
/*
 * void calculateSF(
 *    const UnitCell &unitCell,
 *    const std::vector<sf_engine_data_types::HC_WfnParam> &wfn_parameters,
 *    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
 *    const std::vector<int> &atom_to_wfn_map,
 *    const std::vector<int> &atom_to_type_map,
 *    const std::vector<Vector3<REAL> > &atomicPositions,
 *    const std::vector<std::vector<REAL> > &atomic_displacement_parameters,
 *    const std::vector<REAL> &atomic_occupancy,
 *    const std::vector<REAL> &atomic_multiplicity_factor,
 *    const std::vector<Matrix3<REAL> > &local_coordinate_systems,
 *    const std::vector<sf_engine_data_types::SymmetryOperation>
 * &symmetry_operations, bool centrosymmetric, const Vector3<REAL>
 * &inversionTranslation, const std::vector<Vector3<REAL> > &h_vectors, const
 * std::vector<Vector3i >& hkl_indices, std::vector<std::complex<REAL> > &f,
 *    std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
 *    const std::vector<std::complex<REAL> > &dTarget_df,
 *    const std::vector<bool> &include_atom_contribution,
 *    int nThreads);
 *
 */

void HansenCoppens_SF_Engine4::select_P10P20_atoms(
    const std::vector<sf_engine_data_types::HC_TypeParam> &type_parameters,
    const std::vector<int> &atom_to_type_map,
    std::vector<bool> &atom_selection) {
    int nTypes = type_parameters.size();
    int nAtoms = atom_to_type_map.size();
    vector<bool> pz_dz_type(nTypes, false);
    for (int typeIdx = 0; typeIdx < nTypes; typeIdx++) {
        auto const &plms = type_parameters[typeIdx].p_lm;
        if (plms.size() == 3)
            if (plms[0][0] == 0.0) {
                pz_dz_type[typeIdx] = true;
                for (int l = 1; l <= 2; l++)
                    for (int i = 0; i < 2 * l + 1; i++)
                        if (plms[l][i] != 0.0)
                            if (l != i) pz_dz_type[typeIdx] = false;
            }
    }
    atom_selection.resize(nAtoms);
    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        atom_selection[atomIdx] = pz_dz_type[atom_to_type_map[atomIdx]];
}

}  // namespace discamb
