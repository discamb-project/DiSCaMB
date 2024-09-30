#ifndef _DISCAMB_MATHUTILITIES_MATHUTILITIES_H_
#define _DISCAMB_MATHUTILITIES_MATHUTILITIES_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace discamb {
    /**
    * \defgroup MathUtilities MathUtilities
    \brief Math utilities - 3d algebra, numerical integration, spherical harmonics, interpolation,...
    * @{
    */

namespace math_utilities{

inline int roundInt( double r ) { return (r > 0.0) ? int(r + 0.5) : int(r - 0.5); }

template<typename T>
T factorial(int n);

template<typename T>
T binomialCoefficient(int n, int k);

template<typename T>
T power(const T& base, int power);

template<typename T>
T doubleFactorial(int n);


}
/** @} */

//##################### IMPLEMENTATION ####################

template<typename T>
T math_utilities::factorial(
    int n)
{
    int i;
    T result = 1;
    for (i = 1; i <= n; i++)
        result *= i;
    return result;
}

template <typename T>
T math_utilities::binomialCoefficient(
    int n,
    int k)
{
    return factorial<T>(n) / (factorial<T>(n - k) * factorial<T>(k));
}

template<typename T>
T math_utilities::power(
    const T& base,
    int exponent)
{
    T result;
    int i;
    result = 1;
    for (i = 1; i <= exponent; i++)
        result *= base;
    return result;
}

template<typename T>
T math_utilities::doubleFactorial(
    int n)
{
    int i;
    T result = 1;
    for (i = 2 - n % 2; i <= n; i += 2)
        result *= i;
    return result;
}

}



#endif /*_DISCAMB_MATHUTILITIES_MATHUTILITIES_H_*/
