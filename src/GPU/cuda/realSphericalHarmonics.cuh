// vim: noet: sw=3: ts=3: foldmethod=marker: foldlevel=0
#include<cmath>
//assumes x,y,z are normalized

//REAL SPHERICAL HARMONICS ORIGINAL-----------------------------------------{{{
	template<int l, int m, typename T>
__host__ __device__ __forceinline__ T realSphericalHarmonicsOriginal(
		T const x,
		T const y,
		T const z)
{

	T const oneOverSqrtPi = ((T)1.0)/(T)sqrt(M_PI);
	T const r = x*x + y*y + z*z; //TODO: mam nadzieje ze to sie zoptymalizuje (10.07.14 by szmigacz)

	T value;


	if(l == 0)
		value = ((T)0.5)*oneOverSqrtPi;

	if(l == 1) {
		if(m == -1)
			value = (T)(sqrt((T)0.75)*oneOverSqrtPi*y);
		if(m == 0)
			value = (T)(sqrt((T)0.75)*oneOverSqrtPi*z);
		if(m == 1)
			value = (T)(sqrt((T)0.75)*oneOverSqrtPi*x);
	}

	if(l == 2) {
		if(m == -2)
			value = (T)((T)0.5*(T)sqrt((T)15.0)*oneOverSqrtPi*x*y);
		if(m == -1)
			value = (T)((T)0.5*(T)sqrt((T)15.0)*oneOverSqrtPi*y*z);
		if(m == 0)
			value = (T)((T)0.25*(T)sqrt((T)5.0)*oneOverSqrtPi*((T)2.0*z*z-x*x-y*y));
		if(m == 1)
			value = (T)((T)0.5*(T)sqrt((T)15.0)*oneOverSqrtPi*z*x);
		if(m == 2)
			value = (T)((T)0.25*(T)sqrt((T)15.0)*oneOverSqrtPi*(x*x-y*y));
	}

	if(l == 3) {
		if (m == -3)
			value = ((T)0.25*(T)sqrt(35.0/2.0)*oneOverSqrtPi* y*( (T)3.0*x*x  - y*y));
		if (m == -2)
			value = ((T)0.50*(T)sqrt(105.0)*oneOverSqrtPi*x*y*z);
		if (m == -1)
			value = ((T)0.25*(T)sqrt(21.0/2.0)*oneOverSqrtPi*y*((T)4.0*z*z - x*x - y*y));
		if (m == 0)
			value = ((T)0.25*(T)sqrt(7.0)*oneOverSqrtPi*z*( (T)2.0*z*z - (T)3.0*x*x - (T)3.0*y*y));
		if (m == 1)
			value = ((T)0.25*(T)sqrt(21.0/2.0)*oneOverSqrtPi*x*((T)4.0*z*z - x*x - y*y));
		if (m == 2)
			value = ((T)0.25*(T)sqrt(105.0)*oneOverSqrtPi*z*(x*x - y*y));
		if (m == 3)
			value = ((T)0.25*(T)sqrt(35.0/2.0)*oneOverSqrtPi* x*( x*x  - (T)3.0*y*y));
	}

	if(l == 4) {
		if (m == -4)
			value = ( (T)0.75*(T)sqrt(35.0)*oneOverSqrtPi*x*y*(x*x - y*y));
		if (m == -3)
			value = ( (T)0.75*(T)sqrt(35.0/2.0)*oneOverSqrtPi*y*z*((T)3.0*x*x - y*y));
		if (m == -2)
			value = ( (T)0.75*(T)sqrt(5.0)*oneOverSqrtPi*x*y*((T)7.0*z*z - r*r));
		if (m == -1)
			value = ( (T)0.75*(T)sqrt(5.0/2.0)*oneOverSqrtPi*y*z*((T)7.0*z*z - (T)3.0*r*r));
		if (m == 0)
			value = ( (T)(3.0/16.0)*oneOverSqrtPi*((T)35.0*z*z*z*z -(T)30.0*z*z*r*r + (T)3.0*r*r*r*r));
		if (m == 1)
			value = ( (T)0.75*(T)sqrt(5.0/2.0)*oneOverSqrtPi*x*z*((T)7.0*z*z - (T)3.0*r*r));
		if (m == 2)
			value = ( (T)(3.0/8.0)*(T)sqrt(5.0)*oneOverSqrtPi*(x*x - y*y)*((T)7.0*z*z - r*r));
		if (m == 3)
			value = ( (T)0.75*(T)sqrt(35.0/2.0)*oneOverSqrtPi*x*z*(x*x - (T)3.0*y*y));
		if (m == 4)
			value = ( (T)(3.0/16.0)*(T)sqrt(35.0)*oneOverSqrtPi*(x*x*(x*x - (T)3.0*y*y) - y*y*((T)3.0*x*x - y*y)));
	}
	return value;
}
//END-REAL SPHERICAL HARMONICS ORIGINAL-------------------------------------}}}

//REAL SPHERICAL HARMONICS OPT----------------------------------------------{{{
	template<int l, int m, typename T>
__host__ __device__ __forceinline__ T realSphericalHarmonicsOpt(
		T const x,
		T const y,
		T const z)
{

	T const oneOverSqrtPi = ((T)1.0)/(T)sqrt(M_PI);
	T const r = x*x + y*y + z*z; //TODO: mam nadzieje ze to sie zoptymalizuje (10.07.14 by szmigacz). // raczej sie zoptymalizowalo (4.08 by kkiewicz)

	T value;

	T const x_pow2 = x*x;
	T const y_pow2 = y*y;
	T const z_pow2 = z*z;
	T const r_pow2 = r*r;

	T const z_pow4 = z_pow2 * z_pow2;
	T const r_pow4 = r_pow2 * r_pow2;

	if(l == 0)
		value = ((T)0.5)*oneOverSqrtPi;

	if(l == 1) {
		if(m == -1)
			value = (T)(sqrt((T)0.75)*oneOverSqrtPi*y);
		if(m == 0)
			value = (T)(sqrt((T)0.75)*oneOverSqrtPi*z);
		if(m == 1)
			value = (T)(sqrt((T)0.75)*oneOverSqrtPi*x);
	}

	if(l == 2) {
		if(m == -2)
			value = (T)((T)0.5*(T)sqrt(15.0)*oneOverSqrtPi*x*y);
		if(m == -1)
			value = (T)((T)0.5*(T)sqrt(15.0)*oneOverSqrtPi*y*z);
		if(m == 0)
			value = (T)((T)0.25*(T)sqrt(5.0)*oneOverSqrtPi*((T)2.0*z_pow2-x_pow2-y_pow2));
		if(m == 1)
			value = (T)((T)0.5*(T)sqrt(15.0)*oneOverSqrtPi*z*x);
		if(m == 2)
			value = (T)((T)0.25*(T)sqrt(15.0)*oneOverSqrtPi*(x_pow2-y_pow2));
	}

	if(l == 3) {
		if (m == -3)
			value = ((T)0.25*(T)sqrt(17.5)*oneOverSqrtPi* y*( (T)3.0*x_pow2  - y_pow2));
		if (m == -2)
			value = ((T)0.50*(T)sqrt(105.0)*oneOverSqrtPi*x*y*z);
		if (m == -1)
			value = ((T)0.25*(T)sqrt(10.5)*oneOverSqrtPi*y*((T)4.0*z_pow2 - x_pow2 - y_pow2));
		if (m == 0)
			value = ((T)0.25*(T)sqrt(7.0)*oneOverSqrtPi*z*( (T)2.0*z_pow2 - (T)3.0*x_pow2 - (T)3.0*y_pow2));
		if (m == 1)
			value = ((T)0.25*(T)sqrt(10.5)*oneOverSqrtPi*x*((T)4.0*z_pow2 - x_pow2 - y_pow2));
		if (m == 2)
			value = ((T)0.25*(T)sqrt(105.0)*oneOverSqrtPi*z*(x_pow2 - y_pow2));
		if (m == 3)
			value = ((T)0.25*(T)sqrt(17.5)*oneOverSqrtPi* x*( x_pow2  - (T)3.0*y_pow2));
	}

	if(l == 4) {
		if (m == -4)
			value = ( (T)0.75*(T)sqrt(35.0)*oneOverSqrtPi*x*y*(x_pow2 - y_pow2));
		if (m == -3)
			value = ( (T)0.75*(T)sqrt(17.5)*oneOverSqrtPi*y*z*((T)3.0*x_pow2 - y_pow2));
		if (m == -2)
			value = ( (T)0.75*(T)sqrt(5.0)*oneOverSqrtPi*x*y*((T)7.0*z_pow2 - r_pow2));
		if (m == -1)
			value = ( (T)0.75*(T)sqrt(2.5)*oneOverSqrtPi*y*z*((T)7.0*z_pow2 - (T)3.0*r_pow2));
		if (m == 0)
			value = ( (T)(0.1875)*oneOverSqrtPi*((T)35.0*z_pow4 -(T)30.0*z_pow2*r_pow2 + (T)3.0*r_pow4));
		if (m == 1)
			value = ( (T)0.75*(T)sqrt(2.5)*oneOverSqrtPi*x*z*((T)7.0*z_pow2 - (T)3.0*r_pow2));
		if (m == 2)
			value = ( (T)(0.375)*(T)sqrt(5.0)*oneOverSqrtPi*(x_pow2 - y_pow2)*((T)7.0*z_pow2 - r_pow2));
		if (m == 3)
			value = ( (T)0.75*(T)sqrt(17.5)*oneOverSqrtPi*x*z*(x_pow2 - (T)3.0*y_pow2));
		if (m == 4)
			value = ( (T)(0.1875)*(T)sqrt(35.0)*oneOverSqrtPi*(x_pow2*(x_pow2 - (T)3.0*y_pow2) - y_pow2*((T)3.0*x_pow2 - y_pow2)));
	}
	return value;
}
//END-REAL SPHERICAL HARMONICS OPT------------------------------------------}}}

