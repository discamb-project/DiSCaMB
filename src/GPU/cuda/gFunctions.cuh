// vim: noet: sw=3: ts=3: foldmethod=marker: foldlevel=0
#include<cmath>
#include<cassert>

//G FUNCTION OPT------------------------------------------------------------{{{
template<int l, typename T>
__host__ __device__ __forceinline__ T gFunctionOpt(
		const int n,
		T const h,
		T const Z) 
{
	const T K = (T)(2.0*M_PI)*h; // K and Z symbols like in Coppens book

	const T K_pow2 = K*K;
	const T K_pow3 = K_pow2 * K;
	const T K_pow4 = K_pow2 * K_pow2;
	const T K_pow5 = K_pow2 * K_pow3;
	const T K_pow6 = K_pow4 * K_pow2;
	const T K_pow7 = K_pow6 * K;
	const T K_pow8 = K_pow4 * K_pow4;
	const T K_pow9 = K_pow8 * K;

	const T Z_pow2 = Z*Z;
	const T Z_pow3 = Z_pow2 * Z;
	const T Z_pow4 = Z_pow2 * Z_pow2;
	const T Z_pow5 = Z_pow2 * Z_pow3;
	const T Z_pow6 = Z_pow4 * Z_pow2;
	const T Z_pow7 = Z_pow6 * Z;
	const T Z_pow8 = Z_pow4 * Z_pow4;
	const T Z_pow9 = Z_pow8 * Z;

	const T d = K_pow2 + Z_pow2;

	const T d_inv = (T)1.0/d;
	const T d_inv_pow2 = d_inv * d_inv;
	const T d_inv_pow4 = d_inv_pow2 * d_inv_pow2;
	const T d_inv_pow8 = d_inv_pow4 * d_inv_pow4;
	const T d_inv_pow3 = d_inv_pow2 * d_inv;
	const T d_inv_pow5 = d_inv_pow4 * d_inv;
	const T d_inv_pow6 = d_inv_pow4 * d_inv_pow2;
	const T d_inv_pow7 = d_inv_pow6 * d_inv;
	const T d_inv_pow9 = d_inv_pow8 * d_inv;
	const T d_inv_pow10 = d_inv_pow4 * d_inv_pow6;

	T value = (T)0;

	if(l==0){
		if(n==2)
			value = 2*Z*d_inv_pow2;
		else if(n==3)
			value = 2*(3*Z_pow2 - K_pow2)*d_inv_pow3;
		else if(n==4)
			value = 24*Z*(Z_pow2 - K_pow2)*d_inv_pow4;
		else if(n==5)
			value = 24*(5*Z_pow4 - 10*K_pow2*Z_pow2 + K_pow4)*d_inv_pow5;//poprawionywspolczynnk!
		else if(n==6)
			value = 240*Z*(K_pow2 - 3*Z_pow2)*(3*K_pow2 - Z_pow2)*d_inv_pow6;
		else if(n==7)
			value = 720*(7*Z_pow6 - 35*K_pow2*Z_pow4 + 21*K_pow4*Z_pow2 - K_pow6)*d_inv_pow7;
		else if(n==8)
			value = 40320*(Z_pow7 - 7*K_pow2*Z_pow5 + 7*K_pow4*Z_pow3 - K_pow6*Z)*d_inv_pow8;
		else if(n==9)
			value = (362880*Z_pow8 - 3386880*K_pow2*Z_pow6 + 5080320*K_pow4*Z_pow4 - 1451520*K_pow6*Z_pow2 + 40320*K_pow8)*d_inv_pow9;
		else if(n==10)
			value = (3628800*Z_pow9 - 43545600*K_pow2*Z_pow7 + 91445760*K_pow4*Z_pow5 - 43545600*K_pow6*Z_pow3 + 3628800*K_pow8*Z)*d_inv_pow10;
	}

	if(l==1){
		if(n==3)
			value = 8*K*Z*d_inv_pow3;
		else if(n==4)
			value = 8*K*(5*Z_pow2 - K_pow2)*d_inv_pow4;
		else if(n==5)
			value = 48*K*Z*(5*Z_pow2 - 3*K_pow2)*d_inv_pow5;
		else if(n==6)
			value = 48*K*(35*Z_pow4 - 42*K_pow2*Z_pow2 + 3*K_pow4)*d_inv_pow6;
		else if(n==7)
			value = 1920*K*Z*(7*Z_pow4 - 14*K_pow2*Z_pow2 + 3*K_pow4)*d_inv_pow7;
		else if(n==8)
			value = 5760*K*(21*Z_pow6 - 63*K_pow2*Z_pow4 + 27*K_pow4*Z_pow2 - K_pow6)*d_inv_pow8;
		else if(n==9)
			value = (1209600*K*Z_pow7 - 5080320*K_pow3*Z_pow5 + 3628800*K_pow5*Z_pow3 - 403200*K_pow7*Z)*d_inv_pow9;
		else if(n==10)
			value = (13305600*K*Z_pow8 - 74511360*K_pow3*Z_pow6 + 79833600*K_pow5*Z_pow4 - 17740800*K_pow7*Z_pow2 + 403200*K_pow9)*d_inv_pow10;

	}

	if(l==2){
		if(n==4)
			value = 48*K_pow2*Z*d_inv_pow4;
		else if(n==5)
			value = 48*K_pow2*(7*Z_pow2 - K_pow2)*d_inv_pow5;
		else if(n==6)
			value = 384*K_pow2*Z*(7*Z_pow2 - 3*K_pow2)*d_inv_pow6;
		else if(n==7)
			value = 1152*K_pow2*(21*Z_pow4 - 18*K_pow2*Z_pow2 + K_pow4)*d_inv_pow7;
		else if(n==8)
			value = 11520*K_pow2*Z*(21*Z_pow4 - 30*K_pow2*Z_pow2 + 5*K_pow4)*d_inv_pow8;
		else if(n==9)
			value = (2661120*K_pow2*Z_pow6 - 5702400*K_pow4*Z_pow4 + 1900800*K_pow6*Z_pow2 - 57600*K_pow8)*d_inv_pow9;
		else if(n==10)
			value = (31933440*K_pow2*Z_pow7 - 95800320*K_pow4*Z_pow5 + 53222400*K_pow6*Z_pow3 - 4838400*K_pow8*Z)*d_inv_pow10;
	}
	if(l==3){
		if(n==5)
			value = (384*K_pow3*Z)*d_inv_pow5;
		else if(n==6)
			value = (3456*K_pow3*Z_pow2 - 384*K_pow5)*d_inv_pow6;
		else if(n==7)
			value = (34560*K_pow3*Z_pow3 - 11520*K_pow5*Z)*d_inv_pow7;
		else if(n==8)
			value = (380160*K_pow3*Z_pow4 - 253440*K_pow5*Z_pow2 + 11520*K_pow7)*d_inv_pow8;
		else if(n==9)
			value = (4561920*K_pow3*Z_pow5 - 5068800*K_pow5*Z_pow3 + 691200*K_pow7*Z)*d_inv_pow9;
		else if(n==10)
			value = (59304960*K_pow3*Z_pow6 - 98841600*K_pow5*Z_pow4 + 26956800*K_pow7*Z_pow2 - 691200*K_pow9)*d_inv_pow10;
	}
	if(l==4){
		if(n==6)
			value = (3840*K_pow4*Z)*d_inv_pow6;
		else if(n==7)
			value = (42240*K_pow4*Z_pow2 - 3840*K_pow6)*d_inv_pow7;
		else if(n==8)
			value = (506880*K_pow4*Z_pow3 - 138240*K_pow6*Z)*d_inv_pow8;
		else if(n==9)
			value = (6589440*K_pow4*Z_pow4 - 3594240*K_pow6*Z_pow2 + 138240*K_pow8)*d_inv_pow9;
		else if(n==10)
			value = (92252160*K_pow4*Z_pow5 - 83865600*K_pow6*Z_pow3 + 9676800*K_pow8*Z)*d_inv_pow10;
	}

	return value;
}
//END-G FUNCTION OPT--------------------------------------------------------}}}
