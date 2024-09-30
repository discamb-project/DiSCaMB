// vim: noet: sw=3: ts=3
#ifndef FAKECOMPLEX_CUH_0SYHM1GK
#define FAKECOMPLEX_CUH_0SYHM1GK

#ifdef __CUDACC__
#define ALIGN(x)  __align__(x)
#else
#ifdef __GNUC__
#define ALIGN(x)  __attribute__ ((aligned (x)))
#endif
#endif

#define CUDA_CALLABLE __host__ __device__

template<typename T, int alignment>
class ALIGN(alignment) cudaComplex {
	public:
		T real;
		T imag;

		CUDA_CALLABLE cudaComplex(T _real, T _imag) : real(_real), imag(_imag){}
		/*CUDA_CALLABLE cudaComplex() : real((T)0), imag((T)0){}*/
		CUDA_CALLABLE cudaComplex(){}

		CUDA_CALLABLE __forceinline__ cudaComplex& operator+=(const cudaComplex &rhs) {
			this->real += rhs.real;
			this->imag += rhs.imag;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ volatile cudaComplex& operator+=(const volatile cudaComplex &rhs) volatile {
			this->real += rhs.real;
			this->imag += rhs.imag;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ volatile cudaComplex& operator=(const volatile cudaComplex &rhs) volatile {
			this->real = rhs.real;
			this->imag = rhs.imag;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ cudaComplex& operator=(const cudaComplex &rhs) {
			this->real = rhs.real;
			this->imag = rhs.imag;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ volatile cudaComplex& operator=(const cudaComplex &rhs) volatile{
			this->real = rhs.real;
			this->imag = rhs.imag;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ cudaComplex& operator-=(const cudaComplex &rhs) {
			this->real -= rhs.real;
			this->imag -= rhs.imag;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ cudaComplex& operator*=(const cudaComplex &rhs) {
			T tempReal = this->real;
			T tempImag = this->imag;
			this->real = tempReal*rhs.real - tempImag*rhs.imag;
			this->imag = tempReal*rhs.imag + tempImag*rhs.real;
			return *this;
		}
		CUDA_CALLABLE __forceinline__ cudaComplex& operator*=(const T &rhs) {
			this->real *= rhs;
			this->imag *= rhs;
			return *this;
		}

		CUDA_CALLABLE __forceinline__ const cudaComplex operator*(const T &other) const {
			cudaComplex result = *this;
			result *= other;
			return result;
		}

		CUDA_CALLABLE __forceinline__ const cudaComplex operator+(const cudaComplex &other) const {
			cudaComplex result = *this;
			result += other;
			return result;
		}
		CUDA_CALLABLE __forceinline__ const cudaComplex operator-(const cudaComplex &other) const {
			cudaComplex result = *this;
			result -= other;
			return result;
		}
		CUDA_CALLABLE __forceinline__ const cudaComplex operator*(const cudaComplex &other) const {
			cudaComplex result = *this;
			result *= other;
			return result;
		}
};

typedef cudaComplex<double, 16> cudaComplexDouble;
typedef cudaComplex<float, 8>   cudaComplexFloat;

CUDA_CALLABLE __forceinline__ float fabsf(const cudaComplexFloat &arg) {
	return sqrtf( arg.real*arg.real + arg.imag*arg.imag);
}
CUDA_CALLABLE __forceinline__ double fabs(const cudaComplexDouble &arg) {
	return sqrt( arg.real*arg.real + arg.imag*arg.imag);
}
CUDA_CALLABLE __forceinline__ cudaComplexFloat exp(const cudaComplexFloat &arg) {
	float absolute = expf(arg.real);
	return cudaComplexFloat( absolute*cosf(arg.imag), absolute*sinf(arg.imag));
}
CUDA_CALLABLE __forceinline__ cudaComplexDouble exp(const cudaComplexDouble &arg) {
	double absolute = exp(arg.real);
	return cudaComplexDouble( absolute*cos(arg.imag), absolute*sin(arg.imag));
}

#undef CUDA_CALLABLE
#endif /* end of include guard: FAKECOMPLEX_CUH_0SYHM1GK */
