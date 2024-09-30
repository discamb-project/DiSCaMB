// vim: noet: sw=3: ts=3
#ifndef CUSTOMLDG_CUH_LF8M4OAH
#define CUSTOMLDG_CUH_LF8M4OAH

#define LDG(x) __ldg(&x)
//#define LDG(x) x

#if (defined(_MSC_VER) && defined(_WIN64)) || defined(__LP64__)
#define __LDG_PTR   "l"
#else
#define __LDG_PTR   "r"
#endif

static __device__ __inline__
cudaComplexType __ldg(const cudaComplexType *ptr) {
	cudaComplexType ret;
	asm volatile ("ld.global.nc.v2.f64 {%0,%1}, [%2];"  : "=d"(ret.real), "=d"(ret.imag) : __LDG_PTR (ptr));
	return ret;
}
static __device__ __inline__
double3 __ldg(const double3 *ptr) {
	double3 ret;
	//asm volatile ("ld.global.nc.v2.f64 {%0,%1}, [%2];"  : "=d"(ret.x), "=d"(ret.y) : __LDG_PTR (ptr));

	double * pomocniczy = (double*)ptr;
	asm volatile ("ld.global.nc.f64 %0, [%1];"  : "=d"(ret.x) : __LDG_PTR (pomocniczy));
	asm volatile ("ld.global.nc.f64 %0, [%1];"  : "=d"(ret.y) : __LDG_PTR (pomocniczy+1));
	asm volatile ("ld.global.nc.f64 %0, [%1];"  : "=d"(ret.z) : __LDG_PTR (pomocniczy+2));
	return ret;
}
static __device__ __inline__
float3 __ldg(const float3 *ptr) {
	float3 ret;

	float * pomocniczy = (float*)ptr;
	asm volatile ("ld.global.nc.f32 %0, [%1];"  : "=f"(ret.x) : __LDG_PTR (pomocniczy));
	asm volatile ("ld.global.nc.f32 %0, [%1];"  : "=f"(ret.y) : __LDG_PTR (pomocniczy+1));
	asm volatile ("ld.global.nc.f32 %0, [%1];"  : "=f"(ret.z) : __LDG_PTR (pomocniczy+2));
	return ret;
}
//TODO: nigdzie z tego nie korzystam (06.08.14 by szmigacz)
//static __device__ __inline__
//wfnCoreParams __ldg(const wfnCoreParams* ptr) {
	//wfnCoreParams ret;

	//unsigned long long * pomocniczy = (unsigned long long *)ptr;
	//asm volatile ("ld.global.nc.u64 %0, [%1];"  : "=l"(ret.coeff) : __LDG_PTR (pomocniczy));
	//asm volatile ("ld.global.nc.u64 %0, [%1];"  : "=l"(ret.exp) : __LDG_PTR (pomocniczy+1));
	//asm volatile ("ld.global.nc.u64 %0, [%1];"  : "=l"(ret.pow) : __LDG_PTR (pomocniczy+2));
	//asm volatile ("ld.global.nc.s32 %0, [%1];"  : "=r"(ret.n) : __LDG_PTR (pomocniczy+3));
	//return ret;

//}
#undef __LDG_PTR

#endif /* end of include guard: CUSTOMLDG_CUH_LF8M4OAH */

