// vim: noet: sw=3: ts=3: foldmethod=marker: foldlevel=0
#ifndef CUDAREDUCTIONS_CUH_6VQDTQDA
#define CUDAREDUCTIONS_CUH_6VQDTQDA

//REDUCE IN WARP------------------------------------------------------------{{{
template<typename T>
__forceinline__ __device__ void reduceInWarp(T & value) {
	for (int i=16; i>=1; i/=2)
		value += __shfl_xor(value, i, 32);
}

template<typename T>
__forceinline__ __device__ void reduceInWarpShared(volatile T* const tab) {
	tab[threadIdx.x] += tab[threadIdx.x+16];
	tab[threadIdx.x] += tab[threadIdx.x+8];
	tab[threadIdx.x] += tab[threadIdx.x+4];
	tab[threadIdx.x] += tab[threadIdx.x+2];
	tab[threadIdx.x] += tab[threadIdx.x+1];
}
//END-REDUCE IN WARP--------------------------------------------------------}}}

//REDUCE IN SHARED (OLD CODE)-----------------------------------------------{{{
//Assumes: size is a power of 2 in range [64,1024]
template<int size, typename T>
__forceinline__ __device__ void reduceArrayAndSave(volatile T * const v, T * const output) {
	if ( size ==1024 ) {
		if(threadIdx.x<512)
			v[threadIdx.x] += v[threadIdx.x+512];
		__syncthreads();
	}
	if ( size >=512 ) {
		if(threadIdx.x<256)
			v[threadIdx.x] += v[threadIdx.x+256];
		__syncthreads();
	}
	if ( size >=256 ) {
		if(threadIdx.x<128)
			v[threadIdx.x] += v[threadIdx.x+128];
		__syncthreads();
	}
	if ( size >=128 ) {
		if(threadIdx.x<64)
			v[threadIdx.x] += v[threadIdx.x+64];
		__syncthreads();
	}
	if ( size >=64) {
		if(threadIdx.x<32) {
			v[threadIdx.x] += v[threadIdx.x+32];
			v[threadIdx.x] += v[threadIdx.x+16];
			v[threadIdx.x] += v[threadIdx.x+8];
			v[threadIdx.x] += v[threadIdx.x+4];
			v[threadIdx.x] += v[threadIdx.x+2];
			v[threadIdx.x] += v[threadIdx.x+1];
			if(threadIdx.x==0)
				*output = v[0];
		}
	}
}

template<int size, typename T>
__forceinline__ __device__ void reduceArray(volatile T * const v) {
	if ( size ==1024 ) {
		if(threadIdx.x<512)
			v[threadIdx.x] += v[threadIdx.x+512];
		__syncthreads();
	}
	if ( size >=512 ) {
		if(threadIdx.x<256)
			v[threadIdx.x] += v[threadIdx.x+256];
		__syncthreads();
	}
	if ( size >=256 ) {
		if(threadIdx.x<128)
			v[threadIdx.x] += v[threadIdx.x+128];
		__syncthreads();
	}
	if ( size >=128 ) {
		if(threadIdx.x<64)
			v[threadIdx.x] += v[threadIdx.x+64];
		__syncthreads();
	}
	if ( size >=64) {
		if(threadIdx.x<32) {
			v[threadIdx.x] += v[threadIdx.x+32];
			v[threadIdx.x] += v[threadIdx.x+16];
			v[threadIdx.x] += v[threadIdx.x+8];
			v[threadIdx.x] += v[threadIdx.x+4];
			v[threadIdx.x] += v[threadIdx.x+2];
			v[threadIdx.x] += v[threadIdx.x+1];
		}
	}
}
//END-REDUCE IN SHARED (OLD CODE)-------------------------------------------}}}

#endif /* end of include guard: CUDAREDUCTIONS_CUH_6VQDTQDA */
