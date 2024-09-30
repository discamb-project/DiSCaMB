#include "HC_SF_GPU.h"

// vim: noet: sw=3: ts=3: foldmethod=marker: foldlevel=0

//INCLUDES------------------------------------------------------------------{{{
#include<cstdio>
#include"globals.h"
#include"fakeComplex.cuh"
#include"types.h"
#include <fstream>
#include<string>
#include<iostream>
#include<vector>
#include"gpu.cuh"
#include<cassert>
#include "cudaConstant.cuh"
#include"memoryAlloc2.cuh"
#include<sys/time.h>
#include<cuda.h>
//END-INCLUDES--------------------------------------------------------------}}}


#define align_pointer(ptr) (void *)(((uintptr_t)ptr+127) & ~ (uintptr_t)0x7F) // align to 128bytes

//SECOND (TIME MEASURE)-----------------------------------------------------{{{
double second(void) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double)tv.tv_sec + (double)tv.tv_usec / 1000000.0;
}
//END-SECOND (TIME MEASURE)-------------------------------------------------}}}

//COPY RESULTS TO CPU-------------------------------------------------------{{{
void klasa::copyResultsToCpu() {
	cudaSafeCall(cudaMemcpy(arrayF, devArrayF, nH_vec*sizeof(arrayF[0]), cudaMemcpyDeviceToHost));

	// kopiowanie wszystkich wynikow o pochodnych
	size_t n_bytes = countPointerDiff(occupancyDerivatives, rotationsRow1);
	cudaSafeCall(cudaMemcpy(occupancyDerivatives, devOccupancyDerivatives, n_bytes, cudaMemcpyDeviceToHost));

	// przeniesione mnozenie pochodnych 3-5 przez 2, z kernela
	for (int i = 0; i < nAtoms; ++i) {
		adpDerivatives3[i].real *= 2, adpDerivatives3[i].imag *= 2;
		adpDerivatives4[i].real *= 2, adpDerivatives4[i].imag *= 2;
		adpDerivatives5[i].real *= 2, adpDerivatives5[i].imag *= 2;
	}
	//cudaSafeCall(cudaMemcpy(adpDerivatives0, devAdpDerivatives0, nAtoms*sizeof(adpDerivatives0[0]), cudaMemcpyDeviceToHost));
	//cudaSafeCall(cudaMemcpy(adpDerivatives1, devAdpDerivatives1, nAtoms*sizeof(adpDerivatives1[0]), cudaMemcpyDeviceToHost));
	//cudaSafeCall(cudaMemcpy(adpDerivatives2, devAdpDerivatives2, nAtoms*sizeof(adpDerivatives2[0]), cudaMemcpyDeviceToHost));
	//cudaSafeCall(cudaMemcpy(adpDerivatives3, devAdpDerivatives3, nAtoms*sizeof(adpDerivatives3[0]), cudaMemcpyDeviceToHost));
	//cudaSafeCall(cudaMemcpy(adpDerivatives4, devAdpDerivatives4, nAtoms*sizeof(adpDerivatives4[0]), cudaMemcpyDeviceToHost));
	//cudaSafeCall(cudaMemcpy(adpDerivatives5, devAdpDerivatives5, nAtoms*sizeof(adpDerivatives5[0]), cudaMemcpyDeviceToHost));
}
//END-COPY RESULTS TO CPU---------------------------------------------------}}}

//RUN GPU-------------------------------------------------------------------{{{
void klasa::runGpu() {
	//lepiej jest dac L1 bo redukcje mam shuflami i nie potrzebuje tego tutaj, a
	//celowo robie spille do L1
	cudaSafeCall(cudaDeviceSetCacheConfig(cudaFuncCachePreferL1));
	//cudaSafeCall(cudaDeviceSetLimit(cudaLimitStackSize, 150*1024));
	//cudaSafeCall(cudaDeviceSetLimit(cudaLimitPrintfFifoSize,2*2048*1024));
	cudaSafeCall(cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte));
	//size_t stacksize;
	//cudaSafeCall(cudaDeviceGetLimit(&stacksize, cudaLimitStackSize));

	cudaSafeCall(cudaMemset((void*)devArrayF, 0, nH_vec * sizeof(devArrayF[0])));

	calculateSFKernel<<<nBlocks, nThreads>>>(
			nH_vec                              ,
			devHVectors                         ,
			devAtomToWfnMap                     ,
			devAtomToTypeMap                    ,
			devAtomic_occupancy                 ,
			devAtomic_multiplicity_factor       ,
			devWfnCore                          ,
			devWfnValence                       ,
			devTypeParams                       ,
			devLocalCoordsRow1                  ,
			devLocalCoordsRow2                  ,
			devLocalCoordsRow3                  ,
			devRotationsRow1                    ,
			devRotationsRow2                    ,
			devRotationsRow3                    ,
			devAtomicPositions                  ,
			devTranslations                     ,
			nAtoms0                             ,
			nAtoms2                             ,
			nAtoms4                             ,
			nSymmOp                             ,
			devAtomicDisplacementParametersSize ,
			devAtomicDisplacementParameters     ,
			devAtomicDisplacementParametersVec1 ,
			devAtomicDisplacementParametersVec2 ,
			devArrayF                           ,
			devWfnParams                        ,
			devDTargetDf                        ,
			devReduceOccupancyDerivatives       ,
			devReducePositionDerivativesX       ,
			devReducePositionDerivativesY       ,
			devReducePositionDerivativesZ       ,
			devReduceAdpDerivatives0            ,
			devReduceAdpDerivatives1            ,
			devReduceAdpDerivatives2            ,
			devReduceAdpDerivatives3            ,
			devReduceAdpDerivatives4            ,
			devReduceAdpDerivatives5
				);

	reduceDerivatives<<<nAtoms,THREADSREDUCEKERNEL>>>(
			devReduceOccupancyDerivatives ,
			devOccupancyDerivatives       ,
			devReducePositionDerivativesX ,
			devPositionDerivativesX       ,
			devReducePositionDerivativesY ,
			devPositionDerivativesY       ,
			devReducePositionDerivativesZ ,
			devPositionDerivativesZ       ,
			devReduceAdpDerivatives0      ,
			devAdpDerivatives0            ,
			devReduceAdpDerivatives1      ,
			devAdpDerivatives1            ,
			devReduceAdpDerivatives2      ,
			devAdpDerivatives2            ,
			devReduceAdpDerivatives3      ,
			devAdpDerivatives3            ,
			devReduceAdpDerivatives4      ,
			devAdpDerivatives4            ,
			devReduceAdpDerivatives5      ,
			devAdpDerivatives5            ,
			(THREADSMAINKERNEL/WARPSIZE)*nBlocks
				);

	//cudaDeviceSynchronize(); //TODO WYWALIC TO!!!!!!
}
//END-RUN GPU---------------------------------------------------------------}}}

//INITIALIZE CONST----------------------------------------------------------{{{
//zapisuje odwrotnosc silni bo tyle potrzebuje na gpu.
void klasa::initializeConst() {
	double tablica[11];
	tablica[0]=1.0;

	for (int i = 1; i < 11; ++i) {
		tablica[i] = i*tablica[i-1];
	}
	for (int i = 1; i < 11; ++i) {
		tablica[i] = 1.0/tablica[i];
	}
	cudaSafeCall(cudaMemcpyToSymbol(constInverseFactorial, tablica, 11*sizeof(double)));
}
//END-INITIALIZE CONST------------------------------------------------------}}}

//FREE MEMORY---------------------------------------------------------------{{{
void klasa::freeHostMemory() {
	for(unsigned int i=0;i<allocatedHostPointers.size();i++)
		cudaSafeCall(cudaFreeHost(allocatedHostPointers[i]));
}

void klasa::freeGpuMemory() {
//    std::cout<< "n alloc pointers = " << allocatedGpuPointers.size() << endl;
	for(unsigned int i=0;i<allocatedGpuPointers.size();i++)
		cudaSafeCall(cudaFree(allocatedGpuPointers[i]));
 //   cout<< __FILE__ << " " << __LINE__ << endl;
}
//END-FREE MEMORY-----------------------------------------------------------}}}

//WRITE OUTPUT--------------------------------------------------------------{{{
#define PRINTDER_OLD(file, pointer) for(int j=0;j<nAtoms;j++) fprintf(file, "%.5e %.5e\n", pointer[j].real, pointer[j].imag);
#define PRINTDER(file, pointer) \
	for(int j=0;j<nAtoms;j++) tmpArr[atomNrOnInput[j]]=pointer[j]; \
PRINTDER_OLD(file, tmpArr);

void klasa::writeOutput(
    std::vector<double> &f_real,
    std::vector<double> &f_imag,
    std::vector<double> &occDerivatives,
    std::vector<double> &xyzDerivatives,
    std::vector<double> &adpDerivatives) {
	//FILE *fp;
	//fp = fopen(outputFileName, "w");
	//for(int i=0;i<nH_vec;i++) {
	//	fprintf(fp,"%.5e %.5e\n",
	//			arrayF[i].real,
	//			arrayF[i].imag);
	//}
	//cudaComplexType tmpArr[nAtoms];
    for(int i=0;i<nH_vec;i++)
    {
        f_real[i] = arrayF[i].real;
        f_imag[i] = arrayF[i].imag;
    }
    
    for(int i=0;i<nAtoms;i++)
    {
        occDerivatives[i] = occupancyDerivatives[i].real;
        xyzDerivatives[3*i] = positionDerivativesX[i].real;
        xyzDerivatives[3*i+1] = positionDerivativesY[i].real;
        xyzDerivatives[3*i+2] = positionDerivativesZ[i].real;
        adpDerivatives[6*i] = adpDerivatives0[i].real;
        adpDerivatives[6*i+1] = adpDerivatives1[i].real;
        adpDerivatives[6*i+2] = adpDerivatives2[i].real;
        adpDerivatives[6*i+3] = adpDerivatives3[i].real;
        adpDerivatives[6*i+4] = adpDerivatives4[i].real;
        adpDerivatives[6*i+5] = adpDerivatives5[i].real;
    }
   
	//fprintf(fp,"occupancy derivatives\n"); PRINTDER(fp, occupancyDerivatives);
	//fprintf(fp,"position derivativesX\n"); PRINTDER(fp, positionDerivativesX);
	//fprintf(fp,"position derivativesY\n"); PRINTDER(fp, positionDerivativesY);
	//fprintf(fp,"position derivativesZ\n"); PRINTDER(fp, positionDerivativesZ);
	//fprintf(fp,"adp derivatives0\n");     PRINTDER(fp, adpDerivatives0);
	//fprintf(fp,"adp derivatives1\n");     PRINTDER(fp, adpDerivatives1);
	//fprintf(fp,"adp derivatives2\n");     PRINTDER(fp, adpDerivatives2);
	//fprintf(fp,"adp derivatives3\n");     PRINTDER(fp, adpDerivatives3);
	//fprintf(fp,"adp derivatives4\n");     PRINTDER(fp, adpDerivatives4);
	//fprintf(fp,"adp derivatives5\n");     PRINTDER(fp, adpDerivatives5);
	//fclose(fp);
    
}
#undef PRINTDER_OLD
#undef PRINTDER
//END-WRITE OUTPUT----------------------------------------------------------}}}

//TIME MACRO----------------------------------------------------------------{{{
#define TIME(x, wypis) do{ \
	double start = second(); \
	x; \
	printf("%f\t" #wypis "\n", second()-start);\
} while(0)
//END-TIME MACRO------------------------------------------------------------}}}


namespace hc_sf_gpu{

    void gpuInfo(
        bool &hasGpu, 
        size_t &memory, 
        std::string &gpuName)
    {
        int nDevices;
        cudaGetDeviceCount(&nDevices);
        if(nDevices==0)
        {
            hasGpu=false;
            return;
        }
        else
            hasGpu=true;
            
        cudaDeviceProp prop;
        cudaGetDeviceProperties(&prop, 0);
        memory = prop.totalGlobalMem;
        gpuName = prop.name;
        
    }
    
    size_t gpuOnlyMemory(
        size_t n_hkl, 
        size_t n_atoms)
    {
        
        size_t n_threads = THREADSMAINKERNEL;
        size_t n_blocks = (n_hkl + n_threads - 1)/n_threads;//~n_hkl/64
        size_t m = (THREADSMAINKERNEL/WARPSIZE);//2
        // 10 = 1(occ) + 3(xyz) + 6(adps)
        // 2 - two components of complex number
        //INCREASESIZE(s, devReduceOccupancyDerivatives, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
        return 10 * 2 * sizeof(real) * m * n_blocks * n_atoms; // (10 * 2 * 8 * 2 /64) * n_hkl * n_atoms = 5 * n_hkl * n_atoms
    }


    void calculate_multipolar_structure_factors_gpu(
        std::stringstream &in,
        const std::vector<double> &dTarget_dF_real,
        const std::vector<double> &dTarget_dF_imag,
        std::vector<double> &structureFactors_real,
        std::vector<double> &structureFactors_imag,
        std::vector<double> &xyzDerivatives,
        std::vector<double> &adpDerivatives,
        std::vector<double> &occupancyDerivatives
    )
    {
        klasa k = klasa();
        //TIME( cudaFree(0), CREATECONTEXT);
        cudaFree(0);
        k.readInputAllocateMemoryOpt(in, dTarget_dF_real, dTarget_dF_imag);
        k.initializeConst();
        //TIME( k.runGpu(), RUNGPU);
        k.runGpu();
        //TIME( k.copyResultsToCpu(), COPYTOCPU);
        k.copyResultsToCpu();
        k.writeOutput(structureFactors_real, structureFactors_imag, occupancyDerivatives, xyzDerivatives, adpDerivatives);
        //TIME( k.freeHostMemory(), FREEHOSTMEM);
        k.freeHostMemory();
        //TIME( k.freeGpuMemory(), FREEGPUMEM);   
        k.freeGpuMemory();   
    }
    
    void calculate_multipolar_structure_factors_gpu2(
    // input
    const std::vector<std::vector<double> > &symmOps,
    const std::vector<std::vector<double> > &core_coeff, 
    const std::vector<std::vector<double> > &core_exp,
    const std::vector<std::vector<size_t> > &core_pow,
    const std::vector<std::vector<double> > &sph_val_coeff, 
    const std::vector<std::vector<double> > &sph_val_exp,
    const std::vector<std::vector<size_t> > &sph_val_pow,
    const std::vector<double> &def_val_exponent,
    const std::vector<std::vector<size_t> > &def_val_n_l,
    const std::vector<double> &p_val,
    const std::vector<double> &kappa_spherical,
    const std::vector<double> &kappa_def_valence,
    const std::vector<std::vector<std::vector<double> > > &p_lm,
    const std::vector<int> &atomToWfnTypeMap,
    const std::vector<int> &atomToAtomTypeMap,
    const std::vector<std::vector<double> > &atomicPositions,
    const std::vector<std::vector<double> > &adps,
    const std::vector<double> &occupancy,
    const std::vector<double> &multiplicity_weigths,
    const std::vector<std::vector<std::vector<double> > > &lcs,
    const std::vector<vector<double> > &hkl,
    const std::vector<double> &dTarget_dF_real,
    const std::vector<double> &dTarget_dF_imag,
    //output
    std::vector<double> &structureFactors_real,
    std::vector<double> &structureFactors_imag,
    std::vector<double> &xyzDerivatives,
    std::vector<double> &adpDerivatives,
    std::vector<double> &occupancyDerivatives)
    {
        klasa k = klasa();
        TIME( cudaFree(0), CREATECONTEXT);
        //k.readInputAllocateMemoryOpt(in, dTarget_dF_real, dTarget_dF_imag);
        
        k.passInputDataAllocateMemoryOpt(
            symmOps,
            core_coeff, 
            core_exp,
            core_pow,
            sph_val_coeff, 
            sph_val_exp,
            sph_val_pow,
            def_val_exponent,
            def_val_n_l,
            p_val,
            kappa_spherical,
            kappa_def_valence,
            p_lm,
            atomToWfnTypeMap,
            atomToAtomTypeMap,
            atomicPositions,
            adps,
            occupancy,
            multiplicity_weigths,
            lcs,
            hkl,
            dTarget_dF_real,
            dTarget_dF_imag);

        k.initializeConst();
        
        TIME( k.runGpu(), RUNGPU);
        
        TIME( k.copyResultsToCpu(), COPYTOCPU);
        
        k.writeOutput(structureFactors_real, structureFactors_imag, occupancyDerivatives, xyzDerivatives, adpDerivatives);
        
        TIME( k.freeHostMemory(), FREEHOSTMEM);
        TIME( k.freeGpuMemory(), FREEGPUMEM);       
    }

    
}


#undef TIME

