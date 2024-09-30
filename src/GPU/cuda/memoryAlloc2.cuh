// vim: noet: sw=3: ts=3: foldmethod=marker: foldlevel=0
//INCLUDES------------------------------------------------------------------{{{
#include<cstdio>
#include"fakeComplex.cuh"
#include"types.h"
#include <fstream>
#include<string>
#include<iostream>
#include<vector>
#include<cassert>
#include<stdint.h>
#include "cudaConstant.cuh"
//END-INCLUDES--------------------------------------------------------------}}}

using namespace std;

//MACROS--------------------------------------------------------------------{{{
#define INCREASESIZE(size, pointer, amount) \
	size += amount * sizeof(pointer[0]); \
size = (size + 127) & ~ 0x7F;

#define align_pointer(ptr) (void *)(((uintptr_t)ptr+127) & ~ (uintptr_t)0x7F) // align to 128bytes

#define PINPOINTER(pointer, amount) \
	pointer = (typeof(pointer)) curPtr; \
curPtr = (void *)(((typeof(pointer))curPtr) + amount); \
curPtr = (void *)(((uintptr_t)curPtr+127) & ~ (uintptr_t)0x7F);

#define ALLOCATEHOSTPOINTER(pointer,amount) do { \
	cudaSafeCall(cudaMallocHost((void**)&pointer, amount*sizeof(pointer[0]))); \
	allocatedHostPointers.push_back((void*)pointer); \
} while(0)

#define ALLOCATEGPUPOINTER(pointer, amount) do{ \
	cudaSafeCall(cudaMalloc((void**)&(pointer), amount*sizeof((pointer)[0]))); \
	allocatedGpuPointers.push_back((void*)pointer);\
} while(0)
//END-MACROS----------------------------------------------------------------}}}


size_t countPointerDiff(void * a, void * b) {
	return size_t((uintptr_t)b - (uintptr_t)a);
}

//STRUCT ATOM---------------------------------------------------------------{{{
struct Atom { // pomocniczy struct do odpowiedniego sortowania atomow
	int   atomNrOnInput;
	int   atomToWfnMap;
	int   atomToTypeMap;
	real  atomic_occupancy;
	real  atomic_multiplicity_factor;
	real3 atomicPositions;
	int   atomicDisplacementParametersSize;
	real3 atomicDisplacementParametersVec1;
	real3 atomicDisplacementParametersVec2;
	real  atomicDisplacementParameters;
	real3 localCoordsRow1;
	real3 localCoordsRow2;
	real3 localCoordsRow3;
	real3 rotationsRow1;
	real3 rotationsRow2;
	real3 rotationsRow3;
	real3 translations;
};
//END-STRUCT ATOM-----------------------------------------------------------}}}

//COUNT COMMON MEMORY-------------------------------------------------------{{{
size_t klasa::countCommonMemory() {
	size_t s = 0;

	INCREASESIZE(s, dTargetDf, nH_vec);
	INCREASESIZE(s, occupancyDerivatives, nAtoms);
	INCREASESIZE(s, positionDerivativesX, nAtoms);
	INCREASESIZE(s, positionDerivativesY, nAtoms);
	INCREASESIZE(s, positionDerivativesZ, nAtoms);
	for (int i = 0; i < N_ADP_DERIVATIVES; ++i) { // pochodne od 0 do 5
		INCREASESIZE(s, adpDerivatives0, nAtoms);
	}

	INCREASESIZE(s, rotationsRow1, nSymmOp);
	INCREASESIZE(s, rotationsRow2, nSymmOp);
	INCREASESIZE(s, rotationsRow3, nSymmOp);
	INCREASESIZE(s, translations,  nSymmOp);

	INCREASESIZE(s, wfnCore               , nWfnParams);
	INCREASESIZE(s, wfnValence            , nWfnParams);
	INCREASESIZE(s, wfnParams             , nWfnParams);

	INCREASESIZE(s, typeParams, nTypes);

	INCREASESIZE(s, atomicDisplacementParameters    , nAtoms);
	INCREASESIZE(s, atomicDisplacementParametersSize, nAtoms);
	INCREASESIZE(s, atomicDisplacementParametersVec1, nAtoms);
	INCREASESIZE(s, atomicDisplacementParametersVec2, nAtoms);

	INCREASESIZE(s, atomicPositions                 , nAtoms);
	INCREASESIZE(s, atomic_occupancy                , nAtoms);
	INCREASESIZE(s, atomToTypeMap                   , nAtoms);
	INCREASESIZE(s, atomToWfnMap                    , nAtoms);
	INCREASESIZE(s, atomic_multiplicity_factor      , nAtoms);

	INCREASESIZE(s, localCoordsRow1                 , nAtoms);
	INCREASESIZE(s, localCoordsRow2                 , nAtoms);
	INCREASESIZE(s, localCoordsRow3                 , nAtoms);

	INCREASESIZE(s, hVectors, nH_vec);
	INCREASESIZE(s, arrayF  , nH_vec);

	return (s+127) & ~0x7F; // trzeba zwrocic wyrownany wynik
}
//END-COUNT COMMON MEMORY---------------------------------------------------}}}

//COUNT HOST ONLY MEMORY----------------------------------------------------{{{
size_t klasa::countHostonlyMemory() {
	size_t s = 0;
	INCREASESIZE(s, atomNrOnInput, nAtoms);
	INCREASESIZE(s, hostConstDefValenceExp, nWfnParams);
	INCREASESIZE(s, hostConstDefValencePow, nWfnParams);
	INCREASESIZE(s, hostConstPLM, nTypes * SPHHARMCOEF);
	return s;
}
//END-COUNT HOST ONLY MEMORY------------------------------------------------}}}

//COUNT GPU ONLY MEMORY-----------------------------------------------------{{{
size_t klasa::countGPUonlyMemory() {
	size_t s = 0;
	INCREASESIZE(s, devReduceOccupancyDerivatives, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	INCREASESIZE(s, devReducePositionDerivativesX, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	INCREASESIZE(s, devReducePositionDerivativesY, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	INCREASESIZE(s, devReducePositionDerivativesZ, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	for (int i = 0; i < N_ADP_DERIVATIVES; ++i) {
		INCREASESIZE(s, devReduceAdpDerivatives0,(THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	}
	return s;
}
//END-COUNT GPU ONLY MEMORY-------------------------------------------------}}}

//PIN HOST POINTERS TO UNI MEMORY-------------------------------------------{{{
void klasa::pinHostPointersToUnitedMemory() {
	void * curPtr = unitedHostMemory;

	PINPOINTER(dTargetDf, nH_vec);
	PINPOINTER(occupancyDerivatives,  nAtoms);
	PINPOINTER(positionDerivativesX,  nAtoms);
	PINPOINTER(positionDerivativesY,  nAtoms);
	PINPOINTER(positionDerivativesZ,  nAtoms);
	PINPOINTER(adpDerivatives0, nAtoms);
	PINPOINTER(adpDerivatives1, nAtoms);
	PINPOINTER(adpDerivatives2, nAtoms);
	PINPOINTER(adpDerivatives3, nAtoms);
	PINPOINTER(adpDerivatives4, nAtoms);
	PINPOINTER(adpDerivatives5, nAtoms);

	PINPOINTER(rotationsRow1, nSymmOp);
	PINPOINTER(rotationsRow2, nSymmOp);
	PINPOINTER(rotationsRow3, nSymmOp);
	PINPOINTER(translations,  nSymmOp);

	PINPOINTER(wfnCore               , nWfnParams);
	PINPOINTER(wfnValence            , nWfnParams);
	PINPOINTER(wfnParams             , nWfnParams);

	PINPOINTER(typeParams, nTypes);

	PINPOINTER(atomNrOnInput                   , nAtoms);
	PINPOINTER(atomicDisplacementParameters    , nAtoms);
	PINPOINTER(atomicDisplacementParametersSize, nAtoms);
	PINPOINTER(atomicDisplacementParametersVec1, nAtoms);
	PINPOINTER(atomicDisplacementParametersVec2, nAtoms);

	PINPOINTER(atomicPositions                 , nAtoms);
	PINPOINTER(atomic_occupancy                , nAtoms);
	PINPOINTER(atomToTypeMap                   , nAtoms);
	PINPOINTER(atomToWfnMap                    , nAtoms);
	PINPOINTER(atomic_multiplicity_factor      , nAtoms);

	PINPOINTER(localCoordsRow1                 , nAtoms);
	PINPOINTER(localCoordsRow2                 , nAtoms);
	PINPOINTER(localCoordsRow3                 , nAtoms);

	PINPOINTER(hVectors, nH_vec);
	PINPOINTER(arrayF  , nH_vec);

	PINPOINTER(hostConstDefValenceExp, nWfnParams);
	PINPOINTER(hostConstDefValencePow, nWfnParams);
	PINPOINTER(hostConstPLM, nTypes * SPHHARMCOEF);
}
//END-PIN HOST POINTERS TO UNI MEMORY---------------------------------------}}}

//PIN GPU POINTERS TO UNI MEMORY--------------------------------------------{{{
void klasa::pinGPUPointersToUnitedMemory() {
	void * curPtr = unitedGPUMemory;

	PINPOINTER(devDTargetDf, nH_vec);
	PINPOINTER(devOccupancyDerivatives, nAtoms);
	PINPOINTER(devPositionDerivativesX, nAtoms);
	PINPOINTER(devPositionDerivativesY, nAtoms);
	PINPOINTER(devPositionDerivativesZ, nAtoms);
	PINPOINTER(devAdpDerivatives0, nAtoms);
	PINPOINTER(devAdpDerivatives1, nAtoms);
	PINPOINTER(devAdpDerivatives2, nAtoms);
	PINPOINTER(devAdpDerivatives3, nAtoms);
	PINPOINTER(devAdpDerivatives4, nAtoms);
	PINPOINTER(devAdpDerivatives5, nAtoms);

	PINPOINTER(devReduceOccupancyDerivatives, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReducePositionDerivativesX, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReducePositionDerivativesY, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReducePositionDerivativesZ, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReduceAdpDerivatives0, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReduceAdpDerivatives1, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReduceAdpDerivatives2, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReduceAdpDerivatives3, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReduceAdpDerivatives4, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);
	PINPOINTER(devReduceAdpDerivatives5, (THREADSMAINKERNEL/WARPSIZE)*nBlocks*nAtoms);


	PINPOINTER(devRotationsRow1, nSymmOp);
	PINPOINTER(devRotationsRow2, nSymmOp);
	PINPOINTER(devRotationsRow3, nSymmOp);
	PINPOINTER(devTranslations , nSymmOp);

	PINPOINTER(devWfnCore   , nWfnParams);
	PINPOINTER(devWfnValence, nWfnParams);
	PINPOINTER(devWfnParams , nWfnParams);

	PINPOINTER(devTypeParams, nTypes);

	PINPOINTER(devAtomicDisplacementParameters    , nAtoms);
	PINPOINTER(devAtomicDisplacementParametersSize, nAtoms);
	PINPOINTER(devAtomicDisplacementParametersVec1, nAtoms);
	PINPOINTER(devAtomicDisplacementParametersVec2, nAtoms);

	PINPOINTER(devAtomicPositions           , nAtoms);
	PINPOINTER(devAtomic_occupancy          , nAtoms);
	PINPOINTER(devAtomToTypeMap             , nAtoms);
	PINPOINTER(devAtomToWfnMap              , nAtoms);
	PINPOINTER(devAtomic_multiplicity_factor, nAtoms);

	PINPOINTER(devLocalCoordsRow1, nAtoms);
	PINPOINTER(devLocalCoordsRow2, nAtoms);
	PINPOINTER(devLocalCoordsRow3, nAtoms);

	PINPOINTER(devHVectors, nH_vec);
	PINPOINTER(devArrayF  , nH_vec);
}
//END-PIN GPU POINTERS TO UNI MEMORY----------------------------------------}}}

//READ INPUT ALLOCATE MEMORY------------------------------------------------{{{

void klasa::readInputAllocateMemoryOpt(
    std::stringstream &in,
    const std::vector<double> &dTarget_dF_real,
    const std::vector<double> &dTarget_dF_imag)
{
	//printf("czytam dane\n");
	string auxString;
	//ifstream in(inputFileName);
	//if(!in.good()) {
	//	cerr<< "can not read input file " << inputFileName << endl;
	//	exit(1);
	//}
	in>> auxString >> nSymmOp;
	in>> auxString >> nWfnParams;
	in>> auxString >> nTypes;
	in>> auxString >> nAtoms;
	in>> auxString >> nH_vec;

	if (nTypes > MAXATOMTYPES) {
		printf("za duzo typow atomow!!!!\n");
		exit(0);
	}

	if (nWfnParams > MAXWFNPARAMS ) {
		printf("za duzo wfn paramsow\n");
		exit(0);
	}

	if (nSymmOp > MAXSYMMETRYOPERATIONS) {
		printf("za duzo symetrii\n");
		exit(0);
	}

	in>> auxString >> auxString; // "symmetry operations"

	nThreads= THREADSMAINKERNEL;
	nBlocks = (nH_vec + nThreads - 1)/nThreads;

	size_t memorySizeHost = 0;
	size_t memorySizeGPU  = 0;

	// liczymy, ile musimy zaalokowac
	memorySizeHost += countCommonMemory();
	memorySizeGPU  += memorySizeHost;

	memorySizeHost += countHostonlyMemory();
	memorySizeGPU  += countGPUonlyMemory();

	// malloc duzego bloku pamieci
	cudaSafeCall(cudaMallocHost((void**)&unitedHostMemory, memorySizeHost));
	allocatedHostPointers.push_back(unitedHostMemory);
	cudaSafeCall(cudaMalloc((void**)&unitedGPUMemory, memorySizeGPU));
	allocatedGpuPointers.push_back(unitedGPUMemory);

	// i podpiecie odpowiednich wskaznikow
	pinHostPointersToUnitedMemory();
	pinGPUPointersToUnitedMemory();

	// zerowanie sporej ilosci tablic wynikowych
	size_t n_bytes = countPointerDiff(unitedGPUMemory, devRotationsRow1);
	cudaSafeCall(cudaMemset(unitedGPUMemory, 0, n_bytes));

	for(int i=0;i<nSymmOp;i++) {
		in >> rotationsRow1[i].x  >> rotationsRow1[i].y >> rotationsRow1[i].z;
		in >> rotationsRow2[i].x  >> rotationsRow2[i].y >> rotationsRow2[i].z;
		in >> rotationsRow3[i].x  >> rotationsRow3[i].y >> rotationsRow3[i].z;
		in >> translations[i].x   >> translations[i].y  >> translations[i].z;
	}
	// (chwilowo) z tego nie korzystamy
	//cudaSafeCall(cudaMemcpyToSymbol(constRotationsRow1 , rotationsRow1 , nSymmOp*sizeof(rotationsRow1[0])));
	//cudaSafeCall(cudaMemcpyToSymbol(constRotationsRow2 , rotationsRow2 , nSymmOp*sizeof(rotationsRow1[0])));
	//cudaSafeCall(cudaMemcpyToSymbol(constRotationsRow3 , rotationsRow3 , nSymmOp*sizeof(rotationsRow1[0])));
	//cudaSafeCall(cudaMemcpyToSymbol(constTranslations  , translations  , nSymmOp*sizeof(translations[0])));

	// memcpy tablic rotationsRow1,2,3 i translations za pomoca 1 operacji
	n_bytes = countPointerDiff((void *)rotationsRow1, (void *)(translations+nSymmOp));
	cudaSafeCall(cudaMemcpy(devRotationsRow1 , rotationsRow1 , n_bytes , cudaMemcpyHostToDevice));

	in>> auxString >> auxString >> auxString; // "wfn parameter sets "

	for(int i=0;i<nWfnParams;i++) {
		in >> wfnCore[i].n;
		ALLOCATEHOSTPOINTER(wfnCore[i].coeff , wfnCore[i].n);
		ALLOCATEHOSTPOINTER(wfnCore[i].exp   , wfnCore[i].n);
		ALLOCATEHOSTPOINTER(wfnCore[i].pow   , wfnCore[i].n);

		for(int j=0;j<wfnCore[i].n;j++) {
			in >> wfnCore[i].coeff[j] >> wfnCore[i].pow[j] >> wfnCore[i].exp[j];
			wfnCore[i].pow[j] += 2;
		}



		// spherical valence
		in >> wfnValence[i].n;
		ALLOCATEHOSTPOINTER(wfnValence[i].coeff , wfnValence[i].n);
		ALLOCATEHOSTPOINTER(wfnValence[i].exp   , wfnValence[i].n);
		ALLOCATEHOSTPOINTER(wfnValence[i].pow   , wfnValence[i].n);

		for(int j=0;j<wfnValence[i].n;j++) {
			in >> wfnValence[i].coeff[j] >> wfnValence[i].pow[j] >> wfnValence[i].exp[j];
			wfnValence[i].pow[j] += 2;
		}

		in >> wfnParams[i].defValenceExp;
		in >> wfnParams[i].n;

		hostConstDefValenceExp[i] = wfnParams[i].defValenceExp;

		ALLOCATEHOSTPOINTER(wfnParams[i].defValencePow, wfnParams[i].n);

		//printf("ZONK %d\n", wfnParams[i].n);
		//w kodzie zawsze brane jest +2
		for(int j=0;j<wfnParams[i].n;j++) {
			int tempPow;
			in>> tempPow;
			wfnParams[i].defValencePow[j] = tempPow + 2;
		}
	}
	cudaSafeCall(cudaMemcpyToSymbol(constDefValenceExp, hostConstDefValenceExp, nWfnParams*sizeof(hostConstDefValenceExp[0])));


	for(int i=0;i<nWfnParams;i++) {
		for(int j=0;j<wfnParams[i].n;j++) {
			hostConstDefValencePow[i*5+j] = wfnParams[i].defValencePow[j];
		}
	}
	cudaSafeCall(cudaMemcpyToSymbol(constDefValencePow, hostConstDefValencePow, nWfnParams*5*sizeof(int)));


	// memcpy tablic wfnCore/Valence/Params za pomoca 1 operacji
	n_bytes = countPointerDiff((void *)wfnCore, (void *)(wfnParams+nWfnParams));
	cudaSafeCall(cudaMemcpy(devWfnCore, wfnCore, n_bytes, cudaMemcpyHostToDevice));


	in>> auxString >> auxString >> auxString >> auxString; // "atom type parameter sets"


	for(int i=0;i<nTypes;i++) {
		real kappaDefValence;
		real kappaSpherical;
		in >> typeParams[i].pVal >> kappaSpherical;
		in >> kappaDefValence >> typeParams[i].maxL;

		typeParams[i].invKappaDefValence = 1.0/kappaDefValence;
		typeParams[i].invKappaSpherical = 1.0/kappaSpherical;

		ALLOCATEHOSTPOINTER(typeParams[i].pLM, (typeParams[i].maxL+1)*(2*(typeParams[i].maxL+1)+1));

		//printf("MAX L %d\n", typeParams[i].maxL);
		int counter = 0 ;
		for(int j=0;j<=typeParams[i].maxL;j++) {
			for(int k=0;k<2*j+1;k++) {
				in>>typeParams[i].pLM[j*(2*typeParams[i].maxL+1)+k];
				hostConstPLM[i*SPHHARMCOEF+counter] = typeParams[i].pLM[j*(2*typeParams[i].maxL+1)+k]; //TODO: ZONK dla maxL roznego od 4 (20.07.14 by szmigacz)
				counter++;
			}
		}

	}

	cudaSafeCall(cudaMemcpyToSymbol(constPLM, hostConstPLM, nTypes*SPHHARMCOEF*sizeof(hostConstPLM[0])));

	cudaSafeCall(cudaMemcpy(devTypeParams, typeParams, nTypes*sizeof(typeParams[0]), cudaMemcpyHostToDevice));

	in>> auxString >> auxString; // "atomic data"

	std::vector<Atom> atoms[MAX_L]; // koszyki do sortowania po liczbie harmonijek
	for(int i=0;i<nAtoms;i++) {
		Atom a;
		a.atomNrOnInput = i;
		in >> a.atomToWfnMap >> a.atomToTypeMap;
		//printf("zonk %d\n",atomToWfnMap[i]);

		in >> a.atomicPositions.x >> a.atomicPositions.y >> a.atomicPositions.z;
		in >> a.atomic_occupancy  >> a.atomic_multiplicity_factor;

		in >> a.atomicDisplacementParametersSize;

		if ( a.atomicDisplacementParametersSize == 1)
			in >> a.atomicDisplacementParameters;
		else {
			in >> a.atomicDisplacementParametersVec1.x >> a.atomicDisplacementParametersVec1.y >> a.atomicDisplacementParametersVec1.z;
			in >> a.atomicDisplacementParametersVec2.x >> a.atomicDisplacementParametersVec2.y >> a.atomicDisplacementParametersVec2.z;

			a.atomicDisplacementParametersVec2.x *= 2;
			a.atomicDisplacementParametersVec2.y *= 2;
			a.atomicDisplacementParametersVec2.z *= 2;
		}

		in >> a.localCoordsRow1.x >> a.localCoordsRow1.y >> a.localCoordsRow1.z;
		in >> a.localCoordsRow2.x >> a.localCoordsRow2.y >> a.localCoordsRow2.z;
		in >> a.localCoordsRow3.x >> a.localCoordsRow3.y >> a.localCoordsRow3.z;

		atoms[wfnParams[a.atomToWfnMap].n].push_back(a); // dodanie do odpowiedniego koszyka
	}
	// kopiowanie posortowanych danych o atomach
	for (int l = 0, i = 0; l < MAX_L; ++l) {
		for (unsigned int j = 0; j < atoms[l].size(); ++j, ++i) {
			atomNrOnInput[i] = atoms[l][j].atomNrOnInput;
			atomToWfnMap[i]  = atoms[l][j].atomToWfnMap;
			atomToTypeMap[i] = atoms[l][j].atomToTypeMap;

			atomicPositions[i]                  = atoms[l][j].atomicPositions;
			atomic_occupancy[i]                 = atoms[l][j].atomic_occupancy;
			atomic_multiplicity_factor[i]       = atoms[l][j].atomic_multiplicity_factor;
			atomicDisplacementParametersSize[i] = atoms[l][j].atomicDisplacementParametersSize;

			if ( atoms[l][j].atomicDisplacementParametersSize == 1)
				atomicDisplacementParameters[i] = atoms[l][j].atomicDisplacementParameters;
			else {
				atomicDisplacementParametersVec1[i] = atoms[l][j].atomicDisplacementParametersVec1;
				atomicDisplacementParametersVec2[i] = atoms[l][j].atomicDisplacementParametersVec2;
			}

			localCoordsRow1[i] = atoms[l][j].localCoordsRow1;
			localCoordsRow2[i] = atoms[l][j].localCoordsRow2;
			localCoordsRow3[i] = atoms[l][j].localCoordsRow3;
		}
	}
	nAtoms0 = atoms[0].size();
	nAtoms2 = atoms[3].size();
	nAtoms4 = atoms[5].size();

	// skondensowany memcpy kilkunastu tablic
	n_bytes = countPointerDiff((void *)atomicDisplacementParameters, (void *)(localCoordsRow3+nAtoms));
	cudaSafeCall(cudaMemcpy(devAtomicDisplacementParameters, atomicDisplacementParameters, n_bytes, cudaMemcpyHostToDevice));

	in>> auxString >> auxString; // "h vectors"
	for(int i=0;i<nH_vec;i++)
		in>>hVectors[i].x >>hVectors[i].y >>hVectors[i].z;

	cudaSafeCall(cudaMemcpy(devHVectors, hVectors, nH_vec*sizeof(hVectors[0]), cudaMemcpyHostToDevice));
    //dTargetDf;
    //devDTargetDf;
    //dTarget_dF_real
    //dTarget_dF_imag
    for(int i=0;i<nH_vec;i++)
    {
        dTargetDf[i].real = dTarget_dF_real[i];
        dTargetDf[i].imag = dTarget_dF_imag[i];
    }
    cudaSafeCall(cudaMemcpy(devDTargetDf, dTargetDf, nH_vec*sizeof(dTargetDf[0]), cudaMemcpyHostToDevice));
}

//END-READ INPUT ALLOCATE MEMORY--------------------------------------------}}}


void klasa::passInputDataAllocateMemoryOpt(
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
    const std::vector<std::vector<double> > &atomic_positions,
    const std::vector<std::vector<double> > &adps,
    const std::vector<double> &occupancy,
    const std::vector<double> &multiplicity_weigths,
    const std::vector<std::vector<std::vector<double> > > &lcs,
    const std::vector<vector<double> > &hkl,
    const std::vector<double> &dTarget_dF_real,
    const std::vector<double> &dTarget_dF_imag)
{

    nSymmOp = symmOps.size();
	nWfnParams = core_coeff.size();
	nTypes = p_lm.size();
	nAtoms = atomic_positions.size();
	nH_vec = hkl.size();

	if (nTypes > MAXATOMTYPES) {
		printf("too many atom types for GPU calulations\n");
		exit(0);
	}

	if (nWfnParams > MAXWFNPARAMS ) {
		printf("to many atomic wave-function types for GPU calulations\n");
		exit(0);
	}

	if (nSymmOp > MAXSYMMETRYOPERATIONS) {
		printf("too many symmetry operations for GPU calulations\n");
		exit(0);
	}

	

	nThreads= THREADSMAINKERNEL;
	nBlocks = (nH_vec + nThreads - 1)/nThreads;

	size_t memorySizeHost = 0;
	size_t memorySizeGPU  = 0;

	// liczymy, ile musimy zaalokowac
	memorySizeHost += countCommonMemory();
	memorySizeGPU  += memorySizeHost;

	memorySizeHost += countHostonlyMemory();
	memorySizeGPU  += countGPUonlyMemory();

	// malloc duzego bloku pamieci
	cudaSafeCall(cudaMallocHost((void**)&unitedHostMemory, memorySizeHost));
	allocatedHostPointers.push_back(unitedHostMemory);
	cudaSafeCall(cudaMalloc((void**)&unitedGPUMemory, memorySizeGPU));
	allocatedGpuPointers.push_back(unitedGPUMemory);

	// i podpiecie odpowiednich wskaznikow
	pinHostPointersToUnitedMemory();
	pinGPUPointersToUnitedMemory();

	// zerowanie sporej ilosci tablic wynikowych
	size_t n_bytes = countPointerDiff(unitedGPUMemory, devRotationsRow1);
	cudaSafeCall(cudaMemset(unitedGPUMemory, 0, n_bytes));


    
	for(int i=0;i<nSymmOp;i++) {
        rotationsRow1[i].x = symmOps[i][0]; 
        rotationsRow1[i].y = symmOps[i][1]; 
        rotationsRow1[i].z = symmOps[i][2];
        rotationsRow2[i].x = symmOps[i][3];
        rotationsRow2[i].y = symmOps[i][4]; 
        rotationsRow2[i].z = symmOps[i][5];
		rotationsRow3[i].x = symmOps[i][6];
        rotationsRow3[i].y = symmOps[i][7];
        rotationsRow3[i].z = symmOps[i][8];
		translations[i].x = symmOps[i][9];
        translations[i].y = symmOps[i][10];
        translations[i].z = symmOps[i][11];
	}


    
	// memcpy tablic rotationsRow1,2,3 i translations za pomoca 1 operacji
	n_bytes = countPointerDiff((void *)rotationsRow1, (void *)(translations+nSymmOp));
	cudaSafeCall(cudaMemcpy(devRotationsRow1 , rotationsRow1 , n_bytes , cudaMemcpyHostToDevice));

	//###### wfn type parameters

	for(int i=0;i<nWfnParams;i++) {
		wfnCore[i].n = core_coeff[i].size();
		ALLOCATEHOSTPOINTER(wfnCore[i].coeff , wfnCore[i].n);
		ALLOCATEHOSTPOINTER(wfnCore[i].exp   , wfnCore[i].n);
		ALLOCATEHOSTPOINTER(wfnCore[i].pow   , wfnCore[i].n);

		for(int j=0;j<wfnCore[i].n;j++) {
			wfnCore[i].coeff[j] = core_coeff[i][j];
            wfnCore[i].pow[j] = core_pow[i][j];
            wfnCore[i].exp[j] = core_exp[i][j];
			wfnCore[i].pow[j] += 2;
		}



		// spherical valence
		wfnValence[i].n = sph_val_coeff.size();
		ALLOCATEHOSTPOINTER(wfnValence[i].coeff , wfnValence[i].n);
		ALLOCATEHOSTPOINTER(wfnValence[i].exp   , wfnValence[i].n);
		ALLOCATEHOSTPOINTER(wfnValence[i].pow   , wfnValence[i].n);

		for(int j=0;j<wfnValence[i].n;j++) {
			wfnValence[i].coeff[j] = sph_val_coeff[i][j];
            wfnValence[i].pow[j] = sph_val_pow[i][j];
            wfnValence[i].exp[j] = sph_val_exp[i][j];
			wfnValence[i].pow[j] += 2;
		}
        
		wfnParams[i].defValenceExp = def_val_exponent[i];
		wfnParams[i].n = def_val_n_l.size();

		hostConstDefValenceExp[i] = wfnParams[i].defValenceExp;

		ALLOCATEHOSTPOINTER(wfnParams[i].defValencePow, wfnParams[i].n);

		//w kodzie zawsze brane jest +2
		for(int j=0;j<wfnParams[i].n;j++) 
			wfnParams[i].defValencePow[j] = def_val_n_l[i][j] + 2;
		
	}
	cudaSafeCall(cudaMemcpyToSymbol(constDefValenceExp, hostConstDefValenceExp, nWfnParams*sizeof(hostConstDefValenceExp[0])));

    

	for(int i=0;i<nWfnParams;i++) {
		for(int j=0;j<wfnParams[i].n;j++) {
			hostConstDefValencePow[i*5+j] = wfnParams[i].defValencePow[j];
		}
	}
	cudaSafeCall(cudaMemcpyToSymbol(constDefValencePow, hostConstDefValencePow, nWfnParams*5*sizeof(int)));


	// memcpy tablic wfnCore/Valence/Params za pomoca 1 operacji
	n_bytes = countPointerDiff((void *)wfnCore, (void *)(wfnParams+nWfnParams));
	cudaSafeCall(cudaMemcpy(devWfnCore, wfnCore, n_bytes, cudaMemcpyHostToDevice));

    //###### atom type parameters



	for(int i=0;i<nTypes;i++) {

		real kappaDefValence = kappa_def_valence[i];
		real kappaSpherical = kappa_spherical[i];

        typeParams[i].pVal = p_val[i];
		typeParams[i].maxL = p_lm[i].size()-1;

		typeParams[i].invKappaDefValence = 1.0/kappaDefValence;
		typeParams[i].invKappaSpherical = 1.0/kappaSpherical;

		ALLOCATEHOSTPOINTER(typeParams[i].pLM, (typeParams[i].maxL+1)*(2*(typeParams[i].maxL+1)+1));

		//printf("MAX L %d\n", typeParams[i].maxL);
		int counter = 0 ;
		for(int j=0;j<=typeParams[i].maxL;j++) {

			for(int k=0;k<2*j+1;k++) {
				typeParams[i].pLM[j*(2*typeParams[i].maxL+1)+k] = p_lm[i][j][k];
				hostConstPLM[i*SPHHARMCOEF+counter] = typeParams[i].pLM[j*(2*typeParams[i].maxL+1)+k]; //TODO: ZONK dla maxL roznego od 4 (20.07.14 by szmigacz)
				counter++;
			}
		}

	}

	cudaSafeCall(cudaMemcpyToSymbol(constPLM, hostConstPLM, nTypes*SPHHARMCOEF*sizeof(hostConstPLM[0])));

	cudaSafeCall(cudaMemcpy(devTypeParams, typeParams, nTypes*sizeof(typeParams[0]), cudaMemcpyHostToDevice));

    //###### atomic parameters

	std::vector<Atom> atoms[MAX_L]; // koszyki do sortowania po liczbie harmonijek
	for(int i=0;i<nAtoms;i++) {
		Atom a;
		a.atomNrOnInput = i;
		a.atomToWfnMap = atomToWfnTypeMap[i];
        a.atomToTypeMap = atomToAtomTypeMap[i];

		a.atomicPositions.x = atomic_positions[i][0];
        a.atomicPositions.y = atomic_positions[i][1]; 
        a.atomicPositions.z = atomic_positions[i][2];
		a.atomic_occupancy = occupancy[i];
        a.atomic_multiplicity_factor = multiplicity_weigths[i];

		a.atomicDisplacementParametersSize = adps[i].size();

		if ( a.atomicDisplacementParametersSize == 1)
			a.atomicDisplacementParameters = adps[i][0];
		else {
			a.atomicDisplacementParametersVec1.x = adps[i][0];
            a.atomicDisplacementParametersVec1.y = adps[i][1];
            a.atomicDisplacementParametersVec1.z = adps[i][2];
			a.atomicDisplacementParametersVec2.x = adps[i][3];
            a.atomicDisplacementParametersVec2.y = adps[i][4];
            a.atomicDisplacementParametersVec2.z = adps[i][5];

			a.atomicDisplacementParametersVec2.x *= 2;
			a.atomicDisplacementParametersVec2.y *= 2;
			a.atomicDisplacementParametersVec2.z *= 2;
		}

		a.localCoordsRow1.x = lcs[i][0][0];
        a.localCoordsRow1.y = lcs[i][0][1];
        a.localCoordsRow1.z = lcs[i][0][2];
		a.localCoordsRow2.x = lcs[i][1][0];
        a.localCoordsRow2.y = lcs[i][1][1];
        a.localCoordsRow2.z = lcs[i][1][2];
		a.localCoordsRow3.x = lcs[i][2][0];
        a.localCoordsRow3.y = lcs[i][2][1];
        a.localCoordsRow3.z = lcs[i][2][2];

		atoms[wfnParams[a.atomToWfnMap].n].push_back(a); // dodanie do odpowiedniego koszyka
	}
	// kopiowanie posortowanych danych o atomach

	for (int l = 0, i = 0; l < MAX_L; ++l) {
		for (unsigned int j = 0; j < atoms[l].size(); ++j, ++i) {
			atomNrOnInput[i] = atoms[l][j].atomNrOnInput;
			atomToWfnMap[i]  = atoms[l][j].atomToWfnMap;
			atomToTypeMap[i] = atoms[l][j].atomToTypeMap;

			atomicPositions[i]                  = atoms[l][j].atomicPositions;
			atomic_occupancy[i]                 = atoms[l][j].atomic_occupancy;
			atomic_multiplicity_factor[i]       = atoms[l][j].atomic_multiplicity_factor;
			atomicDisplacementParametersSize[i] = atoms[l][j].atomicDisplacementParametersSize;

			if ( atoms[l][j].atomicDisplacementParametersSize == 1)
				atomicDisplacementParameters[i] = atoms[l][j].atomicDisplacementParameters;
			else {
				atomicDisplacementParametersVec1[i] = atoms[l][j].atomicDisplacementParametersVec1;
				atomicDisplacementParametersVec2[i] = atoms[l][j].atomicDisplacementParametersVec2;
			}

			localCoordsRow1[i] = atoms[l][j].localCoordsRow1;
			localCoordsRow2[i] = atoms[l][j].localCoordsRow2;
			localCoordsRow3[i] = atoms[l][j].localCoordsRow3;
		}
	}
	nAtoms0 = atoms[0].size();
	nAtoms2 = atoms[3].size();
	nAtoms4 = atoms[5].size();

	// skondensowany memcpy kilkunastu tablic

	n_bytes = countPointerDiff((void *)atomicDisplacementParameters, (void *)(localCoordsRow3+nAtoms));
	cudaSafeCall(cudaMemcpy(devAtomicDisplacementParameters, atomicDisplacementParameters, n_bytes, cudaMemcpyHostToDevice));
    
    //##### HKL i dTarget_dF
	
	for(int i=0;i<nH_vec;i++)
    {
        hVectors[i].x = hkl[i][0];
        hVectors[i].y = hkl[i][1]; 
        hVectors[i].z = hkl[i][2];
    }
	cudaSafeCall(cudaMemcpy(devHVectors, hVectors, nH_vec*sizeof(hVectors[0]), cudaMemcpyHostToDevice));

    for(int i=0;i<nH_vec;i++)
    {
        dTargetDf[i].real = dTarget_dF_real[i];
        dTargetDf[i].imag = dTarget_dF_imag[i];
    }
    cudaSafeCall(cudaMemcpy(devDTargetDf, dTargetDf, nH_vec*sizeof(dTargetDf[0]), cudaMemcpyHostToDevice));

}