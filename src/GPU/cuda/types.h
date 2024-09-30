// vim: noet: sw=3: ts=3
#ifndef TYPES_H_D7XT0Z5L
#define TYPES_H_D7XT0Z5L
//LICENSE-------------------------------------------------------------------{{{
/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
//END-LICENSE---------------------------------------------------------------}}}

#include <vector>
//#include <sstream>
#include "globals.h"

//typedef cudaComplexFloat cudaComplexType;
//typedef float real;
//typedef float3 real3;

typedef cudaComplexDouble cudaComplexType;
typedef double real;
typedef double3 real3;

struct atomParameters {
	int * atomToWfnMap;
	int * atomToTypeMap;
	real * atomic_occupancy;
	real * atomic_multiplicity_factor;

};

struct wfnCoreParams{
	real *coeff;
	real *exp;
	int *pow;
	int n;
};

struct wfnValenceParams{
	real *coeff;
	real *exp;
	int *pow;
	int n;
};


struct typeParamsValues{
	real invKappaSpherical;
	real invKappaDefValence;
	real pVal;
	real * pLM;
	int maxL;
};

struct wfnParameters{
	cudaComplexType anomalousScattering;
	real defValenceExp;
	int * defValencePow;
	int n;
};

class klasa{
	public:
		char * inputFileName;
		char * outputFileName;

		int nBlocks;
		int nThreads;

		int nSymmOp;
		int nWfnParams;
		int nTypes;
		int nAtoms;
		int nAtoms0, nAtoms2, nAtoms4;
		int nH_vec;

		real * hostConstPLM;
		int  * hostConstDefValencePow;
		real * hostConstDefValenceExp;

		void * unitedHostMemory;
		void * unitedGPUMemory;

		//to jest bardzo brzydkie TODO: poprawic
		cudaComplexType * occupancyDerivatives;
		cudaComplexType * devOccupancyDerivatives;
		cudaComplexType * devReduceOccupancyDerivatives;

		cudaComplexType * positionDerivativesX;
		cudaComplexType * devPositionDerivativesX;
		cudaComplexType * devReducePositionDerivativesX;

		cudaComplexType * positionDerivativesY;
		cudaComplexType * devPositionDerivativesY;
		cudaComplexType * devReducePositionDerivativesY;

		cudaComplexType * positionDerivativesZ;
		cudaComplexType * devPositionDerivativesZ;
		cudaComplexType * devReducePositionDerivativesZ;

		cudaComplexType * adpDerivatives0;
		cudaComplexType * devAdpDerivatives0;
		cudaComplexType * devReduceAdpDerivatives0;

		cudaComplexType * adpDerivatives1;
		cudaComplexType * devAdpDerivatives1;
		cudaComplexType * devReduceAdpDerivatives1;

		cudaComplexType * adpDerivatives2;
		cudaComplexType * devAdpDerivatives2;
		cudaComplexType * devReduceAdpDerivatives2;

		cudaComplexType * adpDerivatives3;
		cudaComplexType * devAdpDerivatives3;
		cudaComplexType * devReduceAdpDerivatives3;

		cudaComplexType * adpDerivatives4;
		cudaComplexType * devAdpDerivatives4;
		cudaComplexType * devReduceAdpDerivatives4;

		cudaComplexType * adpDerivatives5;
		cudaComplexType * devAdpDerivatives5;
		cudaComplexType * devReduceAdpDerivatives5;

		cudaComplexType * dTargetDf;
		cudaComplexType * devDTargetDf;


		wfnCoreParams    * wfnCore;
		wfnCoreParams    * devWfnCore;
		wfnCoreParams    * helperWfnCore;

		wfnValenceParams * wfnValence;
		wfnValenceParams * devWfnValence;
		wfnValenceParams * helperWfnValence;

		wfnParameters    * wfnParams;
		wfnParameters    * devWfnParams;
		wfnParameters    * helperWfnParams;

		typeParamsValues * typeParams;
		typeParamsValues * devTypeParams;
		typeParamsValues * helperTypeParams;


		real3 * hVectors;
		real3 * devHVectors;

		int * atomNrOnInput;
		int * atomToWfnMap;
		int   * devAtomToWfnMap;

		int * atomToTypeMap;
		int   * devAtomToTypeMap;

		real * atomic_occupancy;
		real  * devAtomic_occupancy;

		real * atomic_multiplicity_factor;
		real  * devAtomic_multiplicity_factor;

		real3 * atomicPositions;
		real3 * devAtomicPositions;

		int * atomicDisplacementParametersSize;
		int * devAtomicDisplacementParametersSize;


		real3 * atomicDisplacementParametersVec1;
		real3 * atomicDisplacementParametersVec2;
		real  * atomicDisplacementParameters;
		real3 * devAtomicDisplacementParametersVec1;
		real3 * devAtomicDisplacementParametersVec2;
		real *  devAtomicDisplacementParameters;

		real3 * devLocalCoordsRow1;
		real3 * devLocalCoordsRow2;
		real3 * devLocalCoordsRow3;

		real3 * localCoordsRow1;
		real3 * localCoordsRow2;
		real3 * localCoordsRow3;

		real3 * rotationsRow1;
		real3 * rotationsRow2;
		real3 * rotationsRow3;
		real3 * translations;

		real3 * devRotationsRow1;
		real3 * devRotationsRow2;
		real3 * devRotationsRow3;
		real3 * devTranslations;

		cudaComplexType * arrayF;
		cudaComplexType * devArrayF;

		std::vector<void*> allocatedHostPointers;
		std::vector<void*> allocatedGpuPointers;

		//klasa(char * _inputFileName, char * _outputFileName): inputFileName(_inputFileName), outputFileName(_outputFileName) {}
        klasa(){}
		void readInputAllocateMemory();
		void readInputAllocateMemoryOpt(
                 std::stringstream &in,
                 const std::vector<double> &dTarget_dF_real,
                 const std::vector<double> &dTarget_dF_imag);
        
        void passInputDataAllocateMemoryOpt(
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
            const std::vector<std::vector<double> > &hkl,
            const std::vector<double> &dTarget_dF_real,
            const std::vector<double> &dTarget_dF_imag);

		size_t countCommonMemory();
		size_t countHostonlyMemory();
		size_t countGPUonlyMemory();
		void pinHostPointersToUnitedMemory();
		void pinGPUPointersToUnitedMemory();

		void initializeConst();

		void freeHostMemory();
		void freeGpuMemory();
		void writeOutput(std::vector<double> &f_real,
                         std::vector<double> &f_imag,
                         std::vector<double> &occupancyDerivatives,
                         std::vector<double> &xyzDerivatives,
                         std::vector<double> &adpDerivatives);
		void runGpu();
		void copyResultsToCpu();
};

#endif /* end of include guard: TYPES_H_D7XT0Z5L */
