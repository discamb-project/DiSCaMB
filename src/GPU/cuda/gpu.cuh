// vim: noet: sw=3: ts=3: foldmethod=marker: foldlevel=0
//INCLUDES------------------------------------------------------------------{{{
#include"realSphericalHarmonics.cuh"
#include"fakeComplex.cuh"
#include"gFunctions.cuh"
#include<cmath>
#include<cstdio>
#include<vector_types.h>
#include"types.h"
#include"cudaReductions.cuh"
#include"cudaConstant.cuh"
#include"customLDG.cuh"
//END-INCLUDES--------------------------------------------------------------}}}
//PRINT DEBUGS--------------------------------------------------------------{{{
__device__ __forceinline__ void
debugReal(real liczba) {
	if (threadIdx.x==0)
		printf("%.10e\n", liczba);
}
__device__ __forceinline__ void
debugComplex(cudaComplexType liczba) {
	if (threadIdx.x==0)
		printf("%.10e %.10e\n", liczba.real, liczba.imag);
}
//END-PRINT DEBUGS----------------------------------------------------------}}}
//SHFL XOR CUDACOMPLEX------------------------------------------------------{{{
__forceinline__ __device__ cudaComplexType __shfl_xor(cudaComplexType value, int laneMask, int width=WARPSIZE) {
	cudaComplexType result;
	result.real = __hiloint2double(
			__shfl_xor(__double2hiint(value.real), laneMask, width),
			__shfl_xor(__double2loint(value.real), laneMask, width)
			);
	result.imag = __hiloint2double(
			__shfl_xor(__double2hiint(value.imag), laneMask, width),
			__shfl_xor(__double2loint(value.imag), laneMask, width)
			);
	return result;

}
//END-SHFL XOR CUDACOMPLEX--------------------------------------------------}}}
//UNROLL SPHERICAL HARMONICS------------------------------------------------{{{
//template zeby ladnie rozwijac kod na harmoniki, dla kazdego L mamy pelen
//zakres m od -L do L, bez template trzeba wklejac kod recznie
//TODO: zostawiam stary kod na LDG zeby zrobic pozniej wersje jak brakuje const
//memory
template<int l, int maxM, typename T>
__device__ __forceinline__ T unrollSphericalHarmonics(
		//T const * const pm,
		int const atomTypeIndex,
		T const x,
		T const y,
		T const z)
{
	int offset; //offsety do tablicy ze wspolczynnikami do harmonik sferycznych
	if (l==0) offset =0;
	if (l==1) offset =1;
	if (l==2) offset =4;
	if (l==3) offset =9;
	if (l==4) offset =16;

	if (maxM == 0) {
		//T temp3 = LDG(pm[l]) * realSphericalHarmonicsOpt<l,0,T>(x,y,z);
		T temp3 = constPLM[atomTypeIndex * (SPHHARMCOEF) + offset +l] * realSphericalHarmonicsOpt<l,0,T>(x,y,z);
		return temp3;
	}
	else {
		//TODO: nie wiem jak to zopymalizowac zeby nie liczyc jak jest zero (20.07.14 by szmigacz)
		//T temp1 = LDG(pm[l-maxM]) * realSphericalHarmonicsOpt<l,-maxM, T>(x,y,z);
		//T temp2 = LDG(pm[l+maxM]) * realSphericalHarmonicsOpt<l, maxM, T>(x,y,z);
		T temp1 = constPLM[atomTypeIndex * (SPHHARMCOEF) + offset +l -maxM] * realSphericalHarmonicsOpt<l,-maxM, T>(x,y,z);
		T temp2 = constPLM[atomTypeIndex * (SPHHARMCOEF) + offset +l +maxM] * realSphericalHarmonicsOpt<l, maxM, T>(x,y,z);

		T temp3 = unrollSphericalHarmonics<l, (maxM>0?maxM-1:0) , T>(atomTypeIndex, x, y, z);
		return temp3 + temp1 + temp2;
	}
}
//END-UNROLL SPHERICAL HARMONICS--------------------------------------------}}}
//SCALAR PRODUCT------------------------------------------------------------{{{
__host__ __device__ __forceinline__
real scalarProduct(real3 a, real3 b) {
	real value = a.x*b.x + a.y*b.y + a.z*b.z;
	return value;
}
//END-SCALAR PRODUCT--------------------------------------------------------}}}
//WFN SPHERICAL SF----------------------------------------------------------{{{
__device__ __noinline__ //marker1
real wfnSphericalSF(
		real const * const coefficients,
		real const * const exponents,
		int const * const powers,
		const int nCoefficients,
		real const hLength
		) {
	real result = 0;

	//#pragma unroll 4 //tymczasowo tyle ustawiam, to pomaga ale nie chce tego tutaj miec wlaczonego
	for(int k=0;k<nCoefficients;k++) {
		int localPowers = LDG(powers[k]);
		real localExponents = LDG(exponents[k]);
		real localCoefficient = LDG(coefficients[k]);
		real tempVal = gFunctionOpt<0, real>(localPowers, hLength, localExponents);
		result += localCoefficient*tempVal;
	}
	return result;
}
//END-WFN SPHERICAL SF------------------------------------------------------}}}
//SLATER NORMALIZATION------------------------------------------------------{{{
__device__ __forceinline__ real slaterNormalization(real exponent, int powValence) {
	return pow(exponent, powValence +1) * constInverseFactorial[powValence];
}
//END-SLATER NORMALIZATION--------------------------------------------------}}}
//CALCULATE MULTIPOLAR SF---------------------------------------------------{{{
__device__ __noinline__ //marker3
cudaComplexType calculateMultipolarSF(
		real hVectorLength,
		real3 hVectorRotated,
		real const exponent,
		//int * powerR,
		int const atomWfnIndex,
		//real const * const pLM,
		int atomTypeIndex,
		int n
		) {
	cudaComplexType result;
	if (n == 0) {
		result = cudaComplexType(0,0);
	}
	if (n == 2 || n == 4) {
		real x = hVectorRotated.x;
		real y = hVectorRotated.y;
		real z = hVectorRotated.z;

		real rnorm = rsqrt(x*x + y*y + z*z);
		x *= rnorm;
		y *= rnorm;
		z *= rnorm;

		real radialF;
		int powValence;

		//TODO: zastanowic sie czy mozna pisac przez = czy trzeba przez += (29.07.14 by szmigacz)
		//for L = 0
		//int localPowerR = LDG(powerR[0]);
		powValence = constDefValencePow[atomWfnIndex*5];
		radialF =  gFunctionOpt<0,real>(powValence, hVectorLength, exponent);
		real temp = slaterNormalization(exponent, powValence);;


		radialF *= temp;
		result.real = radialF * unrollSphericalHarmonics<0,0,real>(atomTypeIndex,x,y,z);
		//for L = 1
		//int localpowValenceerR = LDG(powValenceerR[1]);
		powValence = constDefValencePow[atomWfnIndex*5+1];
		radialF =  gFunctionOpt<1,real>(powValence, hVectorLength, exponent);
		radialF *= slaterNormalization(exponent, powValence);
		result.imag = radialF * unrollSphericalHarmonics<1,1,real>(atomTypeIndex,x,y,z);
		//for L = 2
		//int localpowValenceerR = LDG(powValenceerR[2]);
		powValence = constDefValencePow[atomWfnIndex*5+2];
		radialF =  gFunctionOpt<2,real>(powValence, hVectorLength, exponent);
		radialF *= slaterNormalization(exponent, powValence);
		result.real -= radialF * unrollSphericalHarmonics<2,2,real>(atomTypeIndex,x,y,z);
		if (n == 4) {
			powValence = constDefValencePow[atomWfnIndex*5+3];
			radialF =  gFunctionOpt<3,real>(powValence, hVectorLength, exponent);
			radialF *= slaterNormalization(exponent, powValence);
			result.imag -= radialF * unrollSphericalHarmonics<3,3,real>(atomTypeIndex,x,y,z);
			powValence = constDefValencePow[atomWfnIndex*5+4];
			radialF =  gFunctionOpt<4,real>(powValence, hVectorLength, exponent);
			radialF *= slaterNormalization(exponent, powValence);
			result.real += radialF * unrollSphericalHarmonics<4,4,real>(atomTypeIndex,x,y,z);
		}

		result.real *= (4.0*M_PI);
		result.imag *= (4.0*M_PI);
	}
	return result;
}
//END-CALCULATE MULTIPOLAR SF-----------------------------------------------}}}
//MAIN LOOP-----------------------------------------------------------------{{{
__device__ __forceinline__ //marker2
void mainLoop( 
		const int nH_vec,
		real3 const * const hVectors,
		int const * const atomToWfnMap,
		int const * const atomToTypeMap,
		real const * const atomicOccupancy,
		real const * const atomicMultiplictyFactor,
		wfnCoreParams const * const wfnCore,
		wfnValenceParams const * const wfnValence,
		typeParamsValues const * const typeParams,
		real3 const * const localCoordsRow1,
		real3 const * const localCoordsRow2,
		real3 const * const localCoordsRow3,
		real3 const * const rotationsRow1,
		real3 const * const rotationsRow2,
		real3 const * const rotationsRow3,
		real3 const * const atomicPositions,
		real3 const * const translations,
		int const atomPositionsSize0,
		int const atomPositionsSize2,
		int const atomPositionsSize4,
		int const symmetryIndexSize,
		int const * const atomicDisplacementParametersSize,
		real const * const atomicDisplacementParameters,
		real3 const * const atomicDisplacementParametersVec1,
		real3 const * const atomicDisplacementParametersVec2,
		cudaComplexType * const arrayF,
		wfnParameters const * const wfnParams,
		cudaComplexType const * const dTarget_dF,
		cudaComplexType * const occupancyDerivatives,
		cudaComplexType * const positionDerivativesX,
		cudaComplexType * const positionDerivativesY,
		cudaComplexType * const positionDerivativesZ,
		cudaComplexType * const adpDerivatives0,
		cudaComplexType * const adpDerivatives1,
		cudaComplexType * const adpDerivatives2,
		cudaComplexType * const adpDerivatives3,
		cudaComplexType * const adpDerivatives4,
		cudaComplexType * const adpDerivatives5,
		int const hkl_index,
		real3 const hVector,
		real const hVectorLength,
		cudaComplexType & atomicFContribution,
		int const atomIndex,
		int const n_index)
{
	int atomWfnIndex =  LDG(atomToWfnMap[atomIndex]);
	int atomTypeIndex = LDG(atomToTypeMap[atomIndex]);
	cudaComplexType anomalousScattering = LDG(wfnParams[atomWfnIndex].anomalousScattering);

	real atomicMultiplicityRegister = LDG(atomicMultiplictyFactor[atomIndex]);
	real atomWeight = LDG(atomicOccupancy[atomIndex]) * atomicMultiplicityRegister;

	real atomFCore = wfnSphericalSF( 
			wfnCore[atomWfnIndex].coeff,
			wfnCore[atomWfnIndex].exp,
			wfnCore[atomWfnIndex].pow,
			wfnCore[atomWfnIndex].n,
			hVectorLength);

	real atomFSphVal = wfnSphericalSF(
			wfnValence[atomWfnIndex].coeff,
			wfnValence[atomWfnIndex].exp,
			wfnValence[atomWfnIndex].pow,
			wfnValence[atomWfnIndex].n,
			hVectorLength * typeParams[atomTypeIndex].invKappaSpherical);

	atomFSphVal *= typeParams[ atomTypeIndex ].pVal;
	real atomFSph = atomFCore + atomFSphVal;


	cudaComplexType accumulatorOccDer  = cudaComplexType(0,0);
	cudaComplexType accumulatorPosDerX = cudaComplexType(0,0);
	cudaComplexType accumulatorPosDerY = cudaComplexType(0,0);
	cudaComplexType accumulatorPosDerZ = cudaComplexType(0,0);

	cudaComplexType accumulatorAdpDer[N_ADP_DERIVATIVES];
	for (int i = 0; i < N_ADP_DERIVATIVES; ++i)
		accumulatorAdpDer[i] = cudaComplexType(0,0);

	if (hkl_index < nH_vec) {
		for(int symmetryIndex = 0; symmetryIndex < symmetryIndexSize; symmetryIndex++) {
			real invAtomKappaPrime = typeParams[ atomTypeIndex ].invKappaDefValence;
			real3 atomicPosition = atomicPositions[atomIndex]; /*TODO: zoptymalizowac (19.07.14 by szmigacz)*/ 

			real3 localCoordRow1 = localCoordsRow1[atomIndex];
			real3 localCoordRow2 = localCoordsRow2[atomIndex];
			real3 localCoordRow3 = localCoordsRow3[atomIndex];

			real3 rotationRow1 = LDG(rotationsRow1[symmetryIndex]);
			real3 rotationRow2 = LDG(rotationsRow2[symmetryIndex]);
			real3 rotationRow3 = LDG(rotationsRow3[symmetryIndex]);

			real3 translation  = LDG(translations[symmetryIndex]);

			real3 atomPosition;
			atomPosition.x = scalarProduct(rotationRow1, atomicPosition) + translation.x;
			atomPosition.y = scalarProduct(rotationRow2, atomicPosition) + translation.y;
			atomPosition.z = scalarProduct(rotationRow3, atomicPosition) + translation.z;

			cudaComplexType twoPiI = cudaComplexType(0, 2.0*M_PI);
			//cudaComplexType atomicPhaseFactor = exp( twoPiI * scalarProduct( hVector, atomPosition));
			real argument = 2 * scalarProduct(hVector, atomPosition);
			real cosArg;
			real sinArg;
			sincospi(argument, &sinArg, &cosArg);
			cudaComplexType atomicPhaseFactor = cudaComplexType(cosArg, sinArg);

			real3 hRotated;
			hRotated.x =
				hVector.x * rotationRow1.x +
				hVector.y * rotationRow2.x +
				hVector.z * rotationRow3.x;
			hRotated.y =
				hVector.x * rotationRow1.y +
				hVector.y * rotationRow2.y +
				hVector.z * rotationRow3.y;
			hRotated.z =
				hVector.x * rotationRow1.z +
				hVector.y * rotationRow2.z +
				hVector.z * rotationRow3.z;

			real3 hScaledRotated;
			hScaledRotated.x = hRotated.x * invAtomKappaPrime;
			hScaledRotated.y = hRotated.y * invAtomKappaPrime;
			hScaledRotated.z = hRotated.z * invAtomKappaPrime;

			real3 hVectorRotated;
			hVectorRotated.x = 
				hScaledRotated.x * localCoordRow1.x +
				hScaledRotated.y * localCoordRow2.x +
				hScaledRotated.z * localCoordRow3.x;
			hVectorRotated.y = 
				hScaledRotated.x * localCoordRow1.y +
				hScaledRotated.y * localCoordRow2.y +
				hScaledRotated.z * localCoordRow3.y;
			hVectorRotated.z = 
				hScaledRotated.x * localCoordRow1.z +
				hScaledRotated.y * localCoordRow2.z +
				hScaledRotated.z * localCoordRow3.z;

			real hVectorRotatedLength = sqrt(scalarProduct(hScaledRotated,hScaledRotated));

			cudaComplexType atomFDefVal;

			atomFDefVal = calculateMultipolarSF(
					hVectorRotatedLength,
					hVectorRotated,
					constDefValenceExp[atomWfnIndex],
					//wfnParams[atomWfnIndex].defValenceExp,
					//wfnParams[atomWfnIndex].defValencePow,
					atomWfnIndex,
					//typeParams[atomTypeIndex].pLM,
					atomTypeIndex,
					n_index
					);

			int atomicDispParamsSizeRegister = atomicDisplacementParametersSize[atomIndex];
			real atomicDispParamsRegister    = atomicDisplacementParameters[atomIndex];
			real3 adpVec1 = atomicDisplacementParametersVec1[atomIndex];
			real3 adpVec2 = atomicDisplacementParametersVec2[atomIndex];
			real temperatureFactor;
			if (atomicDispParamsSizeRegister == 1) { /* tego if-else'a mozna zoptymalizowac, wybierajac zawsze dobra galaz, ale nie robic tego, bo to spowalnia program ;) (5.08.14 by kkiewicz) */ 
				temperatureFactor = exp( -hVectorLength*hVectorLength * atomicDispParamsRegister);
			}
			else {
				real product = 
					hRotated.x*(hRotated.x*adpVec1.x +
							hRotated.y*adpVec2.x+
							hRotated.z*adpVec2.y) +
					hRotated.y*(hRotated.y*adpVec1.y +
							hRotated.z*adpVec2.z) +
					hRotated.z*hRotated.z*adpVec1.z;
				temperatureFactor = exp(-product);
			}

			cudaComplexType atomFSphComplex = cudaComplexType(atomFSph,0);
			cudaComplexType transformedAtomF = (atomFSphComplex + atomFDefVal + anomalousScattering) *atomicPhaseFactor*temperatureFactor;

			cudaComplexType temporary1 = transformedAtomF * atomWeight;
			atomicFContribution += temporary1;

			//cudaComplexType one = cudaComplexType(1,0);
            cudaComplexType one = dTarget_dF[hkl_index];
			accumulatorOccDer += one * transformedAtomF * atomicMultiplicityRegister;
			//dTarget_dF[hkl_index] * transformedAtomF * atomicMultiplictyFactor[atomIndex];

			accumulatorPosDerX += one * twoPiI * hRotated.x * temporary1; 
			accumulatorPosDerY += one * twoPiI * hRotated.y * temporary1;
			accumulatorPosDerZ += one * twoPiI * hRotated.z * temporary1;

			accumulatorAdpDer[0] -= one * hRotated.x * hRotated.x * temporary1;
			accumulatorAdpDer[1] -= one * hRotated.y * hRotated.y * temporary1;
			accumulatorAdpDer[2] -= one * hRotated.z * hRotated.z * temporary1;
			accumulatorAdpDer[3] -= one * hRotated.x * hRotated.y * temporary1;
			accumulatorAdpDer[4] -= one * hRotated.x * hRotated.z * temporary1;
			accumulatorAdpDer[5] -= one * hRotated.y * hRotated.z * temporary1;

		}
	}
	int warpId = threadIdx.x/WARPSIZE;
	int laneId;
	asm("mov.u32 %0, %laneid;" : "=r"(laneId));

	//performs a reduction across the warp
	reduceInWarp(accumulatorOccDer );
	reduceInWarp(accumulatorPosDerX);
	reduceInWarp(accumulatorPosDerY);
	reduceInWarp(accumulatorPosDerZ);
	for (int i = 0; i < N_ADP_DERIVATIVES; ++i)
		reduceInWarp(accumulatorAdpDer[i]);

	//one thread per warp saves the results to global
	if (laneId == 0 ) {
		int idx = (THREADSMAINKERNEL/WARPSIZE)*(gridDim.x * atomIndex + blockIdx.x)+warpId;
		occupancyDerivatives[idx] = accumulatorOccDer ;
		positionDerivativesX[idx] = accumulatorPosDerX;
		positionDerivativesY[idx] = accumulatorPosDerY;
		positionDerivativesZ[idx] = accumulatorPosDerZ;
		adpDerivatives0[idx] = accumulatorAdpDer[0];
		adpDerivatives1[idx] = accumulatorAdpDer[1];
		adpDerivatives2[idx] = accumulatorAdpDer[2];
		adpDerivatives3[idx] = accumulatorAdpDer[3];
		adpDerivatives4[idx] = accumulatorAdpDer[4];
		adpDerivatives5[idx] = accumulatorAdpDer[5];
	}
}
//END-MAIN LOOP-------------------------------------------------------------}}}
//CALCULATE SF KERNEL-------------------------------------------------------{{{
	__global__
//__launch_bounds__(THREADSMAINKERNEL)
__launch_bounds__(THREADSMAINKERNEL,8)
void calculateSFKernel (
		const int nH_vec,
		real3 const * const hVectors,
		int const * const atomToWfnMap,
		int const * const atomToTypeMap,
		real const * const atomicOccupancy,
		real const * const atomicMultiplictyFactor,
		wfnCoreParams const * const wfnCore,
		wfnValenceParams const * const wfnValence,
		typeParamsValues const * const typeParams,
		real3 const * const localCoordsRow1,
		real3 const * const localCoordsRow2,
		real3 const * const localCoordsRow3,
		real3 const * const rotationsRow1,
		real3 const * const rotationsRow2,
		real3 const * const rotationsRow3,
		real3 const * const atomicPositions,
		real3 const * const translations,
		int const atomPositionsSize0,
		int const atomPositionsSize2,
		int const atomPositionsSize4,
		int const symmetryIndexSize,
		int const * const atomicDisplacementParametersSize,
		real const * const atomicDisplacementParameters,
		real3 const * const atomicDisplacementParametersVec1,
		real3 const * const atomicDisplacementParametersVec2,
		cudaComplexType * const arrayF,
		wfnParameters const * const wfnParams,
		cudaComplexType const * const dTarget_dF,
		cudaComplexType * const occupancyDerivatives,
		cudaComplexType * const positionDerivativesX,
		cudaComplexType * const positionDerivativesY,
		cudaComplexType * const positionDerivativesZ,
		cudaComplexType * const adpDerivatives0,
		cudaComplexType * const adpDerivatives1,
		cudaComplexType * const adpDerivatives2,
		cudaComplexType * const adpDerivatives3,
		cudaComplexType * const adpDerivatives4,
		cudaComplexType * const adpDerivatives5
		)
{
	int hkl_index = threadIdx.x + blockDim.x * blockIdx.x;

	real3 hVector = hVectors[hkl_index];
	real hVectorLength = sqrt( hVector.x*hVector.x + hVector.y*hVector.y + hVector.z*hVector.z);
	cudaComplexType atomicFContribution = cudaComplexType(0,0);

	int atomIndex = 0;
	for(; atomIndex < atomPositionsSize0; atomIndex++) {
		mainLoop(
				nH_vec,
				hVectors,
				atomToWfnMap,
				atomToTypeMap,
				atomicOccupancy,
				atomicMultiplictyFactor,
				wfnCore,
				wfnValence,
				typeParams,
				localCoordsRow1,
				localCoordsRow2,
				localCoordsRow3,
				rotationsRow1,
				rotationsRow2,
				rotationsRow3,
				atomicPositions,
				translations,
				atomPositionsSize0,
				atomPositionsSize2,
				atomPositionsSize4,
				symmetryIndexSize,
				atomicDisplacementParametersSize,
				atomicDisplacementParameters,
				atomicDisplacementParametersVec1,
				atomicDisplacementParametersVec2,
				arrayF,
				wfnParams,
				dTarget_dF,
				occupancyDerivatives,
				positionDerivativesX,
				positionDerivativesY,
				positionDerivativesZ,
				adpDerivatives0,
				adpDerivatives1,
				adpDerivatives2,
				adpDerivatives3,
				adpDerivatives4,
				adpDerivatives5,
				hkl_index,
				hVector,
				hVectorLength,
				atomicFContribution,
				atomIndex,
				0);
	}
	for(; atomIndex < atomPositionsSize0+atomPositionsSize2; atomIndex++) {
		mainLoop(
				nH_vec,
				hVectors,
				atomToWfnMap,
				atomToTypeMap,
				atomicOccupancy,
				atomicMultiplictyFactor,
				wfnCore,
				wfnValence,
				typeParams,
				localCoordsRow1,
				localCoordsRow2,
				localCoordsRow3,
				rotationsRow1,
				rotationsRow2,
				rotationsRow3,
				atomicPositions,
				translations,
				atomPositionsSize0,
				atomPositionsSize2,
				atomPositionsSize4,
				symmetryIndexSize,
				atomicDisplacementParametersSize,
				atomicDisplacementParameters,
				atomicDisplacementParametersVec1,
				atomicDisplacementParametersVec2,
				arrayF,
				wfnParams,
				dTarget_dF,
				occupancyDerivatives,
				positionDerivativesX,
				positionDerivativesY,
				positionDerivativesZ,
				adpDerivatives0,
				adpDerivatives1,
				adpDerivatives2,
				adpDerivatives3,
				adpDerivatives4,
				adpDerivatives5,
				hkl_index,
				hVector,
				hVectorLength,
				atomicFContribution,
				atomIndex,
				2);
	}
	for(; atomIndex < atomPositionsSize0+atomPositionsSize2+atomPositionsSize4; atomIndex++) {
		mainLoop(
				nH_vec,
				hVectors,
				atomToWfnMap,
				atomToTypeMap,
				atomicOccupancy,
				atomicMultiplictyFactor,
				wfnCore,
				wfnValence,
				typeParams,
				localCoordsRow1,
				localCoordsRow2,
				localCoordsRow3,
				rotationsRow1,
				rotationsRow2,
				rotationsRow3,
				atomicPositions,
				translations,
				atomPositionsSize0,
				atomPositionsSize2,
				atomPositionsSize4,
				symmetryIndexSize,
				atomicDisplacementParametersSize,
				atomicDisplacementParameters,
				atomicDisplacementParametersVec1,
				atomicDisplacementParametersVec2,
				arrayF,
				wfnParams,
				dTarget_dF,
				occupancyDerivatives,
				positionDerivativesX,
				positionDerivativesY,
				positionDerivativesZ,
				adpDerivatives0,
				adpDerivatives1,
				adpDerivatives2,
				adpDerivatives3,
				adpDerivatives4,
				adpDerivatives5,
				hkl_index,
				hVector,
				hVectorLength,
				atomicFContribution,
				atomIndex,
				4);
	}
	arrayF [ hkl_index ] = atomicFContribution;
}
//END-CALCULATE SF KERNEL---------------------------------------------------}}}
//REDUCE DERIVATIVES KERNEL-------------------------------------------------{{{
//tego kernela nie ma sensu optymalizowac, dziala i tyle
__global__ void reduceDerivatives (
		cudaComplexType * const reduceOccupancyDerivatives,
		cudaComplexType * const occupancyDerivatives,
		cudaComplexType * const reducePositionDerivativesX,
		cudaComplexType * const positionDerivativesX,
		cudaComplexType * const reducePositionDerivativesY,
		cudaComplexType * const positionDerivativesY,
		cudaComplexType * const reducePositionDerivativesZ,
		cudaComplexType * const positionDerivativesZ,
		cudaComplexType * const reduceAdpDerivative0,
		cudaComplexType * const adpDerivatives0,
		cudaComplexType * const reduceAdpDerivative1,
		cudaComplexType * const adpDerivatives1,
		cudaComplexType * const reduceAdpDerivative2,
		cudaComplexType * const adpDerivatives2,
		cudaComplexType * const reduceAdpDerivative3,
		cudaComplexType * const adpDerivatives3,
		cudaComplexType * const reduceAdpDerivative4,
		cudaComplexType * const adpDerivatives4,
		cudaComplexType * const reduceAdpDerivative5,
		cudaComplexType * const adpDerivatives5,
		int N
		)
{
	__shared__ volatile cudaComplexType shOccDer[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shPosDerX[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shPosDerY[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shPosDerZ[THREADSREDUCEKERNEL];

	__shared__ volatile cudaComplexType shAdpDer0[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shAdpDer1[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shAdpDer2[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shAdpDer3[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shAdpDer4[THREADSREDUCEKERNEL];
	__shared__ volatile cudaComplexType shAdpDer5[THREADSREDUCEKERNEL];

	cudaComplexType twoPiI = cudaComplexType(0, 2.0*M_PI);

	if (threadIdx.x < N) {
		int idx = N*blockIdx.x + threadIdx.x;
		shOccDer[threadIdx.x]  = reduceOccupancyDerivatives[idx];
		shPosDerX[threadIdx.x] = reducePositionDerivativesX[idx];
		shPosDerY[threadIdx.x] = reducePositionDerivativesY[idx];
		shPosDerZ[threadIdx.x] = reducePositionDerivativesZ[idx];
		shAdpDer0[threadIdx.x] = reduceAdpDerivative0[idx];
		shAdpDer1[threadIdx.x] = reduceAdpDerivative1[idx];
		shAdpDer2[threadIdx.x] = reduceAdpDerivative2[idx];
		shAdpDer3[threadIdx.x] = reduceAdpDerivative3[idx];
		shAdpDer4[threadIdx.x] = reduceAdpDerivative4[idx];
		shAdpDer5[threadIdx.x] = reduceAdpDerivative5[idx];
		for(int i = threadIdx.x + blockDim.x; i < N; i += blockDim.x) {
			int idx = N*blockIdx.x+i;
			shOccDer[threadIdx.x]  += reduceOccupancyDerivatives[idx];
			shPosDerX[threadIdx.x] += reducePositionDerivativesX[idx];
			shPosDerY[threadIdx.x] += reducePositionDerivativesY[idx];
			shPosDerZ[threadIdx.x] += reducePositionDerivativesZ[idx];
			shAdpDer0[threadIdx.x] += reduceAdpDerivative0[idx];
			shAdpDer1[threadIdx.x] += reduceAdpDerivative1[idx];
			shAdpDer2[threadIdx.x] += reduceAdpDerivative2[idx];
			shAdpDer3[threadIdx.x] += reduceAdpDerivative3[idx];
			shAdpDer4[threadIdx.x] += reduceAdpDerivative4[idx];
			shAdpDer5[threadIdx.x] += reduceAdpDerivative5[idx];
		}
	}
	else {
		shOccDer[threadIdx.x]  = cudaComplexType(0,0);
		shPosDerX[threadIdx.x] = cudaComplexType(0,0);
		shPosDerY[threadIdx.x] = cudaComplexType(0,0);
		shPosDerZ[threadIdx.x] = cudaComplexType(0,0);
		shAdpDer0[threadIdx.x] = cudaComplexType(0,0);
		shAdpDer1[threadIdx.x] = cudaComplexType(0,0);
		shAdpDer2[threadIdx.x] = cudaComplexType(0,0);
		shAdpDer3[threadIdx.x] = cudaComplexType(0,0);
		shAdpDer4[threadIdx.x] = cudaComplexType(0,0);
		shAdpDer5[threadIdx.x] = cudaComplexType(0,0);
	}
	__syncthreads();
	//TODO: tutaj trzeba pozmieniac (19.07.14 by szmigacz) albo nie ma sensu bo i po co?
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shOccDer  , occupancyDerivatives + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shPosDerX , positionDerivativesX + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shPosDerY , positionDerivativesY + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shPosDerZ , positionDerivativesZ + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shAdpDer0 , adpDerivatives0 + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shAdpDer1 , adpDerivatives1 + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shAdpDer2 , adpDerivatives2 + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shAdpDer3 , adpDerivatives3 + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shAdpDer4 , adpDerivatives4 + blockIdx.x);
	reduceArrayAndSave<THREADSREDUCEKERNEL>( shAdpDer5 , adpDerivatives5 + blockIdx.x);
}
//END-REDUCE DERIVATIVES KERNEL---------------------------------------------}}}
