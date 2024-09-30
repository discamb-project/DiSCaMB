// vim: noet: sw=3: ts=3
#ifndef CUDACONSTANT_CUH_ELPVRG41
#define CUDACONSTANT_CUH_ELPVRG41

#include "types.h"
#include "globals.h"

__constant__ real constPLM[MAXATOMTYPES*SPHHARMCOEF];
__constant__ int constDefValencePow[MAXWFNPARAMS * 5];
__constant__ real constDefValenceExp[MAXWFNPARAMS];

__constant__ real3 constRotationsRow1[MAXSYMMETRYOPERATIONS];
__constant__ real3 constRotationsRow2[MAXSYMMETRYOPERATIONS];
__constant__ real3 constRotationsRow3[MAXSYMMETRYOPERATIONS];
__constant__ real3 constTranslations[MAXSYMMETRYOPERATIONS];

__constant__ double constInverseFactorial[11];

#endif /* end of include guard: CUDACONSTANT_CUH_ELPVRG41 */

