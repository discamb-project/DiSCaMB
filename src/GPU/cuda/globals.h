// vim: noet: sw=3: ts=3
#ifndef GLOBALS_H_VABG8MLH
#define GLOBALS_H_VABG8MLH

#define MAXATOMTYPES 250
#define MAXWFNPARAMS 10
#define MAXSYMMETRYOPERATIONS 48
#define N_ADP_DERIVATIVES 6

// maksymalna liczba L przy petli for (l=0 to l=L) for (m=-l to m=l) //tu nie ma buga
#define MAX_L 6 
//liczba watkow w glownym kernelu (dowolna wielokrotnosc 32)
#define THREADSMAINKERNEL 64
//rozmiar cudowego warpa (nie zmieniac!!)
#define WARPSIZE 32
//liczba wspolczynnikow do harmonik sferycznych
// == 1 + 3 + 5 + 7 + 9 = 25
#define SPHHARMCOEF 25
//liczba watkow wykorzystywanych w kernelu do redukcji
#define THREADSREDUCEKERNEL 128

#define cudaSafeCall(x) do { if((x)!=cudaSuccess) { \
	printf("Error at %s:%d\n",__FILE__,__LINE__);\
	exit(EXIT_FAILURE);}} while(0)

#endif /* end of include guard: GLOBALS_H_VABG8MLH */

