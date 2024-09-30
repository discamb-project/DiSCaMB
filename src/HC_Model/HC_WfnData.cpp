#include "discamb/HC_Model/HC_WfnData.h"
#include "discamb/MathUtilities/math_utilities.h"

using namespace std;

namespace discamb {

namespace sto_atomic_wfn
{

void atomicStoWfnToSphericalDensity(
    const HC_WfnBankEntry &wfn,
    vector<double> &coefficients,
    vector<double> &exponents,
    vector<int> &powerR,
    ElectronDensityType densityType)
{
    int orbitalIdx,nOrbitals,i,j,nSlaters;
    double coefficient, occupancy;
    vector<int> orbitalIndices; 
    
    // clear containers for the results

    coefficients.clear();
    exponents.clear();
    powerR.clear();

    // figure out which orbitals to use depending on densityType
    if(densityType == TOTAL_DENSITY)
    {
        orbitalIndices = wfn.core_orbitals_indices;
        orbitalIndices.insert(orbitalIndices.end(),wfn.valence_orbitals_indices.begin(),
                                                   wfn.valence_orbitals_indices.end());
    }
    else
        densityType == CORE_DENSITY ? orbitalIndices = wfn.core_orbitals_indices : 
                              orbitalIndices = wfn.valence_orbitals_indices;
    
    // add contribution from orbitals 


    nOrbitals = orbitalIndices.size();
    for(int k = 0 ; k < nOrbitals ; k++)
    {
        orbitalIdx = orbitalIndices[k];
        occupancy = (double) wfn.orbital_occupancy[orbitalIdx];
        if( occupancy != 0.0 )
        {
            nSlaters = wfn.orbitals[orbitalIdx].coefficients.size();

            for( i = 0 ; i < nSlaters ; i++ )
                for( j = 0 ; j <= i ; j++ )
                {
                    coefficient = wfn.orbitals[orbitalIdx].coefficients[i] *
                                  wfn.orbitals[orbitalIdx].coefficients[j] * 
                                  occupancy / ( 4 * M_PI );
                    if(i != j)
                        coefficient *= 2;
                    
                    coefficients.push_back(coefficient);

                    exponents.push_back(wfn.orbitals[orbitalIdx].exponents[i] +
                                        wfn.orbitals[orbitalIdx].exponents[j]);

                    powerR.push_back(wfn.orbitals[orbitalIdx].power_r[i] + 
                                     wfn.orbitals[orbitalIdx].power_r[j]);
                }
        }
    }
}

} // namepsace sto_atomic_wfn
} //namespace discamb
