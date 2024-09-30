#ifndef _DISCAMB_SCATTERING_SCATTERING_UTILITIES_H_
#define _DISCAMB_SCATTERING_SCATTERING_UTILITIES_H_

#include "discamb/MathUtilities/Vector3.h"
#include "SF_CalcDataTypes.h"
#include "SfCalculator.h"
#include "Real.h"

#include <vector>
#include <complex>
#include <memory>

namespace discamb {

    /**
    * \defgroup Scattering Scattering
    \brief Form factor and scattering factor calculation.
    * @{
    */




namespace scattering_utilities
{
    std::unique_ptr<SfCalculator> scatteringFactorCalculatorFromJsonFile(const Crystal& crystal, const std::string& jsonFile);

    void centeringSfMultipliers(char centering,const std::vector<Vector3i> &hkl, std::vector<double> &multipliers, bool obverse);

    template<typename T>
    void divideSet(const std::vector<T> &v,int nBatches,std::vector<std::vector<T> > &subsets);

    void combineScatteringFactorsSets(
        const std::vector< std::vector< std::complex< REAL > > > &scatteringFactorsSets,
        std::vector< std::complex< REAL > > &scatteringFactors);

    void merge_dTarget_dParameterSets(
        const std::vector<std::vector<TargetFunctionAtomicParamDerivatives> > &dTarget_dParamSets,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dParam);

    // 
    double scaleFactorFcalc(const std::vector<double>& i_obs, const std::vector<double>& i_calc, const std::vector<double>& i_sigma,
        double a = 0, double b = 0);

}


//----- TEMPLATE IMPLEMENTATION --------------

template<typename T>
void scattering_utilities::divideSet(
    const std::vector<T> &v,
    int nBatches, 
    std::vector<std::vector<T> > &subsets)
{
    
    int maxN_ElementsPerSubset, nElements = v.size();

    if (nElements % nBatches == 0)
        maxN_ElementsPerSubset = nElements / nBatches;
    else
        maxN_ElementsPerSubset = nElements / nBatches + 1;

    subsets.clear();
    subsets.resize(nBatches);


    for (int batchIdx = 0; batchIdx<nBatches; batchIdx++)
    {
        if (batchIdx != nBatches - 1)
            subsets[batchIdx].insert(subsets[batchIdx].end(),
                                      v.begin() + batchIdx * maxN_ElementsPerSubset,
                                      v.begin() + (batchIdx + 1) * maxN_ElementsPerSubset);
        else
            subsets[batchIdx].insert(subsets[batchIdx].end(),
                                      v.begin() + batchIdx * maxN_ElementsPerSubset,
                                      v.end());
    }
    

}
/** @}*/
} // namespace discamb

#endif /*_DISCAMB_SCATTERING_SCATTERING_UTILITIES_H_*/
