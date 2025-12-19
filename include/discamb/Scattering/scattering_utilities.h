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
    void generate_symmetry_equivalent_hkls(
        const SpaceGroup& spaceGroup,
        const std::vector<Vector3i>& hkl,
        std::vector<Vector3i>& symmEquivalentHkls);

    std::unique_ptr<SfCalculator> scatteringFactorCalculatorFromJsonFile(const Crystal& crystal, const std::string& jsonFile);

    void centeringSfMultipliers(char centering,const std::vector<Vector3i> &hkl, std::vector<double> &multipliers, bool obverse);

    template<typename T>
    void divideSet(const std::vector<T> &v,int nBatches,std::vector<std::vector<T> > &subsets);

    template<typename T>
    void divideSet(const std::vector<T>& v, int nBatches, std::vector<std::vector<T> >& subsets, std::vector<std::vector<int> > &subset_indices);

    void combineScatteringFactorsSets(
        const std::vector< std::vector< std::complex< REAL > > > &scatteringFactorsSets,
        std::vector< std::complex< REAL > > &scatteringFactors);

    void merge_dTarget_dParameterSets(
        const std::vector<std::vector<TargetFunctionAtomicParamDerivatives> > &dTarget_dParamSets,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dParam);

    // 
    double scaleFactorFcalc(const std::vector<double>& i_obs, const std::vector<double>& i_calc, const std::vector<double>& i_sigma,
        double a = 0, double b = 0);

    int findPreferredHklOrderingDirection(
        const std::vector<Vector3i>& hkl,
        std::vector<std::vector<Vector3i> >& orderedHklLines,
        std::vector<std::vector<int> >& mapToOriginalSetIndices);



    // splits hkl groups into sets with similar number of hkl
    // https://en.m.wikipedia.org/wiki/Multiway_number_partitioning
    void splitHklLines(
        int nSets,
        const std::vector<std::vector<Vector3i> >& orderedHklLines,
        std::vector<std::vector<int> >& lineGroups);
    
    // 
    void splitHklLines(
        int nSets,
        const std::vector<std::vector<Vector3i> >& orderedHklLines,
        const std::vector<std::vector<int> >& hklIdxInOriginalSet,
        std::vector< std::vector<std::vector<Vector3i> > > & newHklLines,
        std::vector < std::vector <std::vector<std::pair<int, int> > > >& hklMap,
        std::vector < std::vector<std::vector<int> > >& subsetDataHklIdxInOryginalSet);

    void init_line_multipliers(
        const Vector3d& lineStepCart,
        const std::vector< std::vector<Vector3d> >& r_atom_symm,
        const std::vector<std::vector<double> > &adps,
        const std::vector<Matrix3d>& symmOpRotationCart, 
        std::vector<std::vector<std::complex<double> > > &phase_factor_multiplier,
        std::vector<std::vector<double> > &temp_factor_multiplier_multiplier);

    void calculate_line_temperature_factors(
        //in:
        const std::vector<Vector3d>& coordinates,
        const std::vector< std::vector<double> >& atoms_adps,
        const std::vector<Vector3d>& hkl_cart,
        int hkl_idx,
        const std::vector<Vector3i>& hkl_line,
        int idx_in_line,
        int line_direction,
        const std::vector<Vector3<double> >& rotated_h,
        std::vector< std::vector<double> >& temp_factor_multiplier,
        const std::vector< std::vector<double> >& temp_factor_multiplier_multiplier,
        std::vector<std::vector<double> >& adpMultipliers,
        //out:
        std::vector< std::vector<double> >& line_temperature_factors);

    void calculate_line_phase_factors(
        //in:
        const std::vector<std::vector<Vector3d> >& r_at,
        const std::vector<Vector3d>& hkl_cart,
        int hkl_idx,
        const std::vector<Vector3i>& hkl_line,
        int idx_in_line,
        int line_direction,
        const std::vector<Vector3<double> >& rotated_h,
        const std::vector< std::vector<std::complex<double > > >& phase_factor_multiplier,
        //std::vector<Vector3d> &rotated_h_2pi,
        //std::vector<double> &translation_factor_2pi,

        //out:
        std::vector< std::vector<std::complex<double > > >& line_phase_factor);

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

template<typename T>
void scattering_utilities::divideSet(
    const std::vector<T>& v,
    int nBatches,
    std::vector<std::vector<T> >& subsets,
    std::vector<std::vector<int> >& subset_indices)
{
    int maxN_ElementsPerSubset, nElements = v.size();

    if (nElements % nBatches == 0)
        maxN_ElementsPerSubset = nElements / nBatches;
    else
        maxN_ElementsPerSubset = nElements / nBatches + 1;

    subsets.clear();
    subsets.resize(nBatches);
    subset_indices.resize(nBatches);


    for (int batchIdx = 0; batchIdx < nBatches; batchIdx++)
    {
        if (batchIdx != nBatches - 1)
        {
            subsets[batchIdx].insert(subsets[batchIdx].end(),
                v.begin() + batchIdx * maxN_ElementsPerSubset,
                v.begin() + (batchIdx + 1) * maxN_ElementsPerSubset);
            for (int i = 0; i < maxN_ElementsPerSubset; i++)
                subset_indices[batchIdx].push_back(i + batchIdx * maxN_ElementsPerSubset);
        }
        else
        {
            subsets[batchIdx].insert(subsets[batchIdx].end(),
                v.begin() + batchIdx * maxN_ElementsPerSubset,
                v.end());
            int n = nElements - (nBatches-1) * maxN_ElementsPerSubset;
            for (int i = 0; i < n; i++)
                subset_indices[batchIdx].push_back(i + batchIdx * maxN_ElementsPerSubset);
        }
    }
}



/** @}*/
} // namespace discamb

#endif /*_DISCAMB_SCATTERING_SCATTERING_UTILITIES_H_*/
