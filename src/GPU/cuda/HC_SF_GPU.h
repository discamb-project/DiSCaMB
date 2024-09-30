#ifndef _DISCAMB_GPU_CUDA_HC_SF_GPU_
#define _DISCAMB_GPU_CUDA_HC_SF_GPU_

#include <string>
#include <sstream>
#include <vector>


namespace hc_sf_gpu{
    void gpuInfo(bool &hasGpu, size_t &memory, std::string &gpuName);
    size_t gpuOnlyMemory(size_t nHkl, size_t nAtoms);
    
    void calculate_multipolar_structure_factors_gpu(
        //######  I N P U T 
        std::stringstream &in,
        const std::vector<double> &dTarget_dF_real,
        const std::vector<double> &dTarget_dF_imag,
        //######  O U T P U T
        std::vector<double> &structureFactors_real,
        std::vector<double> &structureFactors_imag,
        std::vector<double> &xyzDerivatives,
        std::vector<double> &adpDerivatives,
        std::vector<double> &occupancyDerivatives
    );

    void calculate_multipolar_structure_factors_gpu2(
    //######  I N P U T 
  
        //---- SYMM OPS ---
    /* individual symmOp as list m:
                  |m0 m1 m2| | m9 |
    { M | t } = { |m3 m4 m5||| m10| }
                  |m6 m7 m8| | m11|    */
    const std::vector<std::vector<double> > &symmOps,
    //--- WFN TYPES -----    
    // core 
    const std::vector<std::vector<double> > &core_coeff, 
    const std::vector<std::vector<double> > &core_exp,
    const std::vector<std::vector<size_t> > &core_pow,
    // spherical valence
    const std::vector<std::vector<double> > &sph_val_coeff, 
    const std::vector<std::vector<double> > &sph_val_exp,
    const std::vector<std::vector<size_t> > &sph_val_pow,
    // deformation valence
    const std::vector<double> &def_val_exponent,
    const std::vector<std::vector<size_t> > &def_val_n_l,
    //--- ATOM TYPES ----------
    const std::vector<double> &p_val,
    const std::vector<double> &kappa_spherical,
    const std::vector<double> &kappa_def_valence,
    const std::vector<std::vector<std::vector<double> > > &p_lm,
    //--- ATOMIC PARAMETERS ---
    const std::vector<int> &atomToWfnTypeMap,
    const std::vector<int> &atomToAtomTypeMap,
    const std::vector<std::vector<double> > &atomicPositions,
    const std::vector<std::vector<double> > &adps,
    const std::vector<double> &occupancy,
    const std::vector<double> &multiplicity_weigths,
    const std::vector<std::vector<std::vector<double> > > &lcs,
    //--- HKL dTarget_dF---
    const std::vector<std::vector<double> > &hkl,
    const std::vector<double> &dTarget_dF_real,
    const std::vector<double> &dTarget_dF_imag,
    //######  O U T P U T
    std::vector<double> &structureFactors_real,
    std::vector<double> &structureFactors_imag,
    std::vector<double> &xyzDerivatives,
    std::vector<double> &adpDerivatives,
    std::vector<double> &occupancyDerivatives);

    
}

#endif