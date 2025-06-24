#pragma once

#include "discamb/Scattering/SF_CalcDataTypes.h"

#include <vector>

namespace discamb {
    namespace agreement_factors {

        enum AgreementFactor {Relative_L1_Percent};

        double value(
            const std::vector<double>& v1, 
            const std::vector<double>& v2,
            AgreementFactor af = Relative_L1_Percent);

        double value(
            const std::vector<std::complex<double> > & v1,
            const std::vector<std::complex<double> > & v2,
            AgreementFactor af = Relative_L1_Percent);


        void for_derivatives(
            const std::vector<TargetFunctionAtomicParamDerivatives>& dT_dp1,
            const std::vector<TargetFunctionAtomicParamDerivatives>& dT_dp2,
            bool& size_match,
            double& d_xyz_agreement_factor,
            double& d_adp_agreement_factor,
            double& d_occ_agreement_factor,
            double& d_fpfdp_agreement_factor,
            AgreementFactor af = Relative_L1_Percent);

        // ff1,ff2[hkl_idx][atom_idx] - form factors
        void for_atomic_form_factors(
            std::vector<std::vector<std::complex<double> > > &ff1,
            std::vector<std::vector<std::complex<double> > >& ff2,
            std::vector<double> &agreement_factors,
            AgreementFactor af = Relative_L1_Percent);

    }
}
