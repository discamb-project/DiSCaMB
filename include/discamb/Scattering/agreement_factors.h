#pragma once

#include "discamb/Scattering/SF_CalcDataTypes.h"

#include <vector>
#include <variant>

namespace discamb {
    namespace agreement_factors {

        enum AgreementFactor {Relative_L1_Percent, R2, R2_target};

        /**
        returns a scale factor for i_calc, the value optimizes wR2 agreement factor 
        */

        double scale_for_wR2(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights);

        double wR2(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights);

        double wR2(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights,
            double &scale);

        double wR2_predefined_scale(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights,
            double scale = 1.0);

        double value(
            const std::vector<double>& v_calc, 
            const std::vector<double>& v_ref,
            AgreementFactor af = Relative_L1_Percent);

        double value(
            const std::vector<double>& v_calc,
            const std::vector<double>& v_ref,
            const std::vector<double>& weights,
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
