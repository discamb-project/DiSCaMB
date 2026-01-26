#include "discamb/Scattering/agreement_factors.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb {
    namespace agreement_factors {

        namespace {
            template<typename T>
            double value_relative_L1_percent(
                const std::vector<T>& v1,
                const std::vector<T>& v2,
                const std::vector<double> &weights)
            {
                double numerator, denominator;
                int i, n = v1.size();
                if (n != v2.size())
                    on_error::throwException("vectors size do not match when trying to calculate (w)L1 agreement factor", __FILE__, __LINE__);
                if (n != weights.size())
                    on_error::throwException("weights size do not match when trying to calculate (w)L1 agreement factor", __FILE__, __LINE__);
                numerator = 0.0;
                denominator = 0.0;

                for (i = 0; i < n; i++)
                {
                    numerator += weights[i] * abs(v1[i] - v2[i]);
                    denominator += weights[i] * (abs(v1[i]) + abs(v2[i]));
                }
                if (denominator == 0)
                    on_error::throwException("denominator 0 when trying to calculate L1 agreement factor", __FILE__, __LINE__);

                return numerator / denominator * 200;
            }

            template<typename T>
            double value_r2_target(
                const std::vector<T>& v_calc,
                const std::vector<T>& v_ref,
                const std::vector<double>& weights)
            {
                
                int i, n = v_calc.size();
                if (n != v_ref.size())
                    on_error::throwException("vectors size do not match when trying to calculate (w)R2 agreement factor", __FILE__, __LINE__);
                if (n != weights.size())
                    on_error::throwException("weights size do not match when trying to calculate (w)R2 agreement factor", __FILE__, __LINE__);
                double result = 0.0;

                for (i = 0; i < n; i++)
                {
                    double diff = v_calc[i] - v_ref[i];
                    result += weights[i] * norm(diff);
                }

                return result;
            }

            template<typename T>
            double value_r2(
                const std::vector<T>& v_calc,
                const std::vector<T>& v_ref,
                const std::vector<double>& weights)
            {
                double numerator, denominator;
                int i, n = v_calc.size();
                if (n != v_ref.size())
                    on_error::throwException("vectors size do not match when trying to calculate (w)R2 agreement factor", __FILE__, __LINE__);
                if (n != weights.size())
                    on_error::throwException("weights size do not match when trying to calculate (w)R2 agreement factor", __FILE__, __LINE__);
                numerator = value_r2_target(v_calc, v_ref, weights);
                denominator = 0.0;

                for (i = 0; i < n; i++)
                {
                    double diff = v_calc[i] - v_ref[i];
                    denominator += weights[i] * v_ref[i] * v_ref[i];
                }
                if (denominator == 0)
                    on_error::throwException("denominator 0 when trying to calculate R2 agreement factor", __FILE__, __LINE__);

                return sqrt(numerator / denominator) * 100;
            }


        }
        
        double scale_for_wR2(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights)
        {
            double scale = 1.0;
            double numerator = 0.0;
            double denominator = 0.0;

            for (int i = 0; i < i_obs.size(); i++)
            {
                numerator += weights[i] * i_calc[i] * i_obs[i];
                denominator += weights[i] * i_calc[i] * i_calc[i];
            }

            scale = numerator / denominator;
            return scale;
        }

        double wR2(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights)
        {
            double scale;
            return wR2(i_obs, i_calc, weights, scale);
        }

        double wR2(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights,
            double& scale)
        {
            scale = scale_for_wR2(i_obs, i_calc, weights);
            return wR2_predefined_scale(i_obs, i_calc, weights, scale);
        }

        double wR2_predefined_scale(
            const std::vector<double>& i_obs,
            const std::vector<double>& i_calc,
            const std::vector<double>& weights,
            double scale)
        {
            double numerator = 0.0;
            double denominator = 0.0;

            for (int i = 0; i < i_obs.size(); i++)
            {
                double diff = i_obs[i] - i_calc[i] * scale;
                numerator += weights[i] * diff * diff;
                denominator += weights[i] * i_obs[i] * i_obs[i];
            }

            if (numerator < 0)
                on_error::throwException("invalid numerator when calculating wR2 agreement factor", __FILE__, __LINE__);

            if (denominator <= 0)
                on_error::throwException("invalid denominator when calculating wR2 agreement factor", __FILE__, __LINE__);

            return sqrt(numerator / denominator);

        }


        double value(
            const std::vector<double>& v1, 
            const std::vector<double>& v2,
            AgreementFactor af)
        {
            vector<double> weights(v1.size(), 1.0);
            if(af == Relative_L1_Percent)
                return value_relative_L1_percent(v1, v2, weights);
            else if(af == R2)
                return value_r2(v1, v2, weights);
            else
                on_error::throwException("unknown agreement factor type", __FILE__, __LINE__);
            return -1.0;
        }
        double value(
            const std::vector<double>& v1,
            const std::vector<double>& v2,
            const std::vector<double> &weights,
            AgreementFactor af)
        {
            if (af == Relative_L1_Percent)
                return value_relative_L1_percent(v1, v2, weights);
            else if (af == R2)
                return value_r2(v1, v2, weights);
            else
                on_error::throwException("unknown agreement factor type", __FILE__, __LINE__);
            return -1.0;
        }

        double value(
            const std::vector<std::complex<double> >& v1,
            const std::vector<std::complex<double> >& v2,
            AgreementFactor af)
        {
            vector<double> weights(v1.size(), 1.0);
            return value_relative_L1_percent(v1, v2, weights);
        }


        void for_derivatives(
            const std::vector<TargetFunctionAtomicParamDerivatives>& dT_dp1,
            const std::vector<TargetFunctionAtomicParamDerivatives>& dT_dp2,
            bool& size_match,
            double& d_xyz_agreement_factor,
            double& d_adp_agreement_factor,
            double& d_occ_agreement_factor,
            double& d_fpfdp_agreement_factor,
            AgreementFactor af)
        {
            int atomIdx, nAtoms = dT_dp1.size();
            size_match = (nAtoms == dT_dp2.size());
            vector<double> xyz1, xyz2, adps1, adps2, occ1, occ2;
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                if (dT_dp1[atomIdx].adp_derivatives.size() != dT_dp2[atomIdx].adp_derivatives.size())
                {
                    size_match = false;
                    return;
                }
                for (int i = 0; i < 3; i++)
                {
                    xyz1.push_back(dT_dp1[atomIdx].atomic_position_derivatives[i]);
                    xyz2.push_back(dT_dp2[atomIdx].atomic_position_derivatives[i]);
                }
                for (int i = 0; i < dT_dp1[atomIdx].adp_derivatives.size(); i++)
                {
                    adps1.push_back(dT_dp1[atomIdx].adp_derivatives[i]);
                    adps2.push_back(dT_dp2[atomIdx].adp_derivatives[i]);
                }
                occ1.push_back(dT_dp1[atomIdx].occupancy_derivatives);
                occ2.push_back(dT_dp2[atomIdx].occupancy_derivatives);
            }

            d_xyz_agreement_factor = agreement_factors::value(xyz1, xyz2, af);
            d_adp_agreement_factor = agreement_factors::value(adps1, adps2, af);
            d_occ_agreement_factor = agreement_factors::value(occ1, occ2, af);

        }

        void for_atomic_form_factors(
            std::vector<std::vector<std::complex<double> > >& ff1,
            std::vector<std::vector<std::complex<double> > >& ff2,
            std::vector<double>& agreement_factors,
            AgreementFactor af)
        {
            int nHkl = ff1.size();
            if (nHkl != ff2.size())
                on_error::throwException("cannot calculate agreement factors for atomic form factors - array size mismatch", __FILE__, __LINE__);
            
            if (nHkl == 0)
            {
                agreement_factors.clear();
                return;
            }
            
            int nAtoms = ff1[0].size();
            if(nAtoms != ff2[0].size())
                on_error::throwException("cannot calculate agreement factors for atomic form factors - number of atoms does not match", __FILE__, __LINE__);
            
            agreement_factors.resize(nAtoms);
            vector<complex<double> > v1(nHkl), v2(nHkl);

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
                {
                    v1[hklIdx] = ff1[hklIdx][atomIdx];
                    v2[hklIdx] = ff2[hklIdx][atomIdx];
                }
                agreement_factors[atomIdx] = agreement_factors::value(v1, v2, af);
            }
        }


    }
}
