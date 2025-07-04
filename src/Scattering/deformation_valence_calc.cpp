#include "discamb/Scattering/deformation_valence_calc.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/MathUtilities/math_utilities.h"


namespace discamb {
    namespace deformation_valence_calc {

        deformation_valence_evaluator getDeformationValenceEvaluator(
            const std::string& point_grouop)
        {
            if (point_grouop == "1")
                return &calculateDefVal_pg_1;
            else if (point_grouop == "2")
                return &calculateDefVal_pg_2;
            else
                return &calculateDefVal_general;
        }

        void calculateDefVal_general(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& g_functions_and_slater_normalization,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
            std::vector<std::complex<double> >& defVal)
        {
            on_error::not_implemented(__FILE__, __LINE__);

            //int nSymmetryOps = sphericalHarmonics.size();
            //for (int i = 0; i < nSymmetryOps; i++)
            //    defVal[i] = HansenCoppens_SF_Engine4::calculateDeformationValence(
            //        p_lm, g_functions_and_slater_normalization, maxL, sphericalHarmonics[i]);
        }

        void calculateDefVal_pg_1(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& g_functions_and_slater_normalization,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
            std::vector<std::complex<double> >& defVal)
        { 
            on_error::not_implemented(__FILE__, __LINE__);
            //defVal[0] = HansenCoppens_SF_Engine4::calculateDeformationValence(
              //  p_lm, g_functions_and_slater_normalization, maxL, sphericalHarmonics[0]);
        }

        void calculateDefVal_pg_2(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& radial,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& _sphericalHarmonics,
            std::vector<std::complex<double> >& defVal)
        {
            //{"2", {"x,y,x", "-x,y,-z"}},
            double a_re = 0;
            double a_im = 0;
            double b_re = 0;
            double b_im = 0;
            auto const& sph = _sphericalHarmonics[0];

            a_re += radial[0] * p_lm[0][0] * sph[0][0];

            if (maxL > 0)
            {
                a_im += radial[1] * p_lm[1][0] * sph[1][0];
                b_im += radial[1] * (p_lm[1][1] * sph[1][1]+ p_lm[1][2] * sph[1][2]);


                if (maxL > 1)
                {
                    a_re -= radial[2] * (
                        p_lm[2][2] * sph[2][2] +
                        p_lm[2][3] * sph[2][3] +
                        p_lm[2][4] * sph[2][4]);

                    b_re -= radial[2] * (p_lm[2][0] * sph[2][0] +
                        p_lm[2][1] * sph[2][1]);

                    if (maxL > 2)
                    {
                        a_im -= radial[3] * (
                            p_lm[3][1] * sph[3][1] +
                            p_lm[3][3] * sph[3][3] +
                            p_lm[3][5] * sph[3][5] );
                        b_im -= radial[3] * (p_lm[3][0] * sph[3][0] +
                            p_lm[3][2] * sph[3][2] +
                            p_lm[3][4] * sph[3][4] + 
                            p_lm[3][6] * sph[3][6]);

                        if (maxL > 3)
                        {
                            a_re += radial[4] * (
                                p_lm[4][4] * sph[4][4] +
                                p_lm[4][5] * sph[4][5] +
                                p_lm[4][6] * sph[4][6] +
                                p_lm[4][7] * sph[4][7] +
                                p_lm[4][8] * sph[4][8]);
                            b_re += radial[4] * (p_lm[4][0] * sph[4][0] +
                                p_lm[4][1] * sph[4][1] +
                                p_lm[4][2] * sph[4][2] +
                                p_lm[4][3] * sph[4][3]);
                        }
                    }

                }

            }
            defVal[0] = 4 * M_PI * std::complex<double>(a_re + b_re, a_im + b_im);
            defVal[1] = 4 * M_PI * std::complex<double>(a_re - b_re, a_im - b_im);
        }


    }
}

