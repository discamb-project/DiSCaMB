#pragma once

#include<vector>
#include<complex>

namespace discamb {

    namespace deformation_valence_calc {
        /* Function pointer type for deformation valence evaluation
        defVal[i] - deformation valence for i-th symmetry operation,
        order of symmetry operations follow the order in the 
        crystallographic_point_group_tables::getPointGroupOperationsStr
        are defined in the symmetry_operations vector

        */
        typedef void (*deformation_valence_evaluator)(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& g_functions_and_slater_normalization,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
            std::vector<std::complex<double> >& defVal);

        // Function to get the appropriate deformation valence evaluator based on point group
        deformation_valence_evaluator getDeformationValenceEvaluator(const std::string& point_grouop);

        void calculateDefVal_general(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& g_functions_and_slater_normalization,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
            std::vector<std::complex<double> >& defVal);

        void calculateDefVal_pg_1(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& g_functions_and_slater_normalization,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
            std::vector<std::complex<double> >& defVal);

        void calculateDefVal_pg_2(
            const std::vector<std::vector<double> >& p_lm, // coefficients for multipolar terms (with wavefunction normalization of spherical harmonics)
            const std::vector<double>& g_functions_and_slater_normalization,
            int maxL,
            const std::vector< std::vector <std::vector<double> > >& sphericalHarmonics,
            std::vector<std::complex<double> >& defVal);

    }
}

