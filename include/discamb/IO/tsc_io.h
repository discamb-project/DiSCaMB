#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <string>
#include <vector>
#include <complex>

namespace discamb {

    /**
     * \addtogroup IO
     * @{
     */

    /**
     * \brief Utilities for reading and writing TSC (textual atomic form-factor) files.
     *
     * The functions in this namespace provide simple read/write helpers for TSC-like
     * files where atomic form factors are stored per reflection (h,k,l) and per atom.
     */
    namespace tsc_io {

        /**
         * @brief Read a TSC file with integer Miller indices.
         *
         * Reads atom labels, Miller indices (h,k,l) and per-atom complex form-factors.
         * The form-factors matrix is organised as [hkl_index][atom_index].
         *
         * @param[in] fileName Path to the input TSC file.
         * @param[out] atomLabels Vector that will receive atom labels (column order).
         * @param[out] hkl Vector that will receive integer Miller indices (Vector3i).
         * @param[out] atomicFormFactors Matrix of atomic form factors indexed as [hklIdx][atomIdx].
         */
        void read_tsc(
            const std::string& fileName,
            std::vector<std::string>& atomLabels,
            std::vector<Vector3i>& hkl,
            //[hklIdx][atomIdx]
            std::vector<std::vector<std::complex<double> > >& atomicFormFactors);

        /**
         * @brief Alternative TSC reader (variant 2). 
         * Behaviour is equivalent to \ref read_tsc but speed might be different.
         *
         * @param[in] fileName Path to the input TSC file.
         * @param[out] atomLabels Vector that will receive atom labels (column order).
         * @param[out] hkl Vector that will receive integer Miller indices (Vector3i).
         * @param[out] atomicFormFactors Matrix of atomic form factors indexed as [hklIdx][atomIdx].
         */
        void read_tsc2(
            const std::string& fileName,
            std::vector<std::string>& atomLabels,
            std::vector<Vector3i>& hkl,
            //[hklIdx][atomIdx]
            std::vector<std::vector<std::complex<double> > >& atomicFormFactors);


        /**
         * @brief Write TSC file with integer Miller indices.
         *
         * Writes atom labels, Miller indices and per-atom complex form-factors.
         * An optional \p additionalText can be appended to the file, e.g. information on
         * form factors model, probably should be of format key: value, e.g.:
         * 
         * FORM FACTORS SOURCE:
         *     SCATTERING MODEL - TAAM
         *     TAAM DATABANK - MATTS (Kumar, P., Gruza, B., Bojarowski, S. A. & Dominiak, P. M. (2019). Acta Cryst. A75, 398-408.)
         *
         * @param[in] fileName Path to the output TSC file.
         * @param[in] atomLabels Atom labels corresponding to columns.
         * @param[in] hkl Vector of integer Miller indices (Vector3i).
         * @param[in] atomicFormFactors Matrix of atomic form factors indexed as [hklIdx][atomIdx].
         * @param[in] additionalText Optional additional text to put be befor DATA section in the file.
         */
        void write_tsc(
            const std::string& fileName,
            const std::vector<std::string>& atomLabels,
            const std::vector<Vector3i>& hkl,
            //[hklIdx][atomIdx]
            const std::vector<std::vector<std::complex<double> > >& atomicFormFactors,
            const std::string &additionalText = std::string());

        /**
         * @brief Write TSC(D) file with form factors at reciprocal space 
         * points with real-valued components (i.e. not necessarily integer 
         * Miller indices).
         *
         * Similar to \ref write_tsc but accepts 'Miller indices' with non-integers (as Vector3d).
         *
         * @param[in] fileName Path to the output TSC(D) file.
         * @param[in] atomLabels Atom labels corresponding to columns.
         * @param[in] hkl Vector of Miller indices represented as Vector3d.
         * @param[in] atomicFormFactors Matrix of atomic form factors indexed as [hklIdx][atomIdx].
         * @param[in] additionalText Optional additional text to include in the file.
         */
        void write_tscd(
            const std::string& fileName,
            const std::vector<std::string>& atomLabels,
            const std::vector<Vector3d>& hkl,
            //[hklIdx][atomIdx]
            const std::vector<std::vector<std::complex<double> > >& atomicFormFactors,
            const std::string& additionalText = std::string());


    }
    /**@}*/
}

