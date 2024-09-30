#include "discamb/MathUtilities/Vector3.h"

#include <string>
#include <vector>
#include <complex>

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


    namespace tsc_io {
        void read_tsc(
            const std::string& fileName,
            std::vector<std::string>& atomLabels,
            std::vector<Vector3i>& hkl,
            //[hklIdx][atomIdx]
            std::vector<std::vector<std::complex<double> > >& atomicFormFactors);

        void read_tsc2(
            const std::string& fileName,
            std::vector<std::string>& atomLabels,
            std::vector<Vector3i>& hkl,
            //[hklIdx][atomIdx]
            std::vector<std::vector<std::complex<double> > >& atomicFormFactors);


        void write_tsc(
            const std::string& fileName,
            const std::vector<std::string>& atomLabels,
            const std::vector<Vector3i>& hkl,
            //[hklIdx][atomIdx]
            const std::vector<std::vector<std::complex<double> > >& atomicFormFactors,
            const std::string &additionalText = std::string());

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

