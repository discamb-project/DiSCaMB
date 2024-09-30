#pragma once

#include <vector>
#include <complex>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


	// from test_utilities.h in discamb

	namespace fit_statistics {

		/**
		denominator<1e-10 returns -5, referenceData.empty() - returns -10
		*/

		double rFactor(const std::vector<std::complex<double> > &scattering_factors,
			const std::vector<std::complex<double> > &reference_scattering_factors,
			bool scale);

		double rFactor(const std::vector<double> &data,
			const std::vector<double> &referenceData,
			bool scale);
        /**
        minimizes (referenceData - s*data)^2
        */
        double scaleFactor(const std::vector<double>& data,
            const std::vector<double>& referenceData);
	}
    /** @}*/
}

