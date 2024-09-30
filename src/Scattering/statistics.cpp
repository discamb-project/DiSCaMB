#include "discamb/Scattering/statistics.h"
#include "discamb/BasicUtilities/OnError.h"

using namespace std;

namespace discamb {

	// from test_utilities.h in discamb

	namespace fit_statistics {

		double rFactor(
			const std::vector<std::complex<double> > &data,
			const std::vector<std::complex<double> > &referenceData,
			bool scaleFactor)
		{
			double oo, co, cc;
			int i, n = data.size();
			vector<complex<double> > diff;

			if (data.size() != referenceData.size())
				on_error::throwException("an attempt to calculate R-factors between data sets of differtent size", __FILE__, __LINE__);

			if (data.empty())
				return -10;

			oo = 0;
			co = 0;
			cc = 0;

			diff.resize(n);

			for (i = 0; i < n; i++)
			{
				oo += std::norm(referenceData[i]);
				co += (referenceData[i] * conj(data[i])).real();
				cc += std::norm(data[i]);
				diff[i] = referenceData[i] - data[i];
			}

			if (oo == 0)
				return -5;
			//on_error::throwException("denominator 0 in R-factor", __FILE__, __LINE__);

			double numerator, denominator, s;

			scaleFactor ? s = co / cc : s = 1.0;

			numerator = oo - 2 * s*co + s * s*cc;
			denominator = oo;

			if (denominator < 1e-10)
				return -5;
			if (referenceData.empty())
				return -10;

			return sqrt(abs(numerator / denominator));

		}

		double rFactor(
			const std::vector<double> &data,
			const std::vector<double> &referenceData,
			bool scale)
		{
			double oo;//, co, cc;
			int i, n = data.size();
			double numerator, denominator;


			if (data.size() != referenceData.size())
				on_error::throwException("an attempt to calculate R-factors between data sets of differtent size", __FILE__, __LINE__);


			oo = 0;

			if (scale)
			{
				double co, cc, s;

				co = 0;
				cc = 0;

				for (i = 0; i < n; i++)
				{
					oo += referenceData[i] * referenceData[i];
					co += data[i] * referenceData[i];
					cc += data[i] * data[i];
				}
				s = s = co / cc;
				numerator = oo - 2 * s*co + s * s*cc;
			}
			else
			{
				double diff, diff2 = 0;
				for (i = 0; i < n; i++)
				{
					diff = data[i] - referenceData[i];
					diff2 += diff * diff;
					oo += referenceData[i] * referenceData[i];
				}
				numerator = diff2;
			}

			if (oo == 0)
				return -5;
			//on_error::throwException("denominator 0 in R-factor", __FILE__, __LINE__); 

			denominator = oo;

			if (denominator < 1e-10)
				return -5;
			if (referenceData.empty())
				return -10;

			return sqrt(abs(numerator / denominator));
		}

        double scaleFactor(
            const std::vector<double>& data,
            const std::vector<double>& referenceData)
        {
            double co, cc;
          
            if (data.size() != referenceData.size())
            {
                on_error::throwException("an attempt to calculate scale factor between data sets of differtent size", __FILE__, __LINE__);
                return 0.0;
            }

            if (data.empty())
            {
                on_error::throwException("an attempt to calculate scale factor for empty set of intensities", __FILE__, __LINE__);
                return 0.0;
            }

            co = 0;
            cc = 0;

            for (int i = 0; i < data.size(); i++)
            {
                co += (referenceData[i] * conj(data[i])).real();
                cc += std::norm(data[i]);
            }

            if (cc == 0)
            {
                on_error::throwException("an attempt to calculate scale factor for 0 intensity data", __FILE__, __LINE__);
                return 0.0;
            }

            return co / cc;

        }

	}
}

