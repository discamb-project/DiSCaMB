#include "discamb/HC_Model/hc_model_utilities.h"
#include "discamb/MathUtilities/SphConverter.h"
#include "discamb/MathUtilities/RealSphericalHarmonics.h"

using namespace std;

namespace discamb{

	namespace hc_model_utilities {

		void rotate_plm(
			std::vector<std::vector<double> > &plm,
			const std::vector<Vector3d> &oldXYZ,
			const std::vector<Vector3d> &newXYZ)
		{
			SphConverter converter;
			vector<vector<double> > rotated_plm, plm_wfn, rotated_plm_wfn;
			vector<vector<double> > lcsOld(3, vector<double>(3)), lcsNew(3, vector<double>(3));
			vector<vector<vector<double> > > rotation;
			int i, j, n, l_max, l;

			if (plm.empty())
				return;

			l_max = plm.size() - 1;

			converter.setMaxL(l_max);

			rotation.resize(l_max + 1);

			for (i = 0; i <= l_max; i++)
				rotation[i].resize(2 * i + 1);

			//columns of local_coordinate_systems corresponds to basis vectors
			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					lcsOld[i][j] = oldXYZ[i][j];

			for (i = 0; i < 3; i++)
				for (j = 0; j < 3; j++)
					lcsNew[i][j] = newXYZ[i][j];

			//

			plm_wfn = plm;

			for (l = 0; l <= l_max; l++)
			{
				n = 2 * l + 1;
				for (i = 0; i < n; i++)
					plm_wfn[l][i] /= real_spherical_harmonics::wfnNormalizationMultipliers[l][i] /
								     real_spherical_harmonics::densityNormalizationMultipliers[l][i];
			}

			//rotated_plm_wfn

			converter.convert(lcsOld, lcsNew, rotation);

			rotated_plm_wfn.resize(l_max + 1);

			for (l = 0; l <= l_max; l++)
			{
				n = 2 * l + 1;
				rotated_plm_wfn[l].resize(n);
				for (i = 0; i < n; i++)
					rotated_plm_wfn[l][i] = 0;

				for (i = 0; i < n; i++)
					for (j = 0; j < n; j++)
						rotated_plm_wfn[l][i] += rotation[l][j][i] * plm_wfn[l][j];
			}

			//

			rotated_plm = rotated_plm_wfn;

			for (l = 0; l <= l_max; l++)
			{
				n = 2 * l + 1;
				for (i = 0; i < n; i++)
					rotated_plm[l][i] *= real_spherical_harmonics::wfnNormalizationMultipliers[l][i] /
					real_spherical_harmonics::densityNormalizationMultipliers[l][i];
			}

			
			plm = rotated_plm;

		}

	}
}
