#include "discamb/MathUtilities/algebra3d.h"

#include <cmath>

using namespace std;

namespace discamb {
	namespace algebra3d
	{
		void powerRealSymm(const Matrix3d& m, double p, Matrix3d& result)
		{
			Vector3d v1, v2, v3;
			double e1, e2, e3;
			algebra3d::eigensystemRealSymm(m, v1, v2, v3, e1, e2, e3);

			Matrix3d c{ v1[0], v2[0], v3[0],
						v1[1], v2[1], v3[1],
						v1[2], v2[2], v3[2] };

			Matrix3d l{ pow(e1, p),          0.0,         0.0,
							   0.0,   pow(e2, p),         0.0,
							   0.0,          0.0,   pow(e3, p) };

			Matrix3d c_t = c;

			c_t.transpose();
			result = c * l * c_t;
		}
	}
}