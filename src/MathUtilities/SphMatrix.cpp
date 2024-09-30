#include "discamb/MathUtilities/SphMatrix.h"
#include <cstdlib>

using namespace std;

namespace discamb {

	SphMatrix::SphMatrix(void)
	{
		setL(0);
	}

	SphMatrix::~SphMatrix(void)
	{
	}

	double &SphMatrix::operator()(
		int i,
		int j)
	{
		mZero = 0.0;
		if (abs(i) > mL || abs(j) > mL)
			return mZero;
		return mData[i + mL][j + mL];
	}

	void SphMatrix::setL(
		int l)
	{
		mL = static_cast<int>(l);
		mData.assign(2 * l + 1, vector<double>(2 * l + 1, 0.0));
	}

}
