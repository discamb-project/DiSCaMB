#pragma once

#include <cstdlib>
#include <vector>
#include <cstring>

namespace discamb {
    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */

	class SphMatrix
	{
	public:
		SphMatrix(void);
		~SphMatrix(void);
		double &operator()(int i, int j);
		void setL(int l);


		int mL;
		double mZero;
		std::vector<std::vector<double> > mData;
	};
    /** @}*/
}
