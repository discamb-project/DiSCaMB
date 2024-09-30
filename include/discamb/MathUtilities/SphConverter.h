#pragma once

#include <vector>
#include "SphMatrix.h"


namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


	/**
	to convert some property \f$ p_l=\sum_m P_{lm} Y_{lm}=\mathbf{P}_l \mathbf{Y}_l\f$ then
	\f$ p_l \f$ in new coordinates \f$ p_l=\mathbf{M}^T \mathbf{P}_l\f$
	*/


	class SphConverter
	{
	public:
		SphConverter();
		~SphConverter();
		void setMaxL(int maxL);
		void convert(const std::vector<std::vector<double> > &oldCoordinates,
			const std::vector<std::vector<double> > &newCoordinates,
			std::vector<std::vector<std::vector<double> > > &conversionMatrices);
	private:
		int mMaxL;
		std::vector<SphMatrix> R, U, V, W, u, v, w;
		SphMatrix r;
		double &P(int i, int l, int m, int m_prime);
		std::vector<std::vector<SphMatrix> > P_data;
		void set_r(const std::vector<std::vector<double> > &oldCoordinates,
			const std::vector<std::vector<double> > &newCoordinates);

		double delta(int i, int j);
	};
    /** @}*/
}
