#include "discamb/MathUtilities/SphConverter.h"
#include <cmath>
#include <cstdlib>

using namespace std;

namespace discamb {

	SphConverter::SphConverter(void)
	{
		r.setL(1);
		mMaxL = 0;
	}

	SphConverter::~SphConverter(void)
	{
	}

	void SphConverter::convert(
		const std::vector<std::vector<double> > &oldCoordinates,
		const std::vector<std::vector<double> > &newCoordinates,
		std::vector<std::vector<std::vector<double> > > &conversionMatrices)
	{
		int l, maxL = mMaxL;
		int m, m_prime, intL, i;
		double denominator;

		set_r(oldCoordinates, newCoordinates);

		R[0](0, 0) = 1;
		R[1].mData = r.mData;

		conversionMatrices.resize(maxL + 1);
		conversionMatrices[0] = R[0].mData;
		conversionMatrices[1] = R[1].mData;
		for (l = 2; l <= maxL; l++)
		{
			intL = static_cast<int>(l);
			// calculate u,v,w,P
			for (m = -intL; m <= intL; m++)
				for (m_prime = -intL; m_prime <= intL; m_prime++)
				{
					abs(m_prime) < intL ? denominator = double((l + m_prime)*(l - m_prime)) :
						denominator = double(2 * l*(2 * l - 1));

					u[l](m, m_prime) = sqrt(double((l + m)*(l - m)) / denominator);
					v[l](m, m_prime) = 0.5*sqrt((1.0 + delta(m, 0))*double((l + abs(m) - 1)*(l + abs(m))) / denominator)*(1 - 2 * delta(m, 0));
					w[l](m, m_prime) = -0.5*sqrt(double((l - abs(m) - 1)*(l - abs(m))) / denominator)*(1 - delta(m, 0));

					for (i = -1; i <= 1; i++)
					{
						if (abs(m_prime) < intL)
							P(i, l, m, m_prime) = r(i, 0)*R[l - 1](m, m_prime);
						if (m_prime == intL)
							P(i, l, m, m_prime) = r(i, 1)*R[l - 1](m, intL - 1) - r(i, -1)*R[l - 1](m, -intL + 1);
						if (m_prime == -intL)
							P(i, l, m, m_prime) = r(i, 1)*R[l - 1](m, -intL + 1) + r(i, -1)*R[l - 1](m, intL - 1);
					}
				}
			//calculate U V W
			for (m = -intL; m <= intL; m++)
				for (m_prime = -intL; m_prime <= intL; m_prime++)
				{
					if (m == 0)
					{
						U[l](m, m_prime) = P(0, l, 0, m_prime);
						V[l](m, m_prime) = P(1, l, 1, m_prime) + P(-1, l, -1, m_prime);
						W[l](m, m_prime) = 0;
					}
					if (m > 0)
					{
						U[l](m, m_prime) = P(0, l, m, m_prime);
						V[l](m, m_prime) = P(1, l, m - 1, m_prime)*sqrt(1.0 + delta(m, 1)) - P(-1, l, -m + 1, m_prime)*(1 - delta(m, 1));
						W[l](m, m_prime) = P(1, l, m + 1, m_prime) + P(-1, l, -m - 1, m_prime);
					}
					if (m < 0)
					{
						U[l](m, m_prime) = P(0, l, m, m_prime);
						V[l](m, m_prime) = P(1, l, m + 1, m_prime)*(1 - delta(m, -1)) + P(-1, l, -m - 1, m_prime)*sqrt(1 + delta(m, -1));
						W[l](m, m_prime) = P(1, l, m - 1, m_prime) - P(-1, l, -m + 1, m_prime);
					}
				}
			// calculate R

			for (m = -intL; m <= intL; m++)
				for (m_prime = -intL; m_prime <= intL; m_prime++)
					R[l](m, m_prime) = u[l](m, m_prime)*U[l](m, m_prime)
					+ v[l](m, m_prime)*V[l](m, m_prime)
					+ w[l](m, m_prime)*W[l](m, m_prime);


			conversionMatrices[l] = R[l].mData;
		}

	}

	void SphConverter::set_r(
		const std::vector<std::vector<double> > &oldCoordinates,
		const std::vector<std::vector<double> > &newCoordinates)
	{
		int i, j, k;
		int xyz[] = { 1,-1,0 };


		for (i = 0; i < 3; i++)
			for (j = 0; j < 3; j++)
			{
				r(xyz[i], xyz[j]) = 0;
				for (k = 0; k < 3; k++)
					r(xyz[i], xyz[j]) += oldCoordinates[i][k] * newCoordinates[j][k];
			}
	}


	void SphConverter::setMaxL(
		int maxL)
	{
		int l;
		mMaxL = maxL;

		R.resize(maxL + 1);
		U.resize(maxL + 1);
		V.resize(maxL + 1);
		W.resize(maxL + 1);
		P_data.resize(3, vector<SphMatrix>(maxL + 1));
		u.resize(maxL + 1);
		v.resize(maxL + 1);
		w.resize(maxL + 1);

		for (l = 0; l <= maxL; l++)
		{
			R[l].setL(l);
			U[l].setL(l);
			V[l].setL(l);
			W[l].setL(l);
			P_data[0][l].setL(l);
			P_data[1][l].setL(l);
			P_data[2][l].setL(l);
			u[l].setL(l);
			v[l].setL(l);
			w[l].setL(l);
		}


	}

	double SphConverter::delta(
		int i,
		int j)
	{
		if (i == j)
			return 1.0;
		return 0.0;
	}

	double &SphConverter::P(
		int i,
		int l,
		int m,
		int m_prime)
	{
		return P_data[i + 1][l](m, m_prime);
	}

}