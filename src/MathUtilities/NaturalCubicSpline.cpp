#include "discamb/MathUtilities/NaturalCubicSpline.h"
#include <vector>

using namespace std;

namespace {

	void solveTridiagonal(
		const vector<double>& a,
		const vector<double>& _b,
		const vector<double>& c,
		const vector<double>& _d,
		vector<double>& x)
	{
		std::vector<double> b = _b;
		std::vector<double> d = _d;
		int i, n = a.size();
		double w;
		x.resize(n);
		for (i = 1; i < n; i++)
		{
			w = a[i] / b[i - 1];
			b[i] -= w * c[i - 1];
			d[i] -= w * d[i - 1];
		}

		x[n - 1] = d[i - 1] / b[i - 1];
		for (int k = n - 2; k >= 0; k--)
		{
			i = k;
			x[i] = (d[i] - c[i] * x[i + 1]) / b[i];
		}
	}

}

namespace discamb {

/*


void calcSplajn(const vector<double>& y, vector<vector<double> >& coeff)
{
	coeff.clear();

	if (y.size() < 2)
		return;

	// number of points
	int n = y.size();
	// tridiagonal
	vector<double> a(n,1.0), b(n,4.0), c(n,1.0), d(n), x(n);
	d[0] = 3 * (y[1] - y[0]);

	for(int i=1; i < n-1 ; i++)
		d[i] = 3 * (y[i+1] - y[i-1]);
	d[n - 1] = 3 * (y[n - 1] - y[n - 2]);
	a[0] = 0;
	c[n - 1] = 0;

	solveTridiagonal(a, b, c, d, x);
	coeff.resize(n, vector<double>(4));

	for (int i = 0; i < n-1; i++)
	{
		coeff[i][0] = y[i];
		coeff[i][1] = x[i];
		coeff[i][2] = 3 * (y[i + 1]-y[i]) - 2 * x[i] - x[i + 1];
		coeff[i][3] = 2 * (y[i] - y[i + 1]) + x[i] + x[i + 1];
	}


}

int main(int argc, char* argv[])
{
	vector<double> y, interp;
	vector<vector<double> > coeff;
	int n = 7;
	double t, d0, d1, x;
	y.resize(n);

	d0 = 2 * M_PI / (n - 1);

	for (int i = 0; i < n; i++)
	{
		x = d0 * i;
		y[i] = sin(x);
	}


	calcSplajn(y, coeff);

	n = 70;
	d1 = 2 * M_PI / (n - 1);
	y.resize(n);
	interp.resize(n);
	int k;

	for (int i = 0; i < n; i++)
	{
		x = d1 * i;
		y[i] = sin(x);
		k = int(x / d0);
		t = (x - k * d0) / d0;
		interp[i] = coeff[k][0] + coeff[k][1] * t + coeff[k][2] * t * t + coeff[k][3] * t * t * t;
	}



	ofstream out("out");
	for (int i = 0; i < n; i++)
		out << d1 * i << " " << y[i] << " " << interp[i] << endl;
	out.close();

	cout << "OK" << endl;

}


*/

	NaturalCubicSpline::NaturalCubicSpline()
	{
		mStart = mStep = 0.0;
	}

	NaturalCubicSpline::~NaturalCubicSpline()
	{
	}

	
	NaturalCubicSpline::NaturalCubicSpline(
		const std::vector<double> &values,
		double step,
		double start)
	{
		set(values, step, start);
	}

	void NaturalCubicSpline::set(
		const std::vector<double> &values,
		double step, 
		double start)
	{
		mStart = start;
		mStep = step;


		mCoefficients.clear();

		if (values.size() < 2)
			return;

		// number of points
		int n = values.size();
		// tridiagonal
		vector<double> a(n, 1.0), b(n, 4.0), c(n, 1.0), d(n), x(n);
		d[0] = 3 * (values[1] - values[0]);

		for (int i = 1; i < n - 1; i++)
			d[i] = 3 * (values[i + 1] - values[i - 1]);
		d[n - 1] = 3 * (values[n - 1] - values[n - 2]);
		a[0] = 0;
		c[n - 1] = 0;

		solveTridiagonal(a, b, c, d, x);
		mCoefficients.resize(n, vector<double>(4));

		for (int i = 0; i < n - 1; i++)
		{
			mCoefficients[i][0] = values[i];
			mCoefficients[i][1] = x[i];
			mCoefficients[i][2] = 3 * (values[i + 1] - values[i]) - 2 * x[i] - x[i + 1];
			mCoefficients[i][3] = 2 * (values[i] - values[i + 1]) + x[i] + x[i + 1];
		}

	}
	double NaturalCubicSpline::evaluate(double x) const
	{
		int segmentIndex;
		double segmentIndexAsReal = (x - mStart) / mStep;
		if (segmentIndexAsReal < 0)
			return 0;
		segmentIndex = int(segmentIndexAsReal);
		if (segmentIndex >= mCoefficients.size())
			return 0;
		
		double t = (x - mStart - segmentIndex * mStep) / mStep;

		return mCoefficients[segmentIndex][0] + t * mCoefficients[segmentIndex][1] + t * t * mCoefficients[segmentIndex][2] + t * t * t * mCoefficients[segmentIndex][3];
	}

	double NaturalCubicSpline::operator()(double x) const
	{
		return evaluate(x);
	}


}


