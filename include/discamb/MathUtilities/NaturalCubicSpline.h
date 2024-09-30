#ifndef _DISCAMB_MATHUTILITIES_NATURALCUBICSPLINE_H_
#define _DISCAMB_MATHUTILITIES_NATURALCUBICSPLINE_H_

#include <vector>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */



class NaturalCubicSpline{
public:
	NaturalCubicSpline();
	~NaturalCubicSpline();
	NaturalCubicSpline(const std::vector<double> &values,double step, double start);
	void set(const std::vector<double> &values,double step, double start);
	double evaluate(double x) const;
	double operator()(double x) const;

private:
	/**
	the interpolating function f_i(t) for i-th segment is given by:
	f_i(t) = mCoefficients[i][0] + t * mCoefficients[i][1] + t * t * mCoefficients[i][2] + t * t * t * mCoefficients[i][3];
	assuming that the segment covers t in [0,1]
	*/
	std::vector<std::vector<double> > mCoefficients;
	double mStart, mStep;
};


/** @} */
}


#endif /*_DISCAMB_MATHUTILITIES_MATHUTILITIES_H_*/ 