#include <vector>
#include <string>
#include <cmath>

#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/BasicUtilities/on_error.h"

namespace discamb{

NGaussianFormFactor::NGaussianFormFactor()
{
    //mA.assign(5, 0.0);
    //mB.assign(5, 0.0);
    //mC = 0.0;
    //mLabel = "";
    mN = 0;
}

NGaussianFormFactor::NGaussianFormFactor(
    const std::string &label, 
    std::vector<double> &a, 
    const std::vector<double> &b, 
    double c)
{
    set_parameters(a,b,c);
    mLabel = label;
}

NGaussianFormFactor::NGaussianFormFactor(
    const std::string &label, 
    const double *a,
    const double *b,
    double c,
    int n_gaussians)
{
    set_parameters(label.c_str(),a,b,c,n_gaussians);
}


void NGaussianFormFactor::set_parameters(
    const char *label, 
    const double *a,
    const double *b,
    double c,
    int n_gaussians)
{
    mN = n_gaussians;
    mA.resize(n_gaussians);
    mB.resize(n_gaussians);
    for(int i=0;i<n_gaussians;i++)
    {
        mA[i] = a[i];
        mB[i] = b[i];
    }
    mC = c;
    mLabel = label;
}


void NGaussianFormFactor::set_parameters(
    const std::vector<double> &a,
    const std::vector<double> &b,
    double c)
{
    //if (a.size() != 5 || a.size() != 5)
      //  on_error::throwException("wrong size of form factor parameters list", __FILE__, __LINE__);
    mN = mA.size();
    mA = a;
    mB = b;
    mC = c;
}

void NGaussianFormFactor::set_label(
    const std::string &label)
{
    mLabel = label;
}

void NGaussianFormFactor::get_parameters(
    std::vector<double> &a,
    std::vector<double> &b,
    double &c) 
const
{
    a=mA;
    b=mB;
    c=mC;
}

std::string NGaussianFormFactor::get_label() const
{
    return mLabel;
}

double NGaussianFormFactor::calculate_sinth(
    double x) // x=sin_theta_over_lambda
const
{
    double result=0;
    int  n = mA.size();
    for(int i=0; i < n; i++)
        result += mA[i] * exp( -mB[i] * x * x );

    return result + mC;
}

double NGaussianFormFactor::calculate_h(
double h_length) 
const
{
    return calculate_sinth( h_length / 2.0 );
}

} // namespace discamb
