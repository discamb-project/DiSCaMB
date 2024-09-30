#ifndef _DISCAMB_SCATTERING_NGAUSSIANFORMFACTOR_H_
#define _DISCAMB_SCATTERING_NGAUSSIANFORMFACTOR_H_

#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


class NGaussianFormFactor
{
public:
    NGaussianFormFactor();
    NGaussianFormFactor(const std::string &label, std::vector<double> &a, const std::vector<double> &b, double c);
    NGaussianFormFactor(const std::string &label, const double *a, const double *b, double c, int n_gaussians);
    //NGaussianFormFactor(const char *label,const double *a,const double *b, double c,int n_gaussians);

    void set_parameters(const std::vector<double> &a, const std::vector<double> &b, double c);
    void set_parameters(const char *label,const double *a,const double *b, double c,int n_gaussians);
    void set_label(const std::string &label);

    void get_parameters(std::vector<double> &a, std::vector<double> &b, double &c) const;
    std::string get_label() const;

    double calculate_sinth(double sin_theta_over_lambda) const;
    double calculate_h(double h_length) const;

private:
    std::vector<double> mA,mB;
    double mC;
    std::string mLabel;
    int mN;
};
/** @}*/
} // namespace discamb

#endif /*_DISCAMB_SCATTERING_NGAUSSIANFORMFACTOR_H_*/
