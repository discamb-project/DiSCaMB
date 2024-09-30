#ifndef _DISCAMB_MATHUTILITIES_REALSPHERICALHARMONICS_H_
#define _DISCAMB_MATHUTILITIES_REALSPHERICALHARMONICS_H_

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/MathUtilities/MathUtilities.h"

#include <cstddef>
#include <vector>
#include <string>
#include <cmath>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


namespace real_spherical_harmonics
{

    static const double densityNormalizationMultipliers[5][9] = 
    {
     { 0.0795774715459482},
     { 0.3183098861837907, 0.3183098861837907, 0.3183098861837907},
     { 0.7500000000000036, 0.7500000000000036, 0.2067483357831728, 0.7500000000000036, 0.3750000000000018},
     { 0.4244131815783849, 1.999999999999998, 0.3203330895838574, 0.2448537586029175, 0.3203330895838574, 0.999999999999999, 0.4244131815783849},
     { 1.874999999999999, 1.249999999999994, 0.6611826845788106, 0.4740025188689597, 0.06941752438437611, 0.4740025188689597, 0.3305913422894053, 1.249999999999994, 0.4687499999999998}
    };

    static const double wfnNormalizationMultipliers[8][15] =
    {
     { 0.2820947917738781},
     { 0.4886025119029199, 0.4886025119029199, 0.4886025119029199},
     { 1.092548430592079, 1.092548430592079, 0.31539156525252, 1.092548430592079, 0.5462742152960396},
     { 0.5900435899266435, 2.890611442640554, 0.4570457994644657, 0.3731763325901154, 0.4570457994644657, 1.445305721320277, 0.5900435899266435},
     { 2.503342941796705, 1.77013076977993, 0.9461746957575601, 0.6690465435572892, 0.1057855469152043, 0.6690465435572892, 0.47308734787878, 1.77013076977993, 0.6258357354491761},
     { sqrt(693.0/(512*M_PI)),sqrt(3465.0 / (256 * M_PI)) ,sqrt(385.0 / (512 * M_PI)), sqrt(1155.0 / (64.0 * M_PI)), sqrt(165.0 / (256.0 * M_PI)), sqrt(11.0 / (256.0 * M_PI)), sqrt(165.0 / (256.0 * M_PI)), sqrt(1155.0 / (64.0 * M_PI)),sqrt(385.0 / (512 * M_PI)) ,sqrt(3465.0 / (256 * M_PI)), sqrt(693.0 / (512 * M_PI))},
     { sqrt(3003.0 / (2048.0 * M_PI)),sqrt(9009.0 / (512.0 * M_PI)),sqrt(819.0 / (1024.0 * M_PI)),sqrt(1365.0 / (512.0 * M_PI)),sqrt(1365.0 / (2048.0 * M_PI)),sqrt(273.0 / (256.0 * M_PI)),sqrt(13.0 / (1024.0 * M_PI)),sqrt(273.0 / (256.0 * M_PI)),sqrt(1365.0 / (2048.0 * M_PI)),sqrt(1365.0 / (512.0 * M_PI)),sqrt(819.0 / (1024.0 * M_PI)),sqrt(9009.0 / (512.0 * M_PI)),sqrt(3003.0 / (2048.0 * M_PI))},
        {sqrt(6435.0 / (4096.0 * M_PI)),sqrt(45045.0 / (2048.0 * M_PI)),sqrt(3465.0 / (4096.0 * M_PI)),sqrt(3465.0 / (1024.0 * M_PI)),sqrt(315.0 / (4096.0 * M_PI)),sqrt(315.0 / (2048.0 * M_PI)),sqrt(105.0 / (4096.0 * M_PI)),sqrt(15.0 / (1024.0 * M_PI)),sqrt(105.0 / (4096.0 * M_PI)),sqrt(315.0 / (2048.0 * M_PI)),sqrt(315.0 / (4096.0 * M_PI)),sqrt(3465.0 / (1024.0 * M_PI)),sqrt(3465.0 / (4096.0 * M_PI)),sqrt(45045.0 / (2048.0 * M_PI)),sqrt(6435.0 / (4096.0 * M_PI))}
    };

    /** Calculaes values of density normalized real spherical harmonics, \p point has to be normalized.*/
    template<int MAX_L>
    void getDensityNormalized(const Vector3d &point,std::vector<std::vector<double> > &values);


    /** Calculaes value of density normalized real spherical harmonics, \p point has to be normalized.*/
    double densityNormalized(const Vector3d &point, int l,int m);

    /** Calculaes values of 'wave-function' normalized ('normal' ones, commonly used in math in physics) 
    real spherical harmonics, \p point has to be normalized.*/
    template<int MAX_L>
    void getWfnNormalized(const Vector3d &point,std::vector<std::vector<double> > &values);

    /** Calculaes value of 'wave-function' normalized ('normal' ones, commonly used in math in physics)
    real spherical harmonics, \p point has to be normalized.*/

    double wfnNormalized(const Vector3d &normalizedVector3D, int l, int m);
    /**
    main rotation axis parallel to z, if absent mirror plane perpendicular to z, in addition: 

    \f$\bar{4}2m - C_2 \parallel x \f$; 

    32 \f$ C_2 \parallel y \f$; 

    \f$\bar{3}m - m \bar y \f$;

    \f$\bar{4}m2 - m_2 \bar y \f$.

       
    */
    bool symmetryInvariant(const std::string &pointGroupLabel, std::vector<std::pair<int, int> > &plms);
    
    double polynomial(const Vector3d &noralizedVector3D, int l, int m);

    /** spherical harmonic y_lm can be written as:
        y_lm = r^(-l) * a_lm * w_lm(x,y,z)
        where a_lm is a const. multiplier and w_lm(x,y,z) is a polynomial with integer coefficients

        normalizedPoint - it is expected to be normalized
        values - are expected to be allocated
    */
    template<int MAX_L>
    void getPolynomials(const Vector3d &normalizedPoint, std::vector<std::vector<double> > &values);
    void getMultipliers_DensityNormalization(int maxL, std::vector<std::vector<double> > &values);
    void getMultipliers_WfnNormalization(int maxL, std::vector<std::vector<double> > &values);
    void getDensityToWfnMultipliers(int maxL,std::vector<std::vector<double> > &values);
    //void getWfnToDensityMultipliers(int maxL,std::vector<std::vector<double> > &values);


    // TEMPLATE IMPLEMENTATION


    template<int MAX_L>
    inline void getDensityNormalized(
        const Vector3d &point, 
        std::vector<std::vector<double> > &values)
    {
        const int nL = MAX_L + 1;
        //static const vector<vector<double> > wfnNormalizationFactors = getWfnNormalizationFactors();

        getPolynomials<MAX_L>(point, values);

        for (int l = 0; l<nL; l++)
            for (int m = -l; m <= l; m++)
                values[l][l + m] *= densityNormalizationMultipliers[l][l + m];
    }

    template<int MAX_L>
    inline void getWfnNormalized(
        const Vector3d &point, 
        std::vector<std::vector<double> > &values)
    {
        const int nL = MAX_L + 1;
        //static const vector<vector<double> > wfnNormalizationFactors = getWfnNormalizationFactors();

        getPolynomials<MAX_L>(point, values);

        for (int l = 0; l < nL; l++)
            for (int m = -l; m <= l; m++)
            {
                int idx = l + m;
                values[l][idx] *= wfnNormalizationMultipliers[l][idx];
            }

    }

    template<int maxL>
    inline void getPolynomials(
        const Vector3d &normalizedPoint, 
        std::vector<std::vector<double> > &values)
    {
        if (maxL<0)
            return;

        values[0][0] = 1.0;

        if (maxL > 0)
        {

            const double x = normalizedPoint[0];
            const double y = normalizedPoint[1];
            const double z = normalizedPoint[2];

            values[1][0] = y;
            values[1][1] = z;
            values[1][2] = x;


            if (maxL > 1)
            {

                double xx = x*x;
                double yy = y*y;
                double zz = z*z;
                double xy = x*y;
                double xz = x*z;
                double yz = y*z;
                double xx_p_yy = xx + yy;
                double xx_m_yy = xx - yy;
                double two_zz = 2 * zz;

                values[2][0] = xy;
                values[2][1] = yz;
                values[2][2] = two_zz - xx_p_yy;
                values[2][3] = xz;
                values[2][4] = xx_m_yy;


                if (maxL > 2)
                {
                    double _3xx_m_yy = 3 * xx - yy;
                    double xx_m_3yy = xx - 3 * yy;
                    double four_zz = 4 * zz;

                    values[3][0] = _3xx_m_yy * y;
                    values[3][1] = xy*z;
                    values[3][2] = y * (four_zz - xx_p_yy);
                    values[3][3] = z*(two_zz - 3 * xx_p_yy);
                    values[3][4] = x * (four_zz - xx_p_yy);
                    values[3][5] = (xx_m_yy)*z;
                    values[3][6] = xx_m_3yy * x;


                    if (maxL > 3)
                    {
                        double _7zz = 7 * zz;

                        values[4][0] = xy * xx_m_yy;
                        values[4][1] = _3xx_m_yy * yz;
                        values[4][2] = xy * (_7zz - 1);
                        values[4][3] = yz * (_7zz - 3);
                        values[4][4] = (_7zz - two_zz)*(_7zz - 6) + 3;
                        values[4][5] = xz * (_7zz - 3);
                        values[4][6] = (xx_m_yy)*(_7zz - 1);
                        values[4][7] = xx_m_3yy*xz;
                        values[4][8] = xx * xx_m_3yy - yy * _3xx_m_yy;

                        if (maxL > 4)
                        {
                            double x3 = xx * x;
                            double x4 = xx * xx;
                            double x5 = x4 * x;
                            double y3 = yy * y;
                            double y4 = yy * yy;
                            double y5 = y4 * y;
                            double z3 = zz * z;
                            double z4 = zz * zz;
                            double z5 = z4 * z;

                            values[5][0] = 5*x4*y-10*xx*y3+y5;
                            values[5][1] = z*(4*x3*y-4*x*y3);
                            values[5][2] = (9.0*zz-1.0)*(3.0*xx*y-y3);
                            values[5][3] = 2.0*xy*(3*z3-z);
                            values[5][4] = (21.0*z4-14.0*zz+1)*y;
                            values[5][5] = 63*z5-70*z3+15*z;
                            values[5][6] = (21.0 * z4 - 14.0 * zz + 1) * x;
                            values[5][7] = (3*z3-z)*xx_m_yy;
                            values[5][8] = (9*zz-1)*(x3-3*x*yy);
                            values[5][9] = z*(x4-6*xx*yy+y4);
                            values[5][10] = x5-10*x3*yy+5*x*y4;

                            if (maxL > 5)
                            {
                                double x6 = x * x5;
                                double y6 = y * y5;
                                double z6 = z * z5;
                                values[6][0] = 6*x5*y-20*x3*y3+6*x*y5;
                                values[6][1] = z*(5*x4*y-10*xx*y3+y5);
                                values[6][2] = (11*zz-1)*(4*x3*y-4*x*y3);
                                values[6][3] = (11*z3-3*z)*(3*xx*y-y3);
                                values[6][4] = 2*xy*(33*z4-18*zz+1);
                                values[6][5] = (33*z5-30*z3+5*z)*y;
                                values[6][6] = 231*z6-315*z4+105*zz-5;
                                values[6][7] = (33 * z5 - 30 * z3 + 5*z) * x;
                                values[6][8] = (33*z4-18*zz+1)*xx_m_yy;
                                values[6][9] = (11*z3-3*z)*(x3-3*x*yy);
                                values[6][10] = (11*zz-1)*(x4-6*xx*yy+y4);
                                values[6][11] = z*(x5-10*x3*yy+5*x*y4);
                                values[6][12] = x6-15*x4*yy+15*xx*y4-y6;

                                if (maxL > 6)
                                {
                                    double x7 = x * x6;
                                    double y7 = y * y6;
                                    double z7 = z * z6;
                                    values[7][0] = 7 * x6 * y - 35 * x4 * y3 + 21 * xx * y5 - y7;
                                    values[7][1] = z * (6 * x5 * y - 20 * x3 * y3 + 6 * x * y5);
                                    values[7][2] = (13 * zz - 1) * (5 * x4 * y - 10 * xx * y3 + y5);
                                    values[7][3] = (13 * z3 - 3 * z) * (4 * x3 * y - 4 * x * y3);
                                    values[7][4] = (143 * z4 - 66 * zz + 3) * (3 * xx * y - y3);
                                    values[7][5] = 2 * xy * (143 * z5 - 110 * z3 + 15 * z);
                                    values[7][6] = (429 * z6 - 495 * z4 + 135 * zz - 5) * y;
                                    values[7][7] = 429 * z7 - 693 * z5 + 315 * z3 - 35 * z;
                                    values[7][8] = (429 * z6 - 495 * z4 + 135 * zz - 5) * x;
                                    values[7][9] = (143 * z5 - 110 * z3 + 15 * z) * xx_m_yy;
                                    values[7][10] = (143 * z4 - 66 * zz + 3) * (x3 - 3 * yy * x);
                                    values[7][11] = (13 * z3 - 3 * z) * (x4 - 6 * xx * yy + y4);
                                    values[7][12] = (13 * zz - 1) * (5 * y4 * x - 10 * yy * x3 + x5);
                                    values[7][13] = z * (x6 - 15 * x4 * yy + 15 * y4 * xx - y6);
                                    values[7][14] = -7 * y6 * x + 35 * y4 * x3 - 21 * yy * x5 + x7;
                                }

                            }

                        }
                    } // maxL>3
                } // maxL>2
            } // maxL>1
        } // maxL>0

    } //getPolynomials

}
/** @} */
} // namespace discamb

#endif /*_DISCAMB_MATHUTILITIES_REALSPHERICALHARMONICS_H_*/
