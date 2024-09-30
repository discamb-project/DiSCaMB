#ifndef _DISCAMB_CRYSTALSTRUCTURE_UNITCELL_HPP_
#define _DISCAMB_CRYSTALSTRUCTURE_UNITCELL_HPP_

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/MathUtilities/Matrix3.h"


namespace discamb {

/**
 * \addtogroup CrystalStructure
 * @{
 */


/** 
    Represents unit cell.
*/

class UnitCell
{
public:
    /** Constructs UnitCell with parameters (1, 1, 1, 90, 90 ,90). */
    UnitCell();
    /** Constructs UnitCell, angle parameters in degrees.*/
    UnitCell(double a,double b,double c,double alpha,double beta,double gamma);
    ~UnitCell();
    /** Sets unit cell parameters, angle parameters should be given in degrees */
    void set(double a,double b,double c,double alpha,double beta,double gamma);
    void get(double &a,double &b,double &c,double &alpha,double &beta,double &gamma) const;

    double a() const;
    double b() const;
    double c() const;
    double alpha() const;
    double beta() const;
    double gamma() const;

    void setA(double a);
    void setB(double b);
    void setC(double c);
    void setAlpha(double alpha);
    void setBeta(double beta);
    void setGamma(double gamma);

    /**
     Converts from fractional coordinates to Cartesian.
     The a-axis is collinear with the x-axis, and the c*-axis parallel to z.
    */
    void fractionalToCartesian(const Vector3d &fractional,Vector3d &cartesian) const;

    /**
    Converts from Cartesian coordinates to fractional.
    The a-axis is collinear with the x-axis, and the c*-axis parallel to z.
    */

    void cartesianToFractional(const Vector3d &cartesian,Vector3d &fractional) const;

    /** Returns fractional to Cartesian coordinates conversion matrix \f$ \mathbf{M} \f$ such that:
    \f{equation}{
    \mathbf{r}_{cart} = \mathbf{M r}_{frac}
    \f}
    */

    const Matrix3d & getFractionalToCartesianMatrix() const;

    /** Returns Cartesian to fractional coordinates conversion matrix \f$ \mathbf{M} \f$ such that:
    \f{equation}{
    \mathbf{r}_{frac} = \mathbf{M r}_{cart}
    \f}
    */
    const Matrix3d & getCartesianToFractionalMatrix() const;

private:
    Matrix3d mCartesianToFractional,mFractionaToCartesian;
    double mA,mB,mC,mAlpha,mBeta,mGamma;
    void calculateTransformationMatrices();
};

/**@}*/

}//namespace discamb




#endif /*_DISCAMB_CRYSTALSTRUCTURE_UNITCELL_HPP_*/
