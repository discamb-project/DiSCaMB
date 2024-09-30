
#ifndef _DISCAMB_CRYSTALSTRUCTURE_RECIPROCALLATTICE_H_
#define _DISCAMB_CRYSTALSTRUCTURE_RECIPROCALLATTICE_H_

#include "UnitCell.h"

namespace discamb {

/**
* \addtogroup CrystalStructure
* @{
*/


/** \brief Provides tools for conversion reciprocal space vectors between Cartesian and fractional coordinates.

The a-axis is collinear with the x-axis, and the c*-axis parallel to z.
*/

class ReciprocalLatticeUnitCell
{
public:
    ReciprocalLatticeUnitCell();
    ReciprocalLatticeUnitCell(const UnitCell &uc);
    /** Arguments are direct space unit cell parameters, angles in degrees.*/
    ReciprocalLatticeUnitCell(double a,double b,double c,double alpha,double beta,double gamma);
    ~ReciprocalLatticeUnitCell();

    void set(const UnitCell &uc);
    /** Arguments are direct space unit cell parameters, angles in degrees.*/
    void set(double a,double b,double c,double alpha,double beta,double gamma);

    /** Converts vector in reciprocal space from Cartesian to fractional coordinates. 
    The a-axis is collinear with the x-axis, and the c*-axis parallel to z.*/
    void cartesianToFractional(const Vector3d &cart,Vector3d &frac) const;
    
    /** Converts vector in reciprocal space from fractional to Cartesian coordinates. 
    The a-axis is collinear with the x-axis, and the c*-axis parallel to z.*/
    void fractionalToCartesian(const Vector3d &frac,Vector3d &cart) const;

    /** Returns Cartesian to fractional coordinates conversion matrix \f$ \mathbf{M} \f$ such that:
    \f{equation}{
    \mathbf{r}_{frac} = \mathbf{M r}_{cart}
    \f}
    */
    const Matrix3d &getCartesianToFractionalMatrix() const;

    /** Returns fractional to Cartesian coordinates conversion matrix \f$ \mathbf{M} \f$ such that:
    \f{equation}{
    \mathbf{r}_{cart} = \mathbf{M r}_{frac}
    \f}
    */
    const Matrix3d &getFractionalToCartesianMatrix() const;

private:
    Matrix3d mF2C,mC2F;
};

/**@}*/

}//namespace discamb

#endif /*_DISCAMB_CRYSTALSTRUCTURE_RECIPROCALLATTICE_H_*/


