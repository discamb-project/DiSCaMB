#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
#include "discamb/MathUtilities/algebra3d.h"

namespace discamb{


ReciprocalLatticeUnitCell::ReciprocalLatticeUnitCell()
{
    set(1,1,1,90,90,90);
}

ReciprocalLatticeUnitCell::ReciprocalLatticeUnitCell(const UnitCell &uc)
{
    set(uc);
}

ReciprocalLatticeUnitCell::ReciprocalLatticeUnitCell(double a,double b,double c,double alpha,double beta,double gamma)
{
    set(a,b,c,alpha,beta,gamma);
}

ReciprocalLatticeUnitCell::~ReciprocalLatticeUnitCell()
{
}


void ReciprocalLatticeUnitCell::set(const UnitCell &uc)
{

    Vector3d a,b,c,a_star,b_star,c_star;
    double v;
    uc.fractionalToCartesian(Vector3d(1,0,0),a);
    uc.fractionalToCartesian(Vector3d(0,1,0),b);
    uc.fractionalToCartesian(Vector3d(0,0,1),c);

    v = cross_product(a,b)*c;

    a_star = cross_product(b,c)/v;
    b_star = cross_product(c,a)/v;
    c_star = cross_product(a,b)/v;

    for(int i=0;i<3;i++)
    {
        mF2C(i, 0) = a_star[i];
        mF2C(i, 1) = b_star[i];
        mF2C(i, 2) = c_star[i];
    }

    algebra3d::invert3d(mF2C,mC2F);

}

void ReciprocalLatticeUnitCell::set(double a,double b,double c,double alpha,double beta,double gamma)
{
    UnitCell uc(a,b,c,alpha,beta,gamma);
    set(uc);
}



void ReciprocalLatticeUnitCell::cartesianToFractional(
const Vector3d &cart,
Vector3d &frac) 
const
{
    frac = mC2F*cart;
}

void ReciprocalLatticeUnitCell::fractionalToCartesian(
const Vector3d &frac,
Vector3d &cart) 
const
{
    cart = mF2C*frac;
}



const Matrix3d &ReciprocalLatticeUnitCell::getCartesianToFractionalMatrix() const
{
    return mC2F;
}

const Matrix3d &ReciprocalLatticeUnitCell::getFractionalToCartesianMatrix() const
{
    return mF2C;
}

}//namespace discamb

