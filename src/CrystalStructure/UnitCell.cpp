#include "discamb/CrystalStructure/UnitCell.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#include <cmath>

using namespace std;

namespace discamb {

UnitCell::UnitCell()
{
    set(1,1,1,90,90,90);
}

UnitCell::UnitCell(
double a,
double b,
double c,
double alpha,
double beta,
double gamma)
{
    set(a,b,c,alpha,beta,gamma);
}


UnitCell::~UnitCell()
{
}

void UnitCell::set(
double a,
double b,
double c,
double alpha,
double beta,
double gamma)
{
    mA = a;
    mB = b;
    mC = c;
    mAlpha = alpha;
    mBeta = beta;
    mGamma = gamma;

    calculateTransformationMatrices();
}

void UnitCell::get(
double &a,
double &b,
double &c,
double &alpha,
double &beta,
double &gamma)
    const
{
    a = mA;
    b = mB;
    c = mC;
    alpha = mAlpha;
    beta = mBeta;
    gamma = mGamma;
}

double UnitCell::a()
const
{
    return mA;
}

double UnitCell::b()
const
{
    return mB;
}

double UnitCell::c()
const
{
    return mC;
}

double UnitCell::alpha()
const
{
    return mAlpha;
}

double UnitCell::beta()
const
{
    return mBeta;
}

double UnitCell::gamma()
const
{
    return mGamma;
}


void UnitCell::setA(
double a)
{
    mA = a;
    calculateTransformationMatrices();
}

void UnitCell::setB(
double b)
{
    mB = b;
    calculateTransformationMatrices();
}

void UnitCell::setC(
double c)
{
    mC = c;
    calculateTransformationMatrices();
}

void UnitCell::setAlpha(
double alpha)
{
    mAlpha = alpha;
    calculateTransformationMatrices();
}

void UnitCell::setBeta(
double beta)
{
    mBeta = beta;
    calculateTransformationMatrices();
}

void UnitCell::setGamma(
double gamma)
{
    mGamma = gamma;
    calculateTransformationMatrices();
}

void UnitCell::fractionalToCartesian(
const Vector3d &fractional,
Vector3d &cartesian) 
const
{
    cartesian = mFractionaToCartesian*fractional;
    
}


void UnitCell::cartesianToFractional(
const Vector3d &cartesian,
Vector3d &fractional) 
const
{
    fractional = mCartesianToFractional*cartesian;
}

const Matrix3d & UnitCell::getFractionalToCartesianMatrix() 
const
{
    return mFractionaToCartesian;
}

const Matrix3d & UnitCell::getCartesianToFractionalMatrix() 
const
{
    return mCartesianToFractional;
}


/*
 After International Tables for Crystallography (2006). Vol. B, ch. 3.3, pp. 360-384
*/

void UnitCell::calculateTransformationMatrices()
{
    double cos_alpha,cos_beta,cos_gamma,sin_gamma,phi;
    const double f = M_PI/180;

    cos_alpha = cos(mAlpha*f);
    cos_beta = cos(mBeta*f);
    cos_gamma = cos(mGamma*f);
    sin_gamma = sin(mGamma*f);

    phi = sqrt(1-cos_alpha*cos_alpha-cos_beta*cos_beta-cos_gamma*cos_gamma+2*cos_alpha*cos_beta*cos_gamma);

    mFractionaToCartesian(0,0) = mA;
    mFractionaToCartesian(0,1) = mB * cos_gamma;
    mFractionaToCartesian(0,2) = mC * cos_beta;

    mFractionaToCartesian(1,1) = mB * sin_gamma;
    mFractionaToCartesian(1,2) =( mC / sin_gamma ) * ( cos_alpha - cos_beta * cos_gamma );

    mFractionaToCartesian(2,2) = mC * phi / sin_gamma;

    mCartesianToFractional(0,0) = 1.0 / mA;
    mCartesianToFractional(0,1) = - cos_gamma / ( mA * sin_gamma );
    mCartesianToFractional(0,2) = ( cos_alpha * cos_gamma - cos_beta ) / (mA * phi * sin_gamma);
    
    mCartesianToFractional(1,1) = 1.0 / ( mB * sin_gamma );
    mCartesianToFractional(1,2) = ( cos_beta * cos_gamma - cos_alpha ) / ( mB * phi * sin_gamma);
    
    mCartesianToFractional(2,2) = sin_gamma / ( mC * phi );
}

} //namespace discamb
