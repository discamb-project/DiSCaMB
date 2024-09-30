#include "discamb/CrystalStructure/StructuralParametersConverter.h"
#include "discamb/MathUtilities/algebra3d.h"

using namespace std;


namespace discamb {

namespace spc = structural_parameters_convention;

StructuralParametersConverter::StructuralParametersConverter()
{
    set(UnitCell(1,1,1,90,90,90));
}

StructuralParametersConverter::StructuralParametersConverter(
    const UnitCell & unitCell)
{
    set(unitCell);
}

StructuralParametersConverter::~StructuralParametersConverter()
{
}


void StructuralParametersConverter::set(
    const UnitCell &unitCell)
{   
    ReciprocalLatticeUnitCell reciprocalLattice(unitCell);
    Matrix3d m, m_inv, lambda, lambda_inv;
    Vector3d a_star, b_star, c_star;

    mCartToFrac = unitCell.getCartesianToFractionalMatrix();
    mFracToCart = unitCell.getFractionalToCartesianMatrix();


    reciprocalLattice.fractionalToCartesian(Vector3d(1, 0, 0), a_star);
    reciprocalLattice.fractionalToCartesian(Vector3d(0, 1, 0), b_star);
    reciprocalLattice.fractionalToCartesian(Vector3d(0, 0, 1), c_star);


    mRS_Edges[0] = lambda(0, 0) = sqrt(a_star*a_star);
    mRS_Edges[1] = lambda(1, 1) = sqrt(b_star*b_star);
    mRS_Edges[2] = lambda(2, 2) = sqrt(c_star*c_star);



    for(int i = 0 ; i < 3; i++)
        mRS_EdgesInv[i] = lambda_inv(i,i) = 1.0/lambda(i,i);

    mCifToCart = mFracToCart * lambda;
    mCifToCartT = algebra3d::transpose3d(mCifToCart);
    mCartToCif = lambda_inv * mCartToFrac;
    mCartToCifT = algebra3d::transpose3d(mCartToCif);
    mCartToStar = mCartToFrac;
    mCartToStarT = algebra3d::transpose3d(mCartToStar);
    mStartToCart = mFracToCart;
    mStarToCartT = algebra3d::transpose3d(mStartToCart);
    
    mCifToStar[0] = mRS_Edges[0] * mRS_Edges[0];
    mCifToStar[1] = mRS_Edges[1] * mRS_Edges[1];
    mCifToStar[2] = mRS_Edges[2] * mRS_Edges[2];
    mCifToStar[3] = mRS_Edges[0] * mRS_Edges[1];
    mCifToStar[4] = mRS_Edges[0] * mRS_Edges[2];
    mCifToStar[5] = mRS_Edges[1] * mRS_Edges[2];
    
    for(int i=0;i<6;i++)
        mStarToCif[i] = 1.0 / mCifToStar[i];

    m = reciprocalLattice.getCartesianToFractionalMatrix();
    m_inv = reciprocalLattice.getFractionalToCartesianMatrix();

    mD_XyzCartToFrac = m;
    mD_XyzFracToCart = m_inv;

    mD_UcartToUcif = lambda*m;
    mD_UcartToUcifTransposed = algebra3d::transpose3d(mD_UcartToUcif);
    mD_UcartToUstar = m;
    mD_UcartToUstarTransposed = algebra3d::transpose3d(mD_UcartToUstar);
    mD_UstarToUcart = m_inv;
    mD_UstarToUcartTransposed = algebra3d::transpose3d(mD_UstarToUcart);
    mD_UstarToUcif = lambda;
    mD_UstarToUcifTransposed = mD_UstarToUcif;
    mD_UcifToUcart = m_inv*lambda_inv;
    mD_UcifToUcartTransposed = algebra3d::transpose3d(mD_UcifToUcart);
    mD_UcifToUstar = lambda_inv;
    mD_UcifToUstarTransposed = lambda_inv;
}

void StructuralParametersConverter::xyzFractionalToCartesian(
    const Vector3d &fractional, 
    Vector3d &cartesian) 
const
{
    cartesian = mFracToCart * fractional;
}

void StructuralParametersConverter::xyzCartesianToFractional(
    const Vector3d &cartesian, 
    Vector3d &fractional) 
const
{
    fractional = mCartToFrac * cartesian;
}

void StructuralParametersConverter::uCifToU_star(
    const std::vector<double> &uCif, 
    std::vector<double> &uStar) 
const
{
    uStar = uCif;
    for(int i=0;i<6;i++)
        uStar[i] *= mCifToStar[i];
}

void StructuralParametersConverter::uStarToU_cif(
    const std::vector<double> &uStar, 
    std::vector<double> &uCif) 
const
{
    uCif = uStar;
    for(int i=0;i<6; i++)
        uCif[i] *= mStarToCif[i];
}

void StructuralParametersConverter::uCifToU_cart(
    const std::vector<double> &uCif, 
    std::vector<double> &uCart) 
const
{
    adpTransform(uCif,uCart,mCifToCart,mCifToCartT);
}

void StructuralParametersConverter::uCartToU_cif(
    const std::vector<double> &uCart, 
    std::vector<double> &uCif) 
const
{
    adpTransform(uCart,uCif,mCartToCif,mCartToCifT);
}

void StructuralParametersConverter::uCartToU_star(
    const std::vector<double> &uCart, 
    std::vector<double> &uStar) 
const
{
    adpTransform(uCart,uStar,mCartToStar,mCartToStarT);
}

void StructuralParametersConverter::uStarToU_cart(
    const std::vector<double> &uStar, 
    std::vector<double> &uCart) 
const
{
    adpTransform(uStar, uCart, mStartToCart,mStarToCartT);
}


void StructuralParametersConverter::derivativesXyzFractionalToCartesian(
    const Vector3<std::complex<double> >& derivativesWrtXyzFractional, 
    Vector3<std::complex<double> >& derivativesWrtXyzCartesian)
const
{
    derivativesWrtXyzCartesian = product<complex<double> >(mD_XyzFracToCart, derivativesWrtXyzFractional);
}

void StructuralParametersConverter::derivativesXyzCartesianToFractional(
    const Vector3<std::complex<double> >& derivativesWrtXyzCartesian,
    Vector3<std::complex<double> >& derivativesWrtXyzFractional)
const
{
    derivativesWrtXyzFractional = product<complex<double> >(mD_XyzCartToFrac, derivativesWrtXyzCartesian);
}


void StructuralParametersConverter::derivativesU_CifToU_Star(
    const std::vector<std::complex<double> >& derivativeWrtU_Cif,
    std::vector<std::complex<double> >& derivativeWrtU_Star)
const
{
    adpDerivativeTransform(derivativeWrtU_Cif, derivativeWrtU_Star,
                           mD_UcifToUstar, mD_UcifToUstarTransposed);
}

void StructuralParametersConverter::derivativesU_CifToU_Cartesian(
    const std::vector<std::complex<double> >& derivativeWrtU_Cif,
    std::vector<std::complex<double> >& derivativeWrtU_Cartesian) 
const
{
    adpDerivativeTransform(derivativeWrtU_Cif, derivativeWrtU_Cartesian,
                           mD_UcifToUcart, mD_UcifToUcartTransposed);
}

void StructuralParametersConverter::derivativesU_StarToU_Cif(
    const std::vector<std::complex<double> >& derivativeWrtU_Star,
    std::vector<std::complex<double> >& derivativeWrtU_Cif)
const
{
    adpDerivativeTransform(derivativeWrtU_Star, derivativeWrtU_Cif,
                           mD_UstarToUcif, mD_UstarToUcifTransposed);
}

void StructuralParametersConverter::derivativesU_StarToU_Cartesian(
    const std::vector<std::complex<double> >& derivativeWrtU_Star, 
    std::vector<std::complex<double> >& derivativeWrtU_Cartesian)
const
{
    adpDerivativeTransform(derivativeWrtU_Star, derivativeWrtU_Cartesian,
                           mD_UstarToUcart, mD_UstarToUcartTransposed);
}

void StructuralParametersConverter::derivativesU_CartesianToU_Star(
    const std::vector<std::complex<double> >& derivativeWrtU_Cartesian,
    std::vector<std::complex<double> >& derivativeWrtU_Star)
const
{
    adpDerivativeTransform(derivativeWrtU_Cartesian, derivativeWrtU_Star,
                           mD_UcartToUstar, mD_UcartToUstarTransposed);

}

void StructuralParametersConverter::derivativesU_CartesianToU_Cif(
    const std::vector<std::complex<double> >& derivativeWrtU_Cartesian,
    std::vector<std::complex<double> >& derivativeWrtU_Cif)
const
{
    adpDerivativeTransform(derivativeWrtU_Cartesian, derivativeWrtU_Cif,
                           mD_UcartToUcif, mD_UcartToUcifTransposed);

}

void StructuralParametersConverter::adpTransform(
    const std::vector<double> &adpIn,
    std::vector<double> &adpOut,
    const Matrix3d & trnasformMatrix,
    const Matrix3d & trnasformMatrixTransposed) 
const
{
    mAdpIn.set(adpIn[0], adpIn[3], adpIn[4],
               adpIn[3], adpIn[1], adpIn[5], 
               adpIn[4], adpIn[5], adpIn[2]);
    mAdpOut = trnasformMatrix * mAdpIn * trnasformMatrixTransposed;
    adpOut[0] = mAdpOut(0, 0);
    adpOut[1] = mAdpOut(1, 1);
    adpOut[2] = mAdpOut(2, 2);
    adpOut[3] = mAdpOut(0, 1);
    adpOut[4] = mAdpOut(0, 2);
    adpOut[5] = mAdpOut(1, 2);
}

Matrix3d StructuralParametersConverter::xyzFractionalToCartesianMatrix()
    const
{
    return mFracToCart;
}

Matrix3d StructuralParametersConverter::xyzCartesianToFractionalMatrix()
    const
{
    return mCartToFrac;
}


/*
Returns a matrix M which transforms vector of ADP components u1
into to vector expressed in other convention u2:
  u2 = M u1
*/
void StructuralParametersConverter::getAdpConversionMatrix(
    spc::AdpConvention adpConventionIn,
    spc::AdpConvention adpConventionOut,
    std::vector<std::vector<double> >& w) 
    const
{
    w.assign(6, vector<double>(6, 0));
    int i, j, k, p, q, idx;

    if (adpConventionIn == adpConventionOut)
    {
        for (i = 0; i < 6; i++)
            w[i][i] = 1.0;
        return;
    }

    if (adpConventionIn == spc::AdpConvention::U_cif && adpConventionOut == spc::AdpConvention::U_star)
    {
        for (i = 0; i < 6; i++)
            w[i][i] = mCifToStar[i];
        return;
    }

    if (adpConventionIn == spc::AdpConvention::U_star && adpConventionOut == spc::AdpConvention::U_cif)
    {
        for (i = 0; i < 6; i++)
            w[i][i] = mStarToCif[i];
        return;
    }

    Matrix3d m;

    if ((adpConventionIn == spc::AdpConvention::U_star || adpConventionIn == spc::AdpConvention::U_cif) && adpConventionOut == spc::AdpConvention::U_cart)
        m = mFracToCart;
    else
        m = mCartToFrac;

    
    vector<pair<int, int> > indices{ {0,0}, {1,1}, {2,2}, {0,1}, {0,2}, {1,2} };

    for (idx = 0; idx < 6; idx++)
    {
        i = indices[idx].first;
        j = indices[idx].second;

        for (k = 0; k < 3; k++)
            w[idx][k] = m(i, k) * m(j, k);

        for (k = 3; k < 6; k++)
        {
            p = indices[k].first;
            q = indices[k].second;
            w[idx][k] = m(i, p) * m(j, q) + m(i, q) * m(j, p);
        }
    }


    if (adpConventionIn == spc::AdpConvention::U_cif && adpConventionOut == spc::AdpConvention::U_cart)
        for(i=0;i<6;i++)
            for(j=0;j<6;j++)
                w[i][j] *= mCifToStar[j];

    if (adpConventionIn == spc::AdpConvention::U_cart && adpConventionOut == spc::AdpConvention::U_cif)
        for (i = 0; i < 6; i++)
            for (j = 0; j < 6; j++)
                w[i][j] *= mStarToCif[i];
}


void StructuralParametersConverter::adpDerivativeTransform(
    const std::vector<std::complex<double> >& adpIn,
    std::vector<std::complex<double> >& adpOut, 
    const Matrix3d & trnasformMatrix,
    const Matrix3d & trnasformMatrixTransposed)
const
{
    complex<double> d12 = 0.5 * adpIn[3];
    complex<double> d13 = 0.5 * adpIn[4];
    complex<double> d23 = 0.5 * adpIn[5];

    mDerivativeWrtAdpIn.set(adpIn[0], d12     , d13,
                            d12     , adpIn[1], d23,
                            d13     , d23     , adpIn[2]);

    mDerivativeWrtAdpOut = product<complex<double> >( 
                               trnasformMatrix ,
                               product<complex<double> >(mDerivativeWrtAdpIn,trnasformMatrixTransposed));
    adpOut.resize(6);
    for(int i=0;i<3;i++)
        adpOut[i] = mDerivativeWrtAdpOut(i, i);
    
    adpOut[3] = 2.0 * mDerivativeWrtAdpOut(0, 1);
    adpOut[4] = 2.0 * mDerivativeWrtAdpOut(0, 2);
    adpOut[5] = 2.0 * mDerivativeWrtAdpOut(1, 2);
}

void StructuralParametersConverter::convertXyz(
    const Vector3d &xyzIn, 
    Vector3d &xyzOut, 
    spc::XyzCoordinateSystem coordinateSystemIn,
    spc::XyzCoordinateSystem coordinateSystemOut)
const
{
    if( coordinateSystemIn == coordinateSystemOut)
        xyzOut = xyzIn;
    else
    {   
        if( coordinateSystemIn == spc::XyzCoordinateSystem::cartesian)
            xyzCartesianToFractional(xyzIn,xyzOut);
        else
            xyzFractionalToCartesian(xyzIn, xyzOut);
    }
}

void StructuralParametersConverter::convertADP(
    const std::vector<double> &adpIn, 
    std::vector<double> &adpOut, 
    spc::AdpConvention adpConventionIn,
    spc::AdpConvention adpConventionOut)
const
{
    if(adpConventionIn == adpConventionOut)
        adpOut = adpIn;
    else
    {
        switch(adpConventionIn)
        {
            case spc::AdpConvention::U_cart:
                if(adpConventionOut == spc::AdpConvention::U_cif)
                    uCartToU_cif(adpIn, adpOut);
                else
                    uCartToU_star(adpIn, adpOut);
                break;
            case spc::AdpConvention::U_cif:
                if(adpConventionOut == spc::AdpConvention::U_cart)
                    uCifToU_cart(adpIn, adpOut);
                else
                    uCifToU_star(adpIn, adpOut);
                break;
            case spc::AdpConvention::U_star:
                if(adpConventionOut == spc::AdpConvention::U_cart)
                    uStarToU_cart(adpIn, adpOut);
                else
                    uStarToU_cif(adpIn, adpOut);
        }
        
    }
}


void StructuralParametersConverter::convertDerivativesXyz(
    const Vector3<std::complex<double> > &dXyzIn,
    Vector3<std::complex<double> > &dXyzOut,
    structural_parameters_convention::XyzCoordinateSystem coordinateSystemIn,
    structural_parameters_convention::XyzCoordinateSystem coordinateSystemOut)
const
{
    if (coordinateSystemIn == coordinateSystemOut)
        dXyzOut = dXyzIn;
    else
    {
        if (coordinateSystemIn == spc::XyzCoordinateSystem::cartesian)
            derivativesXyzCartesianToFractional(dXyzIn,dXyzOut);
        else
            derivativesXyzFractionalToCartesian(dXyzIn, dXyzOut);
    }
}

void StructuralParametersConverter::convertDerivativesADP(
    const std::vector<std::complex<double> > &dAdpIn, 
    std::vector<std::complex<double> > &dAdpOut,
    structural_parameters_convention::AdpConvention adpConventionIn,
    structural_parameters_convention::AdpConvention adpConventionOut) 
const
{
    if (adpConventionIn == adpConventionOut)
        dAdpOut = dAdpIn;
    else
    {
        switch (adpConventionIn)
        {
        case spc::AdpConvention::U_cart:
            if (adpConventionOut == spc::AdpConvention::U_cif)
                derivativesU_CartesianToU_Cif(dAdpIn, dAdpOut);
            else
                derivativesU_CartesianToU_Star(dAdpIn, dAdpOut);
            break;
        case spc::AdpConvention::U_cif:
            if (adpConventionOut == spc::AdpConvention::U_cart)
                derivativesU_CifToU_Cartesian(dAdpIn, dAdpOut);
            else
                derivativesU_CifToU_Star(dAdpIn, dAdpOut);
            break;
        case spc::AdpConvention::U_star:
            if (adpConventionOut == spc::AdpConvention::U_cart)
                derivativesU_StarToU_Cartesian(dAdpIn, dAdpOut);
            else
                derivativesU_StarToU_Cif(dAdpIn, dAdpOut);
        }

    }
}

} //namespace discamb
