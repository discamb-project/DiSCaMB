#include "discamb/CrystalStructure/CrystallographicRational.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include <cmath>

using namespace std;

namespace discamb {

// 0 1/2 1/3 2/3 1/4 3/4 1/6 5/6
const double CrystallographicRational::mFractionsReal[8] = { 0.0 , 0.5 , 1.0 / 3.0 , 2.0 / 3.0 , 0.25 , 0.75 , 1.0 / 6.0 , 5.0 / 6.0 };
const int CrystallographicRational::mFractionsDenominator[8] = { 1,2,3,3,4,4,6,6 };
const int CrystallographicRational::mFractionsNumerator[8] = { 0,1,1,2,1,3,1,5 };


CrystallographicRational::CrystallographicRational()
{
   mNumerator = 0;
   mDenominator = 1;
}

CrystallographicRational::CrystallographicRational(
int n)
{
    mNumerator = n;
    mDenominator = 1;
}

CrystallographicRational::CrystallographicRational(
    int numerator, 
    int denominator)
{
    assign(numerator,denominator);
}

CrystallographicRational::CrystallographicRational(
    double value)
{
    assign(value);
}

CrystallographicRational::CrystallographicRational(
    const std::string &s)
{
    assign(s);
}

CrystallographicRational::~CrystallographicRational()
{
}

void CrystallographicRational::assign(
    const std::string &s)
{
    vector<string> words;
    int numerator,denominator;
    string_utilities::split(s,words,'/');
    numerator = atoi(words[0].c_str());
    words.size()==2 ? denominator = atoi(words[1].c_str()) : denominator = 1;
    assign(numerator,denominator);
}

void CrystallographicRational::assign(int n)
{
    mNumerator = n;
}

void CrystallographicRational::assign(
    int numerator, 
    int denominator)
{
    if(denominator<0)
    {
        numerator = -numerator;
        denominator = -denominator;
    }

    bool has_crystallographic_denominator = ( denominator == 1 || denominator == 2 || 
                                              denominator == 3 || denominator == 4 ||
                                              denominator == 6);

    if( ! has_crystallographic_denominator )
    {
        string error_message = string("incorrect denominator in definition of crystallographic translation: ") +
            string_utilities::convertToString(denominator);
        on_error::throwException(error_message,__FILE__,__LINE__);       
    }

    mNumerator = numerator;
    mDenominator = denominator;
    gcdToOne();
}

void CrystallographicRational::gcdToOne()
{
    // check if denominator and numerator have common divisor
    // if yes dividem them by it

    if (mNumerator % mDenominator == 0)
    {
        mNumerator /= mDenominator;
        mDenominator = 1;
    }
    else
    {
        if (mDenominator == 4)
            if (mNumerator % 2 == 0)
            {
                mNumerator /= 2;
                mDenominator = 2;
            }
        if (mDenominator == 6)
        {
            if (mNumerator % 3 == 0)
            {
                mNumerator /= 3;
                mDenominator /= 3;
            }
            if (mNumerator % 2 == 0)
            {
                mNumerator /= 2;
                mDenominator /= 2;
            }
        }
    }

}

void CrystallographicRational::assign(
    double value)
{
    // separate into int and fraction
    double fraction;
    int asInteger;

    asInteger = int(value);
    fraction = value - asInteger; // value = asInteger + fraction;

    // find 'crystallographic' rational number for the fraction
    int numerator, denominator;
    if (!isCrystallographicFraction(fraction, numerator, denominator))
    {
        string error_message = string("an attempt to define 'non-crystallographic fraction' ") + 
            string_utilities::convertToString(fraction) + string(" as symmetry operation translation");
        on_error::throwException(error_message,__FILE__,__LINE__);
    }
    
    // combine integral part with fraction as rational number
    mNumerator = asInteger*denominator + numerator;
    mDenominator = denominator;
    
    gcdToOne();
}

std::string CrystallographicRational::asString()
const
{
    if (mDenominator == 1)
        return string_utilities::convertToString(mNumerator);
    return string_utilities::convertToString(mNumerator)+string("/")+string_utilities::convertToString(mDenominator);
}

bool CrystallographicRational::isCrystallographicFraction(
    double fraction, 
    int &numerator, 
    int &denominator)
{
    bool non_negative = fraction>=0;
    double abs_fraction = fabs(fraction);
    bool found = false;
    numerator = 0;
    denominator = 1;

    for( int i=0;i<8;++i)
        if( fabs(mFractionsReal[i]-abs_fraction)<0.01 )
        {
            found=true;
            numerator = static_cast<int>(mFractionsNumerator[i]);
            denominator = static_cast<int>(mFractionsDenominator[i]);
        }
    
    if(!non_negative)
        numerator *= -1;
    return found;
}

void CrystallographicRational::get(
    int &numerator, 
    int &denominator)
const
{
    numerator = mNumerator;
    denominator = mDenominator;
}

int CrystallographicRational::denominator()
const
{
    return mDenominator;
}

int CrystallographicRational::numerator()
const
{
    return mNumerator;
}


CrystallographicRational &CrystallographicRational::operator*=(
    int n)
{
    mNumerator*=n;
    gcdToOne();
    return *this;
}

CrystallographicRational &CrystallographicRational::operator+=(
const CrystallographicRational &r)
{
    int d1, d2, n1, n2;
    d1 = mDenominator;
    d2 = r.mDenominator;
    n1 = mNumerator;
    n2 = r.mNumerator;

    // looks for least common multiple of denominators
    int lcm;

    if ((d1 == 4 && (d2 == 3 || d2 == 6)) || (d2 == 4 && (d1 == 3 || d1 == 6)))
        lcm = 12;
    else
        (d1 == 2 && d2 == 3) || (d2 == 2 && d1 == 3) ? lcm = 6 : lcm = max(d1, d2);

    if (lcm > 6)
    {
        string error_message = string("an attempt to add crystallographic translations: ") +
            string_utilities::convertToString(n1) + string("/") + string_utilities::convertToString(d1) + string(" and ") +
            string_utilities::convertToString(n2) + string("/") + string_utilities::convertToString(d2) +
            string(" resulting in non crystallographic translation");
        on_error::throwException(error_message, __FILE__, __LINE__);
    }

    n1 *= lcm / d1;
    n2 *= lcm / d2;

    assign(n1 + n2, lcm);
    return *this;
}

CrystallographicRational &CrystallographicRational::operator-=(
    const CrystallographicRational &r)
{
    CrystallographicRational r_minus = -1*r;
    operator+=(r_minus);
    return *this;
}

double CrystallographicRational::get() 
const
{
    return double(mNumerator)/double(mDenominator);
}


CrystallographicRational operator*(const CrystallographicRational &cr, int n)
{
    CrystallographicRational result = cr;
    result *= n;
    return result;
}

CrystallographicRational operator*(int n, const CrystallographicRational &cr)
{
    CrystallographicRational result = cr;
    result *= n;
    return result;
}

CrystallographicRational operator+(const CrystallographicRational &cr1, const CrystallographicRational &cr2)
{
    CrystallographicRational result = cr1;
    result += cr2;
    return result;
}

CrystallographicRational operator-(const CrystallographicRational &cr1, const CrystallographicRational &cr2)
{
    CrystallographicRational result = cr1;
    result -= cr2;
    return result;
}

bool operator<(
    const CrystallographicRational &v1, 
    const CrystallographicRational &v2) 
{ 
    return v1.get()<v2.get(); 
}

bool operator>(
    const CrystallographicRational &v1,
    const CrystallographicRational &v2) 
{ 
    return v1.get()>v2.get(); 
}

bool operator==(
    const CrystallographicRational &v1, 
    const CrystallographicRational &v2)
{
    return ((v1.denominator()==v2.denominator()) && (v1.numerator() == v2.numerator()));
}

}//namespace discamb
