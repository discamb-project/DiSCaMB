
#ifndef _DISCAMB_CRYSTALSTRUCTURE_CRYSTALLOGRAPHICRATIONAL_H_
#define _DISCAMB_CRYSTALSTRUCTURE_CRYSTALLOGRAPHICRATIONAL_H_

#include <string>

namespace discamb {

    /**
    * \addtogroup CrystalStructure
    * @{
    */


/** \brief Rational numbers limited to those which can represent fractional part of translation in space group operations. 
    
    Correct denominators are 1, 2, 3, 4 and 6. 
    In general common divider for numerator and denominator is not searched for at assignment, 
    therefore on attempt to assign CrystalographicRational to e.g. 6/12 simplification to 1/2
    is not performed and exception is thrown.
    */

class CrystallographicRational{
public:
    CrystallographicRational();
    CrystallographicRational(int);
    CrystallographicRational(int numerator, int denominator);
    /** Creates CrystallographicRational, the argument should be within 0.01 of the allowed values of rational numbers
    (otherwie exception is thrown).*/
    CrystallographicRational(double value);
    /** \param s string representing integer (e.g. 5) or rational number (e.g. 12/3)*/
    CrystallographicRational(const std::string &s);
    ~CrystallographicRational();
    /** \brief Assigns to integer. */
    void assign(int numerator);
    /** \brief Assigns to \f$ \frac{numerator}{denominator}\f$.*/
    void assign(int numerator,int denominator);    
    /** \brief assigns to rational number which differs no more than 0.01 from the value.
    If there is no such a number an exception is thrown.*/
    void assign(double value);
    /** \param s string representing integer (e.g. 5) or rational number (e.g. 12/3)*/
    void assign(const std::string &s);
    void get(int &numerator, int &denominator) const;
    int denominator() const;
    int numerator() const;
    double get() const;
    operator double() const {return get();}
    std::string asString() const;
    /**\brief Represents the fraction as rational number (numerator/denominator) which differs no more than 0.01 from the fraction.
    If there is no such a number an exception is thrown.*/
    static bool isCrystallographicFraction(double fraction,int &numerator,int &denominator);
    CrystallographicRational &operator*=(int n);
    CrystallographicRational &operator+=(const CrystallographicRational &r);
    CrystallographicRational &operator-=(const CrystallographicRational &r);
private:
    int mNumerator;
    int mDenominator;
    static const double mFractionsReal[8];
    static const int mFractionsDenominator[8]; 
    static const int mFractionsNumerator[8]; 
    void gcdToOne();
};

CrystallographicRational operator*(const CrystallographicRational &,int);
CrystallographicRational operator*(int, const CrystallographicRational &);
CrystallographicRational operator+(const CrystallographicRational &, const CrystallographicRational &);
CrystallographicRational operator-(const CrystallographicRational &, const CrystallographicRational &);

bool operator<(const CrystallographicRational &v1, const CrystallographicRational &v2);
bool operator>(const CrystallographicRational &v1, const CrystallographicRational &v2);
bool operator==(const CrystallographicRational &v1, const CrystallographicRational &v2);

/**@}*/

}//namespace discamb


#endif /*_DISCAMB_CRYSTALSTRUCTURE_CRYSTALLOGRAPHICRATIONAL_H_*/


