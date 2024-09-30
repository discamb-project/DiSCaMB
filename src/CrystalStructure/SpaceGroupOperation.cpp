#include "discamb/CrystalStructure/SpaceGroupOperation.h"

#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/MathUtilities/math_utilities.h"

#include <cmath>

using namespace std;


namespace {

const int anonymous_namespace_denominators[7] = {2,3,3,4,4,6,6};
const int anonymous_namespace_numerators[7] = {1,1,2,1,3,1,5};
const double anonymous_namespace_fractionsReal[7] = { 0.5 , 1.0/3.0 , 2.0/3.0 , 0.25 , 0.75 , 1.0/6.0 , 5.0/6.0 };
}

namespace discamb {

bool SpaceGroupOperation::isCrystallographicDenominator(
int denominator)
{
    int abs_d = abs(denominator);
    if(abs_d<7)
        if( abs_d != 0 && abs_d != 5 )
            return true;
    return false;

}

bool SpaceGroupOperation::crystallographicFraction(
double fraction,
int &numerator,
int &denominator)
{
    int i,n=7;
    bool found=false;
    
    if(fraction<0.01) 
    {
        numerator = 0;
        denominator = 1;
        found = true;
    }


    for(i=0;i<n;i++)
    {
        if(fabs(fraction - anonymous_namespace_fractionsReal[i])<0.01)
        {
            denominator = anonymous_namespace_denominators[i];
            numerator = anonymous_namespace_numerators[i];
            found = true;
        }
    }
    
    return found;

    //if(!found)
    //    on_error::throwException(string("attempt to define non-crystallographic translation in symmetry operation (translation by : '")+
    //                             string_utilities::convertToString(fraction)+string("' )"),__FILE__,__LINE__);
}

bool SpaceGroupOperation::crystallographicTranslationComponent(
double d,
CrystallographicRational &r)
{
    double fraction;
    int asInteger;
    int numerator,denominator;
    
    asInteger = int(d);
    fraction = d - asInteger; // d = asInteger + fraction;
    if(!crystallographicFraction(fabs(fraction),numerator,denominator))
        return false;
    r = abs(asInteger);
    r += CrystallographicRational(static_cast<int>(numerator),static_cast<int>(denominator));//boost::rational<int>(numerator,denominator);
    if(d<0)
        r *= -1;
    return true;
}

Vector3<CrystallographicRational> SpaceGroupOperation::crystallographicTranslation(
const Vector3d &translation_real,
bool &successful)
{
    Vector3<CrystallographicRational> result;

    successful = crystallographicTranslationComponent(translation_real[0], result[0]) &&
                 crystallographicTranslationComponent(translation_real[1], result[1]) &&
                 crystallographicTranslationComponent(translation_real[2], result[2]);
    return result;

}




SpaceGroupOperation::SpaceGroupOperation()
{
    set(std::string("X,Y,Z"));
}

SpaceGroupOperation::SpaceGroupOperation(
const std::string &operationAsString)
{
    set(operationAsString);

}


SpaceGroupOperation::SpaceGroupOperation(
const Matrix3<int> &rotation,
const Vector3<CrystallographicRational> &translation)
{
    set(rotation,translation);
}


SpaceGroupOperation::~SpaceGroupOperation()
{
}

void SpaceGroupOperation::setToTranslation(
const Vector3d &t)
{
    mRotation.setToIdentity();
    for(int i=0;i<3;i++)
         if(!crystallographicTranslationComponent(t[i],mTranslation[i]))
             on_error::throwException(std::string("attempt to define non-crystallographic translation in symmetry operation (translation by : '") +
                                      string_utilities::convertToString(t[i])+std::string("' )"),__FILE__,__LINE__);

}

void SpaceGroupOperation::setToTranslation(
const Vector3<CrystallographicRational> &t)
{
    Matrix3<int> rotation(1, 0, 0,
                          0, 1, 0,
                          0, 0, 1);
    set(rotation, t);
}

void SpaceGroupOperation::setToTranslation(
int a_num,
int a_denom,
int b_num, 
int b_denom, 
int c_num, 
int c_denom)
{
    setToTranslation(Vector3<CrystallographicRational>(CrystallographicRational(a_num,a_denom),
                                                       CrystallographicRational(b_num,b_denom), 
                                                       CrystallographicRational(c_num,c_denom)));
}

void SpaceGroupOperation::set(
const std::string &operationAsString)
{
    mOperationAsString = operationAsString;
    parse(mOperationAsString,mTranslation,mRotation);
	convertSymmetryOperationToString(mTranslation, mRotation, mOperationAsString);
}

bool SpaceGroupOperation::isIdentity()
const
{
	Matrix3i m(1, 0, 0,
		       0, 1, 0,
		       0, 0, 1);
	Vector3<CrystallographicRational> t(0, 0, 0);
	return (mRotation == m && mTranslation == t);
}

void SpaceGroupOperation::get(
std::string &operationAsString) 
const
{
    operationAsString = mOperationAsString;
}

std::string SpaceGroupOperation::string()
const
{
    return mOperationAsString;
}

void SpaceGroupOperation::apply(
const Vector3d &original,
Vector3d &transformed) 
const
{
    transformed = product<double>(mRotation,original);
    for(int i=0;i<3;i++)
        transformed[i] += mTranslation[i].get();
}

//bool SpaceGroupOperation::isTranslation()
//const
//{
//}

// returns true if this = symmetryOperation + latticeTranslation 
bool SpaceGroupOperation::isLatticeTranslationEquivalent(
const SpaceGroupOperation &symmetryOperation,
Vector3i &latticeTranslation) 
const
{
    int i,j;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            if( mRotation(i,j)!=symmetryOperation.mRotation(i,j) )
                return false;
    
    Vector3<CrystallographicRational> diff = mTranslation - symmetryOperation.mTranslation;

    for(i=0;i<3;i++)
    {
        if( diff[i].denominator() != 1 )
            return false;
        latticeTranslation(i) = diff(i).numerator();
    }
    return true;
}
    

bool SpaceGroupOperation::parse(
const std::string &_symmetryOperationAsString,
Vector3<CrystallographicRational> &translationVector,
Matrix3<int> &rotationMatrix)
{
    vector<std::string> componentsAsString,words;

    // remove white spaces

    std::string symmetryOperationAsString;
    string_utilities::split(_symmetryOperationAsString,words,string_utilities::CharacterType::WHITE_SPACE);
    for(int i=0;i<words.size();i++)
        symmetryOperationAsString += words[i];

    //

    string_utilities::split(symmetryOperationAsString,componentsAsString,',');
    if(componentsAsString.size()!=3)
        return false;

    for(int i=0;i<3;i++)
        if(!parseStringComponent(componentsAsString[i],rotationMatrix(i,0),rotationMatrix(i,1),rotationMatrix(i,2),translationVector(i)))
            return false;
    
    return true;
}

void SpaceGroupOperation::convertSymmetryOperationToString(
const Vector3<CrystallographicRational> &translationVector,
const Matrix3<int> &rotationMatrix,
std::string &symmOpAsString)
{
    int i,componentIndex;
    bool firstTerm;
    char xyz[] = {'X','Y','Z'};
    // check if translationVector is crystallographic
    for(i=0;i<3;i++)
        if(!isCrystallographicDenominator(translationVector[i].denominator()))
            on_error::throwException(std::string("attempt to define non-crystallographic translation in symmetry operation (translation by : '") +
                                 translationVector[i].asString()+std::string("' )"),__FILE__,__LINE__);
    
    symmOpAsString.clear();

    for( componentIndex=0 ; componentIndex<3 ; componentIndex++ )
    {
        firstTerm=true;

        for(i=0;i<3;i++)
            if( rotationMatrix(componentIndex,i) !=0 )
            {
                if( rotationMatrix(componentIndex,i) == 1 )
                    if(!firstTerm)
                        symmOpAsString += '+';
                firstTerm = false;
                if( rotationMatrix(componentIndex,i) == -1 )
                    symmOpAsString += '-';

                symmOpAsString += xyz[i];

            }

        if(translationVector[componentIndex].numerator()!=0)
            {
                if(translationVector[componentIndex].get() > 0)
                    symmOpAsString += std::string("+");
                symmOpAsString += translationVector[componentIndex].asString();
            }

        
        if(componentIndex != 2)
            symmOpAsString += ',';
        
    }

}

void SpaceGroupOperation::invert()
{
    Matrix3i identity(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1);
    
    Matrix3i r = mRotation;

    int n = 1;
    while (!(r == identity) && n < 7)
    {
        r *= mRotation;
        n++;
    }

    if (!(r == identity))
        on_error::throwException("crystallographic point operation inversion not found", __FILE__, __LINE__);

    Matrix3i rInv = mRotation;
    for (int i = 1; i < n - 1; i++)
        rInv *= mRotation;

    //Vector3<CrystallographicRational> tInv = -Matrix3<CrystallographicRational>(rInv) * mTranslation;
    Vector3<CrystallographicRational> tInv(0, 0, 0);
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            tInv[i] -= rInv(i,j) * mTranslation(j);

    set(rInv, tInv);
}

bool SpaceGroupOperation::parseStringComponent(
const std::string &component,
int &x,
int &y,
int &z,
CrystallographicRational &translation)
{
    //int multiplier=1;
    int i,nTerms,nCharacters=component.size();
    // already processed character is x-0,y-1,z-2,somthing else - -1
    //int xyz=-1;
    std::string number;
    double translationAsReal;
    //bool numberIsRational=false;
    vector<std::string> rationalNumberComponents,termsStr,numberComponents;
    vector<int> termSign;
    vector<int> nTermsOfGivenType(3,0);
    vector<char> delimiters(2);

    // set initial values

    x=0;
    y=0;
    z=0;
    translation = 0;

    //

    if(nCharacters==0)
        return false;

    delimiters[0] = '-';
    delimiters[1] = '+';

    string_utilities::split(component,termsStr,delimiters);

    for(i=0;i<nCharacters;i++)
        if( component[i] == '+' )
            termSign.push_back(1);
        else
        {
            if(component[i] == '-')
                termSign.push_back(-1);
            else
            {
                if (i == 0)
                    termSign.push_back(1);
            }
        }
    nTerms = termSign.size();

    if( nTerms != termsStr.size() || nTerms == 0 )
        return false;

    for(i=0;i<nTerms;i++)
    {
        if(termsStr[i] == std::string("X") || termsStr[i] == std::string("x"))
        {
            x = termSign[i];
            ++nTermsOfGivenType[0];
            continue;
        }
        if(termsStr[i] == std::string("Y") || termsStr[i] == std::string("y"))
        {
            y = termSign[i];
            ++nTermsOfGivenType[1];
            continue;
        }
        if(termsStr[i] == std::string("Z") || termsStr[i] == std::string("z"))
        {
            z = termSign[i];
            ++nTermsOfGivenType[2];
            continue;
        }

        // has dot e.g. 0.5
        if( termsStr[i].find('.') != string::npos )
        {
            translationAsReal = string_utilities::convertFromString<double>(termsStr[i]);
            if(!crystallographicTranslationComponent(termSign[i]*translationAsReal,translation))
                return false;
        }
        else 
        {
            translation.assign(termsStr[i]);
            translation *= termSign[i];
        }
    }

    for(i=0;i<3;i++)
        if(nTermsOfGivenType[i]>1)
            return false;

    return true;
}

bool SpaceGroupOperation::isSymmetryOperation(
const std::string &s)
{
    Vector3<CrystallographicRational> translationVector;
    Matrix3<int> rotationMatrix;

    return parse(s,translationVector,rotationMatrix);
}

void SpaceGroupOperation::get(
Matrix3<double> &rotation,
Vector3d &translation)
const
{
    for(int i=0;i<3;i++)
    {
        translation[i] = mTranslation[i].get();
        for(int j=0;j<3;j++)
            rotation(i,j) = mRotation(i,j);
    }
    
}

void SpaceGroupOperation::getRotation(
    Matrix3<int> &rotation) 
const
{
    rotation = mRotation;
}

void SpaceGroupOperation::getRotation(
    Matrix3<double> &rotation) 
const
{
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            rotation(i,j) = mRotation(i, j);
}

void SpaceGroupOperation::getTranslation(
    Vector3<CrystallographicRational> &translation) 
const
{
    translation = mTranslation;
}

void SpaceGroupOperation::getTranslation(
    Vector3d &translation) 
const
{
    translation = mTranslation;
}


SpaceGroupOperation operator*(
const SpaceGroupOperation &op1,
const SpaceGroupOperation &op2)
 {
    Matrix3i rotation_1,rotation_2;
    Vector3<CrystallographicRational> translation_1,translation_2;
    SpaceGroupOperation result;

    op1.get(rotation_1,translation_1);
    op2.get(rotation_2,translation_2);

    result.set(rotation_1*rotation_2,product<CrystallographicRational>(rotation_1,translation_2)+translation_1);
    return result;
}



void SpaceGroupOperation::set(
const Matrix3<int> &rotation,
const Vector3<CrystallographicRational> &translation)
{
    mRotation = rotation;
    mTranslation = translation;
    convertSymmetryOperationToString(translation,rotation,mOperationAsString);
}


void SpaceGroupOperation::get(
Matrix3<int> &rotation,
Vector3<CrystallographicRational> &translation)
const
{
    rotation = mRotation;
    translation = mTranslation;
}

bool SpaceGroupOperation::isCrystallographicTranslation(
const std::string &s, CrystallographicRational &translation)
{
    int i,nChar = s.size();
    int n_dividers,n_dots;

    n_dividers = 0;
    n_dots = 0;

    for(i=0;i<nChar;++i)
        if(!isdigit(s[i]))
        {
            if(s[i] == '.')
                n_dots++;
            else
            {
                if(s[i] == '/')
                    ++n_dividers;
                else
                    return false;
            }
        }

    if(n_dividers+n_dots>1)
        return false;
    if(n_dots == 1)
    {
        if(!crystallographicTranslationComponent(atof(s.c_str()),translation))
            return false;
    }
    else
    {
        translation.assign(s);
        if(!isCrystallographicDenominator(translation.denominator()))
            return false;
    }

    return true;
}

}//namespace discamb
