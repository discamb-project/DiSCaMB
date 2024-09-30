#ifndef _DISCAMB_CRYSTALSTRUCTURE_SPACEGROUPOPERATION_HPP_
#define _DISCAMB_CRYSTALSTRUCTURE_SPACEGROUPOPERATION_HPP_

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/MathUtilities/Matrix3.h"
#include "CrystallographicRational.h"



#include <string>

namespace discamb {

/**
 * \addtogroup CrystalStructure
 * @{
 */


/**
Represents symmetry operation in a space group. 
*/

class SpaceGroupOperation
{
public:
    /**
    \brief Constructs identity operation (X,Y,Z).
    */
    SpaceGroupOperation();
    /**
    \brief Constructs operation given by the string. 

    The operation should be of the
    form component_1,component_2,component_3 e.g. X,Y,Z. The operation string 
    is not case sensitive e.g. x,Y,Z works. Translation can be given as integer, 
    rational number or real number e.g. x+5/4,y-x+1.25,-1+z. 
    */
    SpaceGroupOperation(const std::string &operationAsString);
    /**
    \brief Constructs operation as given by the 'rotation' matrix and translation vector.
    */
    SpaceGroupOperation(const Matrix3<int> &rotation, const Vector3<CrystallographicRational> &translation = 
                                                            Vector3<CrystallographicRational>());
    
    ~SpaceGroupOperation();
    /** \brief Sets the symmetry operation to translation.
    
    \param t the translation to which the symmetry operation is set to ({1|t})
    */
    void setToTranslation(const Vector3<CrystallographicRational> &t);
    /** \brief Sets the symmetry operation to translation.
    
    \f$ \mathbf{T} = [\frac{a_{num}}{a_{denom}} , \frac{b_{num}}{b_{denom}}, \frac{c_{num}}{c_{denom}}] \f$
    */
    void setToTranslation(int a_num,int a_denom,int b_num, int b_denom, int c_num, int c_denom);
    /** \brief Sets the symmetry operation to translation. 
    
    The components of the parameter t should not differ by more than 0.01
    from values allowed for space group translation.*/
    void setToTranslation(const Vector3d &t);

    void invert();
    /**
        \brief Sets operation to the one defined by the argument 'operationAsString'.
    
    The operation should be of the
        form component_1, component_2, component_3 e.g.X, Y, Z.The operation string
        is not case sensitive e.g.x, Y, Z works.Translation can be given as integer,
        rational number or real number e.g.x + 5 / 4, y - x + 1.25, -1 + z.*/
    void set(const std::string &operationAsString);
    /** \brief Sets operation.*/ 
    void set(const Matrix3<int> &rotation, const Vector3<CrystallographicRational> &translation = Vector3<CrystallographicRational>());
    /** \brief Gets string representation of the symmetry operation. */
    void get(std::string &operationAsString) const;
    /** \brief Gets string representation of the symmetry operation. */
    std::string string() const;
    /** \brief Gets algebraic representation of the symmetry operation.*/
    void get(Matrix3<double> &rotation,Vector3d &translation) const;

    /** \brief Gets algebraic representation of the symmetry operation.*/
    void get(Matrix3<int> &rotation, Vector3<CrystallographicRational> &translation) const;


    /** \brief Gets rotational part symmetry operation.*/
    void getRotation(Matrix3<int> &rotation) const;
    /** \brief Gets rotational part symmetry operation.*/
    void getRotation(Matrix3<double> &rotation) const;
    /** \brief Gets translational part symmetry operation.*/
    void getTranslation(Vector3<CrystallographicRational> &translation) const;
    /** \brief Gets translational part symmetry operation.*/
    void getTranslation(Vector3d &translation) const;
    /** \brief Apply the operation to vector.*/
    void apply(const Vector3d &original,Vector3d &transformed) const;
    
    /** \brief Returns true if this symmetry operation (\f$\mathbf{S}\f$) is related to the symmetry operation provided as argument (\f$\mathbf{U}\f$)
    by lattice translation:

    \f$\mathbf{S} = \{1|\mathbf{t}\}\mathbf{U}\f$, \f$\mathbf{t}=[x,y,z]\f$, \f$x,y,z\in \mathbb{Z} \f$
    */
    bool isLatticeTranslationEquivalent(const SpaceGroupOperation &symmetryOperation,Vector3i &latticeTranslation) const;

	bool isIdentity() const;

    /** \brief Parse string representing symmetry operation and converts to algebraic form. 

    Returns true if parsing was successful and false otherwise. */
    static bool parse(const std::string &symmetryOperationAsString,
                      Vector3<CrystallographicRational> &translationVector,Matrix3<int> &rotationMatrix);

    /** \brief Converts symmetry operation represented in algebraic form into string notation.*/
    static void convertSymmetryOperationToString(const Vector3<CrystallographicRational> &translationVector,
                                                 const Matrix3<int> &rotationMatrix,std::string &symmOpAsString);
    /** \brief Checks if provided string represents space group symmetry operation.*/
    static bool isSymmetryOperation(const std::string &s);
    /** \brief Checks if real number can correspond to valid component of translation vector of space group operation.
    
    Returns false if it differs by more than 0.01 from such a value. 
        */
    static bool crystallographicTranslationComponent(double d, CrystallographicRational &result);
    /** \brief Checks if 3D vector of real numbers can correspond to valid translation vector of space group operation.

    Returns false if any of the components of the vector differs by more than 0.01 from such allowed values.
    */

    static Vector3<CrystallographicRational> crystallographicTranslation(const Vector3d &translation_real, bool &successful);
    
private:

    std::string mOperationAsString;
    
    Vector3<CrystallographicRational> mTranslation;
    
    Matrix3<int> mRotation;
    bool mIsIdentity;

    
    static bool parseStringComponent(const std::string &component, int &x, int &y, int &z, CrystallographicRational &translation);
    
    static bool isCrystallographicTranslation(const std::string &s, CrystallographicRational &translation);
    
    static bool isCrystallographicDenominator(int denominator);

    static bool crystallographicFraction(double fraction,int &numerator,int &denominator);
};

SpaceGroupOperation operator*(const SpaceGroupOperation &op1,const SpaceGroupOperation &op2);

/**@}*/

} //namespace discamb




#endif /*_DISCAMB_CRYSTALSTRUCTURE_SPACEGROUPOPERATION_HPP_*/
