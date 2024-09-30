 #ifndef _DISCAMB_MATHUTILITIES_MATRIX3_H_
#define _DISCAMB_MATHUTILITIES_MATRIX3_H_

#include "Vector3.h"

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


/** \brief Represents a matrix in 3 dimensional space.*/

template<typename T>
class Matrix3{
public:

    // constructors

    /** \brief Constructs zero matrix.*/

    Matrix3();

    /** \brief Constructs zero with specified components.
    
    \f[
    \begin{bmatrix}
    a11 & a12 & a13 \\
    a21 & a22 & a23 \\
    a31 & a32 & a33 
    \end{bmatrix}
    \f]
    */

    Matrix3(const T &a11,const T &a12,const T &a13,const T &a21,const T &a22,const T &a23,const T &a31,const T &a32,const T &a33);
    Matrix3(const Matrix3 &m);


    /* Indexable needs to have element access via operator()(int, int) */
    template<typename INDEXABLE>
    Matrix3(const INDEXABLE &m);


    // assignment operator

    Matrix3<T> & operator=(const Matrix3 &);

    /** descructor */

    ~Matrix3();

    /** sets the matrix to:
     / a11 a12 a13 \
    |  a21 a22 a23  |
     \ a31 a32 a33 /
    */

    void set(const T &a11,const T &a12,const T &a13,const T &a21,const T &a22,const T &a23,const T &a31,const T &a32,const T &a33);

    // accessors
    
    T &operator()(int row,int column);
    T operator()(int row,int column) const;

    template<typename INDEX_TYPE>
    T& operator()(INDEX_TYPE row, INDEX_TYPE column);

    template<typename INDEX_TYPE>
    T operator()(INDEX_TYPE row, INDEX_TYPE column) const;


    /** unary minus */

    Matrix3<T> operator-() const;

    // scalar operations

    Matrix3<T> & operator*=(const T &scalar); 
    Matrix3<T> & operator/=(const T &scalar);

    // matrix operations

    Matrix3<T> & operator+=(const Matrix3<T> &);
    Matrix3<T> & operator-=(const Matrix3<T> &);
    Matrix3<T> & operator*=(const Matrix3<T> &); 

    // utilities

    void setToIdentity();

    bool isZero() const;

    void transpose();

private:
    T a[9]; // row major order
};	

typedef Matrix3<double> Matrix3d; 
typedef Matrix3<float> Matrix3f;
typedef Matrix3<int> Matrix3i;


///////////////////////////
//
//  NON-MEMBER OPERATORS
//
///////////////////////////


// matrix - matrix

template<typename T>
Matrix3<T> operator*(const Matrix3<T> &m1,const Matrix3<T> &m2);

template<typename T>  	  
const Matrix3<T> operator+(const Matrix3<T> &matrix1,const Matrix3<T> &matrix2);

template<typename T>
const Matrix3<T> operator-(const Matrix3<T> &matrix1,const Matrix3<T> &matrix2);

template<typename T>
bool operator==(const Matrix3<T> &matrix1, const Matrix3<T> &matrix2);

// matrix - vector

template<typename T>
const Vector3<T> operator*(const Matrix3<T> &matrix,const Vector3<T> &vector);

template<typename T>
const Vector3<T> operator*(const Vector3<T> &vector,const Matrix3<T> &matrix);

// matrix - scalar

template<typename T>
const Matrix3<T> operator*(const T &scalar,const Matrix3<T> &matrix);

template<typename T>
const Matrix3<T> operator*(const Matrix3<T> &matrix,const T &scalar);

template<typename T>
const Matrix3<T> operator/(const Matrix3<T> &matrixIn,const T &scalar);

// mixed type products

template<typename T,typename U,typename W>
const Matrix3<T> product(const Matrix3<U> &m1,const Matrix3<W> &m2);

template<typename T,typename U,typename W>
const Vector3<T> product(const Matrix3<U> &matrix,const Vector3<W> &vector);

template<typename T,typename U,typename W>
const Vector3<T> product(const Vector3<U> &vector,const Matrix3<W> &matrix);

template<typename T,typename U,typename W>
const Matrix3<T> product(const U &scalar,const Matrix3<W> &matrix);

template<typename T,typename U,typename W>
const Matrix3<T> product(const Matrix3<U> &matrix,const W &scalar);



///////////////////////////////////////////////////////////////////
//
//                        IMPLEMENTATION
//
///////////////////////////////////////////////////////////////////


// constructors


template<typename T>
Matrix3<T>::Matrix3()
{
    for(int i=0;i<9;i++)
        a[i] = 0;
}

template<typename T>
inline Matrix3<T>::Matrix3(
const T &a11,
const T &a12,
const T &a13,
const T &a21,
const T &a22,
const T &a23,
const T &a31,
const T &a32,
const T &a33)
{
    set(a11,a12,a13,a21,a22,a23,a31,a32,a33);
}

template<typename T>
inline Matrix3<T>::Matrix3(
const Matrix3<T> &m)
{
    for(int i=0;i<9;++i)
        a[i] = m.a[i];
}

template<typename T>
template<typename INDEXABLE>
Matrix3<T>::Matrix3(const INDEXABLE &m)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            operator()(i, j) = m(i, j);
}



// assignment operator


template<typename T>
inline Matrix3<T> & Matrix3<T>::operator=(
const Matrix3<T> &m) 
{
    for(int i=0;i<9;++i)
        a[i] = m.a[i];
    return *this;
}


// destructor


template<typename T>
inline Matrix3<T>::~Matrix3()
{
}


// setter


template<typename T>
inline void Matrix3<T>::set(
const T &a11,
const T &a12,
const T &a13,
const T &a21,
const T &a22,
const T &a23,
const T &a31,
const T &a32,
const T &a33)
{
    a[0] = a11;
    a[1] = a12;
    a[2] = a13;
    a[3] = a21;
    a[4] = a22;
    a[5] = a23;
    a[6] = a31;
    a[7] = a32;
    a[8] = a33;
}


// accessors


template<typename T>
inline T &Matrix3<T>::operator()(
int row,
int column) 
{
    return a[3*row+column];
}

template<typename T>
T Matrix3<T>::operator()(
int row,
int column)
const
{
    return a[3*row+column];
}

template<typename T>
template<typename INDEX_TYPE>
T& Matrix3<T>::operator()(
    INDEX_TYPE row,
    INDEX_TYPE column)
{
    return Matrix3<T>::operator()(static_cast<int>(row), static_cast<int>(column));
}

template<typename T>
template<typename INDEX_TYPE>
T Matrix3<T>::operator()(
    INDEX_TYPE row,
    INDEX_TYPE column)
    const
{
    return Matrix3<T>::operator()(static_cast<int>(row), static_cast<int>(column));
}


// unary minus


template<typename T>
inline Matrix3<T> Matrix3<T>::operator-()
const
{
    Matrix3<T> m(*this);
    m*=-1;
    return m;
}


// scalar operations


template<typename T>
inline Matrix3<T> & Matrix3<T>::operator*=(
const T &scalar)
{
    for(int i=0 ; i<9 ; ++i)
        a[i] *= scalar;
    return *this;
}

template<typename T>
inline Matrix3<T> & Matrix3<T>::operator/=(
const T &scalar)
{
    for(int i=0;i<9;++i)
        a[i] /= scalar;
    return *this;
}


// matrix operations


template<typename T>
inline Matrix3<T> & Matrix3<T>::operator*=(
const Matrix3<T> &m)
{
    //double b[9];
    T b[9];
    b[0] = a[0]*m.a[0] + a[1]*m.a[3] + a[2]*m.a[6];
    b[1] = a[0]*m.a[1] + a[1]*m.a[4] + a[2]*m.a[7];
    b[2] = a[0]*m.a[2] + a[1]*m.a[5] + a[2]*m.a[8];
    b[3] = a[3]*m.a[0] + a[4]*m.a[3] + a[5]*m.a[6];
    b[4] = a[3]*m.a[1] + a[4]*m.a[4] + a[5]*m.a[7];
    b[5] = a[3]*m.a[2] + a[4]*m.a[5] + a[5]*m.a[8];
    b[6] = a[6]*m.a[0] + a[7]*m.a[3] + a[8]*m.a[6];
    b[7] = a[6]*m.a[1] + a[7]*m.a[4] + a[8]*m.a[7];
    b[8] = a[6]*m.a[2] + a[7]*m.a[5] + a[8]*m.a[8];

    set(b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8]);

    return *this;
}

template<typename T>  	  
inline Matrix3<T> & Matrix3<T>::operator+=(
const Matrix3<T> &m)
{
    for(int i=0;i<9;++i)
        a[i] += m.a[i];
    return *this;
}

template<typename T>
inline Matrix3<T> & Matrix3<T>::operator-=(
const Matrix3<T> &m)
{
    for(int i=0;i<9;++i)
        a[i] -= m.a[i];
    return *this;
}

template<typename T>
inline void Matrix3<T>::setToIdentity()
{
    a[0]=a[4]=a[8]=1;
    a[1]=a[2]=a[3]=a[5]=a[6]=a[7]=0;
}

template<typename T>
inline bool Matrix3<T>::isZero()
const
{
    return (a[0] == 0 && a[1] == 0 && a[2] == 0 &&
            a[3] == 0 && a[4] == 0 && a[5] == 0 &&
            a[6] == 0 && a[7] == 0 && a[8] == 0);
}

template<typename T>
inline void Matrix3<T>::transpose()
{
    Matrix3<T> m;// = *this;
    int i,j;
    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            m(i,j) = operator()(j,i);
    *this = m;
}


//////////////////////////////////////////////////////
//
//    non member operators implementation
//
//////////////////////////////////////////////////////



// matrix - matrix

template<typename T>
bool operator==(
const Matrix3<T> &m1,
const Matrix3<T> &m2)
{
    return m1(0,0) == m2(0,0) && m1(0,1) == m2(0,1) && m1(0,2) == m2(0,2) &&
           m1(1,0) == m2(1,0) && m1(1,1) == m2(1,1) && m1(1,2) == m2(1,2) &&
           m1(2,0) == m2(2,0) && m1(2,1) == m2(2,1) && m1(2,2) == m2(2,2);
}


template<typename T>
inline Matrix3<T> operator*(
const Matrix3<T> &m1,
const Matrix3<T> &m2)
{
    Matrix3<T> product;
    int i,j;

    for(i=0;i<3;++i)
        for(j=0;j<3;++j)
            product(i,j) = m1(i,0)*m2(0,j) + m1(i,1)*m2(1,j) + m1(i,2)*m2(2,j);
    
    return product;
}


template<typename T>  	  
inline const Matrix3<T> operator+(
const Matrix3<T> &matrix1,
const Matrix3<T> &matrix2)
{
    Matrix3<T> product;

    for(int i=0; i<3 ; i++)
        for(int j=0 ; j<3 ; j++)
            product(i,j) = matrix1(i,j) + matrix2(i,j);
    return product;
}


template<typename T>
inline const Matrix3<T> operator-(
const Matrix3<T> &matrix1,
const Matrix3<T> &matrix2)
{
    Matrix3<T> difference;

    for(int i=0 ; i<3 ; i++)
        for(int j=0 ; j<3 ; j++)
            difference(i,j) = matrix1(i,j) - matrix2(i,j);

    return difference;
}



// matrix - vector



template<typename T>
inline const Vector3<T> operator*(
const Matrix3<T> &matrix,
const Vector3<T> &vector)
{
    Vector3<T> product;
    int i,j;

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            product(i) += matrix(i,j) * vector(j);
    
    return product;
}


template<typename T>
const Vector3<T> operator*(
const Vector3<T> &vector,
const Matrix3<T> &matrix)
{
    Vector3<T> product;
    int i,j;

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            product(i) += matrix(j,i) * vector(j);

    return product;
}



// matrix - scalar



template<typename T>
inline const Matrix3<T> operator*(
const T &scalar,
const Matrix3<T> &matrix)
{
    Matrix3<T> product(matrix);
    
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            product(i,j) *= scalar;
    
    return product;
}


template<typename T>
inline const Matrix3<T> operator*(
const Matrix3<T> &matrix,
const T &scalar)
{
    Matrix3<T> product(matrix);
    
    for(int i=0 ; i<3 ; i++)
        for(int j=0 ; j<3 ; j++)
            product(i,j) *= scalar;
    
    return product;
}



template<typename T>
inline const Matrix3<T> operator/(
const Matrix3<T> &matrixIn,
const T &scalar)
{
    Matrix3<T> matrixOut(matrixIn);
    int i,j;

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            matrixOut(i,j) /= scalar;

    return matrixOut;
}




// mixed type products

template<typename T,typename U,typename W>
const Matrix3<T> product(const Matrix3<U> &m1,const Matrix3<W> &m2)
{
    Matrix3<T> result;
    int i,j;

    for(i=0;i<3;++i)
        for(j=0;j<3;++j)
            result(i,j) = static_cast<T>(m1(i,0)*m2(0,j) + m1(i,1)*m2(1,j) + m1(i,2)*m2(2,j));
    
    return result;
}

template<typename T,typename U,typename W>
const Vector3<T> product(const Matrix3<U> &matrix,const Vector3<W> &vector)
{
    Vector3<T> result;
    int i,j;

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            result(i) += static_cast<T>(matrix(i,j) * vector(j));
    
    return result;
}

template<typename T,typename U,typename W>
const Vector3<T> product(const Vector3<U> &vector,const Matrix3<W> &matrix)
{
    Vector3<T> result;
    int i,j;

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            result(i) += static_cast<T>(matrix(i,j) * vector(j));

    return result;
}

template<typename T,typename U,typename W>
const Matrix3<T> product(const U &scalar,const Matrix3<W> &matrix)
{
    Matrix3<T> result;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            result(i,j) = static_cast<T>(matrix(i,j) * scalar);
}

template<typename T,typename U,typename W>
const Matrix3<T> product(const Matrix3<U> &matrix,const W &scalar)
{
    Matrix3<T> result;
    for(int i=0;i<3;i++)
        for(int j=0;j<3;j++)
            result(i,j) = static_cast<T>(matrix(i,j) * scalar);
}

/**@}*/

} // namespace discamb







#endif /*_DISCAMB_MATHUTILITIES_MATRIX3_H_*/
