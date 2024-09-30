#ifndef _DISCAMB_MATHUTILITIES_VECTOR3_H_
#define _DISCAMB_MATHUTILITIES_VECTOR3_H_

#include <iostream>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


/** Represents a vector in 3D space.*/


template<typename T>
class Vector3{
public:

    T x,y,z;

    Vector3();
    Vector3(const T &x,const T &y,const T &z);
    Vector3(const Vector3<T> &v);

    Vector3<T> & operator=(const Vector3<T> &v);

    template<typename INDEXABLE>
    Vector3(const INDEXABLE &v);

    template<typename INDEXABLE>
    Vector3<T> & operator=(const INDEXABLE &v);


    ~Vector3();

    template<typename INDEX_TYPE>
    const T& operator()(INDEX_TYPE) const;

    template<typename INDEX_TYPE>
    T& operator()(INDEX_TYPE);

    template<typename INDEX_TYPE>
    const T& operator[](INDEX_TYPE) const;

    template<typename INDEX_TYPE>
    T& operator[](INDEX_TYPE);


    const T &operator()(int) const;
    T &operator()(int);
    const T &operator[](int) const;
    T &operator[](int);



    void set(const T &x,const T &y,const T &z);
    void get(T &x,T &y,T &z) const;
  
    Vector3<T> & operator+=(const Vector3<T> &v);
    Vector3<T> & operator-=(const Vector3<T> &v);

    Vector3<T> & operator*=(const T &scalar); 
    Vector3<T> & operator/=(const T &scalar);  	  

};	

template<typename T>
std::ostream& operator<<(std::ostream& out, const Vector3<T>& v);

template<typename T>
std::istream& operator>>(std::istream& out, Vector3<T>& v);

// 

typedef Vector3<double> Vector3d; 
typedef Vector3<float> Vector3f;
typedef Vector3<int> Vector3i;

// non-member operators

template <typename T>
Vector3<T> operator-(const Vector3<T> &v);

template <typename T>
const Vector3<T> operator+(const Vector3<T> &u,const Vector3<T> &v);

template <typename T>
const Vector3<T> operator-(const Vector3<T> &u,const Vector3<T> &v);

template <typename T>
const T operator*(const Vector3<T> &u,const Vector3<T> &v);

template <typename T>
const Vector3<T> operator*(const T &scalar,const Vector3<T> &v);

template <typename T>
const Vector3<T> operator*(const Vector3<T> &v,const T &scalar);

/** cross product */
template <typename T>
const Vector3<T> operator^(const Vector3<T> &u,const Vector3<T> &v);

template <typename T>
const Vector3<T> operator/(const Vector3<T> &v1,const T &scalar);

template <typename T>
bool operator<(const Vector3<T> &v1,const Vector3<T> &v2);

template <typename T>
bool operator<=(const Vector3<T> &v1,const Vector3<T> &v2);

template <typename T>
bool operator==(const Vector3<T> &v1,const Vector3<T> &v2);

template <typename T>
bool operator!=(const Vector3<T> &v1,const Vector3<T> &v2);

template <typename T>
bool operator>=(const Vector3<T> &v1,const Vector3<T> &v2);

template <typename T>
bool operator>(const Vector3<T> &v1,const Vector3<T> &v2);

template <typename T>
const Vector3<T> cross_product(const Vector3<T> &v1,const Vector3<T> &v2);

///////////////////////////////////////////////////////////////////
//                      IMPLEMENTATION                           //
///////////////////////////////////////////////////////////////////


template<typename T>
inline Vector3<T>::Vector3() : x(0),y(0),z(0) 
{
}

template<typename T>
Vector3<T>::Vector3(
const T &_x,
const T &_y,
const T &_z) : x(_x),y(_y),z(_z)
{
}

template<typename T>
inline Vector3<T>::Vector3(
const Vector3<T> &vec) : x(vec.x),y(vec.y),z(vec.z)
{
}


template<typename T>
inline Vector3<T> & Vector3<T>::operator=(
const Vector3<T> &vec)
{
    x = vec.x;
    y = vec.y;
    z = vec.z;
    return *this;
}

template<typename T>
template<typename INDEXABLE>
Vector3<T>::Vector3(const INDEXABLE &v)
{
    x = v[0];
    y = v[1];
    z = v[2];
}

template<typename T>
template<typename INDEXABLE>
Vector3<T> & Vector3<T>::operator=(const INDEXABLE &v)
{
    x = static_cast<T>(v[0]);
    y = static_cast<T>(v[1]);
    z = static_cast<T>(v[2]);
    return *this;
}


template<typename T>
inline Vector3<T>::~Vector3(){}

// member operators


template<typename T>
const T &Vector3<T>::operator()(
int i)
const
{
    if(i==0)
        return x;
    else
        if(i==1)
            return y;
        else
            return z;
}

template<typename T>
T &Vector3<T>::operator()(
int i)
{
    if(i==0)
        return x;
    else
        if(i==1)
            return y;
        else
            return z;
}

template<typename T>
const T &Vector3<T>::operator[](
int i)
const
{
    if(i==0)
        return x;
    else
        if(i==1)
            return y;
        else
            return z;
}

template<typename T>
T &Vector3<T>::operator[](
int i)
{
    if(i==0)
        return x;
    else
        if(i==1)
            return y;
        else
            return z;
}

template<typename T>
template<typename INDEX_TYPE>
const T& Vector3<T>::operator()(INDEX_TYPE i)
const
{
    return operator()(static_cast<int>(i));
}

template<typename T>
template<typename INDEX_TYPE>
T& Vector3<T>::operator()(INDEX_TYPE i)
{
    return operator()(static_cast<int>(i));
}

template<typename T>
template<typename INDEX_TYPE>
const T& Vector3<T>::operator[](INDEX_TYPE i)
const
{
    return operator[](static_cast<int>(i));
}

template<typename T>
template<typename INDEX_TYPE>
T& Vector3<T>::operator[](INDEX_TYPE i)
{
    return operator[](static_cast<int>(i));
}


template<typename T>  	  
inline Vector3<T> & Vector3<T>::operator+=(
const Vector3<T> &v)
{
    x += v.x;
    y += v.y;
    z += v.z;

    return *this;
}

template<typename T>
inline Vector3<T> & Vector3<T>::operator-=(
const Vector3<T> &v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;

    return *this;
}

template<typename T>
inline Vector3<T> & Vector3<T>::operator*=(
const T &scalar)
{
    x *= scalar;
    y *= scalar;
    z *= scalar;

    return *this;
}

template<typename T>
inline Vector3<T> & Vector3<T>::operator/=(
const T &scalar)
{
    x /= scalar;
    y /= scalar;
    z /= scalar;
    return *this;
}

template<typename T>
inline void Vector3<T>::set(
const T &_x,
const T &_y,
const T &_z)
{
    x = _x;
    y = _y;
    z = _z;
}

template<typename T>
inline void Vector3<T>::get(
T &_x,
T &_y,
T &_z) 
const
{
    _x = x;
    _y = y;
    _z = z;
}


// non-member operators

template<typename T>
inline Vector3<T> operator-(const Vector3<T> &v)
{
    return Vector3<T>(-v.x,-v.y,-v.z);
}

template<typename T>
inline const Vector3<T> operator+(
const Vector3<T> &u,
const Vector3<T> &v)
{
    return Vector3<T>( u.x + v.x , u.y + v.y , u.z + v.z);
}

template<typename T>
inline const Vector3<T> operator-(
const Vector3<T> &u,
const Vector3<T> &v)
{
return Vector3<T>( u.x - v.x , u.y - v.y , u.z - v.z);
}



template<typename T>
inline const T operator*(
const Vector3<T> &u,
const Vector3<T> &v)
{
    return v.x*u.x + v.y*u.y + v.z*u.z;
}

template<typename T>
inline const Vector3<T> operator^(
const Vector3<T> &u,
const Vector3<T> &v)
{
    return Vector3<T>(
        u.y*v.z - u.z*v.y ,
        u.z*v.x - u.x*v.z ,
        u.x*v.y - u.y*v.x );
}



template<typename T>
inline const Vector3<T> operator*(
const T &scalar,
const Vector3<T> &v)
{
    return Vector3<T>(scalar*v.x,scalar*v.y,scalar*v.z);
}

template<typename T>
inline const Vector3<T> operator*(
const Vector3<T> &v,
const T &scalar)
{
    return Vector3<T>( v.x*scalar , v.y*scalar , v.z*scalar);
}



template<typename T>
inline const Vector3<T> operator/(
const Vector3<T> &v,
const T &scalar)
{
    return Vector3<T>( v.x/scalar , v.y/scalar , v.z/scalar );
}



template <typename T>
inline bool operator<(const Vector3<T> &v1,const Vector3<T> &v2)
{
    if(v1(0)<v2(0))
        return true;
    if(v1(0)>v2(0))
        return false;
    if(v1(1)<v2(1))
        return true;
    if(v1(1)>v2(1))
        return false;
    if(v1(2)<v2(2))
        return true;
    return false;
}

template <typename T>
inline bool operator<=(const Vector3<T> &v1,const Vector3<T> &v2)
{
    return (v1<v2)||(v1==v2);
}

template <typename T>
inline bool operator==(const Vector3<T> &v1,const Vector3<T> &v2)
{
    return ( v1(0)==v2(0) && v1(1)==v2(1) && v1(2)==v2(2) );
}

template <typename T>
inline bool operator!=(const Vector3<T> &v1,const Vector3<T> &v2)
{
    return !(v1==v2);
}

template <typename T>
inline bool operator>=(const Vector3<T> &v1,const Vector3<T> &v2)
{
    return !(v1<v2);
}

template <typename T>
inline bool operator>(const Vector3<T> &v1,const Vector3<T> &v2)
{
    return !(v1<=v2);
}

template <typename T>
inline const Vector3<T> cross_product(
const Vector3<T> &u,
const Vector3<T> &v)
{
    return Vector3<T>(u[1]*v[2]-u[2]*v[1],
                      u[2]*v[0]-u[0]*v[2],
                      u[0]*v[1]-u[1]*v[0]);
}

template<typename T>
std::ostream& operator<<(
    std::ostream& out,
    const Vector3<T>& v)
{
    out << v.x << " " << v.y << " " << v.z;
    return out;
}

template<typename T>
std::istream& operator>>(
    std::istream& inp,
    Vector3<T>& v)
{
    inp >> v.x >> v.y >> v.z;
    return inp;
}

/** @}*/
} // namespace discamb







#endif /*_DISCAMB_MATHUTILITIES_VECTOR3_H_*/
