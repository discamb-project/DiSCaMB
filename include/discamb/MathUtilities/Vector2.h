#ifndef _DISCAMB_MATHUTILITIES_VECTOR2_H_
#define _DISCAMB_MATHUTILITIES_VECTOR2_H_

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


/** Represents a vector in 2D space.*/


template<typename T>
class Vector2{
public:

    T x;
    T y;

    Vector2();
    Vector2(const T &x,const T &y);
    Vector2(const Vector2<T> &v);

    Vector2<T> & operator=(const Vector2<T> &v);

    template<typename INDEXABLE>
    Vector2(const INDEXABLE &v);

    template<typename INDEXABLE>
    Vector2<T> & operator=(const INDEXABLE &v);


    ~Vector2();

    const T &operator()(int) const;
    T &operator()(int);
    const T &operator[](int) const;
    T &operator[](int);


    void set(const T &x,const T &y);
    void get(T &x,T &y) const;
  
    Vector2<T> & operator+=(const Vector2<T> &v);
    Vector2<T> & operator-=(const Vector2<T> &v);

    Vector2<T> & operator*=(const T &scalar); 
    Vector2<T> & operator/=(const T &scalar);  	  

//    auto operator<=>(const Point&) const = default;
};

// 

typedef Vector2<double> Vector2d; 
typedef Vector2<float> Vector2f;
typedef Vector2<int> Vector2i;

// non-member operators

template <typename T>
Vector2<T> operator-(const Vector2<T> &v);

template <typename T>
const Vector2<T> operator+(const Vector2<T> &u,const Vector2<T> &v);

template <typename T>
const Vector2<T> operator-(const Vector2<T> &u,const Vector2<T> &v);

template <typename T>
const T operator*(const Vector2<T> &u,const Vector2<T> &v);

template <typename T>
const Vector2<T> operator*(const T &scalar,const Vector2<T> &v);

template <typename T>
const Vector2<T> operator*(const Vector2<T> &v,const T &scalar);


template <typename T>
const Vector2<T> operator/(const Vector2<T> &v1,const T &scalar);


///////////////////////////////////////////////////////////////////
//                      IMPLEMENTATION                           //
///////////////////////////////////////////////////////////////////


template<typename T>
inline Vector2<T>::Vector2() : x(0),y(0)
{
}

template<typename T>
Vector2<T>::Vector2(
const T &_x,
const T &_y) : x(_x),y(_y)
{
}

template<typename T>
inline Vector2<T>::Vector2(
const Vector2<T> &vec) : x(vec.x),y(vec.y)
{
}


template<typename T>
inline Vector2<T> & Vector2<T>::operator=(
const Vector2<T> &vec)
{
    x = vec.x;
    y = vec.y;
    return *this;
}

template<typename T>
template<typename INDEXABLE>
Vector2<T>::Vector2(const INDEXABLE &v)
{
    x = v[0];
    y = v[1];
}

template<typename T>
template<typename INDEXABLE>
Vector2<T> & Vector2<T>::operator=(const INDEXABLE &v)
{
    x = v[0];
    y = v[1];
    return *this;
}


template<typename T>
inline Vector2<T>::~Vector2(){}

// member operators

template<typename T>
const T &Vector2<T>::operator()(
int i)
const
{
    if(i==0)
        return x;
    else
        return y;
}

template<typename T>
T &Vector2<T>::operator()(
int i)
{
    if(i==0)
        return x;
    else
        return y;
}

template<typename T>
const T &Vector2<T>::operator[](
int i)
const
{
    if(i==0)
        return x;
    else
        return y;
}

template<typename T>
T &Vector2<T>::operator[](
int i)
{
    if(i==0)
        return x;
    else
        return y;
}


template<typename T>  	  
inline Vector2<T> & Vector2<T>::operator+=(
const Vector2<T> &v)
{
    x += v.x;
    y += v.y;

    return *this;
}

template<typename T>
inline Vector2<T> & Vector2<T>::operator-=(
const Vector2<T> &v)
{
    x -= v.x;
    y -= v.y;

    return *this;
}

template<typename T>
inline Vector2<T> & Vector2<T>::operator*=(
const T &scalar)
{
    x *= scalar;
    y *= scalar;

    return *this;
}

template<typename T>
inline Vector2<T> & Vector2<T>::operator/=(
const T &scalar)
{
    x /= scalar;
    y /= scalar;
    return *this;
}

template<typename T>
inline void Vector2<T>::set(
const T &_x,
const T &_y)
{
    x = _x;
    y = _y;
}

template<typename T>
inline void Vector2<T>::get(
T &_x,
T &_y) 
const
{
    _x = x;
    _y = y;
}


// non-member operators

template<typename T>
inline Vector2<T> operator-(const Vector2<T> &v)
{
    return Vector2<T>(-v.x,-v.y);
}

template<typename T>
inline const Vector2<T> operator+(
const Vector2<T> &u,
const Vector2<T> &v)
{
    return Vector2<T>( u.x + v.x , u.y + v.y);
}

template<typename T>
inline const Vector2<T> operator-(
const Vector2<T> &u,
const Vector2<T> &v)
{
return Vector2<T>( u.x - v.x , u.y - v.y);
}



template<typename T>
inline const T operator*(
const Vector2<T> &u,
const Vector2<T> &v)
{
    return v.x*u.x + v.y*u.y;
}



template<typename T>
inline const Vector2<T> operator*(
const T &scalar,
const Vector2<T> &v)
{
    return Vector2<T>(scalar*v.x,scalar*v.y);
}

template<typename T>
inline const Vector2<T> operator*(
const Vector2<T> &v,
const T &scalar)
{
    return Vector2<T>( v.x*scalar , v.y*scalar);
}



template<typename T>
inline const Vector2<T> operator/(
const Vector2<T> &v,
const T &scalar)
{
    return Vector2<T>( v.x/scalar , v.y/scalar);
}


/** @}*/
} // namespace discamb







#endif /*_DISCAMB_MATHUTILITIES_VECTOR2_H_*/
