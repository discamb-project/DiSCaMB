#pragma once

#include "discamb/MathUtilities/Matrix3.h"
#include <cmath>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    namespace algebra3d {

        template<typename T>
        void eigensystemRealSymm(const Matrix3<T>& matrix, Vector3<T>& v1, Vector3<T>& v2, Vector3<T>& v3, T& l1, T& l2, T& l3);

        template<typename T>
        void invert3d(const Matrix3<T>& matrix, Matrix3<T>& invertedMatrix);

        template<typename T>
        Matrix3<T> transpose3d(const Matrix3<T>& matrix);


        template<typename T>
        void transpose3d(const Matrix3<T>& matrix, Matrix3<T>& transposed);

        void powerRealSymm(const Matrix3d& m, double p, Matrix3d& result);

        template<typename T>
        T det3d(const Matrix3<T>& matrix);

        template<typename T>
        T trace(const Matrix3<T>& matrix);

        template<typename T>
        Matrix3<T> outer(const Vector3<T>& v1, const Vector3<T>& v2);

        namespace auxiliary {
            template<typename T> void eigen_decomposition(T A[3][3], T V[3][3], T d[3]);
        }

    }
// IMPLEMENTATION

template<typename T>
T algebra3d::det3d(const Matrix3<T> &m)
{
    return    m(0,0)*m(1,1)*m(2,2) + m(0,1)*m(1,2)*m(2,0) + m(0,2)*m(1,0)*m(2,1)
            - m(0,0)*m(1,2)*m(2,1) - m(0,1)*m(1,0)*m(2,2) - m(0,2)*m(1,1)*m(2,0);
}

template<typename T>
T algebra3d::trace(const Matrix3<T>& m)
{
    return m(0, 0) + m(1, 1) + m(2, 2);
}


template<typename T>
void algebra3d::invert3d(const Matrix3<T> &m,Matrix3<T> &im)
{
 
     im(0,0) = m(1,1)*m(2,2) - m(1,2)*m(2,1);
     im(0,1) = m(0,2)*m(2,1) - m(0,1)*m(2,2);
     im(0,2) = m(0,1)*m(1,2) - m(0,2)*m(1,1);
     im(1,0) = m(1,2)*m(2,0) - m(1,0)*m(2,2);
     im(1,1) = m(0,0)*m(2,2) - m(0,2)*m(2,0);
     im(1,2) = m(0,2)*m(1,0) - m(0,0)*m(1,2);
     im(2,0) = m(1,0)*m(2,1) - m(1,1)*m(2,0);
     im(2,1) = m(0,1)*m(2,0) - m(0,0)*m(2,1);
     im(2,2) = m(0,0)*m(1,1) - m(0,1)*m(1,0);

     im /= det3d(m);
}


template<typename T>
void algebra3d::eigensystemRealSymm(
const Matrix3<T> &matrix,
Vector3<T> &v1,
Vector3<T> &v2,
Vector3<T> &v3,
T &l1,
T &l2,
T &l3)
{
    T A[3][3];
    T V[3][3];
    T d[3];

    int i,j;

    for(i=0;i<3;i++)
        for(j=0;j<3;j++)
            A[i][j] = matrix(i,j);

    auxiliary::eigen_decomposition(A,V,d);

    for(i=0;i<3;i++)
    {
        v1(i) = V[i][0];
        v2(i) = V[i][1];
        v3(i) = V[i][2];
    }

    l1 = d[0];
    l2 = d[1];
    l3 = d[2];
}


// ------ COPIED FROM EIG3 MODIFIED FOR WORK WITH TEMPLATES

namespace algebra3d
{
    namespace auxiliary
    {
    
    
        /* Eigen decomposition code for symmetric 3x3 matrices, copied from the public
           domain Java Matrix library JAMA. */



        //#ifdef MAX
        //#undef MAX
        //#endif

        //#define MAX(a, b) ((a)>(b)?(a):(b))

        //#define n 3
        template<typename T>
        static T hypot2(T x, T y) {
          return sqrt(x*x+y*y);
        }

        // Symmetric Householder reduction to tridiagonal form.
        template<typename T>
        static void tred2(T V[3][3], T d[3], T e[3]) {

        //  This is derived from the Algol procedures tred2 by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.

          for (int j = 0; j < 3; j++) {
            d[j] = V[3-1][j];
          }

          // Householder reduction to tridiagonal form.

          for (int i = 3-1; i > 0; i--) {

            // Scale to avoid under/overflow.

            T scale = 0.0;
            T h = 0.0;
            for (int k = 0; k < i; k++) {
              scale = scale + fabs(d[k]);
            }
            if (scale == 0.0) {
              e[i] = d[i-1];
              for (int j = 0; j < i; j++) {
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
                V[j][i] = 0.0;
              }
            } else {

              // Generate Householder vector.

              for (int k = 0; k < i; k++) {
                d[k] /= scale;
                h += d[k] * d[k];
              }
              T f = d[i-1];
              T g = sqrt(h);
              if (f > 0) {
                g = -g;
              }
              e[i] = scale * g;
              h = h - f * g;
              d[i-1] = f - g;
              for (int j = 0; j < i; j++) {
                e[j] = 0.0;
              }

              // Apply similarity transformation to remaining columns.

              for (int j = 0; j < i; j++) {
                f = d[j];
                V[j][i] = f;
                g = e[j] + V[j][j] * f;
                for (int k = j+1; k <= i-1; k++) {
                  g += V[k][j] * d[k];
                  e[k] += V[k][j] * f;
                }
                e[j] = g;
              }
              f = 0.0;
              for (int j = 0; j < i; j++) {
                e[j] /= h;
                f += e[j] * d[j];
              }
              T hh = f / (h + h);
              for (int j = 0; j < i; j++) {
                e[j] -= hh * d[j];
              }
              for (int j = 0; j < i; j++) {
                f = d[j];
                g = e[j];
                for (int k = j; k <= i-1; k++) {
                  V[k][j] -= (f * e[k] + g * d[k]);
                }
                d[j] = V[i-1][j];
                V[i][j] = 0.0;
              }
            }
            d[i] = h;
          }

          // Accumulate transformations.

          for (int i = 0; i < 3-1; i++) {
            V[3-1][i] = V[i][i];
            V[i][i] = 1.0;
            double h = d[i+1];
            if (h != 0.0) {
              for (int k = 0; k <= i; k++) {
                d[k] = V[k][i+1] / h;
              }
              for (int j = 0; j <= i; j++) {
                T g = 0.0;
                for (int k = 0; k <= i; k++) {
                  g += V[k][i+1] * V[k][j];
                }
                for (int k = 0; k <= i; k++) {
                  V[k][j] -= g * d[k];
                }
              }
            }
            for (int k = 0; k <= i; k++) {
              V[k][i+1] = 0.0;
            }
          }
          for (int j = 0; j < 3; j++) {
            d[j] = V[3-1][j];
            V[3-1][j] = 0.0;
          }
          V[3-1][3-1] = 1.0;
          e[0] = 0.0;
        } 

        // Symmetric tridiagonal QL algorithm.
        template<typename T>
        static void tql2(T V[3][3], double d[3], double e[3]) {

        //  This is derived from the Algol procedures tql2, by
        //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
        //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
        //  Fortran subroutine in EISPACK.

          for (int i = 1; i < 3; i++) {
            e[i-1] = e[i];
          }
          e[3-1] = 0.0;

          T f = 0.0;
          T tst1 = 0.0;
          T eps = pow(2.0,-52.0);
          for (int l = 0; l < 3; l++) {

            // Find small subdiagonal element

              //MAX(a, b) ((a)>(b)?(a):(b))

            tst1 = tst1 > fabs(d[l]) + fabs(e[l]) ? tst1 : fabs(d[l]) + fabs(e[l]);
            int m = l;
            while (m < 3) {
              if (fabs(e[m]) <= eps*tst1) {
                break;
              }
              m++;
            }

            // If m == l, d[l] is an eigenvalue,
            // otherwise, iterate.

            if (m > l) {
              int iter = 0;
              do {
                iter = iter + 1;  // (Could check iteration count here.)

                // Compute implicit shift

                T g = d[l];
                T p = (d[l+1] - g) / (2.0 * e[l]);
                T r = hypot2(p,1.0);
                if (p < 0) {
                  r = -r;
                }
                d[l] = e[l] / (p + r);
                d[l+1] = e[l] * (p + r);
                T dl1 = d[l+1];
                T h = g - d[l];
                for (int i = l+2; i < 3; i++) {
                  d[i] -= h;
                }
                f = f + h;

                // Implicit QL transformation.

                p = d[m];
                T c = 1.0;
                T c2 = c;
                T c3 = c;
                T el1 = e[l+1];
                T s = 0.0;
                T s2 = 0.0;
                for (int i = m-1; i >= l; i--) {
                  c3 = c2;
                  c2 = c;
                  s2 = s;
                  g = c * e[i];
                  h = c * p;
                  r = hypot2(p,e[i]);
                  e[i+1] = s * r;
                  s = e[i] / r;
                  c = p / r;
                  p = c * d[i] - s * g;
                  d[i+1] = h + s * (c * g + s * d[i]);

                  // Accumulate transformation.

                  for (int k = 0; k < 3; k++) {
                    h = V[k][i+1];
                    V[k][i+1] = s * V[k][i] + c * h;
                    V[k][i] = c * V[k][i] - s * h;
                  }
                }
                p = -s * s2 * c3 * el1 * e[l] / dl1;
                e[l] = s * p;
                d[l] = c * p;

                // Check for convergence.

              } while (fabs(e[l]) > eps*tst1);
            }
            d[l] = d[l] + f;
            e[l] = 0.0;
          }
  
          // Sort eigenvalues and corresponding vectors.

          for (int i = 0; i < 3-1; i++) {
            int k = i;
            T p = d[i];
            for (int j = i+1; j < 3; j++) {
              if (d[j] < p) {
                k = j;
                p = d[j];
              }
            }
            if (k != i) {
              d[k] = d[i];
              d[i] = p;
              for (int j = 0; j < 3; j++) {
                p = V[j][i];
                V[j][i] = V[j][k];
                V[j][k] = p;
              }
            }
          }
        }

        template<typename T>
        void eigen_decomposition(T A[3][3], T V[3][3], T d[3]) {
          T e[3];
          for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
              V[i][j] = A[i][j];
            }
          }
          tred2(V, d, e);
          tql2(V, d, e);
        }


    }
}

/** @}*/

template<typename T>
Matrix3<T> algebra3d::transpose3d(
const Matrix3<T> &m)
{
    // Matrix3(const T &a11,const T &a12,const T &a13,const T &a21,const T &a22,const T &a23,const T &a31,const T &a32,const T &a33);
    return Matrix3<T>(m(0,0), m(1,0), m(2,0),
                      m(0,1), m(1,1), m(2,1),
                      m(0,2), m(1,2), m(2,2));
}

template<typename T>
void algebra3d::transpose3d(const Matrix3<T> &m,Matrix3<T> &t)
{
    t.set( m(0,0), m(1,0), m(2,0),
           m(0,1), m(1,1), m(2,1),
           m(0,2), m(1,2), m(2,2));
}

template<typename T>
Matrix3<T> algebra3d::outer(
    const Vector3<T>& v1,
    const Vector3<T>& v2)
{
    Matrix3<T> result;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            result(i, j) = v1(i) * v2(j);
    return result;
}


}
