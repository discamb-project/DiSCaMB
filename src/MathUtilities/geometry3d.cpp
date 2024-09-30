#include "discamb/MathUtilities/geometry3d.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/BasicUtilities/on_error.h"

//#include "glm/glm.hpp"
//#include "glm/gtx/vector_angle.hpp"
//#include "glm/gtx/transform.hpp"

#include <cmath>
#include "discamb/MathUtilities/MathUtilities.h"

using namespace std;

namespace discamb {

    namespace geometry3d {

        void normalize(
            Vector3d &v)
        {
            v /= sqrt(v*v);
        }

        /*Vector3d rotate(
            const Vector3d &v,
            double angle, 
            const Vector3d &axis)
        {
            glm::tmat4x4<double> ;
            glm::dmat4 rotation;
            x = glm::rotate<double>(angle, glm::dvec3(axis.x, axis.y, axis.z));

        }*/
        
        //Matrix3d rotationMatrix(
        //    double angle,
        //    const Vector3d &axis)
        //{
        //    glm::dmat4 r;
        //    Vector3d a(normalizeReturn(axis));
        //    r = glm::rotate<double, glm::precision::highp>(angle, glm::dvec3(a.x, a.y, a.z));
        //    return Matrix3d(r[0][0], r[0][1], r[0][2], 
        //                    r[1][0], r[1][1], r[1][2], 
        //                    r[2][0], r[2][1], r[2][2] );
        //}

        Vector3d  normalizeReturn(
            const Vector3d &v)
        {
            return v / sqrt(v*v);
        }



        double distance(
            const Vector3d &v1, 
            const Vector3d &v2)
        {
            return 0.0;
        }

        double length(
            const Vector3d &v)
        {
            return std::sqrt(v*v);
        }

        double angle(
            const Vector3d& v1,
            const Vector3d& v2)
        {
            double v1_2 = v1 * v1;
            double v2_2 = v2 * v2;
            double denominator = sqrt(v1_2 * v2_2);
            if(denominator==0)
                on_error::throwException("An attempt to calcule angle between vetor(s) of length zero",
                                         __FILE__, __LINE__);

            double cos_alpha = v1 * v2 / denominator;
            if (cos_alpha > 1.0)
                cos_alpha = 1.0;
            if (cos_alpha < -1.0)
                cos_alpha = -1.0;
            double alpha = 180.0 / M_PI * acos(cos_alpha);
            return alpha;
        }

        double angle(
            const Vector3d& p1,
            const Vector3d& p2,
            const Vector3d& p3)
        {
            Vector3d v1, v2;
            v1 = p1 - p2;
            v2 = p3 - p2;
            return angle(v1, v2);
        }


        //double angle(
        //    const Vector3d &_v1, 
        //    const Vector3d &_v2)
        //{
        //    Vector3d v1(normalizeReturn(_v1)), v2(normalizeReturn(_v2));

        //    glm::dvec3 gv1(v1.x, v1.y, v1.z);
        //    glm::dvec3 gv2(v2.x, v2.y, v2.z);


        //    if (glm::angle(gv1, gv2) < 1e-8)
        //        return 0;
        //    
        //    Vector3d ref = normalizeReturn(v1^v2);
        //    glm::dvec3 gref(ref.x, ref.y, ref.z);
        //    return glm::orientedAngle(gv1, gv2, gref);
        //}

        //double angle(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3)
        //{
        //    return angle(v1 - v2, v3 - v2);
        //}

        //double angleAroundAxis(const Vector3d &v1, const Vector3d &v2, const Vector3d &axis)
        //{
        //    Vector3d x, y, a;
        //    
        //    x = normalizeReturn(v1);
        //    y = normalizeReturn(v2);
        //    a = normalizeReturn(axis);

        //    return  glm::orientedAngle(glm::dvec3(x.x, x.y, x.z),
        //                               glm::dvec3(y.x, y.y, y.z),
        //                               glm::dvec3(a.x, a.y, a.z));
        //}

        double angleAroundAxis2(const Vector3d &v1, const Vector3d &v2, const Vector3d &axis)
        {
            Vector3d v1n, v2n, x,y,z, v1_2d, v2_2d;

            v1n = normalizeReturn(v1);
            v2n = normalizeReturn(v2);
            z = normalizeReturn(axis);

            //

            y = normalizeReturn(z ^ v1n);
            x = y^z;

            v1_2d[0] = v1n*x;
            v1_2d[1] = v1n*y;
            normalize(v1_2d);
            v2_2d[0] = v2n*x;
            v2_2d[1] = v2n*y;
            normalize(v2_2d);

            return  atan2(v2_2d.y, v2_2d.x);
        }


        //void findTransformation(
        //    const std::vector<Vector3d> &_g1,
        //    const std::vector<Vector3d> &_g2,
        //    Matrix3d &m,
        //    Vector3d &t)
        //{
        //    if ((_g1.size() != 3) || (_g2.size() != 3))
        //        on_error::throwException("at least one input vector has invalid size", __FILE__, __LINE__);

        //    vector<Vector3d> g1 = _g1;
        //    vector<Vector3d> g2 = _g2;

        //    Vector3d t01 = -g1[0];
        //    Vector3d t02 = -g2[0];

        //    Matrix3d rot1, rot2;

        //    for (auto &v : g1)
        //        v += t01;

        //    for (auto &v : g2)
        //        v += t02;

        //    // already overlap?

        //    // skipping the check

        //    // 

        //    double rotAngle1 = geometry3d::angle(g1[1], g2[1]);

        //    rot1 = geometry3d::rotationMatrix(rotAngle1, g1[1] ^ g2[1]);

        //    //

        //    g2[1] = rot1*g2[1];
        //    g2[2] = rot1*g2[2];

        //    //

        //    //double rotAngle2 = geometry3d::angle(g1[2], g2[2]);


        //    double rotAngle2 = geometry3d::angleAroundAxis(g1[2], g2[2], g2[1]);
        //    rotAngle2 = geometry3d::angleAroundAxis2(g1[2], g2[2], g2[1]);

        //    rot2 = geometry3d::rotationMatrix(rotAngle2, g2[1]);

        //    //

        //    g2[2] = rot2*g2[2];

        //    m = rot2*rot1;
        //    t = m*t02 - t01;

        //}


    }

}