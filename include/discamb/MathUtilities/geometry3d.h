#include "discamb/MathUtilities/Vector3.h"
#include "discamb/MathUtilities/Matrix3.h"

#include<vector>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */


    namespace geometry3d {
        void normalize(Vector3d &v);
        Vector3d normalizeReturn(const Vector3d &in);
        //Vector3d rotate(const Vector3d &v, double angle, const Vector3d &axis);
        //Matrix3d rotationMatrix(double angle, const Vector3d &axis);
        double distance(const Vector3d &v1, const Vector3d &v2);
        double length(const Vector3d &v);

        double angle(const Vector3d &v1, const Vector3d &v2);
        double angle(const Vector3d& p1, const Vector3d& p2, const Vector3d& p3);
        
        // angle between two vectors v1-v2 and v3-v2, or an angle between three atoms at 
        // positions v1, v2 and v3, where v2 is a position of the middle atom
        //double angle(const Vector3d &v1, const Vector3d &v2, const Vector3d &v3);

        //double angleAroundAxis(const Vector3d &v1, const Vector3d &v2, const Vector3d &axis);
        double angleAroundAxis2(const Vector3d &v1, const Vector3d &v2, const Vector3d &axis);
        // finds operation {m|t] transforming vectors v2 into v1 
        //void findTransformation(const std::vector<Vector3d> &v1, const std::vector<Vector3d> &v2,
        //                        Matrix3d &m, Vector3d &t);
    }
    /** @}*/
}