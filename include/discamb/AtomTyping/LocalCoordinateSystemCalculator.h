#pragma once

#include "LocalCoordinateSystem.h"


#include "discamb/CrystalStructure/SpaceGroupOperation.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"
#include "discamb/MathUtilities/CrossProductLcs.h"

#include <vector>
#include <memory>
#include <string>

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

	class LocalCoordinateSystemCalculator: public LocalCoordinateSystemInCrystal
	{
	public:
		LocalCoordinateSystemCalculator();
		LocalCoordinateSystemCalculator(const LocalCoordinateSystem<AtomInCrystalID> &lcs, const Crystal &c);
		~LocalCoordinateSystemCalculator();
		// C(1) C(2),C(3),C(4)[-x,-y,z] X any_orthogonal Y R
		// C(1) average_direction:C(2)[1-x,y,z],C(3) Z C(4) Y L
		virtual void set(const std::string &definition, const Crystal &c);
		void set(const LocalCoordinateSystem<AtomInCrystalID> &lcs, const Crystal &c);
		virtual void calculate(Vector3d &x, Vector3d &y, Vector3d &z, const Crystal &c) const;
        void calculate(Vector3d &x, Vector3d &y, Vector3d &z, const Crystal &c, bool &sameChirality) const;
        // cosAngle is a cosine of the angle between the two directions defining lcs
        void calculate(Vector3d& x, Vector3d& y, Vector3d& z, const Crystal& c, bool& sameChirality, double &cosAngle) const;

        /*
              | x[0] y[0] z[0] |
          m = | x[1] y[1] z[1] |
              | x[2] y[2] z[2] |
          x - new x direction, y - new y direction, z - new z direction
        */
        void calculate(Matrix3d &m, const Crystal &c, bool &sameChirality) const;
		using LocalCoordinateSystemInCrystal::calculate;
    private:
        mutable bool mIsCartesian;
        
        std::shared_ptr<CrossProductLcs> mCrossProductLcs;
        LocalCoordinateSystem<AtomInCrystalID> mLcs;
        int mThirdCoordinate;
        void calcAtomicPositionsBasedDirection(const std::vector<AtomInCrystalID> &refPoint, const AtomInCrystalID &centralAtom, const Crystal &c,
                                               LcsDirectionType type, Vector3d &direction) const;
        static void calcAtomPosition(const AtomInCrystalID &atom, const Crystal &c, Vector3d &r);
        static void calculateAnyOrthogonal(const Vector3d &r0, Vector3d &r);
        void takeChiralityIntoAccount(Vector3d &x, Vector3d &y, Vector3d &z, const Crystal &c, bool &sameChirality) const;
	};

    /**@}*/
}

