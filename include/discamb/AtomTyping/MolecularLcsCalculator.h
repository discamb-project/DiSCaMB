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

	class MolecularLcsCalculator
	{
	public:
        MolecularLcsCalculator();
        MolecularLcsCalculator(const LocalCoordinateSystem<int> &lcs);
		~MolecularLcsCalculator();
		void set(const LocalCoordinateSystem<int>& lcs);
		void calculate(Vector3d &x, Vector3d &y, Vector3d &z,const std::vector<Vector3d> &r) const;
        void calculate(Vector3d &x, Vector3d &y, Vector3d &z, const std::vector<Vector3d>& r, bool &sameChirality) const;
        /*
              | x[0] y[0] z[0] |
          m = | x[1] y[1] z[1] |
              | x[2] y[2] z[2] |
          x - new x direction, y - new y direction, z - new z direction
        */
        void calculate(Matrix3d &m, const std::vector<Vector3d>& r, bool &sameChirality) const;
		
    private:
        mutable bool mIsCartesian;
        
        std::shared_ptr<CrossProductLcs> mCrossProductLcs;
        LocalCoordinateSystem<int> mLcs;
        int mThirdCoordinate;
        void calcAtomicPositionsBasedDirection(const std::vector<int> &refPoint, const int &centralAtom,
            const std::vector<Vector3d>& r, LcsDirectionType type, Vector3d &direction) const;
        
        static void calculateAnyOrthogonal(const Vector3d &r0, Vector3d &r);
        void takeChiralityIntoAccount(Vector3d &x, Vector3d &y, Vector3d &z, const std::vector<Vector3d>& r, bool &sameChirality) const;
	};

    /**@}*/
}

