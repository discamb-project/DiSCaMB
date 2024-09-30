#pragma once

#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"

namespace discamb {

    /**
    * \addtogroup CrystalStructure
    * @{
    */


	class ConstantLocalCoordinateSystem : public LocalCoordinateSystemInCrystal
	{
		Vector3d mX, mY, mZ;
	public:
		ConstantLocalCoordinateSystem();
		ConstantLocalCoordinateSystem(const Vector3d &x, const Vector3d &y, const Vector3d &z);
		virtual ~ConstantLocalCoordinateSystem();
		// [1,1,0],[-0.5,0.5,0],[0,0,0.3] - xyz are normalized
		virtual void set(const std::string &definition, const Crystal &c);
		void set(const Vector3d &x, const Vector3d &y, const Vector3d &z);
		virtual void calculate(Vector3d &x, Vector3d &y, Vector3d &z, const Crystal &c) const;

	};
    /**@}*/
}