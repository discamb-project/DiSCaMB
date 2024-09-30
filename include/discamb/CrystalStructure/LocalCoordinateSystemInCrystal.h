
#ifndef _DISCAMB_CRYSTALSTRUCTURE_LOCALCOORDINATESYSTEMINCRYSTAL_HPP_
#define _DISCAMB_CRYSTALSTRUCTURE_LOCALCOORDINATESYSTEMINCRYSTAL_HPP_

#include "discamb/CrystalStructure/Crystal.h"

namespace discamb {

    /**
    * \addtogroup CrystalStructure
    * @{
    */


	class LocalCoordinateSystemInCrystal
	{
	public:
		/** \brief Constructs the coordinate system as the Cartesian one.*/
		//LocalCoordinateSystemInCrystal();
		virtual ~LocalCoordinateSystemInCrystal() = 0;

		/**
		\brief Define the local coordinate system.
		*/
		virtual void set(const std::string &definition, const Crystal &c) = 0;

		/** \brief Calculates coordinate system vectors for particular geometry of crystal structure.

		   The orthonormal vectors (\p x , \p y , \p z ) of coordinate system corresponds to geometry
		   of the structure provided in the argument \p c .*/
		virtual void calculate(Vector3d &x, Vector3d &y, Vector3d &z, const Crystal &c) const = 0;

		/** \brief Calculates coordinate system vectors for particular geometry of crystal structure.

			Columns of the matrix \p correspond to the new base vectors.
			Geometry of the structure provided in the argument \p c is used in the calculation.*/
		virtual void calculate(Matrix3d &m, const Crystal &c) const
		{
			Vector3d x, y, z;
			this->calculate(x, y, z, c);
			m.set(x[0], y[0], z[0],
				x[1], y[1], z[1],
				x[2], y[2], z[2]);
		}

	};

    /**@}*/
}



#endif
 