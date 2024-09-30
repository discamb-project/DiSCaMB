#pragma once

#include "discamb/CrystalStructure/Crystal.h"


#include <vector>
#include <optional>

namespace discamb {
    /**
    * \addtogroup CrystalStructure
    * @{
    */

	/**
	identifies atom in crystal c with its index in asymmetric unit (i.e. index of atom in c.atoms) and symmetry operation
	*/
	class AtomInCrystalID {
	public:
		AtomInCrystalID();
		AtomInCrystalID(int indexInAsu, const SpaceGroupOperation &symmetryOperation = SpaceGroupOperation());
		AtomInCrystalID(const std::string &label, const Crystal &c, const SpaceGroupOperation &symmetryOperation);
				
		~AtomInCrystalID();

		void set(int indexInAsu, const SpaceGroupOperation &symmetryOperation = SpaceGroupOperation());
		void set(const std::string &label, const Crystal &c, const SpaceGroupOperation &symmetryOperation);

		static bool uniqueLabel(const std::string &label, const Crystal &crystal);
		static int findIndex(const std::string &label, const Crystal &crystal);

        int index() const;
		const SpaceGroupOperation &getSymmetryOperation() const;
		
	private:
		SpaceGroupOperation mSpaceGroupOperation;
		int mIndex = 0;	
		
	};

    /**@}*/
}
