#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

class ChiralityComparator {
public:
    ChiralityComparator(
        AtomInCrystalID &cenraltAtom, 
        AtomInCrystalID &direction_1_atom, 
        AtomInCrystalID &direction_2_atom, 
        AtomInCrystalID &reference_atom);

	~ChiralityComparator();
	/**
	central_atom, atom_1, atom_2, atom_3
	when checking chirality	three vectors are created - from central atom to i-th atom:
	v1 = atom_1 - central_atom
	v2 = atom_2 - central_atom
	v3 = atom_3 - central_atom
	then cross product of the first two vectors is calculated and it is checked
	if the angle between the crossproduct and the third vector is below 90, if yes then true is returned
	*/
	void set(AtomInCrystalID &centralAtom, AtomInCrystalID &direction_1_atom, AtomInCrystalID &direction_2_atom, AtomInCrystalID &reference_atom);
	bool sameChilarity(const Crystal &);
private:
	AtomInCrystalID mCentralAtom, mDirection1Atom, mDirection2Atom, mReferenceAtom;
};
/**@}*/
}