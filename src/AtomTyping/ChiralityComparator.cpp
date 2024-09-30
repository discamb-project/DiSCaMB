#include "discamb/AtomTyping/ChiralityComparator.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb {

        ChiralityComparator::ChiralityComparator(
			AtomInCrystalID &centralAtom, 
			AtomInCrystalID &direction_1_atom, 
			AtomInCrystalID &direction_2_atom, 
			AtomInCrystalID &reference_atom)
		{
			set(centralAtom, direction_1_atom, direction_2_atom, reference_atom);
		}

        ChiralityComparator::~ChiralityComparator()
		{
		}

		/**
		central_atom, atom_1, atom_2, atom_3
		when checking chirality	three vectors are created - from central atom to i-th atom:
		v1 = atom_1 - central_atom
		v2 = atom_2 - central_atom
		v3 = atom_3 - central_atom
		then cross product of the first two vectors is calculated and it is checked
		if the angle between the crossproduct and the third vector is below 90, if yes then true is returned
		*/
		void ChiralityComparator::set(
			AtomInCrystalID &cenraltAtom,
			AtomInCrystalID &direction_1_atom,
			AtomInCrystalID &direction_2_atom,
			AtomInCrystalID &reference_atom)
		{
			mCentralAtom = cenraltAtom;
			mDirection1Atom = direction_1_atom;
			mDirection2Atom = direction_2_atom;
			mReferenceAtom = reference_atom;
		}

		bool ChiralityComparator::sameChilarity(
			const Crystal &crystal)
		{
			on_error::not_implemented(__FILE__, __LINE__);
			return true;
		}


}
