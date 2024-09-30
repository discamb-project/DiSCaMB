#include "discamb/CrystalStructure/AtomInCrystalID.h"
#include "discamb/BasicUtilities/on_error.h"

#include <algorithm>

using namespace std;

namespace discamb {
	

	AtomInCrystalID::AtomInCrystalID()
	{}

	AtomInCrystalID::AtomInCrystalID(
		int indexInAsu, 
		const SpaceGroupOperation &symmetryOperation) 
	{
		set(indexInAsu, symmetryOperation);
	}

	AtomInCrystalID::AtomInCrystalID(
		const std::string &label,
		const Crystal &c,
		const SpaceGroupOperation &symmetryOperation)
	{
		set(label, c, symmetryOperation);
	}



	AtomInCrystalID::~AtomInCrystalID()
	{
	}

	void AtomInCrystalID::set(
		int indexInAsu, 
		const SpaceGroupOperation &symmetryOperation) 
	{
		mIndex = indexInAsu;
		mSpaceGroupOperation = symmetryOperation;
	}

	void AtomInCrystalID::set(
		const std::string &label,
		const Crystal &c,
		const SpaceGroupOperation &symmetryOperation)
	{
		if (uniqueLabel(label, c))
		{
			mIndex = findIndex(label, c);
			mSpaceGroupOperation = symmetryOperation;
		}
		else
			mIndex = -1;
	}

	
	bool AtomInCrystalID::uniqueLabel(
		const std::string &label, 
		const Crystal &c)
	{
		int count = 0;
		for (auto atom : c.atoms)
			if (atom.label == label)
				count++;
		return (count == 1);
	}

	int AtomInCrystalID::findIndex(
		const std::string &label,
		const Crystal &c)
	{
		for (int i = 0; i < c.atoms.size(); i++)
			if (c.atoms[i].label == label)
				return int(i);
		return -1;
	}

	int AtomInCrystalID::index()
		const 
	{ 
		return mIndex;
	}

	const SpaceGroupOperation &AtomInCrystalID::getSymmetryOperation()
		const
	{
		return mSpaceGroupOperation;
	}

	
}