#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/MathUtilities/math_utilities.h"

using namespace std;

namespace discamb {

	ConstFormFactorCalculationsManager::ConstFormFactorCalculationsManager(
		const UnitCell& uc,
		std::map<Vector3i, std::vector<std::complex<double> > >& formFactors) 
	{
		mFormFactors = formFactors;
		mReciprocalLatticeUnitCell.set(uc);
	}

	ConstFormFactorCalculationsManager::~ConstFormFactorCalculationsManager()
	{
		int i = 1;
	}

	void ConstFormFactorCalculationsManager::update(
		const std::vector<AtomInCrystal>& atoms) 
	{
	}


	std::complex<double> ConstFormFactorCalculationsManager::calculateFrac(
		int atomIdx, 
		const Vector3i& hkl) 
		const 
	{
		auto it = mFormFactors.find(hkl);
		if (it == mFormFactors.end())
		{
			string hklString = string("[") + to_string(hkl[0]) + string(", ") + to_string(hkl[1]) + string(", ") + to_string(hkl[2]) + string("]");
			on_error::throwException(string("asking for form factors for hkl with no data, hkl = ") + hklString, __FILE__, __LINE__);
			return { 0,0 };
		}
		return it->second[atomIdx];
	}

	std::complex<double> ConstFormFactorCalculationsManager::calculateCart(
		int atomIdx, 
		const Vector3d& hkl) 
		const 
	{
		Vector3d frac;
		Vector3i hklFrac;
		mReciprocalLatticeUnitCell.cartesianToFractional(hkl, frac);

		hklFrac[0] = math_utilities::roundInt(frac[0]);
		hklFrac[1] = math_utilities::roundInt(frac[1]);
		hklFrac[2] = math_utilities::roundInt(frac[2]);

		return calculateFrac(atomIdx, hklFrac);

	}

	void ConstFormFactorCalculationsManager::calculateFrac(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom) 
		const
	{
		auto it = mFormFactors.find(hkl);
		if (it == mFormFactors.end())
		{
			string hklString = string("[") + to_string(hkl[0]) + string(", ") + to_string(hkl[1]) + string(", ") + to_string(hkl[2]) + string("]");
			on_error::throwException(string("asking for form factors for hkl with no data, hkl = ") + hklString, __FILE__, __LINE__);
			
		}
		int atomIdx, nAtoms = includeAtom.size();
		formFactors.assign(nAtoms, 0);
		for(atomIdx=0; atomIdx<nAtoms;atomIdx++)
			if(includeAtom[atomIdx])
				formFactors[atomIdx] = it->second[atomIdx];
	}

	void ConstFormFactorCalculationsManager::calculateCart(
		const Vector3d& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom) 
		const 
	{
		Vector3d frac;
		Vector3i hklFrac;
		mReciprocalLatticeUnitCell.cartesianToFractional(hkl, frac);

		hklFrac[0] = math_utilities::roundInt(frac[0]);
		hklFrac[1] = math_utilities::roundInt(frac[1]);
		hklFrac[2] = math_utilities::roundInt(frac[2]);

		calculateFrac(hklFrac, formFactors, includeAtom);

	}
}
