#include "discamb/Scattering/HcFormFactorCalculationsManager.h"

using namespace std;

namespace {

    struct SfCollector {

        // [hkl][atom]
        complex<double> atomicFF;
        int atomIndex;

        void onPerHklCalculation(
            int hkl_idx,
            std::complex<double> &structureFactor,
            discamb::SfDerivativesAtHkl &_derivatives)
        {
            atomicFF = structureFactor;
        }

    };

}

namespace discamb {

//    HcFormaFactorCalculationsManager::HcFormaFactorCalculationsManager() {}

    HcFormFactorCalculationsManager::HcFormFactorCalculationsManager(
        const Crystal &crystal,
        const HC_ModelParameters &params,
		const std::vector < std::shared_ptr <LocalCoordinateSystemInCrystal> > &lcs,
        bool frozen_lcs
		//bool newImplementation
        //const std::vector<XdLocalCoordinateSystem> &lcs
	): mLcs(lcs)
    {
        mFrozenLcs = frozen_lcs;
		//mNewImplementation = newImplementation;
        int i, nAtoms = crystal.atoms.size();
        mLcsMatrices.resize(nAtoms);
        mReciprocalLatticeUnitCell.set(crystal.unitCell);
        
        mCrystal = crystal;      


        for (i = 0; i < nAtoms; i++)
        {
            //mLcs[i].calculate(mLcsMatrices[i], crystal);
			mLcs[i]->calculate(mLcsMatrices[i], crystal);
            mCrystal.atoms[i].adp.clear();
            mCrystal.atoms[i].coordinates = Vector3d(0, 0, 0);
            mCrystal.atoms[i].multiplicity = 1;
        }

        mCrystal.spaceGroup.set({ string("x,y,z") });

        mHcCalculator = new HansenCoppensStructureFactorCalculator(mCrystal, params);

        mAuxCrystal = crystal;
    }

    HcFormFactorCalculationsManager::~HcFormFactorCalculationsManager()
    {
        delete mHcCalculator;
    }

	//void HcFormFactorCalculationsManager::useNewImplementation(
	//	bool newImplementation)
	//{
	//	//mNewImplementation = newImplementation;
	//}

    void HcFormFactorCalculationsManager::update(
        const std::vector<AtomInCrystal> &atoms)
    {
        int i, nAtoms = atoms.size();

        mAuxCrystal.atoms = atoms;
        if(!mFrozenLcs)
            for (i = 0; i < nAtoms; i++)
                mLcs[i]->calculate(mLcsMatrices[i], mAuxCrystal);
    }

    std::complex<double> HcFormFactorCalculationsManager::calculateFrac(
        int atomIdx, 
        const Vector3i &hkl)
    const 
    {
        SfCollector collector;
        vector<complex<double> > f(1);
        vector<bool> countAtomContribution(mCrystal.atoms.size(), false);
        countAtomContribution[atomIdx] = true;
        mHcCalculator->calculateStructureFactors(mCrystal.atoms, mLcsMatrices, { hkl }, f, countAtomContribution);
        return f[0]; 
    }

    std::complex<double> HcFormFactorCalculationsManager::calculateCart(
        int atomIdx, 
        const Vector3d &hkl) 
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
 
	void HcFormFactorCalculationsManager::calculateFrac(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
		/*if (!mNewImplementation)
		{
			AtomicFormFactorCalculationsManager::calculateCart(hkl, formFactors, includeAtom);
			return;
		}*/

		mHcCalculator->calculateFormFactors(mLcsMatrices, hkl, formFactors, includeAtom);
	}

	void HcFormFactorCalculationsManager::calculateCart(
		const Vector3d& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
        mHcCalculator->calculateFormFactorsCart(mLcsMatrices, hkl, formFactors, includeAtom);
		/*if (!mNewImplementation)
		{
			AtomicFormFactorCalculationsManager::calculateCart(hkl, formFactors, includeAtom);
			return;
		}*/
        /*
		Vector3d frac;
		Vector3i hklFrac;
		mReciprocalLatticeUnitCell.cartesianToFractional(hkl, frac);

		hklFrac[0] = math_utilities::roundInt(frac[0]);
		hklFrac[1] = math_utilities::roundInt(frac[1]);
		hklFrac[2] = math_utilities::roundInt(frac[2]);

		calculateFrac(hklFrac, formFactors, includeAtom);
        */
	}

}
