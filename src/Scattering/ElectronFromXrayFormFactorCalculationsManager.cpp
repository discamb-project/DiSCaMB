#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/BasicChemistry/PeriodicTable.h"

using namespace std;

namespace discamb {

    ElectronFromXrayFormFactorCalculationManager::ElectronFromXrayFormFactorCalculationManager(
        const UnitCell &unitCell,
        const std::vector<double> &nuclearCharge,
		std::shared_ptr<AtomicFormFactorCalculationsManager> &manager)
    {
        mManager = manager;
        mNuclearCharge = nuclearCharge;
        mRUnitCell.set(unitCell);
		int nAtoms = nuclearCharge.size();
		mFormFactor_x.resize(nAtoms);
        setFormFactorsAtHkl000();
    }

    ElectronFromXrayFormFactorCalculationManager::~ElectronFromXrayFormFactorCalculationManager()
    {
    }

    
    void ElectronFromXrayFormFactorCalculationManager::setFormFactorsAtHkl000()
    {
        //std::vector<std::complex<double> > mFormFactorsAtHkl000;
        int atomIdx, nAtoms = mNuclearCharge.size();
        Crystal c;
        c.atoms.resize(nAtoms);
        vector<int> z;
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            string symbol, label;
            symbol = periodic_table::symbol(lrint(mNuclearCharge[atomIdx]));
            label = symbol + to_string(atomIdx + 1);
            c.atoms[atomIdx].type = symbol;
            c.atoms[atomIdx].label = label;
        }
        IamFormFactorCalculationsManager m(c, "electron-IT");
        mFormFactorsAtHkl000.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            mFormFactorsAtHkl000[atomIdx] = m.calculateFrac(atomIdx, { 0,0,0 });
    }


    void ElectronFromXrayFormFactorCalculationManager::update(
        const std::vector<AtomInCrystal> &atoms)
    {
        mManager->update(atoms);
    }

    std::complex<double> ElectronFromXrayFormFactorCalculationManager::calculateFrac(
        int atomIdx, 
        const Vector3i &hkl)
        const
    {
        if (hkl == Vector3i(0, 0, 0))
            return mFormFactorsAtHkl000[atomIdx];

        Vector3d hCart;
        mRUnitCell.fractionalToCartesian(hkl, hCart);

        return calculateCart(atomIdx, hCart);
    }

    std::complex<double> ElectronFromXrayFormFactorCalculationManager::calculateCart(
        int atomIdx,
        const Vector3d &hklCart)
        const
    {
        if (hklCart == Vector3d(0, 0, 0))
            return mFormFactorsAtHkl000[atomIdx];


        double h = sqrt(hklCart*hklCart);
        complex<double> fx = mManager->calculateCart(atomIdx, hklCart);

		complex<double> nuclear = 0.023934 * mNuclearCharge[atomIdx] / (0.25 * h * h);
		complex<double> electron = -0.023934 * fx / (0.25 * h * h);
        return nuclear + electron;
        //return 0.023934 * (mNuclearCharge[atomIdx] - fx) / (0.25 * h * h);
    }

	void ElectronFromXrayFormFactorCalculationManager::calculateFrac(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom) 
		const
	{
        if (hkl == Vector3i(0, 0, 0))
        {
            formFactors = mFormFactorsAtHkl000;
            return;
        }


		Vector3d hCart;
		mRUnitCell.fractionalToCartesian(hkl, hCart);
		//calculateCart(hkl, formFactors, includeAtom);
        calculateCart(hCart, formFactors, includeAtom);
	}

	void ElectronFromXrayFormFactorCalculationManager::calculateCart(
		const Vector3d& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
        if (hkl == Vector3d(0, 0, 0))
        {
            formFactors = mFormFactorsAtHkl000;
            return;
        }


		double h2 = hkl * hkl;
		mManager->calculateCart(hkl, mFormFactor_x , includeAtom);

		//complex<double> fx = mManager->calculateCart(atomIdx, hklCart);

		int atomIdx, nAtoms = includeAtom.size();
		
		if (formFactors.size() != nAtoms)
			formFactors.resize(nAtoms);

		complex<double> nuclear_ff, ed_ff;

		for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
		{
			nuclear_ff = 0.023934 * mNuclearCharge[atomIdx] / (0.25 * h2);
			ed_ff = -0.023934 * mFormFactor_x[atomIdx] / (0.25 * h2);
			//complex<double> electron = -0.023934 * fx / (0.25 * h * h);
			//return nuclear + electron;
			formFactors[atomIdx] = nuclear_ff + ed_ff;
		}

		//complex<double> nuclear = 0.023934 * mNuclearCharge[atomIdx] / (0.25 * h * h);
		//complex<double> electron = -0.023934 * fx / (0.25 * h * h);
		//return nuclear + electron;
	}

    /**
    formFactors[hkl idx][atom idx]
    */
    void ElectronFromXrayFormFactorCalculationManager::calculateFrac(
        const std::vector<Vector3i>& hkl,
        vector < vector<complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        Vector3d hCart;
        int nHkl = hkl.size();
        mManager->calculateFrac(hkl, formFactors, includeAtom);

        int nAtoms;
        if (!formFactors.empty())
            nAtoms = formFactors[0].size();
        else
            return;

        complex<double> nuclear_ff, ed_ff;

        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            mRUnitCell.fractionalToCartesian(hkl[hklIdx], hCart);
            double h2 = hCart * hCart;
            
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                nuclear_ff = 0.023934 * mNuclearCharge[atomIdx] / (0.25 * h2);
                ed_ff = -0.023934 * formFactors[hklIdx][atomIdx] / (0.25 * h2);
                formFactors[hklIdx][atomIdx] = nuclear_ff + ed_ff;
            }
        }

    }

    /**
    formFactors[hkl idx][atom idx]
    */
    void ElectronFromXrayFormFactorCalculationManager::calculateCart(
        const std::vector <Vector3d>& hkl,
        std::vector < std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        
        int nHkl = hkl.size();
        mManager->calculateCart(hkl, formFactors, includeAtom);

        int nAtoms;
        if (!formFactors.empty())
            nAtoms = formFactors[0].size();
        else
            return;

        complex<double> nuclear_ff, ed_ff;

        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            
            double h2 = hkl[hklIdx] * hkl[hklIdx];

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                nuclear_ff = 0.023934 * mNuclearCharge[atomIdx] / (0.25 * h2);
                ed_ff = -0.023934 * formFactors[hklIdx][atomIdx] / (0.25 * h2);
                formFactors[hklIdx][atomIdx] = nuclear_ff + ed_ff;
            }
        }

    }


    /*
    private:
        ElectronFromXrayFormFactorCalculationManager();
        SharedPointer<AtomicFormFactorCalculationsManager>::type &mManager;

    };*/

}
