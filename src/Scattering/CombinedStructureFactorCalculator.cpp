#include "discamb/Scattering/CombinedStructureFactorCalculator.h"

#include <memory>

using namespace std;

namespace discamb {
    CombinedStructureFactorCalculator::CombinedStructureFactorCalculator(
        std::vector<std::shared_ptr<SfCalculator> > &calculators,
        std::vector<std::vector<int> > &atoms)
    {
        set(calculators, atoms);
    }

    CombinedStructureFactorCalculator::CombinedStructureFactorCalculator(const Crystal &crystal, const nlohmann::json &data) {}
    CombinedStructureFactorCalculator::~CombinedStructureFactorCalculator() {};

    void CombinedStructureFactorCalculator::set(
        std::vector<std::shared_ptr<SfCalculator> > &calculators,
        std::vector<std::vector<int> > &atoms)
    {
        string errorMessage;

        if (calculators.empty())
        {
            errorMessage =
                "combined structure factor calculator should be composed of at least one structure factor calculator";
            on_error::throwException(errorMessage, __FILE__, __LINE__);
        }


        mCalculatorsAtoms = atoms;
        mCalculators = calculators;
        int nCalc, calcIdx, nAtoms;
        nCalc = mCalculators.size();
        mPartial_dTarget_dparam.resize(nCalc);
        mPartialDerivatives.resize(nCalc);
        mPartialScatteringFactorOneHkl.resize(nCalc); 
        mPartialScatteringFactorMultiHkl.resize(nCalc);

        nAtoms = 0;
        for (auto const &atom_list : atoms)
            nAtoms += atom_list.size();
        
        setAnomalous(vector<complex<double> >(nAtoms, 0));

        mCalculatorsAtomsInclusion.resize(nCalc,vector<bool>(nAtoms,false));
        for (calcIdx = 0; calcIdx < nCalc; calcIdx++)
            for (int atomIdx : atoms[calcIdx])
                mCalculatorsAtomsInclusion[calcIdx][atomIdx] = true;
    }


    void CombinedStructureFactorCalculator::setAnomalous(
        const std::vector<std::complex<double> > & anomalous) 
    {
        mAnomalous = anomalous;
        for (auto calculator : mCalculators)
            calculator->setAnomalous(mAnomalous);
    }

	void CombinedStructureFactorCalculator::calculateFormFactors(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
		int calculatorIdx, nCalculators = mCalculators.size();
		int atomIdx, nAtoms = includeAtom.size();
		vector<bool> atomContribution;
//		SfDerivativesAtHkl partialDerivatives;
		std::vector<std::complex<double> > calculatorFormFactors;
		formFactors.assign(nAtoms, 0.0);
		atomContribution = includeAtom;

		for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
		{
			for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
				atomContribution[atomIdx] = includeAtom[atomIdx] && mCalculatorsAtomsInclusion[calculatorIdx][atomIdx];
			calculatorFormFactors.assign(nAtoms, 0.0);

			mCalculators[calculatorIdx]->calculateFormFactors(
					hkl,
				    calculatorFormFactors,
					atomContribution);

			for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
				if (mCalculatorsAtomsInclusion[calculatorIdx][atomIdx])
					formFactors[atomIdx] = calculatorFormFactors[atomIdx];
		}

	}


    void CombinedStructureFactorCalculator::calculateFormFactorsCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        int calculatorIdx, nCalculators = mCalculators.size();
        int atomIdx, nAtoms = includeAtom.size();
        vector<bool> atomContribution;
        //		SfDerivativesAtHkl partialDerivatives;
        std::vector<std::complex<double> > calculatorFormFactors;
        formFactors.assign(nAtoms, 0.0);
        atomContribution = includeAtom;

        for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
        {
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                atomContribution[atomIdx] = includeAtom[atomIdx] && mCalculatorsAtomsInclusion[calculatorIdx][atomIdx];
            calculatorFormFactors.assign(nAtoms, 0.0);

            mCalculators[calculatorIdx]->calculateFormFactorsCart(
                hkl,
                calculatorFormFactors,
                atomContribution);

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (mCalculatorsAtomsInclusion[calculatorIdx][atomIdx])
                    formFactors[atomIdx] = calculatorFormFactors[atomIdx];
        }

    }

    void CombinedStructureFactorCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        const std::vector<std::complex<double> > &dTarget_df,
        const std::vector<bool> &countAtomContribution) 
    {
        int calculatorIdx, nCalculators = mCalculators.size();
        int atomIdx, nAtoms = countAtomContribution.size();
        vector<bool> atomContribution;

        atomContribution = countAtomContribution;

        for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
        {
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                atomContribution[atomIdx] = countAtomContribution[atomIdx] && mCalculatorsAtomsInclusion[calculatorIdx][atomIdx];

            mCalculators[calculatorIdx]->
                calculateStructureFactorsAndDerivatives(
                    atoms,
                    hkl,
                    mPartialScatteringFactorMultiHkl[calculatorIdx],
                    mPartial_dTarget_dparam[calculatorIdx],
                    dTarget_df,
                    atomContribution);
        }

        // merge results for scattering factors
        
        int hklIdx, nHkl = hkl.size();
        
        f = mPartialScatteringFactorMultiHkl[0];
        
        for (calculatorIdx = 1; calculatorIdx < nCalculators; calculatorIdx++)
            for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
                f[hklIdx] += mPartialScatteringFactorMultiHkl[calculatorIdx][hklIdx];

        // merge results for scattering factors derivatives

        dTarget_dparam.resize(atoms.size());
        for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
            for (int atomIdx : mCalculatorsAtoms[calculatorIdx])
                dTarget_dparam[atomIdx] = mPartial_dTarget_dparam[calculatorIdx][atomIdx];
    }

    void CombinedStructureFactorCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        const std::vector<bool> &countAtomContribution) 
    {
        int calculatorIdx, nCalculators = mCalculators.size();
        int atomIdx, nAtoms = countAtomContribution.size();
        vector<bool> atomContribution;

        atomContribution = countAtomContribution;

        for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
        {

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                atomContribution[atomIdx] = countAtomContribution[atomIdx] && mCalculatorsAtomsInclusion[calculatorIdx][atomIdx];


            mCalculators[calculatorIdx]->
                calculateStructureFactors(atoms, hkl, mPartialScatteringFactorMultiHkl[calculatorIdx], atomContribution);
        }

        // merge scattering factors

        f = mPartialScatteringFactorMultiHkl[0];
        int hklIdx, nHkl = f.size();
        for (calculatorIdx = 1; calculatorIdx < nCalculators; calculatorIdx++)
            for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
                f[hklIdx] += mPartialScatteringFactorMultiHkl[calculatorIdx][hklIdx];

    }

    void CombinedStructureFactorCalculator::update(
        const std::vector<AtomInCrystal> &atoms) 
    {
        int calculatorIdx, nCalculators = mCalculators.size();

        for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
            mCalculators[calculatorIdx]->update(atoms);

    }

    void CombinedStructureFactorCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i &hkl,
        std::complex<double> &scatteringFactor,
        SfDerivativesAtHkl &derivatives,
        const std::vector<bool> &countAtomContribution)
    {
        int calculatorIdx, nCalculators = mCalculators.size();
        int atomIdx, nAtoms = countAtomContribution.size();
        vector<bool> atomContribution;
        SfDerivativesAtHkl partialDerivatives;

        atomContribution = countAtomContribution;

        for (calculatorIdx = 0; calculatorIdx < nCalculators; calculatorIdx++)
        {
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                atomContribution[atomIdx] = countAtomContribution[atomIdx] && mCalculatorsAtomsInclusion[calculatorIdx][atomIdx];

            mCalculators[calculatorIdx]->
                calculateStructureFactorsAndDerivatives(
                    hkl,
                    mPartialScatteringFactorOneHkl[calculatorIdx],
                    partialDerivatives,
                    atomContribution);
            if (calculatorIdx == 0)
                derivatives = partialDerivatives;
            else
                for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    int adpComponentIdx, nAdp = partialDerivatives.adpDerivatives[atomIdx].size();
                    for (adpComponentIdx = 0; adpComponentIdx < nAdp; adpComponentIdx++)
                        derivatives.adpDerivatives[atomIdx][adpComponentIdx] += 
                            partialDerivatives.adpDerivatives[atomIdx][adpComponentIdx];
                    for (int i = 0; i < 3; i++)
                        derivatives.atomicPostionDerivatives[atomIdx][i] += partialDerivatives.atomicPostionDerivatives[atomIdx][i];
                    derivatives.occupancyDerivatives[atomIdx] += partialDerivatives.occupancyDerivatives[atomIdx];
                }
        }
        scatteringFactor = mPartialScatteringFactorOneHkl[0];
        for (calculatorIdx = 1; calculatorIdx < nCalculators; calculatorIdx++)
            scatteringFactor += mPartialScatteringFactorOneHkl[calculatorIdx];
    }

}

