#include "discamb/Scattering/AnyHcCalculator.h"
#include "discamb/IO/cif_io.h"
#include "discamb/IO/xd_io.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/BasicUtilities/Timer.h"

using namespace std;

namespace discamb {

    AnyHcCalculator::AnyHcCalculator(
        const Crystal &crystal,
        const HC_ModelParameters &parameters,
        const std::vector<std::shared_ptr<LocalCoordinateSystemInCrystal> > &lcs,
        bool electronScattering,
        bool spherical,
		bool implementationForLargeMolecules,
		int nThreads,
        bool frozenLcs)
    {
        set(crystal, parameters, lcs, electronScattering, spherical, implementationForLargeMolecules, nThreads, frozenLcs);
    }

    void AnyHcCalculator::getModelInformation(
        std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        modelInfo.clear();
        modelInfo.push_back({ "SCATTERING MODEL", "Hansen-Coppens multipole model"});
    }

    AnyHcCalculator::AnyHcCalculator(
        const Crystal &crystal,
        const nlohmann::json &data)
    {
        //const Crystal &crystal,
        //  const HC_ModelParameters &parameters,
        //const std::vector<XdLocalCoordinateSystem> &lcs,
        //bool electronScattering)

        bool electronScattering = false;
        bool hcIam = false;
		

        if (data.find("electron_scattering") != data.end())
            electronScattering = data.find("electron_scattering")->get<bool>();

        if (data.find("spherical neutral") != data.end())
            hcIam = data.find("spherical neutral")->get<bool>();

        string xd_inp("xd.inp"), xd_mas("xd.mas");

        if (data.find("xd mas") != data.end())
            xd_mas = data.find("xd mas")->get<string>();

        if (data.find("xd inp") != data.end())
            xd_inp = data.find("xd inp")->get<string>();

		bool useMacromolecularImplementation = false;
		if (data.find("use macromolecular implementation") != data.end())
			useMacromolecularImplementation = data.find("use macromolecular implementation")->get<bool>();

		int nThreads = 1;
		if (data.find("number of threads") != data.end())
			nThreads = data.find("number of threads")->get<int>();

        bool frozenLcs = data.value("frozen lcs", false);


        HC_ModelParameters hcModelParameters;
        vector<XdLocalCoordinateSystem> localCoordinateSystems;
		vector<shared_ptr<LocalCoordinateSystemInCrystal> > lcs;
        Crystal crystal0;

        xd_io::read(xd_mas, xd_inp, hcModelParameters, crystal0, localCoordinateSystems, true);
        
		for (auto &_lcs : localCoordinateSystems)
			lcs.push_back(shared_ptr<LocalCoordinateSystemInCrystal>(new XdLocalCoordinateSystem(_lcs)));
        
		set(crystal, hcModelParameters, lcs, electronScattering, hcIam, useMacromolecularImplementation, nThreads, frozenLcs);

        string multipoleCifFile = data.value("multipole cif", string());

        //if(!multipoleCifFile.empty())
          //  cif_io::saveCif(multipoleCifFile, crystal, hcModelParameters, lcs);
    }

    void AnyHcCalculator::setAnomalous(
        const std::vector<std::complex<double> > & anomalous)
    {
        mCalculator->setAnoumalous(anomalous);
    }

	void AnyHcCalculator::setN_threads(
		int n)
	{
		if (mUseImplementationForLargeMolecules)
			mHcCalculator->setN_Threads(int(n));
		mN_threads = n;
	}

    void AnyHcCalculator::set(
        const Crystal &crystal,
        const HC_ModelParameters &_parameters,
        const std::vector<std::shared_ptr<LocalCoordinateSystemInCrystal> > &lcs,
        bool electronScattering,
        bool integerChargeSpherical,
		bool implementationForLargeMolecules,
		int nThreads,
        bool frozenLcs)
    {
        //AnyScattererStructureFactorCalculator 
        mCrystal = crystal;
        
        HC_ModelParameters parameters = _parameters;
		mUseImplementationForLargeMolecules = implementationForLargeMolecules;

        if (integerChargeSpherical)
        {
            int nWfnTypes = parameters.wfn_parameters.size();
            vector<int> nCoreElectrons(nWfnTypes, 0);

            int i = 0;
            for (auto &wfnType : parameters.wfn_parameters)
            {
                for (auto index : wfnType.core_orbitals_indices)
                    nCoreElectrons[i] += wfnType.orbital_occupancy[index];

                i++;
            }

            vector<bool> typeProcessed(parameters.type_parameters.size(), false);
            int nAtoms = parameters.atom_to_type_map.size();
            for (i = 0; i < nAtoms; i++)
            {
                int typeIdx = parameters.atom_to_type_map[i];

                if (!typeProcessed[typeIdx])
                {
                    int wfnType = parameters.atom_to_wfn_map[i];
                    parameters.type_parameters[typeIdx].kappa_deformation_valence = 1.0;
                    parameters.type_parameters[typeIdx].kappa_spherical_valence = 1.0;
                    parameters.type_parameters[typeIdx].p_lm.clear();
                    parameters.type_parameters[typeIdx].p_lm.resize(1, vector<double>(1, 0));
                    parameters.type_parameters[typeIdx].p_val = -parameters.wfn_parameters[wfnType].charge +
                        int(parameters.wfn_parameters[wfnType].atomic_number) -
                        int(nCoreElectrons[wfnType]);
                    typeProcessed[typeIdx] = true;
                }

            }

        }

        if (electronScattering)
        {
            vector<double> nuclearCharge(crystal.atoms.size());
            for (int i = 0; i < nuclearCharge.size(); i++)
                nuclearCharge[i] = static_cast<double>(basic_chemistry_utilities::atomicNumberFromLabel(crystal.atoms[i].type));

			shared_ptr<AtomicFormFactorCalculationsManager> hcManager =
				shared_ptr<AtomicFormFactorCalculationsManager>(
                    new HcFormFactorCalculationsManager(crystal, parameters, lcs, frozenLcs));


            mManager = shared_ptr<AtomicFormFactorCalculationsManager>(
                new ElectronFromXrayFormFactorCalculationManager(
                    crystal.unitCell,
                    nuclearCharge,
                    hcManager));
            // SharedPointer<AtomicFormFactorCalculationsManager>::type(
            //   new HcFormFactorCalculationsManager(crystal, parameters, lcs))));
        }
        else
            mManager = shared_ptr<AtomicFormFactorCalculationsManager>(
                new HcFormFactorCalculationsManager(crystal, parameters, lcs, frozenLcs));
        mCalculator = new AnyScattererStructureFactorCalculator(crystal);
        mCalculator->setAtomicFormfactorManager(mManager);

		if (implementationForLargeMolecules)
		{
			mHcCalculator = new HansenCoppensStructureFactorCalculator(mCrystal, parameters);
            if (electronScattering)
                mHcCalculator->setElectronScattering(true);
			mLcs = lcs;
			int nLcs = mLcs.size();
			mLcsMatrices.resize(nLcs);
			for (int i = 0; i < nLcs; i++)
				mLcs[i]->calculate(mLcsMatrices[i], crystal);
			setN_threads(nThreads);
		}
    }


    AnyHcCalculator::~AnyHcCalculator()
    {
        delete mCalculator;
        if (mUseImplementationForLargeMolecules)
            delete mHcCalculator;
    }

    void AnyHcCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        const std::vector<std::complex<double> > &dTarget_df,
        const std::vector<bool> &countAtomContribution)
    {
        
		if (!mUseImplementationForLargeMolecules)
		{
			mCalculator->update(atoms);
			mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
		}
		else
		{
            mHcCalculator->calculateStructureFactorsAndDerivatives(
                atoms, mLcsMatrices, hkl,
                f,
                dTarget_dparam,
                dTarget_df,
                countAtomContribution);
		}
    }

    void AnyHcCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution,
        const DerivativesSelector& derivativesSelector)
    {

        if (!mUseImplementationForLargeMolecules)
        {
            mCalculator->update(atoms);
            mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution, derivativesSelector);
        }
        else
        {
            mHcCalculator->calculateStructureFactorsAndDerivatives(
                atoms, mLcsMatrices, hkl,
                f,
                dTarget_dparam,
                dTarget_df,
                countAtomContribution, derivativesSelector);
        }
    }


	void AnyHcCalculator::calculateFormFactors(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
		mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
	}


    void AnyHcCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        const std::vector<bool> &countAtomContribution)
    {
		update(atoms);

		WallClockTimer timer;

		timer.start();
		if (mUseImplementationForLargeMolecules)
			mHcCalculator->calculateStructureFactors(atoms, mLcsMatrices, hkl, f);
		else
			mCalculator->calculateStructureFactors(hkl, f, countAtomContribution);
		//cout << "time of structure factor calculation in HC model, excluding model setup, in ms " << timer.stop() << endl;
        //vector<complex<double> > fake_dTarget_df(hkl.size());
        //std::vector<TargetFunctionAtomicParamDerivatives> dTarget_dparam;
        //mCalculator->update(atoms);
        //mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, fake_dTarget_df, countAtomContribution);

    }

    void AnyHcCalculator::update(const std::vector<AtomInCrystal> &atoms)
    {
        mCalculator->update(atoms);
		if (mUseImplementationForLargeMolecules)
		{
			mCrystal.atoms = atoms;
			int nLcs = mLcs.size();
			for (int i = 0; i < nLcs; i++)
				mLcs[i]->calculate(mLcsMatrices[i], mCrystal);
		}
        //on_error::not_implemented(__FILE__, __LINE__);
    }

    void AnyHcCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i &hkl,
        std::complex<double> &scatteringFactor,
        discamb::SfDerivativesAtHkl &derivatives,
        const std::vector<bool> &countAtomContribution)
    {
        //on_error::not_implemented(__FILE__, __LINE__);
        mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives, countAtomContribution);
    }


}