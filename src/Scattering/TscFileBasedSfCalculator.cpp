#include "discamb/Scattering/TscFileBasedSfCalculator.h"

#include <fstream>
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/Timer.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"


using namespace std;

namespace discamb {
	

	TscFileBasedSfCalculator::TscFileBasedSfCalculator(
		const Crystal& crystal, 
		const nlohmann::json& data)
	{
		string tscFileName;
		if (data.find("tsc file") != data.end())
			tscFileName = data.find("tsc file")->get<string>();
        else
        {
            tscFileName = file_system_utilities::find_newest_file("tsc", false);
            if(tscFileName.empty())
                on_error::throwException("olex tsc file name not given and no file with tsc extension exists in current directory", __FILE__, __LINE__);
        }

        //cout << "reads tsc file " << tscFileName;
        WallClockTimer timer;
        timer.start();

		ifstream in(tscFileName);

		if (!in.good())
			on_error::throwException("can not read olex tsc file", __FILE__, __LINE__);

		bool dataSection = false;
		string line;
		int nAtoms = crystal.atoms.size();
		int expectedN_columns = nAtoms + 3;
		map<Vector3i, vector<complex<double> > > formFactors;
		Vector3i hkl;
		vector<complex<double> > ff(nAtoms);
		vector<string> words, ff_words;
		while (in.good() && !dataSection)
		{
			getline(in, line);
			string_utilities::split(line, words);
			if (words[0] == string("DATA:"))
				dataSection = true;
		}
		if (dataSection)
		{
			//getline(in, line);
			while (in.good())
			{
				getline(in, line);
				string_utilities::split(line, words);
                
				if (words.size() == expectedN_columns)
				{
					hkl = { stoi(words[0]), stoi(words[1]), stoi(words[2]) };
                    
					for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
					{
						string_utilities::split(words[atomIdx + 3], ff_words, ',');
						ff[atomIdx] = { stod(ff_words[0]), stod(ff_words[1]) };
					}
					formFactors[hkl] = ff;
				}
			}
		}

        //cout << " " << timer.stop() << " ms" << endl;
		mCalculator = new AnyScattererStructureFactorCalculator(crystal);
		shared_ptr<AtomicFormFactorCalculationsManager> manager = shared_ptr<AtomicFormFactorCalculationsManager>(
			new ConstFormFactorCalculationsManager(crystal.unitCell, formFactors));

        bool convertXrayToElectrons = data.value<bool>("x2e", false);

        if (convertXrayToElectrons)
        {
            cout<< "converts X-ray to electron form factors ";
            timer.start();

            vector<double> nuclearCharges(crystal.atoms.size());
            for (int i = 0; i < nuclearCharges.size(); i++)
                nuclearCharges[i] = basic_chemistry_utilities::atomicNumberFromLabel(crystal.atoms[i].type);


            std::shared_ptr<AtomicFormFactorCalculationsManager> _manager =
                std::shared_ptr<AtomicFormFactorCalculationsManager>(
                    new ElectronFromXrayFormFactorCalculationManager(crystal.unitCell, nuclearCharges, manager));
            mCalculator->setAtomicFormfactorManager(_manager);
            cout << " " << timer.stop() << " ms" << endl;
            
        }
        else
		    mCalculator->setAtomicFormfactorManager(manager);
	}

	

	void TscFileBasedSfCalculator::setAnomalous(
		const std::vector<std::complex<double> >& anomalous) 
	{
		mCalculator->setAnoumalous(anomalous);
	}

	void TscFileBasedSfCalculator::calculateStructureFactorsAndDerivatives(
		const std::vector<AtomInCrystal>& atoms,
		const std::vector<Vector3i>& hkl,
		std::vector<std::complex<double> >& f,
		std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
		const std::vector<std::complex<double> >& dTarget_df,
		const std::vector<bool>& countAtomContribution)
	{
		mCalculator->update(atoms);
		mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
	}

	void TscFileBasedSfCalculator::calculateStructureFactors(
		const std::vector<AtomInCrystal>& atoms,
		const std::vector<Vector3i>& hkl,
		std::vector<std::complex<double> >& f,
		const std::vector<bool>& countAtomContribution) 
	{
		update(atoms);
		mCalculator->calculateStructureFactors(hkl, f, countAtomContribution);

	}

	void TscFileBasedSfCalculator::update(
		const std::vector<AtomInCrystal>& atoms) 
	{
		mCalculator->update(atoms);
	}

	void TscFileBasedSfCalculator::calculateStructureFactorsAndDerivatives(
		const Vector3i& hkl,
		std::complex<double>& scatteringFactor,
		discamb::SfDerivativesAtHkl& derivatives,
		const std::vector<bool>& countAtomContribution) 
	{
		mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives, countAtomContribution);
	}

	void TscFileBasedSfCalculator::calculateFormFactors(
		const Vector3i& hkl, 
		std::vector<std::complex<double> >& formFactors, 
		const std::vector<bool>& includeAtom) 
		const
	{
		mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
	}

	TscFileBasedSfCalculator::~TscFileBasedSfCalculator(){
		delete mCalculator;
	}
	
}