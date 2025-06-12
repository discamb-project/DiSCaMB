#include "discamb/Scattering/AnyIamCalculator.h"

#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"



using namespace std;

namespace discamb {

    AnyIamCalculator::AnyIamCalculator(
        const Crystal &crystal,
        bool electronScattering,
        const std::string& table)
    {
        set(crystal, electronScattering, table, "standard");
    }

    AnyIamCalculator::AnyIamCalculator(
        const Crystal &crystal,
        const nlohmann::json &data)
    {
        bool electronScattering = false;
        string table;

        if (data.find("electron_scattering") != data.end())
            electronScattering = data.find("electron_scattering")->get<bool>();

        if (data.find("electron scattering") != data.end())
            electronScattering = data.find("electron scattering")->get<bool>();


        if (data.find("table") != data.end())
            table = data.find("table")->get<string>();

        string algorithm = data.value("algorithm", "standard");

        set(crystal, electronScattering, table, algorithm);
    }

    void AnyIamCalculator::setAnomalous(
        const std::vector<std::complex<double> > & anomalous)
    {
        mCalculator->setAnoumalous(anomalous);
    }

    void AnyIamCalculator::set(
        const Crystal &crystal,
        bool electronScattering,
        const string &table,
        const std::string& algorithm)
    {
        mModelInfo.clear();
        mModelInfo.push_back({ "SCATTERING MODEL", "IAM" });
      /*Waasmeier_Kirfel
        International Tables for Crystallography
            Volume C
            Mathematical, Physicaland Chemical Tables
            Edited by A.J.C.Wilson
            Kluwer Academic Publishers
            Dordrecht / Boston / London
            1992
            //Mott–Bethe
             for W-K
                H entry taken from cctbx r23554 (5 Gaussian fit to SDS)
                using the following code:
                    import cctbx.eltbx.xray_scattering
                    x = cctbx.eltbx.xray_scattering.wk1995("H").fetch()
                    print x.array_of_a(), x.array_of_b(), x.c()
            
            Table 4.3.2.3 of International Tables for Crystallography Vol.C(Cowley et al., 2006[Cowley, J.M., Peng, L.M., Ren, G., Dudarev, S.L. & Whelan, M.J. (2006).International Tables for Crystallography, Vol.C, edited by E.Prince, ch. 4.3.2.Chester:IUCr.]
                Tables 4.2.6.8 and 6.1.1.4 in of International Tables for Crystallography Vol.C(Wilson & Geist, 1993[Wilson, A.J.C. & Geist, V. (1993).Cryst.Res.Technol. 28, 110.])
                */
        string tableName;
        mUseLineAlgorithm = (algorithm == "macromol");

        mCrystal = crystal;
        string message = string(", setting IAM sf calculator with ");
        if (electronScattering)
        {
            if (table == string("") || table == string("electron-IT"))
            {
                mManager = shared_ptr<AtomicFormFactorCalculationsManager>(
                    new IamFormFactorCalculationsManager(crystal, "electron-IT"));
                message += string("electron-IT table");
                tableName = "Table 4.3.2.3 in Cowley, J.M., Peng, L.M., Ren, G., Dudarev, S.L. & Whelan, M.J. (2006).International Tables for Crystallography, Vol.C, edited by E.Prince, ch. 4.3.2.Chester:IUCr.";
            }
            else
            {
                vector<double> nuclearCharge(crystal.atoms.size());
                vector<int> atomicNumbers;
                crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);
                for (int i = 0; i < nuclearCharge.size(); i++)
                    nuclearCharge[i] = double(atomicNumbers[i]);

                shared_ptr<AtomicFormFactorCalculationsManager> xIamManager =
                    shared_ptr<AtomicFormFactorCalculationsManager>(
                        new IamFormFactorCalculationsManager(crystal, table));


                mManager = shared_ptr<AtomicFormFactorCalculationsManager>(
                    new ElectronFromXrayFormFactorCalculationManager(
                        crystal.unitCell,
                        nuclearCharge,
                        xIamManager));
                
                if ("Waasmeier-Kirfel" == table)
                    tableName = "Transformed with Mott–Bethe from X-ray for factors based on parameters from Waasmaier, D.& Kirfel, A. (1995).Acta Cryst.A51, 416 - 431. SDS model for hydrogen atom.";
                else
                    tableName = "Transformed with Mott–Bethe from X-ray for factors based on parameters from Tables 4.2.6.8 and 6.1.1.4 in of International Tables for Crystallography Vol.C (Wilson, A.J.C. & Geist, V. (1993).Cryst.Res.Technol. 28, 110.)";
            }
        }
        else
        {
            if (table == string("") || table == string("IT92"))
            {
                mManager = shared_ptr<AtomicFormFactorCalculationsManager>(
                    new IamFormFactorCalculationsManager(crystal, "IT92"));
                message += string("IT92 table");
                tableName = "Tables 4.2.6.8 and 6.1.1.4 in of International Tables for Crystallography Vol.C (Wilson, A.J.C. & Geist, V. (1993).Cryst.Res.Technol. 28, 110.)";
            }
            else
            {
                tableName = "Waasmaier, D.& Kirfel, A. (1995).Acta Cryst.A51, 416 - 431. SDS model for hydrogen atom.";
                mManager = shared_ptr<AtomicFormFactorCalculationsManager>(
                    new IamFormFactorCalculationsManager(crystal, table));
                message += string(table + string(" table"));
            }
        }

        add_to_log(__LINE__, string(__FILE__) + message);

        mCalculator = new AnyScattererStructureFactorCalculator(crystal);
        mCalculator2 = new AnyScattererStructureFactorCalculator2(crystal);
        /*
std::shared_ptr < AnyScattererStructureFactorCalculator2> mCalculator2;
        std::shared_ptr < AnyScattererStructureFactorCalculator> mCalculator;*/
        //mCalculator = new AnyScattererStructureFactorCalculator(crystal);
        //mCalculator = make_shared<AnyScattererStructureFactorCalculator>(crystal);
        //mCalculator2 = new AnyScattererStructureFactorCalculator2(crystal);
        //mCalculator2 = make_shared<AnyScattererStructureFactorCalculator2>(crystal);
        mCalculator->setAtomicFormfactorManager(mManager);
        mCalculator2->setAtomicFormfactorManager(mManager);

        mModelInfo.push_back({ "FORM FACTORS PARAMETERISATION", tableName });

        //if (!electronScattering)
        //    mModelInfo.push_back({ "RADIATION", "X-ray" });
        //else
        //    mModelInfo.push_back({ "RADIATION", "electron" });

    }

    AnyIamCalculator::~AnyIamCalculator() {
        //delete mCalculator;
        delete mCalculator2;
    }

    void AnyIamCalculator::getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        modelInfo = mModelInfo;
    }

    

    void AnyIamCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
        const std::vector<std::complex<double> > &dTarget_df,
        const std::vector<bool> &countAtomContribution)
    {
        if (mUseLineAlgorithm)
        {
            mCalculator2->update(atoms);
            mCalculator2->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
        }
        else
        {
            mCalculator->update(atoms);
            mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
        }

    }

    void AnyIamCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        const std::vector<bool> &countAtomContribution)
    {
        bool skip = true;
        for (bool count : countAtomContribution)
            if (count)
                skip = false;
        if (skip)
        {
            f.assign(hkl.size(), 0.0);
            return;
        }
        //vector<complex<double> > fake_dTarget_df(hkl.size());
        //std::vector<TargetFunctionAtomicParamDerivatives> dTarget_dparam;
        //mCalculator->update(atoms);
        //mCalculator->calculateStructureFactorsAndDerivatives(hkl, f, dTarget_dparam, fake_dTarget_df, countAtomContribution);
        if (mUseLineAlgorithm)
        {
            cout << "use line algorithm\n";
            mCalculator2->update(atoms);
            mCalculator2->calculateStructureFactors(hkl, f, countAtomContribution);
        }
        else
        {
            mCalculator->update(atoms);
            mCalculator->calculateStructureFactors(hkl, f, countAtomContribution);
        }

    }

	void AnyIamCalculator::calculateFormFactors(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
        if (mUseLineAlgorithm)
		    mCalculator2->calculateFormFactors(hkl, formFactors, includeAtom);
        else
            mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
	}


    void AnyIamCalculator::update(const std::vector<AtomInCrystal> &atoms)
    {
        //add_to_log(__LINE__, "AnyIamCalculator::update called");
        mCalculator->update(atoms);
    }

    void AnyIamCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i &hkl,
        std::complex<double> &scatteringFactor,
        discamb::SfDerivativesAtHkl &derivatives,
        const std::vector<bool> &countAtomContribution)
    {
        if (mUseLineAlgorithm)
            mCalculator2->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives,
                countAtomContribution);
        else
            mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives,
                countAtomContribution);

    }


}