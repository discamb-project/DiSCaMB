#include "discamb/Scattering/StockholderAtomSfCalculator.h"

#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/discamb_env.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/mol2_io.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/QuantumChemistry/distributed_multipoles.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/QuantumChemistry/ProatomDB.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/gar_utilities.h"
#include "discamb/StructuralProperties/SimpleAIF_ScoreCalculator.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include "omp.h"

#include <fstream>
#include <utility>
#include <set>
#include <filesystem>


using namespace std;


namespace discamb {






    void StockholderAtomSfCalculator::set(
        const Crystal& crystal,
        const HirshfeldAtomModelSettings& settings,
        bool electronScattering,
        const std::string& jobName)
    {
        mHamSettings = settings;
        mCrystal = crystal;

        // Hirshfeld partition specific
        string atomFile;

        if (mHamSettings.electronDensityPartition.type == ElectronDensityPartitionType::Hirshfeld)
        {
            auto& partitionData = mHamSettings.electronDensityPartition.partitionSpecificData;
            if (partitionData.find("atoms file") != partitionData.end())
            {
                atomFile = partitionData["atoms file"].get<string>();

                if (!atomFile.empty())
                    if (!filesystem::exists(atomFile))
                        on_error::throwException("atomic spherical densities file '" + atomFile +
                            "' does not exists", __FILE__, __LINE__);
            }
        }
        //------------------
        
        mModelInfo.clear();
        mModelInfo.push_back({ "SCATTERING MODEL", "HAR" });
        mModelInfo.push_back({ "QUANTUM CHEMISTRY METHOD", mHamSettings.wfnCalculation.qmSettings.qmMethod });
        if(!mHamSettings.wfnCalculation.qmSettings.relativisticMethod.empty())
            mModelInfo.push_back({ "QUANTUM CHEMISTRY RELATIVISTIC METHOD", mHamSettings.wfnCalculation.qmSettings.relativisticMethod });
        mModelInfo.push_back({ "BASIS SET", mHamSettings.wfnCalculation.qmSettings.basisSet });
            


        // --------------------------

        mJobName = jobName;



        mFormFactorManager = std::shared_ptr<StockholderAtomFormFactorCalcManager>(new StockholderAtomFormFactorCalcManager(crystal, mHamSettings));

        mManager = static_pointer_cast<AtomicFormFactorCalculationsManager>(mFormFactorManager);
        
        mCalculator = new AnyScattererStructureFactorCalculator(crystal);

        if (electronScattering)
        {
            vector<double> nuclearCharges;
            vector<int> atomicNumbers;
            crystal_structure_utilities::atomicNumbers(mCrystal, atomicNumbers);
            for (int z : atomicNumbers)
                nuclearCharges.push_back(double(z));
            //nuclearCharges = atomicNumbers;
            //for (int i = 0; i < nuclearCharges.size(); i++)
            //    nuclearCharges[i] = double(basic_chemistry_utilities::atomicNumberFromLabel(mCrystal.atoms[i].type));


            std::shared_ptr<AtomicFormFactorCalculationsManager> manager =
                std::shared_ptr<AtomicFormFactorCalculationsManager>(
                    new ElectronFromXrayFormFactorCalculationManager(crystal.unitCell, nuclearCharges, mManager));
            mCalculator->setAtomicFormfactorManager(manager);
        }
        else
            mCalculator->setAtomicFormfactorManager(mManager);

        if (mHamSettings.diagnostics.prerun)
            exit(0);
        cout << "done\n";
        cout << "calculate atomic densities\n";
        
        mFormFactorManager->calculateAtomicDensities();
        cout << "done\n";
    }

    void StockholderAtomSfCalculator::getModelInformation(
        vector< pair< string, string> >& modelInfo) 
        const 
    {
        modelInfo = mModelInfo;
    }

    void StockholderAtomSfCalculator::setAnomalous(const std::vector<std::complex<double> >& anomalous)
    {
        mCalculator->setAnoumalous(anomalous);
    }


    StockholderAtomSfCalculator::StockholderAtomSfCalculator(
        const Crystal& crystal, const HirshfeldAtomModelSettings& settings, 
        bool electronScattering, const std::string& jobname)
    {
        set(crystal, settings, electronScattering, jobname);        
    }

    StockholderAtomSfCalculator::StockholderAtomSfCalculator(
        const Crystal& crystal,
        const nlohmann::json& data)
    {
        HirshfeldAtomModelSettings settings;
        settings.set(data, crystal);
        string jobName = data.value("job name", "noname");
        bool electronScattering = data.value("electron_scattering",false);
        electronScattering = data.value("electron scattering", electronScattering);
        
        set(crystal, settings, electronScattering, jobName);

    }



    StockholderAtomSfCalculator::~StockholderAtomSfCalculator()
    {
        delete mCalculator;
    }


    void StockholderAtomSfCalculator::calculateStructureFactorsAndDerivatives(
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

    void StockholderAtomSfCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        const std::vector<bool>& countAtomContribution)
    {
        mCalculator->update(atoms);
        mCalculator->calculateStructureFactors(hkl, f, countAtomContribution);
    }

    void StockholderAtomSfCalculator::update(const std::vector<AtomInCrystal>& atoms)
    {
        mCalculator->update(atoms);
        mFormFactorManager->calculateAtomicDensities();
    }

    void StockholderAtomSfCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i& hkl,
        std::complex<double>& scatteringFactor,
        discamb::SfDerivativesAtHkl& derivatives,
        const std::vector<bool>& countAtomContribution)
    {
        mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives, countAtomContribution);
    }

    void StockholderAtomSfCalculator::calculateFormFactors(
        const Vector3i& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    void StockholderAtomSfCalculator::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }
                                                    


}

