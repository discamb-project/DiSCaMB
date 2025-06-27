#include "discamb/Scattering/TaamSfCalculator.h"

#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/Scattering/disordered_structure_fragments.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator2.h"
#include "discamb/Scattering/TaamSfCalculatorMultiOrderedImpl.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace discamb {

    
    void TaamSfCalculatorSettings::set(
        const Crystal &crystal,
        const nlohmann::json& data,
        const std::string& bankText)
    {
        electronScattering = data.value("electron_scattering", false);
        electronScattering = data.value("electron scattering", false);        
        string bankFilePath = data.value("bank path",string());
        assignmentInfoFile = data.value("assignment info",string());
        assignmentCsvFile = data.value("assignment csv", string());
        parametersInfoFile = data.value("parameters info", string());
        multipolarCif = data.value("multipole cif", string());
        unitCellCharge = data.value("unit cell charge", 0.0);
        scaleToMatchCharge = data.value("scale", true);
        nThreads = data.value("n cores", 1);
        iamTable = data.value("table", string());
        iamElectronScattering = data.value("iam electron scattering", false); 
        frozen_lcs = data.value("frozen lcs", false);
        algorithm = data.value("algorithm", "standard");
        
        MATTS_BankReader bankReader;
        //vector<AtomType> types;
        //vector<AtomTypeHC_Parameters> hcParameters;
        BankSettings bankSettings;



        if (bankFilePath.empty())
        {
            if (!bankText.empty())
            {
                stringstream bankStream;
                bankStream << bankText;
                bankReader.read(bankStream, atomTypes, parameters, bankSettings, true);
                descriptorsSettings = bankSettings.descriptorsSettings;
            }
            else
            {
                // look for *.bnk file

                std::set<string> bnkFiles;
                string extension;
                for (auto it : filesystem::directory_iterator(filesystem::current_path()))
                    if (filesystem::is_regular_file(it.status()))
                    {
                        string_utilities::toLower(it.path().extension().string(), extension);
                        if (extension == string(".bnk"))
                            bnkFiles.insert(it.path().filename().string());
                    }

                if (bnkFiles.empty())
                    on_error::throwException("no bank file specified or found in working directory", __FILE__, __LINE__);
                else
                {
                    bankReader.read(*bnkFiles.begin(), atomTypes, parameters, bankSettings, true);
                    descriptorsSettings = bankSettings.descriptorsSettings;
                }

            }
        }
        else
        {
            bankReader.read(bankFilePath, atomTypes, parameters, bankSettings, true);
            descriptorsSettings = bankSettings.descriptorsSettings;
        }

        
        splitWithLabels = data.value("split with labels", false);
        splitWithInternalAltlocLabels = data.value("split with internal altloc labels", false);

        if (data.find("substructures file") != data.end())
            disordered_structure_fragments::from_file(data["substructures  file"].get<string>(), orderedSubcrystalAtoms);
        // to fix for gcc
        if (data.find("substructures") != data.end())
        {
            //auto & substructures = data.find("substructures");
            disordered_structure_fragments::substructures_from_json(data["substructures"], orderedSubcrystalAtoms);
        }

        if(data.find("substructures from disorder groups") != data.end())
        {
            auto& disorder_groups_json = data["substructures from disorder groups"];
            disordered_structure_fragments::split_structure_json_str(
                crystal, disorder_groups_json, orderedSubcrystalAtoms);
        }

        //vector<disordered_structure_fragments::Fragment> taamFragments;
        string fragmentsFile = data.value("fragments file", string());
        if (!fragmentsFile.empty())
            disordered_structure_fragments::from_file(crystal, fragmentsFile, taamFragments);



    }

    void TaamSfCalculatorSettings::set(
        const Crystal& crystal,
        const std::string &fileName,
        const std::string& bankText)
    {
        nlohmann::json json_data;
        ifstream jsonFileStream("aspher.json");

        if (jsonFileStream.good())
            jsonFileStream >> json_data;
        jsonFileStream.close();

        set(crystal, json_data, bankText);
    }


    TaamSfCalculator::TaamSfCalculator(
        const Crystal& crystal,
        const TaamSfCalculatorSettings& settings)
    {
        set(crystal, settings);
    }

    TaamSfCalculator::TaamSfCalculator(
        const Crystal& crystal, 
        const nlohmann::json& data)
    {
        TaamSfCalculatorSettings settings;
        settings.set(crystal, data);
        set(crystal, settings);
    }

    TaamSfCalculator::TaamSfCalculator(
        const Crystal& crystal, 
        const std::string& jsonFileName)
    {
        TaamSfCalculatorSettings settings;
        settings.set(crystal, jsonFileName);
        set(crystal, settings);
    }


    TaamSfCalculator::~TaamSfCalculator() {}

    void TaamSfCalculator::getModelInformation(
        std::vector<std::pair<std::string, std::string> >& modelInfo) 
        const
    {
        modelInfo.clear();
        modelInfo.push_back({ "SCATTERING MODEL", "TAAM" });
        modelInfo.push_back({ "TAAM DATABANK", "MATTS (Kumar, P., Gruza, B., Bojarowski, S. A. & Dominiak, P. M. (2019). Acta Cryst. A75, 398-408.)" });
    }

    void TaamSfCalculator::set(
        const Crystal& crystal,
        const TaamSfCalculatorSettings& settings)
    {
        vector<vector<pair<string, double> > > orderedSubcrystalAtoms = settings.orderedSubcrystalAtoms;

        if (settings.splitWithLabels && settings.orderedSubcrystalAtoms.empty())
            disordered_structure_fragments::split_with_labels(crystal, orderedSubcrystalAtoms);
        
        if (settings.splitWithInternalAltlocLabels && settings.orderedSubcrystalAtoms.empty())
            disordered_structure_fragments::split_with_labels_internal_altloc(crystal, orderedSubcrystalAtoms);

        if (!orderedSubcrystalAtoms.empty())
        {
            //cout << "use TaamSfCalculatorMultiOrderedImpl" << endl;
            auto impl = make_shared<TaamSfCalculatorMultiOrderedImpl>(
                crystal,
                orderedSubcrystalAtoms,
                settings.atomTypes,
                settings.parameters,
                settings.electronScattering,
                settings.descriptorsSettings,
                settings.assignmentInfoFile,
                settings.assignmentCsvFile,
                settings.parametersInfoFile,
                settings.multipolarCif,
                settings.nThreads,
                settings.unitCellCharge,
                settings.scaleToMatchCharge,
                settings.iamTable,
                settings.iamElectronScattering,
                settings.frozen_lcs,
                settings.algorithm);
            mImplementation = impl;
            return;
        }
        if (!settings.taamFragments.empty())
        {
            //cout << "use HcAtomBankStructureFactorCalculator2" << endl;
            mImplementation = make_shared<HcAtomBankStructureFactorCalculator2>(
                crystal,
                settings.atomTypes,
                settings.parameters,
                settings.electronScattering,
                settings.descriptorsSettings,
                settings.assignmentInfoFile,
                settings.assignmentCsvFile,
                settings.parametersInfoFile,
                settings.multipolarCif,
                settings.nThreads,
                settings.unitCellCharge,
                settings.scaleToMatchCharge,
                settings.iamTable,
                settings.iamElectronScattering,
                settings.frozen_lcs,
                settings.algorithm,
                settings.taamFragments);
        }
        else
        {
            //cout << "use HcAtomBankStructureFactorCalculator" << endl;
            mImplementation = make_shared<HcAtomBankStructureFactorCalculator>(
                crystal,
                settings.atomTypes,
                settings.parameters,
                settings.electronScattering,
                settings.descriptorsSettings,
                settings.assignmentInfoFile,
                settings.assignmentCsvFile,
                settings.parametersInfoFile,
                settings.multipolarCif,
                settings.nThreads,
                settings.unitCellCharge,
                settings.scaleToMatchCharge,
                settings.iamTable,
                settings.iamElectronScattering,
                settings.frozen_lcs,
                settings.algorithm);
        }
    }

    
    void TaamSfCalculator::setAnomalous(
        const std::vector<std::complex<double> >& anomalous)
    {
        mImplementation->setAnomalous(anomalous);
    }

    void TaamSfCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution)
    {
        mImplementation->calculateStructureFactorsAndDerivatives(atoms, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
    }

    void TaamSfCalculator::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution,
        const DerivativesSelector& selector) 
    {
        mImplementation->calculateStructureFactorsAndDerivatives(atoms, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution, selector);
    }


    void TaamSfCalculator::calculateStructureFactors(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        const std::vector<bool>& countAtomContribution) 
    {
        mImplementation->calculateStructureFactors(atoms, hkl, f, countAtomContribution);
    }

    void TaamSfCalculator::update(
        const std::vector<AtomInCrystal>& atoms) 
    {
        mImplementation->update(atoms);
    }

    void TaamSfCalculator::calculateStructureFactorsAndDerivatives(
        const Vector3i& hkl,
        std::complex<double>& scatteringFactor,
        discamb::SfDerivativesAtHkl& derivatives,
        const std::vector<bool>& countAtomContribution)
    {
        mImplementation->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives, countAtomContribution);
    }


    void TaamSfCalculator::calculateFormFactors(
        const Vector3i& hkl, 
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom) const
    {
        mImplementation->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    void TaamSfCalculator::calculateFormFactors(
        const std::vector<Vector3i>& hkl, 
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom) 
        const
    {
        mImplementation->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    void TaamSfCalculator::calculateFormFactorsCart(
        const Vector3d& hkl, 
        std::vector<std::complex<double> >& formFactors, 
        const std::vector<bool>& includeAtom) 
        const
    {
        mImplementation->calculateFormFactorsCart(hkl, formFactors, includeAtom);
    }


}
