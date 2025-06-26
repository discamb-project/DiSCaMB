#include "discamb/Scattering/TaamSfCalculatorMultiOrderedImpl.h"


#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/HC_Model/DeformationValenceParameters.h"

#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/Scattering/AnyHcCalculator.h"
#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/IO/cif_io.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/Scattering/disordered_structure_fragments.h"
#include "discamb/Scattering/taam_utilities.h"


#include <omp.h>

#include <fstream>
#include <sstream>
#include <set>
#include <iomanip>
#include <filesystem>

using namespace std;

namespace discamb {



    

    TaamSfCalculatorMultiOrderedImpl::TaamSfCalculatorMultiOrderedImpl(
        const Crystal &crystal,
        const std::vector < std::vector <std::pair<std::string, double> > >& atomList,
        const std::vector<AtomType> &atomTypes,
        const std::vector<AtomTypeHC_Parameters> &parameters,
        bool electronScattering,
        const DescriptorsSettings &settings,
        const std::string &assignemntInfoFile,
        const std::string &assignmentCsvFile,
        const std::string& parametersInfoFile,
        const std::string& multipolarCif,
        int nThreads,
        double unitCellCharge,
        bool scaleToMatchCharge,
        const string& iamTable,
        bool iamElectronScattering,
        bool frozen_lcs,
        const std::string& algorithm)
    {
        set(crystal, atomList, atomTypes, parameters, electronScattering, settings, assignemntInfoFile, assignmentCsvFile,
            parametersInfoFile, multipolarCif, nThreads, unitCellCharge, scaleToMatchCharge, iamTable, iamElectronScattering, frozen_lcs, algorithm);

    }

    TaamSfCalculatorMultiOrderedImpl::TaamSfCalculatorMultiOrderedImpl(
        const Crystal& crystal,
        const nlohmann::json& data,
        const std::string& bankString/*,
        std::string& assignemntInfo,
        bool generateAssignementInfo*/)
    {
        set(crystal, data, bankString);
    }

    void TaamSfCalculatorMultiOrderedImpl::setReconstructionData(
        const Crystal& crystal,
        const std::vector < std::vector <std::pair<std::string, double> > >& atomList)
    {
        map<string, int> label2idx;
        int atomIdx, nAtoms = crystal.atoms.size();
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            label2idx[crystal.atoms[atomIdx].label] = atomIdx;

        mReconstructionData.resize(nAtoms);
        int nSubcrystals = atomList.size();
        
        //std::vector<std::vector<int> > mSubcrystal2CrystalAtomsMap;
        mSubcrystal2CrystalAtomsMap.clear();
        mSubcrystal2CrystalAtomsMap.resize(nSubcrystals);
        mSubcrystals.assign(nSubcrystals, crystal);
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
        {
            mSubcrystals[subcrystalIdx].atoms.clear();
            for (int i = 0; i < atomList[subcrystalIdx].size(); i++)
            {
                string atomLabel = atomList[subcrystalIdx][i].first;
                if(label2idx.count(atomLabel)==0)
                    on_error::throwException("incorrect atom label in subcrystal definition for Taam Multi-Order structure factor calculator ", __FILE__, __LINE__);

                int idxInCrystal = label2idx[atomLabel];
                mSubcrystal2CrystalAtomsMap[subcrystalIdx].push_back(idxInCrystal);
                double weight = atomList[subcrystalIdx][i].second;
                mSubcrystals[subcrystalIdx].atoms.push_back(crystal.atoms[idxInCrystal]);
                mReconstructionData[idxInCrystal].atoms.push_back({ subcrystalIdx, i });
                mReconstructionData[idxInCrystal].weights.push_back(weight);
            }
        }
        // normalize weights
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            double sum = 0;
            int nInstances = mReconstructionData[atomIdx].weights.size();
            if (nInstances == 0)
                on_error::throwException("atom " + crystal.atoms[atomIdx].label + " is absent in TAAM subcrystals ", __FILE__, __LINE__);
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
                sum += mReconstructionData[atomIdx].weights[instanceIdx];
            if (sum != 0)
                for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
                    mReconstructionData[atomIdx].weights[instanceIdx] /= sum;
            // weight subcrystal atoms
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                int subcrystalIdx = mReconstructionData[atomIdx].atoms[instanceIdx].first;
                int atomInSubcrystal = mReconstructionData[atomIdx].atoms[instanceIdx].second;
                mSubcrystals[subcrystalIdx].atoms[atomInSubcrystal].occupancy *= mReconstructionData[atomIdx].weights[instanceIdx];
            }

        }
        mSubcrystalAtomsWeight.resize(nSubcrystals);
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            mSubcrystalAtomsWeight[subcrystalIdx].assign(mSubcrystals[subcrystalIdx].atoms.size(), 0.0);

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int nInstances = mReconstructionData[atomIdx].weights.size();
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                int subcrystalIdx = mReconstructionData[atomIdx].atoms[instanceIdx].first;
                int atomInSubcrystal = mReconstructionData[atomIdx].atoms[instanceIdx].second;
                mSubcrystalAtomsWeight[subcrystalIdx][atomInSubcrystal] = mReconstructionData[atomIdx].weights[instanceIdx];
            }

        }

    }


    void TaamSfCalculatorMultiOrderedImpl::set(
        const Crystal& crystal,
        const nlohmann::json& data,
        const std::string& bankString)
    {

        bool electronScattering = false;
        bool writeToDiscamb2tscLog = false;
        if (data.find("electron_scattering") != data.end())
            electronScattering = data.find("electron_scattering")->get<bool>();

        if (data.find("electron scattering") != data.end())
            electronScattering = data.find("electron scattering")->get<bool>();


        string bankPath;

        if (data.find("bank path") != data.end())
            bankPath = data.find("bank path")->get<string>();

        string assignmentInfoFile;
        if (data.find("assignment info") != data.end())
            assignmentInfoFile = data.find("assignment info")->get<string>();

        string assignmentCsvFile = data.value("assignment csv", string());

        string parametersInfoFile;
        if (data.find("parameters info") != data.end())
            parametersInfoFile = data.find("parameters info")->get<string>();

        string multipolarCif;
        if (data.find("multipole cif") != data.end())
            multipolarCif = data.find("multipole cif")->get<string>();

        double unitCellCharge = 0;
        if (data.find("unit cell charge") != data.end())
            unitCellCharge = data.find("unit cell charge")->get<double>();

        bool scaleHcParameters = true;
        if (data.find("scale") != data.end())
            scaleHcParameters = data.find("scale")->get<bool>();

        int nCores = 1;
        if (data.find("n cores") != data.end())
            nCores = data.find("n cores")->get<int>();

        string iamTable;
        if (data.find("table") != data.end())
            iamTable = data.find("table")->get<string>();

        bool iamElectronScattering = false;
        if (data.find("iam electron scattering") != data.end())
            iamElectronScattering = data.find("iam electron scattering")->get<bool>();
        bool frozen_lcs = data.value("frozen lcs", false);
        string algorithm = data.value("algorithm", "standard");

        MATTS_BankReader bankReader;
        vector<AtomType> types;
        vector<AtomTypeHC_Parameters> hcParameters;
        BankSettings bankSettings;



        if (bankPath.empty())
        {
            if (!bankString.empty())
            {
                stringstream bankStream;
                bankStream << bankString;
                bankReader.read(bankStream, types, hcParameters, bankSettings, true);
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
                    bankReader.read(*bnkFiles.begin(), types, hcParameters, bankSettings, true);

            }
        }
        else
            bankReader.read(bankPath, types, hcParameters, bankSettings, true);

        vector < vector <pair<string, double> > > atomList;
        if (data.find("substructures file") != data.end())
            disordered_structure_fragments::from_file(data["substructures file"].get<string>(), atomList);
        else
        {
            atomList.resize(1);
            for (auto& atom : crystal.atoms)
                atomList[0].push_back({ atom.label, 1.0 });
        }


        set(crystal, atomList, types, hcParameters, electronScattering, DescriptorsSettings(), assignmentInfoFile, assignmentCsvFile,
            parametersInfoFile, multipolarCif, nCores, unitCellCharge, scaleHcParameters, iamTable, iamElectronScattering, frozen_lcs, algorithm);

    }

    TaamSfCalculatorMultiOrderedImpl::TaamSfCalculatorMultiOrderedImpl(
        const Crystal& crystal,
        const std::string& jsonFileName)
    {
        nlohmann::json json_data;
        ifstream jsonFileStream(jsonFileName);

        if (jsonFileStream.good())
            jsonFileStream >> json_data;
        jsonFileStream.close();

        set(crystal, json_data, string());
    }

    TaamSfCalculatorMultiOrderedImpl::TaamSfCalculatorMultiOrderedImpl(
        const Crystal &crystal, 
        const nlohmann::json &data)
    {
        set(crystal, data, string());
        
    }


    TaamSfCalculatorMultiOrderedImpl::~TaamSfCalculatorMultiOrderedImpl()
    {
    }
    
    void TaamSfCalculatorMultiOrderedImpl::getModelInformation(
        std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        modelInfo.clear();
        modelInfo.push_back({ "SCATTERING MODEL", "TAAM multiple ordered components" });
        modelInfo.push_back({ "TAAM DATABANK", "MATTS (Kumar, P., Gruza, B., Bojarowski, S. A. & Dominiak, P. M. (2019). Acta Cryst. A75, 398-408.)" });
    }



    void TaamSfCalculatorMultiOrderedImpl::set(
        const Crystal &crystal,
        const std::vector < std::vector <std::pair<std::string, double> > >& atomList,
        const std::vector<AtomType> &atomTypes,
        const std::vector<AtomTypeHC_Parameters> &bankParameters,
        bool electronScattering,
        const DescriptorsSettings &settings,
        const std::string &_assignemntInfoFile,
        const std::string &_assignmentCsvFile,
        const std::string& _parametersInfoFile,
        const std::string& _multipolarCif,
        int nThreads,
        double unitCellCharge,
        bool scaleToMatchCharge,
        const string& iamTable,
        bool iamElectronScattering,
        bool frozen_lcs,
        const std::string &algorithm
        /*bool generateAssignmentInfo*/ )
    {
        mCrystal = crystal;
        mReciprocalLatticeUnitCell.set(crystal.unitCell);
        setReconstructionData(crystal, atomList);
        mCalculators.clear();
        int idx = 1;
        for (auto& subcrystal : mSubcrystals)
        {
            string idxStr = "_" + to_string(idx++);
            string assignemntInfoFile = _assignemntInfoFile;
            if (!assignemntInfoFile.empty())
                assignemntInfoFile += idxStr;
            string assignmentCsvFile = _assignmentCsvFile;
            if (!assignmentCsvFile.empty())
                assignmentCsvFile += idxStr;
            string parametersInfoFile = _parametersInfoFile;
            if (!parametersInfoFile.empty())
                parametersInfoFile += idxStr;
            string multipolarCif = _multipolarCif;
            if (!multipolarCif.empty())
                multipolarCif += idxStr;

            mCalculators.push_back(make_shared<HcAtomBankStructureFactorCalculator>(
                subcrystal, atomTypes, bankParameters, electronScattering, settings, assignemntInfoFile,
                assignmentCsvFile, parametersInfoFile, multipolarCif, nThreads, unitCellCharge, scaleToMatchCharge, iamTable,
                iamElectronScattering, frozen_lcs, algorithm));
        }
    }

    void TaamSfCalculatorMultiOrderedImpl::setAnomalous(
        const std::vector<std::complex<double> > & anomalous) 
    {
        on_error::not_implemented(__FILE__, __LINE__);
        //mCalculator->setAnomalous(anomalous);
    }
 

    void TaamSfCalculatorMultiOrderedImpl::calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution) 
    {
        DerivativesSelector selector;
        calculateStructureFactorsAndDerivatives(
            atoms,
            hkl,
            f,
            dTarget_dparam,
            dTarget_df,
            countAtomContribution,
            selector);

    }

    void TaamSfCalculatorMultiOrderedImpl::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution,
        const DerivativesSelector& selector)
    {
        int subcrystalIdx, nSubcrystals = mSubcrystals.size();
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            for (int atomInSubcrystalIdx = 0; atomInSubcrystalIdx < mSubcrystal2CrystalAtomsMap[subcrystalIdx].size(); atomInSubcrystalIdx++)
                mSubcrystals[subcrystalIdx].atoms[atomInSubcrystalIdx] = atoms[mSubcrystal2CrystalAtomsMap[subcrystalIdx][atomInSubcrystalIdx]];
        int atomIdx, nAtoms = atoms.size();
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int nInstances = mReconstructionData[atomIdx].atoms.size();
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                auto& atom = mReconstructionData[atomIdx].atoms[instanceIdx];
                mSubcrystals[atom.first].atoms[atom.second].occupancy *= mReconstructionData[atomIdx].weights[instanceIdx];
            }
        }
        vector<vector<bool> > useAtomContribution(nSubcrystals);
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            useAtomContribution[subcrystalIdx].assign(mSubcrystals[subcrystalIdx].atoms.size(), false);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int nInstances = mReconstructionData[atomIdx].atoms.size();
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                auto& atom = mReconstructionData[atomIdx].atoms[instanceIdx];
                useAtomContribution[atom.first][atom.second] = countAtomContribution[atomIdx];
            }

        }

        // calculate structure factors and derivatives for subcrystals

        vector< vector<TargetFunctionAtomicParamDerivatives> > dTarget_dparam_subcrystal(nSubcrystals);
        vector< vector<complex<double> > > f_subcrystal(nSubcrystals);
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            mCalculators[subcrystalIdx]->calculateStructureFactorsAndDerivatives(
                mSubcrystals[subcrystalIdx].atoms, hkl, f_subcrystal[subcrystalIdx], dTarget_dparam_subcrystal[subcrystalIdx], dTarget_df, useAtomContribution[subcrystalIdx], selector);

        // merge structure factors
        int nHkl = hkl.size();
        f.assign(nHkl, 0.0);
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
                f[hklIdx] += f_subcrystal[subcrystalIdx][hklIdx];

        // merge structure factor derivatives
        dTarget_dparam.resize(nAtoms);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int adpComponent, nAdpComponents = atoms[atomIdx].adp.size();
            dTarget_dparam[atomIdx].adp_derivatives.assign(nAdpComponents, 0.0);
            dTarget_dparam[atomIdx].atomic_position_derivatives.set(0, 0, 0);
            dTarget_dparam[atomIdx].occupancy_derivatives = 0.0;

            int nInstances = mReconstructionData[atomIdx].atoms.size();
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                auto& atom = mReconstructionData[atomIdx].atoms[instanceIdx];
                int subcrystalIdx = atom.first;
                int atomInSubcrystal = atom.second;
                for (int adpComponent = 0; adpComponent < nAdpComponents; adpComponent++)
                    dTarget_dparam[atomIdx].adp_derivatives[adpComponent] +=
                    dTarget_dparam_subcrystal[subcrystalIdx][atomInSubcrystal].adp_derivatives[adpComponent];
                dTarget_dparam[atomIdx].atomic_position_derivatives +=
                    dTarget_dparam_subcrystal[subcrystalIdx][atomInSubcrystal].atomic_position_derivatives;
                dTarget_dparam[atomIdx].occupancy_derivatives +=
                    dTarget_dparam_subcrystal[subcrystalIdx][atomInSubcrystal].occupancy_derivatives * mReconstructionData[atomIdx].weights[instanceIdx];

            }

        }
    }


	void TaamSfCalculatorMultiOrderedImpl::calculateFormFactors(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
        Vector3d hklCart;
        mReciprocalLatticeUnitCell.fractionalToCartesian(hkl, hklCart);
        calculateFormFactorsCart(hklCart, formFactors, includeAtom);
	}

    void TaamSfCalculatorMultiOrderedImpl::findSubstructureContributingAtoms(
        const std::vector<bool>& countAtomContribution,
        std::vector< std::vector<bool> >& substructureAtomContribution)
        const
    {
        int subcrystalIdx, nSubcrystals = mSubcrystals.size();
        int atomIdx, nAtoms = countAtomContribution.size();
        substructureAtomContribution.resize(nSubcrystals);

        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            substructureAtomContribution[subcrystalIdx].assign(mSubcrystals[subcrystalIdx].atoms.size(), false);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int nInstances = mReconstructionData[atomIdx].atoms.size();
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                auto& atom = mReconstructionData[atomIdx].atoms[instanceIdx];
                substructureAtomContribution[atom.first][atom.second] = countAtomContribution[atomIdx];
            }
        }

    }


    void TaamSfCalculatorMultiOrderedImpl::calculateFormFactorsCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        int subcrystalIdx, nSubcrystals = mSubcrystals.size();
        int atomIdx, nAtoms = includeAtom.size();
        vector<vector<bool> > useAtomContribution;

        findSubstructureContributingAtoms(includeAtom, useAtomContribution);
        formFactors.resize(nAtoms);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (includeAtom[atomIdx])
                formFactors[atomIdx] = 0.0;

        // calculate form factors
        // [hkl idx][atom idx]
        vector< vector<complex<double> > > ff_subcrystal(nSubcrystals);
        for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            mCalculators[subcrystalIdx]->calculateFormFactorsCart(
                hkl, ff_subcrystal[subcrystalIdx], useAtomContribution[subcrystalIdx]);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {

            int nInstances = mReconstructionData[atomIdx].atoms.size();
            for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
            {
                auto& atom = mReconstructionData[atomIdx].atoms[instanceIdx];
                int subcrystalIdx = atom.first;
                int atomInSubcrystal = atom.second;

                formFactors[atomIdx] +=
                    ff_subcrystal[subcrystalIdx][atomInSubcrystal] * mReconstructionData[atomIdx].weights[instanceIdx];
            }
        }

    }


    void TaamSfCalculatorMultiOrderedImpl::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        int hklIdx, nHkl = hkl.size();
        formFactors.resize(nHkl);
        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
            calculateFormFactors(hkl[hklIdx], formFactors[hklIdx], includeAtom);
    }



    void TaamSfCalculatorMultiOrderedImpl::calculateStructureFactors(
        const std::vector<AtomInCrystal> &atoms,
        const std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &f,
        const std::vector<bool> &countAtomContribution) 
    {
        vector<TargetFunctionAtomicParamDerivatives> dTarget_dparam;
        vector<complex<double> > dTarget_df(hkl.size(),0.0);

        TaamSfCalculatorMultiOrderedImpl::calculateStructureFactorsAndDerivatives(
            atoms,
            hkl,
            f,
            dTarget_dparam,
            dTarget_df,
            countAtomContribution);
    }

        void TaamSfCalculatorMultiOrderedImpl::update(
            const std::vector<AtomInCrystal> &atoms) 
        {
            int subcrystalIdx, nSubcrystals = mSubcrystalsForUpdate.size();

            for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            {
                for (int atomInSubcrystalIdx = 0; atomInSubcrystalIdx < mSubcrystal2CrystalAtomsMap[subcrystalIdx].size(); atomInSubcrystalIdx++)
                    mSubcrystalsForUpdate[subcrystalIdx].atoms[atomInSubcrystalIdx] = atoms[mSubcrystal2CrystalAtomsMap[subcrystalIdx][atomInSubcrystalIdx]];
                mCalculators[subcrystalIdx]->update(mSubcrystalsForUpdate[subcrystalIdx].atoms);
            }

        }

        void TaamSfCalculatorMultiOrderedImpl::calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &scatteringFactor,
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution) 
        {
            on_error::not_implemented(__FILE__, __LINE__);
            
            int subcrystalIdx, nSubcrystals = mSubcrystals.size();
            int atomIdx, nAtoms = countAtomContribution.size();
            vector<vector<bool> > useAtomContribution;

            findSubstructureContributingAtoms(countAtomContribution, useAtomContribution);

            
            scatteringFactor = 0.0;

            derivatives.adpDerivatives.resize(nAtoms);
            derivatives.anomalousScatteringDerivatives.resize(nAtoms);
            derivatives.atomicPostionDerivatives.resize(nAtoms);
            derivatives.occupancyDerivatives.resize(nAtoms);

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (countAtomContribution[atomIdx])
                {
                    derivatives.adpDerivatives[atomIdx].assign(mCrystal.atoms[atomIdx].adp.size(),0.0);
                    derivatives.anomalousScatteringDerivatives[atomIdx] = 0.0;
                    derivatives.atomicPostionDerivatives[atomIdx] = Vector3d(0.0,0.0,0.0);
                    derivatives.occupancyDerivatives[atomIdx] = 0.0;
                }
                    

            // calculate form factors
            // [hkl idx][atom idx]
            vector< vector<complex<double> > > ff_subcrystal(nSubcrystals);
            SfDerivativesAtHkl subcrystalDerivatives;
            complex<double> subcrystalScatteringFactor;
            for (int subcrystalIdx = 0; subcrystalIdx < nSubcrystals; subcrystalIdx++)
            {
                mCalculators[subcrystalIdx]->calculateStructureFactorsAndDerivatives(
                    hkl, subcrystalScatteringFactor, subcrystalDerivatives, useAtomContribution[subcrystalIdx]);

                scatteringFactor += subcrystalScatteringFactor;
                
                int nAtomsSubcrystal = mSubcrystals[subcrystalIdx].atoms.size();
                for (int atomIdxSubcrystal = 0; atomIdxSubcrystal < nAtomsSubcrystal; atomIdxSubcrystal++)
                {
                    int atomIdx = mSubcrystal2CrystalAtomsMap[subcrystalIdx][atomIdxSubcrystal];
                    int nAdpComponents = derivatives.adpDerivatives[atomIdx].size();
                    for(int adpComonent=0; adpComonent<nAdpComponents; adpComonent++)
                        derivatives.adpDerivatives[atomIdx][adpComonent] += subcrystalDerivatives.adpDerivatives[atomIdxSubcrystal][adpComonent];
                    
                    derivatives.atomicPostionDerivatives[atomIdx] += subcrystalDerivatives.atomicPostionDerivatives[atomIdxSubcrystal];
                    derivatives.occupancyDerivatives[atomIdx] += subcrystalDerivatives.occupancyDerivatives[atomIdxSubcrystal] * 
                        mSubcrystalAtomsWeight[subcrystalIdx][atomIdxSubcrystal];
                }
            }
        }


}
 