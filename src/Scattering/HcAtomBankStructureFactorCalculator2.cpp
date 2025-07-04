#include "discamb/Scattering/HcAtomBankStructureFactorCalculator2.h"


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
#include "discamb/Scattering/taam_utilities.h"


#include <omp.h>

#include <fstream>
#include <sstream>
#include <set>
#include <iomanip>
#include <filesystem>

using namespace std;

namespace discamb {

    namespace {
        bool isSphericalType(const string& s)
        {
            int nZeros = 0;
            int nDigits = 0;
            for (auto c : s)
            {
                if (c == '0')
                    nZeros++;
                if (isdigit(c))
                    nDigits++;
            }
            if (nDigits == nZeros && nDigits == 3)
                return true;
            return false;
        }
    }

    

    HcAtomBankStructureFactorCalculator2::HcAtomBankStructureFactorCalculator2(
        const Crystal &crystal,
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
        const std::string& algorithm,
        const std::vector<disordered_structure_fragments::Fragment>& taamFragments)
    {
        set(crystal, atomTypes, parameters, electronScattering, settings, assignemntInfoFile, assignmentCsvFile,
            parametersInfoFile, multipolarCif, nThreads, unitCellCharge, scaleToMatchCharge, iamTable, iamElectronScattering, frozen_lcs, algorithm, taamFragments);

    }

    HcAtomBankStructureFactorCalculator2::HcAtomBankStructureFactorCalculator2(
        const Crystal& crystal,
        const nlohmann::json& data,
        const std::string& bankString/*,
        std::string& assignemntInfo,
        bool generateAssignementInfo*/)
    {
        set(crystal, data, bankString);
    }

    void HcAtomBankStructureFactorCalculator2::set(
        const Crystal& crystal,
        const nlohmann::json& data,
        const std::string& bankString/*,
        bool generateAssignementInfo*/)
    {
        bool electronScattering = false;
        //bool writeToDiscamb2tscLog = false;
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

        vector<disordered_structure_fragments::Fragment> taamFragments;
        string fragmentsFile = data.value("fragments file", string());
        if (!fragmentsFile.empty())
            disordered_structure_fragments::from_file(crystal, fragmentsFile, taamFragments);

        set(crystal, types, hcParameters, electronScattering, DescriptorsSettings(), assignmentInfoFile, assignmentCsvFile,
            parametersInfoFile, multipolarCif, nCores, unitCellCharge, scaleHcParameters, iamTable, iamElectronScattering, 
            frozen_lcs, algorithm, taamFragments);

    }

    
    HcAtomBankStructureFactorCalculator2::HcAtomBankStructureFactorCalculator2(
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

    HcAtomBankStructureFactorCalculator2::HcAtomBankStructureFactorCalculator2(
        const Crystal &crystal, 
        const nlohmann::json &data)
    {
        set(crystal, data, string());
        
  //      bool electronScattering = false;
  //      bool writeToDiscamb2tscLog = false;
  //      if (data.find("electron_scattering") != data.end())
  //          electronScattering = data.find("electron_scattering")->get<bool>();

  //      if (data.find("electron scattering") != data.end())
  //          electronScattering = data.find("electron scattering")->get<bool>();


  //      string bankPath;

  //      if (data.find("bank path") != data.end())
  //          bankPath = data.find("bank path")->get<string>();

  //      string assignmentInfoFile;
  //      if (data.find("assignment info") != data.end())
  //          assignmentInfoFile = data.find("assignment info")->get<string>();

  //      string assignmentCsvFile = data.value("assignment csv", string());

  //      string parametersInfoFile;
  //      if (data.find("parameters info") != data.end())
  //          parametersInfoFile = data.find("parameters info")->get<string>();

  //      string multipolarCif;
  //      if (data.find("multipole cif") != data.end())
  //          multipolarCif = data.find("multipole cif")->get<string>();

  //      double unitCellCharge = 0;
  //      if (data.find("unit cell charge") != data.end())
  //          unitCellCharge = data.find("unit cell charge")->get<double>();

  //      bool scaleHcParameters = true;
  //      if (data.find("scale") != data.end())
  //          scaleHcParameters = data.find("scale")->get<bool>();

  //      int nCores=1;
  //      if (data.find("n cores") != data.end())
  //          nCores = data.find("n cores")->get<int>();

  //      string iamTable;
  //      if (data.find("table") != data.end())
  //          iamTable = data.find("table")->get<string>();
  //      
  //      bool iamElectronScattering = false;
  //      if (data.find("iam electron scattering") != data.end())
  //          iamElectronScattering = data.find("iam electron scattering")->get<bool>();
  //      bool frozen_lcs = data.value("frozen lcs", false);
  //      mAlgorithm = data.value("algorithm", "standard");

  //      MATTS_BankReader bankReader;
  //      vector<AtomType> types;
  //      vector<AtomTypeHC_Parameters> hcParameters;
  //      BankSettings bankSettings;


		//if(bankPath.empty())
		//{
  //          // look for *.bnk file
  //          
  //          std::set<string> bnkFiles;
  //          string extension;
  //          for (auto it : filesystem::directory_iterator(filesystem::current_path()))
  //              if (filesystem::is_regular_file(it.status()))
  //              {
  //                  string_utilities::toLower(it.path().extension().string(), extension);
  //                  if (extension == string(".bnk"))
  //                      bnkFiles.insert(it.path().filename().string());
  //              }

  //          if (bnkFiles.empty())
  //          {
  //              on_error::throwException("no bank file specified or found in working directory", __FILE__, __LINE__);
  //              //string bankString;
  //              //stringstream bankStream;
  //              //default_ubdb_bank_string(bankString);
  //              //bankStream << bankString;
  //              //bankReader.read(bankStream, types, hcParameters, bankSettings, true);
  //          }
  //          else
  //          {
  //              bankReader.read(*bnkFiles.begin(), types, hcParameters, bankSettings, true);
  //          }
		//}
		//else
		//	bankReader.read(bankPath, types, hcParameters, bankSettings, true);
  //      set(crystal,types, hcParameters, electronScattering, DescriptorsSettings(), assignmentInfoFile, assignmentCsvFile,
  //          parametersInfoFile, multipolarCif, nCores, unitCellCharge, scaleHcParameters, iamTable, iamElectronScattering, frozen_lcs);
    }


    HcAtomBankStructureFactorCalculator2::~HcAtomBankStructureFactorCalculator2()
    {
    }
    
    void HcAtomBankStructureFactorCalculator2::getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        modelInfo = mModelInfo;
    }

    void HcAtomBankStructureFactorCalculator2::setExtendedCrystal(
        const Crystal& crystal,
        const std::vector<disordered_structure_fragments::Fragment>& taamFragments)
    {
        mExtended2normal.clear();
        mNormal2extended.clear();
        mExtendedCrystal = crystal;
        int nAtoms = crystal.atoms.size();
        mNormal2extended.resize(nAtoms);
        mAtomContributionWeights.resize(nAtoms);
        

        if (taamFragments.empty())
        {
            for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
            {
                mNormal2extended[atomIdx].push_back(atomIdx);
                mAtomContributionWeights[atomIdx].push_back(1.0);
                mExtended2normal.push_back(atomIdx);
            }
            return;
        }
        //vector<vector<int> > normal2extended;
        //vector<int> extended2normal;mNormal2extended
        mExtendedCrystal.atoms.clear();
        
        for (int fragmentIdx = 0; fragmentIdx < taamFragments.size(); fragmentIdx++)
        {
            int nAtoms = taamFragments[fragmentIdx].atomRelativeWeights.size();
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (taamFragments[fragmentIdx].atomRelativeWeights[atomIdx] > 0)
                {
                    SpaceGroupOperation op(taamFragments[fragmentIdx].atomList[atomIdx].second);
                    if(!op.isIdentity())
                        on_error::throwException("cannot use atoms with non-identity space group operation and non-zero weight in TAAM fragment", __FILE__, __LINE__);
                    int idxInCrystal = taamFragments[fragmentIdx].atomList[atomIdx].first;
                    mExtendedCrystal.atoms.push_back(crystal.atoms[idxInCrystal]);
                    mExtendedCrystal.atoms.back().label = crystal.atoms[idxInCrystal].label + "_" + to_string(fragmentIdx+1);   
                    mNormal2extended[idxInCrystal].push_back(mExtendedCrystal.atoms.size() - 1);
                    int idxExtended = mExtendedCrystal.atoms.size()-1;
                    mAtomContributionWeights[idxInCrystal].push_back(taamFragments[fragmentIdx].atomRelativeWeights[atomIdx]);
                    //mExtendedAtomsContributionWeights
                    mExtended2normal.push_back(idxInCrystal);
                    //mNormal2extendedByFrag[idxInCrystal][fragmentIdx] = mExtendedCrystal.atoms.size() - 1;
                }
        }
        // normalize weights
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            double sum = 0.0;
            for(int i=0; i< mAtomContributionWeights[atomIdx].size();i++)
                sum += mAtomContributionWeights[atomIdx][i];
            
            for (int i = 0; i < mAtomContributionWeights[atomIdx].size(); i++)
                if (sum == 0.0)
                    mAtomContributionWeights[atomIdx][i] = 0;
                else
                    mAtomContributionWeights[atomIdx][i] /= sum;
        }
    }

    AtomInCrystalID HcAtomBankStructureFactorCalculator2::regularToExtendedCrystalAtomId(
        const AtomInCrystalID& atom) const
    {
        return AtomInCrystalID(mNormal2extended[atom.index()][0], atom.getSymmetryOperation().string());
    }

    void HcAtomBankStructureFactorCalculator2::assignAtomTypes(
        const Crystal& crystal,
        const std::vector<disordered_structure_fragments::Fragment>& taamFragments,
        const std::vector<AtomType>& atomTypes,
        const DescriptorsSettings& settings,
        std::vector < LocalCoordinateSystem<AtomInCrystalID> > &lcs,
        std::vector<int> &types,
        const std::string& assignemntInfoFile,
        const std::string& assignmentCsvFile)
    {
        lcs.clear();
        types.clear();

        CrystalAtomTypeAssigner assigner;
        assigner.setAtomTypes(atomTypes);
        assigner.setDescriptorsSettings(settings);
        

        if (taamFragments.empty())
            assigner.assign(crystal, types, lcs);
        else
        {
            vector<vector<pair<int, string> > > fragments;
            vector<vector<int> > atomsToAssign(taamFragments.size());

            for (int fragmentIdx = 0; fragmentIdx < taamFragments.size(); fragmentIdx++)
            {
                fragments.push_back(taamFragments[fragmentIdx].atomList);
                for (int atomIdx = 0; atomIdx < taamFragments[fragmentIdx].atomRelativeWeights.size(); atomIdx++)
                    if (taamFragments[fragmentIdx].atomRelativeWeights[atomIdx] > 0)
                        atomsToAssign[fragmentIdx].push_back(atomIdx);
            }
            vector<vector<int> > typeIdsFrag;
            vector < vector < LocalCoordinateSystem<AtomInCrystalID> > > lcsFrag;
            assigner.assign(crystal, fragments, atomsToAssign, typeIdsFrag, lcsFrag);
            for (int fragmentIdx = 0; fragmentIdx < taamFragments.size(); fragmentIdx++)
            {
                types.insert(types.end(), typeIdsFrag[fragmentIdx].begin(), typeIdsFrag[fragmentIdx].end());
                for(auto &lcsFragAtom: lcsFrag[fragmentIdx])
                {
                    LocalCoordinateSystem<AtomInCrystalID> lcsExtended = lcsFragAtom;
                    int idxExtended = mNormal2extended[lcsFragAtom.centralAtom.index()][0];
                    lcsExtended.centralAtom = regularToExtendedCrystalAtomId(lcsFragAtom.centralAtom);
                    lcsExtended.refPoint_1.clear();
                    lcsExtended.refPoint_2.clear();
                    for(auto const &atom: lcsFragAtom.refPoint_1)
                        lcsExtended.refPoint_1.push_back(regularToExtendedCrystalAtomId(atom));
                    for (auto const& atom : lcsFragAtom.refPoint_2)
                        lcsExtended.refPoint_2.push_back(regularToExtendedCrystalAtomId(atom));
                    lcs.push_back(lcsExtended);
                }
                

            }
        }

        if (!assignemntInfoFile.empty())
        {
            if (assignemntInfoFile == string("print_to_discamb2tsc_log_file") || assignemntInfoFile == string("print_to_discambMATTS2tsc_log_file"))
            {
                vector<string> unassignedAtoms;
                for (int i = 0; i < types.size(); i++)
                {
                    if (types[i] < 0)
                        unassignedAtoms.push_back(mExtendedCrystal.atoms[i].label);
                    else
                    {
                        string typeLabel = atomTypes[types[i]].id;
                        int labelSize = typeLabel.size();
                        if (labelSize > 3)
                        {
                            string back = typeLabel.substr(labelSize - 3);
                            string front = typeLabel.substr(0, labelSize - 3);

                            if (back == string("000"))
                            {
                                bool hasDigit = false;
                                for (auto& c : front)
                                    if (isdigit(c))
                                        hasDigit = true;
                                if (!hasDigit)
                                    unassignedAtoms.push_back(crystal.atoms[i].label);
                            }
                        }
                    }
                }

                //ofstream out("discamb2tsc.log", std::ofstream::out | std::ofstream::app);
                //clog << "atomic form factor for X-ray scattering calculated with\nHansen-Coppens model parameterized with MATTS databank.\n"
                //    << "atom type assigned to " << crystal.atoms.size() - unassignedAtoms.size() << " of " << crystal.atoms.size() << " atoms\n";

                assigner.printAssignment(clog, mExtendedCrystal, types, lcs);

                //if (!unassignedAtoms.empty())
                //{
                //    clog << "atoms with unassigned atom types:\n";
                //    for (auto& label : unassignedAtoms)
                //        clog << label << "\n";
                //}
            }
            else
            {
                ofstream out(assignemntInfoFile);
                assigner.printAssignment(out, mExtendedCrystal, types, lcs);
                out.close();
            }
        }
        if (!assignmentCsvFile.empty())
        {
            ofstream out(assignmentCsvFile);
            assigner.printAssignmentCSV(out, crystal, types, lcs);
            out.close();
        }


    }


    void HcAtomBankStructureFactorCalculator2::set(
        const Crystal &crystal,
        const std::vector<AtomType> &atomTypes,
        const std::vector<AtomTypeHC_Parameters> &bankParameters,
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
        const std::string &algorithm,
        const std::vector<disordered_structure_fragments::Fragment>& taamFragments)
    {
        mFragments = taamFragments;
        mCrystal = crystal;
        setExtendedCrystal(crystal, taamFragments);
        mModelInfo.clear();
        mAlgorithm = algorithm;
        mModelInfo.push_back({ "SCATTERING MODEL", "TAAM" });
        mModelInfo.push_back({ "TAAM DATABANK", "MATTS (Kumar, P., Gruza, B., Bojarowski, S. A. & Dominiak, P. M. (2019). Acta Cryst. A75, 398-408.)" });


        mN_Threads = nThreads;

        vector < LocalCoordinateSystem<AtomInCrystalID> > lcs;
        vector<int> types;
        assignAtomTypes(crystal, taamFragments, atomTypes, settings, lcs, types, assignemntInfoFile, assignmentCsvFile);


        HC_ModelParameters multipoleModelPalameters;


        vector<int> atomicNumbers;
        vector<int> nonMultipolarAtoms;
        crystal_structure_utilities::atomicNumbers(mExtendedCrystal, atomicNumbers);

        vector<double> multiplicityTimesOccupancy;
        for (auto const & atom : mExtendedCrystal.atoms)
            multiplicityTimesOccupancy.push_back(atom.multiplicity * atom.occupancy);
        
        

        if (scaleToMatchCharge)
            taam_utilities::type_assignment_to_HC_parameters(
                bankParameters, types, multiplicityTimesOccupancy, atomicNumbers, unitCellCharge,
                multipoleModelPalameters, true, nonMultipolarAtoms);
        else
            taam_utilities::type_assignment_to_unscaled_HC_parameters(bankParameters, types, atomicNumbers,
                multipoleModelPalameters, true, nonMultipolarAtoms);

        if (!multipolarCif.empty())
            cif_io::saveCif(multipolarCif, mExtendedCrystal, multipoleModelPalameters, lcs);


        //ubdb_utilities::ubdb_type_assignment_to_HC_parameters(
        //    bankParameters, types, atomicNumbers,
        //    multipoleModelPalameters, true, nonMultipolarAtoms);

        vector<shared_ptr<LocalCoordinateSystemInCrystal> > lcaCalculators;
        for (int lcsIdx = 0; lcsIdx < lcs.size(); lcsIdx++)
            lcaCalculators.push_back(
                shared_ptr<LocalCoordinateSystemInCrystal>(
                    types[lcsIdx] >= 0 ?
                    new LocalCoordinateSystemCalculator(lcs[lcsIdx], mExtendedCrystal) :
                    new LocalCoordinateSystemCalculator()
                )
            );
        if (!parametersInfoFile.empty())
            printAssignedMultipolarParameters(
                parametersInfoFile,
                mExtendedCrystal,
                lcs,
                atomTypes,
                bankParameters,
                multipoleModelPalameters,
                types);
        if(mAlgorithm == "standard")
            mHcCalculator = shared_ptr<SfCalculator>(
                new AnyHcCalculator(mExtendedCrystal, multipoleModelPalameters, lcaCalculators, electronScattering, false, false, nThreads, frozen_lcs));
        else
            mHcCalculator = shared_ptr<SfCalculator>(
                                new AnyHcCalculator(mExtendedCrystal, multipoleModelPalameters, lcaCalculators, electronScattering, false, true, nThreads, frozen_lcs));
        mIamCalculator = shared_ptr<SfCalculator>(
            new AnyIamCalculator(mExtendedCrystal, iamElectronScattering, iamTable));
        vector< std::shared_ptr<SfCalculator> > calculators{ mHcCalculator, mIamCalculator };
        vector<vector<int> > calculatorAtoms(2);
        for (int atomIdx = 0; atomIdx < mExtendedCrystal.atoms.size(); atomIdx++)
            if (find(nonMultipolarAtoms.begin(), nonMultipolarAtoms.end(), atomIdx) == nonMultipolarAtoms.end())
                calculatorAtoms[0].push_back(atomIdx);
            else
                calculatorAtoms[1].push_back(atomIdx);
        mCalculator = shared_ptr< CombinedStructureFactorCalculator>(new CombinedStructureFactorCalculator(calculators, calculatorAtoms));
        
        //-------------

        if (nThreads == 1)
            return;

        mHcCalculators.clear();
        mIamCalculators.clear();
        mCalculators.clear();

        for (int i = 0; i < mN_Threads; i++)
        {

            mHcCalculators.push_back(shared_ptr<SfCalculator>(
                new AnyHcCalculator(mExtendedCrystal, multipoleModelPalameters, lcaCalculators, electronScattering, false, false, 1, frozen_lcs)));
            mIamCalculators.push_back(shared_ptr<SfCalculator>(
                new AnyIamCalculator(mExtendedCrystal, iamElectronScattering, iamTable)));

            vector< std::shared_ptr<SfCalculator> > calculators2{ mHcCalculators.back(), mIamCalculators.back() };
            mCalculators.push_back(shared_ptr< CombinedStructureFactorCalculator>(new CombinedStructureFactorCalculator(calculators2, calculatorAtoms)));
        }


    }


    void HcAtomBankStructureFactorCalculator2::printHcParameters(
        std::ofstream& out,
        const HC_ModelParameters& multipoleModelPalameters,
        int atomIdx)
    {
        auto const& typeParams = multipoleModelPalameters.type_parameters[multipoleModelPalameters.atom_to_type_map[atomIdx]];
        out << "kappa " << typeParams.kappa_spherical_valence << "\n"
            << "kappa_prime " << typeParams.kappa_deformation_valence << "\n"
            << "P_val " << typeParams.p_val << "\n";
        int nPlm = 0;
        for (auto x : typeParams.p_lm)
            for (auto y : x)
                if (y != 0)
                    nPlm++;
        out << "number of Plm terms " << nPlm << "\n"
            << "Plms:\n";

        for(int l=0; l < int(typeParams.p_lm.size()); l++)
            for (int m = -l; m <= l; m++)
            {
                int m_idx = l + m;
                if (typeParams.p_lm[l][m_idx] != 0)
                    out << l << " " << m << " " << typeParams.p_lm[l][m_idx] << "\n";
            }
    }

    void HcAtomBankStructureFactorCalculator2::printAssignedMultipolarParameters(
        const std::string& fName,
        const Crystal& crystal,
        const std::vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs,
        const std::vector<AtomType>& atomTypes,
        const std::vector<AtomTypeHC_Parameters>& bankParameters,
        const HC_ModelParameters& multipoleModelPalameters,
        const std::vector<int>& assignedTypeIdx)
    {
        ofstream out(fName);
        
        
        int atomIdx, nAtoms = crystal.atoms.size();
        vector<string> atomLabels;
        vector<int> atomicNumbers;
        int nAssigned = 0;
        
        
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            atomLabels.push_back(crystal.atoms[atomIdx].label);
        out << "number of atoms in asymmetric unit " << nAtoms << "\n";

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            out << "\n";

            out << "atom label " << crystal.atoms[atomIdx].label << "\n"
                << "chemical element " << periodic_table::symbol(atomicNumbers[atomIdx]) << "\n"
                << "atom type ";
            int typeIdx = assignedTypeIdx[atomIdx];
            if (typeIdx >= 0)
                out << atomTypes[typeIdx].id << "\n";
            else
                out << "unknown\n";
            
            if (typeIdx >= 0)
            {
                LocalCoordinateSystemCalculator lcsCalculator(lcs[atomIdx], crystal);
                Vector3d x, y, z;
                bool sameChirality;
                lcsCalculator.calculate(x, y, z, crystal, sameChirality);
                out << "local coordinate system:\n"
                    << "    definition " << ubdbLcsAsString(lcs[atomIdx], atomLabels);

                if (lcs[atomIdx].chirality.empty())
                    sameChirality = true;

                if (lcs[atomIdx].isR)
                    out << " R\n";
                else 
                    out << " L\n";
                out << setprecision(6) << fixed;
                out << "    X          " << setw(10) << x.x << setw(10) << x.y << setw(10) << x.z << "\n"
                    << "    Y          " << setw(10) << y.x << setw(10) << y.y << setw(10) << y.z << "\n"
                    << "    Z          " << setw(10) << z.x << setw(10) << z.y << setw(10) << z.z << "\n";
                out << "    same chirality: ";
                if (lcs[atomIdx].chirality.empty())
                    out << "achiral\n";
                else
                    out << boolalpha << sameChirality << "\n";

                //printHcParameters(out, bankParameters[assignedTypeIdx[atomIdx]]);
                //out << "\n";
                printHcParameters(out, multipoleModelPalameters, atomIdx);
                //out << "\n";
            }
            
        }

        out.close();

    }


    void HcAtomBankStructureFactorCalculator2::setAnomalous(
        const std::vector<std::complex<double> > & anomalous) 
    {
        int nAtomsExtended = mExtended2normal.size();
        vector<complex<double> > anomalousExtended(nAtomsExtended);
        for (int atomIdx = 0; atomIdx < nAtomsExtended; atomIdx++)
            anomalousExtended[atomIdx] = anomalous[mExtended2normal[atomIdx]];
        mCalculator->setAnomalous(anomalousExtended);
    }
 
    //std::vector<bool> mExtendedCrystalAtomsContribution;
    void HcAtomBankStructureFactorCalculator2::updateExtendedCrystalAtomsContribution(
        const std::vector<bool>& countAtomContribution)
        const
    {
        int nAtomsExtended = mExtendedCrystal.atoms.size();
        mExtendedCrystalAtomsContribution.resize(nAtomsExtended);
        for (int i = 0; i < mExtended2normal.size(); i++)
            mExtendedCrystalAtomsContribution[i] = countAtomContribution[mExtended2normal[i]];
    }


    void HcAtomBankStructureFactorCalculator2::updateExtendedCrystalAtoms(
        const std::vector<AtomInCrystal>& atoms)
    {
        for (int i = 0; i < mExtended2normal.size(); i++)
            mExtendedCrystal.atoms[i] = atoms[mExtended2normal[i]];
        int nAtoms = atoms.size();
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            for (int instanceIdx = 0; instanceIdx < mNormal2extended[atomIdx].size(); instanceIdx++)
            {
                int extendedIdx = mNormal2extended[atomIdx][instanceIdx];
                mExtendedCrystal.atoms[extendedIdx].occupancy *= mAtomContributionWeights[atomIdx][instanceIdx];
            }
        }
    }

    void HcAtomBankStructureFactorCalculator2::calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution) 
    {
        updateExtendedCrystalAtoms(atoms);
        updateExtendedCrystalAtomsContribution(countAtomContribution);
        mCalculator->calculateStructureFactorsAndDerivatives(
            mExtendedCrystal.atoms, hkl, f, mExtended_dTarget_dparam, dTarget_df, mExtendedCrystalAtomsContribution);

        mergeExtendedDerivatives(dTarget_dparam, mExtended_dTarget_dparam);
    }

    void HcAtomBankStructureFactorCalculator2::mergeExtendedDerivatives(
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<TargetFunctionAtomicParamDerivatives>& extended_dTarget_dparam)
        const
    {
        int nAtoms = mNormal2extended.size();
        dTarget_dparam.resize(nAtoms);
        int instanceIdx = 0;
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            dTarget_dparam[atomIdx] = extended_dTarget_dparam[mNormal2extended[atomIdx][instanceIdx]];
            dTarget_dparam[atomIdx].occupancy_derivatives *= mAtomContributionWeights[atomIdx][instanceIdx];
        }

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            for (instanceIdx = 1; instanceIdx < mNormal2extended[atomIdx].size(); instanceIdx++)
            {
                int extendedIdx = mNormal2extended[atomIdx][instanceIdx];
                for (int k = 0; k < dTarget_dparam[atomIdx].adp_derivatives.size(); k++)
                    dTarget_dparam[atomIdx].adp_derivatives[k] += 
                        extended_dTarget_dparam[extendedIdx].adp_derivatives[k];

                dTarget_dparam[atomIdx].atomic_position_derivatives += 
                    extended_dTarget_dparam[extendedIdx].atomic_position_derivatives;

                dTarget_dparam[atomIdx].occupancy_derivatives += 
                    extended_dTarget_dparam[extendedIdx].occupancy_derivatives *
                    mAtomContributionWeights[atomIdx][instanceIdx];
            }
            

    }


    void HcAtomBankStructureFactorCalculator2::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution,
        const DerivativesSelector& selector)
    {
        updateExtendedCrystalAtoms(atoms);
        updateExtendedCrystalAtomsContribution(countAtomContribution);
        mCalculator->calculateStructureFactorsAndDerivatives(
            mExtendedCrystal.atoms, hkl, f, mExtended_dTarget_dparam, dTarget_df, mExtendedCrystalAtomsContribution, selector);
        mergeExtendedDerivatives(dTarget_dparam, mExtended_dTarget_dparam);
    }


	void HcAtomBankStructureFactorCalculator2::calculateFormFactors(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
        updateExtendedCrystalAtomsContribution(includeAtom);
        mExtendedFormFactors1Hkl.resize(mExtended2normal.size());
		mCalculator->calculateFormFactors(hkl, mExtendedFormFactors1Hkl, mExtendedCrystalAtomsContribution);
        mergeExtendedFormFactors(mExtendedFormFactors1Hkl, includeAtom, formFactors);
	}

    void HcAtomBankStructureFactorCalculator2::calculateFormFactorsCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        updateExtendedCrystalAtomsContribution(includeAtom);
        mExtendedFormFactors1Hkl.resize(mExtended2normal.size());
        mCalculator->calculateFormFactorsCart(hkl, mExtendedFormFactors1Hkl, mExtendedCrystalAtomsContribution);
        mergeExtendedFormFactors(mExtendedFormFactors1Hkl, includeAtom, formFactors);
    }


    void HcAtomBankStructureFactorCalculator2::mergeExtendedFormFactors(
        const std::vector<std::complex<double> >& extended,
        const std::vector<bool>& includeAtom,// regular (not extended)
        std::vector < std::complex<double> >& merged)
        const
    {
        int atomIdx, nAtoms = mNormal2extended.size();      
        merged.resize(nAtoms);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (includeAtom[atomIdx])
            {
                int nInstances = mNormal2extended[atomIdx].size();
                merged[atomIdx] = 0.0;
                for (int instanceIdx = 0; instanceIdx < mNormal2extended[atomIdx].size(); instanceIdx++)
                    merged[atomIdx] += extended[mNormal2extended[atomIdx][instanceIdx]] * mAtomContributionWeights[atomIdx][instanceIdx];
            }

    }


    void HcAtomBankStructureFactorCalculator2::mergeExtendedFormFactors(
        const std::vector < std::vector<std::complex<double> > >& extended,
        const std::vector<bool>& includeAtom,
        std::vector < std::vector<std::complex<double> > >& merged)
        const
    {
        int hklIdx, nHkl = extended.size();
        merged.resize(nHkl);      

        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
            mergeExtendedFormFactors(extended[hklIdx], includeAtom, merged[hklIdx]);
        
    }


    void HcAtomBankStructureFactorCalculator2::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        updateExtendedCrystalAtomsContribution(includeAtom);


        vector< vector<complex<double> > > formFactorsExtended;

        if (mN_Threads == 1)
        {
            formFactorsExtended.resize(hkl.size());
            for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
                mCalculator->calculateFormFactors(hkl[hklIdx], formFactorsExtended[hklIdx], mExtendedCrystalAtomsContribution);
            //SfCalculator::calculateFormFactors(hkl, formFactorsExtended, mExtendedCrystalAtomsContribution);
        }
        else
        {

            vector< vector<Vector3i> > hklsPerThread(mN_Threads);
            vector < vector< vector<complex<double> > > > formFactorsPerThread(mN_Threads);

            int nHkl = hkl.size();
            int reminder = nHkl % mN_Threads;
            int quotient = nHkl / mN_Threads;

            for (int i = 0; i < mN_Threads; i++)
            {
                int start = i * quotient + min(i, reminder);
                int end = (i + 1) * quotient + min(i + 1, reminder);
                hklsPerThread[i].insert(hklsPerThread[i].begin(), hkl.begin() + start, hkl.begin() + end);
            }

            omp_set_num_threads(static_cast<int>(mN_Threads));
#pragma omp parallel
            {

                int id = omp_get_thread_num();
                cout << "thread " << id << "\n";
                int nHkl = hklsPerThread[id].size();
                formFactorsPerThread[id].resize(nHkl);
                for (int i = 0; i < nHkl; i++)
                    mCalculators[id]->calculateFormFactors(
                        hklsPerThread[id][i],
                        formFactorsPerThread[id][i], mExtendedCrystalAtomsContribution);
            }

            formFactorsExtended.clear();
            for (int i = 0; i < mN_Threads; i++)
                formFactorsExtended.insert(formFactorsExtended.end(), formFactorsPerThread[i].begin(), formFactorsPerThread[i].end());
        }
        mergeExtendedFormFactors(formFactorsExtended, includeAtom, formFactors);

    }



        void HcAtomBankStructureFactorCalculator2::calculateStructureFactors(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            const std::vector<bool> &countAtomContribution) 
        {
            updateExtendedCrystalAtoms(atoms);
            updateExtendedCrystalAtomsContribution(countAtomContribution);
            mCalculator->calculateStructureFactors(mExtendedCrystal.atoms, hkl, f, mExtendedCrystalAtomsContribution);
        }

        void HcAtomBankStructureFactorCalculator2::update(
            const std::vector<AtomInCrystal> &atoms) 
        {
            updateExtendedCrystalAtoms(atoms);
            mCalculator->update(mExtendedCrystal.atoms);
        }

        void HcAtomBankStructureFactorCalculator2::calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &scatteringFactor,
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution) 
        {
            updateExtendedCrystalAtomsContribution(countAtomContribution);
            mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, mExtendedDerivatives1hkl, mExtendedCrystalAtomsContribution);
            /*merge the derivatives*/
            int nAtoms = mNormal2extended.size();
            derivatives.adpDerivatives.resize(nAtoms);
            derivatives.anomalousScatteringDerivatives.resize(nAtoms);
            derivatives.atomicPostionDerivatives.resize(nAtoms);
            derivatives.occupancyDerivatives.resize(nAtoms);

            for(int atomIdx=0; atomIdx<nAtoms; atomIdx++)
                if (countAtomContribution[atomIdx])
                {
                    int nAdps = mCrystal.atoms[atomIdx].adp.size();
                    derivatives.adpDerivatives[atomIdx].assign(nAdps, 0.0);
                    derivatives.anomalousScatteringDerivatives[atomIdx] = 0.0;
                    derivatives.atomicPostionDerivatives[atomIdx] = Vector3d(0.0, 0.0, 0.0);
                    derivatives.occupancyDerivatives[atomIdx] = 0.0;

                    int nInstances = mNormal2extended[atomIdx].size();
                    for (int instanceIdx = 0; instanceIdx < nInstances; instanceIdx++)
                    {
                        int idxExtended = mNormal2extended[atomIdx][instanceIdx];
                        for (int adpIdx = 0; adpIdx < nAdps; adpIdx++)
                            derivatives.adpDerivatives[atomIdx][adpIdx] += mExtendedDerivatives1hkl.adpDerivatives[idxExtended][adpIdx];

                        derivatives.atomicPostionDerivatives[atomIdx] += mExtendedDerivatives1hkl.atomicPostionDerivatives[idxExtended];
                        derivatives.occupancyDerivatives[atomIdx] += mExtendedDerivatives1hkl.occupancyDerivatives[idxExtended] * mAtomContributionWeights[atomIdx][instanceIdx];
                    }
                }
        }


}
 