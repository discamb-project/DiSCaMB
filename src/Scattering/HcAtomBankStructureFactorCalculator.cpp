#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"


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

    

    HcAtomBankStructureFactorCalculator::HcAtomBankStructureFactorCalculator(
        const Crystal &crystal,
        const std::vector<AtomType> &atomTypes,
        const std::vector<AtomTypeHC_Parameters> &parameters,
        bool electronScattering,
        const DescriptorsSettings &settings,
        const std::string &assignemntInfoFile,
        const std::string& parametersInfoFile,
        const std::string& multipolarCif,
        int nThreads,
        double unitCellCharge,
        bool scaleToMatchCharge,
        const string& iamTable,
        bool iamElectronScattering,
        bool frozen_lcs)
    {
        set(crystal, atomTypes, parameters, electronScattering, settings, assignemntInfoFile, parametersInfoFile, multipolarCif, nThreads, unitCellCharge, scaleToMatchCharge, iamTable, iamElectronScattering, frozen_lcs);

    }


    HcAtomBankStructureFactorCalculator::HcAtomBankStructureFactorCalculator(
        const Crystal &crystal, 
        const nlohmann::json &data)
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

        int nCores=1;
        if (data.find("n cores") != data.end())
            nCores = data.find("n cores")->get<int>();

        string iamTable;
        if (data.find("table") != data.end())
            iamTable = data.find("table")->get<string>();
        
        bool iamElectronScattering = false;
        if (data.find("iam electron scattering") != data.end())
            iamElectronScattering = data.find("iam electron scattering")->get<bool>();
        bool frozen_lcs = data.value("frozen lcs", false);

        MATTS_BankReader bankReader;
        vector<AtomType> types;
        vector<AtomTypeHC_Parameters> hcParameters;
        BankSettings bankSettings;


		if(bankPath.empty())
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
            {
                on_error::throwException("no bank file specified or found in working directory", __FILE__, __LINE__);
                //string bankString;
                //stringstream bankStream;
                //default_ubdb_bank_string(bankString);
                //bankStream << bankString;
                //bankReader.read(bankStream, types, hcParameters, bankSettings, true);
            }
            else
            {
                bankReader.read(*bnkFiles.begin(), types, hcParameters, bankSettings, true);
            }
		}
		else
			bankReader.read(bankPath, types, hcParameters, bankSettings, true);
        set(crystal,types, hcParameters, electronScattering, DescriptorsSettings(), assignmentInfoFile,
            parametersInfoFile, multipolarCif, nCores, unitCellCharge, scaleHcParameters, iamTable, iamElectronScattering, frozen_lcs);
    }


    HcAtomBankStructureFactorCalculator::~HcAtomBankStructureFactorCalculator()
    {
    }
    
    void HcAtomBankStructureFactorCalculator::getModelInformation(std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        modelInfo = mModelInfo;
    }

    void HcAtomBankStructureFactorCalculator::set(
        const Crystal &crystal,
        const std::vector<AtomType> &atomTypes,
        const std::vector<AtomTypeHC_Parameters> &bankParameters,
        bool electronScattering,
        const DescriptorsSettings &settings,
        const std::string &assignemntInfoFile,
        const std::string& parametersInfoFile,
        const std::string& multipolarCif,
        int nThreads,
        double unitCellCharge,
        bool scaleToMatchCharge,
        const string& iamTable,
        bool iamElectronScattering,
        bool frozen_lcs)
    {
        mModelInfo.clear();

        mModelInfo.push_back({ "SCATTERING MODEL", "TAAM" });
        mModelInfo.push_back({ "TAAM DATABANK", "MATTS (Kumar, P., Gruza, B., Bojarowski, S. A. & Dominiak, P. M. (2019). Acta Cryst. A75, 398-408.)" });


        mN_Threads = nThreads;

        CrystalAtomTypeAssigner assigner;
        assigner.setAtomTypes(atomTypes);
        assigner.setDescriptorsSettings(settings);
        vector < LocalCoordinateSystem<AtomInCrystalID> > lcs;
        vector<int> types;
        assigner.assign(crystal, types, lcs);
        if (!assignemntInfoFile.empty())
        {
            if (assignemntInfoFile == string("print_to_discamb2tsc_log_file") || assignemntInfoFile == string("print_to_discambMATTS2tsc_log_file"))
            {
                vector<string> unassignedAtoms;
                for (int i = 0; i < types.size(); i++)
                {
                    if (types[i] < 0)
                        unassignedAtoms.push_back(crystal.atoms[i].label);
                    else
                    {
                        string typeLabel = atomTypes[types[i]].id;
                        int labelSize = typeLabel.size();
                        if (labelSize > 3)
                        {
                            string back = typeLabel.substr(labelSize - 3);
                            string front = typeLabel.substr(0, labelSize - 3);
                            
                            if ( back == string("000"))
                            {
                                bool hasDigit = false;
                                for (auto& c : front)
                                    if (isdigit(c))
                                        hasDigit = true;
                                if(!hasDigit)
                                    unassignedAtoms.push_back(crystal.atoms[i].label);
                            }
                        }
                    }
                }

                //ofstream out("discamb2tsc.log", std::ofstream::out | std::ofstream::app);
                //clog << "atomic form factor for X-ray scattering calculated with\nHansen-Coppens model parameterized with MATTS databank.\n"
                //    << "atom type assigned to " << crystal.atoms.size() - unassignedAtoms.size() << " of " << crystal.atoms.size() << " atoms\n";
                
                assigner.printAssignment(clog, crystal, types, lcs);

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
                assigner.printAssignment(out, crystal, types, lcs);
                out.close();
            }
        }

        HC_ModelParameters multipoleModelPalameters;


        vector<int> atomicNumbers;
        vector<int> nonMultipolarAtoms;
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);

        vector<double> multiplicityTimesOccupancy;
        for (auto& atom : crystal.atoms)
            multiplicityTimesOccupancy.push_back(atom.multiplicity * atom.occupancy);
        
        if (scaleToMatchCharge)
            taam_utilities::type_assignment_to_HC_parameters(
                bankParameters, types, multiplicityTimesOccupancy, atomicNumbers, unitCellCharge,
                multipoleModelPalameters, true, nonMultipolarAtoms);
        else
            taam_utilities::type_assignment_to_unscaled_HC_parameters(bankParameters, types, atomicNumbers,
                multipoleModelPalameters, true, nonMultipolarAtoms);

        if (!multipolarCif.empty())
            cif_io::saveCif(multipolarCif, crystal, multipoleModelPalameters, lcs);


        //ubdb_utilities::ubdb_type_assignment_to_HC_parameters(
        //    bankParameters, types, atomicNumbers,
        //    multipoleModelPalameters, true, nonMultipolarAtoms);

        vector<shared_ptr<LocalCoordinateSystemInCrystal> > lcaCalculators;
        for (int lcsIdx = 0; lcsIdx < lcs.size(); lcsIdx++)
            lcaCalculators.push_back(
                shared_ptr<LocalCoordinateSystemInCrystal>(
                    types[lcsIdx] >= 0 ?
                    new LocalCoordinateSystemCalculator(lcs[lcsIdx], crystal) :
                    new LocalCoordinateSystemCalculator()
                )
            );
        if (!parametersInfoFile.empty())
            printAssignedMultipolarParameters(
                parametersInfoFile,
                crystal,
                lcs,
                atomTypes,
                bankParameters,
                multipoleModelPalameters,
                types);

        mHcCalculator = shared_ptr<SfCalculator>(
                            new AnyHcCalculator(crystal, multipoleModelPalameters, lcaCalculators, electronScattering, false, false, 1, frozen_lcs));
        mIamCalculator = shared_ptr<SfCalculator>(
            new AnyIamCalculator(crystal, iamElectronScattering, iamTable));
        vector< std::shared_ptr<SfCalculator> > calculators{ mHcCalculator, mIamCalculator };
        vector<vector<int> > calculatorAtoms(2);
        for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
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
                new AnyHcCalculator(crystal, multipoleModelPalameters, lcaCalculators, electronScattering, false, false, 1, frozen_lcs)));
            mIamCalculators.push_back(shared_ptr<SfCalculator>(
                new AnyIamCalculator(crystal, iamElectronScattering, iamTable)));

            vector< std::shared_ptr<SfCalculator> > calculators2{ mHcCalculators.back(), mIamCalculators.back() };
            mCalculators.push_back(shared_ptr< CombinedStructureFactorCalculator>(new CombinedStructureFactorCalculator(calculators2, calculatorAtoms)));
        }


    }

    //void HcAtomBankStructureFactorCalculator::printHcParameters(
    //    std::ofstream& out,
    //    const AtomTypeHC_Parameters& bankParameters)
    //{
    //    out << "kappa " << bankParameters.kappa << "\n"
    //        << "kappa_prime " << bankParameters.kappa_prime << "\n"
    //        << "P_val " << bankParameters.p_val << "\n"
    //        << "number of Plm terms " << bankParameters.p_lms.size() << "\n"
    //        << "Plms:\n";

    //    for (int i = 0; i < bankParameters.p_lms.size(); i++)
    //        out << bankParameters.p_lm_indices[i].first << " "
    //            << bankParameters.p_lm_indices[i].second << " " << bankParameters.p_lms[i] << "\n";
    //}

    void HcAtomBankStructureFactorCalculator::printHcParameters(
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

    void HcAtomBankStructureFactorCalculator::printAssignedMultipolarParameters(
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


    void HcAtomBankStructureFactorCalculator::setAnomalous(
        const std::vector<std::complex<double> > & anomalous) 
    {
        mCalculator->setAnomalous(anomalous);
    }
    
    void HcAtomBankStructureFactorCalculator::calculateStructureFactorsAndDerivatives(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dparam,
            const std::vector<std::complex<double> > &dTarget_df,
            const std::vector<bool> &countAtomContribution) 
    {
        mCalculator->calculateStructureFactorsAndDerivatives(
            atoms, hkl, f, dTarget_dparam, dTarget_df, countAtomContribution);
    }

	void HcAtomBankStructureFactorCalculator::calculateFormFactors(
		const Vector3i& hkl,
		std::vector<std::complex<double> >& formFactors,
		const std::vector<bool>& includeAtom)
		const
	{
		mCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
	}

    void HcAtomBankStructureFactorCalculator::calculateFormFactorsCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mCalculator->calculateFormFactorsCart(hkl, formFactors, includeAtom);
    }

    //void HcAtomBankStructureFactorCalculator::calculateFormFactors(
    //    const std::vector<Vector3d>& hkl,
    //    std::vector< std::vector<std::complex<double> > >& formFactors,
    //    const std::vector<bool>& includeAtom) const
    //{

    //}

    void HcAtomBankStructureFactorCalculator::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        if (mN_Threads == 1)
        {
            SfCalculator::calculateFormFactors(hkl, formFactors, includeAtom);
            return;
        }

        vector< vector<Vector3i> > hklsPerThread(mN_Threads);
        vector < vector< vector<complex<double> > > > formFactorsPerThread(mN_Threads);

        int nHkl = hkl.size();
        int reminder = nHkl % mN_Threads;
        int quotient = nHkl / mN_Threads;

        for (int i = 0; i < mN_Threads; i++)
        {
            int start = i * quotient + min(i, reminder);
            int end = (i + 1) * quotient + min(i+1, reminder);
            hklsPerThread[i].insert(hklsPerThread[i].begin(), hkl.begin() + start, hkl.begin() + end);
        }

        omp_set_num_threads(static_cast<int>(mN_Threads));
#pragma omp parallel
        {
            
            int id = omp_get_thread_num();
            cout << "thread " << id << "\n";
            int nHkl = hklsPerThread[id].size();
            formFactorsPerThread[id].resize(nHkl);
            for(int i=0; i<nHkl; i++)
                mCalculators[id]->calculateFormFactors(
                    hklsPerThread[id][i], 
                    formFactorsPerThread[id][i], includeAtom);
        }

        formFactors.clear();
        for (int i = 0; i < mN_Threads; i++)
            formFactors.insert(formFactors.end(), formFactorsPerThread[i].begin(), formFactorsPerThread[i].end());

    }



        void HcAtomBankStructureFactorCalculator::calculateStructureFactors(
            const std::vector<AtomInCrystal> &atoms,
            const std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &f,
            const std::vector<bool> &countAtomContribution) 
        {
            mCalculator->calculateStructureFactors(atoms, hkl, f, countAtomContribution);
        }

        void HcAtomBankStructureFactorCalculator::update(
            const std::vector<AtomInCrystal> &atoms) 
        {
            mCalculator->update(atoms);
        }

        void HcAtomBankStructureFactorCalculator::calculateStructureFactorsAndDerivatives(
            const Vector3i &hkl,
            std::complex<double> &scatteringFactor,
            discamb::SfDerivativesAtHkl &derivatives,
            const std::vector<bool> &countAtomContribution) 
        {
            mCalculator->calculateStructureFactorsAndDerivatives(hkl, scatteringFactor, derivatives, countAtomContribution);
        }


}
 