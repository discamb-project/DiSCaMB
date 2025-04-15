#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/cif_io.h"
#include "discamb/IO/discamb_io.h"
#include "discamb/IO/fragmentation_io.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamSfCalculator.h"
#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/statistics.h"
#include "discamb/Scattering/taam_utilities.h"
#include "discamb/Scattering/TscFileBasedSfCalculator.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include <fstream>

using namespace std;
using namespace discamb;

void test_taam_disorder(int argc, char *argv[])
{
    if (argc != 6)
    {
        cout<< "expected 5 arguments:\n"
            " (1) cif/res file\n"
            " (2) bank path\n"
            " (3) fragment definitions\n"
            " (4) representatives file\n"
            " (5) output file\n";
        exit(0);
    }
    string structureFile = argv[1];
    string bankPath = argv[2];
    string fragmentDefinitionsFile = argv[3];
    string representativesFile = argv[4];
    string outputFile = argv[5];

    // structure

    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);

    // atom types

    MATTS_BankReader bankReader;
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> typeParameters;
    BankSettings bankSettings;

    bankReader.read(bankPath, atomTypes, typeParameters, bankSettings);

    // read fragments
    vector<FragmentConstructionData> fragmentConstructionData;
    fragmentation_io::readPlainText(fragmentDefinitionsFile, fragmentConstructionData);
    vector<QmFragmentInCrystal> qmFragments;

    fragmentation::make_qm_fragments(crystal, fragmentConstructionData, qmFragments);

    // convert atoms representation

    vector< vector<AtomInCrystalID> > fragments(qmFragments.size());
    vector<vector<pair<string, string> > > fragmentsStrStr;
    for (int i = 0; i < fragments.size(); i++)
    {
        crystal_structure_utilities::convertAtomList(crystal, qmFragments[i].atoms.atomList, fragments[i]);
        fragmentsStrStr.push_back(qmFragments[i].atoms.atomList);
    }

    // read representatives
    
    vector<vector<AtomRepresentativeInfo> > representatives;
    vector<string> fragmentLabels;
    for (auto& fragment : qmFragments)
        fragmentLabels.push_back(fragment.label);
    discamb_io::read_representatives(representativesFile, crystal, fragmentLabels, fragmentsStrStr, representatives);

     
    // assign atom types
    CrystalAtomTypeAssigner assigner;
    assigner.setAtomTypes(atomTypes);
    assigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
    vector<vector<int> > typeIds;
    vector < vector <LocalCoordinateSystem<AtomInCrystalID> > > lcs;
    assigner.assign(crystal, fragments, representatives, typeIds, lcs);

    // print assignement

    ofstream out(outputFile);
    if (!out.good())
        on_error::throwException("cannot write output file", __FILE__, __LINE__);

    int atomIdx, fragmentIdx, nFragments = fragments.size();
    for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
    {
        int nAtoms = fragments[fragmentIdx].size();
        out << "\nFragment " << fragmentLabels[fragmentIdx] << "\n\n";
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            out << setw(25) << left << fragmentsStrStr[fragmentIdx][atomIdx].first + " " + fragmentsStrStr[fragmentIdx][atomIdx].second;
            if (typeIds[fragmentIdx][atomIdx] < 0)
                out << "---";
            else
                out << atomTypes[typeIds[fragmentIdx][atomIdx]].id;
            out << "\n";
        }
        
    }
    out.close();
}

double agreementFactorAbs(
    const vector<complex<double> >& data1,
    const vector<complex<double> >& data2)
{
    vector<double> abs_sf1, abs_sf2;

    double s12 = 0, s22 = 0;
    int n = data1.size();
    
    for (int i = 0; i < n; i++)
    {
        abs_sf1.push_back(abs(data1[i]));
        abs_sf2.push_back(abs(data2[i]));
        s12 += abs_sf1[i] * abs_sf2[i];
        s22 += abs_sf2[i] * abs_sf2[i];
    }

    double scale = s12 / s22;
    double num = 0, den = 0;

    for (int i = 0; i < n; i++)
    {
        num += abs(abs_sf1[i] - scale * abs_sf2[i]);
        den += abs(abs_sf1[i] + scale * abs_sf2[i]);
    }

    return 200 * num / den;

}

/*

    f_calc_a = flex.abs(f_calc_a)
    f_calc_b = flex.abs(f_calc_b)
    scale = flex.sum(f_calc_a * f_calc_b) / flex.sum(f_calc_b * f_calc_b)
    num = flex.sum(flex.abs(f_calc_a - scale * f_calc_b))
    den = flex.sum(flex.abs(f_calc_a + scale * f_calc_b))
    diff = flex.abs(f_calc_a - f_calc_b)
*/

double agreementFactor(
    const vector<complex<double> >& data1,
    const vector<complex<double> >& data2)
{
    double abs_sf1, abs_sf2;

    double s12 = 0, s22 = 0;
    int n = data1.size();

    for (int i = 0; i < n; i++)
    {
        abs_sf1 = abs(data1[i]);
        abs_sf2 = abs(data2[i]);
        s12 += abs_sf1 * abs_sf2;
        s22 += abs_sf2 * abs_sf2;
    }

    double scale = s12 / s22;
    double num = 0, den = 0;

    for (int i = 0; i < n; i++)
    {
        num += abs(data1[i] - scale * data2[i]);
        den += abs(data2[i] + scale * data2[i]);
    }

    return 200 * num / den;

}



//######

bool read_bank_and_assign_atoms(
    const string bank_filepath,
    const Crystal& crystal,
    vector<AtomType>& atomTypes,
    vector<AtomTypeHC_Parameters>& hcParameters,
    BankSettings& bankSettings,
    CrystalAtomTypeAssigner& crystalAssigner,
    int& nAtomsInAsymmetricUnit,
    StructureWithDescriptors& structure,
    vector<int>& typeIds,
    vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs,
    vector<string>& lcs_strings
) {

    // read bank parameters
    MATTS_BankReader bankReader;

    ifstream bankStream{ bank_filepath };
    bankReader.read(bankStream, atomTypes, hcParameters, bankSettings, true);

    crystalAssigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
    crystalAssigner.setAtomTypes(atomTypes);

    vector<int> atomicNumbers;
    vector<string> labels;
    bool processingCrystalStructure = true;

    nAtomsInAsymmetricUnit = crystal.atoms.size();

    crystalAssigner.assign(crystal, typeIds, lcs, structure);

    for (auto& atom : crystal.atoms)
        labels.push_back(atom.label);

    for (auto& coordinateSystem : lcs)
        lcs_strings.push_back(ubdbLcsAsString(coordinateSystem, labels));

    return true;
}


void TAAM_types(
    Crystal& crystal,
    string bank_filepath) {

    // read bank parameters
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> hcParameters;
    BankSettings bankSettings;
    CrystalAtomTypeAssigner crystalAssigner;
    vector < LocalCoordinateSystem<AtomInCrystalID> > lcs;
    vector<int> typeIds;
    int n_atoms;
    StructureWithDescriptors structure;
    vector<string> lcs_strings;

    //
    MATTS_BankReader bankReader;

    ifstream bankStream{ bank_filepath };
    bankReader.read(bankStream, atomTypes, hcParameters, bankSettings, true);

    crystalAssigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
    crystalAssigner.setAtomTypes(atomTypes);

    vector<int> atomicNumbers;
    vector<string> labels;
    bool processingCrystalStructure = true;

    int nAtomsInAsymmetricUnit = crystal.atoms.size();

    crystalAssigner.assign(crystal, typeIds, lcs, structure);

    for (auto& atom : crystal.atoms)
        labels.push_back(atom.label);

    for (auto& coordinateSystem : lcs)
        lcs_strings.push_back(ubdbLcsAsString(coordinateSystem, labels));

    //

    //read_bank_and_assign_atoms(
    //    bank_filepath,
    //    crystal,
    //    atomTypes,
    //    hcParameters,
    //    bankSettings,
    //    crystalAssigner,
    //    n_atoms,
    //    structure,
    //    typeIds,
    //    lcs,
    //    lcs_strings
    //);

    //if (log_assignment) {
    //    write_assignment_logs(
    //        crystal,
    //        atomTypes,
    //        hcParameters,
    //        bankSettings,
    //        crystalAssigner,
    //        n_atoms,
    //        structure,
    //        typeIds,
    //        lcs,
    //        lcs_strings
    //    );
    //}

    // get TAAM parameters with unit cell charge scaled/shifted to 0
    HC_ModelParameters multipoleModelPalameters;

    
    vector<int> nonMultipolarAtoms;
    crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);

    vector<double> multiplicityTimesOccupancy;
    for (auto& atom : crystal.atoms)
        multiplicityTimesOccupancy.push_back(atom.multiplicity * atom.occupancy);
    double unitCellCharge = 0.0;
    taam_utilities::TaamAtomicChargeInfo chargeInfo;
    taam_utilities::type_assignment_to_HC_parameters(
        hcParameters,
        typeIds,
        multiplicityTimesOccupancy,
        atomicNumbers,
        unitCellCharge,
        multipoleModelPalameters,
        true,
        nonMultipolarAtoms,
        chargeInfo
    );

    for (int i = 0; i < crystal.atoms.size(); i++)
        cout << crystal.atoms[i].label << " " << chargeInfo.atomicChargesBeforeScaling[i] << " " << chargeInfo.atomicChargesAfterScaling[i] << endl;
    vector<shared_ptr<LocalCoordinateSystemInCrystal> > lcaCalculators;
    for (auto coordinateSystem : lcs)
        lcaCalculators.push_back(
            shared_ptr<LocalCoordinateSystemInCrystal>(
                new LocalCoordinateSystemCalculator(coordinateSystem, crystal)));

    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator =
        shared_ptr<AtomicFormFactorCalculationsManager>(
            new HcFormFactorCalculationsManager(crystal, multipoleModelPalameters, lcaCalculators));

    ////calculate electron structure factors 
    //if (is_electron) {
    //    vector<double> nuclearCharge;
    //    for (int z : atomicNumbers)
    //        nuclearCharge.push_back(z);
    //    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
    //        new ElectronFromXrayFormFactorCalculationManager(
    //            crystal.unitCell,
    //            nuclearCharge,
    //            formfactorCalculator));
    //}

    
}

//######

void compare_sf_taam_qm(
    int argc, char* argv[])
{
    vector<string> arguments, options;
    map<string, string> optionsWithValues;
    parse_cmd::get_args_and_options(argc, argv, arguments, options, optionsWithValues);


    if (arguments.size() != 3)
    {
        cout << "expected 3 arguments:\n"
            " (1) crystal structure cif/res file\n"
            " (2) tsc file for model 1\n"
            " (3) tsc file for model 2\n";
        exit(0);
    }

    string crystalFile = arguments[0];
    string tscFile1 = arguments[1];
    string tscFile2 = arguments[2];
    bool electrons = parse_cmd::hasOption(options, "-e");
    bool no000 = parse_cmd::hasOption(options, "-no000");

    // tsc
    vector<string> atomLabels;
    vector<Vector3i> hkl;
    vector<vector<complex<double> > > formFactors1, formFactors2;
    tsc_io::read_tsc(tscFile1, atomLabels, hkl, formFactors1);
    tsc_io::read_tsc(tscFile2, atomLabels, hkl, formFactors2);

    if (no000)
    {
        auto it = find(hkl.begin(), hkl.end(), Vector3i(0, 0, 0));
        if (it != hkl.end())
        {
            int hkl000_idx = distance(hkl.begin(), it);
            hkl.erase(hkl.begin() + hkl000_idx);
            formFactors1.erase(formFactors1.begin() + hkl000_idx);
            formFactors2.erase(formFactors2.begin() + hkl000_idx);
        }
    }
    

    
    // crystal

    Crystal crystal;
    structure_io::read_structure(crystalFile, crystal);
    
    for (auto& atom : crystal.atoms)
        atom.adp.clear();

    //#########################
    // sf calculation


    AnyScattererStructureFactorCalculator sfCalculator(crystal);
    vector<complex<double> > structureFactors1, structureFactors2;
    map < Vector3i, vector<complex<double> > > formFactors;
    //shared_ptr<AtomicFormFactorCalculationsManager> calcManager;
    //
    //// for model 1

    //for (int i = 0; i < hkl.size(); i++)
    //    formFactors[hkl[i]] = formFactors1[i];
    //calcManager = std::make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell,formFactors);
    //    
    //sfCalculator.setAtomicFormfactorManager(calcManager);
    //sfCalculator.calculateStructureFactors(hkl, structureFactors1);

    //// for model 2

    //for (int i = 0; i < hkl.size(); i++)
    //    formFactors[hkl[i]] = formFactors2[i];
    //calcManager = std::make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, formFactors);

    //sfCalculator.setAtomicFormfactorManager(calcManager);
    //sfCalculator.calculateStructureFactors(hkl, structureFactors2);

    // calculation 2
    nlohmann::json data;
    if(electrons)
        data["x2e"] = true;
    data["tsc file"] = tscFile1;
    TscFileBasedSfCalculator tscCalc1(crystal, data);
    data["tsc file"] = tscFile2;
    TscFileBasedSfCalculator tscCalc2(crystal, data);
    vector<bool> countAtomContribution(crystal.atoms.size(), true);
    tscCalc1.calculateStructureFactors(crystal.atoms, hkl, structureFactors1, countAtomContribution);
    tscCalc2.calculateStructureFactors(crystal.atoms, hkl, structureFactors2, countAtomContribution);

    //#########################
    // compare structure factors

    cout << "diff = " << agreementFactor(structureFactors1, structureFactors2) << endl;

    cout << "atom by atom :\n";
    int atomIdx, nAtoms = atomLabels.size();
    int nHkl = hkl.size();
    vector<complex<double> > ff1(nHkl), ff2(nHkl);
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            ff1[hklIdx] = formFactors1[hklIdx][atomIdx];
            ff2[hklIdx] = formFactors2[hklIdx][atomIdx];
        }

        cout << atomLabels[atomIdx] << " " << agreementFactor(ff1, ff2) << endl;
        
    }
}

void nGaussian(const string &fName)
{
    Crystal crystal;
    structure_io::read_structure(fName, crystal);

    vector<double> array_of_a = { 0.0893, 0.2563, 0.7570, 1.0487, 0.3575  };
    vector<double> array_of_b = { 0.2465, 1.7100, 6.4094, 18.6113, 50.2523 };
    map<string, NGaussianFormFactor> scatterers;
    for (int i = 1; i < 100; i++)
        scatterers[periodic_table::symbol(i)] = NGaussianFormFactor(periodic_table::symbol(i), array_of_a, array_of_b, 0.0);
    shared_ptr<AtomicFormFactorCalculationsManager> formfactorCalculator;
    formfactorCalculator = shared_ptr<AtomicFormFactorCalculationsManager>(
        new IamFormFactorCalculationsManager(crystal, scatterers));
    AnyScattererStructureFactorCalculator calculator(crystal);
    calculator.setAtomicFormfactorManager(formfactorCalculator);
    vector<complex<double>> sf;
    vector<Vector3i> hkl = { {1,0,0},{2,3,4} };
    calculator.calculateStructureFactors(hkl, sf);

}

void check_args(
    int argc,
    char* argv[],
    int nExpected,
    const vector<string>& argsDescription)
{
    if (argc - 1 == nExpected)
        return;
    cout << "expected " << nExpected << " arguments\n:";
    
    for (int i = 0; i < argsDescription.size(); i++)
        cout << "   (" << i + 1 << ")  " << argsDescription[i] << "\n";
    exit(0);
}

void testTypeAssignemntLog(const string& structureFile, const string &bnk)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    nlohmann::json data;
    data["bank path"] = bnk;
    HcAtomBankStructureFactorCalculator calculator(crystal, data);
    vector<Vector3i> hkl{ Vector3i(1,0,0),Vector3i(1,1,0) };
    vector<complex<double> > f;
    vector<bool> countAtomContribution(crystal.atoms.size(), true);
    calculator.calculateStructureFactors(crystal.atoms, hkl, f, countAtomContribution);
    for (int i = 0; i < hkl.size(); i++)
        cout << hkl[i] << " " << f[i] << "\n";
}


void test_int_floor()
{
    for (int i = 0; i < 30; i++)
    {
        double d = i * 1.01;
        int j = floor(d);
        cout << j << "\n";
    }
}

void test_connectivity(const string& structFile)
{
    vector<string> structFiles;
    if (!structFile.empty())
        structFiles.push_back(structFile);
    else
        file_system_utilities::find_files("res", structFiles);
    
    for (auto& fileName : structFiles)
    {
        Crystal crystal;
        structure_io::read_structure(fileName, crystal);
        UnitCellContent uc(crystal);
        vector<vector<UnitCellContent::AtomID> > connectivity, connectivity3;
            
        structural_properties::calcUnitCellConnectivity(uc, connectivity, 0.4, "simple");
        structural_properties::calcUnitCellConnectivity(uc, connectivity3, 0.4, "boxes");
        for (auto& entry : connectivity)
            sort(entry.begin(), entry.end());
        for (auto& entry : connectivity3)
            sort(entry.begin(), entry.end());
        if (connectivity != connectivity3)
        {
            cout << fileName << " connectivity differ\n";
            for(int i=0;i< connectivity.size(); i++)
                if (connectivity[i] != connectivity3[i])
                {
                    string label, symmOp;
                    uc.interpreteAtomID(UnitCellContent::AtomID(i), label, symmOp);
                    cout << "difference for " << i << " " << label << " " << symmOp << "\n";
                    if (connectivity[i].size() != connectivity3[i].size())
                    {
                        cout << "  different size\n";
                        cout << "    connectivity:\n";
                        for (int j = 0; j < connectivity[i].size(); j++)
                        {
                            uc.interpreteAtomID(connectivity[i][j], label, symmOp);
                            cout << "      " << label << " " << symmOp 
                                 << " ID " << connectivity[i][j].atomIndex << " " << connectivity[i][j].unitCellPosition << "\n";
                        }
                        cout << "    connectivity3:\n";
                        for (int j = 0; j < connectivity3[i].size(); j++)
                        {
                            uc.interpreteAtomID(connectivity3[i][j], label, symmOp);
                            cout << "      " << label << " " << symmOp 
                                 << " ID " << connectivity3[i][j].atomIndex << " " << connectivity3[i][j].unitCellPosition << "\n";
                        }

                    }
                    else
                    {
                        cout << "  connectivity:\n";
                        for (int j = 0; j < connectivity[i].size(); j++)
                        {
                            uc.interpreteAtomID(connectivity[i][j], label, symmOp);
                            cout << "    " << label << " " << symmOp << "\n"; 
                        }
                        cout << "  connectivity3:\n";
                        for (int j = 0; j < connectivity3[i].size(); j++)
                        {
                            uc.interpreteAtomID(connectivity3[i][j], label, symmOp);
                            cout << "    " << label << " " << symmOp << "\n";
                        }

                    }
                }
            break;
        }
        else
            cout << fileName << " connectivity match\n";
    }
}

void tsc_differences(
    const string& tsc1,
    const string& tsc2)
{
    vector<string> labels1, labels2;
    vector<Vector3i> hkl;
    vector<vector<complex<double> > > ff1, ff2;
    tsc_io::read_tsc(tsc1, labels1, hkl, ff1);
    int nLabels = labels1.size();
    tsc_io::read_tsc(tsc2, labels2, hkl, ff2);

    if (nLabels != labels2.size())
        cout << "file differ by number of scatterers " << "\n";
    if(labels1!=labels2)
        cout << "file differ by scatterer names " << "\n";

    int nHkl = ff1.size();
    if(nHkl!=ff2.size())
        cout << "file differ by number of hkl " << "\n";
    cout << "number of atoms " << nLabels << "\n";
    cout << "number of hkl " << nHkl << "\n";
    for (int atomIdx = 0; atomIdx < nLabels; atomIdx++)
    {
        //cout << atomIdx << " " << ff1[nHkl - 1][atomIdx] << " " << ff2[nHkl - 1][atomIdx] << "\n";
        bool different_atomic_ff = false;
        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
            if (ff1[hklIdx][atomIdx] != ff2[hklIdx][atomIdx])
            {
                if (!different_atomic_ff)
                {
                    cout << "difference for atom " << labels1[atomIdx] << "\n";
                    different_atomic_ff = true;
                }
                cout << "   " << hkl[hklIdx] << " " << ff1[hklIdx][atomIdx] << " " << ff2[hklIdx][atomIdx] << endl;
            }
    }

}

void calc_sf(
    const string& structure_file,
    const string& hkl_file)
{
    Crystal crystal;
    WallClockTimer timer;
    timer.start();
    structure_io::read_structure(structure_file, crystal);
    cout << "structure read in " << timer.stop() << " ms\n";


    vector<Vector3i> hkl;

    hkl_io::readHklIndices(hkl_file.c_str(), hkl);

    vector<complex<double> > sf;
    nlohmann::json json_data;
    ifstream jsonFileStream("aspher.json");

    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    timer.start();
    auto sf_calculator = SfCalculator::create(crystal, json_data);
    cout << "sf calculator created in " << timer.stop() << " ms\n";
    vector<bool> countAtom(crystal.atoms.size(), true);
    timer.start();
    cout << "calculating sf for " << crystal.atoms.size() << " atoms and " << hkl.size() << " reflections\n";
    sf_calculator->calculateStructureFactors(crystal.atoms, hkl, sf, countAtom);
    cout << "sf calculated in " << timer.stop() << " ms\n";
    ofstream out("cals_sf.txt");
    for (int i = 0; i < hkl.size(); i++)
        out << hkl[i] << " " << sf[i] << "\n";
    out.close();
    vector<TargetFunctionAtomicParamDerivatives> dTarget_dParam;
    vector<complex<double> > dTarget_dF;
    dTarget_dF.resize(hkl.size(),0.01);
    cout << "calculating sf and derivatives\n";
    timer.start();
    sf_calculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf, dTarget_dParam, dTarget_dF, countAtom);
    cout << "sf calculated in " << timer.stop() << " ms\n";
}

void read_cctb_sf(
    vector<Vector3i>& hkl,
    vector<complex<double> >& sf)
{
    ifstream in("cctbx_sf");
    string line;
    vector<string> words;
    while (getline(in, line))
    {
        string_utilities::split(line, words);
        if (words.size() == 5)
        {
            Vector3i h;
            complex<double> s;
            h[0] = stoi(words[0]);
            h[1] = stoi(words[1]);
            h[2] = stoi(words[2]);
            s.real(stod(words[3]));
            s.imag(stod(words[4]));
            hkl.push_back(h);
            sf.push_back(s);
        }
    }
    in.close();
}

void setAnomaluous(SfCalculator* sfCalculator, const Crystal crystal)
{
    vector<complex<double> > df{ {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0.002,0.002}, {0.004,0.003}, {0.008,0.006} };
    vector<complex<double> > anomalous;
    vector<int> atomic_numbers;
    crystal_structure_utilities::atomicNumbers(crystal, atomic_numbers);

    for(int z: atomic_numbers)
        anomalous.push_back(df[z]);
    sfCalculator->setAnomalous(anomalous);

}

void read_hkl_and_structure_factors(
    const string & inputFile,
    vector<Vector3i>& hkl,
    vector<complex<double> >& sf)
{
    cout << "input File " << inputFile << endl;
    if (inputFile == string("cctbx_sf"))
        read_cctb_sf(hkl, sf);
    else
    {
        filesystem::path p(inputFile);
        string extension = p.extension().string();
        cout << "extension " << extension << "\n";
        if(extension ==(".res") || extension == (".ins") || extension == (".cif"))
        {
            Crystal crystal;
            structure_io::read_structure(inputFile, crystal);

            nlohmann::json json_data;
            ifstream jsonFileStream("aspher.json");
            if (jsonFileStream.good())
                jsonFileStream >> json_data;
            jsonFileStream.close();
            auto sf_calculator = SfCalculator::create(crystal, json_data);
            //setAnomaluous(sf_calculator, crystal);
            vector<bool> countAtom(crystal.atoms.size(), true);
            sf_calculator->calculateStructureFactors(crystal.atoms, hkl, sf, countAtom);
        }
        else
        {
            SpaceGroup sg;
            UnitCell uc;
            vector<double> refln_F_meas, refln_F_sigma, refln_A_calc, refln_B_calc;
            cif_io::readFcf(inputFile, sg, uc, hkl, refln_F_meas, refln_F_sigma, refln_A_calc, refln_B_calc);
            for (int i = 0; i < hkl.size(); i++)
                sf.push_back(complex<double>(refln_A_calc[i], refln_B_calc[i]));
            
        }
    }

}

void adjust_sf(
    vector<complex<double> >& sf1,
    vector<complex<double> >& sf2,
    vector<Vector3i>& hkl1,
    vector<Vector3i>& hkl2)
{
    vector<int> idx1, idx2;

    vector<Vector3i> hkl;
    for (int i = 0; i < hkl1.size(); i++)
    {
        auto it = find(hkl2.begin(), hkl2.end(), hkl1[i]);
        if (it != hkl2.end())
        {
            idx1.push_back(i);
            idx2.push_back(distance(hkl2.begin(), it));
        }
    }

    vector<complex<double> > _sf01(idx1.size()), _sf02(idx2.size());
    
    for (int i = 0; i < idx1.size(); i++)
    {
        _sf01[i] = sf1[idx1[i]];
        _sf02[i] = sf2[idx2[i]];
        hkl.push_back(hkl1[idx1[i]]);
    }
    sf1 = _sf01;
    sf2 = _sf02;
    hkl1 = hkl; 
    hkl2 = hkl;
}


void compare_sf(
    const std::string &input_file_1,
    const std::string &input_file_2)
{
    vector<Vector3i> hkl, hkl0;

    vector<complex<double> > sf_model_1, sf_model_2;
    string extension1 = filesystem::path(input_file_1).extension().string();

    if (extension1 == ".tsc")
    {
        vector<string> labels;
        vector<vector<complex<double> > > atomic_form_factors;
        tsc_io::read_tsc(input_file_1, labels, hkl, atomic_form_factors);
        Crystal crystal;
        structure_io::read_structure(input_file_2, crystal);

        nlohmann::json json_data;
        json_data["model"] = "tsc";
        json_data["tsc file"] = input_file_1;
        auto sf_calculator = SfCalculator::create(crystal, json_data);
        setAnomaluous(sf_calculator, crystal);
        vector<bool> countAtom(crystal.atoms.size(), true);
        sf_calculator->calculateStructureFactors(crystal.atoms, hkl, sf_model_1, countAtom);
    }
    else
        read_hkl_and_structure_factors(input_file_1, hkl, sf_model_1);
    
    hkl0 = hkl;
    read_hkl_and_structure_factors(input_file_2, hkl, sf_model_2);
    if (!hkl0.empty())
        if (hkl0 != hkl)
        {
            //   on_error::throwException("hkl differ", __FILE__, __LINE__);
            cout << "hkl differ\n";
            adjust_sf(sf_model_1, sf_model_2, hkl0, hkl);
        }

    cout<< "agreementFactor  = " << agreementFactor(sf_model_1, sf_model_2) << endl;
    ofstream out("sf_diff.txt");
    vector<pair<double, int> > diff;
    for (int i = 0; i < hkl.size(); i++)
        diff.push_back(make_pair(abs(sf_model_1[i] - sf_model_2[i]), i));
    sort(diff.begin(), diff.end());
    for (int p = 0; p < hkl.size(); p++)
    {
        int idx = diff[p].second;
        out << hkl[idx] << " " << sf_model_1[idx] << " " << sf_model_2[idx] << "\n";
    }
    out.close();
}

void makeWholeHklSet(
    const vector<Vector3i>& hkl0,
    const SpaceGroup& spaceGroup,
    vector<Vector3i>& hkl)
{
    int nSymm = spaceGroup.nSymmetryOperationsInSubset();
    vector<Matrix3i> rotations(nSymm);

    for (int i = 0; i < nSymm; i++)
        spaceGroup.getSpaceGroupOperation(0, 0, i).getRotation(rotations[i]);
    for (int i = 0; i < nSymm; i++)
        rotations.push_back(-1 * rotations[i]);


    hkl.clear();
    set<Vector3i> uniqueHkl;
    vector<Vector3i> hkl1 = hkl0;

    for (auto const& rotation : rotations)
        for (auto const& h : hkl1)
            uniqueHkl.insert(h * rotation);

    hkl.assign(uniqueHkl.begin(), uniqueHkl.end());
}

void form_factors_calculation_time(
    const string& structureFile,
    const string& hklFile)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<complex<double> > formFactors;
    map<Vector3i, vector<complex<double> > > formFactorsMap;
    vector<Vector3i> hkl, hkl0;
    WallClockTimer timer;

    hkl_io::readHklIndices(hklFile, hkl0);
    makeWholeHklSet(hkl0, crystal.spaceGroup, hkl);
    nlohmann::json json_data;
    ifstream jsonFileStream("aspher.json");
    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    auto sf_calculator = SfCalculator::create(crystal, json_data);
    vector<bool> countAtom(crystal.atoms.size(), true);
    timer.start();
    for(int i=0;i<hkl.size();i++)
        sf_calculator->calculateFormFactors(hkl[i], formFactors, countAtom);
    cout << "form factors calculated in " << timer.stop() << " ms\n";

}


void sf_calculator_new_implementation(
    const string &structureFile,
    const string &tscOrHklFile)
{ 
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<vector<complex<double> > > formFactors;
    map<Vector3i, vector<complex<double> > > formFactorsMap;
    vector<Vector3i> hkl, hkl0;
    WallClockTimer timer;

    if (filesystem::path(tscOrHklFile).extension().string() == ".tsc")
    {
        vector<string> labels;
        tsc_io::read_tsc(tscOrHklFile, labels, hkl0, formFactors);
    }
    else
    {
        hkl_io::readHklIndices(tscOrHklFile,hkl0);
        makeWholeHklSet(hkl0, crystal.spaceGroup, hkl);
        nlohmann::json json_data;
        ifstream jsonFileStream("aspher.json");
        if (jsonFileStream.good())
            jsonFileStream >> json_data;
        jsonFileStream.close();
        auto sf_calculator = SfCalculator::create(crystal, json_data);
        vector<bool> countAtom(crystal.atoms.size(), true);
        timer.start();
        sf_calculator->calculateFormFactors(hkl, formFactors, countAtom);
        cout << "form factors calculated in " << timer.stop() << " ms\n";
    }

    for (int i = 0; i < hkl.size(); i++)
        formFactorsMap[hkl[i]] = formFactors[i];

    auto ffCalc = make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, formFactorsMap);
    AnyScattererStructureFactorCalculator calculator(crystal);
    AnyScattererStructureFactorCalculator2 calculator2(crystal);
    IamSfCalculator iamCalculator(crystal);
    auto _ffCalc = std::static_pointer_cast<AtomicFormFactorCalculationsManager>(ffCalc);
    calculator.setAtomicFormfactorManager(_ffCalc);
    calculator2.setAtomicFormfactorManager(_ffCalc);

    vector<complex<double> > sf1, sf2;
    
    
    timer.start();
    calculator.calculateStructureFactors(hkl0, sf1);
    cout<< "original implementation execution time (ms): " << timer.stop() << "\n";
    timer.start();
    //calculator2.calculateStructureFactors(hkl0, sf2);
    iamCalculator.calculateStructureFactors(hkl0, sf2);
    cout << "new implementation execution time (ms): " << timer.stop() << "\n";
    cout << "agreement factor = " << agreementFactor(sf1, sf2) << "\n";
    ofstream out("sf_diff.txt");
    for (int i = 0; i < hkl0.size(); i++)
        out << hkl0[i] << " " << sf1[i] << " " << sf2[i] << "\n";
    out.close();
     
    return;
 }

void iam_sf_calculator_new_implementation(
    const string& structureFile,
    const string& tscOrHklFile)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<vector<complex<double> > > formFactors;
    vector<Vector3i> hkl, hkl0;
    WallClockTimer timer;

    if (filesystem::path(tscOrHklFile).extension().string() == ".tsc")
    {
        vector<string> labels;
        tsc_io::read_tsc(tscOrHklFile, labels, hkl0, formFactors);
    }
    else
        hkl_io::readHklIndices(tscOrHklFile, hkl0);

    
    nlohmann::json json_data;
    json_data["model"] = "iam";
    auto current_iam_calculator = SfCalculator::create(crystal, json_data);
    timer.start();
    vector<bool> countAtom(crystal.atoms.size(), true);
    current_iam_calculator->calculateFormFactors(hkl, formFactors, countAtom);
    cout << "form factors calculated in " << timer.stop() << " ms\n";

    IamSfCalculator new_iam_calculator(crystal);

    vector<complex<double> > sf1, sf2;


    timer.start();
    current_iam_calculator->calculateStructureFactors(crystal.atoms, hkl0, sf1, countAtom);
    cout << "original implementation execution time (ms): " << timer.stop() << "\n";
    timer.start();
    
    new_iam_calculator.calculateStructureFactors(hkl0, sf2);
    cout << "new implementation execution time (ms): " << timer.stop() << "\n";
    cout << "agreement factor = " << agreementFactor(sf1, sf2) << "\n";
    //ofstream out("sf_diff.txt");
    //for (int i = 0; i < hkl0.size(); i++)
    //    out << hkl0[i] << " " << sf1[i] << " " << sf2[i] << "\n";
    //out.close();

    return;
}


void sort_hkl(const string& hklFile)
{
    vector<Vector3i> hkl;
    hkl_io::readHklIndices(hklFile.c_str(), hkl);
    vector<vector<Vector3i> > orderedHklLines;
    vector<vector<int> > mapToOriginalSetIndices;
    int direction = AnyScattererStructureFactorCalculator2::findPreferredHklOrderingDirection(hkl, orderedHklLines, mapToOriginalSetIndices);
    Vector3i shift(0, 0, 0);
    shift[direction] = 1;
    ofstream out("sorted_hkl.txt");
    for (auto& hklLine : orderedHklLines)
    {
        for (auto& h : hklLine)
            out << h << " ";
        out << "\n";
    }
    out.close();

    int nLineBreaks = 0;
    int breakLenghtTotal = 0;
    int longestBreak = 0;
    for (auto& hklLine : orderedHklLines)
        for (int i = 1; i < hklLine.size(); i++)
            if (hklLine[i] != shift + hklLine[i - 1])
            {
                nLineBreaks++;
                int breakLenght = hklLine[i][direction] - hklLine[i - 1][direction] + 1;
                if (breakLenght > longestBreak)
                    longestBreak = breakLenght;
                breakLenghtTotal += breakLenght;
            }
    cout << "n reflections = " << hkl.size() << "\n";
    cout << "n lines = " << orderedHklLines.size() << "\n";
    cout << "n line breaks = " << nLineBreaks << "\n";
    cout << "break length = " << breakLenghtTotal << "\n";
    cout << "longest break = " << longestBreak << "\n";
}

void testAdpCalc()
{
    Crystal crystal;
    crystal.unitCell.set(10, 10, 10, 90, 90, 90);
    crystal.atoms.resize(1);
    auto& atom = crystal.atoms[0];
    double u_iso = 0.01;
    atom.adp.push_back(u_iso);
    atom.type = "C";
    atom.coordinates = { 0.0, 0.0, 0.0 };

    AnyIamCalculator calc(crystal);
    vector<Vector3i> hkl = { {1,2,3}, {1,2,4}, {1,2,5}, {1,2,6}, {1,2,7}, {1,2,8} };
    vector<complex<double> > sf, sf_0_adps;
    calc.calculateStructureFactors(crystal.atoms, hkl, sf, vector<bool>(1, true));
    crystal.atoms[0].adp.clear();
    calc.calculateStructureFactors(crystal.atoms, hkl, sf_0_adps, vector<bool>(1, true));

    double t;
    ReciprocalLatticeUnitCell reciprocalUc(crystal.unitCell);
    Vector3d hklCart, l_cart;
    reciprocalUc.fractionalToCartesian(hkl[0], hklCart);
    double c_iter, c;
    double u_2_pi2 = u_iso * 2.0 * M_PI * M_PI;
    c = exp(-2 * u_2_pi2);
    reciprocalUc.fractionalToCartesian(hkl[0], hklCart);
    reciprocalUc.fractionalToCartesian(Vector3d(0,0,1), l_cart);
    double beta = l_cart * l_cart;
    c = exp(-2 * beta * u_2_pi2);
    double t0 = exp(-hklCart * hklCart * u_2_pi2);
    reciprocalUc.fractionalToCartesian(hkl[1], hklCart);
    double t1 = exp(-hklCart * hklCart * u_2_pi2);
    c_iter = t1 / t0;

    for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
    {
        reciprocalUc.fractionalToCartesian(hkl[hklIdx], hklCart);
        double h2 = hklCart * hklCart;
        //t = exp(-h2 * u_2_pi2);
        if (hklIdx == 0)
            t = t0;
        if(hklIdx == 1)
            t = t1;
        if (hklIdx > 1)
        {
            c_iter *= c;
            t *= c_iter;
        }
        cout << setw(20) << hkl[hklIdx] << " " << setw(30) << sf[hklIdx] << "  " << sf_0_adps[hklIdx] << " " << sf[hklIdx].real() / sf_0_adps[hklIdx].real() << " " << t << endl;
    }
    //##########################
    //   U anisotropic
    //##########################
    cout << "test anisotropic adp\n";

    crystal.atoms[0].adp = { 0.02200, 0.01630, 0.01257, -0.00795, 0.00475, -0.00413 };
    calc.calculateStructureFactors(crystal.atoms, hkl, sf, vector<bool>(1, true));
    vector<double> u_cart(6);
    StructuralParametersConverter converter;
    converter.convertADP(crystal.atoms[0].adp, u_cart, structural_parameters_convention::AdpConvention::U_cif, structural_parameters_convention::AdpConvention::U_cart);
    Matrix3d U(u_cart[0], u_cart[3], u_cart[4],
               u_cart[3], u_cart[1], u_cart[5],
               u_cart[4], u_cart[5], u_cart[2]);
    Matrix3d m = U* (2.0 * M_PI * M_PI);
    reciprocalUc.fractionalToCartesian(hkl[0], hklCart);
    c = exp(-2.0 * l_cart * m * l_cart);
    t0 = exp(-hklCart * m * hklCart);
    reciprocalUc.fractionalToCartesian(hkl[1], hklCart);
    t1 = exp(-hklCart * m * hklCart);
    c_iter = t1 / t0;

    for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
    {
        reciprocalUc.fractionalToCartesian(hkl[hklIdx], hklCart);
        double h2 = hklCart * hklCart;
        //t = exp(-h2 * u_2_pi2);
        if (hklIdx == 0)
            t = t0;
        if (hklIdx == 1)
            t = t1;
        if (hklIdx > 1)
        {
            c_iter *= c;
            t *= c_iter;
        }
        cout << setw(20) << hkl[hklIdx] << " " << setw(30) << sf[hklIdx] << "  " << sf_0_adps[hklIdx] << " " << sf[hklIdx].real() / sf_0_adps[hklIdx].real() << " " << t << endl;
    }

}

void check_line(int argc, char* argv[])
{
    vector<int> line;
    for(int i=1;i<argc;i++)
        line.push_back(stoi(argv[i]));
    cout << "line points_before_break points_after_break\n";
    for (int idx_in_line = 0; idx_in_line < line.size(); idx_in_line++)
    {
        int points_after_break = 0;
        int points_before_break = 0;
        for (int i = idx_in_line + 1; i < line.size() ; i++)
        {
            if (line[i] == line[i - 1] + 1)
                points_before_break++;
            else
                break;
        }
        for (int i = idx_in_line - 1; i >= 0; i--)
        {
            if (line[i] == line[i + 1] - 1)
                points_after_break++;
            else
                break;
        }
        cout << line[idx_in_line] << " " << points_before_break << " " << points_after_break << "\n";
    }

}

void calc_derivatives_numerically(
    const Crystal &crystal,
    SfCalculator* sfCalculator,
    const Vector3i& hkl,
    SfDerivativesAtHkl& derivatives,
    const vector<bool>& countAtom)
{
    double step=0.0001;
    double stepAdp = 0.001;
    complex<double> sf0, sf_plus, sf_minus;
    SfDerivativesAtHkl derivatives0, derivatives_plus, derivatives_minus;
    sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf0, derivatives0, countAtom);

    int nAtoms = crystal.atoms.size();
    derivatives.atomicPostionDerivatives.resize(nAtoms,Vector3<complex<double> >(0.0,0.0,0.0));
    derivatives.adpDerivatives.resize(nAtoms);
    derivatives.occupancyDerivatives.resize(nAtoms, 0.0);
    for (int coordIdx=0; coordIdx<3; coordIdx++)
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            auto atoms = crystal.atoms;
            atoms[atomIdx].coordinates[coordIdx] += step;
            sfCalculator->update(atoms);
            sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_plus, derivatives_plus, countAtom);
            atoms[atomIdx].coordinates[coordIdx] -= 2*step;
            sfCalculator->update(atoms);
            sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_minus, derivatives_minus, countAtom);
            derivatives.atomicPostionDerivatives[atomIdx][coordIdx] = (sf_plus - sf_minus) / (2 * step);
            atoms[atomIdx].coordinates[coordIdx] -= step;
        }

    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        auto atoms = crystal.atoms;
        atoms[atomIdx].occupancy += step;
        sfCalculator->update(atoms);
        sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_plus, derivatives_plus, countAtom);
        atoms[atomIdx].occupancy -= 2 * step;
        sfCalculator->update(atoms);
        sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_minus, derivatives_minus, countAtom);
        derivatives.occupancyDerivatives[atomIdx] = (sf_plus - sf_minus) / (2 * step);
        atoms[atomIdx].occupancy += step;
        derivatives.adpDerivatives[atomIdx].resize(atoms[atomIdx].adp.size());

        for (int i = 0; i < atoms[atomIdx].adp.size(); i++)
        {
            atoms[atomIdx].adp[i] += stepAdp;
            sfCalculator->update(atoms);
            sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_plus, derivatives_plus, countAtom);
            atoms[atomIdx].adp[i] -= 2 * stepAdp;
            sfCalculator->update(atoms);
            sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_minus, derivatives_minus, countAtom);
            derivatives.adpDerivatives[atomIdx][i] = (sf_plus - sf_minus) / (2 * stepAdp);
            atoms[atomIdx].adp[i] -= stepAdp;
        }
    }

}

void calc_derivatives_numerically_shift_all(
    const Crystal& crystal,
    SfCalculator* sfCalculator,
    const Vector3i& hkl,
    SfDerivativesAtHkl& derivatives,
    const vector<bool>& countAtom)
{
    derivatives.atomicPostionDerivatives.resize(crystal.atoms.size());
    double step = 0.0001;
    complex<double> sf0, sf_plus, sf_minus;
    SfDerivativesAtHkl derivatives0, derivatives_plus, derivatives_minus;
    sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf0, derivatives0, countAtom);



    int nAtoms = crystal.atoms.size();
    vector<double> occupancies;
    for (auto& atom : crystal.atoms)
        occupancies.push_back(atom.occupancy);

    auto atoms = crystal.atoms;
    
    for (int coordinateIdx = 0; coordinateIdx < 3; coordinateIdx++)
    {
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            for (auto& atom : atoms)
            {
                atom.occupancy = 0.0;
                atom.coordinates[coordinateIdx] += step;
            }

            atoms[atomIdx].occupancy = occupancies[atomIdx];
            sfCalculator->update(atoms);
            sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_plus, derivatives_plus, countAtom);

            for (auto& atom : atoms)
                atom.coordinates[coordinateIdx] -= 2 * step;

            sfCalculator->update(atoms);
            sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_minus, derivatives_minus, countAtom);
            
            derivatives.atomicPostionDerivatives[atomIdx][coordinateIdx] = (sf_plus - sf_minus) / (2 * step);
            for (auto& atom : atoms)
                atom.coordinates[coordinateIdx] += step;

        }

        sfCalculator->update(crystal.atoms);
    }

}


void derivatives(
    const string& structureFile,
    const string& hklFile)
{
    
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<vector<complex<double> > > formFactors;
    vector<Vector3i> hkl;
    
    hkl_io::readHklIndices(hklFile, hkl);

    int nHkl = hkl.size();
    int nAtoms = crystal.atoms.size();

    nlohmann::json json_data;
    ifstream jsonFileStream("aspher.json");
    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    auto sfCalculator = SfCalculator::create(crystal, json_data);

    vector<bool> countAtom(crystal.atoms.size(), true);
    vector<complex<double> > sf, sf_plus, sf_minus;
    vector<TargetFunctionAtomicParamDerivatives> dT_dp;
    vector<complex<double> > dt_df(nHkl, 0.0);
    // x y z U occ
    SfDerivativesAtHkl derivatives_ant, derivatives_num, derivatives_num_shift_all;
    //sfCalculator->calculateStructureFactors(crystal.atoms, hkl, sf, countAtom);
    dt_df[0] = 1.0;
    sfCalculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf, dT_dp, dt_df, countAtom);
    for (int i = 0; i < dT_dp.size(); i++)
        cout << i << " " << dT_dp[i].atomic_position_derivatives << "\n";
    
    ofstream out("derivatives");

    for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
    {
        out << "\n HKL = " << hkl[hklIdx] << "\n\n";
        calc_derivatives_numerically(crystal, sfCalculator, hkl[hklIdx], derivatives_num, countAtom);
        calc_derivatives_numerically_shift_all(crystal, sfCalculator, hkl[hklIdx], derivatives_num_shift_all, countAtom);
        sfCalculator->calculateStructureFactorsAndDerivatives(hkl[hklIdx], sf[hklIdx], derivatives_ant, countAtom);
        out << "dF/dxyz\n";
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            out << atomIdx << "\n";
            for (int i = 0; i < 3; i++)
                out << "    " << i << " " << fixed << setprecision(6) << derivatives_ant.atomicPostionDerivatives[atomIdx][i] << "  "
                << fixed << setprecision(6) << derivatives_num.atomicPostionDerivatives[atomIdx][i] << "  "
                << fixed << setprecision(6) << derivatives_num_shift_all.atomicPostionDerivatives[atomIdx][i] << "\n";
        }
        out << "dF/d_occ\n";
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            out << atomIdx << " " << fixed << setprecision(6) << derivatives_ant.occupancyDerivatives[atomIdx] << "  "
            << fixed << setprecision(6) << derivatives_num.occupancyDerivatives[atomIdx] << "\n";
        out << "dF/d_adp\n";
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            out << atomIdx;
            for (int i = 0; i < derivatives_ant.adpDerivatives[atomIdx].size(); i++)
                out << "   " << i << " " << fixed << setprecision(6) << derivatives_ant.adpDerivatives[atomIdx][i] << "  "
                << fixed << setprecision(6) << derivatives_num.adpDerivatives[atomIdx][i] << "\n";
        }
    }

    for (int hklIdx = 0; hklIdx < hkl.size(); hklIdx++)
    {
        out << "\n HKL = " << hkl[hklIdx] << "\n\n";
        calc_derivatives_numerically(crystal, sfCalculator, hkl[hklIdx], derivatives_num, countAtom);
        calc_derivatives_numerically_shift_all(crystal, sfCalculator, hkl[hklIdx], derivatives_num_shift_all, countAtom);
        sfCalculator->calculateStructureFactorsAndDerivatives(hkl[hklIdx], sf[hklIdx], derivatives_ant, countAtom);
        out << "dF/dxyz\n";
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            for (int i = 0; i < 3; i++)
            {
                complex<double> v1 = derivatives_ant.atomicPostionDerivatives[atomIdx][i];
                complex<double> v2 = derivatives_ant.atomicPostionDerivatives[atomIdx][i];
                if (fabs(v1.real()) > 1e-10 || fabs(v2.real()) > 1e-10)
                    if (fabs(v1.real() - v2.real())/ (fabs(v1.real()) + fabs(v2.real())) > 1e-4)
                        out<< atomIdx << " " << i << " real\n";
                if (fabs(v1.imag()) > 1e-10 || fabs(v2.imag()) > 1e-10)
                    if (fabs(v1.imag() - v2.imag()) / (fabs(v1.imag()) + fabs(v2.imag())) > 1e-4)
                        out << atomIdx << " " << i << " imag\n";

            }
        }
        out << "dF/d_occ\n";
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            complex<double> v1 = derivatives_ant.occupancyDerivatives[atomIdx];
            complex<double> v2 = derivatives_ant.occupancyDerivatives[atomIdx];
            if (fabs(v1.real()) > 1e-10 || fabs(v2.real()) > 1e-10)
                if (fabs(v1.real() - v2.real()) / (fabs(v1.real()) + fabs(v2.real())) > 1e-4)
                    out << atomIdx << " " << " real\n";
            if (fabs(v1.imag()) > 1e-10 || fabs(v2.imag()) > 1e-10)
                if (fabs(v1.imag() - v2.imag()) / (fabs(v1.imag()) + fabs(v2.imag())) > 1e-4)
                    out << atomIdx << " " << " imag\n";
        }
        out << "dF/d_adp\n";
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            out << atomIdx;
            for (int i = 0; i < derivatives_ant.adpDerivatives[atomIdx].size(); i++)
            {
                complex<double> v1 = derivatives_ant.adpDerivatives[atomIdx][i];
                complex<double> v2 = derivatives_ant.adpDerivatives[atomIdx][i];
                if (fabs(v1.real()) > 1e-10 || fabs(v2.real()) > 1e-10)
                    if (fabs(v1.real() - v2.real()) / (fabs(v1.real()) + fabs(v2.real())) > 1e-4)
                        out << atomIdx << " " << i << " real\n";
                if (fabs(v1.imag()) > 1e-10 || fabs(v2.imag()) > 1e-10)
                    if (fabs(v1.imag() - v2.imag()) / (fabs(v1.imag()) + fabs(v2.imag())) > 1e-4)
                        out << atomIdx << " " << i << " imag\n";

            }
        }
    }


    out.close();
}

void printFractional(
    const string &structureFile, 
    const string& xyzFile)
{
    Crystal crystal;
    MoleculeData molData;
    structure_io::read_structure(structureFile, crystal);
    xyz_io::readXyz(xyzFile, molData);
    for (int i = 0; i < molData.atomPositions.size(); i++)
    {
        Vector3d r_frac;
        crystal.unitCell.cartesianToFractional(molData.atomPositions[i], r_frac);
        cout << left << setw(8) << periodic_table::symbol(molData.atomicNumbers[i]) << r_frac << endl;
    }

}

int main(int argc, char* argv[])
{
    try {

        printFractional(argv[1], argv[2]);
        return 0;

        iam_sf_calculator_new_implementation(argv[1], argv[2]);
        return 0;


        sf_calculator_new_implementation(argv[1], argv[2]);
        return 0;
        
        form_factors_calculation_time(argv[1], argv[2]);
        return 0;

        derivatives(argv[1], argv[2]);
        return 0;


        check_line(argc, argv);
        return 0;

        testAdpCalc();
        return 0;

        sort_hkl(argv[1]);
        return 0;

        calc_sf(argv[1], argv[2]);
        return 0;

        if (argc != 3)
        {
            cout << "expected 3 arguments:\n"
                " (1) structure file or fcf or cctbx_sf\n"
                " (2) fcf or cctbx_sf file\n"
                " or\n"
                " (1) tsc file\n"
                " (2) structure file\n";
            exit(0);
        }
        compare_sf(argv[1], argv[2]);
        return 0;




        tsc_differences(argv[1], argv[2]);
        return 0;
        
        string structFile = "";
        if (argc == 2)
            structFile = argv[1];
        test_connectivity(structFile);
        return 0;

        test_int_floor();
        return 0;

        check_args(argc, argv, 2, { "structure file", "bank file path" });
        testTypeAssignemntLog(argv[1], argv[2]);
        return 0;

        check_args(argc, argv, 2, { "structure file", "bank file path"});
        Crystal crystal;
        structure_io::read_structure(argv[1], crystal);
        TAAM_types(crystal, argv[2]);
        return 0;

        check_args(argc, argv, 1, { "structure file" });
        nGaussian(argv[1]);
        return 0;
        compare_sf_taam_qm(argc, argv);
        return 0;
        test_taam_disorder(argc, argv);

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

