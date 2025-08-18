#define _CRTDBG_MAP_ALLOC

#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/crystallographic_point_group_tables.h"
#include "discamb/IO/cif_io.h"
#include "discamb/IO/discamb_io.h"
#include "discamb/IO/NativeIAM_Reader.h"
#include "discamb/IO/fragmentation_io.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/IO/xyz_io.h"
#include "discamb/MathUtilities/SphConverter.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/Scattering/agreement_factors.h"
#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/gar_utilities.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator.h"
#include "discamb/Scattering/HcAtomBankStructureFactorCalculator2.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/Scattering/HirshfeldAtomModelSettings.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"
#include "discamb/Scattering/IamSfCalculator.h"
#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/Scattering/statistics.h"
#include "discamb/Scattering/taam_utilities.h"
#include "discamb/Scattering/TaamSfCalculatorMultiOrderedImpl.h"
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
    //cout << "scale = " << scale << "\n";
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
        bankSettings.wfn_databank,
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

    for (int z : atomic_numbers)
        if (z < 9)
            anomalous.push_back(df[z]);
        else
            anomalous.push_back(0.0);
    sfCalculator->setAnomalous(anomalous);

}

void setAnomaluous(shared_ptr<SfCalculator> &sfCalculator, const Crystal crystal)
{
    vector<complex<double> > df{ {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0.002,0.002}, {0.004,0.003}, {0.008,0.006} };
    vector<complex<double> > anomalous;
    vector<int> atomic_numbers;
    crystal_structure_utilities::atomicNumbers(crystal, atomic_numbers);

    for (int z : atomic_numbers)
        if (z < 9)
            anomalous.push_back(df[z]);
        else
            anomalous.push_back(0.0);
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
    int direction = scattering_utilities::findPreferredHklOrderingDirection(hkl, orderedHklLines, mapToOriginalSetIndices);
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

void copy_compare()
{
    auto starting_point = filesystem::current_path();

    for (auto entry : filesystem::directory_iterator(starting_point))
        if(filesystem::is_directory(entry.path()))
        {
            vector<string> words;
            string_utilities::split(entry.path().string(), words, '\\');
            string dirName = words.back();
            if (dirName == "compare")
                continue;

            filesystem::current_path(entry.path());
            string res = file_system_utilities::find_newest_file("res");
            string cif = file_system_utilities::find_newest_file("cif");
            string npy;
            vector<string> npyFiles;
            file_system_utilities::find_files("npy", npyFiles);
            for (auto npyFile : npyFiles)
                if (npyFile.find("discamb") == string::npos)
                    npy = npyFile;
            filesystem::copy(res, starting_point / "compare" / (dirName + ".res"));
            filesystem::copy(cif, starting_point / "compare" / (dirName + ".cif"));
            filesystem::copy(npy, starting_point / "compare" / (dirName + ".npy"));
            filesystem::current_path(starting_point);
        }
}

//double relative_l1_100_agreement_factor(
//    const vector<double> &v1,
//    const vector<double> &v2)
//{
//    double numerator, denominator;
//    int i, n = v1.size();
//    if (n != v2.size())
//        on_error::throwException("vectors size do not match when trying to calculate L1 agreement factor", __FILE__, __LINE__);
//    numerator = 0.0;
//    denominator = 0.0;
//
//    for (i = 0; i < n; i++)
//    {
//        numerator += abs(v1[i] - v2[i]);
//        denominator += abs(v1[i]) + abs(v2[i]);
//    }
//    if(denominator == 0)
//        on_error::throwException("denominator 0 when trying to calculate L1 agreement factor", __FILE__, __LINE__);
//
//    return numerator / denominator * 200;
//}
//
//void compare_derivatives(
//    const vector<TargetFunctionAtomicParamDerivatives>& dT_dp1,
//    const vector<TargetFunctionAtomicParamDerivatives>& dT_dp2,
//    bool &size_match,
//    double &d_xyz_agreement_factor,
//    double& d_adp_agreement_factor,
//    double& d_occ_agreement_factor,
//    double& d_fpfdp_agreement_factor)
//{
//    /*
//    Vector3<REAL> atomic_position_derivatives;
//    std::vector<REAL> adp_derivatives;
//    REAL occupancy_derivatives = 0;
//    */
//    int atomIdx, nAtoms = dT_dp1.size();
//    size_match = (nAtoms == dT_dp2.size());
//    vector<double> xyz1, xyz2, adps1, adps2, occ1, occ2;
//    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
//    {
//        if (dT_dp1[atomIdx].adp_derivatives.size() != dT_dp2[atomIdx].adp_derivatives.size())
//        {
//            size_match = false;
//            return;
//        }
//        for (int i = 0; i < 3; i++)
//        {
//            xyz1.push_back(dT_dp1[atomIdx].atomic_position_derivatives[i]);
//            xyz2.push_back(dT_dp2[atomIdx].atomic_position_derivatives[i]);
//        }
//        for (int i = 0; i < dT_dp1[atomIdx].adp_derivatives.size(); i++)
//        {
//            adps1.push_back(dT_dp1[atomIdx].adp_derivatives[i]);
//            adps2.push_back(dT_dp2[atomIdx].adp_derivatives[i]);
//        }
//        occ1.push_back(dT_dp1[atomIdx].occupancy_derivatives);
//        occ2.push_back(dT_dp2[atomIdx].occupancy_derivatives);
//    }
//    
//    d_xyz_agreement_factor=relative_l1_100_agreement_factor(xyz1, xyz2); 
//    d_adp_agreement_factor=relative_l1_100_agreement_factor(adps1, adps2);
//    d_occ_agreement_factor=relative_l1_100_agreement_factor(occ1, occ2);
//
//    //for(int )
//}

void taam_parallel(
    const string& structureFile,
    const string& hklFile,
    bool calcOnlyNew,
    bool sf_and_sfAndDf,
    bool anomalous)
{
    
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<vector<complex<double> > > formFactors;
    vector<Vector3i> hkl_all, hkl;
    hkl_io::readHklIndices(hklFile, hkl);
    //hkl_io::readHklIndices(hklFile, hkl_all);
    //set<Vector3i> uniqueHkl(hkl_all.begin(), hkl_all.end());
    //hkl.insert(hkl.end(),uniqueHkl.begin(), uniqueHkl.end());

    int nHkl = hkl.size();
    int nAtoms = crystal.atoms.size();

    nlohmann::json json_data;
    ifstream jsonFileStream("aspher.json");
    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    //auto sfCalculator = SfCalculator::create(crystal, json_data);
    shared_ptr<SfCalculator> sfCalculator = shared_ptr<SfCalculator>(SfCalculator::create(crystal, json_data));
    json_data["algorithm"] = "macromol";
    shared_ptr<SfCalculator> sfCalculator2 = shared_ptr<SfCalculator>(SfCalculator::create(crystal, json_data));

    if (anomalous)
    {
        setAnomaluous(sfCalculator, crystal);
        setAnomaluous(sfCalculator2, crystal);
    }
    

    vector<bool> countAtom(crystal.atoms.size(), true);
    vector<complex<double> > sf_new, sf_new2, sf_standard, sf_standard2, sf_plus, sf_minus;
    vector<TargetFunctionAtomicParamDerivatives> dT_dp, dT_dp_new;
    vector<complex<double> > dt_df(nHkl, complex<double>(1.0,0.3));
    // x y z U occ
    SfDerivativesAtHkl derivatives_ant, derivatives_num, derivatives_num_shift_all;
    WallClockTimer timer;

    // structure factors call 'calculateStructureFactors'

    if (sf_and_sfAndDf)
    {
        if (!calcOnlyNew)
        {
            timer.start();
            sfCalculator->calculateStructureFactors(crystal.atoms, hkl, sf_standard, countAtom);
            cout << "time = " << timer.stop() << " ms\n";
        }
        timer.start();
        sfCalculator2->calculateStructureFactors(crystal.atoms, hkl, sf_new, countAtom);
        cout << "time = " << timer.stop() << " ms\n";
        cout << "agreement factor " << agreementFactor(sf_new, sf_standard) << "\n";
    }

    // structure factors and derivatives:  'calculateStructureFactorsAndDerivatives'
    if (!calcOnlyNew)
    {
        //auto sfCalculator = SfCalculator::create(crystal, json_data);
        timer.start();
        sfCalculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_standard2, dT_dp, dt_df, countAtom);
        cout << "time = " << timer.stop() << " ms\n";
    }

    DerivativesSelector derivativesSelector;
    //derivativesSelector.d_adp = false;
    //derivativesSelector.d_xyz = false;
    timer.start();

    sfCalculator2->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_new2, dT_dp_new, dt_df, countAtom, derivativesSelector);
    cout << "time = " << timer.stop() << " ms\n";
    if(!calcOnlyNew)
        cout << "agreement factor " << agreement_factors::value(sf_new2, sf_standard2) << "\n";

    if (sf_and_sfAndDf)
    {
        cout << "agreement factor standar sf vs sf_dsf " << agreementFactor(sf_standard, sf_standard2) << "\n";
        cout << "agreement factor new sf vs sf_dsf " << agreementFactor(sf_new, sf_new2) << "\n";
    }


    //return;

    // check calculateStructureFactorsAndDerivatives
    //cout<< "calculate structure factors ansd derivatives\n";
    //timer.start();
    //sfCalculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_standard, dT_dp, dt_df, countAtom);
    //cout << "time = " << timer.stop() << " ms\n";
    //timer.start();
    //sfCalculator2->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_new, dT_dp_new, dt_df, countAtom);
    //cout << "time = " << timer.stop() << " ms\n";

    bool size_match;
    double d_xyz_agreement_factor, d_adp_agreement_factor, d_occ_agreement_factor, d_fpfdp_agreement_factor;
    agreement_factors::for_derivatives(
        dT_dp,
        dT_dp_new,
        size_match,
        d_xyz_agreement_factor, d_adp_agreement_factor, d_occ_agreement_factor, d_fpfdp_agreement_factor);

    cout << "xyz : " << d_xyz_agreement_factor << "\n"
         << "adp : " << d_adp_agreement_factor << "\n"
         << "occ : " << d_occ_agreement_factor << "\n";


    //cout << "size match " << boolalpha << size_match << "\n"
    //    << "agreement factors\n"
    //    << "xyz : " << d_xyz_agreement_factor << "\n"
    //    << "adp : " << d_adp_agreement_factor << "\n"
    //    << "occ : " << d_occ_agreement_factor << "\n";
    

}


void convert_sph_harmonics_test()
{
    double v = 1.0 / sqrt(2.0);
    vector<vector<double> > newCoordinates{ {v,v,0.0}, {-v,v,0.0},{0.0,0.0,1.0} };
    vector<vector<double> > oldCoordinates{ {1.0,0.0,0.0}, {0.0,1.0,0.0},{0.0,0.0,1.0} };
    SphConverter sphConverter;
    vector<vector<vector<double> > > conversionMatrices;
    sphConverter.setMaxL(4);
    sphConverter.convert(oldCoordinates, newCoordinates, conversionMatrices);
    cout << conversionMatrices.size() << endl;
    vector<double> sph_l1_conv(3), sph_l1{0,0,1};

    for (int i = 0; i < 3; i++)
    {
        sph_l1_conv[i] = 0;
        for (int j = 0; j < 3; j++)
            sph_l1_conv[i] += conversionMatrices[1][i][j] * sph_l1[j];
        cout << sph_l1_conv[i] << " ";
    }

}

void sort_split_hkl(const string& hklFile, int nSets)
{
    vector<Vector3i> hkl;
    hkl_io::readHklIndices(hklFile.c_str(), hkl);
    vector<vector<Vector3i> > orderedHklLines;
    vector<vector<int> > mapToOriginalSetIndices;
    int direction = scattering_utilities::findPreferredHklOrderingDirection(hkl, orderedHklLines, mapToOriginalSetIndices);
    vector<vector<vector<Vector3i> > > splitOrderedLines;
    vector < vector < vector <pair<int, int> > > > splitOrderedIndices;
    vector<vector<vector<int> > > subsetDataHklIdxInOryginalSet;
    
    scattering_utilities::splitHklLines(nSets, orderedHklLines, mapToOriginalSetIndices, splitOrderedLines, splitOrderedIndices, subsetDataHklIdxInOryginalSet);
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


    ofstream out2("sorted_split_hkl.txt");

    for (int setIdx = 0; setIdx < nSets; setIdx++)
    {
        out2 << "\n SET " << setIdx + 1 << "\n\n";

        for (int i=0;i< splitOrderedLines[setIdx].size(); i++)
        {
            for (int j=0;j< splitOrderedLines[setIdx][i].size(); j++)
                out2 << splitOrderedLines[setIdx][i][j] << "  " << splitOrderedIndices[setIdx][i][j].first << " " << splitOrderedIndices[setIdx][i][j].second 
                     << " " << subsetDataHklIdxInOryginalSet[setIdx][i][j] << " , ";
            out2 << "\n";
        }
    }
    out2.close();

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

void compare_olex_taam(
    const string &fcf_file,
    const string &dis_file)
{
    vector<cif_io::DataSet> data_sets;
    cif_io::readCif(fcf_file, data_sets);
    cif_io::DataSet& fcf = data_sets[0];
    vector<Vector3i> hkl;
    vector<complex<double> > sf_fcf, sf_discamb;
    vector<int> tagsIndices;
    int loopIdx;
    bool hassLoop = cif_io::findLoopIndex(fcf, { "_refln_index_h" }, loopIdx, tagsIndices);
    vector<string>& h_str = fcf.loops[loopIdx].values[0];
    vector<string>& k_str = fcf.loops[loopIdx].values[1];
    vector<string>& l_str = fcf.loops[loopIdx].values[2];
    vector<string>& a_str = fcf.loops[loopIdx].values[5];
    vector<string>& b_str = fcf.loops[loopIdx].values[6];

    for (int i = 0; i < h_str.size(); i++)
    {
        hkl.push_back({ stoi(h_str[i]), stoi(k_str[i]) ,stoi(l_str[i]) });
        sf_fcf.push_back({ stod(a_str[i]), stod(b_str[i]) });
    }

    ////

    Crystal crystal;

    structure_io::read_structure(dis_file, crystal);
    //vector<vector<complex<double> > > formFactors;
    //vector<Vector3i> hkl;

    //hkl_io::readHklIndices(hklFile, hkl);

    int nHkl = hkl.size();
    int nAtoms = crystal.atoms.size();

    nlohmann::json json_data;
    ifstream jsonFileStream("aspher.json");
    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    auto sfCalculator = SfCalculator::create(crystal, json_data);

    vector<bool> countAtom(crystal.atoms.size(), true);
    sfCalculator->calculateStructureFactors(crystal.atoms, hkl, sf_discamb, countAtom);
    cout << "agreement factor " << agreementFactor(sf_discamb, sf_fcf);

}
struct Choice {
    bool dx = false;
    bool docc = false;
    bool dU = false;
    bool dfpfdp = false;
    void set(int option)
    {

    }
};

void convert(int& optionInt, Choice& optionStruc)
{

}

void options()
{
    // dx, docc, dU, dfpfdp
    
    int options;
    enum class settings  {dx = 0, docc=2, dU=4, dfpfdp=8};

}

void make_hkl(
    const string& structureFile,
    const string& indicesFile,
    const std::string &_jsonFile)
{
    string jsonFile = _jsonFile;
    if (jsonFile.empty())
        jsonFile = "aspher.json";
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<Vector3i> hkl;
    hkl_io::readHklIndices(indicesFile, hkl);
    vector<complex<double> > sf;
    vector<bool> countAtom(crystal.atoms.size(), true);
    auto sfCalculator = SfCalculator::create(crystal, jsonFile);
    sfCalculator->calculateStructureFactors(crystal.atoms, hkl, sf, countAtom);
    
    
    /*
    void writeShelxHkl(
        const std::string & fileName,
        const std::vector<Vector3i> &hkl,
        const std::vector<double> &intensities,
        const std::vector<double> &sigma,
        const std::vector<int> &batchNumber,
        bool freeFormat);
    */

    vector<double> intensities;
    vector<double> sigma;
    for (int i = 0; i < sf.size(); i++)
    {
        intensities.push_back(abs(sf[i]) * abs(sf[i]));
        sigma.push_back(0.02 * intensities.back() + 0.01);
    }
    vector<int> batchNumber(sf.size(), 1);


    hkl_io::writeShelxHkl("out.hkl", hkl, intensities, sigma, batchNumber, false);
}

void readFragments(
    const string& fragmentsFile,
    const Crystal& crystal,
    vector < vector <pair<int, std::string> > >& fragmentAtoms,
    vector<vector<int> > &atomsToAssign)
{
    fragmentAtoms.clear();
    atomsToAssign.clear();

    ifstream in(fragmentsFile);
    if (!in.good())
        on_error::throwException("cannot open fragments file " + fragmentsFile, __FILE__, __LINE__);

    string line;

    while (getline(in, line))
    {
        if (line.empty())
            continue;
        if (line[0] == '#')
            continue; // skip comments

        vector<int> fragmentAtomsToAssign;
        vector < pair<int, string> > atomListIntStr;
        vector < pair<string, string> > atomListStrStr;
        int nAtoms = stoi(line);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            getline(in, line);
            vector<string> words;
            string_utilities::split(line, words, ' ');
            if (words.size() < 2)
                on_error::throwException("incorrect line in fragments file " + fragmentsFile + ": " + line, __FILE__, __LINE__);
            double weight = stod(words[0]);
            if(weight > 0.0)
                fragmentAtomsToAssign.push_back(atomIdx);
            string label = words[1];
            string symmOp = "x,y,z"; // default
            if (words.size() == 3)
                symmOp = words[2];
            atomListStrStr.push_back({ label, symmOp });
        }
        crystal_structure_utilities::convertAtomList(crystal, atomListStrStr, atomListIntStr);
        fragmentAtoms.push_back(atomListIntStr);
        atomsToAssign.push_back(fragmentAtomsToAssign);
    }
}

void type_assign_disorder(
    const string &structureFile,
    const string &fragmentsFile,
    const string& bankFile)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector < vector <pair<int, std::string> > > fragmentAtoms;
    vector<vector<int> > atomsToAssign;
    readFragments(fragmentsFile, crystal, fragmentAtoms, atomsToAssign);
    CrystalAtomTypeAssigner assigner;
    vector< vector<int> > typeID;
    vector< vector<LocalCoordinateSystem<AtomInCrystalID> > > lcs;

    //--------
    MATTS_BankReader bankReader;
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> typeParameters;
    BankSettings bankSettings;

    bankReader.read(bankFile, atomTypes, typeParameters, bankSettings);

    assigner.setAtomTypes(atomTypes);
 
    assigner.assign(crystal, fragmentAtoms, atomsToAssign, typeID, lcs);
    //--------
    ofstream out("assigned_types.txt");
    for (int i = 0; i < fragmentAtoms.size(); i++)
    {
        out << "Fragment " << i << "\n";
        for (int j = 0; j < atomsToAssign[i].size(); j++)
        {
            int atomIdx = fragmentAtoms[i][atomsToAssign[i][j]].first;
            out << "    " << crystal.atoms[atomIdx].label << " ";
            if(typeID[i][j]>=0)
                out << atomTypes[typeID[i][j]].id << " ";
            else
                out << "no type assigned ";
            out << "\n";
        }
    }
    out.close();

}

void type_assign_order(
    const string& structureFile,
    const string& bankFile)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    
    CrystalAtomTypeAssigner assigner;
    vector<int> typeID;
    vector<LocalCoordinateSystem<AtomInCrystalID> > lcs;

    //--------
    MATTS_BankReader bankReader;
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> typeParameters;
    BankSettings bankSettings;

    bankReader.read(bankFile, atomTypes, typeParameters, bankSettings);

    assigner.setAtomTypes(atomTypes);

    assigner.assign(crystal, typeID, lcs);
    ofstream out("assignment_info.txt");
    assigner.printAssignment(out, crystal, typeID, lcs);

}

void sf_taam_disorder(
    const string &strucureFile,
    std::vector<std::string> substructureFiles,
    const string &jsonFile_taam_multiordered,
    const string &jsonFile_taam_disordered,
    const string &jsonFile_taam_regular,
    const vector<string> &aspherFiles)
{
    Crystal crystal;
    int nSubstructures = substructureFiles.size();

    if (nSubstructures == 0)
        on_error::throwException("zero substructures defined", __FILE__, __LINE__);

    vector<Crystal> substructures(nSubstructures);
    for(int substructureIdx=0; substructureIdx<nSubstructures; substructureIdx++)
        structure_io::read_structure(substructureFiles[substructureIdx], substructures[substructureIdx]);
    structure_io::read_structure(strucureFile, crystal);

    Crystal crystal_buster = crystal;
    for (auto& atom : crystal_buster.atoms)
        atom.label += "   1    A    H105    Z N      1    A X CA     1    A";

    TaamSfCalculatorMultiOrderedImpl tamm_mo_calc(crystal, jsonFile_taam_multiordered);
    HcAtomBankStructureFactorCalculator2 taam_disordered_calc(crystal, jsonFile_taam_disordered);

    nlohmann::json json_data;
    ifstream jsonFileStream(jsonFile_taam_regular);

    if (jsonFileStream.good())
        jsonFileStream >> json_data;
    jsonFileStream.close();
    
    vector<shared_ptr<HcAtomBankStructureFactorCalculator2> > taam_regular_calc(nSubstructures);// (crystal, jsonFile_taam_disordered);
    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
        taam_regular_calc[substructureIdx] = make_shared<HcAtomBankStructureFactorCalculator2>(substructures[substructureIdx], json_data);

    vector<shared_ptr<SfCalculator> > aspherSfCalculators;
    int aspherFileIdx = 1;
    for (auto const& aspherFile : aspherFiles)
    {
        cout << "reading " << aspherFile << "\n";
        nlohmann::json aspherJson;
        ifstream aspherJsonFileStream(aspherFile);
        if (aspherJsonFileStream.good())
            aspherJsonFileStream >> aspherJson;
        aspherJsonFileStream.close();
        if(aspherFileIdx !=3)
            aspherSfCalculators.push_back(shared_ptr<SfCalculator>(SfCalculator::create(crystal, aspherJson)));
        else
        {
            cout << "buster tyoe labels used with " << aspherFile << "\n";
            aspherSfCalculators.push_back(shared_ptr<SfCalculator>(SfCalculator::create(crystal_buster, aspherJson)));
        }
        aspherFileIdx++;
    }

    int nAdditionalCalculators = aspherSfCalculators.size();

   // vector<Vector3i> hkl{ Vector3i(0,0,1), Vector3i(0,0,2), Vector3i(0,2,0) };
    //vector<Vector3i> hkl{ Vector3i(0,2,0) };

    vector<Vector3i> hkl{
        { 0, 0, 1 }, { 0, 0, 2 }, { 0, 0, 3 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 1, 2 },
        { 0, 1, 3 }, { 0, 2, 0 }, { 0, 2, 1 }, { 0, 2, 2 }, { 0, 2, 3 }, { 0, 3, 0 },
        { 0, 3, 1 }, { 0, 3, 2 }, { 0, 3, 3 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 0, 2 },
        { 1, 0, 3 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 1, 2 }, { 1, 1, 3 }, { 1, 2, 0 },
        { 1, 2, 1 }, { 1, 2, 2 }, { 1, 2, 3 }, { 1, 3, 0 }, { 1, 3, 1 }, { 1, 3, 2 },
        { 1, 3, 3 }, { 2, 0, 0 }, { 2, 0, 1 }, { 2, 0, 2 }, { 2, 0, 3 }, { 2, 1, 0 },
        { 2, 1, 1 }, { 2, 1, 2 }, { 2, 1, 3 }, { 2, 2, 0 }, { 2, 2, 1 }, { 2, 2, 2 },
        { 2, 2, 3 }, { 2, 3, 0 }, { 2, 3, 1 }, { 2, 3, 2 }, { 2, 3, 3 }, { 3, 0, 0 },
        { 3, 0, 1 }, { 3, 0, 2 }, { 3, 0, 3 }, { 3, 1, 0 }, { 3, 1, 1 }, { 3, 1, 2 },
        { 3, 1, 3 }, { 3, 2, 0 }, { 3, 2, 1 }, { 3, 2, 2 }, { 3, 2, 3 }, { 3, 3, 0 },
        { 3, 3, 1 }, { 3, 3, 2 }, { 3, 3, 3 }
    };



    //#######################################
    //   structure factors
    //#######################################

    vector<complex<double> > f_mo, f_dis, f_reg;
    vector<vector<complex<double> > > f_part(nSubstructures),f_additional(nAdditionalCalculators);

    cout << "tamm_mo_calc\n";
    tamm_mo_calc.calculateStructureFactors(crystal.atoms, hkl, f_mo, vector<bool>(crystal.atoms.size(), true));
    cout << "taam_disordered_calc\n";
    taam_disordered_calc.calculateStructureFactors(crystal.atoms, hkl, f_dis, vector<bool>(crystal.atoms.size(), true));
    
    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        cout << "taam_regular_calc " << substructureIdx << "\n";
        taam_regular_calc[substructureIdx]->calculateStructureFactors(substructures[substructureIdx].atoms, hkl, f_part[substructureIdx],
            vector<bool>(substructures[substructureIdx].atoms.size(), true));
    }

    //vector<complex
    f_reg = f_part[0];
    int nHkl = hkl.size();
    for (int substructureIdx = 1; substructureIdx < nSubstructures; substructureIdx++)
        for (int i = 0; i < nHkl; i++)
            f_reg[i] += f_part[substructureIdx][i];

    for (int i = 0; i < nAdditionalCalculators; i++)
        aspherSfCalculators[i]->calculateStructureFactors(crystal.atoms, hkl, f_additional[i]);

    cout << "structure factors comparison - methods: reglar vs. disordered\n"
         << agreement_factors::value(f_dis, f_reg) << "\n";
    cout << "structure factors comparison - methods: reglar vs. multi-ordered\n"
        << agreement_factors::value(f_reg, f_mo) << "\n";
    for (int i = 0; i < nAdditionalCalculators; i++)
        cout << "structure factors comparison - methods: reglar vs. method " << i+1 << "\n"
             << agreement_factors::value(f_reg, f_additional[i]) << "\n";

        


    /*
    for (int i = 0; i < nHkl; i++)
    {
        cout << f_mo[i] << " " << f_dis[i] << f_reg[i];
        for (int j = 0; j < nAdditionalCalculators; j++)
            cout << f_additional[j][i] << " ";
        cout << endl;
    }
    */

    //#######################################
    //   structure factor derivatives
    //#######################################

    vector<TargetFunctionAtomicParamDerivatives> dT_dp_reg, dT_dp_mo, dT_dp_dis;
    vector< vector<TargetFunctionAtomicParamDerivatives> > dT_dp_additional(nAdditionalCalculators);
    vector<vector<TargetFunctionAtomicParamDerivatives> > dT_dp_reg_partial(nSubstructures);
    vector<std::complex<double> > dTarget_df(nHkl, 1.0);

    cout << "tamm_mo_calc derivatives\n";
    tamm_mo_calc.calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, f_mo, dT_dp_mo, dTarget_df);
    cout << "taam_disordered_calc derivatives\n";
    taam_disordered_calc.calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, f_dis, dT_dp_dis, dTarget_df);

    for (int i = 0; i < nAdditionalCalculators; i++)
    {
        cout << "derivatives additional calculator " << i+1 << "\n";
        aspherSfCalculators[i]->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, f_dis, dT_dp_additional[i], dTarget_df);
    }


    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        cout << "taam_regular_calc derivatives " << substructureIdx << "\n";
        taam_regular_calc[substructureIdx]->calculateStructureFactorsAndDerivatives(substructures[substructureIdx].atoms, hkl, f_part[substructureIdx],
            dT_dp_reg_partial[substructureIdx], dTarget_df,
            vector<bool>(substructures[substructureIdx].atoms.size(), true));
    }
    
    int atomIdx, nAtoms = crystal.atoms.size();
    TargetFunctionAtomicParamDerivatives dT_dp_0;
    dT_dp_0.atomic_position_derivatives.set(0.0, 0.0, 0.0);
    dT_dp_0.occupancy_derivatives = 0.0;
    dT_dp_reg.assign(nAtoms, dT_dp_0);
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        dT_dp_reg[atomIdx].adp_derivatives.assign(crystal.atoms[atomIdx].adp.size(), 0.0);
    vector<double> total_weight(nAtoms, 0.0);
    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        map<string, int> label2idx_substructure;
        for (int i = 0; i < substructures[substructureIdx].atoms.size(); i++)
            label2idx_substructure[substructures[substructureIdx].atoms[i].label] = i;

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            if (label2idx_substructure.count(crystal.atoms[atomIdx].label) == 0)
                continue;
            int idxPartial = label2idx_substructure[crystal.atoms[atomIdx].label];
            int nAdpComponents = crystal.atoms[atomIdx].adp.size();
            for (int i = 0; i < nAdpComponents; i++)
                dT_dp_reg[atomIdx].adp_derivatives[i] += dT_dp_reg_partial[substructureIdx][idxPartial].adp_derivatives[i];
            dT_dp_reg[atomIdx].atomic_position_derivatives += dT_dp_reg_partial[substructureIdx][idxPartial].atomic_position_derivatives;
            dT_dp_reg[atomIdx].occupancy_derivatives += dT_dp_reg_partial[substructureIdx][idxPartial].occupancy_derivatives *
                substructures[substructureIdx].atoms[idxPartial].occupancy;
            total_weight[atomIdx] += substructures[substructureIdx].atoms[idxPartial].occupancy;
        }
    }
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        dT_dp_reg[atomIdx].occupancy_derivatives /= total_weight[atomIdx];
    // dT_dp_reg, dT_dp_mo, dT_dp_dis;
    bool size_match;
    double rel_l1_100_adp, rel_l1_100_xyz, rel_l1_100_occ, rel_l1_100_fpfdp;
    agreement_factors::for_derivatives(dT_dp_reg, dT_dp_mo, size_match, rel_l1_100_xyz, rel_l1_100_adp, rel_l1_100_occ, rel_l1_100_fpfdp);
    cout<< "derivatives comparison - methods: reglar vs. multi-ordered\n"
        << "  xyz " << rel_l1_100_xyz << "\n"
        << "  adp " << rel_l1_100_adp << "\n"
        << "  occ " << rel_l1_100_occ << "\n";
    agreement_factors::for_derivatives(dT_dp_reg, dT_dp_dis, size_match, rel_l1_100_xyz, rel_l1_100_adp, rel_l1_100_occ, rel_l1_100_fpfdp);
    cout << "derivatives comparison - methods: reglar vs. disordered\n"
        << "  xyz " << rel_l1_100_xyz << "\n"
        << "  adp " << rel_l1_100_adp << "\n"
        << "  occ " << rel_l1_100_occ << "\n";
    for (int i = 0; i < nAdditionalCalculators; i++)
    {
        double af_xyz, af_adp, af_occ, af_fpfdp;
        agreement_factors::for_derivatives(dT_dp_reg, dT_dp_additional[i], size_match, af_xyz, af_adp, af_occ, af_fpfdp);
        
        cout << "derivatives comparison - methods: reglar vs. method " << i + 1 << "\n"
            << "  xyz " << af_xyz << "\n"
            << "  adp " << af_adp << "\n"
            << "  occ " << af_occ << "\n";

    }

    //dT_dp_additional[i]
    //#######################################
    //   atomic form factors
    //#######################################

    vector < vector<complex<double> > > regular_ff, multi_ordered_ff, disordered_ff, substructure_ff;

    tamm_mo_calc.calculateFormFactors(hkl, multi_ordered_ff, vector<bool>(crystal.atoms.size(), true));
    taam_disordered_calc.calculateFormFactors(hkl, disordered_ff, vector<bool>(crystal.atoms.size(), true));
    vector<vector<int> > substructureAtom2StructureAtom(nSubstructures);
    map<string, int> label2idx;
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        label2idx[crystal.atoms[atomIdx].label] = atomIdx;

    regular_ff.assign(nHkl, vector<complex<double> >(nAtoms, 0.0));

    for (int substructureIdx = 0; substructureIdx < nSubstructures; substructureIdx++)
    {
        taam_regular_calc[substructureIdx]->calculateFormFactors(hkl, substructure_ff,
            vector<bool>(substructures[substructureIdx].atoms.size(), true));
        int nAtomsSubstructure = substructures[substructureIdx].atoms.size();
        for (int atomIdxSubstructure = 0; atomIdxSubstructure < nAtomsSubstructure; atomIdxSubstructure++)
        {
            int atomIdx = label2idx[substructures[substructureIdx].atoms[atomIdxSubstructure].label];
            for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
                regular_ff[hklIdx][atomIdx] += substructure_ff[hklIdx][atomIdxSubstructure]*substructures[substructureIdx].atoms[atomIdxSubstructure].occupancy;
        }
    }
    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
            regular_ff[hklIdx][atomIdx] /= total_weight[atomIdx];


    vector<double> atomic_agreement_factors;
    agreement_factors::for_atomic_form_factors(regular_ff, multi_ordered_ff, atomic_agreement_factors);
    auto it = max_element(atomic_agreement_factors.begin(), atomic_agreement_factors.end());
    int idx = distance(atomic_agreement_factors.begin(), it);
    cout << "max agreement factor for atoms (reglar vs. multi-ordered)\n"
        << " " << *it << " for " << crystal.atoms[idx].label << "\n";
    
    agreement_factors::for_atomic_form_factors(regular_ff, disordered_ff, atomic_agreement_factors);
    it = max_element(atomic_agreement_factors.begin(), atomic_agreement_factors.end());
    idx = distance(atomic_agreement_factors.begin(), it);
    cout << "max agreement factor for atoms (reglar vs. disordered)\n"
        << " " << *it << " for " << crystal.atoms[idx].label << "\n";

    for (int i = 0; i < nAdditionalCalculators; i++)
    {
        vector < vector<complex<double> > > additional_ff;
        aspherSfCalculators[i]->calculateFormFactors(hkl, additional_ff, vector<bool>(crystal.atoms.size(), true));

        agreement_factors::for_atomic_form_factors(regular_ff, additional_ff, atomic_agreement_factors);
        it = max_element(atomic_agreement_factors.begin(), atomic_agreement_factors.end());
        idx = distance(atomic_agreement_factors.begin(), it);
        cout << "max agreement factor for atoms (reglar vs. model" << i+1 << ")\n"
            << " " << *it << " for " << crystal.atoms[idx].label << "\n";

    }


}

void test_json()
{
    nlohmann::json data;
    data["list"] = { 1,2,3 };
    for (auto& item : data["list"].items())
        cout << item.key() << " " << item.value() << "\n";
    data["list2"] = R"(
      [ ["C1",0.1],["C2",0.2]]
      )"_json;
    for (auto& item : data["list2"].items())
        cout << item.key() << " " << item.value() << "\n";

    /*
    43, 42, 12
    41, 40, 11

    32, 33, 14
    34, 35, 15

    19, 20, 38, 39
    17, 18, 36, 37
    */
}

string convertBusterLabel(
    const string& label)
{
    string newLabel;
    size_t pos = label.find('.');
    if (pos == string::npos)
        newLabel = label;
    else
    {
        string dotChar = label.substr(pos, 2);
        newLabel = label.substr(0, pos) + label.substr(pos + 2) + dotChar;
    }
    return newLabel;
}

void test()
{
    
    cout<< convertBusterLabel("H2.B   1    A    H105    Z N      1    A X CA     1    A") << "\n";
    cout << convertBusterLabel("H2   1    A    H105    Z N      1    A X CA     1    A") << "\n";
}

void point_group_test()
{
    vector<string> pg_222 = { "x,y,x", "-x,-y,z", "-x,y,-z", "x,-y,-z" };
    vector<string> pg_222_reordered = { "-x,-y,z", "x,y,x", "x,-y,-z", "-x,y,-z" };
    vector<int> canonical_order;
    string pg = crystallographic_point_group_tables::findPointGroup(pg_222, canonical_order);
    cout<< pg << "\n";
    for (int i : canonical_order)
        cout << i << "\n";
    pg = crystallographic_point_group_tables::findPointGroup(pg_222_reordered, canonical_order);
    cout << pg << "\n";
    for (int i : canonical_order)
        cout << i << "\n";

}

template<int N>
int sum_to_n()
{
    int x = 0;
    for (int i = 1; i < N; i++)
        x += i;
    return x;
}

template<>
int sum_to_n<3>()
{
    return 111;
}

void three_way_sf_calculation(
    const Crystal& crystal,
    shared_ptr<SfCalculator>& calculator,
    const  vector<Vector3i>& hkl,
    vector<complex<double> >& sf, 
    vector<complex<double> >& sf_with_derivatives_1,
    vector<complex<double> >& sf_with_derivatives_2)
{
    vector<bool> countAtom(crystal.atoms.size(), true);
    vector<complex<double> > dt_df(hkl.size(), complex<double>(1.0, 0.3));
    vector<TargetFunctionAtomicParamDerivatives> dTarget_dparam;
    discamb::SfDerivativesAtHkl derivatives;

    calculator->calculateStructureFactors(crystal.atoms, hkl, sf);
    calculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkl, sf_with_derivatives_1, dTarget_dparam, dt_df, countAtom);
    sf_with_derivatives_2.resize(hkl.size());
    for (int i = 0; i < hkl.size(); i++)
        calculator->calculateStructureFactorsAndDerivatives(hkl[i], sf_with_derivatives_2[i], derivatives, countAtom);
}

void test_anomalous(
    const string &modelName,
    const Crystal &crystal,
    shared_ptr<SfCalculator> &calculator,
    const  vector<Vector3i> &hkl,
    vector<complex<double> > &d_sf)
{
    vector<complex<double> > sf, sf_anomalous, sf_with_derivatives_1, sf_with_derivatives_1_anomalous, sf_with_derivatives_2, sf_with_derivatives_2_anomalous;

    calculator->setAnomalous(vector<complex<double> >(crystal.atoms.size(), 0.0));

    three_way_sf_calculation(crystal, calculator, hkl, sf, sf_with_derivatives_1, sf_with_derivatives_2);

    vector<complex<double> > df{ {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0.002,0.002}, {0.004,0.003}, {0.008,0.006} };
    vector<complex<double> > anomalous;
    vector<int> atomic_numbers;
    crystal_structure_utilities::atomicNumbers(crystal, atomic_numbers);

    for (int z : atomic_numbers)
        if (z < 9)
            anomalous.push_back(df[z]);
        else
            anomalous.push_back(0.0);
    
    calculator->setAnomalous(anomalous);

    three_way_sf_calculation(crystal, calculator, hkl, sf_anomalous, sf_with_derivatives_1_anomalous, sf_with_derivatives_2_anomalous);

    cout << "test for " << modelName << "\n";
    cout << " no anomalous vs anomalous:\n" 
         << "  " << agreement_factors::value(sf, sf_anomalous) << "\n"
         << "  " << agreement_factors::value(sf_with_derivatives_1, sf_anomalous) << "\n"
         << "  " << agreement_factors::value(sf_with_derivatives_2, sf_anomalous) << "\n"
         << "  " << agreement_factors::value(sf_with_derivatives_2, sf_with_derivatives_1_anomalous) << "\n"
         << "  " << agreement_factors::value(sf_with_derivatives_2, sf_with_derivatives_2_anomalous) << "\n";
    int nHkl = hkl.size();
    d_sf.resize(nHkl);
    for (int i = 0; i < nHkl; i++)
        d_sf[i] = sf[i] - sf_anomalous[i];
}

void test_anomalous(
    const string& structure_file_name)
{
    Crystal crystal;
    structure_io::read_structure(structure_file_name, crystal);

    auto iam = SfCalculator::create_shared_ptr(crystal, {{ "model", "IAM" }});
    auto taam_standard = SfCalculator::create_shared_ptr(crystal, { { "model", "TAAM" }, { "algorithm", "standard"}, {"bank path","MATTS2021databank.txt"}});
    auto taam_macromol = SfCalculator::create_shared_ptr(crystal, { { "model", "TAAM" }, { "algorithm", "macromol"}, {"bank path","MATTS2021databank.txt"} });

    vector<Vector3i> hkl{
    { 0, 0, 1 }, { 0, 0, 2 }, { 0, 0, 3 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 1, 2 },
    { 0, 1, 3 }, { 0, 2, 0 }, { 0, 2, 1 }, { 0, 2, 2 }, { 0, 2, 3 }, { 0, 3, 0 },
    { 0, 3, 1 }, { 0, 3, 2 }, { 0, 3, 3 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 0, 2 },
    { 1, 0, 3 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 1, 2 }, { 1, 1, 3 }, { 1, 2, 0 },
    { 1, 2, 1 }, { 1, 2, 2 }, { 1, 2, 3 }, { 1, 3, 0 }, { 1, 3, 1 }, { 1, 3, 2 },
    { 1, 3, 3 }, { 2, 0, 0 }, { 2, 0, 1 }, { 2, 0, 2 }, { 2, 0, 3 }, { 2, 1, 0 },
    { 2, 1, 1 }, { 2, 1, 2 }, { 2, 1, 3 }, { 2, 2, 0 }, { 2, 2, 1 }, { 2, 2, 2 },
    { 2, 2, 3 }, { 2, 3, 0 }, { 2, 3, 1 }, { 2, 3, 2 }, { 2, 3, 3 }, { 3, 0, 0 },
    { 3, 0, 1 }, { 3, 0, 2 }, { 3, 0, 3 }, { 3, 1, 0 }, { 3, 1, 1 }, { 3, 1, 2 },
    { 3, 1, 3 }, { 3, 2, 0 }, { 3, 2, 1 }, { 3, 2, 2 }, { 3, 2, 3 }, { 3, 3, 0 },
    { 3, 3, 1 }, { 3, 3, 2 }, { 3, 3, 3 }
    };
    int nHkl = hkl.size();
    vector<complex<double> > d_sf_iam, d_sf_taam_standard, d_sf_taam_macromol;

    test_anomalous("IAM", crystal, iam, hkl, d_sf_iam);
    test_anomalous("TAAM standard", crystal, taam_standard, hkl, d_sf_taam_standard);
    test_anomalous("TAAM macromol", crystal, taam_macromol, hkl, d_sf_taam_macromol);

    cout << "anomalous diff - IAM vs. TAAM standard " << agreement_factors::value(d_sf_iam, d_sf_taam_standard) << "\n"
         << "anomalous diff - IAM vs. TAAM macromol " << agreement_factors::value(d_sf_iam, d_sf_taam_macromol) << "\n";
    // no anomalous
    
    //three_way_sf_calculation(crystal, iam, hkl, d_sf_iam);
    //three_way_sf_calculation(crystal, taam_standard, hkl, sf_taam_standard, sf_taam_standard_with_derivatives, sf_taam_standard_with_derivatives2);
    //three_way_sf_calculation(crystal, taam_macromol, hkl, sf_taam_macromol, sf_taam_macromol_with_derivatives, sf_taam_macromol_with_derivatives2);

}

void test_neighbour()
{
    Crystal crystal;
    structure_io::read_structure("urea.res", crystal);
    UnitCellContent unitCellContent(crystal);
    UnitCellContent::AtomID atomId(4, Vector3i(0, 0, -2));
    vector<UnitCellContent::AtomID> centralPart{ atomId };
    vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;
    std::vector< std::vector<discamb::UnitCellContent::AtomID> > molecules;
    structural_properties::splitUnitCellIntoMolecules(unitCellContent, molecules, networkBonds);
    vector<UnitCellContent::AtomID> clusterAtoms;
    structural_properties::makeCluster(unitCellContent, centralPart, molecules, clusterAtoms, 2.1, false);
    for (auto& atom : clusterAtoms)
        cout << atom.atomIndex << " " << atom.unitCellPosition << "n";
}

void test_atom_position()
{
    Crystal crystal;
    structure_io::read_structure("urea.res", crystal);

    for (auto& atom : crystal.atoms)
    {
        Vector3d xyz;
        crystal.unitCell.fractionalToCartesian(atom.coordinates, xyz);
        cout << atom.label << " " << xyz << endl;
    }

    UnitCellContent unitCellContent(crystal);
    int nAtoms = crystal.atoms.size();
    for (int i = 0; i < nAtoms; i++)
    {
        UnitCellContent::AtomID atomId;
        unitCellContent.findAtom(crystal.atoms[i].label, "x,y,z", atomId);
        cout << crystal.atoms[i].label << " " << atomId.atomIndex << " " << atomId.unitCellPosition << endl;
    }
    //UnitCellContent::AtomID atomId(4, Vector3i(0, 0, -2));
    //vector<UnitCellContent::AtomID> centralPart{ atomId };
    //vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;
    //std::vector< std::vector<discamb::UnitCellContent::AtomID> > molecules;
    //structural_properties::splitUnitCellIntoMolecules(unitCellContent, molecules, networkBonds);
    //vector<UnitCellContent::AtomID> clusterAtoms;
    //structural_properties::makeCluster(unitCellContent, centralPart, molecules, clusterAtoms, 2.1, false);
    //for (auto& atom : clusterAtoms)
    //    cout << atom.atomIndex << " " << atom.unitCellPosition << "n";
}

void test_find_neighbour()
{
    Crystal crystal;
    structure_io::read_structure("urea.res", crystal);
    UnitCellContent unitCellContent(crystal);
    UnitCellContent::AtomID atomId(4);
    vector<UnitCellContent::AtomID> centralPart{ atomId }, clusterAtoms;
    structural_properties::makeAtomsCluster(unitCellContent, centralPart, clusterAtoms, 2.1, false);
    for (auto& atom : clusterAtoms)
        cout << atom.atomIndex << " " << atom.unitCellPosition << endl;
}

void get_representatives()
{
    Crystal crystal;
    structure_io::read_structure("urea.res", crystal);
    nlohmann::json data;
    data["model"] = "HAR";
    std::vector<QmFragmentInCrystal> crystalFragments;
    ham_settings::setCrystalFragments(data, crystal, crystalFragments);
    vector<vector<AtomRepresentativeInfo> > repesentatives;
    vector<vector<pair<string, string> > > subsystemAtoms(1);
    subsystemAtoms[0] = {
        {"C","x,y,z"},
        {"O","x,y,z"},
        {"N","x,y,z"},
        {"N","-x,-y+1,z"},
        {"Hb","x,y,z"},
        {"Ha","x,y,z"},
        {"Hb","-x,-y+1,z"},
        {"Ha","-x,-y+1,z"}
    };
    gar_utilities::findDefaultRepresentatives(
        crystal, subsystemAtoms, repesentatives);


    
    ham_settings::setRepresentatives(data, crystal, crystalFragments, repesentatives);
    UnitCellContent unitCellContent(crystal);
}

void frozen_lcs()
{
    
    nlohmann::json aspher_frozen = nlohmann::json::parse(R"(
  {
    "model": "taam",
    "bank path":"MATTS2021databank.txt",
    "algorithm": "macromol",
    "frozen lcs": true
  }
)");

    nlohmann::json aspher_no_frozen = nlohmann::json::parse(R"(
  {
    "model": "taam",
    "bank path":"MATTS2021databank.txt",
    "algorithm": "macromol",
    "frozen lcs": false
  }
)");

    Crystal crystal, crystal_rot;
    structure_io::read_structure("l_ala.res", crystal);
    crystal_rot = crystal;
    Vector3d n, h1, h2, h3, frac, v1, v2, v3, d1, d2, d3;
    crystal.unitCell.fractionalToCartesian(crystal.atoms[0].coordinates, n);
    crystal.unitCell.fractionalToCartesian(crystal.atoms[1].coordinates, h1);
    crystal.unitCell.fractionalToCartesian(crystal.atoms[2].coordinates, h2);
    crystal.unitCell.fractionalToCartesian(crystal.atoms[3].coordinates, h3);
    v1 = h1 - n;
    v2 = h2 - n;
    v3 = h3 - n;
    
    d1 = v1 + v2;
    d1 = 1.02 * d1 / sqrt(d1 * d1);
    d2 = v2 + v3;
    d2 = 1.02 * d2 / sqrt(d2 * d2);
    d3 = v3 + v1;
    d3 = 1.02 * d3 / sqrt(d3 * d3);

    h1 = n + d1;
    h2 = n + d2;
    h3 = n + d3;

    crystal.unitCell.cartesianToFractional(h1, crystal_rot.atoms[1].coordinates);
    crystal.unitCell.cartesianToFractional(h2, crystal_rot.atoms[2].coordinates);
    crystal.unitCell.cartesianToFractional(h3, crystal_rot.atoms[3].coordinates);

    auto fcal_frozen = SfCalculator::create_shared_ptr(crystal, aspher_frozen);
    auto fcal_no_frozen = SfCalculator::create_shared_ptr(crystal, aspher_no_frozen);

    structure_io::write_structure("l_ala_rot.res", crystal_rot);

    vector<Vector3i> hkl{
    { 0, 0, 1 }, { 0, 0, 2 }, { 0, 0, 3 }, { 0, 1, 0 }, { 0, 1, 1 }, { 0, 1, 2 },
    { 0, 1, 3 }, { 0, 2, 0 }, { 0, 2, 1 }, { 0, 2, 2 }, { 0, 2, 3 }, { 0, 3, 0 },
    { 0, 3, 1 }, { 0, 3, 2 }, { 0, 3, 3 }, { 1, 0, 0 }, { 1, 0, 1 }, { 1, 0, 2 },
    { 1, 0, 3 }, { 1, 1, 0 }, { 1, 1, 1 }, { 1, 1, 2 }, { 1, 1, 3 }, { 1, 2, 0 },
    { 1, 2, 1 }, { 1, 2, 2 }, { 1, 2, 3 }, { 1, 3, 0 }, { 1, 3, 1 }, { 1, 3, 2 },
    { 1, 3, 3 }, { 2, 0, 0 }, { 2, 0, 1 }, { 2, 0, 2 }, { 2, 0, 3 }, { 2, 1, 0 },
    { 2, 1, 1 }, { 2, 1, 2 }, { 2, 1, 3 }, { 2, 2, 0 }, { 2, 2, 1 }, { 2, 2, 2 },
    { 2, 2, 3 }, { 2, 3, 0 }, { 2, 3, 1 }, { 2, 3, 2 }, { 2, 3, 3 }, { 3, 0, 0 },
    { 3, 0, 1 }, { 3, 0, 2 }, { 3, 0, 3 }, { 3, 1, 0 }, { 3, 1, 1 }, { 3, 1, 2 },
    { 3, 1, 3 }, { 3, 2, 0 }, { 3, 2, 1 }, { 3, 2, 2 }, { 3, 2, 3 }, { 3, 3, 0 },
    { 3, 3, 1 }, { 3, 3, 2 }, { 3, 3, 3 }
    };
    vector<complex<double> > sf_frozen, sf_not_frozen;
    fcal_frozen->calculateStructureFactors(crystal_rot.atoms, hkl, sf_frozen);
    fcal_no_frozen->calculateStructureFactors(crystal_rot.atoms, hkl, sf_not_frozen);

    cout << "r-factor frozen vs no frozen " << agreement_factors::value(sf_frozen, sf_not_frozen) << endl;
    
    vector<TargetFunctionAtomicParamDerivatives> dt_dp, dt_dp2;
    vector<complex<double> > dt_df(hkl.size(),1.0), f1, f2;

    fcal_frozen->calculateStructureFactorsAndDerivatives(crystal_rot.atoms, hkl, f1, dt_dp, dt_df, vector<bool>(crystal.atoms.size(), true));
    crystal_rot.atoms[0].coordinates += Vector3d(0.000001, 0.0, 0.0);
    fcal_frozen->calculateStructureFactorsAndDerivatives(crystal_rot.atoms, hkl, f2, dt_dp2, dt_df, vector<bool>(crystal.atoms.size(), true));

    
}

void atomTypesDisorder()
{
    CrystalAtomTypeAssigner assigner;
    vector < vector <pair<int, string> > > fragmentAtoms(1);
    fragmentAtoms[0].push_back({ 0, "x,y,z" });
    fragmentAtoms[0].push_back({ 1, "x,y,z" });
    fragmentAtoms[0].push_back({ 1, "-X+1/2,-Y+1/2,Z" });

    vector< vector<int> > atomsToAssign{ {0, 1, 2} };
    vector< vector<int> > typeID(1);
    vector< vector<LocalCoordinateSystem<AtomInCrystalID> > > lcs;

    MATTS_BankReader bankReader;
    vector<AtomType> types;
    vector<AtomTypeHC_Parameters> hcParameters;
    BankSettings bankSettings;
    bankReader.read("MATTS2021databank.txt", types, hcParameters, bankSettings);
    assigner.setAtomTypes(types);

    Crystal crystal;
    structure_io::read_structure("ice.res", crystal);

    assigner.assign(
            crystal,
            fragmentAtoms,
            atomsToAssign,
            typeID,
            lcs);

    for (auto& ids : typeID[0])
        if (ids >= 0)
            cout << types[ids].id << endl;
        else
            cout << "unassigned" << endl;

}

int main(int argc, char* argv[])
{
    try {

        atomTypesDisorder();
        return 0;

        frozen_lcs();
        return 0;

        get_representatives();
        return 0;

        test_atom_position();
        return 0;

        test_neighbour();
        return 0;

        test_anomalous(argv[1]);
        if (argc != 2)
            on_error::throwException("expected structure file an an argument\n", __FILE__, __LINE__);

        return 0;

        vector<string> arguments, options;
        parse_cmd::get_args_and_options(argc, argv, arguments, options);
        
        if (arguments.size() != 2)
            on_error::throwException("expected structure file, hkl indices file and optionally -n - for new version calculation\n", __FILE__, __LINE__);
        
        bool onlyNew = parse_cmd::hasOption(options, "-n");
        bool sf_and_sfAndDf = !parse_cmd::hasOption(options, "-d");
        bool anomalous = parse_cmd::hasOption(options, "-a");
        
        taam_parallel(argv[1], argv[2], onlyNew, sf_and_sfAndDf, anomalous);

#ifdef _MSC_VER
        _CrtSetReportMode(_CRT_WARN, _CRTDBG_MODE_DEBUG);
        _CrtDumpMemoryLeaks();
#endif
        return 0;


        check_args(argc, argv, 3, { "structure file", "indices file", "json file"});
        make_hkl(argv[1], argv[2], argc > 3 ? argv[3] : "");
        return 0;

        check_args(argc, argv, 2, { "fcf file", "structure file" });
        compare_olex_taam(argv[1], argv[2]);
        return 0;

        check_args(argc, argv, 2, { "hkl indices", "n sets" });
        sort_split_hkl(argv[1], atoi(argv[2]));
        return 0;

        
        cout << sum_to_n<0>() << endl;
        cout << sum_to_n<1>() << endl;
        cout << sum_to_n<2>() << endl;
        cout << sum_to_n<3>() << endl;
        cout << sum_to_n<4>() << endl;

        return 0;

        point_group_test();
        return 0;

        //test();
        //return 0;

        string commandLineArgsInfo =
            "(1) structure file\n" 
            "(2) substructure files as file1,file2,.. \n" 
            "(3) json file for TAAM multiordered\n" 
            "(4) json file for TAAM disordered\n" 
            "(5) json file for TAAM regular\n" 
            "(6) optionally json files for model defined in the file as file1,file2,..\n";
        if (argc < 6)
            on_error::throwException("expected arguments:\n" + commandLineArgsInfo, __FILE__, __LINE__);

        int nSubstructureFiles = argc - 5;

        string strucureFile = argv[1];
        string substructureFileList = argv[2];
        vector<string> substructureFiles;
        string_utilities::split(substructureFileList, substructureFiles, ',');
        string jsonFile_taam_multiordered = argv[3];
        string jsonFile_taam_disordered = argv[4];
        string jsonFile_taam_regular = argv[5];
        vector<string> aspherFiles;

        if (argc == 7)
        {
            string aspherFilesList = argv[6];
            string_utilities::split(aspherFilesList, aspherFiles, ',');
        }

        sf_taam_disorder(
            strucureFile,
            substructureFiles,
            jsonFile_taam_multiordered,
            jsonFile_taam_disordered,
            jsonFile_taam_regular,
            aspherFiles);
        return 0;

        test_json();
        return 0;

        check_args(argc, argv, 3, { "structure file", "fragments file or -ordered", "bank file" });
        string arg2 = argv[2];
        if (arg2 == "-ordered")
            type_assign_order(argv[1], argv[3]);
        else
            type_assign_disorder(argv[1], argv[2], argv[3]);
        return 0;


        //convert_sph_harmonics_test();
        //return 0;

        copy_compare();
        return 0;

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

