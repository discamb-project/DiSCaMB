#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/discamb_io.h"
#include "discamb/IO/fragmentation_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/QuantumChemistry/fragmentation.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/TscFileBasedSfCalculator.h"
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

double agreementFactor(
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
            " (3) tsc file for model 1\n";
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

int main(int argc, char* argv[])
{
    try { 
     
        compare_sf_taam_qm(argc, argv);
        return 0;
        test_taam_disorder(argc, argv);

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

