#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/TscFileBasedSfCalculator.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include <fstream>

using namespace std;
using namespace discamb;


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

struct TscComparison {
    double rf = 0.0;
    vector<double> atomRf;
    vector<double> chargeDiff;
};

void compareFfAndSf(
    const vector<vector<complex<double> > > &ff1,
    const vector<vector<complex<double> > >& ff2,
    const vector<complex<double> > &sf1,
    const vector<complex<double> >& sf2,
    TscComparison &comparison)
{
    comparison.rf = agreementFactor(sf1, sf2);
    int nAtoms = ff1[0].size();
    comparison.atomRf.clear();
    for (int i = 0; i < nAtoms; i++)
    {
        vector<complex<double> > ff1_atom, ff2_atom;
        for (int j = 0; j < ff1.size(); j++)
        {
            ff1_atom.push_back(ff1[j][i]);
            ff2_atom.push_back(ff2[j][i]);
        }
        comparison.atomRf.push_back(agreementFactor(ff1_atom, ff2_atom));
    }
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


void compare_tsc(
    const string& structureFile,
    const string& reference_tsc,
    const vector<string> &_otherTsc)
{
    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);
    vector<vector<complex<double> > > formFactors, formFactorsRef;
    vector<complex<double> > structureFactors, structureFactorsRef;
    vector<Vector3i> hkl, hkl0;
    vector<string> otherTsc, atomLabels;
    vector<TscComparison> tscComparison;
    
    if (_otherTsc.empty())
    {
        file_system_utilities::find_files("tsc", otherTsc);
        if (otherTsc.size() < 2)
            on_error::throwException("no tsc files to compare with the reference data", __FILE__, __LINE__);
        auto it = find(otherTsc.begin(), otherTsc.end(), reference_tsc);
        if (it != otherTsc.end())
            otherTsc.erase(it);
    }
    else
        otherTsc = _otherTsc;

    tsc_io::read_tsc(reference_tsc, atomLabels, hkl, formFactorsRef);
    map<Vector3i, vector<complex<double> > > formFactorsMap;
    for(int i = 0; i < hkl.size(); i++)
        formFactorsMap[hkl[i]] = formFactorsRef[i];
    shared_ptr<AtomicFormFactorCalculationsManager> ff_calc = make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, formFactorsMap);
    AnyScattererStructureFactorCalculator sfCalculator(crystal);
    sfCalculator.setAtomicFormfactorManager(ff_calc);
    vector<bool> countAtom(crystal.atoms.size(), true);
    sfCalculator.calculateStructureFactors(hkl, structureFactorsRef);
    

    for (auto& tscFileName : otherTsc)
    {
        tsc_io::read_tsc(tscFileName, atomLabels, hkl, formFactors);
        formFactorsMap.clear();
        for (int i = 0; i < hkl.size(); i++)
            formFactorsMap[hkl[i]] = formFactors[i];
        std::dynamic_pointer_cast<ConstFormFactorCalculationsManager>(ff_calc)->resetFormFactors(formFactorsMap);
        sfCalculator.calculateStructureFactors(hkl, structureFactors);
        TscComparison tscComp;
        compareFfAndSf(formFactorsRef, formFactors, structureFactorsRef, structureFactors, tscComp);
        tscComparison.push_back(tscComp);
    }

    ofstream out("tsc_comparison");
    
    out << "overall R-factor\n";
    for (int i=0;i<otherTsc.size();i++)
    {
        out << left << setw(20) << otherTsc[i] << setprecision(3) << fixed << tscComparison[i].rf << "\n";
    }
    out << "atomic r-factors\n" << setw(20) << " ";
    for (int i = 0; i < crystal.atoms.size(); i++)
        out << setw(14) << crystal.atoms[i].label;
    out << endl;
    for (int i = 0; i < otherTsc.size(); i++)
    {
        out << left << setw(20) << otherTsc[i] << setprecision(3);
        for (int j = 0; j < crystal.atoms.size(); j++)
            out << setw(14) <<  fixed << tscComparison[i].atomRf[j];
        out << endl;

    }

    out << "\n\nSORTED\n\n";

    vector<pair<double, string> > rf_sorted;
    for (int i = 0; i < otherTsc.size(); i++)
        rf_sorted.push_back(make_pair(tscComparison[i].rf, otherTsc[i]));
    sort(rf_sorted.begin(), rf_sorted.end());

    for (auto item: rf_sorted)
        out << left << setw(20) << item.second << setprecision(3) << fixed << item.first << "\n";


    out << "atomic r-factors\n" << setw(20) << " ";
    for (int i = 0; i < crystal.atoms.size(); i++)
        out << setw(14) << crystal.atoms[i].label;
    out << endl;


    vector<pair<double, int> > h_atom_rf_sorted;
    for (int i = 0; i < otherTsc.size(); i++)
    {
        double average = 0;
        int n = 0;
        for (int j = 0; j < crystal.atoms.size(); j++)
            if (crystal.atoms[j].type == "H")
            {
                average += tscComparison[i].atomRf[j];
                n++;
            }
        if (n > 0)
            average /= n;
        h_atom_rf_sorted.push_back(make_pair(average, i));
    }
    sort(h_atom_rf_sorted.begin(), h_atom_rf_sorted.end());
    for (int p = 0; p < otherTsc.size(); p++)
    {
        int i = h_atom_rf_sorted[p].second;
        out << left << setw(20) << otherTsc[i] << setprecision(3);
        for (int j = 0; j < crystal.atoms.size(); j++)
            out << setw(14) << fixed << tscComparison[i].atomRf[j];
        out << endl;

    }

    out.close();
}




int main(int argc, char* argv[])
{
    try {

        if (argc < 3)
        {
            cout << "expected at least 2 arguments:\n"
                "   (1) structure file\n"
                "   (2) ref tsc file\n"
                "and optionally list of tsc to compare\n";
            exit(0);
        }

        vector<string> otherTsc;
        for (int i = 3; i < argc; i++)
            otherTsc.push_back(argv[i]);
        compare_tsc(argv[1], argv[2], otherTsc);
        
        return 0;



    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

