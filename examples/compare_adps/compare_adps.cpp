#include "discamb/BasicChemistry/PeriodicTable.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/StructuralParametersConverter.h"
#include "discamb/IO/structure_io.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include "json.hpp"

#include <fstream>
#include <filesystem>
#include <iomanip>

using namespace std;

using namespace discamb;
struct AdpDescriptors {
    double eta_r;
    double msdCorrelation;
    double ratioUeq;
};

void compareAdps(
    const vector<double>& a,
    const vector<double>& a_ref,
    const UnitCell& unitCell,
    const UnitCell& unitCellRef,
    AdpDescriptors & adpDescriptors)
{
    vector<double> ac(6), ac_ref(6);
    StructuralParametersConverter converter, converterRef;
    converter.set(unitCell);
    converterRef.set(unitCellRef);
    converter.convertADP(a, ac, structural_parameters_convention::AdpConvention::U_cif, structural_parameters_convention::AdpConvention::U_cart);
    converterRef.convertADP(a_ref, ac_ref, structural_parameters_convention::AdpConvention::U_cif, structural_parameters_convention::AdpConvention::U_cart);

    adpDescriptors.eta_r = 100*(1.0 - crystal_structure_utilities::overlappingAdpSimilarityIndex(ac_ref, ac));
    adpDescriptors.msdCorrelation = crystal_structure_utilities::msdCorrelation(ac_ref, ac);
    adpDescriptors.ratioUeq = (ac[0] + ac[1] + ac[2]) / (ac_ref[0] + ac_ref[1] + ac_ref[2]);
}



ostream& operator<<(ostream& os, const AdpDescriptors& adpDescriptors)
{
    int precision = os.precision();
    auto formatFlags = os.flags();
    os << setprecision(7) << fixed;

    os  << setw(14) << adpDescriptors.eta_r
        << setw(14) << adpDescriptors.msdCorrelation
        << setw(14) << adpDescriptors.ratioUeq;

    os.precision(precision);
    os.setf(formatFlags);

    return os;
}



void calculateAdpGroupStats(
    const vector<vector<AdpDescriptors> >& adpDescriptors,
    const vector<int>& atomsInGroup,
    vector<AdpDescriptors>& groupDescriptors)
{
    if (adpDescriptors.empty())
        return;

    int structureIdx, nStructures = adpDescriptors[0].size();
    int nAtomsInGroup = atomsInGroup.size();
    
    groupDescriptors.resize(nStructures);

    for (structureIdx = 0; structureIdx < nStructures; structureIdx++)
    {
        groupDescriptors[structureIdx].ratioUeq = 0;
        groupDescriptors[structureIdx].eta_r = 0;
        groupDescriptors[structureIdx].msdCorrelation = 0;


        for (int atomIdx: atomsInGroup)
        {
            AdpDescriptors const& d = adpDescriptors[atomIdx][structureIdx];
            groupDescriptors[structureIdx].ratioUeq += d.ratioUeq;
            groupDescriptors[structureIdx].eta_r += d.eta_r;
            groupDescriptors[structureIdx].msdCorrelation += d.msdCorrelation;
        }
        groupDescriptors[structureIdx].ratioUeq /= nAtomsInGroup;
        groupDescriptors[structureIdx].eta_r /= nAtomsInGroup;
        groupDescriptors[structureIdx].msdCorrelation /= nAtomsInGroup;
    }

}


void calculateAdpStatsAtomByAtom(
    const Crystal& refCrystal,
    const vector<Crystal>& crystals,
    vector<vector<AdpDescriptors> >& adpDescriptors)
{
    int structureIdx, nStructures, atomIdx, nAtoms;

    nAtoms = refCrystal.atoms.size();
    nStructures = crystals.size();
    adpDescriptors.clear();
    adpDescriptors.resize(nAtoms);
    vector<double> u, u_ref;

    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {

        adpDescriptors[atomIdx].resize(nStructures);

        if (refCrystal.atoms[atomIdx].adp.size() == 6)
        {
            vector<double> const& u_ref = refCrystal.atoms[atomIdx].adp;

            for (structureIdx = 0; structureIdx < nStructures; structureIdx++)
                compareAdps(crystals[structureIdx].atoms[atomIdx].adp, u_ref, 
                            crystals[structureIdx].unitCell, refCrystal.unitCell, 
                            adpDescriptors[atomIdx][structureIdx]);
        }
    }

}


void printStats(
    const Crystal& crystal,
    const string& fileName,
    const vector<string>& structureNames,
    const vector<vector<AdpDescriptors> > &adpDescriptors,
    const vector<AdpDescriptors> &hAdpsStats,
    const vector<AdpDescriptors> &nonhAdpsStats,
    const map<int, vector<AdpDescriptors> > hAdpDescriptorsByNeighbourElement)
{
    ofstream out(fileName);
    size_t i, maxNameLength;
    vector<string> names = structureNames;
    vector<size_t> nameSize;
    for (auto& name : names)
        nameSize.push_back(name.size());


    maxNameLength = *max_element(nameSize.begin(), nameSize.end());

    out << "     compare_adps v. 0.1           \n"
        << "                                   \n"
        << "  The program compares ADPs in terms of         \n"
        << "     eta_r        - rescaled overlapping coefficient\n"
        << "     msd_corr     - correlation of atomic mean square displacements\n"
        << "     Ueq/Ueq(Ref) - ratio of equivalent isotropic atomic displacement parameters, \n"
        << "                    reference value in denominator\n"
        << "                                   \n"
        << "  References:                      \n"
        << "      DiSCaMB:                         \n"
        << "          Chodkiewicz, M.L., Migacz, S., Rudnicki, W., Makal, A., Kalinowski, J.A.,\n"
        << "          Moriarty, N.W., Grosse - Kunstleve, R.W., Afonine, P.V., Adams, P.D.&\n"
        << "          Dominiak, P.M. (2018).J.Appl.Cryst. 51, 193 - 199. doi: 10.1107/S1600576717015825\n"
        << "      rescaled overlapping coefficient: \n"
        << "          Chodkiewicz, M., Patrikeev, L., Pawledzio, S. & Wozniak, K. (2024).\n"
        << "          IUCrJ, 11, 249-259. doi: 10.1107/S2052252524001507\n"
        << "      correlation of atomic mean square displacements: \n"
        << "          Chodkiewicz, M. & Wozniak. K (2024)\n"
        << "                                   \n"
        << "  The program is a part of DiSCaMB library distribution\n"
        << "  copyrights: University of Warsaw\n"
        << "  All right reserved              \n\n\n";

    out << "content:\n\n"
        << "  (1) averaged H atom ADPs statistics\n"
        << "      (1.1) all H atoms\n"
        << "      (1.2) grouped by bonding element\n"
        << "  (2) averaged non-H atom ADPs statistics\n"
        << "  (3) case by case ADPs comparison\n";



    //---------------------------------------------------------------------------------------------

    out << "\n\n(1) averaged H atom ADPs statistics\n\n";
    
    out << "      (1.1) all H atoms\n";

    out << string(maxNameLength + 4, ' ') << "      eta_r        msd_corr    Ueq/Ueq(Ref)\n";
    //                                       "              .             .             .             .
    for (i = 0; i < hAdpsStats.size(); i++)
        out << setw(maxNameLength + 4) << names[i] << hAdpsStats[i] << "\n";
     
    out << "\n      (1.2) staistics grouped by H bonding element\n\n";

    out << string(maxNameLength + 4, ' ') << "      eta_r        msd_corr    Ueq/Ueq(Ref)\n";

    for (auto &item: hAdpDescriptorsByNeighbourElement)
    {
        auto& stats = item.second;
        if(item.first < 0)
            out << "  H-polar\n";
        else
            out << "  H-" << periodic_table::symbol(item.first) << "\n";
        for (i = 0; i < hAdpsStats.size(); i++)
            out << setw(maxNameLength + 4) << names[i] << stats[i] << endl;

    }

    out << "\n\n(2) averaged non-H atom ADPs statistics\n\n";
                                             
    out << string(maxNameLength + 4, ' ') << "      eta_r        msd_corr    Ueq/Ueq(Ref)\n";
    for (i = 0; i < nonhAdpsStats.size(); i++)
        out << setw(maxNameLength + 4) << names[i] << nonhAdpsStats[i] << "\n";


    out << "\n\n(3) case by case ADPs comparison\n\n";

    out << string(maxNameLength + 4, ' ') << "      eta_r        msd_corr    Ueq/Ueq(Ref)\n";

    size_t atomIdxInCrystal, atomIdxInGroup;
    for (int atomIdx = 0; atomIdx < adpDescriptors.size(); atomIdx++)
    {

        out << crystal.atoms[atomIdx].label << "\n";
        for (i = 0; i < adpDescriptors[atomIdx].size(); i++)
            out << setw(maxNameLength+4) << names[i] << adpDescriptors[atomIdx][i] << endl;
    }


    out.close();
}

nlohmann::json adpStats2Json(const vector<AdpDescriptors>& d)
{
    vector<double> u_eq_ratio, eta_r, msd_corr;

    for (int i = 0; i < d.size(); i++)
    {
        eta_r.push_back(d[i].eta_r);
        msd_corr.push_back(d[i].msdCorrelation);
        u_eq_ratio.push_back(d[i].ratioUeq);
    }

    nlohmann::json data = { {"eta_r", eta_r}, {"MSD_corr", msd_corr}, {"Ueq ratio", u_eq_ratio}};

    return data;
}


void printStatsJson(
    const Crystal& crystal,
    const string& fileName,
    const vector<string>& structureNames,
    const vector<vector<AdpDescriptors> >& adpDescriptors,
    const vector<AdpDescriptors>& hAdpsStats,
    const vector<AdpDescriptors>& nonhAdpsStats,
    const map<int, vector<AdpDescriptors> > hAdpDescriptorsByNeighbourElement)
{
    nlohmann::json jsonData;
    
    jsonData["content"] = "structure comparison";
    jsonData["refinement labels"] = structureNames;


    auto& adpData = jsonData["ADPs"];

    auto& hAdpData = adpData["H"];

    hAdpData["average"] = adpStats2Json(hAdpsStats);

    auto& adpsByBondedElement = hAdpData["by neighbour element"];

    for (auto& item : hAdpDescriptorsByNeighbourElement)
    {
        string symbol = (item.first > 0 ? periodic_table::symbol(item.first) : string("polar"));
        adpsByBondedElement[symbol] = adpStats2Json(item.second);
    }

    adpData["non-H"] = adpStats2Json(nonhAdpsStats);

    ofstream out(fileName);
    out << jsonData.dump(4);
    out.close();
}

 
int main(int argc, char* argv[])
{
    try { 
     
        string help = "usage:\n"
                      "  compare_adps reference_structure structure_1 ..  structure_n output_file\n\n"
                      "structure files should have one of the following extensions: cif, res, ins\n";

        if (argc < 4)
        {
            if (argc == 1)
                cout << "error - missing arguents," << help << "\n\n";
            else
                on_error::throwException("Expected at least 3 arguments, " + help, __FILE__, __LINE__);
        }

        string outputFileName = argv[argc-1];

        Crystal refCrystal;
        vector<Crystal> crystals;
        vector<string> structureNames;
        structure_io::read_structure(argv[1], refCrystal);
        for (size_t i = 2; i < argc - 1; i++)
        {
            crystals.resize(crystals.size() + 1);
            structure_io::read_structure(argv[i], crystals.back());
            structureNames.push_back(argv[i]);
        }

        vector<vector<pair<int, string> > > bonds;
        structural_properties::asymmetricUnitConnectivity(refCrystal, bonds, 0.4);

        vector<vector<AdpDescriptors> > adpDescriptors;
        vector<AdpDescriptors> hAdpsStats;
        vector<AdpDescriptors> nonhAdpsStats;
        vector<int> atomicNumbers;
      

        crystal_structure_utilities::atomicNumbers(refCrystal, atomicNumbers);

        calculateAdpStatsAtomByAtom(refCrystal, crystals, adpDescriptors);


        vector<int> hAtoms, nonHAtoms;
        // h - atoms grouped by neighbour element
        map<int, vector<int> > hAtomsByNeighbourElement;

        for (int i=0; i < atomicNumbers.size(); i++)
        {
            if (atomicNumbers[i] == 1)
            {
                hAtoms.push_back(i);
                for (auto& bond : bonds[i])
                    hAtomsByNeighbourElement[atomicNumbers[bond.first]].push_back(i);
            }
            else
                nonHAtoms.push_back(i);
        }

        calculateAdpGroupStats(adpDescriptors, hAtoms, hAdpsStats);
        calculateAdpGroupStats(adpDescriptors, nonHAtoms, nonhAdpsStats);

        map<int, vector<AdpDescriptors> > hAdpDescriptorsByNeighbourElement;
        map<int, int> hAdpByNeighbourElementCount;


        for (auto& item : hAtomsByNeighbourElement)
            calculateAdpGroupStats(adpDescriptors, item.second, hAdpDescriptorsByNeighbourElement[item.first]);

        printStats(refCrystal, outputFileName, structureNames, adpDescriptors,
            hAdpsStats, nonhAdpsStats, hAdpDescriptorsByNeighbourElement);

        printStatsJson(refCrystal, outputFileName + ".json", structureNames,  adpDescriptors,
            hAdpsStats, nonhAdpsStats, hAdpDescriptorsByNeighbourElement);

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

