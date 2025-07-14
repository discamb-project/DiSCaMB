#include "discamb/HC_Model/SlaterOrbitalWfnData.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicChemistry/periodic_table.h"



using namespace std;

namespace wfndata
{

// max subshell D
void stringToConfiguration(
    const string s,
    vector<vector<int> > &configuration)
{
    int i, shell, subshell, nElectrons;
    vector<string> words;

    configuration.clear();

    discamb::string_utilities::split(s, words, ')');

    // K(2 L(8 3S(2 3P(6 4S(0 3D(4 

    for (i = 0; i < words.size(); i++)
    {
        // shell e.g. K(2)
        if (!isdigit(words[i][0]))
        {
            if (words[i][0] == 'K')
                configuration.resize(1, vector<int>(1, 2));
            if (words[i][0] == 'L')
            {
                configuration.resize(2, vector<int>(2, 2));
                configuration[1][1] = 6;
            }
            if (words[i][0] == 'M')
            {
                configuration.resize(3, vector<int>(3, 2));
                configuration[2][1] = 6;
                configuration[2][2] = 10;
            }

        }
        else
        {
            shell = atoi((string() + words[i][0]).c_str());
            if (words[i][1] == 'S')
                subshell = 0;
            else
            {
                if (words[i][1] == 'P')
                    subshell = 1;
                else
                {
                    if (words[i][1] == 'D')
                        subshell = 2;
                    else
                        discamb::on_error::throwException("error in Clementi-Roetti data code", __FILE__, __LINE__);
                    //{
                    //    cout << "ERROR" << endl;
                    //    exit(1);
                    //}
                }

            }

            nElectrons = atoi(words[i].substr(3).c_str());

            if (nElectrons > 0)
            {
                if (shell > configuration.size())
                    configuration.resize(shell);
                if (subshell + 1 > configuration[shell - 1].size())
                    configuration[shell - 1].resize(subshell + 1);
                configuration[shell - 1][subshell] = nElectrons;
            }

        }

    }
}

void setConfigurationData(
    map<string, vector<vector<int> > > &configurations,
    string configuration_data[])
{
    int i;
    vector<string> words;
    vector<vector<int> > configuration;
    //configuration.resize(68);

    //"Mn3+ K(2)L(8)3S(2)3P(6)4S(0)3D(4)"

    for (i = 0; i < 68; i++)
    {
        discamb::string_utilities::split(configuration_data[i], words);
        stringToConfiguration(words[1], configuration);
        configurations[words[0]] = configuration;
    }
}

void valenceOrbitals(
    LocalDef_Wfn &wfn,
    const vector<vector<int> > &configuration,
    vector<pair<int, int> > &val_orb)
{

    int _no_valence[] = { 0, 2, 10, 18, 36, 54 };
    vector<int> no_valence(_no_valence, _no_valence + 6);
    int n_electrons;
    int outerShell = 1;
    int period = discamb::periodic_table::period(wfn.atomicNumber);
    discamb::periodic_table::Block block = discamb::periodic_table::block(wfn.atomicNumber);

    n_electrons = wfn.atomicNumber - wfn.charge; 

    val_orb.clear();

    if (find(no_valence.begin(), no_valence.end(), n_electrons) != no_valence.end())
        return;

    if (wfn.charge > 0 && (n_electrons == 28 || n_electrons == 46))
        return;

    if (block == discamb::periodic_table::Block::S || block == discamb::periodic_table::Block::P)
        val_orb.push_back(make_pair(period, 0));

    if (block == discamb::periodic_table::Block::P)
        val_orb.push_back(make_pair(period, 1));


    // D block elements with 10 d electrons - no valence orbitals

    if (block == discamb::periodic_table::Block::D)
        if (configuration.size()>2)
            if (configuration[2].size() > 2)
            {
                if (configuration[2][2] != 10)
                    val_orb.push_back(make_pair(3, 2));
                else
                    if (configuration.size()>3)
                        val_orb.push_back(make_pair(4, 0));
            }


}

double factorial(int n)
{
    double result = 1;

    for (int i = 1; i <= n; i++)
        result *= i;
    return result;
}

double slaterOverlap(int n1, double a1, int n2, double a2)
{
    int n = n1 + n2 + 2;
    double a = a1 + a2;
    double f = factorial(n);
    double p = pow(a, n + 1);
    return f / p;
}

void renormalizeOrbital(
    discamb::SlaterTypeAtomicOrbitalRdf& orbital)
{
    int nFunctions = orbital.coefficients.size();
    double overlap = 0;
    for (int i = 0; i < nFunctions; i++)
        for (int j = 0; j < nFunctions; j++)
            overlap += orbital.coefficients[i] * orbital.coefficients[j] *
            slaterOverlap(orbital.power_r[i], orbital.exponents[i], orbital.power_r[j], orbital.exponents[j]);

    double multiplier = sqrt(1.0 / overlap);

    for (int i = 0; i < nFunctions; i++)
        orbital.coefficients[i] *= multiplier;
}

}

namespace discamb
{

using namespace wfndata;

SlaterOrbitalWfnData::SlaterOrbitalWfnData(LocalDef_Wfn wfns[], string *configuration_data, const int max_atomic_number)
{
	int entryIdx;
	HC_WfnBankEntry entry;
    // LocalDef_Wfn wfns
    /*
    struct LocalDef_Wfn {
    std::string label;
    int atomicNumber;
    int charge;
    int nOrbs;
    LocalDef_Orbital orbs[20];
    };
    */
    map<string, vector<vector<int> > > configurations;
    vector<pair<int, int> > val_orb;
    pair<int, int> orbitalIndex;
    setConfigurationData(configurations, configuration_data);
    int stoIdx, nSTO;
    double angstrom = 1.0 / 0.52917721092;
    double normalizationFactor;

	for (entryIdx = 0; entryIdx < max_atomic_number; entryIdx++)
	{
        HC_WfnBankEntry entry;
        entry.atomic_number = wfns[entryIdx].atomicNumber;
        entry.charge = wfns[entryIdx].charge;
        valenceOrbitals(wfns[entryIdx], configurations[wfns[entryIdx].label], val_orb);
        
        entry.orbitals.resize(wfns[entryIdx].nOrbs);
        entry.orbital_occupancy.resize(wfns[entryIdx].nOrbs);
        entry.wfnAtomTypeName = wfns[entryIdx].label;


        for (int i = 0; i < wfns[entryIdx].nOrbs; i++)
        {
            orbitalIndex.first = wfns[entryIdx].orbs[i].principal_number;
            orbitalIndex.second = wfns[entryIdx].orbs[i].azimutal_number;
            if (find(val_orb.begin(), val_orb.end(), orbitalIndex) != val_orb.end())
                entry.valence_orbitals_indices.push_back(i);
            else
                entry.core_orbitals_indices.push_back(i);

            entry.orbital_occupancy[i] = configurations[wfns[entryIdx].label][orbitalIndex.first - 1][orbitalIndex.second];

            nSTO = wfns[entryIdx].orbs[i].nTerms;

            entry.orbitals[i].coefficients.resize(nSTO);
            entry.orbitals[i].exponents.resize(nSTO);
            entry.orbitals[i].power_r.resize(nSTO);

            for (stoIdx = 0; stoIdx < nSTO; stoIdx++)
            {
                
                entry.orbitals[i].exponents[stoIdx] = wfns[entryIdx].orbs[i].term[stoIdx].exponent * angstrom;
                entry.orbitals[i].power_r[stoIdx] = wfns[entryIdx].orbs[i].term[stoIdx].powerR;
                normalizationFactor = sto_atomic_wfn::stoNormalizationFactor(entry.orbitals[i].power_r[stoIdx],
                                                                             entry.orbitals[i].exponents[stoIdx]);
                
                entry.orbitals[i].coefficients[stoIdx] = normalizationFactor * 
                                                         wfns[entryIdx].orbs[i].term[stoIdx].coefficient;
            }

            entry.orbitals[i].azimuthal_number = wfns[entryIdx].orbs[i].azimutal_number;
            entry.orbitals[i].principal_number = wfns[entryIdx].orbs[i].principal_number;
            renormalizeOrbital(entry.orbitals[i]);
        }
        
        mEntries[entry.wfnAtomTypeName] = entry;
	}
}

SlaterOrbitalWfnData::~SlaterOrbitalWfnData()
{

}

void SlaterOrbitalWfnData::getEntries(
    std::vector<discamb::HC_WfnBankEntry> &entries)
{
    std::map<std::string, discamb::HC_WfnBankEntry>::const_iterator it = mEntries.begin();
    entries.clear();
    for (; it != mEntries.end(); it++)
        entries.push_back(it->second);
}

const HC_WfnBankEntry &SlaterOrbitalWfnData::getEntry(
    const std::string &atomType, 
    bool &hasType) 
const
{
    hasType = (mEntries.find(atomType) != mEntries.end());

    if (hasType)
        return mEntries.find(atomType)->second;
    return mEmptyEntry;
}

const HC_WfnBankEntry &SlaterOrbitalWfnData::getEntry(
    const std::string &atomType)
const
{
    std::map<std::string, discamb::HC_WfnBankEntry>::const_iterator it;
    it = mEntries.find(atomType);
    if (it == mEntries.end())
        on_error::throwException("Invalid atomType label when attempting to access Clementi-Roetti wavefunction data",
            __FILE__, __LINE__);
    return mEntries.find(atomType)->second;
}

}