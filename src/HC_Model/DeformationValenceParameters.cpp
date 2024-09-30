#include "discamb/HC_Model/DeformationValenceParameters.h"
#include "discamb/HC_Model/ClementiRaimondiData.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/BasicChemistry/periodic_table.h"

#include <fstream>
#include <iostream>

using namespace std;

namespace{

struct Type
{
    void set(string label);
    string label;
    string symbol;
    int atomicNumber = 1;
    int nElectrons = 1;
    int charge = 0;
    int period = 1;
    discamb::periodic_table::Block block = discamb::periodic_table::Block::NONE;
};

void Type::set(string _label)
{
    string nonSymbol;
    int sign;
    symbol.clear();


    for (int i=0; i<_label.size(); i++)
        if (isalpha(_label[i]))
            symbol += _label[i];
        else
            nonSymbol += _label[i];

    atomicNumber = discamb::periodic_table::atomicNumber(symbol);
    
    block = discamb::periodic_table::block(atomicNumber);
    period = discamb::periodic_table::period(atomicNumber);
    label = _label;

    charge = 0;

    if (!nonSymbol.empty())
    {
        sign = (nonSymbol.find('-') != string::npos ? -1 : 1);

        if (nonSymbol.size() == 1)
            charge = sign;
        else
        {
            for (int i=0; i<nonSymbol.size(); i++)
                if (isdigit(nonSymbol[i]))
                    charge = atoi((string() + nonSymbol[i]).c_str());
            charge *= sign;
        }
    }

    int nEl = static_cast<int>(atomicNumber) - charge;
    nElectrons = static_cast<int>(nEl);
}

void assignPowers(
    vector<int> &powers,
    Type &type)
{
    powers.resize(5);
    // first row
    if (type.period == 1)
        for (int i = 0; i < 5; i++)
            powers[i] = i;

    // second row
    if (type.period == 2)
    {
        powers[0] = 2;
        powers[1] = 2;
        powers[2] = 2;
        powers[3] = 3;
        powers[4] = 4;
    }

    // 3-rd period elements and 4 period transition metals
    if (type.period == 3 || (type.period == 4 && type.block == discamb::periodic_table::Block::D))
        powers.assign(5, 4);

    // 4-rth period main group's elements
    if (type.period == 4 && type.block != discamb::periodic_table::Block::D)
        powers.assign(5, 6);

}


void getConfigurationAndValenceOrbitals(
    const discamb::HC_WfnBankEntry &entry,
    vector<vector<int> > &configuration,
    vector<pair<int, int> > &valence)
{
    int i, nOrbs, n, l;
    configuration.clear();
    valence.clear();

    nOrbs = entry.orbitals.size();
    for (i = 0; i < nOrbs; i++)
    {
        n = entry.orbitals[i].principal_number;
        l = entry.orbitals[i].azimuthal_number;
        if (n > configuration.size())
            configuration.resize(n);
        if (l + 1 > configuration[n - 1].size())
            configuration[n - 1].resize(l + 1);

        configuration[n - 1][l] = entry.orbital_occupancy[i];
        if (find(entry.valence_orbitals_indices.begin(), entry.valence_orbitals_indices.end(), i) != entry.valence_orbitals_indices.end())
            valence.push_back(make_pair(n, l));
        
    }
}



double calcExponent(
    Type &type,
    const vector<vector<double> > exponents,
    const vector<vector<int> > &configuration,
    const vector<pair<int, int> > &valence)
{
    vector<double> exps;
    vector<int> occ;

    if (valence.empty())
    {
        if (type.block == discamb::periodic_table::Block::S)
            return 2 * exponents.back()[0];
        if (type.block == discamb::periodic_table::Block::P)
            return (exponents.back()[0] + 3 * exponents.back()[1]) / 2.0;
        if (type.block == discamb::periodic_table::Block::D)
            return 2 * exponents[exponents.size() - 2].back();

        
        discamb::on_error::throwException(
            "Problem encountered when calculating exponent for deformation valence parameters"
            " in multipolar model - scatterer type definition does not have valid specification"
            " - invalid assignment of block in periodic table, expected S, P or D", __FILE__, __LINE__);

        return 0.0;

    }
    else
    {
        for (int i = 0; i < valence.size(); i++)
        {
            exps.push_back(exponents[valence[i].first - 1][valence[i].second]);
            occ.push_back(configuration[valence[i].first - 1][valence[i].second]);
        }
    }

    double totalOcc = 0;
    double result = 0;

    for (int i = 0; i < occ.size(); i++)
    {
        totalOcc += occ[i];
        result += occ[i] * exps[i];
    }

    return 2 * result / totalOcc;
}

}

namespace discamb {

DeformationValenceParameters::DeformationValenceParameters()
{ 
    setStandardParameterization();
}

DeformationValenceParameters::DeformationValenceParameters(
    ParametersType type)
{
    set(type);
}


DeformationValenceParameters::~DeformationValenceParameters()
{
}

void DeformationValenceParameters::set(
    ParametersType type)
{
    type == ParametersType::STANDARD ? setStandardParameterization() : setUbdbParameterization();
}

void DeformationValenceParameters::setParameters(
    const std::vector<std::string> &types,
    const std::vector<double> &zeta,
    const std::vector<std::vector<int> > &powers_r)
{
    if( types.size()!=zeta.size() || types.size() != powers_r.size())
        on_error::throwException("incompatible sizes of arrays defining deformation valence parameters in multipolar model",
                                 __FILE__,__LINE__);
    int i,n = types.size();
    mParameters.clear();
    for(i=0;i<n;i++)
        mParameters[types[i]] = make_pair(zeta[i],powers_r[i]);

}

void DeformationValenceParameters::setParameter(
    const std::string &type, 
    double zeta, 
    const std::vector<int> &powers_r)
{
    mParameters[type] = make_pair(zeta,powers_r);
}


bool DeformationValenceParameters::getParameters(
    const std::string &type, 
    double &zeta, 
    std::vector<int> &powers_r)
    const
{
    std::map<std::string, std::pair<double, std::vector<int> > >::const_iterator it;
    it = mParameters.find(type);
    if(it == mParameters.end())
        return false;
    
    zeta = it->second.first;
    powers_r = it->second.second;
    return true;
}

void DeformationValenceParameters::getParameters(
    std::vector<std::string> &types, 
    std::vector<double> &zeta,
    std::vector<std::vector<int> > &powers_r)
const
{
    std::map<std::string, std::pair<double, std::vector<int> > >::const_iterator it;
    int i,n = mParameters.size();
    types.resize(n);
    zeta.resize(n);
    powers_r.resize(n);
    i=0;
    for(it = mParameters.begin();it!=mParameters.end();it++)
    {
        types[i] = it->first;
        zeta[i] = it->second.first;
        powers_r[i] = it->second.second;
        i++;
    }
}

bool DeformationValenceParameters::getParameters(
    const std::string &typeName, 
    double &zeta, 
    std::vector<int> &powers_r,
    const std::vector<std::pair<int, int> > &valenceOrbitals,
    const std::vector<std::vector<int> > &_configuration,
    const std::vector<std::vector<double> > &_exponents)
    const
{
    

    Type type;    
    vector<vector<double> > exponents;
    vector<vector<int> > configuration;
    vector<pair<int, int> > valence;

    type.set(typeName);

    if (type.atomicNumber > 36)
        return false;

    assignPowers(powers_r, type);

    if (_configuration.empty())
    {
        ClementiRoettiData cr_data;
        HC_WfnBankEntry entry = cr_data.getEntry(typeName);

        if (entry.atomic_number == 0)
            return false;

        getConfigurationAndValenceOrbitals(entry, configuration, valence);
    }
    else
        configuration = _configuration;

    valence = valenceOrbitals;

    if (_exponents.empty())
        clementi_raimondi::getExponents(type.atomicNumber, exponents);
    else
        exponents = _exponents;

    zeta = calcExponent(type, exponents, configuration, valence);
    
    return true;
}


void DeformationValenceParameters::setUbdbParameterization()
{
    setStandardParameterization();

    std::map<std::string, std::pair<double, std::vector<int> > >::iterator it;
    string symbol;

    for (it = mParameters.begin(); it != mParameters.end(); it++)
    {
        symbol.clear();

        for (int i = 0; i < it->first.size(); i++)
            if (isalpha(it->first[i]))
                symbol += it->first[i];

        if (symbol == string("P"))
        {
            it->second.second = vector<int>(5, 6);
            it->second.second[0] = 4;
        }

        if (symbol == string("S"))
        {
            it->second.second[0] = 1;
            for(int i=1; i<5; i++)
                it->second.second[i] = 2*i;
        }
    }

}


void DeformationValenceParameters::setStandardParameterization()
{
    vector<HC_WfnBankEntry> entries;
    ClementiRoettiData cr_data;
    int entryIdx, nEntries;
    cr_data.getEntries(entries);
    Type type;
    vector<int> powers;
    double exponent;
    //double exponent2;
    vector<vector<double> > exponents;
    vector<vector<int> > configuration;
    vector<pair<int, int> > valence;
    
    nEntries = entries.size();
    

    for ( entryIdx = 0; entryIdx < nEntries; entryIdx++)
    {
        HC_WfnBankEntry &entry = entries[entryIdx];
        type.set(entry.wfnAtomTypeName);
        assignPowers(powers, type);
        getConfigurationAndValenceOrbitals(entry, configuration, valence);
        clementi_raimondi::getExponents(entry.atomic_number, exponents);
        exponent = calcExponent(type, exponents, configuration, valence);
        
        mParameters[entry.wfnAtomTypeName] = make_pair(exponent, powers);
        
    }
}

} // namespace discamb 

