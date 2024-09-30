#include "discamb/Scattering/BankHcFunctions.h"

#include "discamb/BasicChemistry/PeriodicTable.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/HC_Model/DeformationValenceParameters.h"
#include "discamb/MathUtilities/RealSphericalHarmonics.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/Scattering/MATTS_Default.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"


#include <sstream>
#include <set>
#include <map>
#include <fstream>

using namespace std;

namespace discamb {

    BankHcFunctions::BankHcFunctions() 
    {
        setDefault();
    }
    
    BankHcFunctions::BankHcFunctions(
        const std::string& hcBankFileName,
        bool addSpherical)
    {
        set(hcBankFileName, addSpherical);
    }

    void BankHcFunctions::set(
        const std::string& hcBankFileName,
        bool addSpherical)
    {
        ifstream in(hcBankFileName);

        if (!in.good())
            on_error::throwException("cannot read file '" + hcBankFileName + "'", __FILE__, __LINE__);
        set(in, addSpherical);
        in.close();
    }

    void BankHcFunctions::setHC_ModelParameters(
        const std::vector<AtomTypeHC_Parameters>& bankMultipoleParameters,
        int atomType,
        int atomicNumber,
        HC_ModelParameters& parameters)
    {
        parameters.atom_to_type_map.assign(1, 0);
        parameters.atom_to_wfn_map.assign(1, 0);

        // set wfn types

        parameters.wfn_parameters.resize(1);
        discamb::ClementiRoettiData clementiRoettiData;
        discamb::DeformationValenceParameters def_val(discamb::DeformationValenceParameters::ParametersType::UBDB);

        string chemicalElementSymbol = periodic_table::symbol(atomicNumber);

        parameters.wfn_parameters[0] = clementiRoettiData.getEntry(chemicalElementSymbol);
        def_val.getParameters(chemicalElementSymbol, parameters.wfn_parameters[0].deformation_valence_exponent,
            parameters.wfn_parameters[0].deformation_valence_power);

        // set bank types

        parameters.type_parameters.resize(1);

        parameters.type_parameters[0].kappa_deformation_valence = bankMultipoleParameters[atomType].kappa_prime;
        parameters.type_parameters[0].kappa_spherical_valence = bankMultipoleParameters[atomType].kappa;
        parameters.type_parameters[0].p_val = bankMultipoleParameters[atomType].p_val;
        int maxL = -1;
        for (auto p : bankMultipoleParameters[atomType].p_lm_indices)
            maxL = std::max(maxL, int(p.first));

        // set all Plm to zero
        vector<vector<double> >& plm = parameters.type_parameters[0].p_lm;
        if (maxL > -1)
        {
            int nL = int(maxL) + 1;
            plm.resize(nL);
            for (int i = 0; i < nL; i++)
                plm[i].assign(2 * i + 1, 0);
        }
        // set values for non-zero Plm parameters
        for (int i = 0; i < bankMultipoleParameters[atomType].p_lm_indices.size(); i++)
        {
            int l = int(bankMultipoleParameters[atomType].p_lm_indices[i].first);
            int m = bankMultipoleParameters[atomType].p_lm_indices[i].second;
            int idx = l + m;
            plm[l][idx] = bankMultipoleParameters[atomType].p_lms[i];
        }

    }


    void BankHcFunctions::set(
        std::istream& bankStream, 
        bool addSpherical)
    {
        mTypeParameters.clear();
        mTypeIdx.clear();
        mTypeLabel.clear();

        MATTS_BankReader bankReader;
        vector<AtomType> types;
        vector<AtomTypeHC_Parameters> typeParams;
        vector<AtomTypeHC_Parameters> bankMultipoleParameters;
        BankSettings bankSettings;
        bankReader.read(bankStream, types, bankMultipoleParameters, bankSettings, addSpherical);
        int typeIdx, nTypes = types.size();

        mTypeParameters.resize(nTypes);
        
        for (typeIdx = 0; typeIdx < nTypes; typeIdx++)
        {
            int z = *types[typeIdx].atoms[0].atomic_number_range.begin();
            setHC_ModelParameters(bankMultipoleParameters, typeIdx, z, mTypeParameters[typeIdx]);
            mTypeIdx[types[typeIdx].id] = typeIdx;
            mTypeLabel.push_back(types[typeIdx].id);
        }

    }


    void BankHcFunctions::setDefault(
        bool addSpherical)
    {
        string bankString;
        stringstream bankStream;
        default_ubdb_bank_string(bankString);
        bankStream << bankString;
        set(bankStream, addSpherical);
    }

    std::complex<double> BankHcFunctions::calculateComplex(
        HC_ModelParameters& hcParams,
        double h)
        const
    {
        int atom_to_type_map = hcParams.atom_to_type_map[0];
        int type_2_wfn_map = hcParams.atom_to_wfn_map[0];
        //int nH = h_vec.size();
        Crystal crystal;
        bool h000 = (h == 0);
        if (!h000)
            crystal.unitCell.set(1.0 / h, 1.0 / h, 1.0 / h, 90, 90, 90);
        //ReciprocalLatticeUnitCell reciprocalLatticeUnitCell(crystal.unitCell);
        //Vector3d hCart;
        //if (h000)
        //    reciprocalLatticeUnitCell.fractionalToCartesian(Vector3d(0, 0, 0), hCart);
        //else
        //    reciprocalLatticeUnitCell.fractionalToCartesian(Vector3d(1, 0, 0), hCart);
        //cout << "hCart = " << hCart[0] << " " << hCart[1] << " " << hCart[2] << "\n";
        crystal.atoms.resize(1);
        AtomInCrystal& atom = crystal.atoms[0];
        atom.adp.resize(1, 0.0);
        atom.coordinates.set(0, 0, 0);
        atom.type = periodic_table::symbol(hcParams.wfn_parameters[0].atomic_number);
        atom.multiplicity = 1;
        atom.occupancy = 1.0;

        vector < shared_ptr <LocalCoordinateSystemInCrystal> > lcs;
        lcs.push_back(shared_ptr<LocalCoordinateSystemInCrystal>(new LocalCoordinateSystemCalculator()));
        HcFormFactorCalculationsManager manager(crystal, hcParams, lcs);

        complex<double> f;
        if (h000)
            f = manager.calculateFrac(0, Vector3i(0, 0, 0));
        else
            f = manager.calculateFrac(0, Vector3i(0, 0, 1));
        return f;
    }

    double BankHcFunctions::calculate(
        HC_ModelParameters& hcParams,
        double h)
        const 
    {
        complex<double> f = calculateComplex(hcParams, h);
        return f.real();
    }


    double BankHcFunctions::f_core(
        int typeIdx,
        double h)
        const
    {
        auto pVal = mTypeParameters[typeIdx].type_parameters[0].p_val;
        vector<vector<double> > plm;
        plm.swap(mTypeParameters[typeIdx].type_parameters[0].p_lm);
        

        mTypeParameters[typeIdx].type_parameters[0].p_val = 0;
        mTypeParameters[typeIdx].type_parameters[0].p_lm.clear();

        double result = calculate(mTypeParameters[typeIdx], h);

        mTypeParameters[typeIdx].type_parameters[0].p_val = pVal;
        mTypeParameters[typeIdx].type_parameters[0].p_lm.swap(plm);
        return result;
    }


    double BankHcFunctions::f_val(
        int typeIdx,
        double h)
        const
    {
        double fCore = f_core(typeIdx, h);
        double pVal = mTypeParameters[typeIdx].type_parameters[0].p_val;
        vector<vector<double> > plm;

        plm.swap(mTypeParameters[typeIdx].type_parameters[0].p_lm);
        mTypeParameters[typeIdx].type_parameters[0].p_val = 1;

        double fVal = calculate(mTypeParameters[typeIdx], h);

        fVal -= fCore;

        mTypeParameters[typeIdx].type_parameters[0].p_val = pVal;
        plm.swap(mTypeParameters[typeIdx].type_parameters[0].p_lm);

        return fVal;
    }

    double BankHcFunctions::p_val(
        int typeIdx)
        const
    {
        return mTypeParameters[typeIdx].type_parameters[0].p_val;
    }

    std::string BankHcFunctions::getTypeLabel(
        int idx)
        const
    {
        return mTypeLabel[idx];
    }

    double BankHcFunctions::fourier_bessel_transform(
        int typeIdx,
        int l,
        double h)
        const
    {
        if (mTypeParameters[typeIdx].type_parameters[0].p_lm.empty())
            return 0.0;

        double fCore = f_core(typeIdx, h);
        double pVal = mTypeParameters[typeIdx].type_parameters[0].p_val;
        vector<vector<double> > plm;

        int maxL = mTypeParameters[typeIdx].type_parameters[0].p_lm.size() - 1;

        if(l>maxL)
            return 0.0;

        plm.resize(maxL + 1);
        for (int i = 0; i <= maxL; i++)
            plm[i].assign(2 * i + 1, 0.0);

        plm[l][l] = 1.0;




        plm.swap(mTypeParameters[typeIdx].type_parameters[0].p_lm);
        mTypeParameters[typeIdx].type_parameters[0].p_val = 0;

        complex<double> f = calculateComplex(mTypeParameters[typeIdx], h);

        f -= complex<double>(fCore, 0.0);

        double d_l0 = real_spherical_harmonics::densityNormalized(Vector3d(0, 0, 1), int(l), 0);
        double result;
        if (l % 2 == 0)
            result = f.real();
        else
            result = f.imag();

        vector<double> m = { 1,1,-1,-1,1 };

        result /= 4.0 * M_PI * d_l0 * m[l];

        

        mTypeParameters[typeIdx].type_parameters[0].p_val = pVal;
        plm.swap(mTypeParameters[typeIdx].type_parameters[0].p_lm);

        return result;
    }
    
    int BankHcFunctions::getTypeIdx(
        const std::string& typeLabel)
        const
    {
        if (hasType(typeLabel))
            return mTypeIdx.find(typeLabel)->second;
        return 0;
    }

    bool BankHcFunctions::hasType(
        const std::string& typeLabel)
        const
    {
        if (mTypeIdx.find(typeLabel) != mTypeIdx.end())
            return true;
        return false;
    }

}

