#define _CRTDBG_MAP_ALLOC

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/agreement_factors.h"
#include "discamb/Scattering/ConstFormFactorCalculationsManager.h"
#include "discamb/Scattering/SfCalculator.h"
#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include <fstream>


using namespace std;
using namespace discamb;




double target(
    const vector<double>& I_obs,
    const vector<double>& I_calc,
    const vector<complex<double> >& F_calc)
{
    double sum = 0;
    for (auto const& f : F_calc)
        sum += f.real() * f.real();
    return sum;
    //double numerator = 0.0;
    //double denominator = 0.0;
    //for (double i_obs : I_obs)
    //    denominator += i_obs * i_obs;

    //for (int k = 0; k < I_obs.size(); k++)
    //{
    //    double diff = I_obs[k] - I_calc[k];
    //    numerator += diff * diff;
    //}

    ////return sqrt(numerator / denominator);
    //return numerator;
}




//void test_derivatives(
//    const string &structure_file, 
//    const string& hkl_file)
//{
//    double step = 1e-5;
//
//    Crystal crystal;
//    structure_io::read_structure(structure_file, crystal);
//    vector<Vector3i> hkls;
//    vector<double> intensities, sigmas;
//    vector<int> batch;
//    hkl_io::readShelxHkl(hkl_file, hkls, intensities, sigmas, batch, false);
//
//    nlohmann::json settings = nlohmann::json::parse(R"(
//  {
//    "model": "taam",
//    "bank path":"MATTS2021databank.txt",
//    "#algorithm": "macromol",
//    "frozen lcs": true
//  }
//)");
//
//
//
//
//    auto sf_calculator = SfCalculator::create_shared_ptr(crystal, settings);
//    //sf_calculator->calculateStructureFactorsAndDerivatives()
//
//    
//    int nAtoms = crystal.atoms.size();
//    vector<TargetFunctionAtomicParamDerivatives> numerical_derivatives(nAtoms), analytical_derivatives, analytical_derivatives_2(nAtoms);
//    
//    vector<complex<double> > f_calc, f_model;
//    int hklIdx, nHkl = hkls.size();
//    vector<double> I_obs, I_calc(hkls.size());
//    make_f_model(crystal, hkls, settings, f_model);
//    for (auto& f : f_model)
//        I_obs.push_back(norm(f));
//    sf_calculator->calculateStructureFactors(crystal.atoms, hkls, f_calc);
//
//    for (int coordIdx = 0; coordIdx < 3; coordIdx++)
//        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
//        {
//            auto atoms = crystal.atoms;
//            //cart[coordIdx] += step;
//            atoms[atomIdx].coordinates[coordIdx] += step;
//            sf_calculator->calculateStructureFactors(atoms, hkls, f_calc);
//            for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//                I_calc[hklIdx] = norm(f_calc[hklIdx]);
//            double target_plus = target(I_obs, I_calc, f_calc);
//
//            atoms[atomIdx].coordinates[coordIdx] -= 2*step;
//            sf_calculator->calculateStructureFactors(atoms, hkls, f_calc);
//            for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//                I_calc[hklIdx] = norm(f_calc[hklIdx]);
//            double target_minus = target(I_obs, I_calc, f_calc);
//
//            numerical_derivatives[atomIdx].atomic_position_derivatives[coordIdx] = (target_plus - target_minus) / (2 * step);
//        }
//
//        //for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
//        //{
//        //    auto atoms = crystal.atoms;
//        //    Vector3d cart, frac;
//        //    crystal.unitCell.fractionalToCartesian(atoms[atomIdx].coordinates, cart);
//        //    cart[coordIdx] += step;
//        //    crystal.unitCell.cartesianToFractional(cart, atoms[atomIdx].coordinates);
//        //    sf_calculator->calculateStructureFactors(atoms, hkls, f_calc);
//        //    for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//        //        I_calc[hklIdx] = norm(f_calc[hklIdx]);
//        //    double target_plus = target(I_obs, I_calc);
//
//        //    cart[coordIdx] -= 2*step;
//        //    crystal.unitCell.cartesianToFractional(cart, atoms[atomIdx].coordinates);
//        //    sf_calculator->calculateStructureFactors(atoms, hkls, f_calc);
//        //    for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//        //        I_calc[hklIdx] = norm(f_calc[hklIdx]);
//        //    double target_minus = target(I_obs, I_calc);
//
//        //    numerical_derivatives[atomIdx].atomic_position_derivatives[coordIdx] = (target_plus - target_minus) / (2 * step);
//        //}
//
//    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
//    {
//        auto atoms = crystal.atoms;
//        atoms[atomIdx].occupancy += step;
//        sf_calculator->calculateStructureFactors(atoms, hkls, f_calc);
//        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//            I_calc[hklIdx] = norm(f_calc[hklIdx]);
//        double target_plus = target(I_obs, I_calc, f_calc);
//        atoms[atomIdx].occupancy -= 2 * step;
//        sf_calculator->calculateStructureFactors(atoms, hkls, f_calc);
//        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//            I_calc[hklIdx] = norm(f_calc[hklIdx]);
//        double target_minus = target(I_obs, I_calc, f_calc);
//        numerical_derivatives[atomIdx].occupancy_derivatives = (target_plus - target_minus) / (2 * step);
//
//        for (int i = 0; i < atoms[atomIdx].adp.size(); i++)
//            numerical_derivatives[atomIdx].adp_derivatives.push_back(0.0);
//
//        //for (int i = 0; i < atoms[atomIdx].adp.size(); i++)
//        //{
//        //    atoms[atomIdx].adp[i] += stepAdp;
//        //    sfCalculator->update(atoms);
//        //    sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_plus, derivatives_plus, countAtom);
//        //    atoms[atomIdx].adp[i] -= 2 * stepAdp;
//        //    sfCalculator->update(atoms);
//        //    sfCalculator->calculateStructureFactorsAndDerivatives(hkl, sf_minus, derivatives_minus, countAtom);
//        //    derivatives.adpDerivatives[atomIdx][i] = (sf_plus - sf_minus) / (2 * stepAdp);
//        //    atoms[atomIdx].adp[i] -= stepAdp;
//        //}
//
//    }
//    vector<complex<double> > dt_df;
//    dTarget_dF(I_obs, f_calc, dt_df);
//    sf_calculator->calculateStructureFactorsAndDerivatives(crystal.atoms, hkls, f_calc, analytical_derivatives, dt_df);
//
//    /*analytical derivatives 2*/
//    vector<bool> count_atom(nAtoms, true);
//    for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
//    {
//        complex<double> sf;
//        SfDerivativesAtHkl df;
//        sf_calculator->calculateStructureFactorsAndDerivatives(hkls[hklIdx], sf, df, count_atom);
//        //analytical_derivatives_2
//        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
//        {
//            for (int i = 0; i < 3; i++)
//                analytical_derivatives_2[atomIdx].atomic_position_derivatives[i] += dt_df[hklIdx].real() * df.atomicPostionDerivatives[atomIdx][i].real() +
//                dt_df[hklIdx].imag() * df.atomicPostionDerivatives[atomIdx][i].imag();
//        }
//        
//    }
//
//    for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
//    {
//        cout << crystal.atoms[atomIdx].label << "\n";
//        cout << "dT/dxyz    numerical    analytical    analytical2 \n";
//        string xyz_str = "xyz";
//        for (int i = 0; i < 3; i++)
//        {
//            cout << "  " << xyz_str[i] << "  " << numerical_derivatives[atomIdx].atomic_position_derivatives[i]
//                << "  " << analytical_derivatives[atomIdx].atomic_position_derivatives[i]  
//                << "  " << analytical_derivatives_2[atomIdx].atomic_position_derivatives[i] << "\n";
//        }
//        cout<< "  occ " << "  " << numerical_derivatives[atomIdx].occupancy_derivatives
//            << "  " << analytical_derivatives[atomIdx].occupancy_derivatives << "\n";
//    }
//    // derivatives dF/dp
//
//
//}












void save_one_hkl(
    const Vector3i& hkl,
    vector<Vector3i>& hkls,
    vector<double> & intensities,
    vector<double> & weights)
{ 
    auto it = find(hkls.begin(), hkls.end(), hkl);
    if (it == hkls.end())
        on_error::throwException("hkl " + hkl.string() + " not found in the list.", __FILE__, __LINE__);
    int idx = distance(hkls.begin(), it);
    double intensity = intensities[idx];
    double weight = weights[idx];
    
    hkls.clear();
    intensities.clear();
    weights.clear();
    
    hkls.push_back(hkl);
    intensities.push_back(intensity);
    weights.push_back(weight);
    
}

void prepare_synethetic_data(
const string &structureFile,
const string &hklFile )
{
    Crystal crystal, crystal_distorted;
    structure_io::read_structure(structureFile, crystal);
    vector<Vector3i> hkls;
    vector<double> int_obs, int_calc, sigmas;
    vector<int> batch;
    hkl_io::readShelxHkl(hklFile, hkls, int_obs, sigmas, batch, false);
    
    double scale = 1.0;

    // skip hkl 000
    if (hkls.back() == Vector3i())
    {
        hkls.pop_back();
        int_obs.pop_back();
        sigmas.pop_back();
    }

    shared_ptr<SfCalculator> sf_calculator = SfCalculator::create_shared_ptr(crystal, string("aspher.json"));
    vector<complex<double> > f_calc;

    sf_calculator->calculateStructureFactors(crystal.atoms, hkls, f_calc);


    vector<double> weights;

    for (auto& f : f_calc)
        int_calc.push_back(norm(f));

    scattering_utilities::scaleFactorFcalc(int_obs, int_calc, sigmas, scale);

    for(auto &sigma: sigmas)
        sigma /= scale;

    hkl_io::writeShelxHkl("synthetic.hkl", hkls, int_calc, sigmas, batch, false);

    crystal_structure_utilities::StructureDistortionParameters distortionParams;
    distortionParams.distort_adps = true;
    distortionParams.distort_occupancies = false;
    distortionParams.distort_positions = true;
    distortionParams.maxAdpChange = 0.04;
    crystal_structure_utilities::distort_structure(crystal, crystal_distorted, distortionParams);

    structure_io::write_structure("distorted.cif", crystal_distorted);

}

void test_xyz_derivatives(
    const Crystal &crystal,
    const vector<Vector3i> &hkls,
    const vector<double> &weights)
{
    shared_ptr<SfCalculator> sf_calculator = SfCalculator::create_shared_ptr(crystal, string("aspher.json"));
    vector<complex<double> > f_calc_oryginal, f_calc;

    sf_calculator->calculateStructureFactors(crystal.atoms, hkls, f_calc_oryginal);

    vector<double> i_obs, i_calc;

    for (auto& f : f_calc_oryginal)
        i_obs.push_back(norm(f));

    Crystal crystal_distorted;
    crystal_structure_utilities::distort_structure(crystal, crystal_distorted);

    sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, f_calc);

    for (auto& f : f_calc)
        i_calc.push_back(norm(f));
    cout << "R2 " << agreement_factors::value(i_calc, i_obs, weights, agreement_factors::R2);

    double numerator = 0.0;
    double denominator = 0.0;

    for (int i = 0; i < i_obs.size(); i++)
    {
        double diff = i_obs[i] - i_calc[i];
        numerator += weights[i] * diff * diff;
        denominator += weights[i] * i_obs[i] * i_obs[i];
    }

    //##############################################
    //       xyz derivative 1-st atom
    //##############################################

    // numerical derivative
    double step = 0.00001;
    auto& atom = crystal_distorted.atoms[0];
    double tartget_plus, target_minus;
    vector<complex<double> > f_calc_plus, f_calc_minus;
    atom.coordinates[0] += step;
    sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, f_calc_plus);
    atom.coordinates[0] -= step;
    i_calc.clear();
    for (auto& f : f_calc_plus)
        i_calc.push_back(norm(f));
    tartget_plus = agreement_factors::value(i_calc, i_obs, weights, agreement_factors::R2);

    atom.coordinates[0] -= step;
    sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, f_calc_minus);
    atom.coordinates[0] += step;
    i_calc.clear();
    for (auto& f : f_calc_minus)
        i_calc.push_back(norm(f));
    target_minus = agreement_factors::value(i_calc, i_obs, weights, agreement_factors::R2);

    double numerical_derivative = (tartget_plus - target_minus) / (2 * step);
    cout << "\nNumerical derivative dR2/dx1: " << numerical_derivative << endl;


    // analytical derivative


    sf_calculator->update(crystal_distorted.atoms);
    vector<bool> count_atom(crystal_distorted.atoms.size(), true);
    SfDerivativesAtHkl derivatives;
    SfDerivativesAtHkl cumulativeDerivatives;
    int atomIdx, nAtoms = crystal_distorted.atoms.size();

    // initialize cumulative derivatives
    cumulativeDerivatives.adpDerivatives.resize(nAtoms);
    cumulativeDerivatives.atomicPostionDerivatives.resize(nAtoms);
    complex<double> complex_zero(0.0, 0.0);
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        cumulativeDerivatives.adpDerivatives[atomIdx].resize(crystal.atoms[atomIdx].adp.size(), complex_zero);
        cumulativeDerivatives.atomicPostionDerivatives[atomIdx].set(complex_zero, complex_zero, complex_zero);
        cumulativeDerivatives.occupancyDerivatives.push_back(complex_zero);
    }

    double d = 0.0;

    for (int i = 0; i < hkls.size(); i++)
    {
        complex<double> scatteringFactor;
        sf_calculator->calculateStructureFactorsAndDerivatives(hkls[i], scatteringFactor, derivatives, count_atom);

        complex<double> dx = derivatives.atomicPostionDerivatives[0].x;
        double diff = (i_obs[i] - norm(scatteringFactor));
        d += weights[i] * diff * (scatteringFactor.real() * dx.real() + scatteringFactor.imag() * dx.imag());
        //d += -2.0 * weights[i] * (i_obs[i] - norm(scatteringFactor)) * dx.imag();
        //denominator += weights[i] * i_obs[i];
    }
    d *= -200.0 / sqrt(numerator * denominator);
    cout << "\nAnalytical derivative dR2/dx1: " << d << endl;


    // analytical derivative with chain rule


    double d_chain = 0.0;
    vector<complex<double> > dt_df(hkls.size(), 0.0), sf;
    sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, sf, vector<bool>(crystal_distorted.atoms.size(), true));
    for (int i = 0; i < hkls.size(); i++)
    {
        double diff = i_obs[i] - norm(sf[i]);
        dt_df[i] = -200.0 * weights[i] * diff * sf[i] / sqrt(numerator * denominator);
    }
    vector<TargetFunctionAtomicParamDerivatives> dt_dp;
    sf_calculator->calculateStructureFactorsAndDerivatives(crystal_distorted.atoms, hkls, sf, dt_dp, dt_df, vector<bool>(crystal_distorted.atoms.size(), true));

    cout << "\nAnalytical derivative v.2 dR2/dx1: " << dt_dp[0].atomic_position_derivatives.x << endl;


    // analytical derivative with chain rule & no form factor recalculation


    d_chain = 0.0;
    vector < vector<complex<double> > > form_factors;
    vector<Vector3i> symmetryEquivalentHkls;
    scattering_utilities::generate_symmetry_equivalent_hkls(
        crystal.spaceGroup,
        hkls,
        symmetryEquivalentHkls);
    sf_calculator->calculateFormFactors(symmetryEquivalentHkls, form_factors, vector<bool>(crystal_distorted.atoms.size(), true));
    map<Vector3i, vector<complex<double> > > form_factors_map;
    for (int i = 0; i < symmetryEquivalentHkls.size(); i++)
        form_factors_map[symmetryEquivalentHkls[i]] = form_factors[i];
    shared_ptr<AtomicFormFactorCalculationsManager> ff_manager = make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, form_factors_map);
    AnyScattererStructureFactorCalculator asf_calculator(crystal_distorted);
    asf_calculator.setAtomicFormfactorManager(ff_manager);

    dt_df.assign(hkls.size(), complex<double>());
    sf.clear();
    asf_calculator.calculateStructureFactors(hkls, sf, vector<bool>(crystal_distorted.atoms.size(), true));
    for (int i = 0; i < hkls.size(); i++)
    {
        complex<double> scatteringFactor;
        asf_calculator.calculateStructureFactorsAndDerivatives(hkls[i], scatteringFactor, derivatives, count_atom);

        double diff = i_obs[i] - norm(sf[i]);
        dt_df[i] = -200.0 * weights[i] * diff * sf[i] / sqrt(numerator * denominator);
    }
    vector<TargetFunctionAtomicParamDerivatives> dt_dp_nff_rec;
    asf_calculator.calculateStructureFactorsAndDerivatives(hkls, sf, dt_dp_nff_rec, dt_df, vector<bool>(crystal_distorted.atoms.size(), true));

    cout << "\nAnalytical derivative v.3 dR2/dx1: " << dt_dp_nff_rec[0].atomic_position_derivatives.x << endl;


    // analytical derivative without chain rule & no form factor recalculation

    vector < vector<complex<double> > > form_factors2;
    sf_calculator->update(crystal_distorted.atoms);
    sf_calculator->calculateFormFactors(symmetryEquivalentHkls, form_factors2, vector<bool>(crystal_distorted.atoms.size(), true));
    form_factors_map.clear();

    for (int i = 0; i < symmetryEquivalentHkls.size(); i++)
        form_factors_map[symmetryEquivalentHkls[i]] = form_factors2[i];

    shared_ptr<AtomicFormFactorCalculationsManager> ff_manager2 = make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, form_factors_map);
    AnyScattererStructureFactorCalculator asf_calculator2(crystal_distorted);
    asf_calculator2.setAtomicFormfactorManager(ff_manager2);

    d = 0.0;
    for (int i = 0; i < hkls.size(); i++)
    {
        complex<double> scatteringFactor;
        asf_calculator2.calculateStructureFactorsAndDerivatives(hkls[i], scatteringFactor, derivatives, vector<bool>(crystal_distorted.atoms.size(), true));

        complex<double> dx = derivatives.atomicPostionDerivatives[0].x;
        double diff = (i_obs[i] - norm(scatteringFactor));
        d += weights[i] * diff * (scatteringFactor.real() * dx.real() + scatteringFactor.imag() * dx.imag());
    }
    d *= -200.0 / sqrt(numerator * denominator);
    cout << "\nAnalytical derivative v.4 dR2/dx1: " << d << endl;

}


void test_adps_derivatives(
    const Crystal& crystal,
    const vector<Vector3i>& hkls,
    const vector<double>& weights)
{
    shared_ptr<SfCalculator> sf_calculator = SfCalculator::create_shared_ptr(crystal, string("aspher.json"));
    vector<complex<double> > f_calc_oryginal, f_calc;

    sf_calculator->calculateStructureFactors(crystal.atoms, hkls, f_calc_oryginal);

    vector<double> i_obs, i_calc;

    for (auto& f : f_calc_oryginal)
        i_obs.push_back(norm(f));

    Crystal crystal_distorted;
    crystal_structure_utilities::distort_structure(crystal, crystal_distorted);

    sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, f_calc);

    for (auto& f : f_calc)
        i_calc.push_back(norm(f));
    cout << "R2 " << agreement_factors::value(i_calc, i_obs, weights, agreement_factors::R2);

    double numerator = 0.0;
    double denominator = 0.0;

    for (int i = 0; i < i_obs.size(); i++)
    {
        double diff = i_obs[i] - i_calc[i];
        numerator += weights[i] * diff * diff;
        denominator += weights[i] * i_obs[i] * i_obs[i];
    }

    //##############################################
    //       adps derivative 1-st atom
    //##############################################
    vector<double> numerical_derivatives;

    for (int adpComponentIdx = 0; adpComponentIdx < 6; adpComponentIdx++)
    {

        // numerical derivative
        double step = 0.00001;
        auto& atom = crystal_distorted.atoms[0];

        if (atom.adp.size() < 6)
            on_error::throwException("Atom has no anisotropic ADPs.", __FILE__, __LINE__);

        double tartget_plus, target_minus;
        vector<complex<double> > f_calc_plus, f_calc_minus;
        atom.adp[adpComponentIdx] += step;
        sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, f_calc_plus);
        atom.adp[adpComponentIdx] -= step;
        i_calc.clear();
        for (auto& f : f_calc_plus)
            i_calc.push_back(norm(f));
        tartget_plus = agreement_factors::value(i_calc, i_obs, weights, agreement_factors::R2);

        atom.adp[adpComponentIdx] -= step;
        sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, f_calc_minus);
        atom.adp[adpComponentIdx] += step;
        i_calc.clear();
        for (auto& f : f_calc_minus)
            i_calc.push_back(norm(f));
        target_minus = agreement_factors::value(i_calc, i_obs, weights, agreement_factors::R2);

        double numerical_derivative = (tartget_plus - target_minus) / (2 * step);
        cout << "\nNumerical derivative dR2/dx1: " << numerical_derivative << endl;
    }

    // analytical derivative


    sf_calculator->update(crystal_distorted.atoms);
    vector<bool> count_atom(crystal_distorted.atoms.size(), true);
    SfDerivativesAtHkl derivatives;
    SfDerivativesAtHkl cumulativeDerivatives;
    int atomIdx, nAtoms = crystal_distorted.atoms.size();

    // initialize cumulative derivatives
    cumulativeDerivatives.adpDerivatives.resize(nAtoms);
    cumulativeDerivatives.atomicPostionDerivatives.resize(nAtoms);
    complex<double> complex_zero(0.0, 0.0);
    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        cumulativeDerivatives.adpDerivatives[atomIdx].resize(crystal.atoms[atomIdx].adp.size(), complex_zero);
        cumulativeDerivatives.atomicPostionDerivatives[atomIdx].set(complex_zero, complex_zero, complex_zero);
        cumulativeDerivatives.occupancyDerivatives.push_back(complex_zero);
    }
    
    double d = 0.0;
    vector<double> analyticalDerivatives(6, 0.0);
    
    for (int i = 0; i < hkls.size(); i++)
    {
        complex<double> scatteringFactor;
        sf_calculator->calculateStructureFactorsAndDerivatives(hkls[i], scatteringFactor, derivatives, count_atom);

        complex<double> dp = derivatives.adpDerivatives[0][0];
        double diff = (i_obs[i] - norm(scatteringFactor));
        d += weights[i] * diff * (scatteringFactor.real() * dp.real() + scatteringFactor.imag() * dp.imag());
        for(int adpCompIdx = 0; adpCompIdx < 6; adpCompIdx++)
        {
            complex<double> dp = derivatives.adpDerivatives[0][adpCompIdx];
            analyticalDerivatives[adpCompIdx] += weights[i] * diff * (scatteringFactor.real() * dp.real() + scatteringFactor.imag() * dp.imag());
        }
        //d += -2.0 * weights[i] * (i_obs[i] - norm(scatteringFactor)) * dx.imag();
        //denominator += weights[i] * i_obs[i];
    }
    d *= -200.0 / sqrt(numerator * denominator);
    for (int adpCompIdx = 0; adpCompIdx < 6; adpCompIdx++)
    {
        analyticalDerivatives[adpCompIdx] *= -200.0 / sqrt(numerator * denominator);
        cout << "\nAnalytical derivative dR2/dAdp" << adpCompIdx + 1 << ": " << analyticalDerivatives[adpCompIdx] << endl;
    }
    //cout << "\nAnalytical derivative dR2/dx1: " << d << endl;
    
    return;
    // analytical derivative with chain rule


    double d_chain = 0.0;
    vector<complex<double> > dt_df(hkls.size(), 0.0), sf;
    sf_calculator->calculateStructureFactors(crystal_distorted.atoms, hkls, sf, vector<bool>(crystal_distorted.atoms.size(), true));
    for (int i = 0; i < hkls.size(); i++)
    {
        double diff = i_obs[i] - norm(sf[i]);
        dt_df[i] = -200.0 * weights[i] * diff * sf[i] / sqrt(numerator * denominator);
    }
    vector<TargetFunctionAtomicParamDerivatives> dt_dp;
    sf_calculator->calculateStructureFactorsAndDerivatives(crystal_distorted.atoms, hkls, sf, dt_dp, dt_df, vector<bool>(crystal_distorted.atoms.size(), true));

    cout << "\nAnalytical derivative v.2 dR2/dx1: " << dt_dp[0].atomic_position_derivatives.x << endl;


    // analytical derivative with chain rule & no form factor recalculation

    
    d_chain = 0.0;
    vector < vector<complex<double> > > form_factors;
    vector<Vector3i> symmetryEquivalentHkls;
    scattering_utilities::generate_symmetry_equivalent_hkls(
        crystal.spaceGroup,
        hkls,
        symmetryEquivalentHkls);
    sf_calculator->calculateFormFactors(symmetryEquivalentHkls, form_factors, vector<bool>(crystal_distorted.atoms.size(), true));
    map<Vector3i, vector<complex<double> > > form_factors_map;
    for (int i = 0; i < symmetryEquivalentHkls.size(); i++)
        form_factors_map[symmetryEquivalentHkls[i]] = form_factors[i];
    shared_ptr<AtomicFormFactorCalculationsManager> ff_manager = make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, form_factors_map);
    AnyScattererStructureFactorCalculator asf_calculator(crystal_distorted);
    asf_calculator.setAtomicFormfactorManager(ff_manager);

    dt_df.assign(hkls.size(), complex<double>());
    sf.clear();
    asf_calculator.calculateStructureFactors(hkls, sf, vector<bool>(crystal_distorted.atoms.size(), true));
    for (int i = 0; i < hkls.size(); i++)
    {
        complex<double> scatteringFactor;
        asf_calculator.calculateStructureFactorsAndDerivatives(hkls[i], scatteringFactor, derivatives, count_atom);

        double diff = i_obs[i] - norm(sf[i]);
        dt_df[i] = -200.0 * weights[i] * diff * sf[i] / sqrt(numerator * denominator);
    }
    vector<TargetFunctionAtomicParamDerivatives> dt_dp_nff_rec;
    asf_calculator.calculateStructureFactorsAndDerivatives(hkls, sf, dt_dp_nff_rec, dt_df, vector<bool>(crystal_distorted.atoms.size(), true));

    cout << "\nAnalytical derivative v.3 dR2/dx1: " << dt_dp_nff_rec[0].atomic_position_derivatives.x << endl;


    // analytical derivative without chain rule & no form factor recalculation

    vector < vector<complex<double> > > form_factors2;
    sf_calculator->update(crystal_distorted.atoms);
    sf_calculator->calculateFormFactors(symmetryEquivalentHkls, form_factors2, vector<bool>(crystal_distorted.atoms.size(), true));
    form_factors_map.clear();

    for (int i = 0; i < symmetryEquivalentHkls.size(); i++)
        form_factors_map[symmetryEquivalentHkls[i]] = form_factors2[i];

    shared_ptr<AtomicFormFactorCalculationsManager> ff_manager2 = make_shared<ConstFormFactorCalculationsManager>(crystal.unitCell, form_factors_map);
    AnyScattererStructureFactorCalculator asf_calculator2(crystal_distorted);
    asf_calculator2.setAtomicFormfactorManager(ff_manager2);

    d = 0.0;
    for (int i = 0; i < hkls.size(); i++)
    {
        complex<double> scatteringFactor;
        asf_calculator2.calculateStructureFactorsAndDerivatives(hkls[i], scatteringFactor, derivatives, vector<bool>(crystal_distorted.atoms.size(), true));

        complex<double> dx = derivatives.atomicPostionDerivatives[0].x;
        double diff = (i_obs[i] - norm(scatteringFactor));
        d += weights[i] * diff * (scatteringFactor.real() * dx.real() + scatteringFactor.imag() * dx.imag());
    }
    d *= -200.0 / sqrt(numerator * denominator);
    cout << "\nAnalytical derivative v.4 dR2/dx1: " << d << endl;

}


int main(int argc, char* argv[])
{

    try {
        
        //crystal_structure_utilities::generate_hkl(
        //    UnitCell(10.0, 10.0, 10.0, 90.0, 90.0, 90.0),
        //    1.0, hkls);
        

        vector<string> arguments, options;
        parse_cmd::get_args_and_options(argc, argv, arguments, options);

        if (arguments.size() < 2)
        {
            cout << "Usage: pre_lsq <structure file> <hkl file>\n";
            return 0;
        }
        
        if(filesystem::exists("aspher.json") == false)
        {
            cout << "aspher.json is missing.\n";
            return 0;
        }

        prepare_synethetic_data(arguments[0],
            arguments[1]);
        return 0;

        Crystal crystal, crystal_distorted;
        structure_io::read_structure(arguments[0], crystal);
        vector<Vector3i> hkls0, hkls;
        vector<double> intensities, sigmas;
        vector<int> batch;
        hkl_io::readShelxHkl(arguments[1], hkls0, intensities, sigmas, batch, false);
        vector<int> indicesOfOriginals;

        //save_one_hkl(
        //    Vector3i(1, 0, 2),
        //    hkls0,
        //    intensities,
        //    sigmas);

        // skip hkl 000
        if (hkls0.back() == Vector3i())
            hkls0.pop_back();

        crystal_structure_utilities::filter_hkl(crystal.unitCell, 0.8, hkls0, hkls, indicesOfOriginals);
        vector<double> weights;
        
        for (int i = 0; i < indicesOfOriginals.size(); i++)
        {
            double s = sigmas[indicesOfOriginals[i]];
            weights.push_back(1.0 / (s * s));
        }

        test_xyz_derivatives(crystal, hkls, weights);
        test_adps_derivatives(crystal, hkls, weights);
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

