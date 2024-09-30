/*! \file example2.cpp
 
 */

#include "discamb/Scattering/HansenCoppensStructureFactorCalculator.h"
#include "discamb/HC_Model/DeformationValenceParameters.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/IO/xd_io.h"

#include <iostream>

using namespace std;

// Example of structure implementing function onPerHklCalculation
// Such structure/type is necessary for data collection in structure
// factor and derivatives calculations for small molecules

struct DataCollector
{
    std::vector<std::complex<double> > structureFactorsData;
    std::vector<discamb::SfDerivativesAtHkl> derivativesData;

    // this function is called by discamb::HansenCoppensStructureFactorCalculator
    // for each hkl index in order to pass results for that hkl to function
    // which in turn saves it in 
    void onPerHklCalculation(
        size_t hkl_idx,
        std::complex<double> &structureFactor,
        discamb::SfDerivativesAtHkl &derivatives)
    {
        // note that structureFactor and derivativesData should be
        // allocated/resized before calling the function 
        // e.g. structureFactorsData.resize(number_of_hkl);
        structureFactorsData[hkl_idx] = structureFactor;
        derivativesData[hkl_idx] = derivatives;
    }
};



int main(int argc, char *argv[])
{
    //######################   define atom    #################################

    discamb::AtomInCrystal c1, o1, n1, h1, h2;
    c1.label = "C(1)";
    c1.coordinates = { 0.00000,  0.50000,  0.32600 };
    c1.adp = { 0.0088,  0.0088,  0.0035,   0.0002,  0.0000,  0.0000 };
    c1.multiplicity = 2;
    c1.occupancy = 1;
    c1.type = "C";

    // we skip the rest of atoms, in practice all information will be read from XD files

    //#############  define unit cell, transform coordinates  #################

    discamb::UnitCell unitCell(5.565, 5.565, 4.684, 90, 90, 90);
    discamb::Vector3d cartesian;
    unitCell.fractionalToCartesian(c1.coordinates, cartesian);
    cout << c1.type << " " << cartesian[0] << " " << cartesian[1]
         << " " << cartesian[2] << std::endl;

    //###################   define space group   ##############################

    // default constructor creates identity operation
    discamb::SpaceGroupOperation identity;

    // space group operations can be set using either matrix - vector notation
    discamb::SpaceGroupOperation fourBar(discamb::Matrix3i{ 0,  1,  0,
                                                           -1,  0,  0,
                                                            0,  0, -1 });

    // or string x,y,z notation
    discamb::SpaceGroupOperation mirrorPlane("1/2-Y, 1/2-X, Z");
    // or alternatively 
    mirrorPlane.set({ 0, -1,  0,
                     -1,  0,  0,
                      0,  0,  1 }, { 1 / 2, 1 / 2, 0 });


    // lets generata all the operations necessary to define space group P -4 21 m in DiSCaMB
    vector< discamb::SpaceGroupOperation > operations =
        { identity                    , mirrorPlane * identity                   ,
          fourBar                     , mirrorPlane * fourBar                    ,
          fourBar * fourBar           , mirrorPlane * fourBar * fourBar          ,
          fourBar * fourBar * fourBar , mirrorPlane * fourBar * fourBar * fourBar };

    discamb::SpaceGroup p4bar21m;
    p4bar21m.set(operations);

    //############################  define crystal  ###########################

    discamb::Crystal crystal;
    crystal.adpConvention = discamb::structural_parameters_convention::U_cif;
    // we have previously omitted atomic parameters specification for atoms other than c1
    // atoms redefined later anyway when XD files are read
    crystal.atoms = { c1, o1, n1, h1, h2 }; 
    crystal.spaceGroup = p4bar21m;
    crystal.unitCell = unitCell;

    crystal.xyzCoordinateSystem = discamb::structural_parameters_convention::fractional; 

    //###########################

    // set wave function related data

    discamb::DeformationValenceParameters deformationValenceParameters;
    discamb::HC_WfnType h_wfn, c_wfn, o_wfn, n_wfn;
    discamb::ClementiRoettiData clementiRoettiData;

    c_wfn = clementiRoettiData.getEntry("C");
    deformationValenceParameters.getParameters("C", c_wfn.deformation_valence_exponent,
    c_wfn.deformation_valence_power);

    // set data specific for atom type

    discamb::HC_AtomTypeParameters c_type, n_type, h_type, o_type;

    c_type.kappa_spherical_valence = 0.99;
    c_type.kappa_deformation_valence = 0.835;
    c_type.p_val = 4.2488;

                  // P_{0,0} 
    c_type.p_lm = { { 0.0000 },
                  // P_{1,-1}, P_{1,0}, P_{1,1}
                    { 0.0000,  0.0203,  0.0000 },
                  // P_{2,-2}, P_{2,-1}, P_{2,0}, P_{2,1}, P_{2,2}
                    { 0.0000,  0.0000,  0.0578,  0.0000,  0.0532 },
                    { 0.0000,  0.0000,  0.0000,  0.1060,  0.0000, -0.0742,  0.0000 },
                    { 0.0000,  0.0000,  0.0000,  0.0000, -0.0120,  0.0000,  0.0275,  0.0000,  0.0043 } };

    // define multipole model parameters

    discamb::HC_ModelParameters multipoleModelParameters;

    multipoleModelParameters.wfn_parameters = { c_wfn, o_wfn, n_wfn, h_wfn };
    multipoleModelParameters.atom_to_wfn_map = { 0, 1, 2, 3, 3 };
    multipoleModelParameters.type_parameters = { c_type, o_type, n_type, h_type };
    multipoleModelParameters.atom_to_type_map = { 0, 1, 2, 3, 3 };

    // alternatively it can be read from XD files

    vector<discamb::XdLocalCoordinateSystem> lcs;
    discamb::xd_io::read("xd.mas", "xd.inp", multipoleModelParameters, crystal, lcs, false, true);

    //####################### structure factor calculation

    // lets define some hkl indices
    vector<discamb::Vector3i> hkl{ { 1, 0, 0 },{ 1, -1, 2 },{ 1, 2, 3 } };


    // set coordinate system for the first atom in crystal.atoms
    // just an example definition not included in further calculations

    discamb::XdLocalCoordinateSystem lcs_0;

    lcs_0.set("C()", "O()", "C()", "DUM1", "Z", "X", true, crystal, false, { { -0.144680, 0.355320, 0.179010 } });

    // and get local coordinates for current geometry

    vector< discamb::Matrix3d > lcsMatrices(5);

    lcs[0].calculate(lcsMatrices[0], crystal); // .. and so on for other atoms

    for(size_t i=1;i<5;i++)
        lcs[i].calculate(lcsMatrices[i], crystal);

    //#################################### macromolecular-like setup of structure factors and derivatives calc

    // the derivatives of target function with respect to structure factor components
    vector<complex<double> > dT_dF{ { 1.2, -0.34 },{ 0.13, 3.1 },{ 1, 2 } };

    // the deravatives to be computed - derivatives of target function
    // with respect to structural parameters
    std::vector< discamb::TargetFunctionAtomicParamDerivatives > derivatives;
    
    // container for structure factors
    vector<complex<double> > structureFactors;
    
    // set up the calculator
    discamb::HansenCoppensStructureFactorCalculator calculator(crystal, multipoleModelParameters);
    
    // calculate
    calculator.calculateStructureFactorsAndDerivatives(crystal.atoms, lcsMatrices,
        hkl, structureFactors, derivatives,
        dT_dF, vector<bool>(5, true));
    
    //################## small molecule setup for structure factor calculation
    
    DataCollector collector;
    
    // get number of hkl vectors
    size_t nHkl = hkl.size();
    
    // and allocate memory for collector structures 
    collector.structureFactorsData.resize(nHkl);
    collector.derivativesData.resize(nHkl);
    
    // the same calculator as for macromolecular-like setup is used
    
    
    // calculate, the calculator already set up in the 'macromolecular' part above
    
    calculator.calculateStructureFactorsAndDerivatives(crystal.atoms, lcsMatrices,
        hkl, collector);
    
    // print out structure factors
    
    cout << "h k l and structure factor value" << endl;
    
    for (size_t i = 0; i<nHkl; i++)
        cout << hkl[i][0] << " " << hkl[i][1] << " " << hkl[i][2]
        << " " << collector.structureFactorsData[i] << endl;
    
}

