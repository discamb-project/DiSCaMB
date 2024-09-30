#include "discamb/IO/NativeIAM_Reader.h"

#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
//#include "discamb/CrystalStructure/crystal_structure_utilities.hpp"
#include <fstream>
#include <algorithm>

using namespace std;

namespace discamb {

NativeIAM_Reader::NativeIAM_Reader()
{
}

NativeIAM_Reader::~NativeIAM_Reader()
{
}


void NativeIAM_Reader::read(
const std::string &fileName, 
Crystal &crystal,
ScatteringParameters &scatteringParameters)
const
{
    ifstream in(fileName.c_str());
    string s,s2,line;

    vector<string> words;
    int format_type; // 1 or 2
    int n_symm;
    vector<SpaceGroupOperation> space_group_operations;
    vector<string> atomScatteringType;

    
    double a,b,c,alpha,beta,gamma;
    Vector3d vec_3d;
    Vector3<CrystallographicRational> inv_center_translation_rational_vec;
    char centering;
    bool hasInversionCenter,correct_translation_vector;

    if(!in.good())
        on_error::throwException(string("can not read file '")+fileName,__FILE__,__LINE__);

    in>> s >> s >> a >> b >> c >> alpha >> beta >> gamma;

    crystal.unitCell.set(a, b, c, alpha, beta, gamma);
    getline(in,line);
    getline(in,line);
    in>> s;
    if(s == string("centering"))
    {
        /*
        centering  P
        inversion centre: yes inversion centre translation:  0.5 0.5 0.0
        */

        format_type = 1;
        in >> centering;
        in >> s >> s >> s2;
        if(s2==string("yes"))
        {
            hasInversionCenter=true;
            in>> s >> s >> s >> vec_3d[0] >> vec_3d[1] >> vec_3d[2];

            inv_center_translation_rational_vec = 
                SpaceGroupOperation::crystallographicTranslation(vec_3d, correct_translation_vector);
        }
        else
            hasInversionCenter = false;
        
        in >> s >> s >> s >> n_symm;
    }
    else
    {
        format_type = 2;
        in >> s >> s >> n_symm;
    }

    // n symmetry operations  8
    
    getline(in,line);
    space_group_operations.resize(n_symm);
    for(int i=0;i<n_symm;i++)
    {
        getline(in, line);
        space_group_operations[i].set(line);
    }

    if(format_type == 1)
        crystal.spaceGroup.set(space_group_operations, centering, true, hasInversionCenter, inv_center_translation_rational_vec);
    else
        crystal.spaceGroup.set(space_group_operations);

    /*
    n scatterers  33
    C1   C      4    0.7805200   -0.2496000    0.8659300 1.000000 u_aniso  0.028900 0.018800 0.024600 -0.000700 0.000400 -0.001300
    H11  H      4    0.8682000   -0.2443000    0.8342000 1.000000 u_iso 0.029400
    */

    int n_scatterers;
    in>> s >> s >> n_scatterers;
    crystal.atoms.resize(n_scatterers);

    double x,y,z,occupancy,u_iso,u11,u22,u33,u12,u13,u23;
    string atom_label,u_type;
    int multiplicity;
    bool is_u_iso;
    atomScatteringType.resize(n_scatterers);

    for(int i=0;i<n_scatterers;++i)
    {
        in>> atom_label >> atomScatteringType[i] >> multiplicity >> x >> y >> z >> occupancy;
        in>> u_type;
        is_u_iso = (u_type == string("u_iso"));
        if(is_u_iso)
        {
            in>> u_iso;
            crystal.atoms[i].adp.assign(1, u_iso);
        }
        else
        {
            in>> u11 >> u22 >> u33 >> u12 >> u13 >> u23;
            //C++11 crystal.atoms[i].adp = {u11,u22,u33,u12,u13,u23};
            crystal.atoms[i].adp.resize(6);
            crystal.atoms[i].adp[0] = u11;
            crystal.atoms[i].adp[1] = u22;
            crystal.atoms[i].adp[2] = u33;
            crystal.atoms[i].adp[3] = u12;
            crystal.atoms[i].adp[4] = u13;
            crystal.atoms[i].adp[5] = u23;
        }
        crystal.atoms[i].coordinates=Vector3d(x,y,z);
        crystal.atoms[i].label = atom_label;
        crystal.atoms[i].multiplicity = multiplicity;
        crystal.atoms[i].occupancy = occupancy;
        crystal.atoms[i].type = atomScatteringType[i];
    }
    // sets ScatteringParameters::atomTypeSymbols and ScatteringParameters::atomToTypeMap

    scatteringParameters.groupAtomTypes(atomScatteringType,scatteringParameters.atomTypeSymbols,scatteringParameters.atomToTypeMap);

    // read anomalous dispersion if available
    int nAtomTypes = scatteringParameters.atomTypeSymbols.size();
    scatteringParameters.atomTypeAnomalousDispersion.assign(nAtomTypes,0.0);
    getline(in, line);
    getline(in, line);
    string_utilities::split(line,words);
    vector<string>::iterator it;
    
    if(!words.empty()) //anomalus dispersion 
    {
        getline(in,line);
        string_utilities::split(line, words);

        while(in.good())
            if(!words.empty())
            {
                if(words.size()!=3)
                    on_error::throwException(string("Error when reading anomaloous dispersion data from file '")+
                                             fileName + string(" '"),__FILE__,__LINE__);
            
                it = find(scatteringParameters.atomTypeSymbols.begin(), scatteringParameters.atomTypeSymbols.end(), words[0]);
                if(it == scatteringParameters.atomTypeSymbols.end())
                    on_error::throwException(string("Error when reading anomaloous dispersion data from file '") +
                                             fileName + string(" ', unknown atom type: '") + words[0] + string("'"),
                                             __FILE__, __LINE__);
                int index = distance(scatteringParameters.atomTypeSymbols.begin(),it);

                scatteringParameters.atomTypeAnomalousDispersion[index]  = complex<double>(atof(words[1].c_str()),atof(words[2].c_str()));
                getline(in, line);
                string_utilities::split(line, words);
            }
        
    }

    in.close();

    StructuralParametersConverter conv(crystal.unitCell);
    vector<double> _adp;
    for (auto& atom : crystal.atoms)
        if (atom.adp.size() == 6)
        {
            _adp = atom.adp;
            conv.uStarToU_cif(_adp, atom.adp);
        }


    crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
    crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;

}

}// namespace discamb
