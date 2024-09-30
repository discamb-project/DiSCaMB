#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/HC_Model/XdLocalCoordinateSystem.h"

namespace discamb{

    /**
    * \addtogroup IO IO
    * @{
    */


struct XdMasData
{
    UnitCell unitCell;
    SpaceGroup spaceGroup;
    std::vector<std::string> atomLabels;
    std::vector<std::vector<std::string> > lcs;
    std::vector<Vector3d> dummy_atoms;
    bool dummyAtomIndexingFrom0 = true;
    std::vector<std::vector<std::string> > scat_tables;
    std::vector<int> atomKappaSetIndex;
    std::vector<int> atomScatTableIndex;
};

struct XdParData
{
    std::vector<Vector3d> positionsFractional;
    std::vector<std::string> atomLabels;
    std::vector<std::vector<double> > adps;
    std::vector<std::vector<double> > kappa_sets;
    std::vector<std::vector<std::vector<double> > > plms;
    std::vector<double> p_val;
    std::vector<int> atomToKappaSetMap;
    std::vector<int> kappaSetToScatTableMap;
    std::vector<Vector3d> dummyAtomsPositionsFractional;
    std::vector<double> occupancy;
};

// error: always crystal.atoms[i].multiplicity = 1;



namespace xd_io{

//  -discamb
    void read(const std::string & masterFileName,const std::string & parameterFileName,Crystal &crystal);
    
    void readXdMas(const std::string & masterFileName,XdMasData &xdMasData);

    void readXdPar(const std::string & parameterFileName,XdParData &xdParData);

    void readXdMasSections(const std::string & masterFileName,std::vector<std::string> &generalSection,
                          std::vector<std::pair<std::string,std::vector<std::string> > > &sections);

    void getCrystal(const XdMasData &xdMasData,const XdParData &xdParData,Crystal &crystal);

    void read(const std::string & masterFileName,const std::string & parameterFileName,
                     HC_ModelParameters &hc_parameters,Crystal &crystal,
                     std::vector<XdLocalCoordinateSystem> &localCoordinateSystems,bool group_types,bool comparePval=true);
 //  -discamb

    void write(const std::string &xdMasName, const std::string &xdParName, 
               const Crystal &c, const HC_ModelParameters &params,
               const std::vector<XdLocalCoordinateSystem> &lcs);

    void writeMasterFile(
        const std::string &xdMasName,
        const Crystal &c,
        const HC_ModelParameters &params,
        const std::vector<XdLocalCoordinateSystem> &lcs);

#ifdef _MSC_VER

    void readXdFou(const std::string &fileName,
                   std::vector<Vector3i> &hkl,
                   std::vector<double> &f_obs,
                   std::vector<double> &sigma,
                   std::vector<double> &phase_obs,
                   std::vector<std::complex<double> > &f_model_1,
                   std::vector<std::complex<double> > &f_model_2);

    void writeXdFou(const std::string &fileName,
        const std::vector<Vector3i> &hkl,
        const std::vector<double> &f_obs,
        const std::vector<double> &sigma,
        const std::vector<double> &phase_obs,
        const std::vector<std::complex<double> > &f_model_1,
        const std::vector<std::complex<double> > &f_model_2);
#endif
}
/**@}*/
}

