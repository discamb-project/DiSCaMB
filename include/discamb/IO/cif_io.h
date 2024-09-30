#ifndef _DISCAMBDEV_IO_CIFIO_H_
#define _DISCAMBDEV_IO_CIFIO_H_

#include "discamb/Scattering/ScatteringData.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/AtomTyping/LocalCoordinateSystem.h"
#include "discamb/BasicUtilities/string_utilities.h"



namespace discamb {

/**
* \defgroup IO IO
\brief Reading and writeng files in various formats.
* @{
*/



    namespace cif_io {

        struct Loop {
            // [key index]
            std::vector<std::string> loop_headers;
            /** [key index][value index] */
            std::vector<std::vector<std::string> > values;
            std::string to_string(int nItemsPerLine=-1) const;
        };

        struct DataSet {
            /** e.g. data_urea i.e. data_ is also included*/
            std::string dataSetName;
            std::map < std::string, std::string> key_value_items;
            std::vector<Loop> loops;
        };


        void saveCif(const std::string &fName, const Crystal &crystal);
        
        void saveCif(
            const std::string& fName, 
            const Crystal& crystal, 
            const HC_ModelParameters& parameters, 
            const std::vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs);

        void readCif(const std::string &fName, std::vector<DataSet> &data);
        /* WARNING! if _atom_site_symmetry_multiplicity is absent atomic multiplicities are set to 0 */
        void cifDataToCrystal(const DataSet &data, Crystal &crystal);
        void readFcf(const std::string &fName, SpaceGroup &sg, UnitCell &uc, std::vector<Vector3i> &hkl,
                     std::vector<double> &_refln_F_meas, std::vector<double> &_refln_F_sigma,
                     std::vector<double> &_refln_A_calc, std::vector<double> &_refln_B_calc);
//        void writeFcf(const std::string& fName, const SpaceGroup& sg, const UnitCell& uc, const std::vector<Vector3i>& hkl,
//            const std::vector<double>& _refln_F_meas, const std::vector<double>& _refln_F_sigma,
//            const std::vector<double>& _refln_A_calc, const std::vector<double>& _refln_B_calc,
//            bool scale = false, double a=0, double b=0);

        void writeFcf(const std::string& fName, const SpaceGroup& sg, const UnitCell& uc, const std::vector<Vector3i>& hkl,
            const std::vector<double>& _refln_F_squared_meas, const std::vector<double>& _refln_F_squared_sigma,
            const std::vector<double>& _refln_F_squared_calc, bool scale = false, double a = 0, double b = 0);

        void extractAtomSiteData(const DataSet &data, std::vector<AtomInCrystal> &atoms);
        void extractSpaceGroupData(const DataSet &data, SpaceGroup &spaceGroup);
        void extractUnitCellData(const DataSet &data, UnitCell &unitCell);

        void extractReflectionData(const DataSet &data, ReflectionsData &reflectionsData, std::map<std::string, double> &reflns_scale_meas_intensity);

        template<typename T>
        void convertVectorOfStrings(const std::vector<std::string> stringData, std::vector<T> &convertedData);

        template<>
        void convertVectorOfStrings<char>(const std::vector<std::string> stringData, std::vector<char> &convertedData);


        /** looks for loop with entry given by \p tag*/
        bool findLoopIndex(const DataSet &data, const std::string &tag, int &loopIdx, int &headerIdx);
        bool findLoopTagIndex(const Loop &loop, const std::string &tag, int &headerIdx);
        /** looks for loop containing entries given by \p tags assuming that there is no tab which appears in two different loops*/
        bool findLoopIndex(const DataSet &data, const std::vector<std::string> &tags, int &loopIdx, 
                           std::vector<int> &tagsIndices);
        

        // Hansen-Coppens multipole model

        void extractMultipoleModel(
            const Crystal &crystal,
            const DataSet& data, 
            HC_ModelParameters& parameters,
            std::vector < LocalCoordinateSystem<AtomInCrystalID> > &lcs);

        void saveMultipoleModel(
            const Crystal& crystal,
            DataSet& data,
            std::vector<Vector3d> &dummyAtomFractionalPosition,
            const HC_ModelParameters& parameters,
            const std::vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs);

    }
    /**@}*/
}

template<typename T>
void discamb::cif_io::convertVectorOfStrings(
    const std::vector<std::string> stringData,
    std::vector<T> &convertedData)
{
    int i, n = stringData.size();
    convertedData.resize(n);
    for (i = 0; i < n; i++)
        discamb::string_utilities::convertFromString(stringData[i], convertedData[i]);
}

#endif

