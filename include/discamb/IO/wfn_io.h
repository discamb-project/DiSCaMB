#pragma once

#include <string>
#include <istream>
#include <map>
#include <vector>
#include "discamb/MathUtilities/Vector3.h"

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */


    namespace wfn_io {

        struct AdditionalElectronDensity
        {
            std::vector<int> primitive_to_center;
            std::vector<int> primitive_type;
            std::vector<double> primitive_exponents;
            std::vector<double> primitive_coefficients;
            std::string name;
        };

        struct WfnFileData {
            std::string label;
            bool isGaussian = true;
            std::vector<int> atomic_numbers;
            std::vector<std::string> center_label;
            std::vector<Vector3d> center_position;
            std::vector<double> center_charge;

            /** numeration starts from 1 */
            std::vector<int> primitive_to_center;
            std::vector<int> primitive_type;
            std::vector<double> primitive_exponents;
            std::vector<double> molecular_orbital_energy;
            std::vector<double> molecular_orbital_occupancy;
            std::vector<std::vector<double> > molecular_orbitals;
            double energy = 0.0;
            double virial = 0.0;
            
            // Additional Electron Density Functions
            std::vector<AdditionalElectronDensity> edfs;
        };

        void read_wfn(const std::string &fileName, WfnFileData &data, bool skipNonAtomicCenters=true);
        void read_wfn(std::istream &input, WfnFileData &data, bool skipNonAtomicCenters=true);
        void read_atomic_wfn_database(const std::string &fileName, std::vector<WfnFileData> &data);

        void read_wfx(const std::string& fileName, WfnFileData& data);
        void read_wfx(std::istream& input, WfnFileData& data);

        void read_wavefunction(const std::string& fileName, WfnFileData& data);
        

        void read_wfx(const std::string& fileName, std::map<std::string, std::vector<std::string> >& data,
                     std::vector< std::map<std::string, std::vector<std::string> > > &edfs);
        void read_wfx(std::istream& input, std::map<std::string, std::vector<std::string> >& data,
                      std::vector< std::map<std::string, std::vector<std::string> > >& edfs);
        // merge additional density functions
        void mergeEdfs(const std::vector<AdditionalElectronDensity>& edfs, AdditionalElectronDensity &edf);

        //void removeCentersWithNoBasisFunctions(WfnFileData& data);
    }
    /**@}*/
}
