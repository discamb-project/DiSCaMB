#include "discamb/QuantumChemistry/ProatomDB.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/IO/proatom_db_io.h"
#include "discamb/BasicUtilities/discamb_env.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "json.hpp"


#include <fstream>

#include <filesystem>
#include <map>

using namespace std;
using namespace nlohmann;
//namespace fsystem = discamb::file_system;

namespace discamb {

    namespace {
        void getAlternativeBasisSetName(
            map<string, string>& alternativeBasisSetName)
        {
            alternativeBasisSetName.clear();
            string discambPathStr = discamb_env::get_discamb_path();
            if (discambPathStr.empty())
                return;
            
            filesystem::path filePath = filesystem::path(discambPathStr) / "data/atomic_densities/txt/_alternative_basis_set_names.json";

            if (!filesystem::exists(filePath))
                return;
                        
            nlohmann::json jsonData;
            string fPath = filePath.string();
            ifstream jsonFileStream(fPath);

            if (jsonFileStream.good())
                jsonFileStream >> jsonData;
            else
                on_error::throwException("can not read file " + fPath, __FILE__, __LINE__);
            vector<string> alternativeNames;
            for (auto [key, value]: jsonData.items()) {
                string_utilities::split(value, alternativeNames, ',');
                for (auto basisSetName : alternativeNames)
                    alternativeBasisSetName[string_utilities::toLower(basisSetName)] = key;
                alternativeBasisSetName[string_utilities::toLower(key)] = key;
            }

        }

        //void getAvailableBasisSets(
        //    const string& method,
        //    vector<string>& basisSetNames)
        //{
        //    basisSetNames.clear();
        //    for (auto it : filesystem::directory_iterator(filesystem::current_path()))
        //        if(filesystem::is_regular_file(it.path()))
        //        {
        //            string s = it.path().filename().string();
        //            if (s.find("atoms_" + method) != string::npos)
        //            {

        //            }
        //        }
        //}

    }

    ProatomDB::ProatomDB()
    {   
    }

    ProatomDB::~ProatomDB()
    {
    }

    bool ProatomDB::hasDiscambDataFile(
        const std::string& method,
        const std::string& basisSet,
        std::string& fileName)
    {
        clog << "attempting to find spherical atomic densities file\n";
        
        string discambPathStr = discamb_env::get_discamb_path();
        if (discambPathStr.empty())
        {
            clog << "  not available, DISCAMB_PATH not set\n";
            return false;
        }

        

        filesystem::path discambPath(discamb_env::get_discamb_path());

        string unifiedBasisSetName = basisSet;
        map<string, string> alternativeBasisSetNames;
        getAlternativeBasisSetName(alternativeBasisSetNames);
        string basisSetLowerCase = string_utilities::toLower(basisSet);
        if (alternativeBasisSetNames.find(basisSetLowerCase) != alternativeBasisSetNames.end())
            unifiedBasisSetName = alternativeBasisSetNames.find(basisSetLowerCase)->second;


        filesystem::path atomFilePath = discambPath /
            ("data/atomic_densities/txt/atoms_" + method + "_" + unifiedBasisSetName + ".txt");
        clog << "  " << "looking for file " << atomFilePath.string() << "\n";
        if (!filesystem::exists(atomFilePath))
        {
            clog << "  file not found\n";
            return false;
        }
        clog << "  found\n";
        fileName = atomFilePath.string();
        return true;
    }


    void ProatomDB::setFromFile(
        const std::string& fileName)
    {
        vector<int> z;
        vector<int> q;
        vector<vector<double> > values;

        proatom_db_io::read_proatom_db(fileName, z, q, values);
        //std::map<std::pair<int, int>, std::vector<double> > mDensities;
        mDensities.clear();
        int i, n = z.size();
        for (i = 0; i < n; i++)
            mDensities[{z[i], q[i]}] = values[i];

    }

    bool ProatomDB::hasEntry(
        int atomicNumber,
        int charge)
        const
    {
        return (mDensities.find({ atomicNumber, charge }) != mDensities.end());
    }

    bool ProatomDB::getSphericalAtom(
        int atomicNumber,
        int charge,
        SphericalAtomicDensity& density)
        const
    {
        if (hasEntry(atomicNumber, charge))
        {
            density.setValues(mDensities.find({ atomicNumber, charge })->second, 0.001);
            return true;
        }
        return false;
    }
    
    
    
}
