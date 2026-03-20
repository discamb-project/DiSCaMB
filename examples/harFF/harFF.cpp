#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/rsp_io.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/StockholderAtomFormFactorCalcManager.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"


#include "json.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>

using namespace discamb;
using namespace std;

void readSettings(
    const std::string &jsonFile,
    const Crystal& crystal,
    nlohmann::json &jsonData
    //HirshfeldAtomModelSettings& harSettings
)
{
    
   
    ifstream jsonFileStream("aspher.json");

    if (jsonFileStream.good())
        jsonFileStream >> jsonData;
    else
        on_error::throwException("cannot read file '" + jsonFile + "'", __FILE__, __LINE__);

    jsonFileStream.close();
    
    //harSettings.set(jsonData, crystal);
}



int main(int argc, char *argv[])
{


    try {
        try {

            Crystal crystal;
            int nAtoms;

            ofstream logFile("harFF.log");
            auto clog_orginal_rdbuf = std::clog.rdbuf();
            std::clog.rdbuf(logFile.rdbuf());
            
            if (argc != 5)
                on_error::throwException("expected arguments: (1) structure file (cif/ins/res) (2) rsp file (3) settings file (4) output file", __FILE__, __LINE__);

            filesystem::path structureFilePath(argv[1]);
            if (!filesystem::exists(structureFilePath))
                on_error::throwException("declared structure file missing", __FILE__, __LINE__);
           
            structure_io::read_structure(argv[1], crystal);

            vector<Vector3d> q_vectors;
            rsp_io::read(argv[2], q_vectors);
                       
            HirshfeldAtomModelSettings hamSettings;
            nlohmann::json settings;


            nAtoms = crystal.atoms.size();
            vector< vector<complex<double> > > formFactors;
            vector<bool> includeAtom(nAtoms, true);
            readSettings(argv[3], crystal, settings);

            if (settings.value("model", string("iam")) == string("iam"))
            {
                string table = settings.value("table", string("Waasmeier-Kirfel"));
                IamFormFactorCalculationsManager iamCalculator(crystal, table);
                iamCalculator.calculateFrac(q_vectors, formFactors, includeAtom);
            }
            else
            {
                if (settings.value("model", string("iam")) == string("taam"))
                {

                }
                else
                {
                    HirshfeldAtomModelSettings harSettings;
                    harSettings.set(settings, crystal);
                    StockholderAtomFormFactorCalcManager harCalculator(crystal, harSettings);
                    harCalculator.calculateAtomicDensities();
                    harCalculator.calculateFrac(q_vectors, formFactors, includeAtom);
                }
            }
            
            
            vector<string> atomLabels;
            for (auto const& atom : crystal.atoms)
                atomLabels.push_back(atom.label);

            tsc_io::write_tscd(argv[4], atomLabels, q_vectors, formFactors);

            std::clog.rdbuf(clog_orginal_rdbuf);
            logFile.close();

        }
        catch (nlohmann::json::parse_error& e)
        {
            std::cout << "message: " << e.what() << '\n'
                << "exception id: " << e.id << '\n'
                << "byte position of error: " << e.byte << std::endl;

            stringstream ss;

            ss << "message: " << e.what() << '\n'
               << "exception id: " << e.id << '\n'
               << "byte position of error: " << e.byte << std::endl;

            string message = "Error when parsing JSON file\n" + ss.str();
            on_error::throwException(message, __FILE__, __LINE__);
        }
    }
    catch (exception &e)
    {
        cout << e.what() << endl;
        ofstream out("discamb_error.log" , ostream::app);
        out << e.what();
        out.close();
    }
}


