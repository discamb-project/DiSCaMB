#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/IO/structure_io.h"

#include "json.hpp"

#include <iostream>
#include <fstream>
#include <filesystem>

using namespace discamb;
using namespace std;

struct TscSource {
    string name;
    double weight;
    vector<string> atoms;
};


int main(int argc, char *argv[])
{


    try {
        try {

                      
            if (argc != 2)
                on_error::throwException("expected settings file as an argument", __FILE__, __LINE__);
                                    
            // read json file
            
            string jsonFile = argv[1];
            nlohmann::json settings;
            ifstream jsonFileStream(jsonFile);


            if (jsonFileStream.good())
                jsonFileStream >> settings;
            else
                on_error::throwException("cannot read file '" + jsonFile + "'", __FILE__, __LINE__);

            jsonFileStream.close();

            // process tsc

            string outputTsc = settings.value("output tsc", string());
            if (outputTsc.empty())
                on_error::throwException("missing \"output tsc\" entry in settings file", __FILE__, __LINE__);
            string structureFile = settings.value("structure", string());
            if (structureFile.empty())
                on_error::throwException("missing \"structure\" entry in settings file", __FILE__, __LINE__);

            vector<TscSource> tscSources;
            
            if(settings.find("form factors") == settings.end())
                on_error::throwException("missing \"form factors\" entry in settings file", __FILE__, __LINE__);
            else
            {
                auto ff = settings.find("form factors");
                for (auto item : *ff)
                {
                    TscSource tscSource;
                    tscSource.name = item.value("file", string());
                    if(tscSource.name.empty())
                        on_error::throwException("missing \"file\" entry in settings file, \"form factors\" part", __FILE__, __LINE__);
                    tscSource.weight = item.value("weight", 1.0);
                    
                    
                    string listStr = item.value("atoms", string());
                    string_utilities::split(listStr, tscSource.atoms, ',');
                    
                    tscSources.push_back(tscSource);
                }
            }

            // 
            Crystal crystal;
            structure_io::read_structure(structureFile, crystal);
            vector<string> atomLabels;
            for (auto& atom : crystal.atoms)
                atomLabels.push_back(atom.label);
            int nAtoms = atomLabels.size();
            vector<double> weights_sum(nAtoms, 0.0);
            //
            vector<Vector3i> hkl;
            bool initialized = false;
            vector<vector<complex<double> > > finalFf;
            int nHkl=0;
            for (auto tscSource : tscSources)
            {
                vector<string> tscAtomLabels;
                vector<vector<complex<double> > > ff;
                tsc_io::read_tsc(tscSource.name, tscAtomLabels, hkl, ff);
                if (!initialized)
                {
                    initialized = true;
                    nHkl = hkl.size();
                    finalFf.resize(nHkl, vector<complex<double> >(nAtoms, { 0,0 }));
                }
                if (tscSource.atoms.empty())
                    tscSource.atoms = tscAtomLabels;

                for (auto& atomLabel : tscSource.atoms)
                {
                    auto it = find(atomLabels.begin(), atomLabels.end(), atomLabel);
                    int idxInCrystal;
                    if (it != atomLabels.end())
                        idxInCrystal = std::distance(atomLabels.begin(), it);
                    else
                        on_error::throwException("atom with label '" + atomLabel + "' present in tsc file but absent in structure file", __FILE__, __LINE__);
                    int idxInTsc;
                    it = find(tscAtomLabels.begin(), tscAtomLabels.end(), atomLabel);

                    if (it != tscAtomLabels.end())
                        idxInTsc = std::distance(tscAtomLabels.begin(), it);
                    else
                        on_error::throwException("atom with label '" + atomLabel + "' missing in tsc file '" + tscSource.name +"'", __FILE__, __LINE__);

                    weights_sum[idxInCrystal] += tscSource.weight;
                    
                    for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
                        finalFf[hklIdx][idxInCrystal] += ff[hklIdx][idxInTsc];
                }
            }

            //

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (weights_sum[atomIdx] == 0)
                    on_error::throwException("atom '" + atomLabels[atomIdx] + "' has no assigned form factors", __FILE__, __LINE__);

            // 

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
                    finalFf[hklIdx][atomIdx] /= weights_sum[atomIdx];

            tsc_io::write_tsc(outputTsc, atomLabels, hkl, finalFf);

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


