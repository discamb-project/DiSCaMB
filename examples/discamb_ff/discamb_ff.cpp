#include "discamb/BasicUtilities/Timer.h"
#include "discamb/IO/hkl_io.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/IO/shelx_io.h"
#include "discamb/IO/tsc_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/Scattering/SfCalculator.h"
#include "discamb/BasicUtilities/parse_cmd.h"

#include "json.hpp"


#include <cstdio>
#include <fstream>
#include <memory>
#include <sstream>
#include <set>
#include <filesystem>

using namespace discamb;
using namespace std;




void makeWholeHklSet(
    const vector<Vector3i>& hkl0, 
    const SpaceGroup &spaceGroup, 
    vector<Vector3i> &hkl)
{
    int nSymm = spaceGroup.nSymmetryOperationsInSubset();
    vector<Matrix3i> rotations(nSymm);

    for (int i = 0; i < nSymm; i++)
        spaceGroup.getSpaceGroupOperation(0, 0, i).getRotation(rotations[i]);
    for (int i = 0; i < nSymm; i++)
        rotations.push_back(-1 * rotations[i]);


	hkl.clear();
	set<Vector3i> uniqueHkl;
    vector<Vector3i> hkl1 = hkl0;

	for (auto const& rotation : rotations)
        for (auto const& h : hkl1)
            uniqueHkl.insert(h * rotation);

	hkl.assign(uniqueHkl.begin(), uniqueHkl.end());
}




std::unique_ptr<SfCalculator> sfCalculatorFromJsonFile(
    const Crystal &crystal)
{
    nlohmann::json jsonData;
    
    if (filesystem::exists(filesystem::path("aspher.json")))
    {
        ifstream jsonFileStream("aspher.json");

        if (jsonFileStream.good())
            jsonFileStream >> jsonData;
        else
        {
            on_error::throwException("can not read aspher.json file, expected to be present in the current directory", __FILE__, __LINE__);
            return std::unique_ptr<SfCalculator>(nullptr);
        }

    }
    else
        on_error::throwException("no aspher.json file", __FILE__, __LINE__);

    return std::unique_ptr<SfCalculator>(SfCalculator::create(crystal, jsonData));
}

void readHkl(
    const string& hklFile,
    vector<Vector3i>& hkls,
    bool shelxFreeFormat)
{
    hkls.clear();
    vector<int> batchNumbers;
    vector<double> intensities, sigmas;

    if (filesystem::path(hklFile).extension().string() == string(".tsc"))
    {
        vector<string> atomLabels;
        vector<vector<complex<double> > > ff;
        tsc_io::read_tsc(hklFile, atomLabels, hkls, ff);
        hkls.push_back(Vector3i(0, 0, 0));
    }
    else
        hkl_io::readShelxHkl(hklFile, hkls, intensities, sigmas, batchNumbers, shelxFreeFormat);

}

int main(int argc, char *argv[])
{
    Crystal crystal;
    
    ofstream out;
    vector<string> arguments; 
    vector<string> options;
    map<string, string> optionsWithValues;
    


    try {
        try {
            WallClockTimer timerAll;
            timerAll.start();
            ofstream logFile("discamb2tsc.log");
            auto clog_orginal_rdbuf = std::clog.rdbuf();
            std::clog.rdbuf(logFile.rdbuf());

            time_t time_now = std::chrono::system_clock::to_time_t(chrono::system_clock::now());

            string header = string("\n  discamb_ff") + string("\n  compiled on ") +
                __DATE__ + string(" , ") + __TIME__ + string(".\n\n");

            cout << header;

            // read crystal data

            string structureFile, hklFile;

            parse_cmd::get_args_and_options(argc, argv, arguments, options, optionsWithValues);

            if (arguments.size() != 2)
            {
                cout << "expected two arguments: structure file (cif/ins) and hkl file (shelx hkl or tsc)\n";
                exit(0);
            }
            else
            {
                structureFile = arguments[0];
                hklFile = arguments[1];
            }

            bool shelxFreeFormat = (find(options.begin(), options.end(), "-f") != options.end());
            structure_io::read_structure(structureFile, crystal);

            vector<Vector3i> hkls, hklAll;

            readHkl(hklFile, hkls, shelxFreeFormat);

            makeWholeHklSet(hkls, crystal.spaceGroup, hklAll);

            WallClockTimer timer;

            timer.start();

            clog << "file created at " << ctime(&time_now) << "\n";
           
            
            std::unique_ptr<SfCalculator> calculator = 
                sfCalculatorFromJsonFile(crystal);
            

            string stem = filesystem::path(structureFile).stem().string();
            string outputFileName = stem + string(".tsc");
            ofstream out(outputFileName);
            
            vector<pair<string, string> > modelInfo;
            calculator->getModelInformation(modelInfo);

            vector<string> atomLabels;
            for (auto& atom : crystal.atoms)
                atomLabels.push_back(atom.label);

            timer.start();

            vector< vector<complex<double> > > formFactors;
            vector<bool> includeAtom(crystal.atoms.size(), true);
            calculator->calculateFormFactors(hklAll, formFactors, includeAtom);

            clog << "calculation of form factors: " << timer.stop() << " ms\n";

            string tsc_comment = "FORM FACTORS SOURCE:\n";
            for (auto& p : modelInfo)
                tsc_comment += "   " + p.first + " - " + p.second + "\n";
            tsc_io::write_tsc(outputFileName, atomLabels, hklAll, formFactors, tsc_comment);

            out.close();

            clog << "total time: " << timerAll.stop() << " ms\n";
            
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


