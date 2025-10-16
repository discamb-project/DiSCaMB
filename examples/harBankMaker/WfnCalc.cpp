#include "WfnCalc.h"
#include "ChooseStructures.h"

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/Timer.h"
#include "discamb/IO/mol2_io.h"
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include <fstream>
#include <memory>

using namespace std;
using namespace discamb;

WfnCalc::WfnCalc()
{

}

WfnCalc::~WfnCalc()
{
}

void WfnCalc::set() 
{

    nlohmann::json data;
    readSettings(data);

    mQmSettings.set(data);

    mQmProgram = data.value("qm program", mQmProgram);

    mRunner = shared_ptr<WaveFunctionDataGeneratorRunner>(WaveFunctionDataGeneratorRunner::create(mQmProgram));

    mRunner->set(data);

    mHardware.set(data);


    if(data.contains("max time minutes"))
        mMaxTimeMinutes = data["max time minutes"].get<int>();

    ChooseStructures::readChosenMol(mStructureName, mMolIdxInStructure);

    //    std::vector<std::string>&structureName,
    //    std::vector<int>&mol_idx)

    //input.open("chosen_mol");
    //vector<string> words;
    //string line;
    //do 
    //{
    //    getline(input, line);
    //    string_utilities::split(line, words);
    //    if (words.size() > 2)
    //    {
    //        mStructureName.push_back(words[0]);
    //        mMolIdxInStructure.push_back(stoi(words[1]));
    //    }
    //}
    //while (!words.empty());
    //input.close();

}

void WfnCalc::run()
{
    WallClockTimer timer;
    timer.start();
    WaveFunctionCalculationData wfnData;
    wfnData.hardware = mHardware;
    wfnData.qmSettings = mQmSettings;
    

    filesystem::path initialPath = filesystem::current_path();
    filesystem::path calcFolder = filesystem::current_path() / "wfnCalc";
    filesystem::create_directory(calcFolder);
    filesystem::current_path(calcFolder);
    vector<string> failedCalc;

    int nPrecalculated = 0;
    int nCalculated = 0;
    bool timebreak = false;
    for (int i = 0; i < mStructureName.size(); i++)
    {
        if (filesystem::exists(filesystem::path("wfnCalc/" + mStructureName[i] + ".wfx")))
        {
            nPrecalculated++;
            continue;
        }

        string mol2FileName = (filesystem::path("../chosen_mol_dir") / (mStructureName[i] + ".mol2")).string();
        mol2_io::Mol2Data data;
        mol2_io::read(mol2FileName, data);


        vector<int> atomicNumbers, z;
        vector<Vector3d> positions;
        int charge = 0;

        mol2_io::atomicNumbers(data, z);
        wfnData.qmSystem.atomicNumbers.clear();
        wfnData.qmSystem.positions.clear();

        for(int j=0;j<data.atomId.size(); j++)
            if (data.substructureIdx[j] == mMolIdxInStructure[i])
            {
                wfnData.qmSystem.atomicNumbers.push_back(z[j]);
                wfnData.qmSystem.positions.push_back(data.atomPosition[j]);
                
                charge += data.atomicCharge[j];
            }
       
        structural_properties::normalizeXhBondLengths(wfnData.qmSystem.positions, wfnData.qmSystem.atomicNumbers);

        wfnData.qmSystem.charge = charge;
        wfnData.qmSystem.spin_multilicity = 1;
        string fileName = mStructureName[i] + "_" + to_string(mMolIdxInStructure[i]) + ".inp";
        wfnData.jobName = mStructureName[i] + "_" + to_string(mMolIdxInStructure[i]);
        if (!mRunner->run(wfnData))
            failedCalc.push_back(wfnData.jobName);
        else
            nCalculated++;

        if (timer.elapsedTime() > mMaxTimeMinutes * 60000)
        {
            timebreak = true;
            break;
        }
    }
    filesystem::current_path(initialPath);

    ofstream out("wfn_calc.log");
    if (mStructureName.empty())
        out << "no systems for wfn calculations\n";
    else
    {
        if(timebreak)
            out << "time limit reached, calculations terminated\n";

        out << "run calculations for " << mStructureName.size() << " systems\n";
        out << nPrecalculated << " systems already had wfn file\n";
        out << nCalculated << " systems calculated\n";
        out << failedCalc.size() << " systems failed\n";
        out << mStructureName.size() - nPrecalculated - nCalculated - failedCalc.size() << " remaining systems not calculated\n";

        if (failedCalc.empty())
            out << "all calculations terminated normally\n";
        else
        {
            out << failedCalc.size() << " calculations not finished succesfully:\n";
            for (string& s : failedCalc)
                out << s << "\n";
        }
    }
    out.close();

}
