#include "discamb/QuantumChemistry/OrcaRunner.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/Sheduler.h"
#include "discamb/BasicUtilities/Task.h"
#include "discamb/BasicUtilities/discamb_env.h"
#include "discamb/BasicChemistry/periodic_table.h"

#include <iomanip>

#include <sstream>
#include <fstream>
#include <filesystem>

using namespace std;

namespace fs = std::filesystem; 

namespace {

    struct OrcaTask : public discamb::Task
    {
        discamb::OrcaRunner runner;
        discamb::WaveFunctionCalculationData data;
        string orcaFolder;
        virtual void run();
    };

    void OrcaTask::run()
    {
        //runner.setHardware(data.hardware.nCores, data.hardware.totalMemoryMB);
        string inputFileName = name + ".inp";
        
        runner.printInputFile(inputFileName, name + ".pc", data);
        //runner.printInputFile(name + ".inp", "", data.qmSystem.positions, data.qmSystem.atomicNumbers, 
          //  data.qmSystem.spin_multilicity, data.qmSystem.charge, vector<double>(0),
            //vector<discamb::Vector3d>(0), map<int, string>(), true);

        string outputFileName = name + ".log";

        string command = orcaFolder + "\\orca " + inputFileName + " > " + outputFileName;
        system(command.c_str());

        if (!runner.succesfulRun(name))
        {
            string errorMessage = "unsuccesful ORCA run, see file '" + outputFileName + "'";
            discamb::on_error::throwException(errorMessage, __FILE__, __LINE__);
        }
        
        system((orcaFolder + "\\orca_2aim " + name).c_str());
        string wfxFile = name + ".wfx";
        if (filesystem::exists(name + ".wfx"))
        {
            if (data.wfnFileName != wfxFile)
            {
                cout << data.wfnFileName << " " << wfxFile << "\n";
                command = "copy " + wfxFile + " " + data.wfnFileName;
                system(command.c_str());
            }
        }
        else
            discamb::on_error::throwException(data.wfnFileName + " was not created", __FILE__, __LINE__);
       
        runner.ecpWfxPostPorcessing(name, wfxFile);

        return;

        //--

        //string script_file = "run_" + name + ".bat";
        //ofstream out(script_file.c_str());
        //string executionFolder = "orca_dir\\" + name;
        //out << "mkdir " + executionFolder + "\n"
        //    << "copy " + inputFileName + " " + executionFolder + "\n"
        //    << "cd " + executionFolder + "\n"
        //    << orcaFolder + "\\orca " + inputFileName + " > " + name + ".log\n"
        //    << orcaFolder + "\\orca_2aim " + name + "\n"
        //    << "copy " << name << ".wfx " << filesystem::current_path().string() + "\\" << data.wfnFileName << "\n";
        //out.close();
        //system(script_file.c_str());
    }


}

namespace discamb {


    OrcaRunner::OrcaRunner()
    {
        setDefaults();
    }

    void OrcaRunner::setDefaults()
    {
        mExecFullPath = "orca";
        nPreviuslyUsedCharges = 0;
        mDefaultFileFormat = WaveFunctionFileFormat::wfn;
        mTryToReadGuess = true;
        //mTemplateHasCharges = false;
        tryToSetPathToOrcaFromSettingsFile();
        tryToSetPathToMolden2AimFromSettingsFile();
    }

    OrcaRunner::OrcaRunner(
        const nlohmann::json& settings)
    {
        set(settings);
    }

    OrcaRunner::~OrcaRunner()
    {
    }


    void OrcaRunner::set(
        const nlohmann::json& settings)
    {
        string exec, folder;       

        if (settings.find("wfn format") != settings.end())
            mDefaultFileFormat = WaveFunctionDataGeneratorRunner::wfnFileFormat(settings.find("wfn format").value().get<string>());

        if (settings.find("try reuse") != settings.end())
            mTryToReadGuess = settings.find("try reuse").value().get<bool>();


        if (settings.find("exec path") != settings.end())
        {
            mExecFullPath = settings.find("exec path").value().get<string>();
            mExecFolder = filesystem::path(mExecFullPath).remove_filename().string();
        }
        else
        {
            if (settings.find("folder") != settings.end() || settings.find("qm folder") != settings.end())
            {
                folder = settings.value("folder", folder);
                folder = settings.value("qm folder", folder);
                setExecFolder(folder);
            }
            else
                tryToSetPathToOrcaFromSettingsFile();
        }
        
        if (settings.find("molden2aim folder") != settings.end())
            setMolden2AimFolder(settings.find("molden2aim folder").value().get<string>());
        else
            tryToSetPathToMolden2AimFromSettingsFile();

    }

    void OrcaRunner::setExecFolder(
        const std::string& execFolder)
    {
        mExecFullPath = (fs::path(execFolder) / fs::path("orca")).string();
        mExecFolder = execFolder;
    }

    std::string OrcaRunner::getExecFolder()
        const
    {
        return mExecFolder;
    }

    void OrcaRunner::tryToSetPathToOrcaFromSettingsFile()
    {
        nlohmann::json defaults = discamb_env::get_discamb_dir_json("settings/settings.json");
        if (defaults.is_object())
        {
            if (defaults.find("orca folder") != defaults.end())
                setExecFolder(defaults.find("orca folder").value().get<string>());
                //mExecFullPath = (fs::path(defaults.find("orca folder").value().get<string>()) / fs::path("orca")).string();
        }
    }

    void OrcaRunner::tryToSetPathToMolden2AimFromSettingsFile()
    {
        nlohmann::json defaults = discamb_env::get_discamb_dir_json("settings/settings.json");
        if (defaults.is_object())
        {
            if (defaults.find("molden2aim folder") != defaults.end())
                setMolden2AimFolder(defaults["molden2aim folder"].get<string>());
        }
    }
    
    void OrcaRunner::ecpWfxPostPorcessing(
        const std::string& jobName,
        const std::string& wfxFile)
    {
        printMolden2AimIniFile();
        vector<string> elementsWithEcp;
        vector<int> nEcpElectrons;

        string orcaOutput = jobName + ".log";

        if (!checkForEcp(jobName + ".log", elementsWithEcp, nEcpElectrons))
            return;
        // make molden file 
        string command = mExecFolder + "\\orca_2mkl " + jobName + " -emolden";
        system(command.c_str());
        string moldenFileNameLong = jobName + ".molden.input";
        string moldenFileName = jobName + ".molden";
        filesystem::rename(moldenFileNameLong, moldenFileName);
        command = mMolden2AimFolder + "\\molden2aim.exe " + " -i " + moldenFileName;
        system(command.c_str());
    }

    void OrcaRunner::setMolden2AimFolder(
        const std::string& folder)
    {
        mMolden2AimFolder = folder;
        string execPath = getMolden2AimExecPath();
    }

    void OrcaRunner::printMolden2AimIniFile()
    {
        //if (filesystem::exists(filesystem::path("m2a.ini")))
          //  return;

        vector<string> lines  
        {   "molden= -1", "wfn= -1", "wfncheck= -1",
            "wfx= 1", "wfxcheck= -1", "nbo= -1",
            "nbocheck= -1", "wbo= -1", "program=1",
            "edftyp=1", "allmo = 0", "unknown = 1",
            "nosupp=1", "clear=1"};

        ofstream out("m2a.ini");
        if (!out.good())
            on_error::throwException("it was not possible to write to m2a.ini file needed when processing ORCA output for calculations with ECP", __FILE__, __LINE__);

        for (auto& line : lines)
            out << line << "\n";
        out.close();
    }


    bool OrcaRunner::checkForEcp(
        const std::string& orcaOutput,
        std::vector<std::string>& elementsWithEcp,
        std::vector<int>& nEcpElectrons)
    {
        elementsWithEcp.clear();
        nEcpElectrons.clear();

        ifstream in(orcaOutput);
        if (!in.good())
            on_error::throwException("cannot read ORCA output file '" + orcaOutput + "'", __FILE__, __LINE__);
        string line;
        vector<string> words;
        bool foundWhatExpected;
        while (in.good())
        {
            getline(in, line);
            if (line.find("#ECP")!=string::npos)
            {
                /*
                NewECP Sb
                N_core 28
                */

                foundWhatExpected = false;
                getline(in, line);
                string_utilities::split(line, words);

                if (words.size() == 2)
                    if (words[0] == string("NewECP"))
                    {
                        foundWhatExpected = true;
                        elementsWithEcp.push_back(words[1]);
                    }
                if (!foundWhatExpected)
                    on_error::throwException("problem with processing ECP information from ORCA output, unexpected format", __FILE__, __LINE__);

                foundWhatExpected = false;
                getline(in, line);
                string_utilities::split(line, words);
                if (words.size() == 2)
                    if (words[0] == string("N_core"))
                    {
                        foundWhatExpected = true;
                        nEcpElectrons.push_back(stoi(words[1]));
                    }
                if (!foundWhatExpected)
                    on_error::throwException("problem with processing ECP information from ORCA output, unexpected format", __FILE__, __LINE__);

            }
        }
        in.close();
        return (!nEcpElectrons.empty());
    }

    std::string OrcaRunner::getMolden2AimExecPath(
        bool required)
        const
    {
        filesystem::path molden2aimFolderPath(mMolden2AimFolder);
        filesystem::path molden2aimPath = molden2aimFolderPath / "molden2aim.exe";

        if (!filesystem::exists(molden2aimPath))
        {
            if (required)
                on_error::throwException("cannot find molden2aim at the path provided:\n"+ molden2aimPath.string(), __FILE__, __LINE__);
            else
                return string();
        }

        return molden2aimPath.string();
    }

    void OrcaRunner::addEcpIndicatorToMoldenFile(
        const string& moldenFileName,
        const string& outputName,
        vector<string>& elementsWithEcp,
        vector<int>& nEcpElectrons)
    {
        ifstream in(moldenFileName);
        ofstream out(outputName);
        string line;
        
        while (in.good())
        {
            getline(in, line);
            if (line == string("[GTO]"))
            {
                out << "[Core]\n";
                for (int i = 0; i < elementsWithEcp.size(); i++)
                    out << elementsWithEcp[i] << ": " << nEcpElectrons[i] << "\n";
            }
            out << line << "\n";
        }
        in.close();
        out.close();
    }

    bool OrcaRunner::succesfulRun(
        const std::string& jobName)
    {
        ifstream in(jobName + ".log");
        if (!in.good())
            return false;
        string line;
        bool terminatedNormally = false;
        while (in.good())
        {
            getline(in, line);
            if (line.find("****ORCA TERMINATED NORMALLY****") != string::npos)
                terminatedNormally = true;
        }

        if (!terminatedNormally)
            return false;
        
        //if(!filesystem::exists(filesystem::current_path()/(jobName+".wfx")))
        //    return false;

        return true;
    }

    bool OrcaRunner::run(const WaveFunctionCalculationData& inputData)
        const
    {

        //system("mkdir orca_dir");

        OrcaTask orcaTask;
        orcaTask.orcaFolder = mExecFolder;

        orcaTask.data = inputData;
        orcaTask.n = inputData.hardware.nCores;
        orcaTask.data.hardware.totalMemoryMB = inputData.hardware.totalMemoryMB;
        orcaTask.name = inputData.jobName;
        orcaTask.orcaFolder = mExecFolder;
        orcaTask.run();

        return succesfulRun(inputData.jobName);
    }

    void OrcaRunner::runMultipleJobs(
        const std::vector<WaveFunctionCalculationData>& inputData,
        int nCores,
        int totalMemory,
        std::vector<bool>& succesful)
        const
    {
        runMultipleJobs(inputData, nCores, totalMemory);
        succesful.clear();
        for (const auto& inp : inputData)
            succesful.push_back(succesfulRun(inp.jobName));
    }

    void OrcaRunner::runMultipleJobs(
        const std::vector<WaveFunctionCalculationData>& inputData, 
        int nCores,
        int memoryPerCore)
        const
    {
        if (nCores == 0)
            nCores = 1;
        if (memoryPerCore == 0)
            memoryPerCore = 200;
        
        system("mkdir orca_dir");
//        OrcaRunnerData orcaRunnerData;
        int fragmentIdx, nFragments = inputData.size();

        Scheduler scheduler;
        scheduler.nFreeCPU = nCores;
        OrcaTask orcaTask;
        orcaTask.orcaFolder = mExecFolder;
        for (fragmentIdx = 0; fragmentIdx < nFragments; fragmentIdx++)
        {
            //orcaRunnerData.atomicNumbers = atomicNumbers[fragmentIdx];
            //orcaRunnerData.charge = charge[fragmentIdx];
            //orcaRunnerData.spin_multilicity = spinMultiplicity[fragmentIdx];
            //orcaRunnerData.nCpu = bestCombination[fragmentIdx];
            //orcaRunnerData.memoryMB = 3000 * orcaRunnerData.nCpu;
            //orcaRunnerData.positions = fragments[fragmentIdx].atomPosition;
            orcaTask.data = inputData[fragmentIdx];
            orcaTask.n = inputData[fragmentIdx].hardware.nCores;
            orcaTask.data.hardware.totalMemoryMB = orcaTask.n * memoryPerCore;
            //orcaTask.name = "frag_" + to_string(fragmentIdx + 1);
            string wfnFile = inputData[fragmentIdx].wfnFileName;
            orcaTask.name = wfnFile.substr(0, wfnFile.find(".wfx"));
            orcaTask.orcaFolder = mExecFolder;
            orcaTask.runner = *this;
            scheduler.tasks.push_back({ orcaTask.n, orcaTask.name, true, shared_ptr<Task>(new OrcaTask(orcaTask)) });
        }
        scheduler.nToFinish = scheduler.tasks.size();
        scheduler.run();

    }

    void OrcaRunner::correctMoldenFile(
        const std::string& wfnFile,
        int nAtoms,
        int nCharges)
        const
    {
        ifstream in(wfnFile);
        stringstream ss;
        if (!in.good())
            on_error::throwException(string("cannot read file '") + wfnFile + string("'"), __FILE__, __LINE__);

        string line;
        bool checked = false;
        bool removedQ = false;
        bool number_and_0_section_to_remove;
        vector<string> words;
        char eol = char(10);
        while (getline(in, line))
        {
            if (checked)
                ss << line << eol;
            else
            {
                if (!line.empty())
                {
                    if (!removedQ)
                    {
                        if (line[0] == 'Q')
                        {
                            for (int i = 1; i < nCharges; i++)
                                getline(in, line);
                            removedQ = true;
                        }
                        else
                            ss << line << eol;
                    }
                    else
                    {
                        string_utilities::split(line, words);
                        number_and_0_section_to_remove = false;
                        if (words.size() == 2)
                        {
                            bool digits = true;
                            for (char c : words[0])
                                if (!isdigit(c))
                                    digits = false;
                            for (char c : words[1])
                                if (!isdigit(c))
                                    digits = false;
                            if (digits)
                                if(atoi(words[0].c_str())>nAtoms)
                                    number_and_0_section_to_remove = true;
                        }

                        if (number_and_0_section_to_remove)
                        {
                            for (int i = 1; i < nCharges; i++)
                            {
                                getline(in, line);
                                getline(in, line);
                            }
                            checked = true;
                        }
                        else
                            ss << line << eol;
                            
                    }
                }
                else
                    ss << eol;

            }
        }
        in.close();
        ofstream out(wfnFile, std::ios_base::binary | std::ios_base::out);
        out << ss.rdbuf();
        out.close();
    }

    void OrcaRunner::makeWfnLookGaussianMade(
        const std::string fileName)
        const
    {
        ifstream in(fileName);
        if (!in.good())
            on_error::throwException(string("cannot read file '")+fileName + string("'"), __FILE__, __LINE__);
        stringstream ss;
        ss << in.rdbuf();
        in.close();
        ofstream out;
        out.open(fileName);
        string line;
        getline(ss, line);
        out << line << "\n";
        getline(ss, line);
        line = string("GAUSSIAN") + line.substr(8);
        out << line << "\n";
        while(getline(ss, line))
        {
            out << line << "\n";
        }
        out.close();

    }

        /*  */
    bool OrcaRunner::supportsFormat(
        WaveFunctionFileFormat format)
        const
    {
        if (format == WaveFunctionFileFormat::wfn || format == WaveFunctionFileFormat::wfx || 
            format == WaveFunctionFileFormat::mkl || format == WaveFunctionFileFormat::molden)
            return true;
        return false;
    }

    void OrcaRunner::supportedFormats(
        std::vector<WaveFunctionFileFormat>& formats) 
        const
    {
        formats.clear();
        formats.push_back(WaveFunctionFileFormat::wfn);
        formats.push_back(WaveFunctionFileFormat::wfx);
        formats.push_back(WaveFunctionFileFormat::mkl);
        formats.push_back(WaveFunctionFileFormat::molden);
    }

    WaveFunctionFileFormat OrcaRunner::defaultFileFormat()
        const
    {
        return mDefaultFileFormat;
    }

    void OrcaRunner::printInputFileFromTemplate(
        const std::string fileName,
        const WaveFunctionCalculationData& inputData)
        const
    {
        map<string, string> dictionary;
        string method;
        string_utilities::toUpper(inputData.qmSettings.qmMethod, method);

        // geometry as string

        string geometryStr;
        stringstream geometryStreamString;
        for (int i = 0; i < inputData.qmSystem.positions.size(); i++)
        {

            geometryStreamString << periodic_table::symbol(inputData.qmSystem.atomicNumbers[i]) << " ";
            geometryStreamString << fixed << setw(14) << setprecision(6) << inputData.qmSystem.positions[i][0]
                << fixed << setw(14) << setprecision(6) << inputData.qmSystem.positions[i][1]
                << fixed << setw(14) << setprecision(6) << inputData.qmSystem.positions[i][2] << "\n";
        }
        geometryStr = geometryStreamString.str();
        geometryStr = geometryStr.substr(0, geometryStr.size() - 1);

        int nCores = (inputData.hardware.nCores == 0) ? 1 : inputData.hardware.nCores;
        int memoryPerCore = (inputData.hardware.totalMemoryMB == 0) ? 200 : inputData.hardware.totalMemoryMB / nCores;


        dictionary["n cpu"] = to_string(nCores);
        dictionary["memory"] = to_string(memoryPerCore);
        dictionary["method"] = method;
        dictionary["basis set"] = inputData.qmSettings.basisSet;
        dictionary["spin"] = to_string(inputData.qmSystem.spin_multilicity);
        dictionary["charge"] = to_string(inputData.qmSystem.charge);
        dictionary["geometry"] = geometryStr;
        //inputFileName
        string chargesFile = "pointcharges.pc";
        if (fileName.find_last_of('.') != string::npos)
            chargesFile = fileName.substr(0, fileName.find_last_of('.')) + ".pc";

        dictionary["point charge file"] = chargesFile;// "pointcharges.pc";

        //
        string fileContent;

        string_utilities::fill_template(dictionary, inputData.qmSettings.inputTemplate, fileContent, '$');

        ofstream out(fileName);
        if (!out.good())
            on_error::throwException(string("cannot print to Orca input file '") + fileName + string("'"), __FILE__, __LINE__);
        out << fileContent;
        out.close();

        if (!inputData.qmSystem.pointChargeValue.empty())
            printPointCharges(chargesFile, inputData.qmSystem.pointChargeValue, inputData.qmSystem.pointChargePosition);


    }


    //void OrcaRunner::printInputFileFromTemplate(
    //    const std::string inputFileName,
    //    const std::vector<Vector3d>& positions,
    //    const std::vector<int>& atomicNumbers,
    //    int spin_multilicity,
    //    int charge,
    //    const std::vector<double>& pointChargeValue,
    //    const std::vector<Vector3d>& pointChargePosition)
    //    const
    //{

    //    map<string, string> dictionary;
    //    string method;
    //    string_utilities::toUpper(mMethod, method);

    //    // geometry as string

    //    string geometryStr;
    //    stringstream geometryStreamString;
    //    for (int i = 0; i < positions.size(); i++)
    //    {

    //        geometryStreamString << periodic_table::symbol(atomicNumbers[i]) << " ";
    //        geometryStreamString << fixed << setw(14) << setprecision(6) << positions[i][0]
    //            << fixed << setw(14) << setprecision(6) << positions[i][1]
    //            << fixed << setw(14) << setprecision(6) << positions[i][2] << "\n";
    //    }
    //    geometryStr = geometryStreamString.str();
    //    geometryStr = geometryStr.substr(0, geometryStr.size() - 1);
    //    // point charges as string

    //    //stringstream chargesStringStream;
    //    //string chargesAsString;
    //    //if (!pointChargeValue.empty())
    //    //{
    //    //    for (int i = 0; i < pointChargeValue.size(); i++)
    //    //        chargesStringStream << "Q "
    //    //        << fixed << setw(14) << setprecision(6) << pointChargeValue[i]
    //    //        << fixed << setw(14) << setprecision(6) << pointChargePosition[i][0]
    //    //        << fixed << setw(14) << setprecision(6) << pointChargePosition[i][1]
    //    //        << fixed << setw(14) << setprecision(6) << pointChargePosition[i][2] << "\n";
    //    //    chargesAsString = chargesStringStream.str();
    //    //}

    //    //

    //    //

    //    dictionary["n cpu"] = to_string(mN_Core);
    //    dictionary["memory"] = mMemoryPerCore;
    //    dictionary["method"] = method;
    //    dictionary["basis set"] = mBasisSet;
    //    dictionary["spin"] = to_string(spin_multilicity);
    //    dictionary["charge"] = to_string(charge);
    //    dictionary["geometry"] = geometryStr;
    //    //inputFileName
    //    string chargesFile = "pointcharges.pc";
    //    if (inputFileName.find_last_of('.') != string::npos)
    //        chargesFile = inputFileName.substr(0, inputFileName.find_last_of('.')) + ".pc";

    //    dictionary["point charge file"] = chargesFile;// "pointcharges.pc";

    //    //
    //    string fileContent;

    //    string_utilities::fill_template(dictionary, mInputTemplate, fileContent, '$');

    //    ofstream out(inputFileName);
    //    if (!out.good())
    //        on_error::throwException(string("cannot print to Orca input file '") + inputFileName + string("'"), __FILE__, __LINE__);
    //    out << fileContent;
    //    out.close();

    //    //printPointCharges("pointcharges.pc", pointChargeValue, pointChargePosition);
    //    if(!pointChargeValue.empty())
    //        printPointCharges(chargesFile, pointChargeValue, pointChargePosition);

    //}

    void OrcaRunner::printPointCharges(
        const std::string& fName,
        const std::vector<double>& pointChargeValue,
        const std::vector<Vector3d>& pointChargePosition)
    {
        ofstream out(fName);
        if (!out.good())
            on_error::throwException("cannot write to '" + fName + "'", __FILE__, __LINE__);
        out << pointChargeValue.size() << "\n";
        for (int i = 0; i < pointChargeValue.size(); i++)
        {
            out << setprecision(6) << fixed << pointChargeValue[i];
            for (int j = 0; j < 3; j++)
                out << " " << setprecision(6) << fixed << pointChargePosition[i][j];
            out << "\n";
        }
        out.close();

    }


    //void OrcaRunner::printInputFile(
    //    const std::string fileName,
    //    const std::string pointChargeFileName,
    //    const std::vector<Vector3d>& positions,
    //    const std::vector<int>& atomicNumbers,
    //    int spin_multilicity,
    //    int charge,
    //    const std::vector<double>& pointChargeValue,
    //    const std::vector<Vector3d>& pointChargePosition,
    //    const std::map<int, std::string>& atomIdx2BasisSetMap,
    //    bool printAimCommand)
    //    const
    //{

    //    if (!mInputTemplate.empty())
    //    {
    //        printInputFileFromTemplate(fileName, positions, atomicNumbers, spin_multilicity,
    //            charge, pointChargeValue, pointChargePosition);
    //        return;
    //    }


    //    ofstream out(fileName);
    //    if (!out.good())
    //        on_error::throwException(string("cannot open file '") + fileName + string("' for writing"), __FILE__, __LINE__);

    //    out << "%pal nprocs " << mN_Core << "\n"
    //        << "end\n\n";


    //    if(mMemoryPerCore != 0)
    //        out << "%maxcore " << mMemoryPerCore << "\n\n";

    //    out << "! " << mMethod;
    //    if (!mRelativisticMethod.empty())
    //        out << " " << mRelativisticMethod;
    //    if (printAimCommand)
    //        out << " AIM";

    //    out << " " << mBasisSet << " PrintBasis\n\n";


    //    if(!pointChargeValue.empty())
    //        out<< "%pointcharges \""+pointChargeFileName+"\"\n";
    //    
    //    // element specific basis set

    //    if (!mAtomicNumber2BasisSetMap.empty())
    //    {
    //        out << "%basis\n";
    //        for (auto& item : mAtomicNumber2BasisSetMap)
    //            out << "  NewGTO " << item.first << "\n"
    //                << "    \"" << item.second << "\"\n"
    //                << "  end\n";
    //        out << "end\n";
    //    }
    //    
    //    // coordinates and atom sepcific basis set

    //    //if (basisSetPerAtom.empty())
    //    //{
    //    //    out << "* xyz " << charge << " " << spin_multilicity << "\n";

    //    //    for (int i = 0; i < positions.size(); i++)
    //    //    {
    //    //        out << periodic_table::symbol(atomicNumbers[i]);
    //    //        for (int j = 0; j < 3; j++)
    //    //            out << " " << setprecision(6) << fixed << positions[i][j];
    //    //        out << "\n";
    //    //    }
    //    //    out << "*\n";
    //    //}
    //        /*
    //        %coords
    //        CTyp xyz # the type of coordinates xyz or internal
    //        Charge -2 # the total charge of the molecule
    //        Mult 2 # the multiplicity = 2S+1
    //        coords
    //        Cu(1) 0 0 0
    //        Cl(2) 2.25 0 0
    //        Cl(2) -2.25 0 0
    //        Cl(2) 0 2.25 0
    //        Cl(2) 0 -2.25 0
    //        end
    //        end
    //        */
    //    out << "%coords\n"
    //        << " CTyp xyz\n"
    //        << " Charge " << charge << "\n"
    //        << " Mult " << spin_multilicity << "\n"
    //        << " coords\n";

    //    for (int i = 0; i < positions.size(); i++)
    //    {
    //        out << "    " << periodic_table::symbol(atomicNumbers[i]);
    //        for (int j = 0; j < 3; j++)
    //            out << " " << setprecision(6) << fixed << positions[i][j];
    //        out << "\n";
    //        auto atomSpecificBasisSetIterator = atomIdx2BasisSetMap.find(i);
    //        if (atomSpecificBasisSetIterator!=atomIdx2BasisSetMap.end())
    //            out << "     NewGTO\n"
    //                << "       \"" << atomSpecificBasisSetIterator->second << "\"\n"
    //                << "     end\n";
    //    }
    //    out << " end\n"
    //        << "end\n";
    //    out.close();

    //    if (!pointChargeValue.empty())
    //        printPointCharges(pointChargeFileName, pointChargeValue, pointChargePosition);
    //    //{
    //    //    
    //    //    ////------ this was in main input file
    //    //    //for (int i = 0; i < pointChargeValue.size(); i++)
    //    //    //{
    //    //    //    out << "Q " << setprecision(6) << fixed << pointChargeValue[i];
    //    //    //    for (int j = 0; j < 3; j++)
    //    //    //        out << " " << setprecision(6) << fixed << pointChargePosition[i][j];
    //    //    //    out << "\n";
    //    //    //}
    //    //    //-------but external point charges input has different format
    //    //    out.open("pointcharges.pc");
    //    //    out << pointChargeValue.size() << "\n";
    //    //    for (int i = 0; i < pointChargeValue.size(); i++)
    //    //    {
    //    //        out << setprecision(6) << fixed << pointChargeValue[i];
    //    //        for (int j = 0; j < 3; j++)
    //    //            out << " " << setprecision(6) << fixed << pointChargePosition[i][j];
    //    //        out << "\n";
    //    //    }
    //    //    out.close();
    //    //}
    //}
    
//private:
//        std::string mExecFullPath;
//        std::string mMethod;
//        int mN_Core;
//        int mMemoryPerCore;
//        std::string mBasisSet;
//        WaveFunctionFileFormat mDefaultFileFormat;
//        bool mTryToReadGuess;
//        void setDefaults();
//        std::string mInputTemplate;
//        bool mTemplateHasCharges;
//
//        void printInputFileFromTemplate(
//            const std::string fileName,
//            const std::vector<Vector3d>& positions,
//            const std::vector<int>& atomicNumbers,
//            int spin_multilicity,
//            int charge,
//            const std::string& wfnFileName,
//            const std::vector<double>& pointChargeValue,
//            const std::vector<Vector3d>& pointChargePosition,
//            WaveFunctionFileFormat format,
//            const std::string& jobName) const;
//
//
//    };

    /*
    void OrcaRunner::printInputFile(
        const std::string fileName,
        const std::string pointChargeFileName,
        const std::vector<Vector3d>& positions,
        const std::vector<int>& atomicNumbers,
        int spin_multilicity,
        int charge,
        const std::vector<double>& pointChargeValue,
        const std::vector<Vector3d>& pointChargePosition,
        const std::map<int, std::string>& atomIdx2BasisSetMap,
        bool printAimCommand)
        const

    */

    void OrcaRunner::printInputFile(
    const std::string fileName,
    const WaveFunctionCalculationData& inputData)
    const
    {
        printInputFile(fileName, fileName + ".pc", inputData);
    }


    void OrcaRunner::printInputFile(
        const std::string fileName,
        const std::string pointChargeFileName,
        const WaveFunctionCalculationData& inputData)
        const
    {
        auto const &qmSettings = inputData.qmSettings;
        auto const &qmSystem = inputData.qmSystem;

        if (! inputData.qmSettings.inputTemplate.empty())
        {
            printInputFileFromTemplate(fileName, inputData);
            //printInputFileFromTemplate(fileName, positions, atomicNumbers, spin_multilicity,
            //    charge, pointChargeValue, pointChargePosition);
            return;
        }


        ofstream out(fileName);
        if (!out.good())
            on_error::throwException(string("cannot open file '") + fileName + string("' for writing"), __FILE__, __LINE__);

        int nCores = (inputData.hardware.nCores == 0) ? 1 : inputData.hardware.nCores;

        out << "%pal nprocs " << nCores << "\n"
            << "end\n\n";

        int memoryPerCore = (inputData.hardware.totalMemoryMB == 0) ? 200 : inputData.hardware.totalMemoryMB / nCores;

        out << "%maxcore " << memoryPerCore << "\n\n";

        out << "! " << qmSettings.qmMethod;
        if (!qmSettings.relativisticMethod.empty())
            out << " " << qmSettings.relativisticMethod;
        //if (printAimCommand)
        //    out << " AIM";

        out << " " << qmSettings.basisSet << " PrintBasis\n\n";


        if (!qmSystem.pointChargeValue.empty())
            out << "%pointcharges \"" + pointChargeFileName + "\"\n";

        // element specific basis set

        if (! qmSettings.atomicNumber2BasisSetMap.empty())
        {
            out << "%basis\n";
            for (auto& item : qmSettings.atomicNumber2BasisSetMap)
                out << "  NewGTO " << item.first << "\n"
                << "    \"" << item.second << "\"\n"
                << "  end\n";
            out << "end\n";
        }

        out << "%coords\n"
            << " CTyp xyz\n"
            << " Charge " << qmSystem.charge << "\n"
            << " Mult " << qmSystem.spin_multilicity << "\n"
            << " coords\n";

        for (int i = 0; i < qmSystem.positions.size(); i++)
        {
            out << "    " << periodic_table::symbol(qmSystem.atomicNumbers[i]);
            for (int j = 0; j < 3; j++)
                out << " " << setprecision(6) << fixed << qmSystem.positions[i][j];
            out << "\n";
            auto atomSpecificBasisSetIterator = inputData.atomIdx2BasisSetMap.find(i);
            if (atomSpecificBasisSetIterator != inputData.atomIdx2BasisSetMap.end())
                out << "     NewGTO\n"
                << "       \"" << atomSpecificBasisSetIterator->second << "\"\n"
                << "     end\n";
        }
        out << " end\n"
            << "end\n";
        out.close();

        if (!qmSystem.pointChargeValue.empty())
            printPointCharges(pointChargeFileName, qmSystem.pointChargeValue, qmSystem.pointChargePosition);

    }


}
