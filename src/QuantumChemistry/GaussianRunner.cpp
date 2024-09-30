#include "discamb/QuantumChemistry/GaussianRunner.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/discamb_env.h"
#include "discamb/config.h"



#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

namespace discamb {

    void GaussianRunner::setDefaults()
    {
        //mDefaultFileFormat = WaveFunctionFileFormat::fchk;
        mDefaultFileFormat = WaveFunctionFileFormat::wfn;
        mChangeDir = false;
        //mMemory = "";
        //mN_Core = 1;
        mGenerateOutput = true;
        
        mScfOptions = "Conver=8";
        //mTemplateHasCharges = false;
        //setExec(string(), string());
    }

    void GaussianRunner::setExecFolder(const std::string& execFolder)
    {
        string execName = findGaussianExecName(execFolder);
        setExec(execFolder, execName);
    }

    std::string GaussianRunner::getExecFolder()
        const
    {
        return mExecFolder;
    }


    std::string GaussianRunner::findGaussianExecName(
        const std::string& folder)
    {
        vector<string> possibleNames = { "g98","g98.exe","g03","g03.exe","g09","g09.exe","g16","g16.exe" };

        for (auto execName : possibleNames)
        {
            filesystem::path p = filesystem::path(folder) / execName;
            if (filesystem::exists(p))
                return execName;
        }
        return string();
    }


    //void GaussianRunner::setHardware(
    //    int nCore,
    //    int totalMemoryMB)
    //{
    //    mN_Core = nCore;
    //    if (totalMemoryMB % 1000 == 0)
    //        mMemory = to_string(totalMemoryMB / 1000) + "GB";
    //    else
    //        mMemory = to_string(totalMemoryMB) + "MB";
    //}

    std::string GaussianRunner::memoryAsString(
        int memoryMB)
    {
        if (memoryMB % 1000 == 0)
            return (to_string(memoryMB / 1000) + "GB");
        
        return (to_string(memoryMB) + "MB");

    }

	GaussianRunner::GaussianRunner()
	{
        setDefaults();
	}

	GaussianRunner::GaussianRunner(
		const nlohmann::json& settings) 
	{
        setDefaults();
		set(settings);
	}

	GaussianRunner::GaussianRunner(
		const std::string& executionFolder,
		const std::string& execName)
	{
        mDefaultFileFormat = WaveFunctionFileFormat::fchk;
        mChangeDir = false;
        //mMemory = "2GB";
        //mN_Core = 1;

		setExec(executionFolder, execName);
	}

	GaussianRunner::~GaussianRunner() 
	{
	}

	// void GaussianRunner::setTheoryLevel(
	//	const std::string& method, 
	//	const std::string& basis_set,
 //       const std::map<int, std::string>& atomicNumber2BasisSetMap,
 //        const std::string& relativisticMethod)
	//{
	//	mMethod = method;
	//	mBasisSet = basis_set;
 //       mAtomicNumber2BasisSetMap = atomicNumber2BasisSetMap;
 //       mRelativisticMethod = relativisticMethod;
	//}

    void GaussianRunner::setScfOptions(const std::string& scfOptions)
    {
        mScfOptions = scfOptions;
    }


	void GaussianRunner::set(const nlohmann::json& settings) 
	{
	 	string exec, folder;

        //if (settings.find("template") != settings.end())
        //    setInputTemplate(settings.find("template").value().get<string>());

   //     if (settings.find("memory") != settings.end())
			//mMemory = settings.find("memory").value().get<string>();

   //     if (settings.find("n cores") != settings.end())
   //         mN_Core = settings.find("n cores").value().get<int>();

        //if (settings.find("wfn format") != settings.end())
        //    mDefaultFileFormat = WaveFunctionDataGeneratorRunner::wfnFileFormat(settings.find("wfn format").value().get<string>());

        //if (settings.find("try reuse") != settings.end())
        //    mTryToReadGuess = settings.find("try reuse").value().get<bool>();

        nlohmann::json settingsJson = discamb_env::get_discamb_dir_json("settings/settings.json");

        if (settings.find("exec") != settings.end())
            exec = settings.find("exec").value().get<string>();
        else
            exec = mExecName;


        if (settings.find("scf options") != settings.end())
            setScfOptions( settings.find("scf options").value().get<string>() );
 
//        "gaussian folder": "c:\\programy\\G16W",
  //          "gaussian exec" : "g16.exe",

        

        if (settings.find("folder") != settings.end())
        {
            folder = settings.find("folder").value().get<string>();
//            mChangeDir = true;
        }
        else
        {
            folder = mExecFolder;
//            mChangeDir = false;
        }

		setExec(folder, exec);
	}

    void GaussianRunner::getExec(
        std::string& executionFolder,
        std::string& execName)
        const
    {
        executionFolder = mExecFolder;
        execName = mExecName;

    }

	void GaussianRunner::setExec(
		const std::string& execFolder,
		const std::string& execName)
        const
	{
        mExecFolder = execFolder;
		mExecName = execName;

        if (!mExecName.empty() && !mExecFolder.empty())
        {
            checkPaths();
            return;
        }

        nlohmann::json settingsJson = discamb_env::get_discamb_dir_json("settings/settings.json");

        if (settingsJson.is_object())
        {

            if (mExecName.empty())
                if (settingsJson.find("gaussian exec") != settingsJson.end())
                    mExecName = settingsJson["gaussian exec"].get<string>();

            if (mExecFolder.empty())
                if (settingsJson.find("gaussian folder") != settingsJson.end())
                    mExecFolder = settingsJson["gaussian folder"].get<string>();
        }


#ifdef  CMAKE_DETECTED_WIN32
        
        if (mExecName.empty())
        {
            string candidate;
            bool execFound=false;
            filesystem::path gaussianFolderPath(mExecFolder);
            for(auto &entry: filesystem::directory_iterator(gaussianFolderPath))
                if (entry.is_regular_file())
                {
                    candidate = entry.path().filename().string();
                    if (candidate[0] == 'g' && candidate.find(".exe") != string::npos)
                        if (candidate.find("w.exe") == string::npos)
                        {
                            for (int i = 1; i < candidate.size(); i++)
                                if (isdigit(candidate[i]) || candidate[i] == '.')
                                {
                                    if (candidate[i] == '.')
                                    {
                                        mExecName = candidate;
                                        execFound = true;
                                    }
                                }
                                else
                                    i = 10000;
                        }
                }
        }
#endif
        checkPaths();
	}

    void GaussianRunner::checkPaths()
        const
    {
        filesystem::path gaussianExecFullPath(mExecFolder);
        gaussianExecFullPath /= mExecName;

        if (!filesystem::exists(gaussianExecFullPath) || !filesystem::is_regular_file(gaussianExecFullPath))
            on_error::throwException("cannot find Gaussian exec file, tested full path: '" + 
                                      gaussianExecFullPath.string()+ "'", __FILE__, __LINE__);

    }

    void GaussianRunner::runMultipleJobs(
        const std::vector<WaveFunctionCalculationData>& inputData,
        int nCores,
        int memoryPerCore)
        const
    {
        bool overwriteHardware = (nCores == 0);
        HardwareResources hardware;
        hardware.nCores = nCores;
        hardware.totalMemoryMB = memoryPerCore * nCores;
        WaveFunctionCalculationData data;

        for (auto const& inpData : inputData)
        {
            data = inpData;
            if(overwriteHardware)
                data.hardware = hardware;
            run(data);
        }
    }

    void GaussianRunner::runMultipleJobs(
        const std::vector<WaveFunctionCalculationData>& inputData,
        int nCores, int totalMemory,
        std::vector<bool>& succesful)
        const
    {
        runMultipleJobs(inputData, nCores, totalMemory);
        succesful.clear();
        for (const auto& inp : inputData)
            succesful.push_back(succesfulRun(inp.jobName, inp.wfnFileName));

    }

    bool GaussianRunner::run(
        const WaveFunctionCalculationData& inputData)
        const
    {
        auto const& qmSettings = inputData.qmSettings;
        auto const& qmSystem = inputData.qmSystem;
        if (inputData.qmSettings.qmMethod.empty() && inputData.qmSettings.inputTemplate.empty())
            on_error::throwException("wave function calculation method not defined", __FILE__, __LINE__);
        if (qmSettings.basisSet.empty() && qmSettings.inputTemplate.empty())
            if(qmSettings.qmMethod != "PM6")
                on_error::throwException("basis set for wave function calculation not defined", __FILE__, __LINE__);


        if (mExecName.empty())
        {
            setExec(mExecFolder, mExecName);
            if (mExecName.empty())
                on_error::throwException("Gaussian exec name (e.g. 'gw16') is not specified", __FILE__, __LINE__);
        }
        //string wfnFileName = inputData.jobName + ".wfx";
        if (fs::exists(inputData.wfnFileName))
            fs::remove(inputData.wfnFileName);

        fs::path originalCurrentPath = fs::current_path();

        string inputFileName = inputData.jobName + string(".gjf");
        printInputFile(inputFileName, inputData);



        string outputName = inputData.jobName + string(".out");

        //

        string command;

        fs::path p = fs::path(mExecFolder) / fs::path(mExecName);

#ifdef CMAKE_DETECTED_WIN32
        command = "set GAUSS_EXEDIR=" + mExecFolder + "& " + p.string() + " " + inputFileName; 
#else
        command = mExecName + string(" ") + inputFileName;
#endif

        //string command = mExecName + string(" ") + inputFileName;
        if (mGenerateOutput)
            command += string(" ") + outputName;
        cout << "calling " + command << "\n";
        system(command.c_str());

        if (inputData.wfnFormat == WaveFunctionFileFormat::fchk)
        {
            command = string("formchk task.chk ") + inputData.wfnFileName;
            cout << "calling : '" << command << "'\n";
            system(command.c_str());

        }

        if (!fs::exists(inputData.wfnFileName))
            on_error::throwException(string("wave function file '") + inputData.wfnFileName + string("' was not created"), __FILE__, __LINE__);
        
        if (mChangeDir)
            if (fs::exists(originalCurrentPath / inputData.wfnFileName))
                fs::remove(originalCurrentPath / inputData.wfnFileName);

        if (mChangeDir)
        {
            fs::copy_file(inputData.wfnFileName, originalCurrentPath / inputData.wfnFileName);
            fs::current_path(originalCurrentPath);
        }

        return succesfulRun(inputData.jobName, inputData.wfnFileName);
    }

    bool GaussianRunner::succesfulRun(
        const std::string& jobName,
        const std::string& wfnName)
    {
        // Normal termination of Gaussian
        ifstream in(jobName + ".log");
        if (!in.good())
            return false;
        string line;
        bool terminatedNormally = false;
        while (in.good())
        {
            getline(in, line);
            if (line.find("Normal termination of Gaussian") != string::npos)
                terminatedNormally = true;
        }

        if (!terminatedNormally)
            return false;

        if (!filesystem::exists(wfnName))
            return false;

        return true;
    }

    void GaussianRunner::printInputFileFromTemplate(
        const std::string fileName,
        const WaveFunctionCalculationData& inputData)
        const
    {
        auto const& qmSettings = inputData.qmSettings;
        auto const& system = inputData.qmSystem;
        map<string, string> dictionary;

        // geometry as string

        string geometryStr;
        stringstream geometryStreamString;
        geometryStreamString << fixed << setprecision(6);
        for (int i = 0; i < system.positions.size(); i++)
        {

            geometryStreamString << periodic_table::symbol(system.atomicNumbers[i]) << " ";
            geometryStreamString << setw(14) << system.positions[i][0]
                                 << setw(14) << system.positions[i][1]
                                 << setw(14) << system.positions[i][2] << "\n";
        }
        geometryStr = geometryStreamString.str();

        // point charges as string

        stringstream chargesStringStream;
        string chargesAsString;
        if (!system.pointChargeValue.empty())
        {
            for (int i = 0; i < system.pointChargeValue.size(); i++)
                chargesStringStream << setw(14) << system.pointChargePosition[i][0]
                                    << setw(14) << system.pointChargePosition[i][1]
                                    << setw(14) << system.pointChargePosition[i][2]
                                    << setw(14) << system.pointChargeValue[i] << "\n";
            chargesAsString = chargesStringStream.str();
        }

        //

        dictionary["n cpu"] = to_string(inputData.hardware.nCores);
        dictionary["memory"] = memoryAsString(inputData.hardware.totalMemoryMB);
        dictionary["method"] = qmSettings.qmMethod;
        dictionary["basis set"] = qmSettings.basisSet;
        dictionary["job name"] = inputData.jobName;
        dictionary["scf options"] = mScfOptions;
        dictionary["spin"] = to_string(system.spin_multilicity);
        dictionary["total charge"] = to_string(system.charge);
        dictionary["geometry"] = geometryStr;
        dictionary["charge"] = chargesAsString;
        dictionary["wfn file"] = inputData.wfnFileName;

        //
        string fileContent;

        string_utilities::fill_template(dictionary, inputData.qmSettings.inputTemplate, fileContent, '$');

        ofstream out(fileName);
        if (!out.good())
            on_error::throwException(string("cannot print to Gaussian input file '") + fileName + string("'"), __FILE__, __LINE__);
        out << fileContent;
        out.close();

    }


    void GaussianRunner::printInputFile(
        const std::string fileName,
        const WaveFunctionCalculationData& inputData)
        const
    {
        auto const& qmSettings = inputData.qmSettings;
        auto const& system = inputData.qmSystem;
        if (!qmSettings.inputTemplate.empty())
        {
            printInputFileFromTemplate(fileName, inputData);
            return;
        }

        ofstream out(fileName);
        if (!out.good())
            on_error::throwException(string("cannot open file task.gjf for writing"), __FILE__, __LINE__);

        if (inputData.hardware.nCores != 1)
            out << "%NProcShared=" << inputData.hardware.nCores << "\n";


        out << "%MEM=" << memoryAsString(inputData.hardware.totalMemoryMB) << "\n";
        


        if (inputData.wfnFormat == WaveFunctionFileFormat::fchk || qmSettings.tryToReadGuess)
            out << "%chk=" << inputData.jobName << ".chk\n";
        

        out << "#p " << qmSettings.qmMethod;
        if (!qmSettings.basisSet.empty())
            out << "/" << qmSettings.basisSet;
        out<< " NoSymm";

        if (!mScfOptions.empty())
            out << " SCF(" << mScfOptions << ")";


        out << " Density=Current";

        if (!qmSettings.relativisticMethod.empty())
            out << " Integral=DKH2 IOP(3/93=1)";

        if (inputData.hardware.diskSpaceGB > 1)
            out << " maxdisk=" << inputData.hardware.diskSpaceGB << "GB" << "\n";

        if (qmSettings.tryToReadGuess)
            if (fs::exists(inputData.jobName + ".chk"))
                out << " Guess=Read";

        if (!system.pointChargeValue.empty())
            out << " Charge";
        if (inputData.wfnFormat == WaveFunctionFileFormat::wfn)
            out << " OUTPUT=WFN";
        if (inputData.wfnFormat == WaveFunctionFileFormat::wfx || inputData.wfnFormat == WaveFunctionFileFormat::unspecified)
            out << " OUTPUT=WFX";

        out << "\n\ntask\n\n";

        out << system.charge << " " << system.spin_multilicity << "\n";
        out << fixed << setprecision(6);
        for (int i = 0; i < system.positions.size(); i++)
        {
            out << periodic_table::symbol(system.atomicNumbers[i]) << " ";
            out << setw(14) << system.positions[i][0]
                << setw(14) << system.positions[i][1]
                << setw(14) << system.positions[i][2] << "\n";
        }

        if (!system.pointChargeValue.empty())
        {
            out << "\n";
            for (int i = 0; i < system.pointChargeValue.size(); i++)
            {
                out << setw(14) << system.pointChargePosition[i][0]
                    << setw(14) << system.pointChargePosition[i][1]
                    << setw(14) << system.pointChargePosition[i][2]
                    << setw(14) << system.pointChargeValue[i] << "\n";
            }
        }

        if (inputData.wfnFormat != WaveFunctionFileFormat::fchk)
            out << "\n" << inputData.wfnFileName << "\n";

        out << "\n\n";
        out.close();

    }


	bool GaussianRunner::supportsFormat(WaveFunctionFileFormat format) 
		const 
	{
		vector<WaveFunctionFileFormat> formats;
		supportedFormats(formats);
		if (find(formats.begin(), formats.end(), format) == formats.end())
			return false;
		return true;
	}

	void GaussianRunner::supportedFormats(
		std::vector<WaveFunctionFileFormat>& formats) 
		const 
	{
		formats = { WaveFunctionFileFormat::fchk,
				    WaveFunctionFileFormat::wfn, WaveFunctionFileFormat::wfx };
	}

	WaveFunctionFileFormat GaussianRunner::defaultFileFormat() 
		const 
	{ 
        return mDefaultFileFormat;// WaveFunctionFileFormat::fchk;
	};



}