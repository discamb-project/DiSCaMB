#pragma once
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"


namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class OrcaRunner : public WaveFunctionDataGeneratorRunner
    {
    public:

        OrcaRunner();
        OrcaRunner(const nlohmann::json& setting);
        
        ~OrcaRunner();

        virtual std::string name() const { return std::string("ORCA"); }

        virtual void set(const nlohmann::json& settings);
        void setMolden2AimFolder(const std::string& folder);

        virtual void setExecFolder(const std::string& execFolder);
        virtual std::string getExecFolder() const;      
        virtual bool run(const WaveFunctionCalculationData& inputData) const;

        /*  */
        virtual bool supportsFormat(WaveFunctionFileFormat format) const;
        virtual void supportedFormats(std::vector<WaveFunctionFileFormat>& formats) const;
        virtual WaveFunctionFileFormat defaultFileFormat() const;
        
        void printInputFile(const std::string fileName, const std::string pointChargeFileName, const WaveFunctionCalculationData& inputData) const;
        virtual void printInputFile(const std::string fileName, const WaveFunctionCalculationData& inputData) const;

        virtual void runMultipleJobs(const std::vector<WaveFunctionCalculationData>& inputData, int nCores=0,int totalMemory=0) const;
        virtual void runMultipleJobs(const std::vector<WaveFunctionCalculationData>& inputData, int nCores, int totalMemory, std::vector<bool>& succesful) const;

        static bool succesfulRun(const std::string& jobName);
        // fixes file jobName.log
        void ecpWfxPostPorcessing(const std::string& jobName, const std::string& wfxFile);
    private:

        bool mUseBuildinEcpElectronDensityLibraries = false;

        void printInputFileFromTemplate(
            const std::string fileName,
            const WaveFunctionCalculationData& inputData) const;

        std::string mExecFullPath;
        std::string mExecFolder;
        int nPreviuslyUsedCharges;
        std::vector<int> mPreviouslyUsedAtoms;
        WaveFunctionFileFormat mDefaultFileFormat;
        bool mTryToReadGuess;
        
        void setDefaults();

        // ECP handling with molden2aim
        
        bool mHasMolden2Aim = false;
        std::string mMolden2AimFolder;
        std::string getMolden2AimExecPath(bool required = true) const;

        static void addEcpIndicatorToMoldenFile(
            const std::string& moldenFileName, 
            const std::string& outputName, 
            std::vector<std::string>& elementsWithEcp, 
            std::vector<int>& nEcpElectrons);

        static void printMolden2AimIniFile();
        static bool checkForEcp(const std::string& orcaOutput, std::vector<std::string>& elementsWithEcp, std::vector<int>& nEcpElectrons);
        
        //

        void tryToSetPathToOrcaFromSettingsFile();
        void tryToSetPathToMolden2AimFromSettingsFile();

        void makeWfnLookGaussianMade(const std::string fileName) const;
        void correctMoldenFile(const std::string &wfnFile, int nAtoms, int nCharges) const;


        // print even if pointChargeValue.empty()
        static void printPointCharges(const std::string& fName, const std::vector<double>& pointChargeValue,
            const std::vector<Vector3d>& pointChargePosition);

    };

    /**@}*/

}

