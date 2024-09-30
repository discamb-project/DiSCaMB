#pragma once
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"


namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


	class GaussianRunner : public WaveFunctionDataGeneratorRunner
	{
	public:

		GaussianRunner();
		GaussianRunner(const nlohmann::json& setting);
		GaussianRunner(const std::string& executionFolder, const std::string& execName);
		~GaussianRunner();

        virtual std::string name() const { return std::string("Gaussian"); }

  //      virtual void setHardware(int nCore, int totalMemoryMB);

		//virtual void setTheoryLevel(
  //          const std::string& method, 
  //          const std::string& basis_set,
  //          const std::map<int, std::string>& atomicNumber2BasisSetMap = std::map<int, std::string>(),
  //          const std::string& relativisticMethod = std::string());

		virtual void set(const nlohmann::json& settings);
        //virtual void setInputTemplate(const std::string& templateFileName);
        virtual void setExecFolder(const std::string& execFolder);
        virtual std::string getExecFolder() const;
		void setExec(const std::string& executionFolder, const std::string& execName) const;
        void getExec(std::string& executionFolder, std::string& execName) const;
        void setScfOptions(const std::string& scfOptions);

		//virtual void run(const std::vector<Vector3d>& positions, const std::vector<int>& atomicNumbers,
		//	int spin_multilicity, int charge, const std::string& fileName,
		//	WaveFunctionFileFormat format,  const std::string& jobName, 
  //          const std::map<int, std::string>& atomIdx2BasisSetMap = std::map<int, std::string>()) const;

		//virtual void run(const std::vector<Vector3d>& positions, const std::vector<int>& atomicNumbers,
		//	int spin_multilicity, int charge, const std::string& fileName,
		//	const std::vector<double>& pointChargeValue, const std::vector<Vector3d>& pointChargePosition,
		//	WaveFunctionFileFormat format, const std::string& jobName,
  //          const std::map<int, std::string>& atomIdx2BasisSetMap = std::map<int, std::string>()) const;

        virtual bool run(const WaveFunctionCalculationData& inputData) const;
        virtual void runMultipleJobs(const std::vector<WaveFunctionCalculationData>& inputData, int nCores, int memoryPerCore) const;
        virtual void runMultipleJobs(const std::vector<WaveFunctionCalculationData>& inputData, int nCores, int totalMemory, std::vector<bool>& succesful) const;
		/*  */
		virtual bool supportsFormat(WaveFunctionFileFormat format) const;
		virtual void supportedFormats(std::vector<WaveFunctionFileFormat>& formats) const;
		virtual WaveFunctionFileFormat defaultFileFormat() const;
        
        virtual void printInputFile(const std::string fileName, const WaveFunctionCalculationData& inputData) const;

        //void printInputFile(const std::string fileName,
        //    const std::vector<Vector3d>& positions,
        //    const std::vector<int>& atomicNumbers,
        //    int spin_multilicity,
        //    int charge,
        //    const std::string& wfnFileName,
        //    const std::vector<double>& pointChargeValue,
        //    const std::vector<Vector3d>& pointChargePosition,
        //    WaveFunctionFileFormat format,
        //    const std::string& jobName) const;

        static std::string memoryAsString(int memoryMB);
        static std::string findGaussianExecName(const std::string &folder);
        static bool succesfulRun(const std::string& jobName, const std::string& wfnName);
	private:
        void printInputFileFromTemplate(
            const std::string fileName,
            const WaveFunctionCalculationData& inputData) const;
        
        //void printInputFileFromTemplate(
        //    const std::string fileName,
        //    const std::vector<Vector3d>& positions,
        //    const std::vector<int>& atomicNumbers,
        //    int spin_multilicity,
        //    int charge,
        //    const std::string& wfnFileName,
        //    const std::vector<double>& pointChargeValue,
        //    const std::vector<Vector3d>& pointChargePosition,
        //    WaveFunctionFileFormat format,
        //    const std::string& jobName) const;
        bool mChangeDir;
        bool mGenerateOutput;
        
		mutable std::string mExecName;
        mutable std::string mExecFolder;
		//std::string mMethod;
  //      std::string mMemory;
        std::string mScfOptions;
  //      int mN_Core;
		//std::string mBasisSet;
  //      std::string mRelativisticMethod;
  //      std::map<int, std::string> mAtomicNumber2BasisSetMap;
        WaveFunctionFileFormat mDefaultFileFormat;
        //bool mTryToReadGuess;
        void setDefaults();
        /** finds out if mExecName*/
        void checkPaths() const;
        //std::string mInputTemplate;
        //bool mTemplateHasCharges;



	};

    /**@}*/

}

