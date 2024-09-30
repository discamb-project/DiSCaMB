#pragma once
#include "discamb/QuantumChemistry/WaveFunctionFileFormat.h"
#include "discamb/QuantumChemistry/WaveFunctionCalculationData.h"
#include "discamb/MathUtilities/Vector3.h"
#include "json.hpp"

#include <vector>
#include <string>

namespace discamb{
    
    /**
    * \defgroup QuantumChemistry QuantumChemistry
    \brief Quantum chemistry.
    * @{
    */




    class WaveFunctionDataGeneratorRunner
    {
        public:
	            
            virtual void set(const nlohmann::json &settings) = 0;
            virtual void setExecFolder(const std::string& execFolder);
            virtual std::string getExecFolder() const;
            virtual bool supportsFormat(WaveFunctionFileFormat format) const = 0;
            virtual void supportedFormats(std::vector<WaveFunctionFileFormat>& formats) const = 0;
            virtual WaveFunctionFileFormat defaultFileFormat() const = 0;
            virtual void printInputFile(const std::string fileName, const WaveFunctionCalculationData& inputData) const=0;

            virtual std::string name() const = 0;

            virtual bool run(const WaveFunctionCalculationData& inputData) const = 0;

            virtual void runMultipleJobs(const std::vector<WaveFunctionCalculationData>& inputData, int nCores,int totalMemory) const = 0;
            virtual void runMultipleJobs(const std::vector<WaveFunctionCalculationData>& inputData, int nCores, int totalMemory, std::vector<bool> &succesful) const = 0;

            static std::string wfnFormatAsString(WaveFunctionFileFormat format) {
                if (format == WaveFunctionFileFormat::fchk)
                    return std::string("fchk");
                if (format == WaveFunctionFileFormat::molden)
                    return std::string("molden");
                if (format == WaveFunctionFileFormat::wfn)
                    return std::string("wfn");
                if (format == WaveFunctionFileFormat::wfx)
                    return std::string("wfx");
                return std::string();
            }
            
            static WaveFunctionFileFormat wfnFileFormat(const std::string& s);
            /*
            types: gaussian, orca, // to do nwchem, psi4, gamess,qchem
            */
            static WaveFunctionDataGeneratorRunner *create(const std::string &type);
            enum class Program {GAUSSIAN, ORCA, UNDEFINED};
            static WaveFunctionDataGeneratorRunner* create(Program program);
            static std::string program2string(Program program);
            static Program string2program(const std::string& type, bool throwOnError = true);

    };
    /**@}*/
}

