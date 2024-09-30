#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"
#include "discamb/QuantumChemistry/GaussianRunner.h"
#include "discamb/QuantumChemistry/OrcaRunner.h"

#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"


using namespace std;

namespace discamb{
    
    WaveFunctionDataGeneratorRunner *WaveFunctionDataGeneratorRunner::create(
    const string &type)
    {
        Program program = string2program(type);
        return create(program);
    }        
    
    WaveFunctionDataGeneratorRunner* WaveFunctionDataGeneratorRunner::create(
        Program program)
    {
        if (program == Program::UNDEFINED)
            on_error::throwException("cannot create WaveFunctionDataGeneratorRunner of type UNDEFINED", __FILE__, __LINE__);

        if (program == Program::GAUSSIAN)
            return new GaussianRunner;
        if (program == Program::ORCA)
            return new OrcaRunner;

        return NULL;
    }

    string WaveFunctionDataGeneratorRunner::program2string(
        WaveFunctionDataGeneratorRunner::Program program)
    {
        if (program == WaveFunctionDataGeneratorRunner::Program::GAUSSIAN)
            return "Gaussian";
        if (program == WaveFunctionDataGeneratorRunner::Program::ORCA)
            return "ORCA";
        return "undefined";
    }

    WaveFunctionDataGeneratorRunner::Program WaveFunctionDataGeneratorRunner::string2program(
        const string& _type, 
        bool throwOnError)
    {
        string type = string_utilities::toLower(_type);

        if (type == string("gaussian"))
            return WaveFunctionDataGeneratorRunner::Program::GAUSSIAN;
        if (type == string("orca"))
            return WaveFunctionDataGeneratorRunner::Program::ORCA;
        
        if(throwOnError)
            on_error::throwException(string("invalid name of wave function generator '") + _type + string("'"), __FILE__, __LINE__);

        return WaveFunctionDataGeneratorRunner::Program::UNDEFINED;
    }


    void WaveFunctionDataGeneratorRunner::setExecFolder(
        const std::string& execFolder)
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }

    std::string WaveFunctionDataGeneratorRunner::getExecFolder() const
    {
        on_error::not_implemented(__FILE__, __LINE__);
        return string();
    }
        
    //void WaveFunctionDataGeneratorRunner::runMultipleJobs(
    //    const std::vector<WaveFunctionCalculationData>& inputData,
    //    int nCores,
    //    int memory)
    //{
    //    on_error::not_implemented(__FILE__, __LINE__);
    //}

    //void WaveFunctionDataGeneratorRunner::set(
    //    const QmCalculationSettings& qmSettings)
    //{
    //    setTheoryLevel(qmSettings.qmMethod, qmSettings.basisSet, qmSettings.atomicNumber2BasisSetMap, qmSettings.relativisticMethod);
    //    setHardware(qmSettings.nCores, qmSettings.totalMemory);

    //    virtual void set(const QmCalculationSettings & qmSettings);
    //    virtual void set(const nlohmann::json & settings) = 0;

    //    virtual bool supportsFormat(WaveFunctionFileFormat format) const = 0;
    //    virtual void supportedFormats(std::vector<WaveFunctionFileFormat>&formats) const = 0;
    //    virtual WaveFunctionFileFormat defaultFileFormat() const = 0;
    //    virtual void setInputTemplate(const std::string & templateFileName);

    //    virtual std::string name() const = 0;

    //    virtual void setHardware(int nCore, int totalMemoryMB) = 0;

    //}

    //void WaveFunctionDataGeneratorRunner::setInputTemplate(const std::string& templateFileName)
    //{
    //    string errorMessage = "input template not implemented for " + this->name();
    //    on_error::throwException(errorMessage, __FILE__, __LINE__);
    //}

    WaveFunctionFileFormat WaveFunctionDataGeneratorRunner::wfnFileFormat(
        const std::string& s)
    {
        string type;
        for (char c : s)
            type += tolower(c);

        if (type == string("wfn"))
            return WaveFunctionFileFormat::wfn;

        if (type == string("wfx"))
            return WaveFunctionFileFormat::wfx;

        if (type == string("fchk"))
            return WaveFunctionFileFormat::fchk;

        if (type == string("mkl"))
            return WaveFunctionFileFormat::mkl;

        if (type == string("molden"))
            return WaveFunctionFileFormat::molden;

        on_error::throwException(string("invalid name of wave function format '") + s + string("'"), __FILE__, __LINE__);

        return WaveFunctionFileFormat::molden;
    }

    //void WaveFunctionDataGeneratorRunner::run(const std::vector<Vector3d>& positions, const std::vector<int>& atomicNumbers,
    //    int spin_multilicity, int charge, const std::string& fileName,
    //    const std::vector<double>& pointChargeValue, const std::vector<Vector3d>& pointDipoleValue, const std::vector<Vector3d>& multipoleCenters,
    //    WaveFunctionFileFormat format, const std::string& jobName,
    //    const std::map<int, std::string>& atomIdx2BasisSetMap) const
    //{
    //    on_error::not_implemented(__FILE__, __LINE__);
    //}

}
