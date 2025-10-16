#include "Program.h"
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"



#include <vector>
#include <string>
#include <memory>
#include <limits>

class WfnCalc : public Program
{
public:
    WfnCalc();
    virtual ~WfnCalc();
    virtual void set();
    virtual void run();
private:
    discamb::QmSettings mQmSettings;
    std::vector<std::string> mStructureName;
    std::vector<int> mMolIdxInStructure;
    std::string mQmProgram = std::string("orca");
    discamb::HardwareResources mHardware;
    std::shared_ptr<discamb::WaveFunctionDataGeneratorRunner> mRunner;
    int mMaxTimeMinutes = std::numeric_limits<int>::max();
};

