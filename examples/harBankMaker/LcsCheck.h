#pragma once

#include "Program.h"
#include "discamb/AtomTyping/AtomType.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"


class LcsCheck : public Program
{
public:
    LcsCheck();
    virtual ~LcsCheck();
    virtual void set();
    virtual void run();
private:
    std::string mOutputFileName = std::string("lcs.log");
    std::vector< discamb::AtomType> mAtomTypes;
    discamb::DescriptorsSettings mDescriptorsSettings;

};
