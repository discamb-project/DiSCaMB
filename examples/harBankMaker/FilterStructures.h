#pragma once
#include "Program.h"

#include <limits>
#include <filesystem>

class FilterStructures : public Program
{
public:
    FilterStructures();
    virtual ~FilterStructures();
    virtual void set();
    virtual void run();
private:
    // percent of structures with the highest APDs
    // which should be disregarded
    double mHighAdpsStructurePercent = 0;
    // max Ueq in structure
    double mMaxUeqivalent = std::numeric_limits<double>::max();
    std::filesystem::path mResFolder;
    std::string mChosenResFolder;
};

