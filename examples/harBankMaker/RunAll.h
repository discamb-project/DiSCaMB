#pragma once
#include "Program.h"

#include <memory>
#include <vector>

class RunAll : public Program
{
public:
    RunAll();
    virtual ~RunAll();
    virtual void set();
    virtual void run();
private:
    std::vector<std::shared_ptr<Program> > mPrograms;
};

