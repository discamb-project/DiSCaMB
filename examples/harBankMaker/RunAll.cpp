#include "RunAll.h"

using namespace std;

RunAll::RunAll()
{

}

RunAll::~RunAll()
{

}

void RunAll::set()
{
    mPrograms.push_back(shared_ptr<Program>(Program::create("filter")));
    mPrograms.push_back(shared_ptr<Program>(Program::create("assign")));
    mPrograms.push_back(shared_ptr<Program>(Program::create("choose")));
}

void RunAll::run()
{
    for (auto& program : mPrograms)
        program->run();
}


