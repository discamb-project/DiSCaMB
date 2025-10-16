#include "Program.h"
#include "AssignAtomTypes.h"
#include "CalcDensities.h"
#include "ChooseStructures.h"
#include "FilterStructures.h"
#include "WfnCalc.h"
#include "RunAll.h"

#include "discamb/BasicUtilities/on_error.h"

#include <fstream>

using namespace std;

Program::~Program()
{}

Program* Program::create(
    const std::string& name)
{
    if (name == "filter")
        return new FilterStructures();
    if (name == "assign")
        return new AssignAtomTypes();
    if (name == "choose")
        return new ChooseStructures();
    if (name == "all")
        return new RunAll();
    if (name == "calcwfn")
        return new WfnCalc();
    if (name == "calcden")
        return new CalcDensities();

    discamb::on_error::throwException("unknown program '" + name + "'", __FILE__, __LINE__);
    return nullptr;
}

void Program::readSettings(
    nlohmann::json& settings)
{
    settings = nlohmann::json();

    if (!filesystem::exists(filesystem::path("settings.json")))
        discamb::on_error::throwException("cannot read file settings.json", __FILE__, __LINE__);

    std::ifstream input("settings.json");
    settings = nlohmann::json::parse(input);
    input.close();

}


