#include <filesystem>
#include <iostream>
#include <discamb/BasicUtilities/string_utilities.h>
#include <discamb/QuantumChemistry/OrcaRunner.h>

using namespace std;

int main(int argc, char* argv[])
{
    string folder;
    if (discamb::OrcaRunner::findOrcaFolder(folder))
        return 0;
    else
        return 1;
}

