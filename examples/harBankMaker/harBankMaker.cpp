#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/BasicUtilities/string_utilities.h"

#include "Program.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <deque>

using namespace discamb;
using namespace std;



int main(int argc, char *argv[])
{

    try 
    {
        if (argc != 1 && argc != 2)
            on_error::throwException("expects no or one argument (execution mode)");

        string programName;
        if (argc == 1)
            programName = "all";
        else
            programName = argv[1];

        shared_ptr<Program> program = shared_ptr<Program>(Program::create(programName));
        program->set();
        program->run();
    }
    catch (exception &e)
    {
        cout << "\n" << e.what() << endl;
    }
}
