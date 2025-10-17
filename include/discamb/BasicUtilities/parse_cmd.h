#pragma once

#include <vector>
#include <string>
#include <set>
#include <map>

namespace discamb{

    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */


    namespace parse_cmd{
        /**
        handles command line arguments and options
        Arguments are strings not starting with '-'
        Options are strings starting with '-'
        Example:    
        program_name arg1 arg2 -opt1 -opt2 -opt3 val3 arg3
        Arguments: arg1, arg2, arg3
        Options: -opt1, -opt2
        */
        void get_args_and_options(
                int nArg,
                char* args[],
                std::vector<std::string>& arguments,
                std::vector<std::string>& options);
        /**
        handles command line arguments and options
        Arguments are strings not starting with '-'
        Options are strings starting with '-'
        Options with values are strings starting with '-' immidiately followed by = and a value
        Example:
        program_name arg1 arg2 -opt1 -opt2 -opt3=x val3 arg3
        Arguments: arg1, arg2, arg3
        Options: -opt1, -opt2
        Options with values: -opt3=x
        */

        void get_args_and_options(
            int nArg,
            char* args[],
            std::vector<std::string>& arguments,
            std::vector<std::string>& options,
            std::map<std::string,std::string>& optionsWithValues);
        /**
        returns true if option is present in options
        Example:
        program_name arg1 arg2 -opt1 -opt2 -opt3=x val3
        Has option -opt2 -> true
        */
        bool hasOption(const std::vector<std::string>& options, const std::string value);

        void get_args_and_options(
            int nArg,
            char* args[],
            std::vector<std::string>& arguments,
            std::set<std::string>& options,
            std::map<std::string, std::string>& optionsWithValues);

    }
    /**@}*/
}
