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
        void get_args_and_options(
                int nArg,
                char* args[],
                std::vector<std::string>& arguments,
                std::vector<std::string>& options);

        void get_args_and_options(
            int nArg,
            char* args[],
            std::vector<std::string>& arguments,
            std::vector<std::string>& options,
            std::map<std::string,std::string>& optionsWithValues);

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
