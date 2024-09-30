#include "discamb/BasicUtilities/parse_cmd.h"
#include "discamb/BasicUtilities/string_utilities.h"

using namespace discamb;
using namespace std;

namespace discamb{
    namespace parse_cmd{
        void get_args_and_options(
        int nArg,
        char* args[],
        std::vector<std::string>& arguments,
        std::vector<std::string>& options)
        {
            std::map<std::string, std::string> optionsWithValues;
            get_args_and_options(nArg, args, arguments, options, optionsWithValues);

            //arguments.clear();
            //options.clear();

            //int i;

            //for (i = 1; i < nArg; i++)
            //    if (args[i][0] == '-')
            //        options.push_back(args[i]);
            //    else
            //        arguments.push_back(args[i]);
        }

        bool hasOption(
            const std::vector<std::string>& options,
            const std::string value)
        {
            return (find(options.begin(), options.end(), value) != options.end());
        }

        void get_args_and_options(
            int nArg,
            char* args[],
            std::vector<std::string>& arguments,
            std::set<std::string>& options,
            std::map<std::string, std::string>& optionsWithValues)
        {
            vector<string> _options;
            get_args_and_options(nArg, args, arguments, _options, optionsWithValues);
            options.clear();
            for (auto const& option : _options)
                options.insert(option);
        }


        void get_args_and_options(
            int nArg,
            char* args[],
            vector<string>& arguments,
            vector<string>& options,
            map<string, string>& optionsWithValues)
        {
            arguments.clear();
            options.clear();
            optionsWithValues.clear();

            vector<string> words;

            int i;

            for (i = 1; i < nArg; i++)
                if (args[i][0] == '-')
                {
                    string_utilities::split(args[i], words, '=');
                    if (words.size() == 2)
                        optionsWithValues[words[0]] = words[1];
                    else
                        options.push_back(args[i]);
                }
                else
                    arguments.push_back(args[i]);
        }


    }
}
