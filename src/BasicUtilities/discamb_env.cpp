#include "discamb/BasicUtilities/discamb_env.h"


#include <cstdlib>
#include <filesystem>
#include <fstream>

using namespace std;
namespace fs = std::filesystem;

namespace discamb
{
    namespace discamb_env
    {
        std::string get_discamb_path()
        {

            
            if (const char* c = getenv("DISCAMB_PATH"))
                return string(c);

            return string();
        }

        nlohmann::json get_discamb_dir_json(
            const std::string& relative_path)
        {
            string discamb_path = discamb_env::get_discamb_path();
            fs::path jsonPath = fs::path(discamb_path) / fs::path(relative_path);
            nlohmann::json result;

            if (fs::exists(jsonPath))
            {
                ifstream in(jsonPath.string());
                if (in.good())
                    in >> result;
                in.close();
            }

            return result;
        }
    }

    
}
