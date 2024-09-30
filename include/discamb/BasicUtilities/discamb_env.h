#include <string>
#include "json.hpp"

namespace discamb
{
    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */

    namespace discamb_env
    {
        std::string get_discamb_path();
        nlohmann::json get_discamb_dir_json(const std::string &relative_path);
    }
    /**@}*/
}
