#include <map>
#include <vector>
#include <string>

namespace discamb
{

    /**
    * \addtogroup IO IO
    * @{
    */

    namespace proatom_db_io
    {
        void read_proatom_db(
            const std::string& fileName,
            std::vector<int>& atomicNumber, 
            std::vector<int>& charge,
            std::vector<std::vector<double> >& data);

        void write_proatom_db(
            const std::string& fileName,
            std::vector<int>& atomicNumber,
            std::vector<int>& charge,
            std::vector<std::vector<double> >& data);

    }
    /**@}*/
}

