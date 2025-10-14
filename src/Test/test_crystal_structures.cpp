#include "discamb/Test/test_crystal_structures.h"
#include "discamb/IO/shelx_io.h"
#include <sstream>

using namespace std;

namespace {
    const char* urea_txt = R"(
TITL UREA_p1
CELL 0.71073 5.578 5.578 4.686 90 90 90
ZERR 2 0.0006 0.0006 0.0007 0 0 0
LATT -1
SYMM +Y,-X,-Z
SYMM -X,-Y,+Z
SYMM -Y,+X,-Z
SYMM 0.5+X,0.5-Y,-Z
SYMM 0.5-X,0.5+Y,-Z
SYMM 0.5-Y,0.5-X,+Z
SYMM 0.5+Y,0.5+X,+Z
SFAC C H N O
UNIT 2 8 4 2
TEMP -150(2)
WGHT 0.021479 0.001484
FVAR 3.148215
C     1     0.00000  0.50000  0.32822  10.25000  0.01465  0.01465  0.00659 =
 0.00000  0.00000 -0.00012 
O     4     0.00000  0.50000  0.59634  10.25000  0.01945  0.01945  0.00660 =
 0.00000  0.00000  0.00171 
N     3     0.14486  0.64486  0.17833  10.50000  0.02881  0.02881  0.00941 =
 0.00007  0.00007 -0.01513 
Hb    2     0.25516  0.75516  0.28543  10.50000  0.04756  0.04756  0.02374 =
 -0.00348 -0.00348 -0.02435 
Ha    2     0.14287  0.64287 -0.03608  10.50000  0.04412  0.04412  0.01955 =
 -0.00152 -0.00152 -0.01660 
HKLF 4
END
)";
    

}

namespace discamb {
    namespace test_crystal_structures {
        
        std::vector<std::string> get_test_crystal_structure_names()
        {
            vector<std::string> names{ "urea" };
            return names;
        }

        void get_test_crystal_structure(
            const std::string& name,
            Crystal& crystal)
        {
            const char* res_txt;
            if("urea" == name) {
                res_txt = urea_txt;
            } else {
                throw std::runtime_error("Unknown test crystal structure name: " + name);
            }
            stringstream input_stream(res_txt);
            discamb::shelx_io::read(input_stream, crystal);
        }
    }
}