#include "discamb/CrystalStructure/structure_library.h"
#include "discamb/IO/shelx_io.h"

#include <map>
#include <sstream>

using namespace std;

namespace {

    
const string urea = R"(TITL UREA
CELL 0.71073 5.578 5.578 4.686 90 90 90
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

SIZE 0.23 0.13 0.09
FVAR 3.151334

C     1     0.00000  0.50000  0.32823  10.25000  0.01456  0.01456  0.00668 =
 0.00000 -0.00000  0.00005 
O     4     0.00000  0.50000  0.59620  10.25000  0.01926  0.01926  0.00662 =
 0.00000 -0.00000  0.00165 
N     3     0.14480  0.64480  0.17849  10.50000  0.02846  0.02846  0.00948 =
 0.00006  0.00006 -0.01461 
H1    2     0.25507  0.75507  0.28463  10.50000  0.05388  0.05388  0.01832 =
 -0.00472 -0.00472 -0.03662 
H2    2     0.14163  0.64163 -0.03415  10.50000  0.04788  0.04788  0.01602 =
 -0.00127 -0.00127 -0.01486 
HKLF 4

END
)";

const std::map<std::string, std::string> label2structure = { {"urea", urea} };

}

namespace discamb {

    namespace structure_library {

        bool getStructure(
            const std::string& name, 
            Crystal& structure)
        { 
            structure = Crystal();

            if (label2structure.find(name) != label2structure.end())
            {
                stringstream ss(label2structure.find(name)->second);    
                shelx_io::read(ss, structure);
                return true;
            }
            return false;
        }

        bool hasStructure(const std::string& name)
        {
            if (label2structure.find(name) != label2structure.end())
                return true;
            return false;
        }

        void structureList(
            std::set<std::string>& structureList)
        {
            structureList.clear();
            for (auto const& item : label2structure)
                structureList.insert(item.first);
        }

    }
}
