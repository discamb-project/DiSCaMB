#include <vector>
#include <map>
#include <string>

namespace discamb
{

/**
* \addtogroup QuantumChemistry
* @{
*/


    namespace basis_set_data {
        /**
        n_basis_functions[z] - number of basis functions for element with atomic number z, size of the vector is
        1 + maximum atomic number described by the basis set
        */
        void n_basis_functions(const std::string& basis_set_name, std::vector<int>& n_basis_functions);
        int n_basis_functions(const std::string& basis_set_name, const std::map<int, int>& formula);
        int n_basis_functions(const std::string& basis_set_name, int atomic_number);
    }
    /**@}*/
}
