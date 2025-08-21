#include "CrystalStructure.h"

namespace taam_utilities {



    std::vector<int> find_atom_types(
        const CrystalStructure& crystal_structure,
        const std::string &bank_path);

    std::vector<std::vector<int> > find_atom_types_for_fragments(
        const CrystalStructure& crystal_structure,
        const std::string& bank_path,
        const std::vector< std::vector<std::vector<int> > >& fragments);

    std::vector < std::string > get_atom_type_names(const std::string& bank_path);

}

