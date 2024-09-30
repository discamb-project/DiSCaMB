#include "AtomType.h"

#include <set>
#include <map>


namespace discamb {

    struct TypeTreeNode
    {
        std::string id;
        std::string parent;
        std::set<std::string> children;
        int depth = 0;
    };

    struct TypeTree {
        std::map<std::string, TypeTreeNode> nodes;
        std::vector<std::string> roots;
        void set(const std::vector<AtomType>& types);
    };



}