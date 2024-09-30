#pragma once

#include <string>

namespace discamb {

    struct BasisSetId {
        enum class Type {name, file, object};
        std::string value = "";
        Type type = Type::name;
    };

}