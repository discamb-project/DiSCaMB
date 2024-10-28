#pragma once

#include <map>
#include <set>
#include <string>
#include <vector>

namespace discamb{
    namespace basic_chemistry_utilities{
        
        int atomicNumberFromString(const std::string& label);
        int atomicNumberFromLabel(const std::string &label);

        // s can be an atomic number or symbol
        int atomicNumberFromString(const std::string& s);

        void getElementsList(const std::string& s, std::set<int>& atomic_numbers);
        void getFormula(const std::vector<int>& atomic_numbers, std::map<int, int>& formula);
        std::string formulaAsString(const std::map<int, int>& formula);
    }
}
