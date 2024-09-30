#pragma once

#include <string>

namespace discamb {

    class ChemicalElement {
    public:
        ChemicalElement();
        ChemicalElement(int atomicNumber);
        ChemicalElement(const std::string &symbol);
        ~ChemicalElement();
        int atomicNumber() const;
        void set(int z);
        std::string symbol() const;
        void set(const std::string &s);
    private:
        std::string mSymbol;
        int mAtomicNumber;
    };

}
