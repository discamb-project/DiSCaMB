#include "discamb/BasicChemistry/ChemicalElement.h"
#include "discamb/BasicChemistry/periodic_table.h"


namespace discamb {

    ChemicalElement::ChemicalElement()
    {
        mAtomicNumber = 0;
        mSymbol = "X";
    }

    ChemicalElement::ChemicalElement(int  atomicNumber)
    {
        set(atomicNumber);
    }

    ChemicalElement::ChemicalElement(
        const std::string& symbol) 
    {
        set(symbol);
    }

    ChemicalElement::~ChemicalElement() 
    {
    }
    
    int ChemicalElement::atomicNumber()
        const 
    { 
        return mAtomicNumber; 
    }

    void ChemicalElement::set(
        int  atomicNumber)
    {
        mAtomicNumber = atomicNumber;
        mSymbol = periodic_table::symbol(atomicNumber);
    }
    
    std::string ChemicalElement::symbol() 
        const 
    {
        return mSymbol;
    }
    
    void ChemicalElement::set(
        const std::string& symbol)
    {
        mSymbol = symbol;
        mAtomicNumber = periodic_table::atomicNumber(symbol);
    }

}

//private:
//    std::string mSymbol;
//    int mAtomicNumber;
//};
