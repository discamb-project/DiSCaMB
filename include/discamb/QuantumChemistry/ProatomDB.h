#pragma once 

#include "json.hpp"

#include <vector>
#include <map>
#include <string>
#include <utility>

#include "SphericalAtomicDensity.h"

namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class ProatomDB
    {
    public:
        ProatomDB();
        ~ProatomDB();
        void setFromFile(const std::string &fileName);
        static bool hasDiscambDataFile(const std::string& method, const std::string& basisSet, std::string& fileName);
        bool hasEntry(int atomicNumber, int charge) const;
        bool getSphericalAtom(int atomicNumber, int charge, SphericalAtomicDensity& density) const;
    private:
        std::map<std::pair<int, int>, std::vector<double> > mDensities;
        double mStep = 0.001;
    };
    /**@}*/
}

