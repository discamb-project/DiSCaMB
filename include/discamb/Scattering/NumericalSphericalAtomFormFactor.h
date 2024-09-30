#pragma once
#include "discamb/QuantumChemistry/SphericalAtomicDensity.h"


namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class NumericalSphericalAtomFormFactor{
    public:
        NumericalSphericalAtomFormFactor(const SphericalAtomicDensity &sphericalDensity);
        virtual ~NumericalSphericalAtomFormFactor();
        void set(const SphericalAtomicDensity& sphericalDensity);
        // h in 1/A
        double calculate(double h) const;

    private:
        std::vector<double> mValues;
        double mStart, mStep;
    };
    /** @}*/
}

