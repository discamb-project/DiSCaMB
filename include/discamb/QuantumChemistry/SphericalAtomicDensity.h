#pragma once

#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"
#include "discamb/QuantumChemistry/WaveFunctionCalculationData.h"

#include "json.hpp"

#include <cstddef>
#include <vector>

namespace discamb
{
    
    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class SphericalAtomicDensity {
    public:
        SphericalAtomicDensity();
        ~SphericalAtomicDensity();

        void setFromWfxFile(const std::string& wfxFileName);

        void set_1rdm(
            const std::vector<int> &primitive_type,
            const std::vector<double> &primitive_exponents,
            const std::vector<double> &orbital_occupancy,
            const std::vector<std::vector<double> > &orbitals);
        
        void set_1rdmAndAdditionalDensity(
            const std::vector<int>& primitive_type,
            const std::vector<double>& primitive_exponents,
            const std::vector<double>& orbital_occupancy,
            const std::vector<std::vector<double> >& orbitals,
            const std::vector<int>& additinal_density_primitive_type,
            const std::vector<double>& additinal_density_primitive_exponents,
            const std::vector<double>& additinal_density_primitive_coefficients);


        /**
        step in a.u.
        */
        void setValues(const std::vector<double>& values, double step);

        //void calculate(int atomicNumber, int charge, int spinMultiplicity, 
        //    const std::string& qmSoftware, const std::string& qm_folder, const std::string& qm_method,
        //    const std::string& basis_set, const std::string& relativistic_method, int nCores, int totalMemory, 
        //    const nlohmann::json& settings = nlohmann::json());

        void calculate(int atomicNumber, int charge, int spinMultiplicity, const QmSettings& qmSettings, 
            const std::string& qmSoftware, const std::string& qm_folder, const HardwareResources &hardware,
            const nlohmann::json& settings = nlohmann::json());


        /**
        r in a.u.
        */
        double calculate(double r) const;
    private:
        std::vector<double> mValues;
        int mN_Values;
        double mStep;
        void calculateSphericallyAveraged(ElectronDensityCalculator &calculator);
    };
    /**@}*/
}
