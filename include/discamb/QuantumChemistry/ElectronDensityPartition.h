#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"

#include "json.hpp"

#include <memory>

namespace discamb 
{

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class ElectronDensityPartition {
    public:
        virtual ~ElectronDensityPartition() = 0;
        virtual double calculate(int atom, const Vector3d& r) const = 0;
        virtual double calculate(const Vector3d& r) const = 0;
        virtual double calculate(int atom, const Vector3d& r, int molecularOrbital) const = 0;
        virtual ElectronDensityPartition* clone() const = 0;
        virtual void setAtom(int atomIdx) = 0;
        virtual void setIncludeRange(double range) = 0;
        virtual void applySettings(const nlohmann::json& data) = 0;
        
        virtual void set(
            const nlohmann::json& data, 
            const std::vector<int> &atomicNumbers,
            const std::vector<Vector3d> &atomPositions,
            std::shared_ptr<ElectronDensityCalculator> &edCalculator) = 0;
        
        static ElectronDensityPartition* create(const std::string& partitionName);// , const nlohmann::json& data);
        //static ElectronDensityPartition* std::shared_ptr<ElectronDensityPartition> create(
        //    const std::string &partitionName, 
        //    const std::string &wfxFileName, 
        //    const nlohmann::json& data);

    };
    /**@}*/
}

