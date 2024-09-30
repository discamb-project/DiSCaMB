#pragma once

#include <string>
#include <map>
#include <vector>

namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class OrcaExecutionTimeEstimator {
    public:
        OrcaExecutionTimeEstimator();
        OrcaExecutionTimeEstimator(const std::string &method, const std::string& basisSet);
        ~OrcaExecutionTimeEstimator();
        double executionTime(const std::map<int, int> &formula, int nCpu) const;
        double executionTime(int nBasisFunctions, int nCpu) const;
        double efficiency(double time_1_core, int nCpu) const;
        void setTheoryLevel(const std::string& method, const std::string& basisSet);
    private:
        std::string mMethod = "PBE";
        std::string mBasisSet = "cc-pVDZ";
        double oneCoreExecutionTime(const std::map<int, int>& formula) const;
        double oneCoreExecutionTime(int nBasisFunctions) const;
        
    };
    /**@}*/
}

