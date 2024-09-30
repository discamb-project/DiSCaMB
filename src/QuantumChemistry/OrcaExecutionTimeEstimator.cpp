#include "discamb/QuantumChemistry/OrcaExecutionTimeEstimator.h"
#include "discamb/QuantumChemistry/basis_set_data.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/MathUtilities/Vector3.h"

using namespace std;

namespace discamb {

    
    
    OrcaExecutionTimeEstimator::OrcaExecutionTimeEstimator()
    {
        on_error::not_implemented(__FILE__, __LINE__);
        
    }

    OrcaExecutionTimeEstimator::OrcaExecutionTimeEstimator(
        const std::string& method,
        const std::string& basisSet)
    {
        mMethod = method;
        mBasisSet = basisSet;
    }

    OrcaExecutionTimeEstimator::~OrcaExecutionTimeEstimator()
    {
    }

    double OrcaExecutionTimeEstimator::executionTime(
        int nBasisFunctions,
        int nCpu)
        const
    {
        double one_core_time = oneCoreExecutionTime(nBasisFunctions);
        double efficiency = this->efficiency(one_core_time, nCpu);

        return one_core_time / (nCpu * efficiency);
    }

    double OrcaExecutionTimeEstimator::executionTime(
        const std::map<int, int>& formula, int nCpu)
        const
    {
        double one_core_time = oneCoreExecutionTime(formula);
        double efficiency = this->efficiency(one_core_time, nCpu);
        
        return one_core_time / (nCpu * efficiency);
    }

    double OrcaExecutionTimeEstimator::efficiency(
        double time_1_core, 
        int nCores) const
    {
        if (nCores == 1)
            return 1.0;

        map<pair<string, int>, Vector3d> coeff{
            {{"B3LYP", 12}, Vector3d(0.50311712, 0.81875758, 66.1123164)},
            {{"B3LYP",  6}, Vector3d(0.39018626, 0.91695235, 32.66534788)},
            {{"B3LYP",  4}, Vector3d(0.33281589, 0.94677067, 20.79433368)},
            {{"B3LYP",  3}, Vector3d(0.31805961, 0.96933359, 15.07810509)},
            {{"B3LYP",  2}, Vector3d(0.31053947, 0.9837747 , 9.29503728)},
            {{"PBE",   12}, Vector3d(0.13394426, 0.53329741, 37.45643303)},
            {{"PBE",    6}, Vector3d(0.14943759, 0.75261377, 23.32396693)},
            {{"PBE",    4}, Vector3d(0.09816327, 0.83781318, 15.85447446)},
            {{"PBE",    3}, Vector3d(0.11996668, 0.88951881, 12.0391344)},
            {{"PBE",    2}, Vector3d(0.0742676 , 0.93583639, 7.39168017)}
        };

        Vector3d c;

        if (nCores == 2)
            c = coeff[{mMethod, 2}];
        if (nCores == 3)
            c = coeff[{mMethod, 3}];


        Vector3d c4 = coeff[{mMethod, 4}];
        Vector3d c6 = coeff[{mMethod, 6}];
        Vector3d c12 = coeff[{mMethod, 12}];

        if (nCores == 4)
            c = c4;
        if (nCores == 5)
            c = 0.5 * (c4 + c6);
        if (nCores >= 6 && nCores <= 12)
            c = (nCores - 6) / 6.0 * c12 + (12 - nCores) / 6.0 * c6;

        if (nCores > 12)
            c = c12;

        return (c[0] + c[1] * time_1_core) / (c[2] + time_1_core);

    }


    void OrcaExecutionTimeEstimator::setTheoryLevel(
        const std::string& method,
        const std::string& basisSet)
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }

    double OrcaExecutionTimeEstimator::oneCoreExecutionTime(
        const std::map<int, int>& formula)
        const
    {
        int nFunctions = basis_set_data::n_basis_functions(mBasisSet, formula);
        return oneCoreExecutionTime(nFunctions);

    }

    double OrcaExecutionTimeEstimator::oneCoreExecutionTime(
        int nBasisFunctions)
        const
    {
        map<pair<string, string>, vector<double> > coeff{
            { {"B3LYP", "cc-pVTZ"  }, { 1.32797487e-03, -1.84930575e-01,  9.85304092e+00} },
            { {"B3LYP", "cc-pVDZ"  }, { 8.26817241e-04,  1.61282675e-01, -1.39722368e+01} },
            { {"B3LYP", "def2-SVP" }, { 2.18054418e-03, -9.76901403e-02,  3.71166960e+00} },
            { {"B3LYP", "def2-TZVP"}, { 1.85550322e-03, -1.89184410e-01,  8.36919400e+00} },
            { {"PBE"  , "cc-pVTZ"  }, { 1.73152922e-04,  9.13867274e-03,  3.10074281e-01} },
            { {"PBE"  , "cc-pVDZ"  }, { 8.26358652e-05,  7.17283102e-02, -2.20815169e+00} },
            { {"PBE"  , "def2-SVP" }, { 3.87777513e-04,  1.62947298e-02,  1.35600678e+00} },
            { {"PBE"  , "def2-TZVP"}, { 2.95049258e-04, -1.29559761e-03,  1.47335488e+00} }
        };

        vector<double>& c = coeff[{mMethod, mBasisSet}];

        return c[0] * nBasisFunctions * nBasisFunctions + c[1] * nBasisFunctions + c[2];
    }
}