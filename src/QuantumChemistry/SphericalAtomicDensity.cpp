#include "discamb/QuantumChemistry/SphericalAtomicDensity.h"

#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/IO/wfn_io.h"
#include "discamb/MathUtilities/LebedevGrid.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"
#include "discamb/QuantumChemistry/WaveFunctionDataGeneratorRunner.h"

#include <memory>
#include <filesystem>

using namespace std;



namespace discamb {

    SphericalAtomicDensity::SphericalAtomicDensity()
    {
        mN_Values = 100000;
        mStep = 0.0001;
        mValues.assign(mN_Values, 0.0);
    }

    SphericalAtomicDensity::~SphericalAtomicDensity() 
    {
    }

    void SphericalAtomicDensity::setFromWfxFile(const std::string& wfxFileName)
    {
        wfn_io::WfnFileData wfnData;
        wfn_io::read_wfx(wfxFileName, wfnData);

        if (!wfnData.edfs.empty())
            set_1rdmAndAdditionalDensity(
                wfnData.primitive_type,
                wfnData.primitive_exponents,
                wfnData.molecular_orbital_occupancy,
                wfnData.molecular_orbitals,
                wfnData.edfs[0].primitive_type,
                wfnData.edfs[0].primitive_exponents,
                wfnData.edfs[0].primitive_coefficients);
        else
            set_1rdm(
                wfnData.primitive_type,
                wfnData.primitive_exponents,
                wfnData.molecular_orbital_occupancy,
                wfnData.molecular_orbitals);
    }

    void SphericalAtomicDensity::setValues(
        const std::vector<double>& values,
        double step)
    {
        mValues = values;
        mN_Values = values.size();
        mStep = step;
    }

    void SphericalAtomicDensity::calculate(int atomicNumber, int charge, int spinMultiplicity, const QmSettings& qmSettings,
        const std::string& qmSoftware, const std::string& qm_folder, const HardwareResources& hardware,
        const nlohmann::json& qmProgramSpecificData)
    {
        if (atomicNumber == 1 && string_utilities::toLower(qmSoftware) == string("orca"))
        {
            mValues.resize(mN_Values);

            for (int i = 0; i < mN_Values; i++)
            {
                double r = i * mStep;

                mValues[i] = 1.0 / M_PI * exp(-2.0 * r);
            }
            return;
        }

        WaveFunctionCalculationData data;
        data.hardware = hardware;
        string coreName = qmSettings.qmMethod  + "_" + qmSettings.basisSet + "_" + qmSettings.relativisticMethod + 
                          "_" + to_string(atomicNumber) + "_" + to_string(spinMultiplicity);
        coreName = string_utilities::replace(coreName, ',', '_');
        data.jobName = coreName;
        data.qmProgramSpecificData = qmProgramSpecificData;
        data.qmSettings = qmSettings;
        data.qmSystem.atomicNumbers.push_back(atomicNumber);
        data.qmSystem.charge = charge;
        data.qmSystem.positions.push_back(Vector3d(0, 0, 0));
        data.qmSystem.spin_multilicity = spinMultiplicity;

        string fName = coreName + ".wfx";

        data.wfnFileName = fName;
        shared_ptr< WaveFunctionDataGeneratorRunner > runner;
        runner = shared_ptr< WaveFunctionDataGeneratorRunner >(WaveFunctionDataGeneratorRunner::create(qmSoftware));
        if (!qm_folder.empty())
            runner->setExecFolder(qm_folder);
        runner->set(qmProgramSpecificData);
        //runner->setHardware(nCores, totalMemory);
        //runner->setTheoryLevel(qm_method, basis_set, map<int, string>(), relativistic_method);
        //runner->setExecFolder(qm_folder);
        //runner->set(settings);
        //string coreName = qm_method + "_" + basis_set + "_" + relativistic_method + "_" + to_string(atomicNumber) + "_" + to_string(spinMultiplicity);
        // replace , by _ in coreName
        
        //runner->run(vector<Vector3d>(1, Vector3d(0, 0, 0)), vector<int>(1, atomicNumber), spinMultiplicity, charge, fName, coreName);

        if (!filesystem::exists(filesystem::current_path() / fName))
            runner->run(data);
            //runner->run(vector<Vector3d>(1, Vector3d(0, 0, 0)), vector<int>(1, atomicNumber), spinMultiplicity, charge, fName,
            //    vector<double>(), vector<Vector3d>(), WaveFunctionFileFormat::wfx, coreName);

        setFromWfxFile(fName);

    }


    //void SphericalAtomicDensity::calculate(int atomicNumber, int charge, int spinMultiplicity,
    //    const std::string& qmSoftware, const std::string& qm_folder, const std::string& qm_method,
    //    const std::string& basis_set, const std::string& relativistic_method, 
    //    int nCores, int totalMemory, const nlohmann::json& settings)
    //{
    //    
    //    //if (atomicNumber == 1 && string_utilities::toLower(qmSoftware) == string("orca"))
    //    //{
    //    //    mValues = hDensity;
    //    //    return;
    //    //}

    //    if (atomicNumber == 1 && string_utilities::toLower(qmSoftware) == string("orca"))
    //    {
    //        mValues.resize(mN_Values);

    //        for (int i = 0; i < mN_Values; i++)
    //        {
    //            double r = i * mStep;

    //            mValues[i] = 1.0/M_PI * exp(-2.0*r);
    //        }



    //        return;
    //    }


    //    shared_ptr< WaveFunctionDataGeneratorRunner > runner;
    //    runner = shared_ptr< WaveFunctionDataGeneratorRunner >(WaveFunctionDataGeneratorRunner::create(qmSoftware));
    //    runner->setHardware(nCores, totalMemory);
    //    runner->setTheoryLevel(qm_method, basis_set, map<int,string>(), relativistic_method);
    //    runner->setExecFolder(qm_folder);
    //    runner->set(settings);
    //    string coreName = qm_method + "_" + basis_set + "_" + relativistic_method + "_" + to_string(atomicNumber) + "_" + to_string(spinMultiplicity);
    //    // replace , by _ in coreName
    //    coreName = string_utilities::replace(coreName, ',', '_');
    //    string fName = coreName + ".wfx";
    //    //runner->run(vector<Vector3d>(1, Vector3d(0, 0, 0)), vector<int>(1, atomicNumber), spinMultiplicity, charge, fName, coreName);

    //    if(!filesystem::exists(filesystem::current_path()/ fName))
    //        runner->run(vector<Vector3d>(1, Vector3d(0, 0, 0)), vector<int>(1, atomicNumber), spinMultiplicity, charge, fName,
    //            vector<double>(), vector<Vector3d>(), WaveFunctionFileFormat::wfx, coreName);

    //    setFromWfxFile(fName);
    //}


    void SphericalAtomicDensity::set_1rdmAndAdditionalDensity(
        const std::vector<int>& primitive_type,
        const std::vector<double>& primitive_exponents,
        const std::vector<double>& orbital_occupancy,
        const std::vector<std::vector<double> >& orbitals,
        const std::vector<int>& additinal_density_primitive_type,
        const std::vector<double>& additinal_density_primitive_exponents,
        const std::vector<double>& additinal_density_primitive_coefficients)
    {
        ElectronDensityCalculator calc;
        vector<int> primitive_to_center(primitive_type.size(), 1);
        
        calc.set_1rdm(
            { {0,0,0} },
            primitive_to_center,
            primitive_type, 
            primitive_exponents,
            orbital_occupancy,
            orbitals);

        primitive_to_center.assign(additinal_density_primitive_type.size(), 1);

        calc.setAdditionalDensity(
            primitive_to_center, 
            additinal_density_primitive_type, 
            additinal_density_primitive_exponents, 
            additinal_density_primitive_coefficients);

        calculateSphericallyAveraged(calc);
    }


    void SphericalAtomicDensity::set_1rdm(
        const std::vector<int> &primitive_type,
        const std::vector<double> &primitive_exponents,
        const std::vector<double> &orbital_occupancy,
        const std::vector<std::vector<double> > &orbitals)
    {
        ElectronDensityCalculator calc;
        vector<int> primitive_to_center(primitive_type.size(), 1);
        calc.set_1rdm({ {0,0,0} }, primitive_to_center, primitive_type, primitive_exponents, orbital_occupancy, orbitals);

        calculateSphericallyAveraged(calc);

        //double r, value;
        //vector<double> w;
        //vector<Vector3d> xyz;
        //int angularGridSize = 194;
        //lebedev_laikov::get_grid(angularGridSize, xyz, w);

        //for (int i = 0; i < mN_Values; i++)
        //{
        //    r = i * mStep;
        //    
        //    value = 0;
        //    for (int j = 0; j < angularGridSize; j++)
        //        value += calc.calculate(xyz[j].x * r, xyz[j].y * r, xyz[j].z * r)*w[j];

        //    mValues[i] = value; // (4 * M_PI);
        //    
        //}
    }

    void SphericalAtomicDensity::calculateSphericallyAveraged(
        ElectronDensityCalculator& calculator)
    {
        double r, value;
        vector<double> w;
        vector<Vector3d> xyz;
        int angularGridSize = 194;
        lebedev_laikov::get_grid(angularGridSize, xyz, w);

        for (int i = 0; i < mN_Values; i++)
        {
            r = i * mStep;

            value = 0;
            for (int j = 0; j < angularGridSize; j++)
                value += calculator.calculate(xyz[j].x * r, xyz[j].y * r, xyz[j].z * r) * w[j];

            mValues[i] = value; // (4 * M_PI);
        }

    }

    double SphericalAtomicDensity::calculate(double r)
        const
    {
        int n = int(r / mStep);
        if (n > mN_Values-2)
            return 0;
        return ((r - n * mStep)*mValues[n] + ((n + 1) * mStep - r)*mValues[n + 1]) / mStep;
    }

    
        
    
}