#pragma once

#include "discamb/QuantumChemistry/ElectronDensityPartition.h"
#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"
#include "discamb/QuantumChemistry/PromoleculeElectronDensity.h"
#include <memory>


namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */

    
    class HirshfeldPartition: public ElectronDensityPartition
    {
    public:
        HirshfeldPartition();
        HirshfeldPartition(std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
                           std::shared_ptr<PromoleculeElectronDensity> &promoleculeDensity);
        HirshfeldPartition(const std::string& wfnFile,
            const std::string& atomicDensitiesFile,
            const std::string& method = std::string(),
            const std::string& basisSet = std::string());

        //HirshfeldPartition(const std::string& wfnFile,
        //    const std::string& method,
        //    const std::string& basisSet);


        HirshfeldPartition(
            const std::string& wfnFile,
            const std::string& method,
            const std::string& basisSet,
            const std::vector<std::string>& basis_sets,
            const std::map<int, int>& atomicNumberToBasisSetMap,
            const std::map<int, int>& atomToBasisSetMap,
            bool calculateSphericalDensityIfMissing = true,
            const std::string& qmProgram = std::string("orca"),
            const std::string& qmProgramFolder = std::string("c:\\Orca"),
            const std::string& relativisticMethod = std::string(),
            int nCores = 1,
            int mamoryMB = 1000);


        virtual ~HirshfeldPartition() override;

        virtual void set(
            const nlohmann::json& data, 
            const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& atomPositions,
            std::shared_ptr<ElectronDensityCalculator>& edCalculator);

        void set(std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
            std::shared_ptr<PromoleculeElectronDensity> &promoleculeDensity);

        void setFromFile(
            const std::string& wfn,
            const std::string& atomicDensitiesFile,
            const std::string& method = std::string(),
            const std::string& basisSet = std::string());


        bool setFromDefaultFiles(
            const std::string& wfn,
            const std::string& method,
            const std::string& basisSet,
            const std::vector<std::string>& basis_sets,
            const std::map<int, int>& atomicNumberToBasisSetMap,
            const std::map<int, int>& atomToBasisSetMap,
            bool calculateSphericalDensityIfMissing = true,
            const std::string &qmProgram = std::string("orca"),
            const std::string& qmProgramFolder = std::string("c:\\Orca"),
            const std::string& relativisticMethod = std::string(),
            int nCores = 1,
            int memoryMB = 1000);
            


        virtual void applySettings(const nlohmann::json& data);

        virtual double calculate(int atom, const Vector3d &r) const;
        double calculateDifferenceDensity(int atom, const Vector3d& r) const;
        virtual double calculate(int atom, const Vector3d& r,int molecularOrbital) const override;
        virtual void setIncludeRange(double range) override;
        void unsetIncludeRange();
        virtual void setAtom(int atomIdx) override;
        virtual double calculate(const Vector3d& r) const { return calculate(mAtomIdx, r); };
        HirshfeldPartition* clone() const override;
        void setPower(double p);
    private:
        int mN_Atoms = 0;
        int mAtomIdx = 0;
        double mPower = 1.0;
        double mIncludeRange = -1;
        std::vector<std::vector<int> > mIncludeAtoms;
        std::vector<Vector3d> mAtomPositions;
        std::shared_ptr<ElectronDensityCalculator> mEdCalculator;
        std::shared_ptr<PromoleculeElectronDensity> mPromoleculeDensity;
        void setEdCalculatorFromWfnFile(const std::string &file, std::vector<Vector3d> &positions, std::vector<int>& atomicNumbers);
    };
    /**@}*/
}