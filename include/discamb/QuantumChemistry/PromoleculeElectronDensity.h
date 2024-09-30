#pragma once
#include "discamb/QuantumChemistry/SphericalAtomicDensity.h"
#include "discamb/QuantumChemistry/ProatomDB.h"
#include "discamb/MathUtilities/Vector3.h"

#include "json.hpp"

#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class PromoleculeElectronDensity
    {
    public:
        PromoleculeElectronDensity(const std::vector<SphericalAtomicDensity> &densities,
                                   const std::vector<Vector3d> &positions,
                                   const std::vector<int> &atomToDensityMap);

        PromoleculeElectronDensity(
            const ProatomDB &db,
            const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& positions);

        PromoleculeElectronDensity(
            const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& positions,
            const nlohmann::json& settings);

        PromoleculeElectronDensity();

        ~PromoleculeElectronDensity();
        
        PromoleculeElectronDensity* clone() const;

        int nAtoms() const { return mPositions.size(); }

        void set(const std::vector<SphericalAtomicDensity> &densities,
                 const std::vector<Vector3d> &positions,
                 const std::vector<int> &atomToDensityMap);

        void set(const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& positions,
            const nlohmann::json& settings);


        void get(std::vector<SphericalAtomicDensity>& densities,
                 std::vector<Vector3d>& positions,
                 std::vector<int>& atomToDensityMap);


        void set(const ProatomDB& db, 
            const std::vector <int>& atomicNumbers,
            const std::vector<Vector3d>& positions);


        void set(const std::string& sphericalDensitiesFile,
            const std::vector <int>& atomicNumbers,
            const std::vector<Vector3d>& positions);

        bool setFromDefaultFiles(
            const std::vector <int>& atomicNumbers,
            const std::vector<Vector3d>& positions,
            const std::string& method,
            const std::string& basisSet,
            const std::vector<std::string>& basis_sets = std::vector<std::string>(),
            const std::map<int, int>& atomicNumberToBasisSetMap = std::map<int, int>(),
            const std::map<int, int>& atomToBasisSetMap = std::map<int, int>(),
            bool calculateIfMissing = true,
            const std::string& qmProgram = std::string("orca"),
            const std::string& qmProgramFolder = std::string("c:\\Orca"),
            const std::string& relativisticMethod = std::string(),
            int nCores = 1,
            int memoryMB = 1000);
        

        void setGeometry(const std::vector<Vector3d> &positions);
        void getGeometry(std::vector<Vector3d>& positions) const;
        double calculate(const Vector3d &r) const;
        double calculate(const Vector3d &r, int atomIdx) const;
    private:
        
        std::vector<SphericalAtomicDensity> mDensities;
        std::vector<Vector3d> mPositions;
        std::vector<int> mAtomToDensityMap;
    };
    /**@}*/
}

