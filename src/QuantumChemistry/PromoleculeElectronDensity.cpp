#include "discamb/QuantumChemistry/PromoleculeElectronDensity.h"
#include "discamb/QuantumChemistry/WaveFunctionCalculationData.h"
#include "discamb/BasicChemistry/chemical_element_data.h"

#include "discamb/IO/proatom_db_io.h"
#include "discamb/BasicUtilities/OnError.h"

#include <cmath>
#include <set>
#include <map>

using namespace std;

namespace discamb {

    PromoleculeElectronDensity::PromoleculeElectronDensity()
    {

    }

    PromoleculeElectronDensity::PromoleculeElectronDensity(
        const std::vector<SphericalAtomicDensity> &densities,
        const std::vector<Vector3d> &positions,
        const std::vector<int> &atomToDensityMap)
    {
        set(densities, positions, atomToDensityMap);
    }

    PromoleculeElectronDensity::PromoleculeElectronDensity(
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& positions,
        const nlohmann::json& settings)
    {
        set(atomicNumbers, positions, settings);
    }


    PromoleculeElectronDensity::PromoleculeElectronDensity(
        const ProatomDB& db,
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& positions)
    {
        set(db, atomicNumbers, positions);
    }


    PromoleculeElectronDensity::~PromoleculeElectronDensity()
    {

    }

    PromoleculeElectronDensity* PromoleculeElectronDensity::clone() const
    {
        return new PromoleculeElectronDensity(*this);
    }

    void PromoleculeElectronDensity::set(
        const std::string& sphericalDensitiesFile,
        const std::vector <int>& atomicNumbers,
        const std::vector<Vector3d>& positions)
    {
        ProatomDB db;
        db.setFromFile(sphericalDensitiesFile);
        set(db, atomicNumbers, positions);
    }

    void PromoleculeElectronDensity::set(
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& positions,
        const nlohmann::json& settings)
    {
        if (settings.is_null())
            on_error::throwException("missing 'settings' when calling PromoleculeElectronDensity::set", __FILE__, __LINE__);


        string sphericalAtomsFile = settings.value("atoms file", string());
        if (!sphericalAtomsFile.empty())
        {
            set(sphericalAtomsFile, atomicNumbers, positions);
            return;
        }

        QmSettings qmSettings;
        qmSettings.set(settings);

        HardwareResources hardware;
        hardware.set(settings);

        string qm_program = settings.value("qm program", string());

        vector<string> basis_sets;
        map<int, int> atomicNumberToBasisSetMap;


        if (!qmSettings.atomicNumber2BasisSetMap.empty())
        {
            for (auto& item : qmSettings.atomicNumber2BasisSetMap)
            {
                auto it = find(basis_sets.begin(), basis_sets.end(), item.second);
                if (it != basis_sets.end())
                    atomicNumberToBasisSetMap[item.first] = distance(basis_sets.begin(), it);
                else
                {
                    atomicNumberToBasisSetMap[item.first] = basis_sets.size();
                    basis_sets.push_back(item.second);
                }
            }
        }

        setFromDefaultFiles(
            atomicNumbers,
            positions,
            qmSettings.qmMethod,
            qmSettings.basisSet,
            basis_sets,
            atomicNumberToBasisSetMap,
            map<int, int>(),
            true,
            qm_program,
            settings.value("qm folder", string()),
            qmSettings.relativisticMethod,
            hardware.nCores,
            hardware.totalMemoryMB);
    }




    bool PromoleculeElectronDensity::setFromDefaultFiles(
        const std::vector <int>& atomicNumbers,
        const std::vector<Vector3d>& positions,
        const std::string& method,
        const std::string& basisSet,
        const std::vector<std::string>& _basis_sets,
        const std::map<int, int>& atomicNumberToBasisSetMap,
        const std::map<int, int>& atomToBasisSetMap,
        bool calculateIfMissing,
        const std::string& qmProgram,
        const std::string& qmProgramFolder,
        const std::string& relativisticMethod,
        int nCores,
        int memoryMB)
    {
        int atomIdx, nAtoms = atomicNumbers.size();
        vector<string> basis_sets = _basis_sets;

        basis_sets.push_back(basisSet);
        int basisSetIdx, nBasisSets = basis_sets.size();
        vector<int> atomIdx2BasisSet(nAtoms, nBasisSets-1);

        for (auto z2idx : atomicNumberToBasisSetMap)
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (atomicNumbers[atomIdx] == z2idx.first)
                    atomIdx2BasisSet[atomIdx] = z2idx.second;

        for (auto atom2base : atomToBasisSetMap)
            atomIdx2BasisSet[atom2base.first] = atom2base.second;

        vector<ProatomDB> proatomDBs(nBasisSets);

        for (basisSetIdx = 0; basisSetIdx < nBasisSets; basisSetIdx++)
        {
            string fileName;
            if (ProatomDB::hasDiscambDataFile(method, basis_sets[basisSetIdx], fileName))
                proatomDBs[basisSetIdx].setFromFile(fileName);
            //else
              //  return false;
        }
        
        vector<int> atomToSphDensityMap(nAtoms);
        map<pair<int, int>, vector<int> > dbEntry2Atoms;
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            dbEntry2Atoms[{atomIdx2BasisSet[atomIdx], atomicNumbers[atomIdx]}].push_back(atomIdx);

        int idx = 0;
        vector<SphericalAtomicDensity> sphericalDensities;
        SphericalAtomicDensity sphericalDensity;
        for (auto item : dbEntry2Atoms)
        {
            int basisSetIdx = item.first.first;
            if (!proatomDBs[item.first.first].getSphericalAtom(item.first.second, 0, sphericalDensity))
            {
                if(!calculateIfMissing)
                    on_error::throwException("not spherical atom density for Z = " + to_string(item.first.second), __FILE__, __LINE__);
                else
                {
                    int z = item.first.second;
                    int spinMultiplicity = chemical_element_data::groundStateSpinMultiplicity(z);
                    QmSettings qmSettings;
                    qmSettings.basisSet = basis_sets[basisSetIdx];
                    qmSettings.qmMethod = method;
                    qmSettings.relativisticMethod = relativisticMethod;
                    HardwareResources hardware;
                    hardware.nCores = nCores;
                    hardware.totalMemoryMB = memoryMB;
                    sphericalDensity.calculate(z, 0, spinMultiplicity, qmSettings, qmProgram, qmProgramFolder, hardware);
                    //sphericalDensity.calculate(z, 0, spinMultiplicity, qmProgram, qmProgramFolder, method, basis_sets[basisSetIdx], relativisticMethod, nCores, memoryMB);
                }
            }
            sphericalDensities.push_back(sphericalDensity);
            for (auto i : item.second)
                atomToSphDensityMap[i] = idx;
            idx++;
        }
        
        set(sphericalDensities, positions, atomToSphDensityMap);

        return true;
    }


    void PromoleculeElectronDensity::set(
        const ProatomDB& db,
        const std::vector <int>& atomicNumbers,
        const std::vector<Vector3d>& positions)
    {
        vector<SphericalAtomicDensity> densities;
        vector<int> atomToDensityMap;

        map<int, int> z2idx;
        std::set<int> uniqueZ(atomicNumbers.begin(), atomicNumbers.end());
        int idx = 0;
        densities.resize(uniqueZ.size());
        atomToDensityMap.resize(atomicNumbers.size());
        for (auto z : uniqueZ)
        {
            if (!db.getSphericalAtom(z, 0, densities[idx]))
                on_error::throwException(
                    string("no spherical atom density for file with atomic number ") + to_string(z),
                    __FILE__, __LINE__);
            z2idx[z] = idx++;
        }
        idx = 0;
        for (auto z : atomicNumbers)
            atomToDensityMap[idx++] = z2idx[z];
        set(densities, positions, atomToDensityMap);
    }


    void PromoleculeElectronDensity::set(
        const vector<SphericalAtomicDensity> &densities,
        const vector<Vector3d> &positions,
        const vector<int> &atomToDensityMap)
    {
        mDensities = densities;
        mPositions = positions;
        mAtomToDensityMap = atomToDensityMap;
    }

    void PromoleculeElectronDensity::get(
        std::vector<SphericalAtomicDensity>& densities,
        std::vector<Vector3d>& positions,
        std::vector<int>& atomToDensityMap)
    {
        densities = mDensities;
        positions = mPositions;
        atomToDensityMap = mAtomToDensityMap;
    }



    //void PromoleculeElectronDensity::set(
    //    const std::string& sphericalDensitiesFile,
    //    const std::vector <int>& atomicNumbers,
    //    const std::vector<Vector3d>& positions)
    //{
    //    int index;
    //    vector<int> dbAtomicNumber, atomToDensityMap;
    //    vector <int> charge;
    //    vector<vector<double> > dbData, subsetData;
    //    std::set<int> uniqueZ;
    //    proatom_db_io::read_proatom_db(sphericalDensitiesFile, dbAtomicNumber, charge, dbData);



    //    for (auto z : atomicNumbers)
    //        uniqueZ.insert(z);
    //    
    //    
    //    vector<SphericalAtomicDensity> densities;
    //    map<int, int> zToDbEntry, zToDensityIdx;
    //    SphericalAtomicDensity sphericalDensity;
    //    

    //    for (index = 0; index < dbAtomicNumber.size(); index++)
    //        if (charge[index] == 0)
    //            zToDbEntry[dbAtomicNumber[index]] == index;

    //    index = 0;
    //    for (auto z : uniqueZ)
    //    {
    //        auto it = zToDbEntry.find(z);
    //        if (it == zToDbEntry.end())
    //            on_error::throwException(string("no data for element with atomic number ") + to_string(z) +
    //                string(" in spherical densities  file '") + sphericalDensitiesFile + string("'"),
    //                __FILE__, __LINE__);
    //        sphericalDensity.setValues(dbData[it->second], 0.001);
    //        densities.push_back(sphericalDensity);
    //        zToDensityIdx[z] = index;
    //        index++;
    //    }

    //    for (auto z : atomicNumbers)
    //        atomToDensityMap.push_back(zToDensityIdx[z]);

    //    set(densities, positions, atomToDensityMap);

    //}


    double PromoleculeElectronDensity::calculate(
        const Vector3d &r,
        int atomIdx) 
        const
    {
        Vector3d rLocal = r - mPositions[atomIdx];
        return mDensities[mAtomToDensityMap[atomIdx]].calculate(sqrt(rLocal*rLocal));
    }

    void PromoleculeElectronDensity::setGeometry(
        const std::vector<Vector3d> &positions)
    {
        mPositions = positions;
    }

    void PromoleculeElectronDensity::getGeometry(
        std::vector<Vector3d>& positions)
        const
    {
        positions = mPositions;
    }

    double PromoleculeElectronDensity::calculate(
        const Vector3d &r) 
        const
    {
        int atomIdx, nAtoms = mPositions.size();
        double value=0.0;
        Vector3d rLocal;

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            rLocal = r - mPositions[atomIdx];
            value += mDensities[mAtomToDensityMap[atomIdx]].calculate(sqrt(rLocal*rLocal));
        }

        return value;
    }
/*
        std::vector<SphericalAtomicDensity> mDensities;
        std::vector<Vector3d> mPositions;
        std::vector<int> mAtomToDensityMap;
  */

}
