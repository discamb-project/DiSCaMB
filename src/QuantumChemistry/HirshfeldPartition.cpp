#include "discamb/QuantumChemistry/HirshfeldPartition.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/IO/wfn_io.h"

#include <cmath>

using namespace std;

namespace discamb {

    HirshfeldPartition::HirshfeldPartition()
    {}

    HirshfeldPartition::HirshfeldPartition(
        std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
        std::shared_ptr<PromoleculeElectronDensity> &promoleculeDensity)
    {

        mEdCalculator = ed_calculator;
        mPromoleculeDensity = promoleculeDensity;
        mN_Atoms = promoleculeDensity->nAtoms();
        promoleculeDensity->getGeometry(mAtomPositions);
        setIncludeRange(-1);
    }

    HirshfeldPartition::HirshfeldPartition(
        const std::string& wfnFile,
        const std::string& atomicDensitiesFile,
        const std::string& method,
        const std::string& basisSet)
    {
        setFromFile(wfnFile, atomicDensitiesFile, method, basisSet);
    }

    



    HirshfeldPartition::HirshfeldPartition(
        const std::string& wfnFile,
        const std::string& method,
        const std::string& basisSet,
        const std::vector<std::string>& basis_sets,
        const std::map<int, int>& atomicNumberToBasisSetMap,
        const std::map<int, int>& atomToBasisSetMap,
        bool calculateSphericalDensityIfMissing,
        const std::string& qmProgram,
        const std::string& qmProgramFolder,
        const std::string& relativisticMethod,
        int nCores,
        int memoryMB)
    {
        if (!setFromDefaultFiles(wfnFile, method, basisSet, basis_sets, atomicNumberToBasisSetMap, atomToBasisSetMap, calculateSphericalDensityIfMissing,
            qmProgram, qmProgramFolder, relativisticMethod,nCores, memoryMB))
            on_error::throwException("cannot find matching discamb data files for spherical atoms when creating Hirshfeld partition", __FILE__, __LINE__);
    }



    HirshfeldPartition::~HirshfeldPartition()
    {
    }

    void HirshfeldPartition::set(
        const nlohmann::json& data,
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& atomPositions,
        std::shared_ptr<ElectronDensityCalculator>& edCalculator)
    {
        shared_ptr<PromoleculeElectronDensity> promoleculeDensity = shared_ptr<PromoleculeElectronDensity>(new PromoleculeElectronDensity());
        promoleculeDensity->set(atomicNumbers, atomPositions, data);

        set(edCalculator, promoleculeDensity);
        
        if(data.find("power")!=data.end())
            setPower(data.find("power")->get<double>());
        if(mIncludeRange>0.0)
            setIncludeRange(mIncludeRange);
    }

    void HirshfeldPartition::applySettings(
        const nlohmann::json& data)
    {

    }

    void HirshfeldPartition::setIncludeRange(double range)
    {
        mIncludeRange = range;
        

        int i, j, nAtoms = mAtomPositions.size();
        if(nAtoms==0)
            return;
        mIncludeAtoms.clear();
        mIncludeAtoms.resize(nAtoms);

        if (range < 0)
        {
            for (i = 0; i < nAtoms; i++)
                mIncludeAtoms[0].push_back(i);
            for (i = 1; i < nAtoms; i++)
                mIncludeAtoms[i] = mIncludeAtoms[0];
            return;
        }

        double length2, range2 = range*range;
        Vector3d r;
        for (i = 0; i < nAtoms; i++)
        {
            mIncludeAtoms[i].push_back(i);
            for (j = 0; j < i; j++)
            {
                r = mAtomPositions[i] - mAtomPositions[j];
                length2 = r * r;
                if (length2 <= range2)
                {
                    mIncludeAtoms[i].push_back(j);
                    mIncludeAtoms[j].push_back(i);
                }
            }
        }

        mEdCalculator->setContributingCenters(mIncludeAtoms[mAtomIdx]);
    }
    
    void HirshfeldPartition::unsetIncludeRange()
    {
        mIncludeRange = -1;
    }

    void HirshfeldPartition::setEdCalculatorFromWfnFile(
        const std::string& wfnFile,
        std::vector<Vector3d>& positions,
        std::vector<int>& atomicNumbers)
    {
        wfn_io::WfnFileData wfn;
        wfn_io::read_wfx(wfnFile, wfn);

        positions = wfn.center_position;
        atomicNumbers = wfn.atomic_numbers;

        mEdCalculator = std::shared_ptr<ElectronDensityCalculator>(new ElectronDensityCalculator);

        mEdCalculator->set_1rdm(
            wfn.center_position,
            wfn.primitive_to_center,
            wfn.primitive_type,
            wfn.primitive_exponents,
            wfn.molecular_orbital_occupancy,
            wfn.molecular_orbitals);

        wfn_io::add_edf_from_library(wfn);

        if (!wfn.edfs.empty())
            mEdCalculator->setAdditionalDensity(
                wfn.edfs[0].primitive_to_center,
                wfn.edfs[0].primitive_type,
                wfn.edfs[0].primitive_exponents,
                wfn.edfs[0].primitive_coefficients);

    }

    void HirshfeldPartition::setFromFile(
        const std::string& wfnFile,
        const std::string& _atomicDensitiesFile,
        const std::string& method,
        const std::string& basisSet)
    {
        //std::shared_ptr<ElectronDensityCalculator> ed_calculator;
        //std::shared_ptr<PromoleculeElectronDensity> promoleculeDensity;
        string atomicDensitiesFile = _atomicDensitiesFile;
        //wfn_io::WfnFileData wfn;
        //wfn_io::read_wfx(wfnFile, wfn);
        //mN_Atoms = wfn.center_position.size();
        //mEdCalculator = std::shared_ptr<ElectronDensityCalculator>(new ElectronDensityCalculator);
        //
        //mEdCalculator->set_1rdm(
        //    wfn.center_position,
        //    wfn.primitive_to_center,
        //    wfn.primitive_type,
        //    wfn.primitive_exponents,
        //    wfn.molecular_orbital_occupancy,
        //    wfn.molecular_orbitals);

        //if (!wfn.edfs.empty())
        //    mEdCalculator->setAdditionalDensity(
        //        wfn.edfs[0].primitive_to_center,
        //        wfn.edfs[0].primitive_type,
        //        wfn.edfs[0].primitive_exponents,
        //        wfn.edfs[0].primitive_coefficients);

        std::vector<Vector3d> positions;
        std::vector<int> atomicNumbers;
        setEdCalculatorFromWfnFile(wfnFile, positions, atomicNumbers);
        mN_Atoms = atomicNumbers.size();
        ProatomDB proatomDB;
        
        if (atomicDensitiesFile.empty())
            if (!ProatomDB::hasDiscambDataFile(method, basisSet, atomicDensitiesFile))
                on_error::throwException("atomic densities file not given neither available from DiSCaMB data",__FILE__,__LINE__);
        
        proatomDB.setFromFile(atomicDensitiesFile);

        //hasDiscambDataFile(mMethod, mBase, mAtomsFile)
        mPromoleculeDensity = shared_ptr<PromoleculeElectronDensity>(
            new PromoleculeElectronDensity(proatomDB, atomicNumbers, positions));
        mPromoleculeDensity->getGeometry(mAtomPositions);
        setIncludeRange(-1);
    }

    bool HirshfeldPartition::setFromDefaultFiles(
        const std::string& wfnFile,
        const std::string& method,
        const std::string& basisSet,
        const std::vector<std::string>& basis_sets,
        const std::map<int, int>& atomicNumberToBasisSetMap,
        const std::map<int, int>& atomToBasisSetMap,
        bool calculateSphericalDensityIfMissing,
        const std::string& qmProgram,
        const std::string& qmProgramFolder,
        const std::string& relativisticMethod,
        int nCores,
        int memoryMB)
    {

        std::vector<Vector3d> positions;
        std::vector<int> atomicNumbers;
        setEdCalculatorFromWfnFile(wfnFile, positions, atomicNumbers);
        mN_Atoms = atomicNumbers.size();

        ProatomDB proatomDB;

        mPromoleculeDensity = shared_ptr<PromoleculeElectronDensity>(new PromoleculeElectronDensity());

        if (!mPromoleculeDensity->setFromDefaultFiles(
            atomicNumbers,
            positions,
            method,
            basisSet,
            basis_sets,
            atomicNumberToBasisSetMap,
            atomToBasisSetMap,
            true,
            qmProgram,
            qmProgramFolder,
            relativisticMethod,nCores,memoryMB))
            return false;
        mPromoleculeDensity->getGeometry(mAtomPositions);
        setIncludeRange(-1);
        return true;
    }


    void HirshfeldPartition::setPower(double p)
    {
        mPower = p;
    }

    HirshfeldPartition* HirshfeldPartition::clone()
        const
    {
        HirshfeldPartition *hirshfeldPartition;
        shared_ptr<ElectronDensityCalculator> ed_calculator(mEdCalculator->clone());
        shared_ptr<PromoleculeElectronDensity> promoleculeDensity(mPromoleculeDensity->clone());
                
        hirshfeldPartition = new HirshfeldPartition();
        hirshfeldPartition->set(ed_calculator, promoleculeDensity);
        hirshfeldPartition->setPower(mPower);
        hirshfeldPartition->setIncludeRange(mIncludeRange);
        hirshfeldPartition->setAtom(mAtomIdx);

        return hirshfeldPartition;
    }

    void HirshfeldPartition::set(
        std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
        std::shared_ptr<PromoleculeElectronDensity> &promoleculeDensity)
    {
        mEdCalculator = ed_calculator;
        mPromoleculeDensity = promoleculeDensity;
        mN_Atoms = promoleculeDensity->nAtoms();
        promoleculeDensity->getGeometry(mAtomPositions);
        setIncludeRange(mIncludeRange);
    }

    void HirshfeldPartition::setAtom(int atomIdx)
    { 
        mAtomIdx = atomIdx; 
        mEdCalculator->setContributingCenters(mIncludeAtoms[atomIdx]);
    };

    double HirshfeldPartition::calculateDifferenceDensity(
        int atom,
        const Vector3d& r)
        const
    {
        double atomInMolDensity = calculate(atom, r);
        double isolatedAtomDensity = mPromoleculeDensity->calculate(r, atom);
        return atomInMolDensity - isolatedAtomDensity;
    }

    double HirshfeldPartition::calculate(
        int atom,
        const Vector3d &r)
        const 
    {
        //mEdCalculator->setContributingCenters(mIncludeAtoms[atom]);

        if (mPower == 1.0)
        {
            double promol = mPromoleculeDensity->calculate(r);
            if (promol > 0)
            {
                double rho = mEdCalculator->calculate2(r.x, r.y, r.z);
                double promol_atom = mPromoleculeDensity->calculate(r, atom);
                return rho * promol_atom / promol;
            }
        }
        else
        {
            int nAtoms = mPromoleculeDensity->nAtoms();
            double weight_denominator = 0;
            double atomicDensity;
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                atomicDensity = mPromoleculeDensity->calculate(r, atomIdx);
                if(atomicDensity>0)
                    weight_denominator += pow(atomicDensity , mPower);
            }
        
            if (weight_denominator > 0)
            {
                double rho = mEdCalculator->calculate2(r.x, r.y, r.z);
                atomicDensity = mPromoleculeDensity->calculate(r, atom);
                if (atomicDensity > 0)
                {
                    double weight_numerator = pow(mPromoleculeDensity->calculate(r, atom), mPower);
                    return rho * weight_numerator / weight_denominator;
                }
            }
        }
            return 0.0;
    }

    double HirshfeldPartition::calculate(
        int atom,
        const Vector3d& r,
        int molecularOrbitalIdx)
        const
    {
        double promol = mPromoleculeDensity->calculate(r);
        if (promol > 0)
        {
            double rho = mEdCalculator->calculate2(r.x, r.y, r.z, molecularOrbitalIdx);
            double promol_atom = mPromoleculeDensity->calculate(r, atom);
            return rho * promol_atom / promol;
        }

        return 0.0;
    }

    
}
