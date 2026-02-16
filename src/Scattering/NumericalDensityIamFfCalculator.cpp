#include "discamb/Scattering/NumericalDensityIamFfCalculator.h"

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/QuantumChemistry/ProatomDB.h"


using namespace std;

namespace discamb {

    NumericalDensityIamFfCalculator::NumericalDensityIamFfCalculator()
    {
    }

    NumericalDensityIamFfCalculator::NumericalDensityIamFfCalculator(
        const Crystal& crystal,
        const nlohmann::json& data)
    {
        string proatomDbFile = data.value("atomic densities", "");
        if (proatomDbFile.empty())
            on_error::throwException("atomic densities file not defined", __FILE__, __LINE__);

        mNGaussianIamCalculator = shared_ptr<IamFormFactorCalculationsManager>(
            new IamFormFactorCalculationsManager(crystal, "Waasmeier-Kirfel"));

        mReciprocalLatticeUnitCell.set(crystal.unitCell);
        mElementFormFactor.resize(121);

        ProatomDB db;
        db.setFromFile(proatomDbFile);
       
        std::set<int> uniqueZ;

        crystal_structure_utilities::atomicNumbers(crystal, mAtomicNumber);
        uniqueZ.insert(mAtomicNumber.begin(), mAtomicNumber.end());

        for (auto& z : uniqueZ)
        {
            SphericalAtomicDensity density;
            db.getSphericalAtom(z, 0, density);
            mElementFormFactor[z].set(density);
        }

    }

    NumericalDensityIamFfCalculator::~NumericalDensityIamFfCalculator()
    {
    }

    void NumericalDensityIamFfCalculator::update(
        const std::vector<AtomInCrystal>& atoms)
    {
    }

    std::complex<double> NumericalDensityIamFfCalculator::calculateFrac(
        int atomIdx,
        const Vector3i& hkl)
        const
    {
        Vector3d cart;
        mReciprocalLatticeUnitCell.fractionalToCartesian(hkl, cart);
        return calculateCart(atomIdx, cart);
    }

    std::complex<double> NumericalDensityIamFfCalculator::calculateCart(
        int atomIdx,
        const Vector3d& hkl)
        const
    {
        int z = mAtomicNumber[atomIdx];
        if (z == 1)
            return mNGaussianIamCalculator->calculateCart(atomIdx, hkl);
        return mElementFormFactor[mAtomicNumber[atomIdx]].calculate(hkl.norm());
    }


}
