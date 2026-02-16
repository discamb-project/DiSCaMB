#include "discamb/Scattering/SfCalculator.h"
#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator2.h"
#include "discamb/Scattering/NumericalSphericalAtomFormFactor.h"
#include "discamb/Scattering/IamFormFactorCalculationsManager.h"

namespace discamb {

    class NumericalDensityIamFfCalculator: public AtomicFormFactorCalculationsManager {
        NumericalDensityIamFfCalculator();
        std::vector<NumericalSphericalAtomFormFactor> mElementFormFactor;
        std::shared_ptr<AnyScattererStructureFactorCalculator2> mSfCalculator;
        std::vector<int> mAtomicNumber;
        std::shared_ptr<IamFormFactorCalculationsManager> mNGaussianIamCalculator;
        ReciprocalLatticeUnitCell mReciprocalLatticeUnitCell;
    public:
        NumericalDensityIamFfCalculator(const Crystal& crystal, const nlohmann::json& data);
        ~NumericalDensityIamFfCalculator();
        virtual void update(const std::vector<AtomInCrystal>& atoms);
        virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i& hkl) const;
        virtual std::complex<double> calculateCart(int atomIdx, const Vector3d& hkl) const;

    };

}
