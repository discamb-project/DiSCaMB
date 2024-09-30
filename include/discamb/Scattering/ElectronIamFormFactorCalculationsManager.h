#ifndef _DISCAMB_SCATTERING_IAMFORMFACTORCALCULATIONMANAGER_H_
#define _DISCAMB_SCATTERING_IAMFORMFACTORCALCULATIONMANAGER_H_


#include "AtomicFormFactorCalculationsManager.h"
#include "discamb/Scattering/NGaussianFormFactor.h"


namespace discamb {

    class IamFormFactorCalculationsManager: public AtomicFormFactorCalculationsManager {
    public:
        IamFormFactorCalculationsManager(const Crystal &crystal);
        virtual ~IamFormFactorCalculationsManager();
        //virtual void setStructure(const Crystal &crystal);
        virtual void update(const std::vector<AtomInCrystal> &atoms);
        virtual std::complex<double> calculateFrac(size_t atomIdx, const Vector3i &hkl) const;
        virtual std::complex<double> calculateCart(size_t atomIdx, const Vector3d &hkl) const;
    private:
        IamFormFactorCalculationsManager();
        std::vector<NGaussianFormFactor> mFormFactors;
        std::vector<size_t> mAtomToFormFactorMap;
        ReciprocalLatticeUnitCell mReciprocalSpaceUnitCell;

    };
}

#endif