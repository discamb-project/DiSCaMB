#ifndef _DISCAMB_SCATTERING_IAMFORMFACTORCALCULATIONMANAGER_H_
#define _DISCAMB_SCATTERING_IAMFORMFACTORCALCULATIONMANAGER_H_


#include "AtomicFormFactorCalculationsManager.h"
#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/Scattering/NGaussianFormFactorsTable.h"

#include <map>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class IamFormFactorCalculationsManager: public AtomicFormFactorCalculationsManager {
    public:
        IamFormFactorCalculationsManager(const Crystal &crystal, const std::string &type = std::string("IT92"));
        IamFormFactorCalculationsManager(const Crystal &crystal, const std::map<std::string, NGaussianFormFactor> &formFactors);
        virtual ~IamFormFactorCalculationsManager();
        //virtual void setStructure(const Crystal &crystal);
        virtual void update(const std::vector<AtomInCrystal> &atoms);
        virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i &hkl) const;
        virtual std::complex<double> calculateCart(int atomIdx, const Vector3d &hkl) const;

        virtual void calculateFrac(
            const Vector3i& hkl,
            std::vector<std::complex<double> >& formFactors,
            const std::vector<bool>& includeAtom) const;

        virtual void calculateCart(
            const Vector3d& hkl,
            std::vector<std::complex<double> >& formFactors,
            const std::vector<bool>& includeAtom) const;


        virtual void calculateFrac(
            const std::vector<Vector3d>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;

    private:
        IamFormFactorCalculationsManager();
        std::vector<NGaussianFormFactor> mFormFactors;
        mutable std::vector<double> mFormFactorValues;
        std::vector<int> mAtomToFormFactorMap;
        ReciprocalLatticeUnitCell mReciprocalSpaceUnitCell;

    };
    /** @}*/
}

#endif