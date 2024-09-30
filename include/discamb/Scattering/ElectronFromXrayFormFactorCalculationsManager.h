 #ifndef _DISCAMB_SCATTERING_ELECTRONFROMXRAYFORMFACTORCALCULATIONMANAGER_H_
#define _DISCAMB_SCATTERING_ELECTRONFROMXRAYFORMFACTORCALCULATIONMANAGER_H_


#include "AtomicFormFactorCalculationsManager.h"
#include "discamb/Scattering/NGaussianFormFactor.h"
#include "discamb/CrystalStructure/ReciprocalLatticeUnitCell.h"


#include <memory>


namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class ElectronFromXrayFormFactorCalculationManager: public AtomicFormFactorCalculationsManager {
    public:
        ElectronFromXrayFormFactorCalculationManager(const UnitCell &unitCell,
                                                     const std::vector<double> &nuclearCharge,
												     std::shared_ptr<AtomicFormFactorCalculationsManager> &manager);
        virtual ~ElectronFromXrayFormFactorCalculationManager();

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

        /**
        formFactors[hkl idx][atom idx]
        */
        virtual void calculateFrac(
            const std::vector<Vector3i>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;

        /**
        formFactors[hkl idx][atom idx]
        */
        virtual void calculateCart(
            const std::vector <Vector3d>& hkl,
            std::vector < std::vector<std::complex<double> > >& formFactors,
            const std::vector<bool>& includeAtom) const;



    private:
        ElectronFromXrayFormFactorCalculationManager();
		std::shared_ptr<AtomicFormFactorCalculationsManager> mManager;
        std::vector<double> mNuclearCharge;
        std::vector<std::complex<double> > mFormFactorsAtHkl000;
        void setFormFactorsAtHkl000();
        ReciprocalLatticeUnitCell mRUnitCell;
		mutable std::vector<std::complex<double> > mFormFactor_x;
    };
    /** @}*/
}

#endif