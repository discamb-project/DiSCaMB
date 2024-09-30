#pragma once


#include "discamb/Scattering/AtomicFormFactorCalculationsManager.h"

#include <map>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class ConstFormFactorCalculationsManager: public AtomicFormFactorCalculationsManager
	{
		std::map<Vector3i, std::vector<std::complex<double> > > mFormFactors;
		ReciprocalLatticeUnitCell mReciprocalLatticeUnitCell;
    
	public:
			ConstFormFactorCalculationsManager(const UnitCell &uc, std::map<Vector3i, std::vector<std::complex<double> > > &formFactors);
            virtual ~ConstFormFactorCalculationsManager();
            virtual void update(const std::vector<AtomInCrystal> &atoms);
            virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i &hkl) const;
            virtual std::complex<double> calculateCart(int atomIdx, const Vector3d &hkl) const;

            virtual void calculateFrac(
                const Vector3i &hkl,
                std::vector<std::complex<double> > &formFactors,
                const std::vector<bool> &includeAtom) const;

            virtual void calculateCart(
                const Vector3d &hkl, 
                std::vector<std::complex<double> > &formFactors,
                const std::vector<bool> &includeAtom) const;
    };
    /** @}*/
}


