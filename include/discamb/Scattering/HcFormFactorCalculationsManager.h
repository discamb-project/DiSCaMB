#ifndef _DISCAMB_SCATTERING_HCFORMFACTORCALCULATIONMANAGER_H_
#define _DISCAMB_SCATTERING_HCFORMFACTORCALCULATIONMANAGER_H_

#include "AtomicFormFactorCalculationsManager.h"
#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/Scattering/HansenCoppensStructureFactorCalculator.h"
//#include "discamb/HC_Model/XdLocalCoordinateSystem.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"

#include  <memory>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class HcFormFactorCalculationsManager: public AtomicFormFactorCalculationsManager
    {
        HcFormFactorCalculationsManager();
        HansenCoppensStructureFactorCalculator *mHcCalculator;
        //std::vector<XdLocalCoordinateSystem> mLcs;
		std::vector < std::shared_ptr <LocalCoordinateSystemInCrystal> > mLcs;
        std::vector<Matrix3d> mLcsMatrices;
        Crystal mCrystal, mAuxCrystal;
        ReciprocalLatticeUnitCell mReciprocalLatticeUnitCell;
		
    public:
		HcFormFactorCalculationsManager(
			const Crystal &crystal,
			const HC_ModelParameters &params,
			//const std::vector<XdLocalCoordinateSystem> &lcs);
			const std::vector < std::shared_ptr <LocalCoordinateSystemInCrystal> > &lcs);
        
        virtual ~HcFormFactorCalculationsManager();
        //virtual void setStructure(const Crystal &crystal);
        virtual void update(const std::vector<AtomInCrystal> &atoms);
        virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i &hkl) const;
        virtual std::complex<double> calculateCart(int atomIdx, const Vector3d &hkl) const;

		//void useNewImplementation(bool newImplementation);

		virtual void calculateFrac(
			const Vector3i& hkl,
			std::vector<std::complex<double> >& formFactors,
			const std::vector<bool>& includeAtom) const;

		virtual void calculateCart(
			const Vector3d& hkl,
			std::vector<std::complex<double> >& formFactors,
			const std::vector<bool>& includeAtom) const;

    };
    /** @}*/
}

#endif

