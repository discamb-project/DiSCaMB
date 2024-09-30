#ifndef _DISCAMB_SCATTERING_ATOMICFORMFACTORCALCULATIONMANAGER_H_
#define _DISCAMB_SCATTERING_ATOMICFORMFACTORCALCULATIONMANAGER_H_

#include <vector>
#include <complex>

#include "discamb/CrystalStructure/Crystal.h"

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class AtomicFormFactorCalculationsManager {
        public:
            virtual ~AtomicFormFactorCalculationsManager() {};
            virtual void update(const std::vector<AtomInCrystal> &atoms)=0;
            virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i &hkl) const=0;
            virtual std::complex<double> calculateCart(int atomIdx, const Vector3d &hkl) const=0;

            virtual void calculateFrac(
                const Vector3i &hkl,
                std::vector<std::complex<double> > &formFactors,
                const std::vector<bool> &includeAtom) const;

            virtual void calculateCart(
                const Vector3d &hkl, 
                std::vector<std::complex<double> > &formFactors,
                const std::vector<bool> &includeAtom) const;
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
                const std::vector <Vector3d> & hkl,
                std::vector < std::vector<std::complex<double> > >& formFactors,
                const std::vector<bool>& includeAtom) const;


    protected:
        static Vector3d convertMillerIndexToCartesian(const ReciprocalLatticeUnitCell &reciprocalUnitCell, const Vector3i &hkl);
    };

    /** @}*/

}

#endif
