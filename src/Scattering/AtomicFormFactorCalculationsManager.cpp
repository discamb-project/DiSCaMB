#include "discamb/Scattering/AtomicFormFactorCalculationsManager.h"

using namespace std;

namespace discamb {

/*
virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i &hkl) const=0;
virtual std::complex<double> calculateCart(int atomIdx, const Vector3d &hkl) const=0;
*/


    void AtomicFormFactorCalculationsManager::calculateFrac(
        const Vector3i &hkl,
        std::vector<std::complex<double> > &formFactors,
        const std::vector<bool> &includeAtom)
        const
    {
        int i, n = includeAtom.size();
		formFactors.resize(n);
        for (i = 0; i < n; i++)
            if (includeAtom[i])
                formFactors[i] = calculateFrac(i, hkl);
    }

    void AtomicFormFactorCalculationsManager::calculateCart(
        const Vector3d &hkl,
        std::vector<std::complex<double> > &formFactors,
        const std::vector<bool> &includeAtom)
        const
    {
        int i, n = formFactors.size();
        for (i = 0; i < n; i++)
            if (includeAtom[i])
                formFactors[i] = calculateCart(i, hkl);
    }

    void AtomicFormFactorCalculationsManager::calculateFrac(
        const std::vector<Vector3i>& hkl,
        std::vector < std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom) const
    {
        int hklIdx, nHkl = hkl.size();
        vector<complex<double> > oneH_FormFactors;
        formFactors.clear();
        formFactors.resize(nHkl);

        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            calculateFrac(hkl[hklIdx], oneH_FormFactors, includeAtom);
            formFactors[hklIdx] = oneH_FormFactors;
        }

    }

    void AtomicFormFactorCalculationsManager::calculateCart(
        const std::vector <Vector3d>& hkl,
        std::vector < std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        int hklIdx, nHkl = hkl.size();
        vector<complex<double> > oneH_FormFactors;
        formFactors.clear();
        formFactors.resize(nHkl);

        for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            calculateCart(hkl[hklIdx], oneH_FormFactors, includeAtom);
            formFactors[hklIdx] = oneH_FormFactors;
        }

    }


    Vector3d AtomicFormFactorCalculationsManager::convertMillerIndexToCartesian(
        const ReciprocalLatticeUnitCell &reciprocalUnitCell,
        const Vector3i &hkl)
    {
        Vector3d result;
        reciprocalUnitCell.fractionalToCartesian(hkl, result);
        return result;
    }
    
 
}

