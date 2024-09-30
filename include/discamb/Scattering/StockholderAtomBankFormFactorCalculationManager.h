#pragma once

#include "AtomicFormFactorCalculationsManager.h"
#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/Scattering/HirshfeldAtomModelSettings.h"

namespace discamb {

class StockholderAtomBankFormFactorCalculationManager : public AtomicFormFactorCalculationsManager {
public:
    StockholderAtomBankFormFactorCalculationManager(
        const Crystal& crystal,
        const std::string &bankFileName, 
        int nThreadas, 
        bool useSphericalHarmonicsExpansion);

    ~StockholderAtomBankFormFactorCalculationManager();

    virtual void update(const std::vector<AtomInCrystal>& atoms);
    virtual std::complex<double> calculateFrac(int atomIdx, const Vector3i& hkl) const;
    virtual std::complex<double> calculateCart(int atomIdx, const Vector3d& hkl) const;

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


    //static void readAtomTypeElectronDensity(
    //    const std::string &fName, 
    //    std::vector<double> &electronDensity,
    //    int & angularGridSize,
    //    int & radialGridSize,
    //    int & atomicNumber);

    static void getGrid(
        int angularGridSize,
        int radialGridSize,
        int atomicNumber,
        std::vector<Vector3d> &gridPoints,
        std::vector<double> &gridWeights);


private:
    int mN_Threads = 1;
    Crystal mCrystal;
    
    ReciprocalLatticeUnitCell mReciprocalLatticeUnitCell;
    std::vector<int> mTypeIdx;
    std::vector<int> mAtomicNumber;
    std::vector<LocalCoordinateSystem<AtomInCrystalID> > mLcs;
    std::vector<std::vector<double> > mTypeElectronDensities;
    std::vector<std::vector<Vector3d> > mElementIntegationGrid;
    std::vector<std::vector<double> > mElementGridWeight;
    std::vector<LocalCoordinateSystemCalculator> mLcsCalculators;
    ham_settings::IntegrationGrid mIntegrationGrid;


    bool mUseSphericalHarmonicsExpansion = true;
    //[atom idx][l][l+m][radialIdx]
    std::vector < std::vector <std::vector < std::vector <double> > > > mSphericalHarmonicsCoefficients;
    ////[atom idx][l][l+m][pointIdx]
    //std::vector<std::vector<std::vector<std::vector<double> > > > mRadialFunctions;
    //[element][idx] 
    std::vector < std::vector <double> > mRadialIntegrationWeights;
    // precalculated mRadialIntegrationWeights[z][i]*mRadialGrid[z][i]*mRadialGrid[z][i]
    std::vector < std::vector <double> > mR2RW;
    std::vector < std::vector <double> > mRadialGrid;
    
    // preallocated vectors for use in calculations
    
    mutable std::vector <double> mRadialIntegrationMultiplier;
    mutable std::vector<std::vector <double> > mRadialIntegrationMultiplier_t;
    mutable std::map<int, std::vector< std::vector <double> > > mRadialIntegrationMultipliers;
    mutable std::vector < std::map<int, std::vector< std::vector <double> > > > mRadialIntegrationMultipliers_t;
    // a vector of size of radial grid
    mutable std::vector <double> mAuxRadial_1, mAuxRadial_2, mAuxRadial_3, mAuxRadial_4, mAuxRadial_5;
    mutable std::vector < std::vector <double> > mAuxRadial_t;
    mutable std::vector < std::vector <double> > mSphBessel;
    
    //--------
    std::set<int> mUniqueAtomicNumbers;
    std::vector <double> mAngularWeights;
    std::vector <Vector3d> mAngularGrid;
    std::complex<double> calculateCartSphExp(int atomIdx, const Vector3d& hkl) const;
    void calculateCartSphExp(const Vector3d& hkl, std::vector< std::complex<double> > &ff, const std::vector<bool>& includeAtom) const;
    void calculateCartSphExp_in_thread(const Vector3d& hkl, std::vector< std::complex<double> >& ff, const std::vector<bool>& includeAtom) const;
    void calculateCartSphExp_RecBessel(const Vector3d& hkl, std::vector< std::complex<double> >& ff, const std::vector<bool>& includeAtom) const;

    void setSphericalHarmonicsExpansionData();

    void calculateCart000(std::vector< std::complex<double> >& ff, const std::vector<bool>& includeAtom) const;
};

}
