#include "discamb/Scattering/StockholderAtomBankFormFactorCalculationManager.h"

#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/Constants.h"
#include "discamb/IO/atom_type_io.h"
#include "discamb/IO/tham_io.h"
#include "discamb/MathUtilities/RadialGrid.h"
#include "discamb/MathUtilities/LebedevGrid.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include "discamb/MathUtilities/RealSphericalHarmonics.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"



#include <set>
#include <cmath>
#include <fstream>
#include <omp.h>

using namespace std;

namespace discamb {

    StockholderAtomBankFormFactorCalculationManager::StockholderAtomBankFormFactorCalculationManager(
        const Crystal& crystal,
        const std::string& bankFileName,
        int nThreadas,
        bool useSphericalHarmonicsExpansion)
    {
        mUseSphericalHarmonicsExpansion = useSphericalHarmonicsExpansion;
        mN_Threads = nThreadas;
        mCrystal = crystal;
        crystal_structure_utilities::atomicNumbers(crystal, mAtomicNumber);

        mReciprocalLatticeUnitCell.set(crystal.unitCell);

        // read bank
        vector<AtomType> types;
        DescriptorsSettings descriptorsSettings;
        atom_type_io::readAtomTypes(bankFileName, types, descriptorsSettings);

        // assign atom types
        CrystalAtomTypeAssigner assigner;
        //vector<int> typeIdx;
        vector<LocalCoordinateSystem<AtomInCrystalID> > lcs;
        assigner.setAtomTypes(types);
        assigner.setDescriptorsSettings(descriptorsSettings);
        assigner.assign(crystal, mTypeIdx, lcs);

        if (*std::min_element(mTypeIdx.begin(), mTypeIdx.end()) < 0)
            on_error::throwException("THAM implementation works only if all atoms have atom type assigned, which is not the case", __FILE__, __LINE__);

        int atomIdx, nAtoms = mTypeIdx.size();
        mLcsCalculators.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            mLcsCalculators[atomIdx].set(lcs[atomIdx], crystal);

        //

        mTypeElectronDensities.resize(types.size());
        mElementIntegationGrid.resize(120);
        mElementGridWeight.resize(120);

        mRadialIntegrationWeights.resize(120);
        mR2RW.resize(120);
        mRadialGrid.resize(120);

        // + mutable std::vector<std::vector <double> > mRadialIntegrationMultiplier_t;
        // + mutable std::vector < std::map<int, std::vector< std::vector <double> > > > mRadialIntegrationMultipliers_t;
        //mutable std::vector < std::vector <double> > mAuxRadial_t;


        // read atom type electron densities

        std::set<int> uniqueTypeIds(mTypeIdx.begin(), mTypeIdx.end());

        

        for (int uniqueTypeIdx : uniqueTypeIds)
        {
            string fileName = "type_" + types[uniqueTypeIdx].id;
        
            int angularGridSize, radialGridSize, atomicNumber;
            tham_io::readAtomTypeElectronDensity(fileName, mTypeElectronDensities[uniqueTypeIdx],
                mIntegrationGrid.angularGridSize, mIntegrationGrid.radialGridSize, atomicNumber);
        }

        //

        if (mUseSphericalHarmonicsExpansion)
        {

            mRadialIntegrationMultiplier.resize(mIntegrationGrid.radialGridSize);
            mRadialIntegrationMultiplier_t.resize(mN_Threads, vector<double>(mIntegrationGrid.radialGridSize));
            mUniqueAtomicNumbers.insert(mAtomicNumber.begin(), mAtomicNumber.end());
            for (int z : mUniqueAtomicNumbers)
                mRadialIntegrationMultipliers[z].assign(8, vector<double>(mIntegrationGrid.radialGridSize));
            mRadialIntegrationMultipliers_t.resize(mN_Threads, mRadialIntegrationMultipliers);

            mAuxRadial_1.resize(mIntegrationGrid.radialGridSize);
            mAuxRadial_2.resize(mIntegrationGrid.radialGridSize);
            mAuxRadial_3.resize(mIntegrationGrid.radialGridSize);
            mAuxRadial_4.resize(mIntegrationGrid.radialGridSize);
            mAuxRadial_5.resize(mIntegrationGrid.radialGridSize);
            mAuxRadial_t.resize(mN_Threads, vector<double>(mIntegrationGrid.radialGridSize));
            //

            lebedev_laikov::get_grid(mIntegrationGrid.angularGridSize, mAngularGrid, mAngularWeights);

            std::set<int> uniqueAtomicNumbers(mAtomicNumber.begin(), mAtomicNumber.end());
            for (int z : uniqueAtomicNumbers)
            {
                radial_grid::mura_knowles(z, mIntegrationGrid.radialGridSize,
                    mRadialGrid[z], mRadialIntegrationWeights[z]);


                int gridSize = mIntegrationGrid.angularGridSize * mIntegrationGrid.radialGridSize;
                mElementGridWeight[z].resize(gridSize);
                mElementIntegationGrid[z].resize(gridSize);
                mR2RW[z].resize(mIntegrationGrid.radialGridSize);

                int idx = 0;

                for (int angularIdx = 0; angularIdx < mIntegrationGrid.angularGridSize; angularIdx++)
                    for (int radialIdx = 0; radialIdx < mIntegrationGrid.radialGridSize; radialIdx++)
                    {
                        mElementIntegationGrid[z][idx] = mRadialGrid[z][radialIdx] * mAngularGrid[angularIdx];
                        mElementGridWeight[z][idx] = mRadialIntegrationWeights[z][radialIdx] * mAngularWeights[angularIdx] * 4 * M_PI * mRadialGrid[z][radialIdx] * mRadialGrid[z][radialIdx];
                        idx++;
                    }
                for (int radialIdx = 0; radialIdx < mIntegrationGrid.radialGridSize; radialIdx++)
                    mR2RW[z][radialIdx] = mRadialIntegrationWeights[z][radialIdx] * mRadialGrid[z][radialIdx] * mRadialGrid[z][radialIdx];
            }
            setSphericalHarmonicsExpansionData();
        }
        else
        {

            lebedev_laikov::get_grid(mIntegrationGrid.angularGridSize, mAngularGrid, mAngularWeights);

            std::set<int> uniqueAtomicNumbers(mAtomicNumber.begin(), mAtomicNumber.end());
            for (int z : uniqueAtomicNumbers)
            {
                radial_grid::mura_knowles(z, mIntegrationGrid.radialGridSize,
                    mRadialGrid[z], mRadialIntegrationWeights[z]);


                int gridSize = mIntegrationGrid.angularGridSize * mIntegrationGrid.radialGridSize;
                mElementGridWeight[z].resize(gridSize);
                mElementIntegationGrid[z].resize(gridSize);
                mR2RW[z].resize(mIntegrationGrid.radialGridSize);

                int idx = 0;

                for (int angularIdx = 0; angularIdx < mIntegrationGrid.angularGridSize; angularIdx++)
                    for (int radialIdx = 0; radialIdx < mIntegrationGrid.radialGridSize; radialIdx++)
                    {
                        mElementIntegationGrid[z][idx] = mRadialGrid[z][radialIdx] * mAngularGrid[angularIdx];
                        mElementGridWeight[z][idx] = mRadialIntegrationWeights[z][radialIdx] * mAngularWeights[angularIdx] * 4 * M_PI * mRadialGrid[z][radialIdx] * mRadialGrid[z][radialIdx];
                        idx++;
                    }
                for (int radialIdx = 0; radialIdx < mIntegrationGrid.radialGridSize; radialIdx++)
                    mR2RW[z][radialIdx] = mRadialIntegrationWeights[z][radialIdx] * mRadialGrid[z][radialIdx] * mRadialGrid[z][radialIdx];
            }
            

        }
    }

    StockholderAtomBankFormFactorCalculationManager::~StockholderAtomBankFormFactorCalculationManager()
    {
    }

    void StockholderAtomBankFormFactorCalculationManager::setSphericalHarmonicsExpansionData()
    {
        //[gridIdx][l][m]
        vector<vector<vector<double> > > sphericalHarmonics(mIntegrationGrid.angularGridSize);
        int atomIdx, nAtoms = mAtomicNumber.size();
        mSphericalHarmonicsCoefficients.resize(nAtoms,vector<vector<vector<double> > >(8));

        int l, m;

        for (auto& v : sphericalHarmonics)
        {
            v.resize(8);
            for (l = 0; l < 8; l++)
                v[l].resize(2 * l + 1);
        }

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            for (l = 0; l < 8; l++)
            {
                mSphericalHarmonicsCoefficients[atomIdx][l].resize(2 * l + 1);
                for (m = -l; m <= l; m++)
                    mSphericalHarmonicsCoefficients[atomIdx][l][l+m].resize(mIntegrationGrid.radialGridSize);
            }

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            Matrix3d lcs;
            mLcsCalculators[atomIdx].calculate(lcs, mCrystal);
            int radialIdx, angularIdx;
            int atomTypeIdx = mTypeIdx[atomIdx];

            for (angularIdx = 0; angularIdx < mIntegrationGrid.angularGridSize; angularIdx++)
            {
                //angularGrid[i] = lcs * mAngularGrid[i];
                Vector3d r = lcs * mAngularGrid[angularIdx];
                real_spherical_harmonics::getWfnNormalized<7>(r, sphericalHarmonics[angularIdx]);
            }
           
            for (l = 0; l < 8; l++)
                for (m = -l; m <= l; m++)
                    for (int radialIdx = 0; radialIdx < mIntegrationGrid.radialGridSize; radialIdx++)
                    {
                        double sum = 0;
                        for (int angularIdx = 0; angularIdx < mIntegrationGrid.angularGridSize; angularIdx++)
                        {
                            double electronDensity = mTypeElectronDensities[atomTypeIdx][mIntegrationGrid.radialGridSize * angularIdx + radialIdx];
                            sum += electronDensity * mAngularWeights[angularIdx] * sphericalHarmonics[angularIdx][l][l+m];
                        }
                        
                        mSphericalHarmonicsCoefficients[atomIdx][l][l + m][radialIdx] = 4 * M_PI * sum;
                    }

        }
    }

    void StockholderAtomBankFormFactorCalculationManager::getGrid(
        int angularGridSize,
        int radialGridSize,
        int atomicNumber,
        std::vector<Vector3d>& gridPoints,
        std::vector<double>& gridWeights)
    {
        vector<double> radialGrid, weightsRadial, weightsAngular;
        std::vector<Vector3d> angularGrid;
        radial_grid::mura_knowles(atomicNumber, radialGridSize, radialGrid, weightsRadial);
        lebedev_laikov::get_grid(angularGridSize, angularGrid, weightsAngular);
        
        int gridSize = angularGridSize * radialGridSize;
        gridWeights.resize(gridSize);
        gridPoints.resize(gridSize);

        int idx = 0;

        for (int angularIdx = 0; angularIdx < angularGridSize; angularIdx++)
            for (int radialIdx = 0; radialIdx < radialGridSize; radialIdx++)
            {
                gridPoints[idx] = radialGrid[radialIdx] * angularGrid[angularIdx];
                gridWeights[idx] = weightsRadial[radialIdx] * weightsAngular[angularIdx] * 4 * M_PI * radialGrid[radialIdx] * radialGrid[radialIdx];
                idx++;
            }

    }


    //void StockholderAtomBankFormFactorCalculationManager::readAtomTypeElectronDensity(
    //    const std::string& fName,
    //    std::vector<double>& electronDensity,
    //    int& angularGridSize,
    //    int& radialGridSize,
    //    int& atomicNumber)
    //{
    //    ifstream in(fName);
    //    if (!in.good())
    //        on_error::throwException("cannot read file '" + fName + "'", __FILE__, __LINE__);
    //    string s;
    //    in >> s >> s >> atomicNumber >> s >> s >> angularGridSize >> s >> s >> radialGridSize;
    //    int gridSize = angularGridSize * radialGridSize;
    //    electronDensity.clear();
    //    electronDensity.resize(gridSize);

    //    for (int i = 0; i < gridSize; i++)
    //        in >> electronDensity[i];

    //    in.close();
    //}


    void StockholderAtomBankFormFactorCalculationManager::update(const std::vector<AtomInCrystal>& atoms)
    {
        mCrystal.atoms = atoms;
    }

    std::complex<double> StockholderAtomBankFormFactorCalculationManager::calculateFrac(int atomIdx, const Vector3i& hkl) 
        const 
    {
        Vector3d h_cart;
        mReciprocalLatticeUnitCell.fractionalToCartesian(hkl, h_cart);
        return calculateCart(atomIdx, h_cart);
    }
    
    std::complex<double> StockholderAtomBankFormFactorCalculationManager::calculateCart(
        int atomIdx, 
        const Vector3d& hkl) 
        const
    {
        //return calculateCartSphExp(atomIdx, hkl);

        Matrix3d lcs, lcs_t;
        mLcsCalculators[atomIdx].calculate(lcs, mCrystal);
        lcs_t = lcs;
        lcs_t.transpose();
        Vector3d h_rot = lcs_t * hkl / constants::Angstrom;

        double real = 0;
        double imaginary = 0;
        const vector<Vector3d>& grid = mElementIntegationGrid[mAtomicNumber[atomIdx]];
        const vector<double> &weights = mElementGridWeight[mAtomicNumber[atomIdx]];
        const vector<double>& densities = mTypeElectronDensities[mTypeIdx[atomIdx]];
        double density;
        double m = 2.0 * M_PI;
        for (int i = 0; i < grid.size(); i++)
        {
            //density = exp(-2.0 * grid[i] * grid[i]);
            density = densities[i];
            double d = m * h_rot * grid[i];
            real += cos(d) * weights[i] * density;
            imaginary += sin(d) * weights[i] * density;
            //real += cos(d) * weights[i] * densities[i];
            //imaginary += sin(d) * weights[i] * densities[i];
        }

        //for (int i = 0; i < grid.size(); i++)
        //{
        //    //density = exp(-2.0 * grid[i] * grid[i]);
        //    density = densities[i];
        //    Vector3d r = lcs * grid[i];
        //    double d = 2.0 * M_PI * hkl * r / constants::Angstrom;
        //    real += cos(d) * weights[i] * density;
        //    imaginary += sin(d) * weights[i] * density;
        //    //real += cos(d) * weights[i] * densities[i];
        //    //imaginary += sin(d) * weights[i] * densities[i];
        //}


        return complex<double>(real, imaginary);
    }


    std::complex<double> StockholderAtomBankFormFactorCalculationManager::calculateCartSphExp(
        int atomIdx,
        const Vector3d& hkl)
        const
    {

        vector<vector<double> > sphericalHarmonics(8);
        int l, m;

        for (l = 0; l < 8; l++)
            sphericalHarmonics[l].resize(2 * l + 1);
        double h_length = sqrt(hkl * hkl);
        double h_length_atomic_units = h_length / constants::Angstrom;
        Vector3d h_normalized(1, 0, 0);

        if (h_length > 0)
            h_normalized = hkl / h_length;

        real_spherical_harmonics::getWfnNormalized<7>(h_normalized, sphericalHarmonics);

        vector<double> iPowL = { 1, 1, -1, -1, 1, 1, -1, -1 };

        double realPart = 0;
        double imagPart = 0;
        for (l = 0; l < 8; l++)
        {

            int z = mAtomicNumber[atomIdx];
            int i, n = mRadialGrid[z].size();

            for (i = 0; i < n; i++)
                mRadialIntegrationMultiplier[i] = 
                sph_bessel(l, h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][i]) *
                mRadialGrid[z][i] * mRadialGrid[z][i] * mRadialIntegrationWeights[mAtomicNumber[atomIdx]][i];

                //sphericalBesselTransform += sph_bessel(l, h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][i]) *
                //mRadialGrid[z][i] * mRadialGrid[z][i] * mSphericalHarmonicsCoefficients[atomIdx][l][l + m][i] *
                //mRadialIntegrationWeights[mAtomicNumber[atomIdx]][i];

            for (m = -l; m <= l; m++)
            {
                double sphericalBesselTransform = 0;
                for (i = 0; i < n; i++)
                    sphericalBesselTransform += mRadialIntegrationMultiplier[i] * mSphericalHarmonicsCoefficients[atomIdx][l][l + m][i];
                    

                if (l % 2 == 0)
                    realPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
                else
                    imagPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
            }
        }
        return complex<double>(realPart, imagPart);
    }

    void StockholderAtomBankFormFactorCalculationManager::calculateFrac(
        const Vector3i& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        Vector3d cart;
        mReciprocalLatticeUnitCell.fractionalToCartesian(hkl, cart);
        calculateCart(cart, formFactors, includeAtom);
    }

    void StockholderAtomBankFormFactorCalculationManager::calculateCart(
        const Vector3d& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom) 
        const
    {
        calculateCartSphExp(hkl, formFactors, includeAtom);
    }



void StockholderAtomBankFormFactorCalculationManager::calculateCart000(
    std::vector< std::complex<double> >& ff, 
    const std::vector<bool>& includeAtom) 
    const
{
    int atomIdx, nAtoms = mAtomicNumber.size();
    ff.assign(nAtoms, 1.0);
}



void StockholderAtomBankFormFactorCalculationManager::calculateCartSphExp_RecBessel(
    const Vector3d& hkl,
    std::vector< std::complex<double> >& ff,
    const std::vector<bool>& includeAtom)
    const
{
    double h_length = sqrt(hkl * hkl);
    if (h_length == 0.0)
        return calculateCart000(ff, includeAtom);
    //cout << "*\n";
    int atomIdx, nAtoms = mAtomicNumber.size();
    ff.resize(nAtoms);
    vector<vector<double> > sphericalHarmonics(8);
    int l, m;

    for (l = 0; l < 8; l++)
        sphericalHarmonics[l].resize(2 * l + 1);

    double h_length_atomic_units = h_length / constants::Angstrom;
    Vector3d h_normalized(1, 0, 0);

    if (h_length > 0)
        h_normalized = hkl / h_length;

    real_spherical_harmonics::getWfnNormalized<7>(h_normalized, sphericalHarmonics);

    vector<double> iPowL = { 1, 1, -1, -1, 1, 1, -1, -1 };

    int gridPointIdx, gridSize = mIntegrationGrid.radialGridSize;

    vector<vector<double> > sph_bess_rec(8, vector<double>(gridSize));
    vector<vector<double> > sph_bess_reg(8, vector<double>(gridSize));
    vector<vector<double> > diff(8, vector<double>(gridSize));
    for (int z : mUniqueAtomicNumbers)
    {
        vector<vector<double> >& v = mRadialIntegrationMultipliers[z];
        vector<double>& two_pi_h_r = mAuxRadial_1;
        vector<double>& sph_bessel_7 = mAuxRadial_2;
        vector<double>& sph_bessel_6 = mAuxRadial_3;

        for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
        {
            double x = h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][gridPointIdx];
            two_pi_h_r[gridPointIdx] = x;
            sph_bessel_7[gridPointIdx] = sph_bessel(7, x); // sin(x) / x;
            sph_bessel_6[gridPointIdx] = sph_bessel(6, x); //sph_bessel_0[gridPointIdx] / x - cos(x) / x;
        }

        vector<double>& sph_bessel_n_plus_2 = sph_bessel_7;
        vector<double>& sph_bessel_n_plus_1 = sph_bessel_6;
        
        const vector<double>& r2rw = mR2RW[z];


        // l=7
        {
            vector<double>& vv = v[7];

            for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
                vv[gridPointIdx] = sph_bessel_7[gridPointIdx] * r2rw[gridPointIdx];


        }

        // l=6
        {
            vector<double>& vv = v[6];

            for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
                vv[gridPointIdx] = sph_bessel_6[gridPointIdx] * r2rw[gridPointIdx];
        }


        for (l = 5; l >= 0; l--)
        {
            vector<double>& vv = v[l];
            double sb_l, x;
            for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
            {
                x = two_pi_h_r[gridPointIdx];
                sb_l = (2.0 * l + 3.0) / x * sph_bessel_n_plus_1[gridPointIdx] - sph_bessel_n_plus_2[gridPointIdx];
                vv[gridPointIdx] = sb_l * r2rw[gridPointIdx];
                sph_bess_rec[l][gridPointIdx] = sb_l;
                sph_bessel_n_plus_2[gridPointIdx] = sph_bessel_n_plus_1[gridPointIdx];
                sph_bessel_n_plus_1[gridPointIdx] = sb_l;
            }
        }


        //        for (l = 2; l < 8; l++)
        //        {
        //            vector<double>& vv = v[l];
        //            //double sb_l, x;
        //            for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
        //            {
        //            //    x = two_pi_h_r[gridPointIdx];
        //            //    sb_l = (2 * l - 1) / x * sph_bessel_n_minus_1[gridPointIdx] - sph_bessel_n_minus_2[gridPointIdx];
        //            //    //vv[gridPointIdx] = sph_bessel(l, two_pi_h_r[gridPointIdx]) * r2rw[gridPointIdx];
        //            //    vv[gridPointIdx] = sb_l * r2rw[gridPointIdx];
        //            //    sph_bessel_n_minus_2[gridPointIdx] = sph_bessel_n_minus_1[gridPointIdx];
        //            //    sph_bessel_n_minus_1[gridPointIdx] = sb_l;
        //                vv[gridPointIdx] = sph_bessel(l, two_pi_h_r[gridPointIdx]) * r2rw[gridPointIdx];
        //            }
        //            
        //            //vv[gridPointIdx] =  mRadialGrid[z][gridPointIdx] * mRadialGrid[z][gridPointIdx] * mRadialIntegrationWeights[z][gridPointIdx];
        ////vv[gridPointIdx] = 1.0;
        ////sph_bessel(l, h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][gridPointIdx]) *
        ////mRadialGrid[z][gridPointIdx] * mRadialGrid[z][gridPointIdx] * mRadialIntegrationWeights[z][gridPointIdx];
        //        }
    }


    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        int z = mAtomicNumber[atomIdx];
        double realPart = 0;
        double imagPart = 0;
        vector<vector<double> >& v = mRadialIntegrationMultipliers[z];
        for (l = 0; l < 8; l++)
        {

            //int i, n = mRadialGrid[z].size();

            //for (i = 0; i < n; i++)
            //    mRadialIntegrationMultiplier[i] =
            //    sph_bessel(l, h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][i]) *
            //    mRadialGrid[z][i] * mRadialGrid[z][i] * mRadialIntegrationWeights[mAtomicNumber[atomIdx]][i];

            //sphericalBesselTransform += sph_bessel(l, h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][i]) *
            //mRadialGrid[z][i] * mRadialGrid[z][i] * mSphericalHarmonicsCoefficients[atomIdx][l][l + m][i] *
            //mRadialIntegrationWeights[mAtomicNumber[atomIdx]][i];
            vector<double>& vv = v[l];
            for (m = -l; m <= l; m++)
            {
                double sphericalBesselTransform = 0;
                const vector<double>& sphCoeff = mSphericalHarmonicsCoefficients[atomIdx][l][l + m];
                for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
                    sphericalBesselTransform += vv[gridPointIdx] * sphCoeff[gridPointIdx];
                //sphericalBesselTransform += mRadialIntegrationMultiplier[gridPointIdx] * mSphericalHarmonicsCoefficients[atomIdx][l][l + m][gridPointIdx];


                if (l % 2 == 0)
                    realPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
                else
                    imagPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
            }
        }
        ff[atomIdx] = complex<double>(realPart, imagPart);
    }

}

void StockholderAtomBankFormFactorCalculationManager::calculateCartSphExp(
    const Vector3d& hkl,
    std::vector< std::complex<double> >& ff,
    const std::vector<bool>& includeAtom)
    const
{
    double h_length = sqrt(hkl * hkl);
    if (h_length == 0.0)
        return calculateCart000(ff, includeAtom);
    //cout << "*\n";
    int atomIdx, nAtoms = mAtomicNumber.size();
    ff.resize(nAtoms);
    vector<vector<double> > sphericalHarmonics(8);
    int l, m;

    for (l = 0; l < 8; l++)
        sphericalHarmonics[l].resize(2 * l + 1);

    double h_length_atomic_units = h_length / constants::Angstrom;
    Vector3d h_normalized(1, 0, 0);

    if (h_length > 0)
        h_normalized = hkl / h_length;

    real_spherical_harmonics::getWfnNormalized<7>(h_normalized, sphericalHarmonics);

    vector<double> iPowL = { 1, 1, -1, -1, 1, 1, -1, -1 };

    int gridPointIdx, gridSize = mIntegrationGrid.radialGridSize;

    for (int z : mUniqueAtomicNumbers)
    {
        vector<vector<double> >& v = mRadialIntegrationMultipliers[z];
        vector<double>& two_pi_h_r = mAuxRadial_1;

        for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
        {
            double x = h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][gridPointIdx];
            two_pi_h_r[gridPointIdx] = x;
        }

        const vector<double>& r2rw = mR2RW[z];

        for (l = 0; l <8; l++)
        {
            vector<double>& vv = v[l];
            double sb_l, x;
            for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
            {
                x = two_pi_h_r[gridPointIdx];
                sb_l = sph_bessel(l, x);
                vv[gridPointIdx] = sb_l * r2rw[gridPointIdx];
            }
        }
    }


    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        int z = mAtomicNumber[atomIdx];
        double realPart = 0;
        double imagPart = 0;
        vector<vector<double> >& v = mRadialIntegrationMultipliers[z];

        for (l = 0; l < 8; l++)
        {

            vector<double>& vv = v[l];
            for (m = -l; m <= l; m++)
            {
                double sphericalBesselTransform = 0;
                const vector<double>& sphCoeff = mSphericalHarmonicsCoefficients[atomIdx][l][l + m];
                for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
                    sphericalBesselTransform += vv[gridPointIdx] * sphCoeff[gridPointIdx];
                

                if (l % 2 == 0)
                    realPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
                else
                    imagPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
            }
        }
        ff[atomIdx] = complex<double>(realPart, imagPart);
    }

}

void StockholderAtomBankFormFactorCalculationManager::calculateCartSphExp_in_thread(
    const Vector3d& hkl,
    std::vector< std::complex<double> >& ff,
    const std::vector<bool>& includeAtom)
    const
{
    int threadIdx = omp_get_thread_num();

    double h_length = sqrt(hkl * hkl);
    if (h_length == 0.0)
        return calculateCart000(ff, includeAtom);
    //cout << "*\n";
    int atomIdx, nAtoms = mAtomicNumber.size();
    ff.resize(nAtoms);
    vector<vector<double> > sphericalHarmonics(8);
    int l, m;

    for (l = 0; l < 8; l++)
        sphericalHarmonics[l].resize(2 * l + 1);

    double h_length_atomic_units = h_length / constants::Angstrom;
    Vector3d h_normalized(1, 0, 0);

    if (h_length > 0)
        h_normalized = hkl / h_length;

    real_spherical_harmonics::getWfnNormalized<7>(h_normalized, sphericalHarmonics);

    vector<double> iPowL = { 1, 1, -1, -1, 1, 1, -1, -1 };

    int gridPointIdx, gridSize = mIntegrationGrid.radialGridSize;

    for (int z : mUniqueAtomicNumbers)
    {
        vector<vector<double> >& v = mRadialIntegrationMultipliers_t[threadIdx][z];
        vector<double>& two_pi_h_r = mAuxRadial_t[threadIdx];

        for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
        {
            double x = h_length_atomic_units * 2.0 * M_PI * mRadialGrid[z][gridPointIdx];
            two_pi_h_r[gridPointIdx] = x;
        }

        const vector<double>& r2rw = mR2RW[z];

        for (l = 0; l < 8; l++)
        {
            vector<double>& vv = v[l];
            double sb_l, x;
            for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
            {
                x = two_pi_h_r[gridPointIdx];
                sb_l = sph_bessel(l, x);
                vv[gridPointIdx] = sb_l * r2rw[gridPointIdx];
            }
        }
    }


    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        int z = mAtomicNumber[atomIdx];
        double realPart = 0;
        double imagPart = 0;
        vector<vector<double> >& v = mRadialIntegrationMultipliers_t[threadIdx][z];

        for (l = 0; l < 8; l++)
        {

            vector<double>& vv = v[l];
            for (m = -l; m <= l; m++)
            {
                double sphericalBesselTransform = 0;
                const vector<double>& sphCoeff = mSphericalHarmonicsCoefficients[atomIdx][l][l + m];
                for (gridPointIdx = 0; gridPointIdx < gridSize; gridPointIdx++)
                    sphericalBesselTransform += vv[gridPointIdx] * sphCoeff[gridPointIdx];


                if (l % 2 == 0)
                    realPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
                else
                    imagPart += sphericalBesselTransform * sphericalHarmonics[l][l + m] * 4.0 * M_PI * iPowL[l];
            }
        }
        ff[atomIdx] = complex<double>(realPart, imagPart);
    }

}


void StockholderAtomBankFormFactorCalculationManager::calculateFrac(
    const vector<Vector3i>& hkl,
    vector < vector<complex<double> > >& formFactors,
    const vector<bool>& includeAtom) const
{
    int hklIdx, nHkl = hkl.size();

    vector<Vector3d> hkl_cart(nHkl);
    
    for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
        mReciprocalLatticeUnitCell.fractionalToCartesian(hkl[hklIdx], hkl_cart[hklIdx]);

    calculateCart(hkl_cart, formFactors, includeAtom);
}

/**
formFactors[hkl idx][atom idx]
*/
void StockholderAtomBankFormFactorCalculationManager::calculateCart(
    const vector <Vector3d>& hkl,
    vector < vector<complex<double> > >& formFactors,
    const vector<bool>& includeAtom) const
{
    int nHkl = hkl.size();
    int nAtoms = mAtomicNumber.size();
    formFactors.resize(nHkl, vector<complex<double> >(nAtoms));

    if (mUseSphericalHarmonicsExpansion)
    {
        cout << "parallel version with " << mN_Threads << " thread(s)\n";
        

        omp_set_num_threads(mN_Threads);
#pragma omp parallel for
        for (int i = 0; i < nHkl; i++)
            calculateCartSphExp_in_thread(hkl[i], formFactors[i], includeAtom);
    }
    else
    {
        for (int i = 0; i < nHkl; i++)
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                formFactors[i][atomIdx] = calculateCart(atomIdx, hkl[i]);
    }

}


}


