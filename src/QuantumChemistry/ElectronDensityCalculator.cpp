#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/IO/wfn_io.h"

#include <cmath>
#include <map>
#include <algorithm>

using namespace std;

namespace discamb
{
    ElectronDensityCalculator::ElectronDensityCalculator(){}
    ElectronDensityCalculator::~ElectronDensityCalculator(){}

    void ElectronDensityCalculator::set_1rdm(
        const std::vector<Vector3d> &center_position,
        const std::vector<int> &primitive_to_center,
        const std::vector<int> &primitive_type,
        const std::vector<double> &primitive_exponents,
        const std::vector<double> &molecular_orbital_occupancy,
        const std::vector<std::vector<double> > &molecular_orbitals)
    {
        mCenterPositions = center_position;
        mPrimitiveToCenter = primitive_to_center;
        mPrimitiveType = primitive_type;
        mPrimitiveExponents = primitive_exponents;
        mMolecularOrbitalOccupancy = molecular_orbital_occupancy;
        mMolecularOrbitals = molecular_orbitals;
        mR2.resize(mCenterPositions.size());
        mR.resize(mCenterPositions.size());
        setGroups();
    }

    void ElectronDensityCalculator::setContributingCenters(
        const std::vector<int>& centers)
    {

        if (mCenterPositions.size() == centers.size())
        {
            mUseSubset = false;
            return;
        }

        mSubsetCalculator = shared_ptr<ElectronDensityCalculator>(new ElectronDensityCalculator);
        mUseSubset = true;

        vector<Vector3d> center_position;
        vector<int> primitive_to_center;
        vector<int> primitive_type;
        vector<double> primitive_exponents;
        vector<double> molecular_orbital_occupancy;
        vector<std::vector<double> > molecular_orbitals;
        vector<bool> countAO, countCenter;
        vector<int> old2newCenterIdx;


        int nCentersOld = mCenterPositions.size();
        old2newCenterIdx.assign(nCentersOld, 0);
        countCenter.assign(nCentersOld, false);
        int nPrimitivesOld = mPrimitiveToCenter.size();
        countAO.assign(nPrimitivesOld, false);

        

        for (auto centerIdx : centers)
        {
            countCenter[centerIdx] = true;
            old2newCenterIdx[centerIdx] = center_position.size();
            center_position.push_back(mCenterPositions[centerIdx]);
        }

        for (int primitiveOldIdx = 0; primitiveOldIdx < nPrimitivesOld; primitiveOldIdx++)
            if (countCenter[mPrimitiveToCenter[primitiveOldIdx]-1])
            {
                countAO[primitiveOldIdx] = true;
                primitive_to_center.push_back(old2newCenterIdx[mPrimitiveToCenter[primitiveOldIdx]-1]+1);
                primitive_type.push_back(mPrimitiveType[primitiveOldIdx]);
                primitive_exponents.push_back(mPrimitiveExponents[primitiveOldIdx]);
            }

        molecular_orbital_occupancy = mMolecularOrbitalOccupancy;
        molecular_orbitals.resize(mMolecularOrbitals.size());
        for (int orbitalIdx = 0; orbitalIdx < mMolecularOrbitals.size(); orbitalIdx++)
            for (int primitiveIdxOld = 0; primitiveIdxOld < nPrimitivesOld; primitiveIdxOld++)
                if (countAO[primitiveIdxOld])
                    molecular_orbitals[orbitalIdx].push_back(mMolecularOrbitals[orbitalIdx][primitiveIdxOld]);

        mSubsetCalculator->set_1rdm(
            center_position,
            primitive_to_center,
            primitive_type,
            primitive_exponents,
            molecular_orbital_occupancy,
            molecular_orbitals);

        //
        if (mEDFsPrimitiveToCenter.empty())
            return;

        vector<int> edfsPrimitiveToCenter;
        vector<int> edfsPrimitiveType;
        vector<double> edfsPrimitiveExponents;
        vector<double> edfsPrimitiveCoefficients;

        int nEdfOld = mEDFsPrimitiveToCenter.size();
        //vector<bool> countEdf(nEdfOld, false);
        
        for (int edfIdx = 0; edfIdx < nEdfOld; edfIdx++)
            if (countCenter[mEDFsPrimitiveToCenter[edfIdx]-1])
            {
                edfsPrimitiveToCenter.push_back(old2newCenterIdx[mEDFsPrimitiveToCenter[edfIdx]-1]+1);
                edfsPrimitiveType.push_back(mEDFsPrimitiveType[edfIdx]);
                edfsPrimitiveExponents.push_back(mEDFsPrimitiveExponents[edfIdx]);
                edfsPrimitiveCoefficients.push_back(mEDFsPrimitiveCoefficients[edfIdx]);
            }

        mSubsetCalculator->setAdditionalDensity(
            edfsPrimitiveToCenter,
            edfsPrimitiveType,
            edfsPrimitiveExponents,
            edfsPrimitiveCoefficients);

    }

    void ElectronDensityCalculator::setFromWavefunctionFile(
        const std::string& fileName)
    {
        wfn_io::WfnFileData wfn;
        wfn_io::read_wavefunction(fileName, wfn);

        set_1rdm(
            wfn.center_position,
            wfn.primitive_to_center,
            wfn.primitive_type,
            wfn.primitive_exponents,
            wfn.molecular_orbital_occupancy,
            wfn.molecular_orbitals);

        if (!wfn.edfs.empty())
            setAdditionalDensity(
                wfn.edfs[0].primitive_to_center,
                wfn.edfs[0].primitive_type,
                wfn.edfs[0].primitive_exponents,
                wfn.edfs[0].primitive_coefficients);

    }

    void ElectronDensityCalculator::setFromWfn(
        const std::string& wfnFile)
    {
        wfn_io::WfnFileData wfn;
        wfn_io::read_wfn(wfnFile, wfn);
        set_1rdm(
            wfn.center_position,
            wfn.primitive_to_center,
            wfn.primitive_type,
            wfn.primitive_exponents,
            wfn.molecular_orbital_occupancy,
            wfn.molecular_orbitals);

    }

    void ElectronDensityCalculator::setFromWfx(
        const std::string& wfxFile)
    {
        wfn_io::WfnFileData wfn;
        wfn_io::read_wfx(wfxFile, wfn);
        
        set_1rdm(
            wfn.center_position,
            wfn.primitive_to_center,
            wfn.primitive_type,
            wfn.primitive_exponents,
            wfn.molecular_orbital_occupancy,
            wfn.molecular_orbitals);

        if (!wfn.edfs.empty())
            setAdditionalDensity(
                wfn.edfs[0].primitive_to_center,
                wfn.edfs[0].primitive_type,
                wfn.edfs[0].primitive_exponents,
                wfn.edfs[0].primitive_coefficients);
    }

    double ElectronDensityCalculator::calculate(
        double x,
        double y,
        double z)
        const
    {
        if (mUseSubset)
            mSubsetCalculator->calculate(x, y, z);
        vector<Vector3d> xyz;
        vector<double> r2;
        vector<double> primitives;
        int centerIdx, nCenters = mCenterPositions.size();
        int nPrimitives, primitiveIdx;
        xyz.resize(nCenters);
        r2.resize(nCenters);
        for (centerIdx = 0; centerIdx < nCenters; centerIdx++)
        {
            xyz[centerIdx] = Vector3d(x, y, z) - mCenterPositions[centerIdx];
            r2[centerIdx] = xyz[centerIdx] * xyz[centerIdx];
        }

        double mo, result=0;

        nPrimitives = mPrimitiveToCenter.size();
        primitives.resize(nPrimitives);
        
        for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
        {
            centerIdx = mPrimitiveToCenter[primitiveIdx] - 1;
            primitives[primitiveIdx] = cartersianPolynomial(mPrimitiveType[primitiveIdx], xyz[centerIdx].x, xyz[centerIdx].y, xyz[centerIdx].z);
            primitives[primitiveIdx] *= exp(-mPrimitiveExponents[primitiveIdx] * r2[centerIdx]);
        }

        int nMo, moIdx;
        nMo = mMolecularOrbitals.size();
        for (moIdx = 0; moIdx < nMo; moIdx++)
        {
            mo = 0;
            for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                mo += primitives[primitiveIdx] * mMolecularOrbitals[moIdx][primitiveIdx];
            result += mo * mo * mMolecularOrbitalOccupancy[moIdx];
        }

        return result + calculateAdditinalDensity(x, y, z);
    }

    void ElectronDensityCalculator::setGroups()
    {
        int primitiveIdx, nPrimitives = mPrimitiveType.size();
        mPrimitiveValue.resize(nPrimitives);

        mExponentGroupPrimitives.clear();
        mPolynomialGroupPrimitives.clear();

        for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
        {
            mExponentGroupPrimitives[{mPrimitiveToCenter[primitiveIdx] - 1, mPrimitiveExponents[primitiveIdx]}].push_back(primitiveIdx);
            mPolynomialGroupPrimitives[{mPrimitiveToCenter[primitiveIdx] - 1, mPrimitiveType[primitiveIdx]}].push_back(primitiveIdx);
        }

        int nCenters = mCenterPositions.size();
        mXyz.resize(nCenters);
        mR2.resize(nCenters);

        //-----------

        mExponentGroupCenter.clear();
        mExponentGroupExponent.clear();
        mExponentGroupValue.clear();
        mPrimitive2ExponentGroup.resize(nPrimitives);
        mPrimitive2PolynomialGroup.resize(nPrimitives);
        mPolynomialGroupType.clear();
        mPolynomialGroupCenter.clear();


        int idx = 0;
        for (auto& x : mExponentGroupPrimitives)
        {
            mExponentGroupCenter.push_back(x.first.first);
            mExponentGroupExponent.push_back(x.first.second);
            mExponentGroupValue.push_back(0);

            for (auto primitiveIdx : x.second)
                mPrimitive2ExponentGroup[primitiveIdx] = idx;

            idx++;
        }

        idx = 0;
        for (auto& x : mPolynomialGroupPrimitives)
        {
            mPolynomialGroupType.push_back(x.first.second);
            mPolynomialGroupCenter.push_back(x.first.first);
            mPolynomialGroupValue.push_back(0);

            for (int primitiveIdx : x.second)
                mPrimitive2PolynomialGroup[primitiveIdx] = idx;
        }

    }

    double ElectronDensityCalculator::calculate3(
        double x,
        double y,
        double z)
        const
    {
        if (mUseSubset)
            mSubsetCalculator->calculate3(x, y, z);


        int centerIdx, nCenters = mCenterPositions.size();
        int primitiveIdx, nPrimitives = mPrimitiveToCenter.size();

        for (centerIdx = 0; centerIdx < nCenters; centerIdx++)
        {
            mXyz[centerIdx] = Vector3d(x, y, z) - mCenterPositions[centerIdx];
            mR2[centerIdx] = mXyz[centerIdx] * mXyz[centerIdx];
        }

        int idx, n = mExponentGroupValue.size();
        for(idx=0; idx<n;idx++)
            mExponentGroupValue[idx] = exp(-mExponentGroupExponent[idx] * mR2[mExponentGroupCenter[idx]]);
        n = mPolynomialGroupValue.size();
        for (idx = 0; idx < n; idx++)
        {
            const Vector3d& r = mXyz[mPolynomialGroupCenter[idx]];
            mPolynomialGroupValue[idx] = cartersianPolynomial(mPolynomialGroupType[idx], r.x, r.y, r.z);
        }

        double mo, result = 0;



        for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
        {
            mPrimitiveValue[primitiveIdx] = mExponentGroupValue[mPrimitive2ExponentGroup[primitiveIdx]];
            mPrimitiveValue[primitiveIdx] *= mPolynomialGroupValue[mPrimitive2PolynomialGroup[primitiveIdx]];
        }

        int nMo, moIdx;
        nMo = mMolecularOrbitals.size();
        for (moIdx = 0; moIdx < nMo; moIdx++)
        {
            mo = 0;
            for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                mo += mPrimitiveValue[primitiveIdx] * mMolecularOrbitals[moIdx][primitiveIdx];
            result += mo * mo * mMolecularOrbitalOccupancy[moIdx];
        }

        return result + calculateAdditinalDensity(x, y, z);
    }

    double ElectronDensityCalculator::calculate4(
        double x,
        double y,
        double z)
        const
    {

        if (mUseSubset)
            mSubsetCalculator->calculate4(x, y, z);

        int centerIdx, nCenters = mCenterPositions.size();
        int primitiveIdx, nPrimitives = mPrimitiveToCenter.size();

        for (centerIdx = 0; centerIdx < nCenters; centerIdx++)
        {
            mXyz[centerIdx] = Vector3d(x, y, z) - mCenterPositions[centerIdx];
            mR2[centerIdx] = mXyz[centerIdx] * mXyz[centerIdx];
        }

        int idx, n = mExponentGroupValue.size();
        for (idx = 0; idx < n; idx++)
            mExponentGroupValue[idx] = exp(-mExponentGroupExponent[idx] * mR2[mExponentGroupCenter[idx]]);
        n = mPolynomialGroupValue.size();
        for (idx = 0; idx < n; idx++)
        {
            const Vector3d& r = mXyz[mPolynomialGroupCenter[idx]];
            mPolynomialGroupValue[idx] = cartersianPolynomial(mPolynomialGroupType[idx], r.x, r.y, r.z);
        }

        double mo, result = 0;



        for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
        {
            mPrimitiveValue[primitiveIdx] = mExponentGroupValue[mPrimitive2ExponentGroup[primitiveIdx]];
            mPrimitiveValue[primitiveIdx] *= mPolynomialGroupValue[mPrimitive2PolynomialGroup[primitiveIdx]];
        }

        int nMo, moIdx;
        nMo = mMolecularOrbitals.size();
        for (moIdx = 0; moIdx < nMo; moIdx++)
        {
            mo = 0;
            for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                mo += mPrimitiveValue[primitiveIdx] * mMolecularOrbitals[moIdx][primitiveIdx];
            result += mo * mo * mMolecularOrbitalOccupancy[moIdx];
        }

        return result + calculateAdditinalDensity(x, y, z);
    }

    double ElectronDensityCalculator::calculate2(
        const Vector3d& r) 
        const
    {
        return calculate2(r.x, r.y, r.z);
    }

    double ElectronDensityCalculator::calculate2(
        double x,
        double y,
        double z)
        const
    {
        if (mUseSubset)
        {
            double d = mSubsetCalculator->calculate2(x, y, z);
            return d;
        }


        int centerIdx, nCenters = mCenterPositions.size();
        int nPrimitives, primitiveIdx;

        for (centerIdx = 0; centerIdx < nCenters; centerIdx++)
        {
            mXyz[centerIdx] = Vector3d(x, y, z) - mCenterPositions[centerIdx];
            mR2[centerIdx] = mXyz[centerIdx] * mXyz[centerIdx];
        }

        nPrimitives = mPrimitiveToCenter.size();
        
        
        double d;
        for (auto& x : mExponentGroupPrimitives)
        {
            d = exp(-x.first.second * mR2[x.first.first]);
            for (int primitiveIdx : x.second)
                mPrimitiveValue[primitiveIdx] = d;
        }   

        for (auto& x : mPolynomialGroupPrimitives)
        {
            centerIdx = x.first.first;
            d = cartersianPolynomial(x.first.second, mXyz[centerIdx].x, mXyz[centerIdx].y, mXyz[centerIdx].z);
            for (int primitiveIdx : x.second)
                mPrimitiveValue[primitiveIdx] *= d;
        }


        double mo, result = 0;

        
        //primitives.resize(nPrimitives);

        //for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
        //{
        //    centerIdx = mPrimitiveToCenter[primitiveIdx] - 1;
        //    primitives[primitiveIdx] = cartersianPolynomial(mPrimitiveType[primitiveIdx], xyz[centerIdx].x, xyz[centerIdx].y, xyz[centerIdx].z);
        //    primitives[primitiveIdx] *= exp(-mPrimitiveExponents[primitiveIdx] * r2[centerIdx]);
        //}

        int nMo, moIdx;
        nMo = mMolecularOrbitals.size();
        for (moIdx = 0; moIdx < nMo; moIdx++)
        {
            mo = 0;
            for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                mo += mPrimitiveValue[primitiveIdx] * mMolecularOrbitals[moIdx][primitiveIdx];
            result += mo * mo * mMolecularOrbitalOccupancy[moIdx];
        }

        return result + calculateAdditinalDensity(x, y, z);
    }

    double ElectronDensityCalculator::calculate2(
        double x,
        double y,
        double z,
        int molecularOrbitalIdx)
        const
    {
        int centerIdx, nCenters = mCenterPositions.size();
        int nPrimitives, primitiveIdx;

        for (centerIdx = 0; centerIdx < nCenters; centerIdx++)
        {
            mXyz[centerIdx] = Vector3d(x, y, z) - mCenterPositions[centerIdx];
            mR2[centerIdx] = mXyz[centerIdx] * mXyz[centerIdx];
        }

        nPrimitives = mPrimitiveToCenter.size();


        double d;
        for (auto& x : mExponentGroupPrimitives)
        {
            d = exp(-x.first.second * mR2[x.first.first]);
            for (int primitiveIdx : x.second)
                mPrimitiveValue[primitiveIdx] = d;
        }

        for (auto& x : mPolynomialGroupPrimitives)
        {
            centerIdx = x.first.first;
            d = cartersianPolynomial(x.first.second, mXyz[centerIdx].x, mXyz[centerIdx].y, mXyz[centerIdx].z);
            for (int primitiveIdx : x.second)
                mPrimitiveValue[primitiveIdx] *= d;
        }


        double mo, result = 0;


        //primitives.resize(nPrimitives);

        //for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
        //{
        //    centerIdx = mPrimitiveToCenter[primitiveIdx] - 1;
        //    primitives[primitiveIdx] = cartersianPolynomial(mPrimitiveType[primitiveIdx], xyz[centerIdx].x, xyz[centerIdx].y, xyz[centerIdx].z);
        //    primitives[primitiveIdx] *= exp(-mPrimitiveExponents[primitiveIdx] * r2[centerIdx]);
        //}

        //int nMo, moIdx;
        //nMo = mMolecularOrbitals.size();
        //for (moIdx = 0; moIdx < nMo; moIdx++)
        //{
            mo = 0;
            for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                mo += mPrimitiveValue[primitiveIdx] * mMolecularOrbitals[molecularOrbitalIdx][primitiveIdx];
            result += mo * mo * mMolecularOrbitalOccupancy[molecularOrbitalIdx];
        //}

            return result + calculateAdditinalDensity(x, y, z);
    }

    double ElectronDensityCalculator::calculateAdditinalDensity(
        double x,
        double y,
        double z)
        const
    {
        //return 0;
        if (mEDFsPrimitiveToCenter.empty())
            return 0;
        double d = 0;
        Vector3d r(x,y,z);
        for (auto center : mEDFsCenters)
        {
            int centerIdx = center - 1;
            mR[centerIdx] = r - mCenterPositions[centerIdx];
            mR2[centerIdx] = mR[centerIdx]* mR[centerIdx];
        }
        /*
        std::vector<int> mEDFsPrimitiveToCenter;
        std::vector<int> mEDFsPrimitiveType;
        std::vector<double> mEDFsPrimitiveExponents;
        std::vector<double> mEDFsPrimitiveCoefficients;
        bool mHasNoSTypeEDFs = false;
        */

        int i, n = mEDFsPrimitiveToCenter.size();

        if (mOnlySTypeEDFs)
            for (i = 0; i < n; i++)
                d += mEDFsPrimitiveCoefficients[i] * exp(-mEDFsPrimitiveExponents[i] * mR2[mEDFsPrimitiveToCenter[i] - 1]);
        else
            for (i = 0; i < n; i++)
            {
                int centerIdx = mEDFsPrimitiveToCenter[i] - 1;
                auto& xyz = mR[centerIdx];
                d += mEDFsPrimitiveCoefficients[i] * exp(-mEDFsPrimitiveExponents[i] * mR2[centerIdx]) *
                    cartersianPolynomial(mEDFsPrimitiveType[i], xyz.x, xyz.y, xyz.z);
            }
        return d;
    }

    ElectronDensityCalculator* ElectronDensityCalculator::clone()
        const
    {

        return new ElectronDensityCalculator(*this);
    }


    double ElectronDensityCalculator::cartersianPolynomial(
        int idx,
        double x,
        double y,
        double z)
    {
        switch (idx)
        {
            // S
            case 1: return 1.0;
            // P
            case 2: return x;
            case 3: return y;
            case 4: return z;
            // D
            case 5: return x * x;
            case 6: return y * y;
            case 7: return z * z;
            case 8: return x * y;
            case 9: return x * z;
            case 10: return y * z;
            // F
            case 11: return x * x*x;
            case 12: return y * y*y;
            case 13: return z * z*z;
            case 14: return x * x*y;
            case 15: return x * x*z;
            case 16: return y * y*z;
            case 17: return x * y*y;
            case 18: return x * z*z;
            case 19: return y * z*z;
            case 20: return x * y*z;
            // G
            case 21: return x * x*x*x;
            case 22: return y * y*y*y;
            case 23: return z * z*z*z;
            case 24: return x * x*x*y;
            case 25: return x * x*x*z;
            case 26: return x * y*y*y;
            case 27: return y * y*y*z;
            case 28: return x * z*z*z;
            case 29: return y * z*z*z;
            case 30: return x * x*y*y;
            case 31: return x * x*z*z;
            case 32: return y * y*z*z;
            case 33: return x * x*y*z;
            case 34: return x * y*y*z;
            case 35: return x * y*z*z;
            // H
            case 36: return z * z * z * z * z;
            case 37: return y * z * z * z * z;
            case 38: return y * y * z * z * z;
            case 39: return y * y * y * z * z;
            case 40: return y * y * y * y * z;
            case 41: return y * y * y * y * y;
            case 42: return x * z * z * z * z;
            case 43: return x * y * z * z * z;
            case 44: return x * y * y * z * z;
            case 45: return x * y * y * y * z;
            case 46: return x * y * y * y * y;
            case 47: return x * x * z * z * z;
            case 48: return x * x * y * z * z;
            case 49: return x * x * y * y * z;
            case 50: return x * x * y * y * y;
            case 51: return x * x * x * z * z;
            case 52: return x * x * x * y * z;
            case 53: return x * x * x * y * y;
            case 54: return x * x * x * x * z;
            case 55: return x * x * x * x * y;
            case 56: return x * x * x * x * x;
            // I
            case 57: return z * z * z * z * z * z;
            case 58: return y * z * z * z * z * z;
            case 59: return y * y * z * z * z * z;
            case 60: return y * y * y * z * z * z;
            case 61: return y * y * y * y * z * z;
            case 62: return y * y * y * y * y * z;
            case 63: return y * y * y * y * y * y;
            case 64: return x * z * z * z * z * z;
            case 65: return x * y * z * z * z * z;
            case 66: return x * y * y * z * z * z;
            case 67: return x * y * y * y * z * z;
            case 68: return x * y * y * y * y * z;
            case 69: return x * y * y * y * y * y;
            case 70: return x * x * z * z * z * z;
            case 71: return x * x * y * z * z * z;
            case 72: return x * x * y * y * z * z;
            case 73: return x * x * y * y * y * z;
            case 74: return x * x * y * y * y * y;
            case 75: return x * x * x * z * z * z;
            case 76: return x * x * x * y * z * z;
            case 77: return x * x * x * y * y * z;
            case 78: return x * x * x * y * y * y;
            case 79: return x * x * x * x * z * z;
            case 80: return x * x * x * x * y * z;
            case 81: return x * x * x * x * y * y;
            case 82: return x * x * x * x * x * z;
            case 83: return x * x * x * x * x * y;
            case 84: return x * x * x * x * x * x;
        }
        on_error::throwException("undefined proimitive function type", __FILE__, __LINE__);
        return 0.0;
    }

    void ElectronDensityCalculator::setAdditionalDensity(
        const std::vector<int>& primitive_to_center,
        const std::vector<int>& primitive_type,
        const std::vector<double>& primitive_exponents,
        const std::vector<double>& primitive_coefficients)
    {
        mEDFsPrimitiveToCenter = primitive_to_center;
        mEDFsPrimitiveType = primitive_type;
        mEDFsPrimitiveExponents = primitive_exponents;
        mEDFsPrimitiveCoefficients = primitive_coefficients;
        
        if (!primitive_to_center.empty())
        {
            mOnlySTypeEDFs = (*max_element(primitive_type.begin(), primitive_type.end()) == 1);
            mEDFsCenters.insert(primitive_to_center.begin(), primitive_to_center.end());
        }

        mR2.resize(mCenterPositions.size());
        mR.resize(mCenterPositions.size());

    }
}