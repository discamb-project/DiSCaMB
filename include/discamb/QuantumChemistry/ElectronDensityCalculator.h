#pragma once
#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <map>
#include <set>
#include <memory>
#include <string>

namespace discamb
{

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class ElectronDensityCalculator 
    {
    public:
        ElectronDensityCalculator();
        ~ElectronDensityCalculator();
        void set_1rdm(
            const std::vector<Vector3d> &center_position,
            const std::vector<int> &primitive_to_center,
            const std::vector<int> &primitive_type,
            const std::vector<double> &primitive_exponents,
            const std::vector<double> &molecular_orbital_occupancy,
            const std::vector<std::vector<double> > &molecular_orbitals);
        void setFromWfn(const std::string& wfnFile);
        void setFromWfx(const std::string& wfxFile);
        void setFromWavefunctionFile(const std::string& fileName);
        void setAdditionalDensity(
            const std::vector<int> &primitive_to_center,
            const std::vector<int> &primitive_type,
            const std::vector<double> &primitive_exponents,
            const std::vector<double> &primitive_coefficients);
        void setContributingCenters(const std::vector<int>& centers);
        double calculateAdditinalDensity(double x, double y, double z) const;
        double calculate(double x, double y, double z) const;
        double calculate2(double x, double y, double z) const;
        double calculate2(const Vector3d &r) const;
        double calculate2(double x, double y, double z, int orbital) const;
        double calculate3(double x, double y, double z) const;
        double calculate4(double x, double y, double z) const;
        ElectronDensityCalculator* clone() const;
    private:
        // all primitives data
        
        std::vector<Vector3d> mCenterPositions;
        std::vector<int> mPrimitiveToCenter;
        std::vector<int> mPrimitiveType;
        std::vector<double> mPrimitiveExponents;
        std::vector<double> mMolecularOrbitalOccupancy;
        std::vector<std::vector<double> > mMolecularOrbitals;

        
        // end of all primitives data

        std::shared_ptr<ElectronDensityCalculator> mSubsetCalculator;
        bool mUseSubset = false;

        // additional density data
        std::vector<int> mEDFsPrimitiveToCenter;
        std::vector<int> mEDFsPrimitiveType;
        std::vector<double> mEDFsPrimitiveExponents;
        std::vector<double> mEDFsPrimitiveCoefficients;
        bool mOnlySTypeEDFs = true;

        std::set<int> mEDFsCenters;
        mutable std::vector<Vector3d> mR;

        // end of additional density data


        static double cartersianPolynomial(int idx, double x, double y, double z);

        //--------- new ------------
        std::map<std::pair<int, double>, std::vector<int> > mExponentGroupPrimitives;
        // key {polynomial type, center}, value - list of primitives
        std::map<std::pair<int, int>, std::vector<int> > mPolynomialGroupPrimitives;
        void setGroups();
        std::vector<int> mExponentGroupCenter;
        std::vector<double> mExponentGroupExponent;
        mutable std::vector<double> mExponentGroupValue;
        std::vector<int> mPrimitive2ExponentGroup;
        std::vector<int> mPrimitive2PolynomialGroup;
        std::vector<int> mPolynomialGroupType;
        std::vector<int> mPolynomialGroupCenter;
        mutable std::vector<double> mPolynomialGroupValue;


        mutable std::vector<double> mPrimitiveValue;
        mutable std::vector<Vector3d> mXyz;
        mutable std::vector<double> mR2;

    };
    /**@}*/
}

