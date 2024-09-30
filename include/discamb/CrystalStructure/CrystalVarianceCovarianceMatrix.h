#pragma once

#include "discamb/CrystalStructure/StructuralParametersConverter.h"
#include "discamb/CrystalStructure/Crystal.h"


#include <vector>
#include <map>

namespace discamb{

    /**
    * \addtogroup CrystalStructure
    * @{
    */


class CrystalVarianceCovarianceMatrix{
public:
        CrystalVarianceCovarianceMatrix();
        ~CrystalVarianceCovarianceMatrix();

        //enum variable_type{x, y, z, u_iso, u_11, u_22, u_33, u_12, u_13, u_23, occupancy};
        enum class variable_type{ xyz, u_iso, u_aniso, occupancy };
        /**
        m11, m21, m22, m31, m32, m33, ..
        */
        void set(const std::vector< std::vector<double> >& vcovMatrix, 
                 const std::vector<std::pair<int, variable_type> > &variableList,
                 structural_parameters_convention::XyzCoordinateSystem xyzConvention = structural_parameters_convention::XyzCoordinateSystem::fractional,
                 structural_parameters_convention::AdpConvention adpConvention = structural_parameters_convention::AdpConvention::U_star);
        void setDiagonal(const Crystal &crystal);
        void changeXyzConvention(structural_parameters_convention::XyzCoordinateSystem xyzConvention, const UnitCell & unitCell);
        void changeAdpConvention(structural_parameters_convention::AdpConvention adpConvention, const UnitCell& unitCell);
        void changeConvention(structural_parameters_convention::XyzCoordinateSystem xyzConvention, 
                              structural_parameters_convention::AdpConvention adpConvention, 
                              const UnitCell& unitCell);
        
        void getData(std::vector< std::vector<double> >& vcovMatrixData,
                     std::vector<std::pair<int, variable_type> >& variableList);

        

        void getConvention(structural_parameters_convention::XyzCoordinateSystem &xyzConvention,
                           structural_parameters_convention::AdpConvention &adpConvention) const;

        bool get(int atom_1_idx, int atom_2_idx, variable_type atom_1_variable, variable_type atom_2_variable, 
                 std::vector<std::vector<double> > &vcov) const;



        void getBlock(int idx1, int idx2, std::vector<std::vector<double> >& vcov) const;

        double getVariance(int atomIdx, variable_type variableType, int variableComponent = 0) const;
        bool hasVariance(int atomIdx, variable_type variableType) const;
        void getVcovForXyz(int atom_1_idx, int atom_2_idx, Matrix3d& m);
        // u11, u22, u33, u12, u13, u23
        void getVcovForAdp(int atom_1_idx, int atom_2_idx, std::vector<std::vector<double> > &m);
        // 3x6, u11, u22, u33, u12, u13, u23
        // atom_1_idx - xyz
        // atom_2_idx - adp
        void getVcovForXyzAdp(int atom_1_idx, int atom_2_idx, std::vector<std::vector<double> >& m);
private:
    structural_parameters_convention::XyzCoordinateSystem mXyzConvention;
    structural_parameters_convention::AdpConvention mAdpConvention;
    std::vector<std::vector<double> > mVcov;
    std::map< std::pair<int, variable_type>, int> mVariableDataStartIdx;
    std::map< std::pair<int, variable_type>, int> mVariableIdx;
    std::vector<std::pair<int, variable_type> > mVariableList;
    std::vector<int> mVariableSize;
    static std::string variableTypeAsString(variable_type v);
    double vcov(int i, int j) const;
    
    
    std::vector<int> mBlockStart;
    std::vector<int> mBlockSize;
    std::vector<int> mBlockAtom;
    void setBlockInformation();
    void setBlock(int idx1, int idx2, const std::vector<std::vector<double> >& vcov);
};
/**@}*/
}

