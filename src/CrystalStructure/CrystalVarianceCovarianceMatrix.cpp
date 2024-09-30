#include "discamb/CrystalStructure/CrystalVarianceCovarianceMatrix.h"


#include "discamb/BasicUtilities/OnError.h"

using namespace std;

namespace {
    void identityMatrix(int n, vector<vector<double> >& m)
    {
        m.clear();
        m.resize(n, vector<double>(n, 0.0));
        for (int i = 0; i < n; i++)
            m[i][i] = 1.0;
    }

    void multiplyMatrices(
        const vector<vector<double> >& m1,
        const vector<vector<double> >& m2,
        bool transposeSecondMatrix, // works only for square matrices
        vector<vector<double> >& m)
    {
        int n_col_1, n_row_1, n_col_2, n_row_2;

        n_col_1 = m1[0].size();
        n_row_1 = m1.size();
        n_col_2 = m2[0].size();
        n_row_2 = m2.size();

        if (n_col_1 != n_row_2)
            discamb::on_error::throwException("incompatible matrices sizes", __FILE__, __LINE__);


        m.assign(n_row_1, vector<double>(n_col_2, 0.0));
        int i, j, k;
        if (transposeSecondMatrix)
            for (i = 0; i < n_row_1; i++)
                for (j = 0; j < n_col_2; j++)
                    for (k = 0; k < n_col_1; k++)
                        m[i][j] += m1[i][k] * m2[j][k];
        else
            for (i = 0; i < n_row_1; i++)
                for (j = 0; j < n_col_2; j++)
                    for (k = 0; k < n_col_1; k++)
                        m[i][j] += m1[i][k] * m2[k][j];
    }


}

namespace discamb {

    CrystalVarianceCovarianceMatrix::CrystalVarianceCovarianceMatrix()
    {
        mAdpConvention = structural_parameters_convention::AdpConvention::U_cif;
        mXyzConvention = structural_parameters_convention::XyzCoordinateSystem::fractional;
    }

    CrystalVarianceCovarianceMatrix::~CrystalVarianceCovarianceMatrix()
    {

    }

//    enum variable_type { x, y, z, u_iso, u_11, u_22, u_33, u_12, u_13, u_23, occupancy };

    /*
        structural_parameters_convention::XyzCoordinateSystem mXyzConvention;
    structural_parameters_convention::AdpConvention mAdpConvention;
    std::vector<std::vector<double> > mVcov;
    std::map< std::pair<int, variable_type>, int> mVariable2Idx;
    std::vector<std::pair<int, variable_type> > mVariableList;
    m11, m21, m22, m31, m32, m33, ..
    */
    void CrystalVarianceCovarianceMatrix::set(
        const std::vector< std::vector<double> >& vcovMatrix,
        const std::vector<std::pair<int, variable_type> >& variableList,
        structural_parameters_convention::XyzCoordinateSystem xyzConvention,
        structural_parameters_convention::AdpConvention adpConvention)
    {
        mXyzConvention = xyzConvention;
        mAdpConvention = adpConvention;

        int i, n;

        mVariableList = variableList; 
        mVariableSize.clear();
        map<variable_type, int> variableTypeSize{ {variable_type::xyz,3},{variable_type::u_iso, 1}, {variable_type::u_aniso, 6}, {variable_type::occupancy, 1} };
        
        mVariableDataStartIdx.clear();
        mVariableIdx.clear();
        n = 0;
        for (i = 0; i < variableList.size(); i++)
        {
            mVariableSize.push_back(variableTypeSize[variableList[i].second]);
            mVariableDataStartIdx[variableList[i]] = n;
            mVariableIdx[variableList[i]] = i;
            n += mVariableSize.back();
        }

        mVcov = vcovMatrix;

    }

    void CrystalVarianceCovarianceMatrix::setDiagonal(const Crystal& crystal)
    {
        vector< vector<double> > data;
        vector< pair<int, CrystalVarianceCovarianceMatrix::variable_type> > variableList;
        vector<double> sigmas;
        int atomIdx, nAtoms = crystal.atoms.size();

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            variableList.push_back({ atomIdx,  CrystalVarianceCovarianceMatrix::variable_type::xyz });
            sigmas.push_back(crystal.atoms[atomIdx].coordinates_sigma[0]);
            sigmas.push_back(crystal.atoms[atomIdx].coordinates_sigma[1]);
            sigmas.push_back(crystal.atoms[atomIdx].coordinates_sigma[2]);
            variableList.push_back({ atomIdx,  CrystalVarianceCovarianceMatrix::variable_type::occupancy });
            sigmas.push_back(crystal.atoms[atomIdx].occupancy_sigma);

            if (crystal.atoms[atomIdx].adp.size() == 6)
            {
                variableList.push_back({ atomIdx,  CrystalVarianceCovarianceMatrix::variable_type::u_aniso });
                sigmas.insert(sigmas.end(), crystal.atoms[atomIdx].adp_sigma.begin(), crystal.atoms[atomIdx].adp_sigma.end());
            }

            if (crystal.atoms[atomIdx].adp.size() == 1)
            {
                variableList.push_back({ atomIdx,  CrystalVarianceCovarianceMatrix::variable_type::u_iso });
                sigmas.push_back(crystal.atoms[atomIdx].adp_sigma[0]);
            }

        }
        int n = sigmas.size();
        data.assign(n, vector<double>(n, 0.0));
        for (size_t i = 0; i < n; i++)
            data[i][i] = sigmas[i] * sigmas[i];
        set(data, variableList, structural_parameters_convention::XyzCoordinateSystem::fractional, structural_parameters_convention::AdpConvention::U_cif);

    }
        //vcov.changeConvention(structural_parameters_convention::XyzCoordinateSystem::cartesian, structural_parameters_convention::AdpConvention::U_cif, c.unitCell);
        //for (atomIdx = 0; atomIdx < c.atoms.size(); atomIdx++)
        //{
        //    n = c.atoms[atomIdx].adp.size();
        //    if (n == 6)
        //        for (i = 0; i < 6; i++)
        //            c.atoms[atomIdx].adp_sigma[i] = sqrt(vcov.getVariance(atomIdx, CrystalVarianceCovarianceMatrix::variable_type::u_aniso, i));
        //    if (n == 1)
        //        c.atoms[atomIdx].adp_sigma[0] = sqrt(vcov.getVariance(atomIdx, CrystalVarianceCovarianceMatrix::variable_type::u_iso));

        //}

    //}

    //enum block_variable_type { xyz, u, occ };
//std::vector<int> mBlockStart;
//std::vector<int> mBlockSize;
//std::vector<int> mBlockAtom;
//std::vector<block_variable_type> mBlockVariableType;

    void CrystalVarianceCovarianceMatrix::setBlockInformation()
    {
        //std::vector<std::pair<int, variable_type> > mVariableList;
    }


    double CrystalVarianceCovarianceMatrix::vcov(
        int i, 
        int j)
        const
    {
        if (i < j)
            return mVcov[j][i];
        return mVcov[i][j];
    }

    void CrystalVarianceCovarianceMatrix::getVcovForXyz(
        int atom_1_idx, 
        int atom_2_idx, 
        Matrix3d & m)
    {
        auto start1 = mVariableDataStartIdx.find({ atom_1_idx, variable_type::xyz });
        auto start2 = mVariableDataStartIdx.find({ atom_2_idx, variable_type::xyz });
        if (start1 == mVariableDataStartIdx.end() || start2 == mVariableDataStartIdx.end())
            on_error::throwException("asking for unavailable xyz covariance data", __FILE__, __LINE__);
        int n1 = start1->second;
        int n2 = start2->second;

        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                m(i, j) = vcov(n1+i, n2+j);

    }

    // u11, u22, u33, u12, u13, u23
    void CrystalVarianceCovarianceMatrix::getVcovForAdp(
        int atom_1_idx,
        int atom_2_idx,
        std::vector<std::vector<double> >& m)
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }

    // 3x6, u11, u22, u33, u12, u13, u23
    // atom_1_idx - xyz
    // atom_2_idx - adp
    void CrystalVarianceCovarianceMatrix::getVcovForXyzAdp(
        int atom_1_idx,
        int atom_2_idx,
        std::vector<std::vector<double> >& m)
    {
        on_error::not_implemented(__FILE__, __LINE__);
    }


    void CrystalVarianceCovarianceMatrix::changeXyzConvention(
        structural_parameters_convention::XyzCoordinateSystem xyzConvention,
        const UnitCell& unitCell)
    {

        on_error::not_implemented(__FILE__, __LINE__);
        Matrix3d m;

        if (mXyzConvention == xyzConvention)
            return;

        if(mXyzConvention == structural_parameters_convention::XyzCoordinateSystem::cartesian)
            m = unitCell.getCartesianToFractionalMatrix();

        if (mXyzConvention == structural_parameters_convention::XyzCoordinateSystem::fractional)
            m = unitCell.getFractionalToCartesianMatrix();


    }
    
    void CrystalVarianceCovarianceMatrix::changeAdpConvention(
        structural_parameters_convention::AdpConvention adpConvention, 
        const UnitCell& unitCell)
    {
        on_error::not_implemented(__FILE__, __LINE__); 
    }

    void CrystalVarianceCovarianceMatrix::changeConvention(
        structural_parameters_convention::XyzCoordinateSystem xyzConvention,
        structural_parameters_convention::AdpConvention adpConvention,
        const UnitCell& unitCell)
    {
        //on_error::not_implemented(__FILE__, __LINE__);
        StructuralParametersConverter converter(unitCell);
        vector<vector<double> > xyzTransform(3,vector<double>(3)), adpTransform(6, vector<double>(6)), one_d_transform(1,vector<double>(1,1)),vcov, vcov_transform_t, vcov_transformed;
        int i, j;

        //converter.xyzFractionalToCartesianMatrix()

        // set xyzTransform

        if (xyzConvention == mXyzConvention)
            identityMatrix(3, xyzTransform);
        else
        {
            Matrix3d m;
            if (xyzConvention == structural_parameters_convention::XyzCoordinateSystem::cartesian)
                m = unitCell.getFractionalToCartesianMatrix();
            else
                m = unitCell.getCartesianToFractionalMatrix();

            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                    xyzTransform[i][j] = m(i, j);
        }

        // set adpTransform

        converter.getAdpConversionMatrix(mAdpConvention, adpConvention, adpTransform);
        int nBlocks = mVariableList.size();

        //vector<vector<vector<double> > > transforms = { xyzTransform, adpTransform, one_d_transform };
        map<variable_type, vector<vector<double> > > type2transfromIdx { {variable_type::xyz,xyzTransform}, {variable_type::u_aniso,adpTransform},
                                                                         {variable_type::u_iso, one_d_transform}, {variable_type::occupancy, one_d_transform} };

        for( i=0; i< nBlocks; i++)
            for (j = 0; j < nBlocks; j++)
            {
                getBlock(i, j, vcov);
                vector<vector<double> >& transform1 = type2transfromIdx[mVariableList[i].second];
                vector<vector<double> >& transform2 = type2transfromIdx[mVariableList[j].second];
                
                multiplyMatrices(vcov, transform2, true, vcov_transform_t);
                multiplyMatrices(transform1, vcov_transform_t, false, vcov_transformed);
                setBlock(i, j, vcov_transformed);
            }
        mXyzConvention = xyzConvention;
        mAdpConvention = adpConvention;

    }


    void CrystalVarianceCovarianceMatrix::getData(
        std::vector< std::vector<double> >& vcov_matrix_data,
        std::vector<std::pair<int, variable_type> >& variable_list)
    {
        vcov_matrix_data = mVcov;
        variable_list = mVariableList;
    }

    void CrystalVarianceCovarianceMatrix::getConvention(
        structural_parameters_convention::XyzCoordinateSystem& xyzConvention,
        structural_parameters_convention::AdpConvention& adpConvention) 
        const
    {
        xyzConvention = mXyzConvention;
        adpConvention = mAdpConvention;
    }
    
    std::string CrystalVarianceCovarianceMatrix::variableTypeAsString(
        variable_type v)
    {
        if (v == variable_type::xyz)
            return "xyz";
        if (v == variable_type::u_iso)
            return "u_iso";
        if (v == variable_type::u_aniso)
            return "u_aniso";
        if (v == variable_type::occupancy)
            return "occupancy";
        return "unknown";
    }

    bool CrystalVarianceCovarianceMatrix::hasVariance(
        int atomIdx,
        variable_type variableType)
        const
    {
        auto idx_ptr = mVariableIdx.find({ atomIdx, variableType });
        return (idx_ptr != mVariableIdx.end());
    }

    double CrystalVarianceCovarianceMatrix::getVariance(
        int atomIdx,
        variable_type variableType,
        int variableComponent)
        const
    {
        auto idx_ptr = mVariableIdx.find({ atomIdx, variableType });
        if (idx_ptr == mVariableIdx.end())
            return 0;
            //on_error::throwException("asking for unavailable covariance data", __FILE__, __LINE__);

        int idx = idx_ptr->second;
        int start = mVariableDataStartIdx.find(mVariableList[idx])->second;

        return mVcov[start + variableComponent][start + variableComponent];
    }

    bool CrystalVarianceCovarianceMatrix::get(
        int atom_1_idx, 
        int atom_2_idx, 
        variable_type atom_1_variable, 
        variable_type atom_2_variable,
        std::vector<std::vector<double> >& vcov)
        const
    {
        auto idx_1_ptr = mVariableIdx.find({ atom_1_idx, atom_1_variable });
        auto idx_2_ptr = mVariableIdx.find({ atom_2_idx, atom_2_variable });

        if (idx_1_ptr == mVariableIdx.end() || idx_2_ptr == mVariableIdx.end())
            return false;
            //on_error::throwException("asking for unavailable covariance data", __FILE__, __LINE__);

        int idx1 = idx_1_ptr->second;
        int idx2 = idx_2_ptr->second;

        getBlock(idx1, idx2, vcov);
        return true;
    }

    void CrystalVarianceCovarianceMatrix::getBlock(
        int idx1,
        int idx2,
        std::vector<std::vector<double> >& vcov)
        const
    {
        int n1 = mVariableSize[idx1];
        int n2 = mVariableSize[idx2];

        int start1 = mVariableDataStartIdx.find(mVariableList[idx1])->second;
        int start2 = mVariableDataStartIdx.find(mVariableList[idx2])->second;

        vcov.clear();
        vcov.resize(n1, vector<double>(n2));

        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                vcov[i][j] = mVcov[start1 + i][start2 + j];
    }
    
    void CrystalVarianceCovarianceMatrix::setBlock(
        int idx1,
        int idx2,
        const std::vector<std::vector<double> >& vcov)
    {
        int n1 = mVariableSize[idx1];
        int n2 = mVariableSize[idx2];

        int start1 = mVariableDataStartIdx[mVariableList[idx1]];
        int start2 = mVariableDataStartIdx[mVariableList[idx2]];


        for (int i = 0; i < n1; i++)
            for (int j = 0; j < n2; j++)
                mVcov[start1 + i][start2 + j] = vcov[i][j];

    }
}

