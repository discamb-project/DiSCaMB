#include "discamb/AtomTyping/MolecularLcsCalculator.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb {

    MolecularLcsCalculator::MolecularLcsCalculator()
	{
        mIsCartesian = true;
        mThirdCoordinate = 2;
	}

    MolecularLcsCalculator::MolecularLcsCalculator(
		const LocalCoordinateSystem<int> &lcs)
	{
		set(lcs);
	}

    MolecularLcsCalculator::~MolecularLcsCalculator()
	{
	}


	void MolecularLcsCalculator::set(
		const LocalCoordinateSystem<int> &lcs)
	{
        vector<string> xyz{ "X","Y","Z" };
        mIsCartesian = false;
        mLcs = lcs;
        mCrossProductLcs = shared_ptr< CrossProductLcs>(new CrossProductLcs(xyz[lcs.coordinate_1], xyz[lcs.coordinate_2], lcs.isR));
        mThirdCoordinate = 3 - lcs.coordinate_1 - lcs.coordinate_2;

	}


    void MolecularLcsCalculator::calcAtomicPositionsBasedDirection(
        const std::vector<int> &refPoint,
        const int &centralAtom,
        const std::vector<Vector3d>& r,
        LcsDirectionType type,
        Vector3d &direction)
        const
    {
        Vector3d centralAtomPosition, atomPosition;

        if (refPoint.size() == 0 || ( (type != LcsDirectionType::AVERAGE_POSITION) && (type != LcsDirectionType::AVERAGE_DIRECTION) ) )
            on_error::throwException("problem when calculating atomic local coordinate system", __FILE__, __LINE__);

        centralAtomPosition = r[centralAtom];

        if (type == LcsDirectionType::AVERAGE_DIRECTION)
        {
            direction.set(0, 0, 0);
            Vector3d directionToNeighbor;
            for (auto atom : refPoint)
            {
                atomPosition = r[atom];
                directionToNeighbor = atomPosition - centralAtomPosition;
                directionToNeighbor /= sqrt(directionToNeighbor*directionToNeighbor);
                direction += directionToNeighbor;
            }
        }
        if (type == LcsDirectionType::AVERAGE_POSITION)
        {
            direction.set(0, 0, 0);
            Vector3d averagePosition;
            for (auto atom : refPoint)
            {
                atomPosition = r[atom];
                averagePosition += atomPosition;
            }
            averagePosition /= static_cast<double>(refPoint.size());
            direction = averagePosition - centralAtomPosition;
        }

        direction /= sqrt(direction*direction);
    }

    void MolecularLcsCalculator::calculateAnyOrthogonal(
        const Vector3d &r0,
        Vector3d &r)
    {
        Vector3d x(1, 0, 0), y(0, 1, 0),rx,ry;
        rx = cross_product(r0, x);
        ry = cross_product(r0, y);

        if (rx*rx > ry*ry)
            r = rx;
        else
            r = ry;
        r /= sqrt(r*r);
    }

    void MolecularLcsCalculator::takeChiralityIntoAccount(
        Vector3d &x,
        Vector3d &y,
        Vector3d &z,
        const std::vector<Vector3d>& r,
        bool &sameChirality)
        const
    {
        sameChirality = true;
        if (mLcs.chirality.empty())
            return;
        Vector3d r1, r2, r3, r0, d1, d2, d3;

        r0 = r[mLcs.centralAtom];
        r1 = r[mLcs.chirality[0]];
        r2 = r[mLcs.chirality[1]];
        r3 = r[mLcs.chirality[2]];

        d1 = r1 - r0;
        d2 = r2 - r0;
        d3 = r3 - r0;

        if (cross_product(d1, d2)*d3 < 0)
        {
            (mThirdCoordinate == 0 ? x : mThirdCoordinate == 1 ? y : z) *= -1;
            sameChirality = false;
        }
    }

    void MolecularLcsCalculator::calculate(
        Matrix3d &m,
        const std::vector<Vector3d>& r,
        bool &sameChirality)
        const
    {
        Vector3d x, y, z;
        this->calculate(x, y, z, r, sameChirality);
        m.set(x[0], y[0], z[0],
            x[1], y[1], z[1],
            x[2], y[2], z[2]);
    }

	void MolecularLcsCalculator::calculate(
		Vector3d &x,
		Vector3d &y,
		Vector3d &z,
        const std::vector<Vector3d>& r)
		const
	{
        bool sameChirality;
        calculate(x, y, z, r, sameChirality);
        /*
        if (mLcs.direction1_type == LcsDirectionType::NOT_SET || mLcs.direction2_type == LcsDirectionType::NOT_SET)
            on_error::throwException("an attempt to calculate local coordinate system with no definition", __FILE__, __LINE__);

        //if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL && mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
          //  on_error::throwException("invalid definition of local coordinate system - both defining vectors defined as ANY_ORTHOGONAL", __FILE__, __LINE__);

        if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL && mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL)
            mIsCartesian = true;
        if (mIsCartesian)
        {
            x.set(1, 0, 0);
            y.set(0, 1, 0);
            z.set(0, 0, 1);
            return;
        }
        Vector3d r1, r2;
        if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL)
        {
            calcAtomicPositionsBasedDirection(mLcs.refPoint_2, mLcs.centralAtom, c, mLcs.direction2_type, r2);
            calculateAnyOrthogonal(r2, r1);
        }
        else
        {
            if (mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
            {
                calcAtomicPositionsBasedDirection(mLcs.refPoint_1, mLcs.centralAtom, c, mLcs.direction1_type, r1);
                calculateAnyOrthogonal(r1, r2);
            }
            else
            {
                calcAtomicPositionsBasedDirection(mLcs.refPoint_1, mLcs.centralAtom, c, mLcs.direction1_type, r1);
                calcAtomicPositionsBasedDirection(mLcs.refPoint_2, mLcs.centralAtom, c, mLcs.direction2_type, r2);
            }
        }

        mCrossProductLcs->calculate(r1, r2, x, y, z);
        bool sameChirality;
        takeChiralityIntoAccount(x, y, z, c, sameChirality);
        */
	}

    void MolecularLcsCalculator::calculate(
        Vector3d& x,
        Vector3d& y,
        Vector3d& z,
        const std::vector<Vector3d>& r,
        bool& sameChirality)
        const
    {
        double cosAngle;
        calculate(x, y, z, r, sameChirality, cosAngle);
    }

    void MolecularLcsCalculator::calculate(
        Vector3d &x, 
        Vector3d &y, 
        Vector3d &z, 
        const std::vector<Vector3d>& r,
        bool &sameChirality,
        double& cosAngle)
        const
    {
        if (mLcs.direction1_type == LcsDirectionType::NOT_SET || mLcs.direction2_type == LcsDirectionType::NOT_SET)
            on_error::throwException("an attempt to calculate local coordinate system with no definition", __FILE__, __LINE__);

        //if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL && mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
          //  on_error::throwException("invalid definition of local coordinate system - both defining vectors defined as ANY_ORTHOGONAL", __FILE__, __LINE__);

        if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL && mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
            mIsCartesian = true;
        if (mIsCartesian)
        {
            x.set(1, 0, 0);
            y.set(0, 1, 0);
            z.set(0, 0, 1);
            return;
        }
        Vector3d r1, r2;
        if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL)
        {
            calcAtomicPositionsBasedDirection(mLcs.refPoint_2, mLcs.centralAtom, r, mLcs.direction2_type, r2);
            calculateAnyOrthogonal(r2, r1);
        }
        else
        {
            if (mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
            {
                calcAtomicPositionsBasedDirection(mLcs.refPoint_1, mLcs.centralAtom, r, mLcs.direction1_type, r1);
                calculateAnyOrthogonal(r1, r2);
            }
            else
            {
                calcAtomicPositionsBasedDirection(mLcs.refPoint_1, mLcs.centralAtom, r, mLcs.direction1_type, r1);
                calcAtomicPositionsBasedDirection(mLcs.refPoint_2, mLcs.centralAtom, r, mLcs.direction2_type, r2);
            }
        }

        mCrossProductLcs->calculate(r1, r2, x, y, z);
        cosAngle = r1 * r2 / (r1.norm() * r2.norm());
        takeChiralityIntoAccount(x, y, z, r, sameChirality);
    }

}
