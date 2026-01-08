#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/BasicUtilities/on_error.h"

using namespace std;

namespace discamb {

	LocalCoordinateSystemCalculator::LocalCoordinateSystemCalculator()
	{
        mIsCartesian = true;
        mThirdCoordinate = 2;
	}

	LocalCoordinateSystemCalculator::LocalCoordinateSystemCalculator(
		const LocalCoordinateSystem<AtomInCrystalID> &lcs,
		const Crystal &c)
	{
		set(lcs, c);
	}

	LocalCoordinateSystemCalculator::~LocalCoordinateSystemCalculator()
	{
	}

	void LocalCoordinateSystemCalculator::set(
		const std::string &definition,
		const Crystal &c)
	{
        mIsCartesian = false;
        on_error::not_implemented(__FILE__, __LINE__);
	}

	void LocalCoordinateSystemCalculator::set(
		const LocalCoordinateSystem<AtomInCrystalID> &lcs,
		const Crystal &c)
	{
        vector<string> xyz{ "X", "Y", "Z" };
        mIsCartesian = false;
        mLcs = lcs;
        mCrossProductLcs = shared_ptr< CrossProductLcs>(new CrossProductLcs(xyz[lcs.coordinate_1], xyz[lcs.coordinate_2], lcs.isR));
        mThirdCoordinate = 3 - lcs.coordinate_1 - lcs.coordinate_2;
	}

    void LocalCoordinateSystemCalculator::calcAtomPosition(
        const AtomInCrystalID &atom,
        const Crystal &c,
        Vector3d &r)
    {
        Vector3d r0_transformed, r0 = c.atoms[atom.index()].coordinates;
        atom.getSymmetryOperation().apply(r0, r0_transformed);
        c.unitCell.fractionalToCartesian(r0_transformed, r);
    }

    void LocalCoordinateSystemCalculator::calcAtomicPositionsBasedDirection(
        const std::vector<AtomInCrystalID> &refPoint,
        const AtomInCrystalID &centralAtom,
        const Crystal &c,
        LcsDirectionType type,
        Vector3d &direction)
        const
    {
        Vector3d centralAtomPosition, atomPosition;

        if (refPoint.size() == 0 || ( (type != LcsDirectionType::AVERAGE_POSITION) && (type != LcsDirectionType::AVERAGE_DIRECTION) ) )
            on_error::throwException("problem when calculating atomic local coordinate system", __FILE__, __LINE__);

        calcAtomPosition(centralAtom, c, centralAtomPosition);

        if (type == LcsDirectionType::AVERAGE_DIRECTION)
        {
            direction.set(0, 0, 0);
            Vector3d directionToNeighbor;
            for (auto atom : refPoint)
            {
                calcAtomPosition(atom, c, atomPosition);
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
                calcAtomPosition(atom, c, atomPosition);
                averagePosition += atomPosition;
            }
            averagePosition /= double(refPoint.size());
            direction = averagePosition - centralAtomPosition;
        }

        direction /= sqrt(direction*direction);
    }

    void LocalCoordinateSystemCalculator::calculateAnyOrthogonal(
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

    void LocalCoordinateSystemCalculator::takeChiralityIntoAccount(
        Vector3d &x,
        Vector3d &y,
        Vector3d &z,
        const Crystal &c,
        bool &sameChirality)
        const
    {
        sameChirality = true;
        if (mLcs.chirality.empty())
            return;
        Vector3d r1, r2, r3, r0, d1, d2, d3;

        calcAtomPosition(mLcs.centralAtom, c, r0);
        calcAtomPosition(mLcs.chirality[0], c, r1);
        calcAtomPosition(mLcs.chirality[1], c, r2);
        calcAtomPosition(mLcs.chirality[2], c, r3);

        d1 = r1 - r0;
        d2 = r2 - r0;
        d3 = r3 - r0;

        if (cross_product(d1, d2)*d3 < 0)
        {
            (mThirdCoordinate == 0 ? x : mThirdCoordinate == 1 ? y : z) *= -1;
            sameChirality = false;
        }
    }

    void LocalCoordinateSystemCalculator::calculate(
        Matrix3d &m,
        const Crystal &c,
        bool &sameChirality)
        const
    {
        Vector3d x, y, z;
        this->calculate(x, y, z, c, sameChirality);
        m.set(x[0], y[0], z[0],
            x[1], y[1], z[1],
            x[2], y[2], z[2]);
    }

	void LocalCoordinateSystemCalculator::calculate(
		Vector3d &x,
		Vector3d &y,
		Vector3d &z,
		const Crystal &c)
		const
	{
        bool sameChirality;
        calculate(x, y, z, c, sameChirality);
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


    void LocalCoordinateSystemCalculator::calculate(
        Vector3d &x, 
        Vector3d &y, 
        Vector3d &z, 
        const Crystal &c, 
        bool &sameChirality) 
        const
    {
        if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL && mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
            mIsCartesian = true;
        if (mIsCartesian)
        {
            x.set(1, 0, 0);
            y.set(0, 1, 0);
            z.set(0, 0, 1);
            return;
        }



        if (mLcs.direction1_type == LcsDirectionType::NOT_SET || mLcs.direction2_type == LcsDirectionType::NOT_SET)
            on_error::throwException("an attempt to calculate local coordinate system with no definition", __FILE__, __LINE__);

        //if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL && mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
          //  on_error::throwException("invalid definition of local coordinate system - both defining vectors defined as ANY_ORTHOGONAL", __FILE__, __LINE__);

        Vector3d r1, r2;
        if (mLcs.direction1_type == LcsDirectionType::ANY_ORTHOGONAL)
        {
            calcAtomicPositionsBasedDirection(mLcs.refPoint_2, mLcs.centralAtom, c, mLcs.direction2_type, r2);
            calculateAnyOrthogonal(r2, r1);
        }
        else
        {
            calcAtomicPositionsBasedDirection(mLcs.refPoint_1, mLcs.centralAtom, c, mLcs.direction1_type, r1);
            if (mLcs.direction2_type == LcsDirectionType::ANY_ORTHOGONAL)
                calculateAnyOrthogonal(r1, r2);
            else 
                calcAtomicPositionsBasedDirection(mLcs.refPoint_2, mLcs.centralAtom, c, mLcs.direction2_type, r2);
            
        }

        mCrossProductLcs->calculate(r1, r2, x, y, z);
        
        takeChiralityIntoAccount(x, y, z, c, sameChirality);
    }

}
