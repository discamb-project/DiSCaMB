#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include "discamb/BasicUtilities/on_error.h"

#include <cmath>

using namespace std;

namespace discamb {

    UnitCellContent::UnitCellContent()
    {

    }

    UnitCellContent::UnitCellContent(
        const Crystal& crystal)
    {
        set(crystal);
    }

    UnitCellContent::~UnitCellContent()
    {
    }

    void UnitCellContent::set(
        const Crystal &crystal)
    {
        int atomIndex, nAtoms, symmOpIndex, nSymmOp, i, n, equivalentPositionIndex, identityOperationIndex;
        int nPreviouslyAddedAtoms;
        Vector3d position, fractionalCoordinates, fractionalCoordinatesTransformed;
        Vector3i shift;
        vector<Vector3d> positions, positionsAsymmetricUnitAtomsShiftedTo01;
        bool foundNewPosition, foundPositionOfAsymmUnitAtom;
        string symmetryOperationAsString;
        SpaceGroupOperation symmetryOperation, translation, shiftedSymmetryOperation, identityOperation;

        mEquivalentAsymmUnitAtom.clear();


        mCrystal = crystal;

        nAtoms = crystal.atoms.size();
        nSymmOp = crystal.spaceGroup.nSymmetryOperations();

        mAtoms = crystal.atoms;
        mGeneratingOperations.resize(nAtoms);
		mIndicesOfSymmetryEquivalentAtoms.clear();
        mIndicesOfSymmetryEquivalentAtoms.resize(nAtoms);
        positionsAsymmetricUnitAtomsShiftedTo01.resize(nAtoms);

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            position = mAtoms[atomIndex].coordinates;

            for (i = 0; i < 3; i++)
            {
                shift(i) = -math_utilities::roundInt(floor(position(i)));
                position(i) = position(i) + shift(i);
            }
            translation.setToTranslation(Vector3<CrystallographicRational>(shift));
                                    
            mAtoms[atomIndex].coordinates = position;
			if (!translation.isIdentity())
			{
				translation.get(symmetryOperationAsString);
				mAtoms[atomIndex].label += string("(") + symmetryOperationAsString + string(")");
			}

            mGeneratingOperations[atomIndex].push_back(translation);
            mIndicesOfSymmetryEquivalentAtoms[atomIndex].push_back(atomIndex);
            mEquivalentAsymmUnitAtom.push_back(atomIndex);
        }

        nPreviouslyAddedAtoms = nAtoms;

        identityOperationIndex = nSymmOp;

        // check which symmetry operation is identity
        for (symmOpIndex = 0; symmOpIndex < nSymmOp; symmOpIndex++)
            if (identityOperation.isLatticeTranslationEquivalent(crystal.spaceGroup.getSpaceGroupOperation(symmOpIndex), shift))
                identityOperationIndex = symmOpIndex;


        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {

            //positions.clear();
            //mIndicesOfSymmetryEquivalentAtoms[atomIndex].clear();

            positions.clear();
            //mIndicesOfSymmetryEquivalentAtoms[atomIndex].clear();

            // add coordinates of an atom from asymmetric unit for comparison
            position = mAtoms[atomIndex].coordinates;
//            if (crystal.xyzCoordinateSystem == structural_parameters_convention::fractional)
//                position = crystal.atoms[atomIndex].coordinates;
//            else
//                crystal.unitCell.cartesianToFractional(crystal.atoms[atomIndex].coordinates, position);
                

            positions.push_back(position);
            //mIndicesOfSymmetryEquivalentAtoms[atomIndex].push_back(atomIndex);

            //----

            for (symmOpIndex = 0; symmOpIndex < nSymmOp; symmOpIndex++)
            {
                if (symmOpIndex == identityOperationIndex) // atoms from asymmetric unit were added at the beginning
                    continue;

                const SpaceGroupOperation &symmetryOperation = crystal.spaceGroup.getSpaceGroupOperation(symmOpIndex);

                if (crystal.xyzCoordinateSystem == structural_parameters_convention::XyzCoordinateSystem::fractional)
                    fractionalCoordinates = crystal.atoms[atomIndex].coordinates;
                else
                    crystal.unitCell.cartesianToFractional(crystal.atoms[atomIndex].coordinates, fractionalCoordinates);

                symmetryOperation.apply(fractionalCoordinates, position);

                //crystal.unitCell.fractionalToCartesian(fractionalCoordinatesTransformed, position);


                n = positions.size();
                foundNewPosition = true;
                for (i = 0; i < n; i++)
                    if (crystal_structure_utilities::equivalentPositions(positions[i], position, mCrystal.unitCell, shift))
                    {
                        equivalentPositionIndex = i;
                        foundNewPosition = false;
                        foundPositionOfAsymmUnitAtom = (equivalentPositionIndex == 0);
                        break;
                    }

                if (foundNewPosition)
                {
                    // shifts atom position to unit cell at origin (i.e. fractional coordinates between 0 and 1)
                    for (i = 0; i < 3; i++)
                    {
                        shift(i) = -math_utilities::roundInt(floor(position(i)));
                        position(i) = position(i) + shift(i);
                    }
                    translation.setToTranslation(Vector3<CrystallographicRational>(shift));
                    positions.push_back(position);
                    mAtoms.push_back(crystal.atoms[atomIndex]);
                    //mAtoms.back().atomicNumber = crystal.atoms[atomIndex].atomicNumber;
                    mAtoms.back().coordinates = position;

                    shiftedSymmetryOperation = translation*symmetryOperation;
                    shiftedSymmetryOperation.get(symmetryOperationAsString);
                    mAtoms.back().label += string("(") + symmetryOperationAsString + string(")");
                    mGeneratingOperations.resize(mGeneratingOperations.size() + 1);
                    mGeneratingOperations.back().push_back(shiftedSymmetryOperation);
                    mIndicesOfSymmetryEquivalentAtoms[atomIndex].push_back(mAtoms.size()-1);
                    mEquivalentAsymmUnitAtom.push_back(atomIndex);
                }
                else
                {
                    translation.setToTranslation(Vector3<CrystallographicRational>(shift));
                    shiftedSymmetryOperation = translation*symmetryOperation;
                    if(foundPositionOfAsymmUnitAtom)
                        mGeneratingOperations[atomIndex].push_back(shiftedSymmetryOperation);
                    else
                        mGeneratingOperations[nPreviouslyAddedAtoms + equivalentPositionIndex - 1].push_back(shiftedSymmetryOperation);
                }

            } // for(symmOpIndex=0;..

            nPreviouslyAddedAtoms += positions.size() - 1;

        } // for( atomIndex=0;..

    }

    Vector3d UnitCellContent::getAtomPositionFrac(
        const AtomID& atomID)
        const
    {
//        label = mCrystal.atoms[mEquivalentAsymmUnitAtom[atomID.atomIndex]].label;
        Vector3d result;
        SpaceGroupOperation t, s = mGeneratingOperations[atomID.atomIndex][0];
        t.setToTranslation(Vector3<CrystallographicRational>(atomID.unitCellPosition));
        (t * s).apply(mCrystal.atoms[mEquivalentAsymmUnitAtom[atomID.atomIndex]].coordinates, result);
        return result;
    }

    Vector3d UnitCellContent::getAtomPositionCart(
        const AtomID& atomID)
        const
    {
        Vector3d result;
        mCrystal.unitCell.fractionalToCartesian(getAtomPositionFrac(atomID), result);
        return result;
    }


    int UnitCellContent::indexOfSymmetryEquivalentAtomInCrystal(
        int atomIndex)
        const
    {
        return mEquivalentAsymmUnitAtom[atomIndex];
    }

    const Crystal &UnitCellContent::getCrystal()
        const
    {
        return mCrystal;
    }

    int UnitCellContent::nGeneratingOperations(int atomIdx) const
    {
        return mGeneratingOperations[atomIdx].size();
    }

    const SpaceGroupOperation &UnitCellContent::getGeneratingOperation(
        int atomIndex, 
        int operationIndex)
        const
    {
        return mGeneratingOperations[atomIndex][operationIndex];
    }


    int UnitCellContent::nAtoms()
        const
    {
        return mAtoms.size();
    }

    const AtomInCrystal &UnitCellContent::getAtom(
        int atomIndex)
        const
    {
        return mAtoms[atomIndex];
    }

    void UnitCellContent::findAtom(
        const std::string& atomLabel,
        const std::string& atomSymmetryCard,
        AtomID& atomID)
        const
    {
        if (!hasAtom(atomLabel, atomSymmetryCard, atomID))
        {
            string errorMessage = "cannot find atom with label '" + atomLabel +
                " and symmetry operation '" + atomSymmetryCard + "'";
            on_error::throwException(errorMessage, __FILE__, __LINE__);
        }
    }


	bool UnitCellContent::hasAtom(
		const std::string& atomLabel,
		const std::string& atomSymmetryCard,
		AtomID& atomID)
		const
	{
		string atomString = atomLabel + string("(") + atomSymmetryCard + string(")");
		return hasAtom(atomString, atomID);
	}

    void UnitCellContent::findAtom(
        const std::string& atomString,
        AtomID& atomID)
        const
    {
        if (!hasAtom(atomString, atomID))
        {
            string errorMessage = "cannot find atom '" + atomString + "'";
            on_error::throwException(errorMessage, __FILE__, __LINE__);
        }
    }

    bool UnitCellContent::hasAtom(
        const std::string &atomString,
        AtomID &atomID)
        const
    {
        string::size_type character_position;
        bool hasSymmetryOperation;
        string label, symmetryOperationAsStr;
        SpaceGroupOperation symmetryOperation;
        Vector3i shift;

        // if atomString contains bracket () at the end i.e. label(symmetryOperationAsStr)
        // then extract label and symmetryOperationAsStr

        if (atomString.empty())
            return false;

        //hasBracket = false;
        hasSymmetryOperation = false;

        if (atomString[atomString.size() - 1] == ')')
        {
            character_position = atomString.find_last_of('(');
            if (character_position != string::npos) // opening bracket found
            {

                label = atomString.substr(0, character_position);
                symmetryOperationAsStr = atomString.substr(character_position + 1, atomString.size() - character_position - 2);
                hasSymmetryOperation = SpaceGroupOperation::isSymmetryOperation(symmetryOperationAsStr);
            }
        }

        if (!hasSymmetryOperation)
            label = atomString;

        // look for atom with given label

        int atomInAsymmetricUnit, atomIndex, nAtoms = mCrystal.atoms.size();
        bool labelFound = false;

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
            if (mCrystal.atoms[atomIndex].label == label)
            {
                atomInAsymmetricUnit = atomIndex;
                labelFound = true;
                break;
            }

        if (!labelFound)
            return false;

        // the case with symmetry operation being identity

		
		bool symmetryOperationIsIdentity;

		if (!hasSymmetryOperation)
			symmetryOperationIsIdentity = true;
		else
			symmetryOperationIsIdentity = SpaceGroupOperation(symmetryOperationAsStr).isIdentity();

        if (symmetryOperationIsIdentity)
        {
            atomID.atomIndex = atomInAsymmetricUnit;
			Vector3d translation;
			mGeneratingOperations[atomInAsymmetricUnit][0].getTranslation(translation);
			atomID.unitCellPosition = -translation;// .set(0, 0, 0);
            return true;
        }

        // the case with symmetry operation appended to atom label
        // - search among atoms genrated by symmetry operation from the atom corresponding to the label

        int counter, i, nGeneratingOperations;
        nAtoms = mIndicesOfSymmetryEquivalentAtoms[atomInAsymmetricUnit].size();

        symmetryOperation.set(symmetryOperationAsStr);

        for (counter = 0; counter < nAtoms; counter++)
        {
            atomIndex = mIndicesOfSymmetryEquivalentAtoms[atomInAsymmetricUnit][counter];
            nGeneratingOperations = mGeneratingOperations[atomIndex].size();
            for (i = 0; i < nGeneratingOperations; i++)
            {
                if (mGeneratingOperations[atomIndex][i].isLatticeTranslationEquivalent(symmetryOperation, shift))
                {
                    // shift = symmetryOperation.translation - mGeneratingOperations[][].translation
                    atomID.atomIndex = atomIndex;
                    atomID.unitCellPosition = -shift;
                    return true;
                }
            }
        }

        // corresponding symmetry operation was not found - atomString invalid
        return false;
    }

    SpaceGroupOperation UnitCellContent::getTransformingOperation(
        const AtomID& atom1Id, 
        const AtomID& atom2Id) 
        const
    {
        SpaceGroupOperation result, translation1, translation2, atom_1_SymmOp, atom_2_SymmOp;

        translation1.setToTranslation(atom1Id.unitCellPosition.x, 1, atom1Id.unitCellPosition.y, 1, atom1Id.unitCellPosition.z, 1);
        translation2.setToTranslation(atom2Id.unitCellPosition.x, 2, atom2Id.unitCellPosition.y, 2, atom2Id.unitCellPosition.z, 2);

        atom_1_SymmOp = getGeneratingOperation(atom1Id.atomIndex, 0);
        atom_1_SymmOp = translation1 * atom_1_SymmOp;
        atom_2_SymmOp = getGeneratingOperation(atom2Id.atomIndex, 0);
        atom_2_SymmOp = translation2 * atom_2_SymmOp;

        atom_1_SymmOp.invert();
        result = atom_2_SymmOp * atom_1_SymmOp;

        return result;
    }



    std::string UnitCellContent::getAtomLabel(
        const AtomID &atomID)
        const
    {
        //if( atomID.unitCellPosition == Vector3i(0,0,0) )
        //    return mAtoms[atomID.atomIndex].label;
        //
        //string symmetryOperationAsString;
        //SpaceGroupOperation t,s = mGeneratingOperations[atomID.atomIndex][0];
        //t.setToTranslation(Vector3<boost::rational<int> >(atomID.unitCellPosition));
        //(t*s).get(symmetryOperationAsString);

        //return mAtoms[mEquivalentAsymmUnitAtom[atomID.atomIndex]].label+string("(")+symmetryOperationAsString+string(")");
        string label, symmetryOperation;
        interpreteAtomID(atomID, label, symmetryOperation);
        return label + string("(") + symmetryOperation + string(")");
    }

    void UnitCellContent::interpreteAtomID(
        const AtomID &atomID,
        std::string &label,
        std::string &symmetryOperation)
        const
    {
        label = mCrystal.atoms[mEquivalentAsymmUnitAtom[atomID.atomIndex]].label;

        //if (mEquivalentAsymmUnitAtom[atomID.atomIndex] == atomID.atomIndex) // is in asymmetric unit
        //    symmetryOperation = "X,Y,Z";
        //else
        //{
        //    SpaceGroupOperation t, s = mGeneratingOperations[atomID.atomIndex][0];
        //    t.setToTranslation(Vector3<CrystallographicRational>(atomID.unitCellPosition));
        //    (t*s).get(symmetryOperation);
        //}

        SpaceGroupOperation t, s = mGeneratingOperations[atomID.atomIndex][0];
        t.setToTranslation(Vector3<CrystallographicRational>(atomID.unitCellPosition));
        (t*s).get(symmetryOperation);

    }

}

