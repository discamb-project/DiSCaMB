#include "discamb/HC_Model/XdLocalCoordinateSystem.h"

#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"


#include <map>
#include <cmath>

using namespace std;


namespace discamb {

    XdLocalCoordinateSystem::XdLocalCoordinateSystem()
    {
        mIsCartesian = true;
        mIsDirection1AtomDummy = mIsDirection2Atom1Dummy = mIsDirection2Atom2Dummy =false;
        mRefAtom = 0;
        mAxesOrderRelatedHandedenessCorrection = 1.0;
        mHandedeness = 1.0;
        mCoordinate1 = 1;
        mCoordinate2 = 2;
        mCoordinate3 = 3;
        //can be 1 or -1
        mCoordinate1VectorMultiplier = mDirection12Multuiplier = 1.0;
    }


    XdLocalCoordinateSystem::~XdLocalCoordinateSystem()
    {
    }


    /*XD style coordinate system definition e.g.:
    N(1A)    C(13A)    Z  N(1A)    H(12A)   X   R
    for dummy atoms, e.g.
    O(1)     DUM6      Z  DUM5     H(15)    Y   R
    DUM1        0.642090    0.607910    0.625000
    ...
    DUM5        0.624450    0.333332    0.540222
    DUM6        0.250000    0.476543    0.500000
    followed the strig followed by the dummy atoms coordinates in the order they appear,
    i.e. in the example "O(1)     DUM6      Z  DUM5     H(15)    Y   R 0.250000    0.476543    0.500000 0.624450    0.333332    0.540222"
    */


    void XdLocalCoordinateSystem::set(
        const std::string &definition,
        const Crystal &c)
    {
        vector<string> words;
        string_utilities::split(definition, words);
        if (words.size() < 7)
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("'"), __FILE__, __LINE__);

        // DUM [1,3,4]
        int nDummy = 0;
        if (words[1].find("DUM") == 0)
        {
            nDummy++;
            words[1] = string("DUM") + string_utilities::convertToString(nDummy);
        }
        if (words[3].find("DUM") == 0)
        {
            nDummy++;
            words[3] = string("DUM") + string_utilities::convertToString(nDummy);
        }
        if (words[4].find("DUM") == 0)
        {
            nDummy++;
            words[4] = string("DUM") + string_utilities::convertToString(nDummy);
        }

        std::vector<Vector3d> dummyAtomsPositions(nDummy);
        if (words.size() != 7 + 3 * nDummy)
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("'"), __FILE__, __LINE__);

        for (int i = 0; i < nDummy; i++)
            for (int j = 0; j < 3; j++)
                dummyAtomsPositions[i][j] = atof(words[7 + 3 * i + j].c_str());

        set(words[0], words[1], words[3], words[4], words[2], words[5], (words[6][0] == 'R'), c, false, dummyAtomsPositions);
    }


    void XdLocalCoordinateSystem::set(
        const std::string &centralAtom,
        const std::string &direction_1_atom,
        const std::string &direction_2_atom_1,
        const std::string &direction_2_atom_2,
        const std::string &vector_1,
        const std::string &vector_2,
        bool rightHanded,
        const Crystal &crystal,
        bool dummyAtomIndexingFrom0,
        const std::vector<Vector3d> &dummyAtomFarctionalPositions)
    {
        mCrossProductLcs = shared_ptr<CrossProductLcs>(new CrossProductLcs(vector_1, vector_2, rightHanded));
        //map<string, int> directionIndexMap;
        //directionIndexMap["X"] = 1;
        //directionIndexMap["Y"] = 2;
        //directionIndexMap["Z"] = 3;
        //directionIndexMap["-X"] = -1;
        //directionIndexMap["-Y"] = -2;
        //directionIndexMap["-Z"] = -3;

        mIsCartesian = false;
        //mCoordinate1VectorMultiplier = directionIndexMap[vector_1] / abs(directionIndexMap[vector_1]);
        //mCoordinate1 = abs(directionIndexMap[vector_1]);
        //mDirection12Multuiplier = directionIndexMap[vector_2] / abs(directionIndexMap[vector_2]);
        //mCoordinate2 = abs(directionIndexMap[vector_2]);
        //mCoordinate3 = 6 - mCoordinate1 - mCoordinate2;

        vector<int> reference_atom;
        vector<SpaceGroupOperation> aux_sym_ops;

        getRefAtomIndex(centralAtom, reference_atom, aux_sym_ops, crystal);
        mRefAtom = reference_atom[0];

        processAtomLabel(direction_1_atom, crystal, dummyAtomIndexingFrom0, dummyAtomFarctionalPositions,
            mIsDirection1AtomDummy, mDummy_1_Position, mDirection1Atoms, mDirection1AtomsSymmOps);

        processAtomLabel(direction_2_atom_1, crystal, dummyAtomIndexingFrom0, dummyAtomFarctionalPositions,
            mIsDirection2Atom1Dummy, mDummy_2_1_Position, mDirection2Atoms1, mDirection2Atoms1SymmOps);

        processAtomLabel(direction_2_atom_2, crystal, dummyAtomIndexingFrom0, dummyAtomFarctionalPositions,
            mIsDirection2Atom2Dummy, mDummy_2_2_Position, mDirection2Atoms2, mDirection2Atoms2SymmOps);


        // 1 for R, -1 for L
        //rightHanded ? mHandedeness = 1 : mHandedeness = -1;

        //mAxesOrderRelatedHandedenessCorrection = 1 - 2 * ((2 + mCoordinate2 - mCoordinate1) % 3);
    }

    void XdLocalCoordinateSystem::processAtomLabel(
        const std::string &atomLabel,
        const Crystal &crystal,
        bool dummyAtomIndexingFrom0,
        const std::vector<Vector3d> &dummyAtomPositionsFract,
        bool &isDummyAtom,
        Vector3d &dummyAtomPositionCartesian,
        std::vector<int> &atomIndices,
        std::vector<SpaceGroupOperation> &symmetryOperations)
    {

        isDummyAtom = isDummyAtomLabel(atomLabel);

        if (isDummyAtom)
            crystal.unitCell.fractionalToCartesian(dummyAtomPositionsFract[dummyAtomIndex(atomLabel, dummyAtomIndexingFrom0)],
                dummyAtomPositionCartesian);
        else
            getRefAtomIndex(atomLabel, atomIndices, symmetryOperations, crystal);
    }


    void XdLocalCoordinateSystem::getRefAtomIndex(
        const std::string &label,
        std::vector<int> &atomIndices,
        std::vector<SpaceGroupOperation> &symmetryOperations,
        const Crystal &c)
    {
        int i, j, nAtomsInDefinition, nAtomsInCrystal = c.atoms.size();
        vector<string> atomsAsStrings;
        int unassignedAtom;
        bool atomAssigned;
        string_utilities::split(label, atomsAsStrings, ';');

        nAtomsInDefinition = atomsAsStrings.size();

        atomIndices.resize(nAtomsInDefinition);
        symmetryOperations.assign(nAtomsInDefinition, SpaceGroupOperation());

        unassignedAtom = -1;
        for (i = 0; i < nAtomsInDefinition; i++)
        {
            atomAssigned = false;
            for (j = 0; j < nAtomsInCrystal; j++)
                if (c.atoms[j].label == atomsAsStrings[i])
                {
                    atomIndices[i] = j;
                    atomAssigned = true;
                }

            if (!atomAssigned)
                unassignedAtom = static_cast<int>(i);
        }

        if (unassignedAtom > -1)
        {
            string errorMessage = string("unknown atom label '");
            errorMessage += atomsAsStrings[unassignedAtom];
            errorMessage += string("' in definition of local coordinate system for aspherical atom in Hansen-Coppens model.");
            on_error::throwException(errorMessage, __FILE__, __LINE__);
        }

    }

    bool XdLocalCoordinateSystem::isDummyAtomLabel(
        const std::string &label)
    {
        if (label.size() < 4)
            return false;
        return (string("DUM") == label.substr(0, 3));
    }

    int XdLocalCoordinateSystem::dummyAtomIndex(
        const std::string &label,
        bool dummyAtomIndexingFrom0)
    {
        int index = atoi(label.substr(3).c_str());
        if (!dummyAtomIndexingFrom0)
            --index;
        return index;
    }



    void XdLocalCoordinateSystem::calculate(
        Vector3d &x,
        Vector3d &y,
        Vector3d &z,
        const Crystal &crystal)
        const
    {
        Vector3d point_1_position, point_2_position, point_3_position, reference_atom_position, v1, v2;

        if (mIsCartesian)
        {
            x = Vector3d(1, 0, 0);
            y = Vector3d(0, 1, 0);
            z = Vector3d(0, 0, 1);
            return;
        }

        if (mIsDirection1AtomDummy)
            point_1_position = mDummy_1_Position;
        else
            point_1_position = averageAtomSetPosition(crystal, mDirection1Atoms, mDirection1AtomsSymmOps);

        if (mIsDirection2Atom1Dummy)
            point_2_position = mDummy_2_1_Position;
        else
            point_2_position = averageAtomSetPosition(crystal, mDirection2Atoms1, mDirection2Atoms1SymmOps);

        if (mIsDirection2Atom2Dummy)
            point_3_position = mDummy_2_2_Position;
        else
            point_3_position = averageAtomSetPosition(crystal, mDirection2Atoms2, mDirection2Atoms2SymmOps);


        if (crystal.xyzCoordinateSystem == structural_parameters_convention::XyzCoordinateSystem::fractional)
            crystal.unitCell.fractionalToCartesian(crystal.atoms[mRefAtom].coordinates, reference_atom_position);
        else
            reference_atom_position = crystal.atoms[mRefAtom].coordinates;


        v1 = point_1_position - reference_atom_position;

        v2 = point_3_position - point_2_position;

        mCrossProductLcs->calculate(v1, v2, x, y, z);
    }

    Vector3d XdLocalCoordinateSystem::averageAtomSetPosition(
        const Crystal &crystal,
        const std::vector<int> &atoms,
        const std::vector<SpaceGroupOperation> &symmetryOperations)
    {
        Vector3d atom_fractional, position_fractional, average_position_fractional, result;
        int atom_index, nAtoms = atoms.size();

        for (atom_index = 0; atom_index < nAtoms; atom_index++)
        {
            // get position in fractional coordinates
            if (crystal.xyzCoordinateSystem == structural_parameters_convention::XyzCoordinateSystem::fractional)
                atom_fractional = crystal.atoms[atoms[atom_index]].coordinates;
            else
                crystal.unitCell.cartesianToFractional(crystal.atoms[atoms[atom_index]].coordinates, atom_fractional);
            symmetryOperations[atom_index].apply(atom_fractional, position_fractional);
            average_position_fractional += position_fractional;
        }

        crystal.unitCell.fractionalToCartesian(average_position_fractional, result);

        return result;
    }

} // namespace discamb



//namespace discamb {
//
//XdLocalCoordinateSystem::XdLocalCoordinateSystem()
//{
//    mIsCartesian = true;
//}
//
//
//XdLocalCoordinateSystem::~XdLocalCoordinateSystem()
//{
//}
//   
//
///*XD style coordinate system definition e.g.:
//N(1A)    C(13A)    Z  N(1A)    H(12A)   X   R
//for dummy atoms, e.g.
//O(1)     DUM6      Z  DUM5     H(15)    Y   R
//DUM1        0.642090    0.607910    0.625000
//...
//DUM5        0.624450    0.333332    0.540222
//DUM6        0.250000    0.476543    0.500000
//followed the strig followed by the dummy atoms coordinates in the order they appear,
//i.e. in the example "O(1)     DUM6      Z  DUM5     H(15)    Y   R 0.250000    0.476543    0.500000 0.624450    0.333332    0.540222"
//*/
//
//
//void XdLocalCoordinateSystem::set(
//	const std::string &definition,
//	const Crystal &c)
//{
//	vector<string> words;
//	string_utilities::split(definition, words);
//	if (words.size() < 7)
//		on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("'"), __FILE__, __LINE__);
//
//	// DUM [1,3,4]
//	int nDummy=0;
//	if (words[1].find("DUM") == 0)
//	{
//		nDummy++;
//		words[1] = string("DUM") + string_utilities::convertToString(nDummy);
//	}
//	if (words[3].find("DUM") == 0)
//	{
//		nDummy++;
//		words[3] = string("DUM") + string_utilities::convertToString(nDummy);
//	}
//	if (words[4].find("DUM") == 0)
//	{
//		nDummy++;
//		words[4] = string("DUM") + string_utilities::convertToString(nDummy);
//	}
//
//	std::vector<Vector3d> dummyAtomsPositions(nDummy);
//	if(words.size()!=7+3*nDummy)
//		on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("'"), __FILE__, __LINE__);
//
//	for (int i = 0; i < nDummy; i++)
//		for (int j = 0; j < 3; j++)
//			dummyAtomsPositions[i][j] = atof(words[7 + 3 * i + j].c_str());
//
//	set(words[0], words[1], words[3], words[4], words[2], words[5], (words[6][0]=='R'), c, false, dummyAtomsPositions);
//}
//
//
//void XdLocalCoordinateSystem::set(
//const std::string &centralAtom,
//const std::string &direction_1_atom,
//const std::string &direction_2_atom_1,
//const std::string &direction_2_atom_2,
//const std::string &vector_1,
//const std::string &vector_2,
//bool rightHanded,
//const Crystal &crystal,
//bool dummyAtomIndexingFrom0,
//const std::vector<Vector3d> &dummyAtomFarctionalPositions)
//{
//    map<string,int> directionIndexMap;
//    directionIndexMap["X"]=1;
//    directionIndexMap["Y"]=2;
//    directionIndexMap["Z"]=3;
//    directionIndexMap["-X"]=-1;
//    directionIndexMap["-Y"]=-2;
//    directionIndexMap["-Z"]=-3;
//
//    mIsCartesian = false;
//    mCoordinate1VectorMultiplier = directionIndexMap[vector_1]/abs(directionIndexMap[vector_1]);
//    mCoordinate1 = abs(directionIndexMap[vector_1]);
//    mDirection12Multuiplier = directionIndexMap[vector_2] / abs(directionIndexMap[vector_2]);
//    mCoordinate2 = abs(directionIndexMap[vector_2]);
//    mCoordinate3 = 6 - mCoordinate1 - mCoordinate2;
//
//    vector<int> reference_atom;
//    vector<SpaceGroupOperation> aux_sym_ops;
//
//    getRefAtomIndex(centralAtom,reference_atom,aux_sym_ops,crystal);
//    mRefAtom = reference_atom[0];
//
//    processAtomLabel(direction_1_atom,crystal,dummyAtomIndexingFrom0,dummyAtomFarctionalPositions,
//                     mIsDirection1AtomDummy,mDummy_1_Position,mDirection1Atoms,mDirection1AtomsSymmOps);
//
//    processAtomLabel(direction_2_atom_1,crystal,dummyAtomIndexingFrom0,dummyAtomFarctionalPositions,
//                     mIsDirection2Atom1Dummy,mDummy_2_1_Position,mDirection2Atoms1,mDirection2Atoms1SymmOps);
//
//    processAtomLabel(direction_2_atom_2,crystal,dummyAtomIndexingFrom0,dummyAtomFarctionalPositions,
//                     mIsDirection2Atom2Dummy,mDummy_2_2_Position,mDirection2Atoms2,mDirection2Atoms2SymmOps);
//
//        
//    // 1 for R, -1 for L
//    rightHanded ? mHandedeness = 1 : mHandedeness = -1;  
//
//    mAxesOrderRelatedHandedenessCorrection = 1 - 2*( (2+mCoordinate2-mCoordinate1) %3);
//}
//
//void XdLocalCoordinateSystem::processAtomLabel(
//const std::string &atomLabel,
//const Crystal &crystal,
//bool dummyAtomIndexingFrom0,
//const std::vector<Vector3d> &dummyAtomPositionsFract,
//bool &isDummyAtom,
//Vector3d &dummyAtomPositionCartesian,
//std::vector<int> &atomIndices,
//std::vector<SpaceGroupOperation> &symmetryOperations)
//{
//
//    isDummyAtom =  isDummyAtomLabel(atomLabel);
//
//    if(isDummyAtom)
//        crystal.unitCell.fractionalToCartesian(dummyAtomPositionsFract[dummyAtomIndex(atomLabel,dummyAtomIndexingFrom0)],
//                                               dummyAtomPositionCartesian);
//    else
//        getRefAtomIndex(atomLabel,atomIndices,symmetryOperations,crystal);
//}
//
//
//void XdLocalCoordinateSystem::getRefAtomIndex(
//const std::string &label,
//std::vector<int> &atomIndices,
//std::vector<SpaceGroupOperation> &symmetryOperations,
//const Crystal &c)
//{
//    int i,j,nAtomsInDefinition,nAtomsInCrystal=c.atoms.size();
//    vector<string> atomsAsStrings;
//    int unassignedAtom;
//    bool atomAssigned;
//    string_utilities::split(label,atomsAsStrings,';');
//
//    nAtomsInDefinition = atomsAsStrings.size();
//
//    atomIndices.resize(nAtomsInDefinition);
//    symmetryOperations.assign(nAtomsInDefinition,SpaceGroupOperation());
//
//    unassignedAtom = -1;
//    for(i=0;i<nAtomsInDefinition;i++)
//    {
//        atomAssigned = false;
//        for(j=0;j<nAtomsInCrystal;j++)
//            if (c.atoms[j].label == atomsAsStrings[i])
//            {
//                atomIndices[i] = j;
//                atomAssigned = true;
//            }
//
//        if(!atomAssigned)
//            unassignedAtom = i;
//    }
//
//    if( unassignedAtom > -1 )
//    {
//        string errorMessage = string("unknown atom label '");
//        errorMessage += atomsAsStrings[unassignedAtom];
//        errorMessage += string("' in definition of local coordinate system for aspherical atom in Hansen-Coppens model.");
//        on_error::throwException(errorMessage, __FILE__, __LINE__);
//    }
//    
//}
//
//bool XdLocalCoordinateSystem::isDummyAtomLabel(
//const std::string &label)
//{
//    if(label.size()<4)
//        return false;
//    return (string("DUM")==label.substr(0,3));
//}
//
//int XdLocalCoordinateSystem::dummyAtomIndex(
//const std::string &label,
//bool dummyAtomIndexingFrom0)
//{
//    int index = atoi(label.substr(3).c_str());
//    if(!dummyAtomIndexingFrom0)
//        --index;
//    return index;
//}
//
//
//
//void XdLocalCoordinateSystem::calculate(
//Vector3d &x,
//Vector3d &y,
//Vector3d &z,
//const Crystal &crystal) 
//const
//{
//    Vector3d point_1_position,point_2_position,point_3_position,reference_atom_position;
//    Vector3d direction[3],direction_12;
//    Vector3d *axes[]={&x,&y,&z};
//    
//    if (mIsCartesian)
//    {
//        x = Vector3d(1, 0, 0);
//        y = Vector3d(0, 1, 0);
//        z = Vector3d(0, 0, 1);
//        return;
//    }
//
//    if(mIsDirection1AtomDummy)
//        point_1_position = mDummy_1_Position;
//    else
//        point_1_position = averageAtomSetPosition(crystal,mDirection1Atoms,mDirection1AtomsSymmOps);
//        
//    //cout << "point_1_position\n";
//    //cout << point_1_position[0] << " " << point_1_position[1] << " " << point_1_position[2] << endl;
//    if(mIsDirection2Atom1Dummy)
//        point_2_position = mDummy_2_1_Position;
//    else
//        point_2_position = averageAtomSetPosition(crystal,mDirection2Atoms1,mDirection2Atoms1SymmOps);
//    //cout << "point_2_position\n";
//    //cout << point_2_position[0] << " " << point_2_position[1] << " " << point_2_position[2] << endl;
//    
//    if(mIsDirection2Atom2Dummy)
//        point_3_position = mDummy_2_2_Position;
//    else
//        point_3_position = averageAtomSetPosition(crystal,mDirection2Atoms2,mDirection2Atoms2SymmOps);
//
//    //cout << "point_3_position\n";
//    //cout << point_3_position[0] << " " << point_3_position[1] << " " << point_3_position[2] << endl;
//
//    if(crystal.xyzCoordinateSystem==structural_parameters_convention::fractional)
//        crystal.unitCell.fractionalToCartesian(crystal.atoms[mRefAtom].coordinates, reference_atom_position);
//    else
//        reference_atom_position = crystal.atoms[mRefAtom].coordinates;
//
//    //cout << "reference_atom_position\n";
//    //cout << reference_atom_position[0] << " " << reference_atom_position[1] << " " << reference_atom_position[2] << endl;
//
//    direction[0] = mCoordinate1VectorMultiplier * (point_1_position - reference_atom_position);
//    direction[0] /= sqrt(direction[0]*direction[0]);
//
//    //cout << "direction[0]\n";
//    //cout << direction[0][0] << " " << direction[0][1] << " " << direction[0][2] << endl;
//     
//
//    direction_12 = mDirection12Multuiplier * (point_3_position - point_2_position);
//    
//    //cout << "direction_12\n";
//    //cout << direction_12[0] << " " << direction_12[1] << " " << direction_12[2] << endl;
//
//    direction[2] = mAxesOrderRelatedHandedenessCorrection * mHandedeness * cross_product(direction[0],direction_12);
//    direction[2] /= sqrt(direction[2]*direction[2]);
//
//    //cout << "direction[2]\n";
//    //cout << direction[2][0] << " " << direction[2][1] << " " << direction[2][2] << endl;
//
//
//    direction[1] = mAxesOrderRelatedHandedenessCorrection * mHandedeness * cross_product(direction[2], direction[0]);
//
//    //cout << "direction[1]\n";
//    //cout << direction[1][0] << " " << direction[1][1] << " " << direction[1][2] << endl;
//
//
//    *axes[mCoordinate1-1] = direction[0];
//    *axes[mCoordinate2-1] = direction[1];
//    *axes[mCoordinate3-1] = direction[2];   
//}
//
//Vector3d XdLocalCoordinateSystem::averageAtomSetPosition(
//const Crystal &crystal,
//const std::vector<int> &atoms,
//const std::vector<SpaceGroupOperation> &symmetryOperations)
//{
//    Vector3d atom_fractional, position_fractional,average_position_fractional,result;
//    int atom_index,nAtoms = atoms.size();
//
//    for(atom_index=0;atom_index<nAtoms;atom_index++)
//    {
//        // get position in fractional coordinates
//        if(crystal.xyzCoordinateSystem == structural_parameters_convention::fractional)
//            atom_fractional = crystal.atoms[atoms[atom_index]].coordinates;
//        else
//            crystal.unitCell.cartesianToFractional(crystal.atoms[atoms[atom_index]].coordinates, atom_fractional);
//        symmetryOperations[atom_index].apply(atom_fractional, position_fractional);
//        average_position_fractional += position_fractional;
//    }
//
//    crystal.unitCell.fractionalToCartesian(average_position_fractional,result);
//
//    return result;
//}
//
//} // namespace discamb
//
