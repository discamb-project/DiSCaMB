#ifndef _DISCAMB_HC_MODEL_XDLOCALCOORDINATESYSTEM_H_
#define _DISCAMB_HC_MODEL_XDLOCALCOORDINATESYSTEM_H_


//#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"
#include "discamb/MathUtilities/CrossProductLcs.h"

#include <memory>

namespace discamb {

/** \ingroup HC_Model */

/** \brief Local coordinate system defined similarily as in XD. */
class XdLocalCoordinateSystem: public LocalCoordinateSystemInCrystal
{
public:
    /** \brief Constructs the coordinate system as the Cartesian one.*/
   XdLocalCoordinateSystem();
   virtual ~XdLocalCoordinateSystem();

   /**
   \brief Define the local coordinate system. 
   
   The first four parameters are used to define four reference
   points. Each of them can be defined in one of the following ways:
   - by atom label (e.g. C1) - corresponds to the position of the corresponding atom 
   - multiple atom labels separater with semicolon ';' (with no spaces e.g. C1;C2;N3) - corresponds to mean ("average")
     position of the corresponding atoms 
   - by dummy atoms - the label should take the following form: DUMn where n is an index of the dummy atom e.g. DUM3 
     corresponds to 3-rd dummy atom in the input vector of dummy atoms (\p dummyAtomFarctionalPositions) or to the 4-th one 
     when the parameter \p dummyAtomIndexingFrom0 is true 

   The labels used have to correspond to atomic labels in the structure specified with the parameter \p c.
   The reference points define 2 directions described by parameters \p coordinate_1 and \p coordinate_2 which can take one 
   of the following values "X", "Y", "Z", "-X", "-Y", "-Z". 
   The direction from the point defined by \p centralAtom to \p direction_1_atom defines the direction (\f$ \mathbf{d}_1 \f$)
   indicated by  the parameter \p coordinate_1. The direction from point defined by \p direction_2_atom_1 to the the one given by 
   \p direction_2_atom_2 defines auxilary vector (\f$ \mathbf{d}_{12} \f$). The direction corresponding to \p coordinate_2 
   (\f$ \mathbf{d}_2 \f$) lays on the same plane as \f$ \mathbf{d}_1 \f$ and \f$ \mathbf{d}_{12} \f$ and it is perpendicular 
   to \f$ \mathbf{d}_1 \f$.
   There are two vectors fulfilling this conditions and the one which is more similar to the direction \f$ \mathbf{d}_{12} \f$ is chosen 
   (i.e. the one which has positive dot product with \f$ \mathbf{d}_{12} \f$). It is given by the following expression:
   \f[
   \mathbf{d}_2 = N ( \mathbf{d}_1 \times \mathbf{d}_{12} ) \times \mathbf{d}_1 
   \f]
   where N is normalization factor. The third remaining direction of coordinate system 
   is assigned according to the definition of the two already defined directions and handedness (specified by the parameter 
   \p rightHanded ).
      
   */ 
    void set(const std::string &centralAtom, const std::string &direction_1_atom, const std::string &direction_2_atom_1,
             const std::string &direction_2_atom_2, const std::string &coordinate_1, const std::string &coordinate_2,
             bool rightHanded, const Crystal &c, bool dummyAtomIndexingFrom0 = true,
             const std::vector<Vector3d> &dummyAtomFarctionalPositions = std::vector<Vector3d>());

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
	if the same dummy atom appears twice its coordinates have to be provided twice
	*/
	virtual void set(const std::string &definition, const Crystal &c);
	using LocalCoordinateSystemInCrystal::calculate;
    /** \brief Calculates coordinate system vectors for particular geometry of crystal structure. 
    
       The orthonormal vectors (\p x , \p y , \p z ) of coordinate system corresponds to geometry
       of the structure provided in the argument \p c .*/
    virtual void calculate(Vector3d &x, Vector3d &y, Vector3d &z, const Crystal &c) const;


private:    
    std::shared_ptr<CrossProductLcs> mCrossProductLcs;
    bool mIsCartesian;
    bool mIsDirection1AtomDummy,mIsDirection2Atom1Dummy,mIsDirection2Atom2Dummy;
    std::vector<int> mDirection1Atoms,mDirection2Atoms1,mDirection2Atoms2;
    int mRefAtom; 
    // 1 or -1
    double mAxesOrderRelatedHandedenessCorrection;
    // 1 for R, -1 for L
    double mHandedeness;  
    // negative if no symmetry operation should be applied
    //std::vector<int> mDirection1AtomSymmOpIndices,mDirection2Atom1SymmOpIndices,mDirection2Atom2SymmOpIndices;
    std::vector<SpaceGroupOperation> mDirection1AtomsSymmOps,mDirection2Atoms1SymmOps,mDirection2Atoms2SymmOps;
    // 1 for X, 2 for Y, 3 for Z
    int mCoordinate1,mCoordinate2,mCoordinate3;
    //can be 1 or -1
    double mCoordinate1VectorMultiplier, mDirection12Multuiplier; 
    Vector3d mDummy_1_Position,mDummy_2_1_Position,mDummy_2_2_Position;
    static bool isDummyAtomLabel(const std::string &label);
    // position of dummy atom (Cartesian coordinates);
    static int dummyAtomIndex(const std::string &label,bool dummyAtomIndexingFrom0);
    /** symmetry operation index is negative if no symmetry operation is specified
    \todo enable finding symmetry operation and check for label correctness
    */
    static void getRefAtomIndex(const std::string &label,std::vector<int> &atomIndices,
                                std::vector<SpaceGroupOperation> &symmetryOperations,const Crystal &c);
    static void processAtomLabel(const std::string &atomLabel,const Crystal &crystal,bool dummyAtomIndexingFrom0,
                                 const std::vector<Vector3d> &dummyAtomPositionsFract,bool &isDummyAtom,
                                 Vector3d &dummyAtomPosition,std::vector<int> &atomIndices,
                                 std::vector<SpaceGroupOperation> &symmetryOperations);
    static Vector3d averageAtomSetPosition(const Crystal &c,const std::vector<int> &atoms,
                                           const std::vector<SpaceGroupOperation> &symmetryOperations);
    
};

}

/**@}*/

#endif /*_DISCAMB_HC_MODEL_XDLOCALCOORDINATESYSTEM_H_*/


