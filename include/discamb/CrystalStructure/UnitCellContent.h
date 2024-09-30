#pragma once


#include "discamb/CrystalStructure/Crystal.h"

/*
atoms from asymmetric unit are at the beginning

requirements:
list of all atoms in [0,1) incl. reference atom in symm unit and the symm operation which makes the ref. atom the given atom 
pointers from assym unit atom to all symmetry equivalent in [0,1)
point group of each atom in [0,1) (?)

*/

namespace discamb {

    /**
     * \addtogroup CrystalStructure
     * @{
     */


    class UnitCellContent
    {
    public:
        UnitCellContent();
        UnitCellContent(const Crystal &crystal);
        ~UnitCellContent();

        void set(const Crystal &crystal);
        const Crystal &getCrystal() const;

        struct AtomID
        {
            //        AtomID(){atomIndex=0;unitCellPosition.set(0,0,0);}
            bool operator< (const AtomID& id) const
            {
                return std::pair<int, Vector3i>(atomIndex, unitCellPosition) < 
                       std::pair<int, Vector3i>(id.atomIndex, id.unitCellPosition);
            };
            int atomIndex;
            Vector3i unitCellPosition;
            AtomID() :atomIndex(0), unitCellPosition(Vector3i(0, 0, 0)) {}
            AtomID(int i, const Vector3i &ucp = Vector3i(0, 0, 0)) :atomIndex(i), unitCellPosition(ucp) {}
        };
        
        int nAtoms() const;
        int nGeneratingOperations(int atomIdx) const;
        const SpaceGroupOperation &getGeneratingOperation(int atomIndex, int operationIndex) const;
        // returns transforation which S transforms atom1 into atom 2
        // S at_1 = at_2
        // atom_1_ucp - atom 1 unit cell position
        SpaceGroupOperation getTransformingOperation(const AtomID& atomId1, const AtomID& atomId2) const;
        /**
        \todo U_ij not transformed
        atoms from asymmetric unit first
        */
        const AtomInCrystal &getAtom(int atomIndex) const;
        
		bool hasAtom(const std::string& atomLabel, const std::string &atomSymmetryCard, AtomID& atomID) const;
        bool hasAtom(const std::string &atomString, AtomID &atomID) const;
        void findAtom(const std::string& atomLabel, const std::string& atomSymmetryCard, AtomID& atomID) const;
        void findAtom(const std::string& atomString, AtomID& atomID) const;
        void interpreteAtomID(const AtomID &atomID, std::string &label, std::string &symmetryOperation) const;
        std::string getAtomLabel(const AtomID &atomID) const;
        Vector3d getAtomPositionFrac(const AtomID& atomID) const;
        Vector3d getAtomPositionCart(const AtomID& atomID) const;
        int indexOfSymmetryEquivalentAtomInCrystal(int atomIndex) const;

    private:
        std::vector<AtomInCrystal> mAtoms;
        std::vector<int> mEquivalentAsymmUnitAtom;

        Crystal mCrystal;
        std::vector<std::vector<int> > mIndicesOfSymmetryEquivalentAtoms;
        std::vector<std::vector<SpaceGroupOperation> > mGeneratingOperations;
    };

    inline bool operator==(const UnitCellContent::AtomID &id1, const UnitCellContent::AtomID &id2)
    {
        return ((id1.atomIndex == id2.atomIndex) && (id1.unitCellPosition == id2.unitCellPosition));
    }
    /**@}*/
}
