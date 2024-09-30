#pragma once

#include "discamb/BasicChemistry/ChemicalElement.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCellContent.h"


namespace discamb {

    /**
    * \addtogroup CrystalStructure
    * @{
    */


    namespace crystal_structure_utilities {

        void stringToAtomList(const std::string& str, std::vector<std::pair<std::string, std::string> >& atomList, char separator);


        /**
        Cartesian position of atom in crystal
        */
        Vector3d atomPosition(int atomIdx, const SpaceGroupOperation& symmOp, const Crystal crystal);

        /**
        Cartesian position of atom in crystal
        */
        Vector3d atomPosition(const std::string &atomLabel, const SpaceGroupOperation& symmOp, const Crystal crystal);


        double u_eq(const UnitCell& uc, const std::vector<double>& uij_cif);
        void u_eq(const Crystal& crystal, std::vector<double>& u_eqivalent);

        /*void translateAtomTo01(const Crystal &crystal, int atomIdx, SpaceGroupOperation &operation,
                               SpaceGroupOperation &operationAfterTranslation, Vector3i &translation);*/
        /**
        "C1,-X,Y+1/2,-Z" => {"C1", "-X,Y+1/2,-Z"}
        if there is no symm op X,Y,Z, will be given
        returns false if encounters problem, true otherwise
        */
        bool splitIntoAtomAndSymmOp(const std::string& s, std::string& atomLabel, std::string& symmOp);
        void splitIntoAtomAndSymmOp(const std::string& s, std::string& atomLabel, std::string& symmOp, bool throwException);
        void splitIntoAtomAndSymmOp(
            const std::vector<std::string>& s,
            std::vector < std::pair<std::string, std::string> >& atoms,
            bool throwException);

        // if labelAndSymmOperation.second is empty then it is trated as identity (i.e. as X,Y,Z)
        void convertToXyzAndElementList(const Crystal& crystal, const std::vector<std::pair<std::string, std::string> >& labelAndSymmOperation,
                                   std::vector<ChemicalElement> &symbol, std::vector<Vector3d> &position);

        void convertToXyzAndElementList(const Crystal& crystal, const UnitCellContent &ucContent, std::vector<UnitCellContent::AtomID>& atoms,
            std::vector<ChemicalElement>& element, std::vector<Vector3d>& position);

        void convertAtomList(
            const UnitCellContent& ucContent, 
            const std::vector<UnitCellContent::AtomID>& atomList, 
            std::vector < std::pair < std::string, std::string > > & atomListConverted);

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector < std::pair < std::string, std::string > >& atomList,
            std::vector<UnitCellContent::AtomID>& atomsListConverted);

        //

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector<UnitCellContent::AtomID>& atomList,
            std::vector < std::pair < int, std::string > >& atomListConverted);

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector < std::pair < int, std::string > >& atomList,
            std::vector<UnitCellContent::AtomID>& atomListConverted);

        void convertAtomList(
            const Crystal & crystal,
            const std::vector < std::pair < int, std::string > >& atomList,
            std::vector < std::pair < std::string, std::string > >& atomListConverted);

        void convertAtomList(
            const Crystal& crystal,
            const std::vector < std::pair < std::string, std::string > >& atomList,
            std::vector < std::pair < int, std::string > >& atomListConverted); 



        double interatomicDistance(
            const Crystal& crystal, 
            const std::string& atomLabel1, 
            const std::string& symmOp1, 
            const std::string& atomLabel2, 
            const std::string& symmOp2);


        bool findAtomSymmetry(const Crystal &crystal, 
                              int atomIdx, 
                              std::vector<std::vector<SpaceGroupOperation> > &pointGroups,
                              double threshold);

        void atomicNumbers(const Crystal &crystal, std::vector<int> &atomicNumbers);

        // returns true if v1 = v2 + shift (shift is vetor with integer components - a lattice vetor)
        bool equivalentPositions(const Vector3d& v1, const Vector3d& v2, const UnitCell& unitCell, Vector3i& shift);


        // r = rInUnitCell + unitCellOrigin, 0 <= rInUnitCell[i] < 1 
        void extractLatticeVector(const Vector3d &r, Vector3d &rInUnitCell, Vector3i &unitCellOrigin);

        void symmetryOperationCartesian(const SpaceGroupOperation &symmOp, const UnitCell & unitCell, Matrix3d &rotationCartesian, Vector3d &translationCartesian);

        void getLatticeVectors(const UnitCell& unitCell, Vector3d& a, Vector3d& b, Vector3d& c);

        double s12AdpSimilarityIndex(const UnitCell& uc, const std::vector<double>& cifAdp1, const std::vector<double>& cifAdp2);
        double s12AdpSimilarityIndex(const std::vector<double>& cartAdp1, const std::vector<double>& cartAdp2);

        double overlappingAdpSimilarityIndex(const std::vector<double>& cartAdp1, const std::vector<double>& cartAdp2);
        
        

        double msdCorrelation(const std::vector<double>& cartAdp1, const std::vector<double>& cartAdp2);
        double rmsdCorrelation(const std::vector<double>& cartAdp1, const std::vector<double>& cartAdp2);

        bool isSuperLatticeNode(const Vector3i& node, const Vector3i& supeLatticeNode, const std::vector<Vector3i>& superLatticeBase);
        void generate_hkl(const UnitCell& uc, double resolution);
    }
    /**@}*/
}