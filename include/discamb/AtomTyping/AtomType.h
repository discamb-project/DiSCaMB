#pragma once

#include "discamb/BasicUtilities/Tribool.h"

#include <vector>
#include <string>
#include <deque>
#include <tuple>
#include <set>
#include <optional>


namespace discamb {

    /**
    * \defgroup AtomTyping AtomTyping
    \brief Atom type assignment.
    * @{
    */




    struct AtomRingInfo
    {
        void set(const std::string &planarRingInfo, const std::string &ring3Info, const std::string &ring4info);

        Tribool inRing = Tribool::Undefined;
        Tribool in3Ring = Tribool::Undefined;
        Tribool in4Ring = Tribool::Undefined;
        /** undefined if negative */
        int n3rings = -1;
        int n4rings = -1;
        bool inAnyAdditionalRing = false;
        /** ring size and label*/
        std::vector<std::pair<int, std::string> > labeledContainingRings;
        std::vector<std::pair<int, std::string> > labeledNonContainingRings;
        std::vector<int> nonLabeledContainingRings;
        std::vector<int> nonLabeledNonContainingRings;
    };


    struct AtomDescriptors
    {
        std::string label;
        
        std::set<int> atomic_number_range;

        bool anyAtomicNumber = false;
        Tribool planar = Tribool::Undefined;
        std::multiset<int> neighborsAtomicNumbers;
        std::vector<std::set<int> > neighborsAtomicNumberRanges;
        int nNeighbours;
        
        AtomRingInfo ringInfo;
        /**
        false if the atom can have more neighbors than these specified in neighborsAtomicNumbers
        and neighborsAtomicNumberRanges
        */
        bool fixedNumberOfNeighbors;
        
    };


    struct AtomTypeLCS
    {
        // MOST_ORTHOGONAL_ATOM is used when one wants to define axis along the direction 
        // from the central atoms to atom ('directing atom')
        // that is most orthogonal to the other axis
        // it requires that the other axis is defined as LABELED_ATOM or ATOM_LIST
        // in the case of central atom with multiple neighbours the directing atom is chosen 
        // from its 1-st neighbours, in the case of central atom with 1 neighbour,
        // the directing atom is chosen from the central atom second neighbours
        enum class LCS_AxisType {
            LABELED_ATOM, NON_LABELED_ATOM,
            NOT_OF_ATOMIC_NUMBER, LABELED_RING,
            NON_LABELED_RING, ATOM_LIST,
            ANY_ATOM, ANY_ORTHOGONAL, AVERAGE_DIRECTION, NOT_SET,
            MOST_ORTHOGONAL_ATOM
        };

        LCS_AxisType lcs_axis_type_1 = LCS_AxisType::ANY_ORTHOGONAL;
        LCS_AxisType lcs_axis_type_2 = LCS_AxisType::ANY_ORTHOGONAL;
        // X-0,Y-1,Z-2
        int lcs_coordinate_1 = 0;
        int lcs_coordinate_2 = 1;
        bool isR = true;
        bool automaticallyDefined = false;
        /**
           Reference points definition:

           In the case of LABELED_ATOM , NON_LABELED_ATOM , NOT_OF_ATOMIC_NUMBER , LABELED_RING , NON_LABELED_RING - single number:

              LABELED_ATOM - index refering to atom in the array atoms
              NON_LABELED_ATOM - atomic number
              NOT_OF_ATOMIC_NUMBER - atomic number
              LABELED_RING - index in the array ringLabels
              NON_LABELED_RING - size of the ring.

           In addition in the case of NON_LABELED_ATOM, NOT_OF_ATOMIC_NUMBER, ANY_ATOM - in braces atom
           to which the atom is bonded or list of atoms to which the atom can be bonded

           In the case of ATOM_LIST
               - component 0 - number (n) of labeled atoms in the list
               - (if there are any labeled atoms) components 1 to n - indices of the corresponding labeled atoms
               - components n+1 and above - atomic numbers of non-labeled atoms

           Empty in the case of ANY_ATOM, ANY_ORTHOGONAL and MOST_ORTHOGONAL_ATOM
        */

        std::vector<int> lcs_axis_1_definition, lcs_axis_2_definition;
        // it is far from being complete
        bool checkCorrectness(std::string &errorMessage) const;
    };

    struct AtomType
    {
        std::vector<std::string> commentLines;
        std::vector<AtomDescriptors> atoms;
        /**
        connectivity matrix
        */
        std::vector<std::vector<int> > connectivity;
        /**
        labels of labeled rings
        */
        std::vector<std::string> ringLabels;

        /**
        sizes of labeled rings (in the same order as ringLabels)
        */
        std::vector<int> ringSizes;
        
        std::string symmetry;
        std::string entryCreationDate, definitionModificationDate, parameterModificationDate;
        std::string id;
        int numberOfInstances = 0;

        /**
    	atom_1, atom_2, atom_3
    	when checking chirality	three vectors are created - from central atom to i-th atom:
    	v1 = atom_1 - central_atom
    	v2 = atom_2 - central_atom
    	v3 = atom_3 - central_atom
    	the angle between cross product of the first two vectors and the third vector is below 90
    	*/
		std::vector<int> chirality;
        //-------------------------------------
        // local coordinate system
        //-------------------------------------

        
        void setLocalCoordinateSystem(const std::string &lcsAsString);
        bool setDefaultLocalCoordinateSystem();
        AtomTypeLCS localCoordinateSystem;
        // move the 4 below later to private (when class is made)
        int getCoordinate(const std::string &coordinate) const;
        void setReferencePoint(const std::string &refPoint, std::vector<int> &definition, AtomTypeLCS::LCS_AxisType &definitionType) const;
        // returns false if there are no such two atoms
        bool chooseNamedNeighboursForLcs(int& neighbour1, int& neighbour2);
        void unnamedNeighboursAtomicNumbers(int atomIdx, std::vector<int>& unnamedNeighborsZ);
        void transformUnnamedAtomsToNamed(int atomicNumber, const AtomDescriptors& defaultDescriptors);
        std::optional<int> getAtomIndex(const std::string &label) const;
        static void sphericalTypes(std::vector<AtomType> &types);
        void transformToCanonicalForm();
    };

    /**@}*/


}
