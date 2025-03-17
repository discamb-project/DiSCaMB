#ifndef DISCAMB_STRUCTURALPROPERTIES_STRUCTURAL_PROPERTIES_H_
#define DISCAMB_STRUCTURALPROPERTIES_STRUCTURAL_PROPERTIES_H_

#include "discamb/BasicUtilities/Tribool.h"
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/CrystalStructure/CrystalVarianceCovarianceMatrix.h"

#include <optional>

namespace discamb {

    /**
    * \defgroup StructuralProperties StructuralProperties
    \brief Structural descriptors.
    * @{
    */


    namespace structural_properties {

        void planarity(const std::vector<int> atomicNumbers, const std::vector<Vector3d>& positinos, const std::vector<std::vector<int> >& connectivity_matrix, std::vector<Tribool> & _planarity, double threshold = 0.1);
        void planarity(const std::vector<int> atomicNumbers, const std::vector<Vector3d>& positinos, const std::vector<int>& atoms, const std::vector<std::vector<int> >& connectivity_matrix, std::vector<Tribool>& _planarity, double threshold = 0.1);
        void planarity(const std::vector<int> atomicNumbers, const std::vector<Vector3d>& positinos, const std::vector<std::vector<int> >& connectivity_matrix, std::vector<Tribool>& _planarity, std::vector<double>& displacement_esd, double threshold = 0.1);
        void planarity(const std::vector<int> atomicNumbers, const std::vector<Vector3d>& positinos, const std::vector<int>& atoms, const std::vector<std::vector<int> >& connectivity_matrix, std::vector<Tribool>& _planarity, std::vector<double>& displacement_esd, double threshold = 0.1);

        Tribool atomPlanarity(int atomIndex, const std::vector<int> atomicNumbers, const std::vector<Vector3d>& positinos, const std::vector<int>& neighbors, double& displacement_esd, double threshold = 0.1, double valence_threshold=180);
        Tribool planarity(const std::vector<Vector3d>& positinos, double threshold, double& displacement_esd);

        

        Vector3d capping_atom_position(const Vector3d& bonded_atom_r, const Vector3d& directing_atom_r, int bonding_atom_atomic_number);

        /**
        Algorith for finding plane with minimal rmsd of points' distance from the plane as described
        in A.G.URZHUMTSEV 'How to Calculate Planarity Restraints' Acta Cryst. (1991). A47,723-727
        */

        void find_plane(const std::vector<Vector3d>& positions, Vector3d &average, Vector3d &normal);

        /**
        Calculates estimated standard deviation of point distance from plane. 
        Algorith for finding the plane is described in A.G.URZHUMTSEV 'How to Calculate Planarity Restraints' Acta Cryst. (1991). A47,723-727
        */
        double planarity_esd(const std::vector<Vector3d>& positions);
        
        void calculateConnectivity(
            const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity,
            double threshold = 0.4);

        void normalizeXhBondLengths(
            std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers);


        void asymmetricUnitConnectivity(const Crystal &c,std::vector<std::vector<std::pair<int, std::string> > > &connectivity, double threshold);

        void assymetricUnitWithNeighbours(const Crystal &c, 
                                          std::vector< std::pair<int, std::string> > &asuWithNeighbours, 
                                          int neighbourRange,
                                          double threshold);

        void assymetricUnitWithNeighbours(const Crystal &c,
            std::vector< std::pair<int, std::string> > &asuWithNeighbours,
            int neighbourRange,
            double threshold,
            std::vector<int> &shellSizes);

        void assymetricUnitWithNeighbours(const Crystal &c,
            std::vector< std::pair<int, std::string> > &asuWithNeighbours,
            std::vector<int> &atomicNumbers,
            std::vector<Vector3d> &positions,
            std::vector<std::string> &labels,
            int neighbourRange,
            double threshold,
            std::vector<int> &shellSizes);


        void assymetricUnitWithNeighbours(
            const Crystal &crystal, 
            std::vector<int> &atomicNumbers,
            std::vector<Vector3d> &positions,
            std::vector<std::string> &labels,
            int neighbourRange,
            double threshold);

        void assymetricUnitWithNeighbours(
            const Crystal &crystal,
            std::vector<int> &atomicNumbers,
            std::vector<Vector3d> &positions,
            std::vector<std::string> &labels,
            int neighbourRange,
            double threshold,
            std::vector<int> &shellSizes);

        void calcUnitCellConnectivity_Boxes(
            const UnitCellContent& uc,
            std::vector<std::vector<UnitCellContent::AtomID> >& connectivity,
            double treshold);

        void calcUnitCellConnectivity_Simple(
            const UnitCellContent& uc,
            std::vector<std::vector<UnitCellContent::AtomID> >& connectivity,
            double treshold);

        void calcUnitCellConnectivity(
            const UnitCellContent &uc, 
            std::vector<std::vector<UnitCellContent::AtomID> > &connectivity,
            double treshold,
            const std::string &method = "boxes");
        
        void graphToNthNeighbour(const UnitCellContent &unitCellContent, const std::vector<UnitCellContent::AtomID> &startingSet,
            std::vector<UnitCellContent::AtomID> &graph, int n, double threshold);

        void graphToNthNeighbour(const UnitCellContent &unitCellContent, const std::vector<UnitCellContent::AtomID> &startingSet,
            std::vector<UnitCellContent::AtomID> &graph, int n, double threshold, std::vector<int> &shellSizes);


        void calculateConnectivity(
            const std::vector<std::vector<UnitCellContent::AtomID> > &unitCellConnectivity,
            const std::vector<UnitCellContent::AtomID> &atoms,
            std::vector<std::vector<int> > &connectivity);

        void findIncludingMolecule(
             int atom, 
             const std::vector<std::vector<UnitCellContent::AtomID> > &connectivity,
             std::vector<UnitCellContent::AtomID> &molecule,
             std::vector<std::pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > &networkBonds);


        // treat atomsToDisconnect as if they wer not connected to other atoms
        // disconnectedBonds[i].first atom in molecule, .second atom not in molecule
        void findIncludingMolecule(
            const UnitCellContent::AtomID &atomId,
            const std::vector<UnitCellContent::AtomID> &atomsToDisconnect,
            const std::vector<std::vector<UnitCellContent::AtomID> >& connectivity,
            std::vector<UnitCellContent::AtomID>& molecule,
            std::vector<std::pair<UnitCellContent::AtomID, UnitCellContent::AtomID> >& disconnectedBonds,
            const UnitCellContent &uc,
            int maxSize = 100000);


		/**
		Finds molecules containing the atoms with index in UnitCellContent given by the list atom.
		Each atom is contained only once and atom present in one of the molecules does not have to
		have UnitCellContent::AtomID::unitCellPosition = [0, 0, 0]
		*/
		
        void findIncludingMolecule(
            const std::vector<int> &atom,
            const UnitCellContent &uc,
            std::vector<std::vector<UnitCellContent::AtomID> > &molecules,
            std::vector< std::vector< std::pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > &networkBonds,
            double threshold);


        void splitIntoMolecules(
            const std::vector<int> &atomicNumbers,
            const std::vector<Vector3d> &positionsInAngstroms,
            std::vector<std::vector<int> > &molecules,
            double threshold);

		/**
		eachAtomPresentOnlyOnce - 
                |        |
        ________|________|________
                |        |
              A-|-B    A-|-B
        ________|________|________
                |        |
                |        |
        if true only one copy of A--B is included
		*/
        void splitUnitCellIntoMolecules(
                 const UnitCellContent &uc,
                 std::vector<std::vector<UnitCellContent::AtomID> > &molecules,
                 std::vector< std::vector< std::pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > &networkBonds,
                 double threshold = 0.4,
			     bool eachAtomPresentOnlyOnce = true);
        //
        void groupSymmetryRelatedMolecules(
            const UnitCellContent& uc,
            std::vector<std::vector<UnitCellContent::AtomID> >& molecules,
            //[group][molecule in group] .first-index .second-symmetry operation transforming representative [0] molecule to the one with 
            // index moleculeGroups[group][molecule in group]
            std::vector<std::vector<std::pair<int, SpaceGroupOperation> > >& moleculeGroups);

        void splitAsymmetricUnitIntoMolecules(
            const UnitCellContent &uc,
            std::vector<std::vector<UnitCellContent::AtomID> > &molecules,
            std::vector< std::vector< std::pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > &networkBonds,
            double threshold);

        /** if vdwThreshold then threshold for including neighbour molecules is 
        sum of van der Waals radii times the threshold, otherwise the threshold 
        is in Angstroms
        */
		void makeCluster(
			const UnitCellContent& uc, 
			const std::vector<UnitCellContent::AtomID> &centralPart,
			std::vector<UnitCellContent::AtomID> &clusterAtoms,
			double threshold, bool vdwThreshold);

        void makeCluster(
            const UnitCellContent& uc,
            const std::vector<UnitCellContent::AtomID>& centralPart,
            const std::vector<std::vector<UnitCellContent::AtomID> > &unitCellMolecules,
            std::vector<UnitCellContent::AtomID>& clusterAtoms,
            double threshold, bool vdwThreshold);

        

        double interatomicDistance(
            const Crystal &crystal,
            int atom1,
            int atom2,
            const SpaceGroupOperation& s1 = SpaceGroupOperation(),
            const SpaceGroupOperation& s2 = SpaceGroupOperation());
        /**
        vcov should be set to structural_parameters_convention::XyzCoordinateSystem::cartesian
        otherwise it will be slow
        */
        double interatomicDistance(
            const Crystal& crystal,
            int atom1,
            int atom2,
            const CrystalVarianceCovarianceMatrix &vcov,
            double &stadardDeviation,
            const SpaceGroupOperation& s1 = SpaceGroupOperation(),
            const SpaceGroupOperation& s2 = SpaceGroupOperation());

        bool molecularGraphIsomorphism(
            const std::vector<int>& z1, 
            const std::vector<std::vector<int> > &connectivity1,
            const std::vector<int>& z2,
            const std::vector<std::vector<int> >& connectivity2);
    }
    /** @}*/
}

#endif
