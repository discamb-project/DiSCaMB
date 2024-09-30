#pragma once

#include <vector>
#include <cstddef>

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    /**

    Describes disorder in set of atoms taken out of average structure of disordered crystal.
    Notation is similar to the one used in CIF format. A group of disordered atoms which are
    always simulateneously present in crystal structure is called 'disorder group' while a group
    of alternative disorder groups is called 'disorder assembly' - as ilustrated on the picture:
    \image html disorder_group_assembly.png

    */

    class MolecularDisorder
    {
    public:
        MolecularDisorder();
        ~MolecularDisorder();

        /**
        Defines number of atoms in the described molecule (causes also clearing all previous info on disorder).
        */
        void setMoleculeSize(int nAtoms);

        /**
        Returns number of atoms in the described molecule.
        */
        int moleculeSize() const;

        // disord group

        /**
        Defines new group of disordered atoms which always appears together, no atoms can belong to more than one disorder group
            \param disorderedAtomsGroup - list of atoms belonging tothe new disorder group (as indices between 0 and moleculeSize()-1)
        */

        void addDisorderGroup(const std::vector<int> &disorderedAtomsGroup);

        /**
        Returns number of groups of disordered atoms.
        */
        int nDisorderAtomGroups() const;

        /**
        Provides list of atoms belonging to given group of disordered atoms.
        \param groupIndex - index of the group
        \param atomIndices - indices of the atoms belonging to the group
        */
        void getDisorderAtomGroup(int groupIndex, std::vector<int> &atomIndices) const;

        /**
        Return index of disorder assembly containing given disorder group (or -1 if none contains it).
        \param disorderAtomGroup - index of the group of disordered atoms
        */
        int groupAssembly(int disorderAtomGroup) const;

        // alternatives

        /**
        Defines new disorder assembly (set of disorder groups which are alternative to each other - i.e. atoms belonging to one of the groups can not appear at
        same time as atoms belonging to the another alternative disorder group).
        \param disorderedAtomsGroups - list of indices of disorder grops belonging to the disorder assembly
        */
        void addDisorderAssembly(const std::vector<int> &disorderedAtomsGroups);

        /**
        Provides list of disorder groups belonging to given disorder assembly.
        \param index - index of the disorder assembly
        \param disorderedAtomsGroups - indices of the disorder atom grops belonging to the assembly
        */
        void getDisorderAssembly(int index, std::vector<int> &disorderedAtomsGroups) const;

        /**
        Returns number of disorder assemblies.
        */
        int nDisorderAssemblies() const;

        /**
        Removes given disorder assembly.
        */
        void removeAssembly(int index);

        // atoms

        /**
        Returns index of the diosrder group containing given atom (and -1 if the atom does not belong to any).
        */
        int atomDisorderGroup(int atomIndex) const;

        /**
        Returns index of the diosrder assembly containing the disorder grouop containing given atom (and -1 if thereis no such a disorder assembly)
        */
        int atomDisorderAssembly(int atomIndex) const;

        void clear();

    private:

        int mN_Atoms;
        std::vector<std::vector<int> > mDisorderGroups;
        std::vector<std::vector<int> > mDisorderAssemblies;
        std::vector<int> mGroupToAssmblyMap;
        std::vector<int> mAtomToGroupMap;


    };
    /**@}*/
} // namespace discamb




