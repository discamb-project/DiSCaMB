#pragma once

#include "discamb/StructuralProperties/ConnectivityAlgorithm.h"

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    template<typename BOND_DETECTOR>
    class GenericConnectivityAlgorithm : public ConnectivityAlgorithm
    {
    public:
        GenericConnectivityAlgorithm();
        virtual ~GenericConnectivityAlgorithm();
        virtual void set(const std::string &bondDetectorSettings);

        virtual void calculateConnectivity(
            const std::vector<Vector3d> &positions,
            const std::vector<int> &atomicNumbers,
            std::vector<std::vector<int> > &connectivity) const;


        virtual void calculateConnectivity(
                        const std::vector<Vector3d> &positions,
                        const std::vector<int> &atomicNumbers,
                        const MolecularDisorder &molecularDisorder, 
                        std::vector<std::vector<int> > &connectivity) const;
        
        virtual void calculateConnectivity(
                         const std::vector<Vector3d> &positions,
                         const std::vector<int> &atomicNumbers,
                         const MolecularDisorder &molecularDisorder, 
                         const std::vector<std::vector<int> > &disjonedAtomGroups, 
                         std::vector<std::vector<int> > &connectivity) const;

    private:
        BOND_DETECTOR mBondDetector;
        void atomGroupConnectivity(
                 const std::vector<int> &atomGroup, 
                 const std::vector<Vector3d> &positions,
                 const std::vector<int> &atomicNumbers,
                 const MolecularDisorder &molecularDisorder, 
                 std::vector<std::vector<int> > &connectivity) const;
    };

    // -------------------- IMPLEMENTATION  ---------------------------------------


    template<typename BOND_DETECTOR>
    GenericConnectivityAlgorithm<BOND_DETECTOR>::GenericConnectivityAlgorithm()
    {
    }


    template<typename BOND_DETECTOR>
    GenericConnectivityAlgorithm<BOND_DETECTOR>::~GenericConnectivityAlgorithm()
    {
    }


    template<typename BOND_DETECTOR>
    void GenericConnectivityAlgorithm<BOND_DETECTOR>::set(
        const std::string &bondDetectorSettings)
    {
        mBondDetector.set(bondDetectorSettings);
    }

    template<typename BOND_DETECTOR>
    void GenericConnectivityAlgorithm<BOND_DETECTOR>::atomGroupConnectivity(
        const std::vector<int> &atomGroup,
        const std::vector<Vector3d> &positions,
        const std::vector<int> &atomicNumbers,
        const MolecularDisorder &molecularDisorder,
        std::vector<std::vector<int> > &connectivity)
        const
    {
        int i, j, nAtoms, atom1_index, atom2_index;
        int disorderAssembly1, disorderAssembly2, disorderGroup1, disorderGroup2;
        int atomic_number_1, atomic_number_2;
        Vector3d atom1_position, atom2_position;

        nAtoms = atomGroup.size();

        for (i = 0; i < nAtoms; i++)
        {
            atom1_index = atomGroup[i];
            connectivity[atom1_index].clear();
            atom1_position = positions[atom1_index];
            atomic_number_1 = atomicNumbers[atom1_index];
            disorderAssembly1 = molecularDisorder.atomDisorderAssembly(atom1_index);
            disorderGroup1 = molecularDisorder.atomDisorderGroup(atom1_index);

            for (j = 0; j < i; j++)
            {
                atom2_index = atomGroup[j];

                disorderAssembly2 = molecularDisorder.atomDisorderAssembly(atom2_index);

                if (disorderAssembly1 >= 0 && disorderAssembly1 == disorderAssembly2) // the same disorder assembly
                {
                    disorderGroup2 = molecularDisorder.atomDisorderGroup(atom2_index);

                    // different disorder group within the same disorder assembly = no bond since they can not be present at the same time
                    if (disorderGroup2 != disorderGroup1)
                        continue;
                }

                atomic_number_2 = atomicNumbers[atom2_index];

                if (mBondDetector.areBonded(atomic_number_1, atom1_position, atomic_number_2, positions[atom2_index]))
                {
                    connectivity[atom1_index].push_back(atom2_index);
                    connectivity[atom2_index].push_back(atom1_index);
                }

            }
        }

    }

    template<typename BOND_DETECTOR>
    void GenericConnectivityAlgorithm<BOND_DETECTOR>::calculateConnectivity(
        const std::vector<Vector3d> &positions,
        const std::vector<int> &atomicNumbers,
        std::vector<std::vector<int> > &connectivity) const
    {
        BOND_DETECTOR bondDetector;
        int nAtoms, i, j;
        

        bondDetector.setThreshold(0.4);//angstroms

        nAtoms = atomicNumbers.size();
        connectivity.clear();
        connectivity.resize(nAtoms);

        for (i = 0; i < nAtoms; i++)
            for (j = 0; j < i; j++)
                if (bondDetector.areBonded(atomicNumbers[i], positions[i], atomicNumbers[j], positions[j]))
                {
                    connectivity[i].push_back(j);
                    connectivity[j].push_back(i);
                }

    }


    template<typename BOND_DETECTOR>
    void GenericConnectivityAlgorithm<BOND_DETECTOR>::calculateConnectivity(
        const std::vector<Vector3d> &positions,
        const std::vector<int> &atomicNumbers,
        const MolecularDisorder &molecularDisorder,
        const std::vector<std::vector<int> > &disjonedAtomGroups,
        std::vector<std::vector<int> > &connectivity)
        const
    {
        int i, nAtoms, nDisjoined, disjoinedIndex, nAtomsInDisjoinedGroup;
        std::vector<int> notDisjoinedIndices;
        std::vector<bool> notDisjoined;

        if (disjonedAtomGroups.empty())
        {
            calculateConnectivity(positions, atomicNumbers, molecularDisorder, connectivity);
            return;
        }

        nAtoms = positions.size();
        connectivity.resize(nAtoms);
        notDisjoined.resize(nAtoms, true);

        nDisjoined = disjonedAtomGroups.size();

        for (disjoinedIndex = 0; disjoinedIndex < nDisjoined; disjoinedIndex++)
        {
            atomGroupConnectivity(disjonedAtomGroups[disjoinedIndex], positions, atomicNumbers, molecularDisorder, connectivity);
            nAtomsInDisjoinedGroup = disjonedAtomGroups[disjoinedIndex].size();
            for (i = 0; i < nAtoms; i++)
                notDisjoined[disjonedAtomGroups[disjoinedIndex][i]] = false;
        }

        notDisjoinedIndices.reserve(nAtoms);

        for (i = 0; i < nAtoms; i++)
            if (notDisjoined[i])
                notDisjoinedIndices.push_back(i);

        atomGroupConnectivity(notDisjoinedIndices, positions, atomicNumbers, molecularDisorder, connectivity);

    }

    template<typename BOND_DETECTOR>
    void GenericConnectivityAlgorithm<BOND_DETECTOR>::calculateConnectivity(
        const std::vector<Vector3d> &positions,
        const std::vector<int> &atomicNumbers,
        const MolecularDisorder &molecularDisorder,
        std::vector<std::vector<int> > &connectivity)
        const
    {
        int i, nAtoms = positions.size();
        std::vector<int> atomIndices(nAtoms);
        connectivity.resize(nAtoms);
        for (i = 0; i < nAtoms; i++)
            atomIndices[i] = i;
        atomGroupConnectivity(atomIndices, positions, atomicNumbers, molecularDisorder, connectivity);
    }

    /** @}*/
}



