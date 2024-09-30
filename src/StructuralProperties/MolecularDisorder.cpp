#include "discamb/StructuralProperties/MolecularDisorder.h"
#include "discamb/BasicUtilities/OnError.h"

using namespace std;
namespace discamb {
    MolecularDisorder::MolecularDisorder() : mN_Atoms(0)
    {
    }

    MolecularDisorder::~MolecularDisorder()
    {
    }

    void MolecularDisorder::setMoleculeSize(int nAtoms)
    {
        mN_Atoms = nAtoms;
        clear();
        mAtomToGroupMap.assign(mN_Atoms, -1);
    }

    int MolecularDisorder::moleculeSize()
        const
    {
        return mN_Atoms;
    }

    void MolecularDisorder::addDisorderGroup(
        const vector<int> &disorderedAtomsGroup)
    {
        int i, nAtoms, disorderGroupIndex;
        nAtoms = disorderedAtomsGroup.size();

        // check if atom indices are correct (in the range 0 - mN_Atoms-1 and do not refere to atoms already assigned to some disorder group)

        for (i = 0; i < nAtoms; i++)
        {
            if (mAtomToGroupMap.at(disorderedAtomsGroup[i]) != -1)
                on_error::throwException("an attempt to assign atom to more than one disorder group", __FILE__, __LINE__);
        }

        // update info on disorder

        disorderGroupIndex = mDisorderGroups.size();
        mDisorderGroups.push_back(disorderedAtomsGroup);
        for (i = 0; i < nAtoms; i++)
            mAtomToGroupMap[disorderedAtomsGroup[i]] = disorderGroupIndex;
        mGroupToAssmblyMap.push_back(-1); // the newly added group is no assigned to any disorder assembly (as indicated by -1)
    }

    int MolecularDisorder::nDisorderAtomGroups()
        const
    {
        return mDisorderGroups.size();
    }

    void MolecularDisorder::getDisorderAtomGroup(
        int groupIndex,
        vector<int> &atomIndices)
        const
    {
        atomIndices = mDisorderGroups.at(groupIndex);
    }



    int MolecularDisorder::groupAssembly(
        int disorderAtomGroup)
        const
    {
        return mGroupToAssmblyMap[disorderAtomGroup];
    }

    void MolecularDisorder::addDisorderAssembly(
        const std::vector<int> &newAssembly)
    {
        int i, newAssemblySize, newAssemblyIndex;

        newAssemblySize = newAssembly.size();

        // check input correctness

        for (i = 0; i < newAssemblySize; i++)
        {
            if (mGroupToAssmblyMap.at(newAssembly[i]) != -1)
                on_error::throwException("an attempt to assign disorder group to more than one disorder assembly", __FILE__, __LINE__);
        }

        // set new info

        newAssemblyIndex = mDisorderAssemblies.size();
        mDisorderAssemblies.push_back(newAssembly);
        for (i = 0; i < newAssemblySize; i++)
            mGroupToAssmblyMap[newAssembly[i]] = newAssemblyIndex;
    }

    void MolecularDisorder::getDisorderAssembly(
        int index,
        std::vector<int> &disorderedAtomsGroups)
        const
    {
        disorderedAtomsGroups = mDisorderAssemblies.at(index);
    }

    int MolecularDisorder::nDisorderAssemblies()
        const
    {
        return mDisorderAssemblies.size();
    }

    void MolecularDisorder::removeAssembly(
        int index)
    {
        mDisorderAssemblies.erase(mDisorderAssemblies.begin() + index);
        for (int i = 0; i < mGroupToAssmblyMap.size(); i++)
            if (mGroupToAssmblyMap[i] == index)
                mGroupToAssmblyMap[i] = -1;
    }


    int MolecularDisorder::atomDisorderGroup(
        int atomIndex)
        const
    {
        return mAtomToGroupMap.at(atomIndex);
    }

    int MolecularDisorder::atomDisorderAssembly(
        int atomIndex)
        const
    {
        int atomGroup = mAtomToGroupMap.at(atomIndex);
        if (atomGroup < 0)
            return -1;
        else
            return mGroupToAssmblyMap[mAtomToGroupMap.at(atomIndex)];
    }

    void MolecularDisorder::clear()
    {
        mDisorderGroups.clear();
        mDisorderAssemblies.clear();
        mGroupToAssmblyMap.clear();
        mAtomToGroupMap.clear();
    }

} // namespcae discamb
