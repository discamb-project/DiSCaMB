#include "discamb/AtomTyping/RingMatchChecker.h"

#include "discamb/BasicUtilities/on_error.h"

#include <set>
#include <algorithm>


using namespace std;

namespace discamb {

    RingMatchChecker::RingMatchChecker()
    {
        mHasLabeledRings = false;
    }

    RingMatchChecker::~RingMatchChecker()
    {
    }

    void RingMatchChecker::setAtomType(
        const AtomType &atomType)
    {
        mN_AtomsInAtomType = atomType.atoms.size();
        mHasLabeledRings = !atomType.ringLabels.empty();

        if (!mHasLabeledRings)
            return;


        ///////////////////////////////////
        //
        //   prepare for the general case
        //
        ///////////////////////////////////

        int i, nRingSizes, nRings;
        int nAtoms, atomIndex, ringIndex;
        map<string, int> ringLabelToIndexMap;
        string ringLabel;

        nRings = atomType.ringLabels.size();
        mAtomTypeRingLabels = atomType.ringLabels;

        for (i = 0; i < nRings; i++)
            ringLabelToIndexMap[atomType.ringLabels[i]] = i;

        nAtoms = atomType.atoms.size();

        mTypeRingsAtStructureAtoms.resize(nAtoms);
        mAtomInTypeRrequiredRings.resize(nAtoms);
        mAtomInTypeForbidenRings.resize(nAtoms);

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            nRings = atomType.atoms[atomIndex].ringInfo.labeledContainingRings.size();

            for (i = 0; i < nRings; i++)
            {
                ringLabel = atomType.atoms[atomIndex].ringInfo.labeledContainingRings[i].second;
                ringIndex = ringLabelToIndexMap[ringLabel];
                mAtomInTypeRrequiredRings[atomIndex].push_back(ringIndex);
            }

            nRings = atomType.atoms[atomIndex].ringInfo.labeledNonContainingRings.size();

            for (i = 0; i < nRings; i++)
            {
                ringLabel = atomType.atoms[atomIndex].ringInfo.labeledNonContainingRings[i].second;
                ringIndex = ringLabelToIndexMap[ringLabel];
                mAtomInTypeForbidenRings[atomIndex].push_back(ringIndex);
            }

        }

        // divide atom type rings according to size

        std::set<int> uniqueRingSizes(atomType.ringSizes.begin(), atomType.ringSizes.end());
        map<int, int> ringSizeToRingGroupMap;
        mRingSizes.clear();
        mRingSizes.insert(mRingSizes.begin(), uniqueRingSizes.begin(), uniqueRingSizes.end());
        nRingSizes = mRingSizes.size();

        for (i = 0; i < nRingSizes; i++)
            ringSizeToRingGroupMap[mRingSizes[i]] = i;

        mSizeGroupedTypeRingIndices.resize(nRingSizes);

        nRings = atomType.ringSizes.size();
        for (i = 0; i < nRings; i++)
            mSizeGroupedTypeRingIndices[ringSizeToRingGroupMap[atomType.ringSizes[i]]].push_back(i);

    }

    void RingMatchChecker::collectRingsOfMatchingStructureAtoms(
        const StructureWithDescriptors &describedStructure,
        const std::vector<int> &substructureToStructureAtomsMap,
        const std::vector<int> &typeAtoms2SubstructureAtomsMap,
        std::vector<int> &ringIndices)
        const
    {
        ringIndices.clear();
        set<int> uniqueRingIndices;
        int atomIndex, nAtoms = typeAtoms2SubstructureAtomsMap.size();
        int atomInStructureIndex, atomInSubstructureIndex;

        for (atomIndex = 0; atomIndex < mN_AtomsInAtomType; atomIndex++)
        {
            atomInSubstructureIndex = typeAtoms2SubstructureAtomsMap[atomIndex];
            atomInStructureIndex = substructureToStructureAtomsMap[atomInSubstructureIndex];
            uniqueRingIndices.insert(describedStructure.atomDescriptors[atomInStructureIndex].planarRingsIndices.begin(),
                describedStructure.atomDescriptors[atomInStructureIndex].planarRingsIndices.end());
        }
        ringIndices.insert(ringIndices.end(), uniqueRingIndices.begin(), uniqueRingIndices.end());
    }

    void RingMatchChecker::sortInvolvedStructureRingsBySize(
        const StructureWithDescriptors &describedStructure,
        const std::vector<int> &involvedStructureRingIndices,
        const std::vector<int> &involvedRingSizes,
        std::vector<std::vector<int> > &sortedInvolvedStructureRings)
        const
    {
        int i, ringSetIndex, ringSize, nRings = involvedStructureRingIndices.size();
        vector<int>::const_iterator ringSizeIterator;

        sortedInvolvedStructureRings.clear();
        sortedInvolvedStructureRings.resize(involvedRingSizes.size());

        for (i = 0; i < nRings; i++)
        {
            ringSize = describedStructure.planarRings[involvedStructureRingIndices[i]].size();
            ringSizeIterator = find(involvedRingSizes.begin(), involvedRingSizes.end(), ringSize);

            if (ringSizeIterator != involvedRingSizes.end())
            {
                ringSetIndex = std::distance(involvedRingSizes.begin(), ringSizeIterator);
                sortedInvolvedStructureRings[ringSetIndex].push_back( i );
            }
        }
    }


    bool RingMatchChecker::match(
        const StructureWithDescriptors &describedStructure,
        const std::vector<int> &substructureToStructureAtomsMap,
        const std::vector<int> &typeAtoms2SubstructureAtomsMap,
        std::map<std::string, int> &labeledRingsMap)
        const
    {
        labeledRingsMap.clear();

        // nothing to check - no labeled rings (non-labeled rings were checked during atom mapping)
        if (!mHasLabeledRings)
            return true;
                


        /////////////////////////
        //
        //    general check
        //
        /////////////////////////

        //get info on rings of the atoms in structure mapped onto atoms in atom type definition

        vector<int> involvedStructureRings;
        int nInvolvedStructureRings;
        /*
        groups rings by size, involvedRingsGroupedBySize[i][j] is an index of substructure rings which is j-th ring 
        of size mRingSizes[i] among involved rings rings. It corresponds to involvedStructureRings[involvedRingsGroupedBySize[i][j]]
        structure ring
        */
        vector<vector<int> > involvedRingsGroupedBySize;
        map<int, int> structure2sustructureRing;
        collectRingsOfMatchingStructureAtoms(describedStructure, substructureToStructureAtomsMap,
                                             typeAtoms2SubstructureAtomsMap, involvedStructureRings);
        sortInvolvedStructureRingsBySize(describedStructure, involvedStructureRings, mRingSizes, involvedRingsGroupedBySize);

        nInvolvedStructureRings = involvedStructureRings.size();

        mRingMap.resize(nInvolvedStructureRings);

        for (int i = 0; i < nInvolvedStructureRings; i++)
            structure2sustructureRing[involvedStructureRings[i]] = i;

        // make surjections sets

        int nRingGroups, ringGroupIndex, nStructureRings, nAtomTypeRings;

        nRingGroups = mRingSizes.size();
        mRingInjectionIterators.resize(nRingGroups);
        for (ringGroupIndex = 0; ringGroupIndex < nRingGroups; ringGroupIndex++)
        {
            nStructureRings = involvedRingsGroupedBySize[ringGroupIndex].size();
            nAtomTypeRings = mSizeGroupedTypeRingIndices[ringGroupIndex].size();
            if (nStructureRings < nAtomTypeRings)
                return false;

            mRingInjectionIterators[ringGroupIndex].set(nStructureRings, nAtomTypeRings);
        }

        // iterate over all maps of the structure rings onto type rings

        bool foundRingMap = false;
        do
        {
            updateRingMap(typeAtoms2SubstructureAtomsMap, involvedRingsGroupedBySize);
            foundRingMap = checkRingMap(typeAtoms2SubstructureAtomsMap, substructureToStructureAtomsMap,
                                        structure2sustructureRing, describedStructure);
            if (foundRingMap)
            {
                for (int i = 0, n = mRingMap.size(); i < n; i++)
                    if (mRingMap[i].has_value())
                        labeledRingsMap[mAtomTypeRingLabels[*mRingMap[i]]] = involvedStructureRings[i];
            }
        } while (nextRingSuriection() && !foundRingMap);

        return foundRingMap;
    }


    bool RingMatchChecker::nextRingSuriection()
        const
    {
        bool foundNextSuriection = false;
        int ringGroupIndex, nRingGroups;

        ringGroupIndex = 0;
        nRingGroups = mSizeGroupedTypeRingIndices.size();

        do {
            foundNextSuriection = mRingInjectionIterators[ringGroupIndex].next();
            if (!foundNextSuriection)
                mRingInjectionIterators[ringGroupIndex].setToStart();
            ringGroupIndex++;
        } while (!foundNextSuriection && ringGroupIndex < nRingGroups);

        return foundNextSuriection;
    }

    void RingMatchChecker::updateRingMap(
        const std::vector<int> &typeAtoms2SubstructureAtomsMap,
        const std::vector<std::vector<int> > &sortedInvolvedStructureRings)
        const
    {
        int nRingGroups, nRingsInGroup, groupIndex, ringIndex, substructureRingIndex, typeRingIndex, groupRingIndex;
        string typeRingLabel;

        nRingGroups = mSizeGroupedTypeRingIndices.size();
        
        mRingMap.assign(mRingMap.size(), optional<int>());

        for (groupIndex = 0; groupIndex < nRingGroups; groupIndex++)
        {
            nRingsInGroup = mSizeGroupedTypeRingIndices[groupIndex].size();

            for (ringIndex = 0; ringIndex < nRingsInGroup; ringIndex++)
            {
                //substructureRingIndex = mSizeGroupedStructureRings[groupIndex][mRingSurjections[groupIndex](ringIndex)];
                groupRingIndex = mRingInjectionIterators[groupIndex](ringIndex);
                //substructureRingIndex = sortedInvolvedStructureRings[groupIndex][ringIndex];
                substructureRingIndex = sortedInvolvedStructureRings[groupIndex][groupRingIndex];
                typeRingIndex = mSizeGroupedTypeRingIndices[groupIndex][ringIndex];
                mRingMap[substructureRingIndex] = typeRingIndex;
            }
        }

    }

    bool RingMatchChecker::checkRingMap(
        const std::vector<int> &typeAtoms2SubstructureAtomsMap,
        const std::vector<int> &substructureToStructureAtomsMap,
        std::map<int, int> structure2sustructureRing,
        const StructureWithDescriptors &describedStructure)
        const
    {
        // its rings <- type_atom -> structure_atom -> its rings -> map to type_atom rings
        int atomIndex, atomInStructure, atomInSubstructure, nAtoms = mAtomInTypeRrequiredRings.size();
        int ringCounter, nRings, structureRingIndex, substructureRingIndex;
        
        // loop over atoms in type definition
        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            mTypeRingsAtStructureAtoms[atomIndex].clear();
            atomInSubstructure = typeAtoms2SubstructureAtomsMap[atomIndex];
            atomInStructure = substructureToStructureAtomsMap[atomInSubstructure];
            nRings = describedStructure.atomDescriptors[atomInStructure].planarRingsIndices.size();

            for (ringCounter = 0; ringCounter < nRings; ringCounter++)
            {
                structureRingIndex = describedStructure.atomDescriptors[atomInStructure].planarRingsIndices[ringCounter];
                substructureRingIndex = structure2sustructureRing[structureRingIndex];
                if (mRingMap[substructureRingIndex].has_value())
                    mTypeRingsAtStructureAtoms[atomIndex].push_back(mRingMap[substructureRingIndex].value());
            }
        }

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            nRings = mAtomInTypeRrequiredRings[atomIndex].size();
            for (ringCounter = 0; ringCounter < nRings; ringCounter++)
                if (find(mTypeRingsAtStructureAtoms[atomIndex].begin(), mTypeRingsAtStructureAtoms[atomIndex].end(), mAtomInTypeRrequiredRings[atomIndex][ringCounter])
                    == mTypeRingsAtStructureAtoms[atomIndex].end())
                    return false;

            nRings = mAtomInTypeForbidenRings[atomIndex].size();
            for (ringCounter = 0; ringCounter < nRings; ringCounter++)
                if (find(mTypeRingsAtStructureAtoms[atomIndex].begin(), mTypeRingsAtStructureAtoms[atomIndex].end(), mAtomInTypeForbidenRings[atomIndex][ringCounter])
                    != mTypeRingsAtStructureAtoms[atomIndex].end())
                    return false;
        }

        return true;
    }

}
