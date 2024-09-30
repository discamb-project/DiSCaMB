#include "discamb/AtomTyping/RingMatchCheckerTypeType.h"

#include "discamb/BasicUtilities/OnError.h"
#include "discamb/MathUtilities/set_theory.h"

#include <set>
#include <algorithm>


using namespace std;

namespace discamb {

    RingMatchCheckerTypeType::RingMatchCheckerTypeType()
    {
        mHasLabeledRings = false;
    }

    RingMatchCheckerTypeType::~RingMatchCheckerTypeType()
    {
    }

    void RingMatchCheckerTypeType::setAtomType(
        const AtomType &atomType)
    {
        //mAtomType = atomType;
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


        mSizeGroupedAtomTypeRingLabels.clear();

        for (i = 0; i < nRings; i++)
        {
            mSizeGroupedAtomTypeRingLabels[atomType.ringSizes[i]].push_back(atomType.ringLabels[i]);
            ringLabelToIndexMap[atomType.ringLabels[i]] = i;
        }

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

    void RingMatchCheckerTypeType::collectRingsOfMatchingSubtypeAtoms(
        const AtomType& subtype,
        const std::vector<int>& typeAtoms2SubtypeAtoms,
        std::map<int, std::set<std::string> > & labeledSubgraphRings)
        const
    {
        labeledSubgraphRings.clear();

        for (int subtypeAtomIdx : typeAtoms2SubtypeAtoms)
        {
            for (auto const& rings : subtype.atoms[subtypeAtomIdx].ringInfo.labeledContainingRings)
                labeledSubgraphRings[rings.first].insert(rings.second);

            for (auto const& rings : subtype.atoms[subtypeAtomIdx].ringInfo.labeledNonContainingRings)
                labeledSubgraphRings[rings.first].insert(rings.second);
        }
    }

    void RingMatchCheckerTypeType::sortInvolvedStructureRingsBySize(
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


    bool RingMatchCheckerTypeType::match(
        const AtomType& subtype,
        const std::vector<int>& typeAtoms2SubtypeAtomsMap)
        const
    {
        
        // nothing to check - no labeled rings (non-labeled rings were checked during atom mapping)
        if (!mHasLabeledRings)
            return true;
                
        // checks if number and size of rings is OK
        // i.e. if for each atom A in type, with reqired labelled rings of sizes 5,6,6 there are 
        // also ring of this size in the corresponging subtype atom, check the same for forbidden rings

        int nAtoms = mAtomInTypeRrequiredRings.size();
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            vector<int> typeRingSizes, subtypeRingSizes;

            int subtypeAtomIdx;

            // check required rings
            int nRequiredRingsAtomInSubtype = subtype.atoms[subtypeAtomIdx].ringInfo.labeledContainingRings.size();
            int nRequiredRingsAtomInType = mAtomInTypeRrequiredRings[atomIdx].size();

            if (nRequiredRingsAtomInType > nRequiredRingsAtomInSubtype)
                return false;

            for (auto ringIdx : mAtomInTypeRrequiredRings[atomIdx])
                typeRingSizes.push_back(mRingSizes[ringIdx]);

            for (auto ring : subtype.atoms[subtypeAtomIdx].ringInfo.labeledContainingRings)
                subtypeRingSizes.push_back(ring.first);

            if (!set_theory::is_subset(subtypeRingSizes, typeRingSizes))
                return false;

            // check forbidden rings

            typeRingSizes.clear();
            subtypeRingSizes.clear();

            int nForbiddenRingsAtomInSubtype = subtype.atoms[subtypeAtomIdx].ringInfo.labeledNonContainingRings.size();
            int nForbiddenRingsAtomInType = mAtomInTypeForbidenRings[atomIdx].size();

            if (nForbiddenRingsAtomInType > nForbiddenRingsAtomInSubtype)
                return false;

            for (auto ringIdx : mAtomInTypeForbidenRings[atomIdx])
                typeRingSizes.push_back(mRingSizes[ringIdx]);

            for (auto ring : subtype.atoms[subtypeAtomIdx].ringInfo.labeledNonContainingRings)
                subtypeRingSizes.push_back(ring.first);

            if (!set_theory::is_subset(subtypeRingSizes, typeRingSizes))
                return false;

        }

        
            

        //

        map<int, set<std::string> > subtypeSubgraphRings;
        map<int, vector<std::string> > subtypeSubgraphRingsVec;
        map<string, string> subtypeTypeRingLabelMap;
        
        collectRingsOfMatchingSubtypeAtoms(subtype, typeAtoms2SubtypeAtomsMap, subtypeSubgraphRings); 
        for (auto& ringsGroup : subtypeSubgraphRings)
            subtypeSubgraphRingsVec[ringsGroup.first].insert(subtypeSubgraphRingsVec[ringsGroup.first].end(), ringsGroup.second.begin(), ringsGroup.second.end());
            
         
        // make injections sets

        int nRingGroups, ringGroupIndex, nSubtypeRings, nAtomTypeRings;

        nRingGroups = mRingSizes.size();
        mRingInjectionIterators.resize(nRingGroups);
        
        for (ringGroupIndex = 0; ringGroupIndex < nRingGroups; ringGroupIndex++)
        {
            nSubtypeRings = 0;
            int ringSize = mRingSizes[ringGroupIndex];
            if (subtypeSubgraphRings.count(ringSize) > 0)
                nSubtypeRings = subtypeSubgraphRings[ringSize].size();
            nAtomTypeRings = mSizeGroupedTypeRingIndices[ringGroupIndex].size();
            if (nSubtypeRings < nAtomTypeRings)
                return false;

            mRingInjectionIterators[ringGroupIndex].set(nSubtypeRings, nAtomTypeRings);
        }

        

        // iterate over all maps of the structure rings onto type rings
        
        bool foundRingMap = false;
        do
        {
            auto sizeGroupedAtomTypeRingLabels = mSizeGroupedAtomTypeRingLabels;
            updateRingMap(sizeGroupedAtomTypeRingLabels, subtypeSubgraphRingsVec, subtypeTypeRingLabelMap);

            foundRingMap = checkRingMap(subtype, subtypeTypeRingLabelMap, typeAtoms2SubtypeAtomsMap);


        } while (nextRingInjection() && !foundRingMap);

        return foundRingMap;



    }


    bool RingMatchCheckerTypeType::nextRingInjection()
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


    void RingMatchCheckerTypeType::updateRingMap(
        std::map<int, std::vector<std::string> > &subtypeSubgraphRingsVec,
        std::map<int, std::vector<std::string> > &typeRingsVec,
        std::map<std::string, std::string> &subtypeTypeRingLabelMap) const
    {
        subtypeTypeRingLabelMap.clear();
        for (const auto& ringGroup : subtypeSubgraphRingsVec)
            for (const string& label : ringGroup.second)
                subtypeTypeRingLabelMap[label] = string();



        int nRingGroups, groupIndex, ringIndex, substructureRingIndex, typeRingIndex, groupRingIndex;
        string typeRingLabel;

        nRingGroups = mSizeGroupedTypeRingIndices.size();

        mRingMap.assign(mRingMap.size(), optional<int>());

        for (groupIndex = 0; groupIndex < nRingGroups; groupIndex++)
        {
            int ringSize = mRingSizes[groupIndex];
            int nRings = typeRingsVec[ringSize].size();

            for (int i = 0; i < nRings; i++)
            {
                string ringType = typeRingsVec[groupIndex][i];
                string ringSubtype = subtypeSubgraphRingsVec[ringSize][mRingInjectionIterators[groupIndex](i)];
                subtypeTypeRingLabelMap[ringSubtype] = ringType;
            }

        }

    }


    //void RingMatchCheckerTypeType::updateRingMap(
    //    const std::vector<int> &typeAtoms2SubstructureAtomsMap,
    //    const std::vector<std::vector<int> > &sortedInvolvedStructureRings)
    //    const
    //{
    //    int nRingGroups, nRingsInGroup, groupIndex, ringIndex, substructureRingIndex, typeRingIndex, groupRingIndex;
    //    string typeRingLabel;

    //    nRingGroups = mSizeGroupedTypeRingIndices.size();
    //    
    //    mRingMap.assign(mRingMap.size(), optional<int>());

    //    for (groupIndex = 0; groupIndex < nRingGroups; groupIndex++)
    //    {
    //        nRingsInGroup = mSizeGroupedTypeRingIndices[groupIndex].size();

    //        for (ringIndex = 0; ringIndex < nRingsInGroup; ringIndex++)
    //        {
    //            //substructureRingIndex = mSizeGroupedStructureRings[groupIndex][mRingSurjections[groupIndex](ringIndex)];
    //            groupRingIndex = mRingInjectionIterators[groupIndex](ringIndex);
    //            //substructureRingIndex = sortedInvolvedStructureRings[groupIndex][ringIndex];
    //            substructureRingIndex = sortedInvolvedStructureRings[groupIndex][groupRingIndex];
    //            typeRingIndex = mSizeGroupedTypeRingIndices[groupIndex][ringIndex];
    //            mRingMap[substructureRingIndex] = typeRingIndex;
    //        }
    //    }

    //}

    bool RingMatchCheckerTypeType::checkRingMap(
        const AtomType& subtype,
        const map<string, string>& subtypeTypeRingLabelMap,
        const std::vector<int>& typeAtoms2SubtypeAtomsMap)
        const
    {
        //std::vector<std::vector<int> > mAtomInTypeRrequiredRings;
        //std::vector<std::vector<int> > mAtomInTypeForbidenRings;
        int nAtoms = mAtomInTypeRrequiredRings.size();
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            //required rings
            int subtypeAtomIdx = typeAtoms2SubtypeAtomsMap[atomIdx];
            vector<string> atomInTypeRingLabels;
            for (int idx : mAtomInTypeRrequiredRings[atomIdx])
                atomInTypeRingLabels.push_back(mAtomTypeRingLabels[idx]);
            vector<string> atomInSubtypeRingLabels;
            for (const auto& item : subtype.atoms[subtypeAtomIdx].ringInfo.labeledContainingRings)
                atomInSubtypeRingLabels.push_back(subtypeTypeRingLabelMap.find(item.second)->second);
            if (!set_theory::is_subset(atomInSubtypeRingLabels, atomInTypeRingLabels))
                return false;
            //forbidden rings
            atomInTypeRingLabels.clear();
            atomInSubtypeRingLabels.clear();

            for (int idx : mAtomInTypeForbidenRings[atomIdx])
                atomInTypeRingLabels.push_back(mAtomTypeRingLabels[idx]);
            
            for (const auto& item : subtype.atoms[subtypeAtomIdx].ringInfo.labeledNonContainingRings)
                atomInSubtypeRingLabels.push_back(subtypeTypeRingLabelMap.find(item.second)->second);
            if (!set_theory::is_subset(atomInSubtypeRingLabels, atomInTypeRingLabels))
                return false;
        }
        return true;
    }

    //bool RingMatchCheckerTypeType::checkRingMap(
    //    const std::vector<int> &typeAtoms2SubstructureAtomsMap,
    //    const std::vector<int> &substructureToStructureAtomsMap,
    //    std::map<int, int> structure2sustructureRing,
    //    const StructureWithDescriptors &describedStructure)
    //    const
    //{
    //    // its rings <- type_atom -> structure_atom -> its rings -> map to type_atom rings
    //    int atomIndex, atomInStructure, atomInSubstructure, nAtoms = mAtomInTypeRrequiredRings.size();
    //    int ringCounter, nRings, structureRingIndex, substructureRingIndex;
    //    
    //    // loop over atoms in type definition
    //    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    //    {
    //        mTypeRingsAtStructureAtoms[atomIndex].clear();
    //        atomInSubstructure = typeAtoms2SubstructureAtomsMap[atomIndex];
    //        atomInStructure = substructureToStructureAtomsMap[atomInSubstructure];
    //        nRings = describedStructure.atomDescriptors[atomInStructure].planarRingsIndices.size();

    //        for (ringCounter = 0; ringCounter < nRings; ringCounter++)
    //        {
    //            structureRingIndex = describedStructure.atomDescriptors[atomInStructure].planarRingsIndices[ringCounter];
    //            substructureRingIndex = structure2sustructureRing[structureRingIndex];
    //            if (mRingMap[substructureRingIndex].has_value())
    //                mTypeRingsAtStructureAtoms[atomIndex].push_back(mRingMap[substructureRingIndex].value());
    //        }
    //    }

    //    for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
    //    {
    //        nRings = mAtomInTypeRrequiredRings[atomIndex].size();
    //        for (ringCounter = 0; ringCounter < nRings; ringCounter++)
    //            if (find(mTypeRingsAtStructureAtoms[atomIndex].begin(), mTypeRingsAtStructureAtoms[atomIndex].end(), mAtomInTypeRrequiredRings[atomIndex][ringCounter])
    //                == mTypeRingsAtStructureAtoms[atomIndex].end())
    //                return false;

    //        nRings = mAtomInTypeForbidenRings[atomIndex].size();
    //        for (ringCounter = 0; ringCounter < nRings; ringCounter++)
    //            if (find(mTypeRingsAtStructureAtoms[atomIndex].begin(), mTypeRingsAtStructureAtoms[atomIndex].end(), mAtomInTypeForbidenRings[atomIndex][ringCounter])
    //                != mTypeRingsAtStructureAtoms[atomIndex].end())
    //                return false;
    //    }

    //    return true;
    //}

}

