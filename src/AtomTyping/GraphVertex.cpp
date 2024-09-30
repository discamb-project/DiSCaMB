#include "discamb/AtomTyping/GraphVertex.h"

#include <algorithm>

namespace discamb {

    GraphVertex::GraphVertex()
    {
    }


    GraphVertex::GraphVertex(
        const AtomDescriptors &descriptors,
        bool allZdefinedNeighborsAsVertices)
    {
        set(descriptors, allZdefinedNeighborsAsVertices);
    }


    GraphVertex::GraphVertex(
        const AtomInStructureDescriptors &descriptors)
    {
        set(descriptors);
    }


    void GraphVertex::set(
        const AtomDescriptors &descriptors, 
        bool allZdefinedNeighborsAsVertices)
    {
        atomicNumberRange = descriptors.atomic_number_range;
        nNeighbors = descriptors.nNeighbours;
        //nNeighbors = descriptors. descriptors.neighborsAtomicNumbers.size();
        additionalNeighboursPossible = !descriptors.fixedNumberOfNeighbors;
        anyAtomicNumber = descriptors.anyAtomicNumber;
        //neighborsAtomicNumbersUniquelyDefined = descriptors.neighborsAtomicNumbersUniquelyDefined;
        neighborsAtomicNumbers = descriptors.neighborsAtomicNumbers;
        neighborsAtomicNumberRanges = descriptors.neighborsAtomicNumberRanges;
        //allNonX_NeighborsNamed = descriptors.allNonX_NeighborsNamed;
        allNonX_NeighborsNamed = allZdefinedNeighborsAsVertices;


        neighborsAtomicNumbersUniquelyDefined = true;
        if (!neighborsAtomicNumberRanges.empty() ||
            additionalNeighboursPossible ||
            nNeighbors != neighborsAtomicNumbers.size())
            neighborsAtomicNumbersUniquelyDefined = false;


        //neighborsAtomicNumbers.clear();
        //neighborsWithUndefinedAtomicNumbers.clear();

        //int neighbor, nNeighbors = descriptors.neighborsAtomicNumbers.size();

        //for (neighbor = 0; neighbor < nNeighbors; neighbor++)
        //    if (descriptors.neighborsAtomicNumbers[neighbor] > 0)
        //        neighborsAtomicNumbers.push_back(descriptors.neighborsAtomicNumbers[neighbor]);
        //    else
        //        neighborsWithUndefinedAtomicNumbers.push_back(descriptors.neighborsAtomicNumbers[neighbor]);

        planar = descriptors.planar;

        ringsSizes.clear();
        int ringIndex, nRings;

        nRings = descriptors.ringInfo.labeledContainingRings.size();

        for (ringIndex = 0; ringIndex < nRings; ringIndex++)
            ringsSizes.push_back(descriptors.ringInfo.labeledContainingRings[ringIndex].first);

        ringsSizes.insert(ringsSizes.end(), descriptors.ringInfo.nonLabeledContainingRings.begin(), descriptors.ringInfo.nonLabeledContainingRings.end());

        std::sort(ringsSizes.begin(), ringsSizes.end());

        possibleMoreRings = descriptors.ringInfo.inAnyAdditionalRing;
        atomInStructure = false;
        isAuxiliaryVertex = false;
        isInRing = descriptors.ringInfo.inRing;
        label = descriptors.label;

        // AtomRingInfo::InRing isIn3Ring, isIn4Ring;
        // int n3Rings, n4Rings;
        isIn3Ring = descriptors.ringInfo.in3Ring;
        isIn4Ring = descriptors.ringInfo.in4Ring;
        n3Rings = descriptors.ringInfo.n3rings;
        n4Rings = descriptors.ringInfo.n4rings;
    }


    void GraphVertex::set(
        const AtomInStructureDescriptors &descriptors)
    {
        label = descriptors.label;
        atomicNumberRange.clear();
        atomicNumberRange .insert(descriptors.atomicNumber);
        nNeighbors = descriptors.nNeighbors;
        additionalNeighboursPossible = false;
        neighborsAtomicNumbers = descriptors.neighborsFormula;
        neighborsAtomicNumbersUniquelyDefined = true;
        planar = descriptors.planar;
        ringsSizes = descriptors.planarRingsSizes;
        possibleMoreRings = false;
        atomInStructure = true;
        isAuxiliaryVertex = false;

        n3Rings = descriptors.n3memberRings;// .threeMemberRingsIndices.size();
        n4Rings = descriptors.n4memberRings;// .fourMemberRingsIndices.size();
        isIn3Ring = (n3Rings == 0 ? Tribool::False : Tribool::True);
        isIn4Ring = (n4Rings == 0 ? Tribool::False : Tribool::True);
    }

}