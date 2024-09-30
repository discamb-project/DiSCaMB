#include "discamb/AtomTyping/GraphVertexComparator.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/MathUtilities/set_theory.h"

#include <algorithm>
#include <vector>

using namespace std;

namespace discamb {

    namespace {
        bool extractSubformula(
            vector<int>& formula, 
            const vector<int>& subformula)
        {
            for (int z : subformula)
                if (find(formula.begin(), formula.end(), z) != formula.end())
                    formula.erase(find(formula.begin(), formula.end(), z));
                else
                    return false;
            return true;
        }

        // map range index into formula index
        bool zRangesFitFormula(
            const vector<int>& formula,
            const vector<set<int> >& ranges,
            //map ranges to formula
            vector<int>& requiredMap)
        {
            
            if (requiredMap.size() == ranges.size())
                return true;
            // true if formula index can be used for new match
            vector<bool> free(formula.size(),true);
            for (int i = 0; i < requiredMap.size(); i++)
                free[requiredMap[i]] = false;

            int currentRangeIdx = requiredMap.size();
            vector<int> requiredMapCandidate;
            for(int i=0;i<formula.size();i++)
                if (free[i])
                    if (ranges[currentRangeIdx].find(formula[i]) != ranges[currentRangeIdx].end())
                    {
                        requiredMapCandidate = requiredMap;
                        requiredMapCandidate.push_back(i);
                        if (zRangesFitFormula(formula, ranges, requiredMapCandidate))
                        {
                            requiredMap = requiredMapCandidate;
                            return true;
                        }
                    }
            return false;

        }

        bool isSubrange(const set<int>& range, const set<int>& subrange)
        {
            
            for (int z : subrange)
                if (find(range.begin(), range.end(), z) == range.end())
                    return false;
            return true;
        }

        bool isNeighbourSubset(
            const multiset<int> &neighborsAtomicNumbersType,
            const vector<set<int> > &neighborsAtomicNumberRangesType,
            int nAtomsType,
            const multiset<int> &neighborsAtomicNumbersSubtype,
            const vector<set<int> > &neighborsAtomicNumberRangesSubtype,
            int nAtomsSubtype,
            bool additionalSubsetNeighboursPossible)
        {
            int nDefType = neighborsAtomicNumbersType.size() + neighborsAtomicNumberRangesType.size();
            int nDefSubtype = neighborsAtomicNumbersSubtype.size() + neighborsAtomicNumberRangesSubtype.size();
            
            if (!additionalSubsetNeighboursPossible)
                if (nAtomsType != nAtomsSubtype)
                    return false;
            vector<set<int> > neighboursType;
            vector<set<int> > neighboursSubtype;
            set<int> set120;

            for (int i = 0; i < 120; i++)
                set120.insert(i);

            neighboursType = neighborsAtomicNumberRangesType;
            for (int z : neighborsAtomicNumbersType)
                neighboursType.push_back(set<int>({ z }));
            for (int i = nDefType; i < nAtomsType; i++)
                neighboursType.push_back(set120);

            neighboursSubtype = neighborsAtomicNumberRangesSubtype;
            for (int z : neighborsAtomicNumbersSubtype)
                neighboursSubtype.push_back(set<int>({ z }));
            for (int i = nDefSubtype; i < nAtomsSubtype; i++)
                neighboursSubtype.push_back(set120);

            vector<vector<bool> > compatible(nAtomsType, vector<bool>(nAtomsSubtype, true));

            for (int tIdx = 0; tIdx < nAtomsType; tIdx++)
            {
                bool hasCompatible = false;
                for (int sIdx = 0; sIdx < nAtomsSubtype; sIdx++)
                {
                    compatible[tIdx][sIdx] = isSubrange(neighboursType[tIdx], neighboursSubtype[sIdx]);
                    if (compatible[tIdx][sIdx])
                        hasCompatible = true;
                }
                if (!hasCompatible)
                    return false;
            }

            vector<int> idx(nAtomsSubtype);
            for (int i = 0; i < nAtomsSubtype; i++)
                idx[i] = i;

            bool found;
            do
            {
                found = true;
                for (int i = 0; i < nAtomsType; i++)
                    if (!compatible[i][idx[i]])
                        found = false;

            } while (next_permutation(idx.begin(), idx.end()) && !found);
            return found;
        }

        //bool isSubset(const vector<int>& set, const vector<int>& subset)
        //{
        //    std::multiset<int> s(set.begin(), set.end());
        //    
        //    for (int element : subset)
        //    {
        //        auto matchingElement = s.find(element);
        //        if (matchingElement != s.end())
        //            s.erase(matchingElement);
        //        else
        //            return false;
        //    }
        //    return true;
        //}
    }

    GraphVertexComparator::GraphVertexComparator()
    {
    }

    GraphVertexComparator::~GraphVertexComparator()
    {
    }

    bool GraphVertexComparator::compatible(
        void *vertex1,
        void *vertex2)
    {
        return compatible(*static_cast<GraphVertex*>(vertex1), *static_cast<GraphVertex*>(vertex2));
    }

    bool GraphVertexComparator::compatible(
        const GraphVertex &v1,
        const GraphVertex &v2)
    {
        if (!v1.atomInStructure && !v2.atomInStructure)
            return compatibleTypeType(v1, v2);

        const GraphVertex *vertex_type;
        const GraphVertex *vertex_structure;


        if (v1.atomInStructure)
        {
            vertex_structure = &v1;
            vertex_type = &v2;
        }
        else
        {
            vertex_structure = &v2;
            vertex_type = &v1;
        }


        if (v1.isAuxiliaryVertex && v2.isAuxiliaryVertex)
            return true;

        if (v1.isAuxiliaryVertex != v2.isAuxiliaryVertex)
            return false;
        if(!vertex_type->anyAtomicNumber)
            if(vertex_type->atomicNumberRange.find(*vertex_structure->atomicNumberRange.begin())== vertex_type->atomicNumberRange.end())
                return false;
        //if ((v1.atomicNumber != 0) && (v2.atomicNumber != 0))
          //  if (v1.atomicNumber != v2.atomicNumber)
            //    return false;

        if (!neighborsCompatible(vertex_type, vertex_structure))
            return false;

        if (!comparePlanarity(vertex_type->planar, vertex_structure->planar))
            return false;

        if (!compareN_Rings(vertex_type->isIn3Ring, vertex_type->n3Rings, vertex_structure->n3Rings.value()))
            return false;

        if (!compareN_Rings(vertex_type->isIn4Ring, vertex_type->n4Rings, vertex_structure->n4Rings.value()))
            return false;
        bool result = comparePlanarRings(vertex_type->ringsSizes, vertex_type->isInRing, vertex_type->possibleMoreRings, vertex_structure->ringsSizes);
        return result;
    }

    bool GraphVertexComparator::compatibleTypeType(
        const GraphVertex& v1,
        const GraphVertex& v2)
    {
     
        if (v1.isAuxiliaryVertex && v2.isAuxiliaryVertex)
            return true;

        if (v1.isAuxiliaryVertex != v2.isAuxiliaryVertex)
            return false;

        const GraphVertex* vertex_type;
        const GraphVertex* vertex_subtype;

        if (v1.subtypeVertex)
        {
            vertex_subtype = &v1;
            vertex_type = &v2;
        }
        else
        {
            vertex_subtype = &v2;
            vertex_type = &v1;
        }

        // otherwise nothing is more general so it is equally  general or subtype
        if (!vertex_type->anyAtomicNumber) 
        {
            if (vertex_subtype->anyAtomicNumber)
                return false;

            if (!isSubrange(vertex_type->atomicNumberRange, vertex_subtype->atomicNumberRange))
                return false;
        }
        if (!neighborsCompatibleTypeType(vertex_type, vertex_subtype))
            return false;

        if (!comparePlanarity(vertex_type->planar, vertex_subtype->planar))
            return false;

        if (!compareN_RingsTypeType(vertex_type->isIn3Ring, vertex_type->n3Rings, vertex_subtype->isIn3Ring, vertex_subtype->n3Rings))
            return false;

        if (!compareN_RingsTypeType(vertex_type->isIn4Ring, vertex_type->n4Rings, vertex_subtype->isIn4Ring, vertex_subtype->n4Rings))
            return false;
        bool result = comparePlanarRingsTypeType(vertex_type->ringsSizes, vertex_type->isInRing, vertex_type->possibleMoreRings,
            vertex_subtype->ringsSizes, vertex_subtype->isInRing, vertex_subtype->possibleMoreRings);
        return result;
        
    }
    
    bool GraphVertexComparator::neighborsCompatibleTypeType(
        const GraphVertex* vertex_type,
        const GraphVertex* vertex_subtype)
        const
    {

        if (vertex_type->nNeighbors > vertex_subtype->nNeighbors)
            return false;

        if (!vertex_type->additionalNeighboursPossible && vertex_subtype->additionalNeighboursPossible)
            return false;

        if (!vertex_type->additionalNeighboursPossible)
            if(vertex_type->nNeighbors != vertex_subtype->nNeighbors)
                return false;


        if (vertex_type->neighborsAtomicNumbersUniquelyDefined)
        {
            if (!vertex_subtype->neighborsAtomicNumbersUniquelyDefined)
                return false;
            else
                return vertex_subtype->neighborsAtomicNumbers == vertex_type->neighborsAtomicNumbers; 
        }



        //vector<vector<bool> > compatible()


        // check if number of neighbors matches 
        
        if (vertex_type->additionalNeighboursPossible)
        {
            if (!isNeighbourSubset(vertex_type->neighborsAtomicNumbers, vertex_type->neighborsAtomicNumberRanges,
                vertex_type->nNeighbors,
                vertex_subtype->neighborsAtomicNumbers, vertex_subtype->neighborsAtomicNumberRanges, 
                vertex_subtype->nNeighbors, true))
                return false;
        }
        else
        {
            if (vertex_subtype->additionalNeighboursPossible)
                return false;
            else
                if (!isNeighbourSubset(vertex_type->neighborsAtomicNumbers, vertex_type->neighborsAtomicNumberRanges,
                    vertex_type->nNeighbors,
                    vertex_subtype->neighborsAtomicNumbers, vertex_subtype->neighborsAtomicNumberRanges, 
                    vertex_subtype->nNeighbors, false))
                    return false;
        }

        return true;
    }

    bool GraphVertexComparator::neighborsCompatible(
        const GraphVertex *vertexInTypeDefinition,
        const GraphVertex *vertexInStructure)
        const
    {
        
        // check if number of neighbors matches 

        if (!vertexInTypeDefinition->additionalNeighboursPossible)
        {
            if (vertexInTypeDefinition->nNeighbors != vertexInStructure->nNeighbors)
                return false;
        }
        else
            if (vertexInTypeDefinition->nNeighbors > vertexInStructure->nNeighbors)
                return false;

        // if all neighbors are named atoms they are represented as vertex of graph
        // and they will be checked anyway, so thre is no need to do it here

        if (vertexInTypeDefinition->allNonX_NeighborsNamed)
            return true;

        // compare neighbors formulas is they are exactly defined

        if(vertexInTypeDefinition->neighborsAtomicNumbersUniquelyDefined)
            return vertexInTypeDefinition->neighborsAtomicNumbers == vertexInStructure->neighborsAtomicNumbers;

        // extract exatly defined part of formula

        vector<int> zSetToMatch(vertexInStructure->neighborsAtomicNumbers.begin(), vertexInStructure->neighborsAtomicNumbers.end());
        if (!extractSubformula(zSetToMatch, 
                               vector<int>(vertexInTypeDefinition->neighborsAtomicNumbers.begin(), 
                                           vertexInTypeDefinition->neighborsAtomicNumbers.end()))
           )
            return false;
        
        // work with part of formula which should be matched with atomic neighbor ranges

        vector<int> requiredMap;
        return zRangesFitFormula(zSetToMatch, vertexInTypeDefinition->neighborsAtomicNumberRanges, requiredMap);


        //if (!vertexInTypeDefinition->additionalNeighboursPossible)
        //{
        //    if (vertexInTypeDefinition->nNeighbors != vertexInStructure->nNeighbors)
        //        return false;

        //    //if (vertexInTypeDefinition->neighborsWithUndefinedAtomicNumbers.empty())
        //    if (vertexInTypeDefinition->neighborsAtomicNumbersUniquelyDefined)
        //        return (vertexInTypeDefinition->neighborsAtomicNumbers == vertexInStructure->neighborsAtomicNumbers);
        //}
        
        //else
        //    if (vertexInTypeDefinition->neighborsWithUndefinedAtomicNumbers.empty())
        //    {
        //        if (vertexInStructure->nNeighbors < vertexInTypeDefinition->nNeighbors)
        //            return false;
        //        vector<int> undefinedZ_Candidates; // atomic numbers of the atoms in structure which should match vertexInTypeDefinition->neighborsWithUndefinedAtomicNumbers
        //        return extractDefinedNeighbors(vertexInStructure->neighborsAtomicNumbers, vertexInTypeDefinition->neighborsAtomicNumbers, undefinedZ_Candidates);
        //    }

        //vector<int> undefinedZ_Candidates; // atomic numbers of the atoms in structure which should match vertexInTypeDefinition->neighborsWithUndefinedAtomicNumbers


        //if (extractDefinedNeighbors(vertexInStructure->neighborsAtomicNumbers, vertexInTypeDefinition->neighborsAtomicNumbers, undefinedZ_Candidates))
        //    return checkUndefinedZ_Neighbors(undefinedZ_Candidates, vertexInTypeDefinition->neighborsWithUndefinedAtomicNumbers, vertexInTypeDefinition->additionalNeighboursPossible);
        //else
        //    return false;

    }

    bool GraphVertexComparator::extractDefinedNeighbors(
        const std::vector<int> &atomInStructureNeighborsZ,
        const std::vector<int> &atomInTypeNeighborsZ,
        std::vector<int> &remainingVerticesZ)
        const
    {

        int nAtomInTypeNeighbors = atomInTypeNeighborsZ.size();
        int atomIndex;
        vector<int>::iterator it;
        remainingVerticesZ = atomInStructureNeighborsZ;

        for (atomIndex = 0; atomIndex < nAtomInTypeNeighbors; atomIndex++)
        {
            it = find(remainingVerticesZ.begin(), remainingVerticesZ.end(), atomInTypeNeighborsZ[atomIndex]);
            if (it != remainingVerticesZ.end())
                remainingVerticesZ.erase(it);
            else
                return false;
        }
        return true;
    }

    bool GraphVertexComparator::checkUndefinedZ_Neighbors(
        const std::vector<int> &remainingVerticesZ,
        const std::vector<int> &neighborsWithUndefinedAtomicNumbers,
        bool additionalNeighborsPossible)
        const
    {
        int i, nRemainingInStructure, nUndefinedZ_InType;
        /**
        mapping[i]-th vertex from remainingVerticesZ is mapped onto i-th vertex of neighborsWithUndefinedAtomicNumbers
        */
        vector<int> mapping;
        bool foundMatch;

        nRemainingInStructure = remainingVerticesZ.size();
        nUndefinedZ_InType = neighborsWithUndefinedAtomicNumbers.size();

        if (!additionalNeighborsPossible)
        {
            if (nUndefinedZ_InType != nRemainingInStructure)
                return false;
        }
        else
            if (nUndefinedZ_InType > nRemainingInStructure)
                return false;


        mapping.resize(nRemainingInStructure);
        for (i = 0; i < nRemainingInStructure; i++)
            mapping[i] = i;

        // do all permutations of remainingVerticesZ to check if its first 
        // nUndefinedZ_InType components in any of the permutations fits to neighborsWithUndefinedAtomicNumbers

        do
        {
            foundMatch = true;

            for (i = 0; i < nUndefinedZ_InType; i++)
                if (neighborsWithUndefinedAtomicNumbers[i] != 0)
                    if (remainingVerticesZ[mapping[i]] == int(-neighborsWithUndefinedAtomicNumbers[i]))
                    {
                        foundMatch = false;
                        break;
                    }

            if (foundMatch)
                break;
        } while (next_permutation(mapping.begin(), mapping.end()));

        return foundMatch;
    }


    bool  GraphVertexComparator::comparePlanarity(
        Tribool planarityOfAtomInTypeDefinition,
        Tribool planarityOfAtomInStructure)
        const
    {
        if (planarityOfAtomInTypeDefinition == Tribool::Undefined || planarityOfAtomInStructure == Tribool::Undefined)
            return true;
        return (planarityOfAtomInStructure == planarityOfAtomInTypeDefinition);
    }

    bool GraphVertexComparator::comparePlanarRingsTypeType(
        const std::vector<int>& ringsOfAtomInTypeDefinition,
        Tribool typeAtomInRing,
        bool possibleMoreRings,
        const std::vector<int>& ringsOfAtomInSubtypeDefinition,
        Tribool subtypeAtomInRing,
        bool possibleMoreRingsSubtype) 
        const
    {
        if (typeAtomInRing == Tribool::Undefined)
            return true;

        if (subtypeAtomInRing == Tribool::Undefined)
            return false;

        if (subtypeAtomInRing != typeAtomInRing)
            return false;

        if (!possibleMoreRings && possibleMoreRingsSubtype)
            return false;

        if (ringsOfAtomInTypeDefinition.empty())
            return true;
        else
        {
            if (possibleMoreRings)
                return set_theory::is_subset(ringsOfAtomInTypeDefinition, ringsOfAtomInSubtypeDefinition);
            else
                return ringsOfAtomInTypeDefinition == ringsOfAtomInSubtypeDefinition;
        }

    }


    bool GraphVertexComparator::comparePlanarRings(
        const std::vector<int> &ringSizesOfAtomInTypeDefinition,
        Tribool typeAtomInRing,
        bool possibleMoreRings,
        const std::vector<int> &ringsOfAtomInStructure)
        const
    {
        if (typeAtomInRing == Tribool::Undefined)
            return true;

        if (typeAtomInRing == Tribool::False)
            return ringsOfAtomInStructure.empty();
        
        // the third option: typeAtomInRing == AtomRingInfo::IN_RING

        // rings required but the atoms does not belong to any
        if (ringsOfAtomInStructure.empty())
            return false;

        // no specific sizes of rings are specified in atom definition
        // any ring is OK, and there is some ring, since if there is 
        // no ring for atom in structure already false was returned
        if (ringSizesOfAtomInTypeDefinition.empty())
            return true;

        if (!possibleMoreRings)
            return (ringSizesOfAtomInTypeDefinition == ringsOfAtomInStructure);

        vector<int> rings(8, 0);
        int i, nRings;

        nRings = ringsOfAtomInStructure.size();
        for (i = 0; i < nRings; i++)
            if(ringsOfAtomInStructure[i] < 8)
                rings[ringsOfAtomInStructure[i]]++;

        nRings = ringSizesOfAtomInTypeDefinition.size();

        for (i = 0; i < nRings; i++)
            if (rings[ringSizesOfAtomInTypeDefinition[i]] == 0)
                return false;
            else
                rings[ringSizesOfAtomInTypeDefinition[i]]--;

        return true;
    }

    bool GraphVertexComparator::compareN_RingsTypeType(
        Tribool typeAtomInRing,
        std::optional<int> typeN_Rings,
        Tribool subtypeAtomInRing,
        std::optional<int> subtypeN_Rings)
        const
    {

        if (typeAtomInRing == Tribool::Undefined)
            return true;

        if (typeAtomInRing == Tribool::True)
        {
            if (subtypeAtomInRing != Tribool::True)
                return false;
            
            if (!typeN_Rings.has_value())
                return true;
            if (typeN_Rings.value() == -1)
                return true;
            // type has defined number of rings
            // subtype also should
            if (!subtypeN_Rings.has_value())
                return false;
            if (typeN_Rings.value() == -1)
                return false;

            return (typeN_Rings.value() == subtypeN_Rings.value());
        }
        // typeAtomInRing == Tribool::False same should be for subtype
        return (subtypeAtomInRing == Tribool::False);
    }

    
    bool GraphVertexComparator::compareN_Rings(
        Tribool typeAtomInRing,
        std::optional<int> typeN_Rings,
        int nN_RingsInStructure)
        const
    {
        
        if (typeAtomInRing == Tribool::Undefined)
            return true;

        if (typeAtomInRing == Tribool::True)
        {
            if (nN_RingsInStructure == 0)
                return false;
            if (!typeN_Rings.has_value())
                return true;
            if (typeN_Rings.value()==-1)
                return true;
            return (typeN_Rings.value() == nN_RingsInStructure);
        }

        return (nN_RingsInStructure==0);
    }


}
