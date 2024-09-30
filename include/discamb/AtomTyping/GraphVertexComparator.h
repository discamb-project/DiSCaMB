#pragma once

#include "GraphVertex.h"

#include "argedit.h"

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    class GraphVertexComparator : public AttrComparator
    {
    public:
        GraphVertexComparator();
        virtual ~GraphVertexComparator();
        virtual bool compatible(void *vertex1, void *vertex2);
        bool compatible(const GraphVertex &v1, const GraphVertex &v2);
    private:
        bool neighborsCompatible(const GraphVertex *vertexInTypeDefinition, const GraphVertex *vertexInStructure) const;
        bool neighborsCompatibleTypeType(const GraphVertex* vertex_type, const GraphVertex* vertex_subtype) const;
        /** try to extract from vertexInStructure neighborsAtomicNumbers part corresponding to vertexInTypeDefinition->neighborsAtomicNumbers
         the remaining part is returned to remainingVerticesZ and can be used to check against vertexInTypeDefinition->neighborsWithUndefinedAtomicNumbers
         */
        
        bool extractDefinedNeighbors(const std::vector<int> &atomInStructureNeighborsZ, const std::vector<int> &atomInTypeNeighborsZ, std::vector<int> &remainingVerticesZ) const;
        bool checkUndefinedZ_Neighbors(const std::vector<int> &remainingVerticesZ, const std::vector<int> &neighborsWithUndefinedAtomicNumbers, bool additionalNeighborsPossible) const;
        bool comparePlanarity(Tribool planarityOfAtomInTypeDefinition, Tribool planarityOfAtomInStructure) const;
        bool comparePlanarRings(const std::vector<int> &ringsOfAtomInTypeDefinition,
                                Tribool typeAtomInRing,
                                bool possibleMoreRings, 
                                const std::vector<int> &ringsOfAtomInStructure) const;
        bool comparePlanarRingsTypeType(const std::vector<int>& ringsOfAtomInTypeDefinition,
            Tribool typeAtomInRing,
            bool possibleMoreRings,
            const std::vector<int>& ringsOfAtomInSubtypeDefinition,
            Tribool subtypeAtomInRing,
            bool possibleMoreSubtype) const;

        // compare 3 member or 4 member rings
        bool compareN_Rings(Tribool typeAtomInRing, std::optional<int> typeN_Rings, int nN_RingsInStructure) const;
        bool compareN_RingsTypeType(Tribool typeAtomInRing, std::optional<int> typeN_Rings, Tribool subtypeAtomInRing, std::optional<int> subtypeN_Rings) const;

        // both vertices are type vertices, used for types comparison
        bool compatibleTypeType(const GraphVertex& v1, const GraphVertex& v2);
    };

    /**@}*/
}
