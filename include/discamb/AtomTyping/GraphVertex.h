#pragma once

//#include "Planarity.h"
#include "AtomInStructureDescriptors.h"
#include "AtomType.h"
#include <vector>
#include <optional>
#include <set>

namespace discamb {
    /**
* \addtogroup AtomTyping
* @{
*/

    struct GraphVertex
    {
        GraphVertex();
        // true if all neighbors with defined atomic number or range of atomic numbers
        // (i.e. not any atomic number, as specified by symbol X) of Z-range
        // are vertices of graph. 
        // Some neigbors my be defined by specifying only its element symbol
        // or symbol of group of elements e.g.
        // C2    CONNECTED_TO  S1,N,X - here N is such a neighbor
        // if there are such neighbors then allZdefinedNeighborsAsVertices is false
        // (there will be no graph vertex for N)
        GraphVertex(const AtomDescriptors& descriptors, bool allZdefinedNeighborsAsVertices);// , bool allZdefinedNeighborsAsVertices);
        GraphVertex(const AtomInStructureDescriptors &descriptors);
        void set(const AtomDescriptors &descriptors, bool allZdefinedNeighborsAsVertices);
        void set(const AtomInStructureDescriptors &descriptors);
        std::string label;// very useful for debugging
        bool anyAtomicNumber = true;
        std::set<int> atomicNumberRange;
        // minimum number of neighbours
        int nNeighbors = 0;
        bool additionalNeighboursPossible = true;
        // is it expressable as an exact formula?
        // false if additionalNeighboursPossible
        // if atomic numbers given by ranges
        // if nNeighbors != neighborsAtomicNumbers.size()

        bool neighborsAtomicNumbersUniquelyDefined = false;

        //std::vector<int> neighborsAtomicNumbers;
        std::multiset<int> neighborsAtomicNumbers;

        std::vector<std::set<int> > neighborsAtomicNumberRanges;
        // applies for type atom, 
        // true if there are no neighbors with specified 
        // atomicNumberRange which are not named
        // e.g. true for neighbors C1, X 
        // false for C1, C

        bool allNonX_NeighborsNamed = true;

        /** 0 - any atom, -Z - any atom with atomic number different than Z*/
        //std::vector<int> neighborsWithUndefinedAtomicNumbers;
        Tribool planar = Tribool::Undefined;
        // concerns planar rings with planar atoms only
        std::vector<int> ringsSizes;
        Tribool isInRing = Tribool::Undefined;
        bool possibleMoreRings = true;
        // concerns 3 and 4 member rings
        Tribool isIn3Ring = Tribool::Undefined;
        Tribool isIn4Ring = Tribool::Undefined;
        std::optional<int> n3Rings;
        std::optional<int> n4Rings;
        // 
        bool atomInStructure = false; // as oposite to atom from atom type definition
        mutable bool subtypeVertex = false; // used for type comparison  
        bool isAuxiliaryVertex = false; // 'auxiliary vertex' - does not correspond to any atom, added to enable making directed graph 'equivalent' to the original undirected one
        
    };

    /**@}*/

}
