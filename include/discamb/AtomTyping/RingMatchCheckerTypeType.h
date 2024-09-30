#pragma once

#include "AtomType.h"
#include "AtomInStructureDescriptors.h"
#include "StructureWithDescriptors.h"
#include "discamb/MathUtilities/InjectionsIterator.h"

#include <map>
#include <string>
#include <optional>

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    class RingMatchCheckerTypeType
    {
    public:
        RingMatchCheckerTypeType();
        ~RingMatchCheckerTypeType();
        void setAtomType(const AtomType &atomType);
        
        bool match(const AtomType &atomType,
                   const std::vector<int> &typeAtoms2SubtypeAtomsMap) const;

        bool uniquelyDefined() const;
    private:
        //AtomType mAtomType;
        int mN_AtomsInAtomType=0;
        std::vector<std::string> mAtomTypeRingLabels;
        std::map<int,std::vector<std::string> > mSizeGroupedAtomTypeRingLabels;
        /**
        mAtomInTypeRrequiredRings[i] - list of indices of rings i-th atom has to be included in.
                                       The indices corresponds to ring names in mRingLabels.
        */
        std::vector<std::vector<int> > mAtomInTypeRrequiredRings;
        std::vector<std::vector<int> > mAtomInTypeForbidenRings;
        std::vector<int> mRingSizes;
        /**
        planarity of ring and its atoms
        */
        std::vector<bool> mRingPlanarity; 
        std::vector<std::vector<int> > mSizeGroupedTypeRingIndices;
        
        /**
        a list of rings for each atom in structure corresponding to atom in type, the rings are rings of 'type'
        corresponding to the current map of rings of structure to rings of type
        */

        
        bool mHasLabeledRings;
        mutable std::vector<InjectionsIterator> mRingInjectionIterators;
        /*
        Indices of the rings in structure which contains atoms corresponding to atoms in atom type definition
        grouped according to ring sizes e.g. if there are 4 such a rings with indices 3 and 5 (five member rings)
        and 4 and 6 (six member rings) then the array is {{3,5},{4,6}}
        */
        //mutable std::vector<std::vector<int> > mSizeGroupedStructureRings;

        /*
        Indices of the rings in structure which contains atoms corresponding to atoms in atom type definition
        grouped according to ring sizes e.g. if there are 4 such a rings with indices 3 and 5 (five member rings)
        and 4 and 6 (six member rings) then the array is {{3,5},{4,6}}
        */
        //mutable std::vector<std::vector<int> > mSizeGroupedStructureRings;



        bool nextRingInjection() const;
        mutable std::vector< std::optional< int > > mRingMap; // maps rings in structure to rings in type definition (-1 if there is no map)
        /**
        i-th entry corresponds to the atom in structure corresponding to i-th atom in type
        the entry is a list of type rings corresponding to structure rings involving the structure atoms
        */
        mutable std::vector<std::vector<int> > mTypeRingsAtStructureAtoms;

        //void updateRingMap(const std::vector<int> &typeAtoms2SubstructureAtomsMap,
        //                   const std::vector<std::vector<int> > &sortedInvolvedStructureRings) const;

        void updateRingMap(std::map<int, std::vector<std::string> > &subtypeSubgraphRingsVec,
            std::map<int, std::vector<std::string> > &typeRingsVec,
            std::map<std::string, std::string> &subtypeTypeRingLabelMap) const;


        //bool checkRingMap(const std::vector<int> &typeAtoms2SubstructureAtomsMap,
        //                  const std::vector<int> &substructureToStructureAtomsMap,
        //                  std::map<int, int> structure2sustructureRing,
        //                  const StructureWithDescriptors &describedStructure) const;

        bool checkRingMap(const AtomType & subtype, const std::map<std::string, std::string> &subtypeTypeRingLabelMap,
            const std::vector<int>& typeAtoms2SubtypeAtomsMap) const;

        /**
        Collects indices of the rings in structure which include atoms corresponding to atoms in atom type definition
        */
        void collectRingsOfMatchingSubtypeAtoms(const AtomType& subtype,
            const std::vector<int> &typeAtoms2SubtypeAtoms,
            std::map<int, std::set<std::string> > &subtypeSubgraphRings) const;

        /**
        sortedInvolvedStructureRings - i-th element of this list is a vector of indices of involvedStructureRingIndices 
                                       pointing to rings of size ringSizes[i]
        */
        void sortInvolvedStructureRingsBySize(const StructureWithDescriptors &describedStructure,
            const std::vector<int> &involvedStructureRingIndices,
            const std::vector<int> &ringSizes,
            std::vector<std::vector<int> > &sortedInvolvedStructureRings) const;

    };
    /**@}*/
}