#pragma once


#include "AtomType.h"
#include "AtomInStructureDescriptors.h"
#include "GraphVertex.h"
#include "RingMatchChecker.h"
#include "RingMatchCheckerTypeType.h"
#include "LocalCoordinateSystemAssigner.h"

template <class Node, class Edge>
class ARGraph;
class ARGEdit;
      
#include <memory>

namespace discamb {

    /**
    * \addtogroup AtomTyping
    * @{
    */

    class TypeMatchAlgorithm
    {
    public:
        TypeMatchAlgorithm();
        TypeMatchAlgorithm(const AtomType &atomType);
        virtual ~TypeMatchAlgorithm();

        virtual void setType(const AtomType &atomType);

        bool generalize(const TypeMatchAlgorithm & typeMatchAlgorithm);

        virtual bool worksForAtomType(const AtomType &atomType) const;

        virtual bool match(int atomIndex, const StructureWithDescriptors &describedStructure) const;
        virtual bool match(int atomIndex, const StructureWithDescriptors &describedStructure, std::vector<MatchMap> &matchMaps) const;
        virtual bool match(int atomIndex,const StructureWithDescriptors &describedStructure,LocalCoordinateSystem<int> &lcs) const;
        virtual bool match(int atomIndex,const StructureWithDescriptors &describedStructure,LocalCoordinateSystem<int> &lcs,std::vector<MatchMap> &matchMaps) const;

        bool checkRingAssignment(const std::vector<int> &typeAtoms2SubstructureAtomsMap);
        std::string typeId() const;
    private:
        std::string mTypeId;
        std::set<int> mAtomicNumberSet;
        LocalCoordinateSystemAssigner mLCS_Assigner;
        bool mHasSimpleNeighbourFormula = false;
        //std::vector<int> mNeighbourFormula;
        std::multiset<int> mNeighbourFormula;

        bool match(int atomIndex, const StructureWithDescriptors &describedStructure, std::vector<MatchMap> &matchMap, bool saveMatchMap) const;

        mutable std::vector<int> mTypeAtoms2SubstructureAtomsMap;
        int mTypeRange;
        bool mTypeDefined;
        bool mMatchOnce;

        //

        mutable std::shared_ptr<ARGraph<GraphVertex, void> > mAtomTypeGraph;
        mutable std::shared_ptr<ARGraph<GraphVertex, void> > mAtomInStructureGraph;

        mutable std::shared_ptr<ARGEdit> mAtomTypeGraphEditor;
        mutable std::shared_ptr<ARGEdit> mAtomInStructureGraphEditor;

        std::vector<GraphVertex> mAtomTypeGraphVertices;

        void makeVflibGraph(const std::vector<std::vector<int> > &connectivity, std::vector<GraphVertex> &verticesDescriptor,
            const std::vector<std::vector<int> > &layers, std::shared_ptr<ARGraph<GraphVertex, void> > &graph,
			std::shared_ptr<ARGEdit> &graphEditor) const;
        mutable RingMatchChecker mRingMatchChecker;
        mutable RingMatchCheckerTypeType mRingMatchCheckerTypeType;
        mutable GraphVertex mAuxiliaryVertexType;
        mutable GraphVertex mAuxiliaryVertexStructure;

        // 

    };


    /**@}*/
}

