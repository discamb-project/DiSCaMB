#include "discamb/AtomTyping/TypeMatchAlgorithm.h"
#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/AtomTyping/GraphVertexComparator.h"


#include "discamb/BasicUtilities/on_error.h"

#include "argedit.h"
#include "argraph.h"
#include "argedit.h"
#include "vf2_sub_state.h"
#include "match.h"


#include <utility>
#include <algorithm>
#include <limits>


using namespace std;

namespace {

struct VisitorData
{
    bool typeTypeComparison = false;
    discamb::RingMatchChecker *ringChecker = nullptr;
    discamb::RingMatchCheckerTypeType* ringCheckerTypeType = nullptr;
    discamb::AtomType *subtype = nullptr;
    bool ringMatchFoundTypeType;
    std::vector<bool> ringMatchFound;
    std::vector<bool> isTypeVertexAugxiliary;
    std::vector<discamb::MatchMap> matchMaps;
    const discamb::StructureWithDescriptors *structure = nullptr;
    std::vector<int> substructureAtomsToStructureAtoms;
};


bool discamb_match_visitor(
int nVerticesToMatch,
node_id *vertices_1, 
node_id *vertices_2,
void *data)
{
    VisitorData *visitorData = static_cast<VisitorData *>(data);

    visitorData->matchMaps.resize(visitorData->matchMaps.size() + 1);
    visitorData->matchMaps.back().atomMatch.resize(nVerticesToMatch);

    for (int i = 0; i < nVerticesToMatch; i++)
        visitorData->matchMaps.back().atomMatch[vertices_1[i]] = vertices_2[i];


    if (visitorData->typeTypeComparison)
    {
        visitorData->ringMatchFoundTypeType = visitorData->ringCheckerTypeType->match(*visitorData->subtype,
            visitorData->matchMaps.back().atomMatch);
        // do not look further if ring match found
        return visitorData->ringMatchFoundTypeType;
    }
    else
    {
        
        visitorData->ringMatchFound.resize(visitorData->ringMatchFound.size() + 1);

        visitorData->ringMatchFound.back() = visitorData->ringChecker->match(*visitorData->structure,
            visitorData->substructureAtomsToStructureAtoms,
            visitorData->matchMaps.back().atomMatch,
            visitorData->matchMaps.back().labeledRingsMap);
    }

	return false;
}


}

namespace discamb {

    TypeMatchAlgorithm::TypeMatchAlgorithm()
    {
        mTypeDefined = false;
        mMatchOnce = true;
        mTypeRange = 0;

        mAuxiliaryVertexType.isAuxiliaryVertex = true;
        mAuxiliaryVertexStructure.isAuxiliaryVertex = true;

        mAuxiliaryVertexType.atomInStructure = false;
        mAuxiliaryVertexStructure.atomInStructure = true;
    }

    TypeMatchAlgorithm::TypeMatchAlgorithm(
        const AtomType& atomType) : TypeMatchAlgorithm()
    {
        setType(atomType);
    }



    TypeMatchAlgorithm::~TypeMatchAlgorithm()
    {
    }

    std::string TypeMatchAlgorithm::typeId()
        const
    {
        return mTypeId;
    }

    void TypeMatchAlgorithm::setType(
        const AtomType &atomType)
    {
        int atomIndex, nAtoms;
        vector<vector<int> > layers;

        mTypeId = atomType.id;
        mAtomicNumberSet = atomType.atoms[0].atomic_number_range;

        mLCS_Assigner.set(atomType);

        mRingMatchChecker.setAtomType(atomType);

        nAtoms = atomType.atoms.size();
        mAtomTypeGraphVertices.resize(nAtoms);

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {

            int nNeighboursWithDefinedZ_orRangeZ = atomType.atoms[atomIndex].neighborsAtomicNumbers.size() + atomType.atoms[atomIndex].neighborsAtomicNumberRanges.size();
            int nNamedNeighboursWithDefinedZ_orRangeZ = 0;
            for (int nIdx : atomType.connectivity[atomIndex])
                if (!atomType.atoms[nIdx].anyAtomicNumber)
                    nNamedNeighboursWithDefinedZ_orRangeZ++;
            
            bool allZdefinedNeighborsAsVertices = (nNamedNeighboursWithDefinedZ_orRangeZ == nNeighboursWithDefinedZ_orRangeZ);
            mAtomTypeGraphVertices[atomIndex].set(atomType.atoms[atomIndex], allZdefinedNeighborsAsVertices);
        }


        graph_algorithms::breadth_first_search(atomType.connectivity, 0, layers);
        makeVflibGraph(atomType.connectivity, mAtomTypeGraphVertices, layers, mAtomTypeGraph, mAtomTypeGraphEditor);

        

        mAtomTypeGraph->SetNodeComparator(new GraphVertexComparator);
        mTypeRange = layers.size() - 1;               

        mHasSimpleNeighbourFormula = mAtomTypeGraphVertices[0].neighborsAtomicNumbersUniquelyDefined;
            
        if (mHasSimpleNeighbourFormula)
        {
            mNeighbourFormula.clear();
            mNeighbourFormula.insert(atomType.atoms[0].neighborsAtomicNumbers.begin(), atomType.atoms[0].neighborsAtomicNumbers.end());
        }

        mTypeDefined = true;
    }

    bool TypeMatchAlgorithm::generalize(
        const TypeMatchAlgorithm& otherType)
    {


        int nNodesOther = otherType.mAtomTypeGraph->NodeCount();
        for (int i = 0; i < nNodesOther; i++)
            static_cast<GraphVertex*>(otherType.mAtomTypeGraph->GetNodeAttr(i))->subtypeVertex = true;
        
        int nNodesThis = mAtomTypeGraph->NodeCount();
        for (int i = 0; i < nNodesThis; i++)
            static_cast<GraphVertex*>(mAtomTypeGraph->GetNodeAttr(i))->subtypeVertex = false;

        if (nNodesThis > nNodesOther)
            return false;
        
            //matchMaps.clear();
        // checks if atomic number of central atom match
        bool atomicNumberMatch = false;
        for (int z1 : mAtomicNumberSet)
            for (int z2 : otherType.mAtomicNumberSet)
                if (z1 == z2)
                    atomicNumberMatch = true;
        if (!atomicNumberMatch)
            return false;

        if (mHasSimpleNeighbourFormula && otherType.mHasSimpleNeighbourFormula)
            if (mNeighbourFormula != otherType.mNeighbourFormula)
                return false;

        if (!mTypeDefined || !otherType.mTypeDefined)
            return false;

        
// XXXXXXXXXXXXXXX STOPPED HERE XXXXXXXXXXXXXX
        
        //    // make a vflib graph for the tested atom in structure and its neighbourhood

        //    vector<vector<int> > subgraphConnectivity;
        //    vector<GraphVertex> vertices;
        //    vector<int> subgraphVertices;
        //    vector<vector<int> > layers;
        //    int vertex, nVertices;

        //    graph_algorithms::bfs_subgraph(describedStructure.connectivity, atomIndex, mTypeRange, subgraphVertices, subgraphConnectivity, layers);
        //    nVertices = subgraphVertices.size();
        //    vertices.resize(nVertices);

        //    for (vertex = 0; vertex < nVertices; vertex++)
        //        vertices[vertex].set(describedStructure.atomDescriptors[subgraphVertices[vertex]]);

        //    makeVflibGraph(subgraphConnectivity, vertices, layers, mAtomInStructureGraph, mAtomInStructureGraphEditor);

        VisitorData visitorData;
        

        //visitorData.ringChecker = &mRingMatchChecker;
        visitorData.ringCheckerTypeType = &mRingMatchCheckerTypeType;
        visitorData.typeTypeComparison = true;
        //    visitorData.structure = &describedStructure;
        //    visitorData.substructureAtomsToStructureAtoms = subgraphVertices;

        //    for (auto& vertex : mAtomTypeGraphVertices)
        //        visitorData.isTypeVertexAugxiliary.push_back(vertex.isAuxiliaryVertex);

       VF2SubState searchState(mAtomTypeGraph.get(), otherType.mAtomTypeGraph.get());

       bool foundMatch = ::match(&searchState, discamb_match_visitor, &visitorData);
       return foundMatch;



    }

    bool TypeMatchAlgorithm::worksForAtomType(
        const AtomType &atomType)
        const
    {
        on_error::not_implemented(__FILE__, __LINE__);
        return false;
    }

    bool TypeMatchAlgorithm::match(
        int atomIndex,
        const StructureWithDescriptors &describedStructure)
        const
    {
		std::vector<MatchMap> matchMaps;
        return match(atomIndex, describedStructure, matchMaps, false);
    }

    bool TypeMatchAlgorithm::match(
        int atomIndex,
        const StructureWithDescriptors &describedStructure,
		std::vector<MatchMap> &matchMaps)
        const
    {
        return match(atomIndex, describedStructure, matchMaps, true);
    }

    bool TypeMatchAlgorithm::match(
        int atomIndex,
        const StructureWithDescriptors &describedStructure,
        LocalCoordinateSystem<int> &lcs)
        const
    {

		std::vector<MatchMap> matchMaps;

        if (match(atomIndex, describedStructure, matchMaps, true))
        {
            mLCS_Assigner.assignLCS(describedStructure, matchMaps, lcs);
            return true;
        }

        return false;
    }

    bool TypeMatchAlgorithm::match(
        int atomIndex,
        const StructureWithDescriptors &describedStructure,
        LocalCoordinateSystem<int> &lcs,
		std::vector<MatchMap> &matchMaps)
        const
    {

        if (match(atomIndex, describedStructure, matchMaps, true))
        {
            mLCS_Assigner.assignLCS(describedStructure, matchMaps, lcs);
            return true;
        }

        return false;
    }

    bool TypeMatchAlgorithm::match(
        int atomIndex,
        const StructureWithDescriptors &describedStructure,
		std::vector<MatchMap> &matchMaps,
        bool saveMatchMap)
        const
    {
 

		matchMaps.clear();

        if (mAtomicNumberSet.find(describedStructure.atomDescriptors[atomIndex].atomicNumber) == mAtomicNumberSet.end())
            return false;
        
        if (mHasSimpleNeighbourFormula)
            if (mNeighbourFormula != describedStructure.atomDescriptors[atomIndex].neighborsFormula)
                return false;

        if (!mTypeDefined)
            return false;

        // make a vflib graph for the tested atom in structure and its neighbourhood

        vector<vector<int> > subgraphConnectivity;
        vector<GraphVertex> vertices;
        vector<int> subgraphVertices;
        vector<vector<int> > layers;
        int vertex, nVertices;

        graph_algorithms::bfs_subgraph(describedStructure.connectivity, atomIndex, mTypeRange, subgraphVertices, subgraphConnectivity, layers);
        nVertices = subgraphVertices.size();
        vertices.resize(nVertices);

        for (vertex = 0; vertex < nVertices; vertex++)
            vertices[vertex].set(describedStructure.atomDescriptors[subgraphVertices[vertex]]);

        makeVflibGraph(subgraphConnectivity, vertices, layers, mAtomInStructureGraph, mAtomInStructureGraphEditor);

        VisitorData visitorData;
        bool foundMatch;

        visitorData.ringChecker = &mRingMatchChecker;
        visitorData.structure = &describedStructure;
        visitorData.substructureAtomsToStructureAtoms = subgraphVertices;

        for (auto &vertex: mAtomTypeGraphVertices)
            visitorData.isTypeVertexAugxiliary.push_back(vertex.isAuxiliaryVertex);

        VF2SubState searchState(mAtomTypeGraph.get(), mAtomInStructureGraph.get());

        
		
		foundMatch = ::match(&searchState, discamb_match_visitor, &visitorData);
		if (!foundMatch)
			return false;

		foundMatch = false;
		for (auto hasMatch : visitorData.ringMatchFound)
			if (hasMatch)
				foundMatch = true;

		if (!foundMatch)
			return false;

		if (saveMatchMap)
			for (int k = 0; k < visitorData.matchMaps.size(); k++)
				if (visitorData.ringMatchFound[k])
				{
					matchMaps.push_back(visitorData.matchMaps[k]);
					// matched atoms numeration in visitorData.second corresponds to molecule subgraph (made of subgraphVertices)
					// the match numeration is corrected bellow to correspod to whole molecule atom numeration
					
                    for (int i = 0; i < mAtomTypeGraphVertices.size(); i++)
                        matchMaps.back().atomMatch[i] = subgraphVertices[visitorData.matchMaps[k].atomMatch[i]];
                    matchMaps.back().atomMatch.resize(mAtomTypeGraphVertices.size());
					


				}
        return true;

    }


    void TypeMatchAlgorithm::makeVflibGraph(
        const std::vector<std::vector<int> > &connectivity,
        std::vector<GraphVertex> &verticesDescriptor,
        const std::vector<std::vector<int> > &layers,
		std::shared_ptr<ARGraph<GraphVertex, void> > &graph,
		std::shared_ptr<ARGEdit> &graphEditor)
        const
    {
        
        
        // makes vertex to layer map (variable vertexLayer)


        int nVertices, layerIndex, nLayers, layerVertexIndex, layerSize;
        vector<int> vertexLayer;

        nVertices = connectivity.size();
        nLayers = layers.size();


        vertexLayer.resize(nVertices);

        for (layerIndex = 0; layerIndex < nLayers; layerIndex++)
        {
            layerSize = layers[layerIndex].size();
            for (layerVertexIndex = 0; layerVertexIndex < layerSize; layerVertexIndex++)
                vertexLayer[layers[layerIndex][layerVertexIndex]] = layerIndex;
        }


        // add nodes to vflib graph loader

        graphEditor = std::shared_ptr<ARGEdit>(new ARGEdit);
        //ARGEdit argEdit;
        int vertexIndex;

        for (vertexIndex = 0; vertexIndex < nVertices; vertexIndex++)
            graphEditor->InsertNode(&verticesDescriptor[vertexIndex]);


        // add edges to vflib graph loader


        int neighbor, neighborIndex, nNeighbors, nAuxVertices;

        nAuxVertices = 0;

        for (vertexIndex = 0; vertexIndex < nVertices; vertexIndex++)
        {
            nNeighbors = connectivity[vertexIndex].size();

            for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
            {
                neighbor = connectivity[vertexIndex][neighborIndex];

                if (neighbor < vertexIndex) // only once per pair
                {

                    if (vertexLayer[neighbor] == vertexLayer[vertexIndex]) // add auxiliary vertex
                    {
                        //auto node = mAuxiliaryVertex
                        if(verticesDescriptor[0].atomInStructure)
                            graphEditor->InsertNode(&mAuxiliaryVertexStructure);
                        else
                            graphEditor->InsertNode(&mAuxiliaryVertexType);

                        if (static_cast<int>(numeric_limits<node_id>::max()) < nAuxVertices + nVertices)
                            on_error::throwException("too many vertices in graph", __FILE__, __LINE__);
                        else
                        {
                            graphEditor->InsertEdge(static_cast<node_id>(vertexIndex), static_cast<node_id>(nAuxVertices + nVertices), NULL);
                            graphEditor->InsertEdge(static_cast<node_id>(neighbor), static_cast<node_id>(nAuxVertices + nVertices), NULL);
                        }
                        nAuxVertices++;
                    }
                    else
                    {
                        if (vertexLayer[neighbor] > vertexLayer[vertexIndex])
                            graphEditor->InsertEdge(static_cast<node_id>(vertexIndex), static_cast<node_id>(neighbor), NULL);
                        else
                            graphEditor->InsertEdge(static_cast<node_id>(neighbor), static_cast<node_id>(vertexIndex), NULL);
                    }
                }
            }
        }

        // make ARGraph

        graph = std::shared_ptr<ARGraph<GraphVertex, void> >(new ARGraph<GraphVertex, void>(graphEditor.get()));
    }




}

