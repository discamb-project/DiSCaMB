#include "discamb/MathUtilities/graph_algorithms.h"

using namespace std;

namespace discamb {

    namespace graph_algorithms
    {

        void breadth_first_search(
            const std::vector<std::vector<int> > &connectivity,
            int start,
            std::vector<std::vector<int> > &neighbors,
            int range)
        {
			neighbors.clear();
            bool added;
            int distance, lastLayerSize, layerVertexIndex, vertexIndex, nNeighbors, neighborIndex, neighbor;
            vector<int> vertexCenterDistance(connectivity.size(), -1);


            neighbors.resize(1, vector<int>(1, start));
            vertexCenterDistance[start] = 0;

            added = true;

            for (distance = 1; (distance <= range) && added; distance++)
            {
                added = false;

                lastLayerSize = neighbors[distance - 1].size();

                for (layerVertexIndex = 0; layerVertexIndex < lastLayerSize; layerVertexIndex++)
                {
                    vertexIndex = neighbors[distance - 1][layerVertexIndex];
                    nNeighbors = connectivity[vertexIndex].size();

                    for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
                    {
                        neighbor = connectivity[vertexIndex][neighborIndex];

                        if (vertexCenterDistance[neighbor] == -1)
                        {
                            if (!added)
                            {
                                neighbors.push_back(vector<int>());
                                added = true;
                            }

                            neighbors.back().push_back(neighbor);
                            vertexCenterDistance[neighbor] = distance;
                        }
                    }

                }
            }

        }

        void breadth_first_search(
            const std::vector<std::vector<int> >& connectivity,
            const std::vector<std::vector<bool> >& breakable_edge,
            int start,
            std::vector<std::vector<int> >& neighbors,
            int range)
        {
            neighbors.clear();
            bool added;
            int distance, lastLayerSize, layerVertexIndex, vertexIndex, nNeighbors, neighborIndex, neighbor;
            vector<int> vertexCenterDistance(connectivity.size(), -1);


            neighbors.resize(1, vector<int>(1, start));
            vertexCenterDistance[start] = 0;

            added = true;

            for (distance = 1; (distance <= range) || added; distance++)
            {
                added = false;
                bool inRange = (distance <= range);

                lastLayerSize = neighbors[distance - 1].size();

                for (layerVertexIndex = 0; layerVertexIndex < lastLayerSize; layerVertexIndex++)
                {
                    vertexIndex = neighbors[distance - 1][layerVertexIndex];
                    nNeighbors = connectivity[vertexIndex].size();

                    for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
                    {
                        if (inRange || !breakable_edge[vertexIndex][neighborIndex])
                        {
                            neighbor = connectivity[vertexIndex][neighborIndex];

                            if (vertexCenterDistance[neighbor] == -1)
                            {
                                if (!added)
                                {
                                    neighbors.push_back(vector<int>());
                                    added = true;
                                }

                                neighbors.back().push_back(neighbor);
                                vertexCenterDistance[neighbor] = distance;
                            }
                        }
                    }

                }
            }

        }


        void bfs_subgraph(
            const std::vector<std::vector<int> > &connectivity,
            int start,
            int range,
            std::vector<int> &subgraphVertices,
            std::vector<std::vector<int> > &subgraphConnectivity)
        {
            vector<vector<int> > layers;

            bfs_subgraph(connectivity, start, range, subgraphVertices, subgraphConnectivity, layers);

        }

        void bfs_subgraph(
            const std::vector<std::vector<int> > &connectivity,
            int start,
            int range,
            std::vector<int> &subgraphVertices,
            std::vector<std::vector<int> > &subgraphConnectivity,
            std::vector<std::vector<int> > &subgraphLayers)
        {
            int layerIndex, nLayers, vertexIndex, nVertices, subgraphSize, neighborIndex, nNeighbors, largeGraphVertexIndex;
            int smallGraphNeighborIndex;
            vector<vector<int> > layers;
            vector<int> indexInSubgraph(connectivity.size(), -1);
            breadth_first_search(connectivity, start, layers, range);

            nLayers = layers.size();
            subgraphSize = 0;
            for (layerIndex = 0; layerIndex < nLayers; layerIndex++)
            {
                nVertices = layers[layerIndex].size();
                for (vertexIndex = 0; vertexIndex < nVertices; vertexIndex++)
                {
                    subgraphVertices.push_back(layers[layerIndex][vertexIndex]);
                    indexInSubgraph[layers[layerIndex][vertexIndex]] = static_cast<int>(subgraphSize);
                    subgraphSize++;
                }
            }

            subgraphConnectivity.resize(subgraphSize);
            for (vertexIndex = 0; vertexIndex < subgraphSize; vertexIndex++)
            {
                largeGraphVertexIndex = subgraphVertices[vertexIndex];
                nNeighbors = connectivity[largeGraphVertexIndex].size();
                for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
                {
                    smallGraphNeighborIndex = indexInSubgraph[connectivity[largeGraphVertexIndex][neighborIndex]];
                    if (smallGraphNeighborIndex != -1)
                        //{
                        subgraphConnectivity[vertexIndex].push_back(smallGraphNeighborIndex);
                    //subgraphConnectivity[smallGraphNeighborIndex].push_back(vertexIndex);
                //}
                }
            }

            // reindex atoms in layers to get local (subgraph) indexing

            subgraphLayers.resize(nLayers);

            for (layerIndex = 0; layerIndex < nLayers; layerIndex++)
            {
                nVertices = layers[layerIndex].size();
                subgraphLayers[layerIndex].resize(nVertices);

                for (vertexIndex = 0; vertexIndex < nVertices; vertexIndex++)
                    subgraphLayers[layerIndex][vertexIndex] = indexInSubgraph[layers[layerIndex][vertexIndex]];
            }

        }


    } // namespace graph_algorithms

}
