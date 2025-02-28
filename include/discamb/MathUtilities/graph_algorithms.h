#pragma once

#include <vector>
#include <cstddef>
#include <functional>
#include <limits>

namespace discamb {

    /**
    * \addtogroup MathUtilities MathUtilities
    * @{
    */

    namespace graph_algorithms
    {
        constexpr int n_pos = std::numeric_limits<int>::max();
        void breadth_first_search(const std::vector<std::vector<int> >& connectivity, int start, std::vector<std::vector<int> >& neighbors, int range = n_pos);

        void breadth_first_search(const std::vector<std::vector<int> >& connectivity, 
            const std::vector<std::vector<bool> > & breakable_edge,int start, std::vector<std::vector<int> >& neighbors, int range = n_pos);

        void bfs_subgraph(const std::vector<std::vector<int> >& connectivity, int start, int range, std::vector<int>& subgraphVertices,
            std::vector<std::vector<int> >& subgraphConnectivity);
        void bfs_subgraph(const std::vector<std::vector<int> >& connectivity, int start, int range, std::vector<int>& subgraphVertices,
            std::vector<std::vector<int> >& subgraphConnectivity, std::vector<std::vector<int> >& subgraphLayers);
        void split(const std::vector<std::vector<int> >& connectivity, std::vector<std::vector<int> >& subGraphs);


    }
    /** @}*/
}




