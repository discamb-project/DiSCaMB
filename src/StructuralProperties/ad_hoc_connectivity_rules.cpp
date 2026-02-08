#include "discamb/StructuralProperties/ad_hoc_connectivity_rules.h"
#include "discamb/BasicUtilities/utilities.h"
#include "discamb/BasicUtilities/on_error.h"

#include <algorithm>

using namespace std;

namespace discamb {

    namespace ad_hoc_connectivity_rules {

        void apply_all(
            const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity)
        {
            disconnect_CC_in_trigonal_bypyramid(positions, atomicNumbers, connectivity);
        }
     
        void apply(const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity,
            Preset preset)
        {
            if (preset == Preset::All)
                apply_all(positions, atomicNumbers, connectivity);
            if (preset == Preset::None)
                return;
        }


        /**
        disconnects possible spurious contact between C1 and C5
           C3
          /  \
       -C1-C2-C5-
          \  /
           C4
        */

        void disconnect_CC_in_trigonal_bypyramid(
            const std::vector<Vector3d>& positions,
            const std::vector<int>& atomicNumbers,
            std::vector<std::vector<int> >& connectivity)
        {
            int atomIdx, nAtoms = positions.size();

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                if(atomicNumbers[atomIdx]==6)
                    if (connectivity[atomIdx].size() > 4)
                    {
                        vector<int> connectedCarbons;
                        for (int neighbour : connectivity[atomIdx])
                            if (atomicNumbers[neighbour] == 6)
                                connectedCarbons.push_back(neighbour);
                        if (connectedCarbons.size() > 3)
                        {
                            vector<pair<double, int> > distanceToNeighbouringC;
                            for (int neighouringC : connectedCarbons)
                            {
                                double dist = (positions[atomIdx] - positions[neighouringC]).norm();
                                distanceToNeighbouringC.push_back({ dist, neighouringC });
                            }
                            
                            int mostDistant = max_element(distanceToNeighbouringC.begin(), distanceToNeighbouringC.end())->second;
                            int nCarbonsConnectingBoth = 0;
                            for (int neighouringC : connectedCarbons)
                                if (utilities::hasElement(connectivity[neighouringC], atomIdx) &&
                                    utilities::hasElement(connectivity[neighouringC], mostDistant))
                                    nCarbonsConnectingBoth++;
                            if (nCarbonsConnectingBoth == 3)
                            {
                                auto it = find(connectivity[atomIdx].begin(), connectivity[atomIdx].end(), mostDistant);
                                if (it != connectivity[atomIdx].end())
                                {
                                    connectivity[atomIdx].erase(it);
                                    it = find(connectivity[mostDistant].begin(), connectivity[mostDistant].end(), atomIdx);
                                    if (it == connectivity[mostDistant].end())
                                        on_error::throwException("nonsymmetric connectivity matrix detected", __FILE__, __LINE__);
                                    connectivity[mostDistant].erase(it);
                                }
                            }

                        }
                    }

            }
        }


    }

}
