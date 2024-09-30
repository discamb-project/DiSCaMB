#include "discamb/AtomTyping/atom_typing_utilities.h"
#include "discamb/AtomTyping/TypeMatchAlgorithm.h"
#include "discamb/AtomTyping/TypeTree.h"
#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/StructuralProperties/RingCalculator.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/HC_Model/DeformationValenceParameters.h"
#include "discamb/BasicChemistry/periodic_table.h"

#include <set>
#include <map>
#include <algorithm>

using namespace std;


namespace discamb {
    namespace atom_typing_utilities {

        int atomTypeRange(
            const AtomType &type,
            int maxPlanarRing,
            int &namedNeighboursRange,
            int &planarRingsRange,
            int &ring34range)
        {
            namedNeighboursRange = 0;
            planarRingsRange = 0;
            ring34range = 0;
            


            int shellIdx, nShells, nAtoms = type.atoms.size();
            vector<vector<int> > graphShells;
            vector<int> n3rings(nAtoms,0), n4rings(nAtoms,0);
            set<set<int> > rings3, rings4;
            graph_algorithms::breadth_first_search(type.connectivity, 0, graphShells);
            int range = 0;
            nShells = graphShells.size();
            RingCalculator::calculate34Rings(type.connectivity, rings3, rings4);

            for (auto &ring : rings3)
                for (int atomIdx : ring)
                    n3rings[atomIdx]++;

            for (auto &ring : rings4)
                for (int atomIdx : ring)
                    n4rings[atomIdx]++;

            range = nShells - 1;

            
            // rings, neighbours, shell

            /*
            1. ringInfo.in3Ring and ringInfo.n3rings
            2. ringInfo.in4Ring and ringInfo.n4rings
            3. planar rings
            4. shell
            5. neighbours
            */
            int atom3ringRange;
            int atom4ringRange;
            int atomPlanarRingRange;

            for(shellIdx=0; shellIdx<nShells; shellIdx++)
                for (int atomIdx : graphShells[shellIdx])
                {
                    AtomRingInfo const &ringInfo = type.atoms[atomIdx].ringInfo;
                    // 3 ring
                    atom3ringRange = 0;
                    // if not in ring we have to check if there is a ring
                    // if undefined we do not check anything
                    // if in ring we check for ring unless nRings = number of rings from connectivity
                    if (ringInfo.in3Ring != Tribool::Undefined &&
                        !(ringInfo.n3rings >= 0 && n3rings[atomIdx] == ringInfo.n3rings))
                        atom3ringRange = shellIdx + 1;
                    // 4 ring
                    atom4ringRange = 0;
                    if (ringInfo.in4Ring != Tribool::Undefined &&
                        !(ringInfo.n4rings >= 0 && n4rings[atomIdx] == ringInfo.n4rings))
                        atom3ringRange = shellIdx + 2;
                    // planar ring
                    atomPlanarRingRange = 0;
                    if (ringInfo.inRing != Tribool::Undefined)
                    {
                        int maxRing;
                        if (ringInfo.inRing == Tribool::False)
                            maxRing = maxPlanarRing;
                        else
                        {
                            if (ringInfo.inAnyAdditionalRing)
                                maxRing = maxPlanarRing;
                            else
                            {
                                int n1, n2;
                                n1 = n2 = 0;
                                if (!ringInfo.labeledContainingRings.empty())
                                    n1 = std::max_element(ringInfo.labeledContainingRings.begin(), ringInfo.labeledContainingRings.end())->first;
                                if (!ringInfo.nonLabeledContainingRings.empty())
                                    n2 = *std::max_element(ringInfo.nonLabeledContainingRings.begin(), ringInfo.nonLabeledContainingRings.end());
                                maxRing = std::max(n1, n2);
                            }
                        }
                        atomPlanarRingRange += (maxRing + 1) / 2 + 1 + shellIdx;
                    }

                    planarRingsRange = max(planarRingsRange, atomPlanarRingRange);
                    ring34range = max( {ring34range, atom3ringRange, atom4ringRange });
                }

            namedNeighboursRange = nShells - 1;
            return max({ planarRingsRange, ring34range, namedNeighboursRange+1 });
        }

        int atomTypesRange(
            const std::vector<AtomType> &types,
            int maxPlanarRing,
            int &namedNeighboursRange,
            int &planarRingsRange,
            int &ring34range)
        {

            int atomRange34ring, atomPlanarRingRange, atomNeighbourRange, atomTotalRange;

            ring34range = 0;
            planarRingsRange = 0;
            namedNeighboursRange = 0;
            int totalRange = 0;

            for (int i = 0; i < types.size(); i++)
            {
                atomTotalRange = atomTypeRange(types[i], maxPlanarRing,
                    atomNeighbourRange, atomPlanarRingRange, atomRange34ring);
                
                namedNeighboursRange = max(atomNeighbourRange, namedNeighboursRange);
                planarRingsRange = max(atomPlanarRingRange, planarRingsRange);
                ring34range = max(atomRange34ring, ring34range);
                totalRange = max(atomTotalRange, totalRange);
            }
            return totalRange;
        }

        void sortTypesByGenarality(
            vector<AtomType>& types)
        {
            TypeTree tree;
            tree.set(types);
            map<int, vector<string> > typesDepthGroups;
            map<string, int> typeLabel2Idx;

            for (int i = 0; i < types.size(); i++)
                typeLabel2Idx[types[i].id] = i;

            for (auto& node : tree.nodes)
                typesDepthGroups[node.second.depth].push_back(node.second.id);

            vector<string> labelsSorted;
            vector<AtomType> typesSorted;
            for (auto& group : typesDepthGroups)
                labelsSorted.insert(labelsSorted.end(), group.second.begin(), group.second.end());
            
            for (auto& label : labelsSorted)
                typesSorted.push_back(types[typeLabel2Idx[label]]);

            types.swap(typesSorted);

            reverse(types.begin(), types.end());
        }

    }

}