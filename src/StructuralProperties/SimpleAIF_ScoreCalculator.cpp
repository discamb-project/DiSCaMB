#include "discamb/StructuralProperties/SimpleAIF_ScoreCalculator.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/MathUtilities/graph_algorithms.h"

using namespace std;

namespace discamb {

	SimpleAIF_ScoreCalculator::SimpleAIF_ScoreCalculator() 
	{
	}

	SimpleAIF_ScoreCalculator::SimpleAIF_ScoreCalculator(
		const nlohmann::json& data)
	{

	}
	
	SimpleAIF_ScoreCalculator::~SimpleAIF_ScoreCalculator()
	{

	}
	
	void SimpleAIF_ScoreCalculator::assignScore(
		const Crystal& crystal,
		const std::vector<std::vector<std::pair<std::string, std::string> > >& clusters,
		std::vector < std::vector < double> >& scores)
	{
		UnitCellContent ucContent;
		ucContent.set(crystal);
		vector<vector<UnitCellContent::AtomID> > ucConnectivity;
		
		structural_properties::calcUnitCellConnectivity(ucContent, ucConnectivity, 0.4);

		scores.clear();
		int clusterIdx, nClusters = clusters.size();
		scores.resize(nClusters);
		for (clusterIdx = 0; clusterIdx < nClusters; clusterIdx++)
			assignScore2(ucContent, ucConnectivity, clusters[clusterIdx], scores[clusterIdx]);
	}

	void SimpleAIF_ScoreCalculator::assignScore(
		UnitCellContent &ucContent,
		std::vector<std::vector<UnitCellContent::AtomID> >& ucConnectivity,
		const std::vector<std::pair<std::string, std::string> >& cluster,
		std::vector < double>& scores)
	{
		
		
		int atomIdx, nAtoms = cluster.size();
		const int maxNeighbour = nAtoms + 1;
		vector<int> distanceToExternalAtom(nAtoms, maxNeighbour);
		vector<vector<int> > connectivity;
		vector<bool> isExternal(nAtoms);
		vector<int> boundedAtomIdx;
		vector<UnitCellContent::AtomID> atoms(nAtoms);


		for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
		{
			if (cluster[atomIdx].first.find('@') == string::npos)
			{
				ucContent.findAtom(cluster[atomIdx].first, cluster[atomIdx].second, atoms[atomIdx]);
				isExternal[atomIdx] = false;
			}
			else
			{
				auto& s = cluster[atomIdx].second;
				string atomLabel, symmOp;
				crystal_structure_utilities::splitIntoAtomAndSymmOp(s, atomLabel, symmOp);
				ucContent.findAtom(atomLabel, symmOp, atoms[atomIdx]);
				isExternal[atomIdx] = true;
			}
		}

		structural_properties::calculateConnectivity(ucConnectivity, atoms, connectivity);
		vector<vector<int> > neighbourShells;
		for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
			if (!isExternal[atomIdx])
			{
				graph_algorithms::breadth_first_search(connectivity, atomIdx, neighbourShells);
				for (int i = 0; i < neighbourShells.size(); i++)
				{
					if (distanceToExternalAtom[atomIdx] < i)
						break;

					for (int neighbourAtomIdx : neighbourShells[i])
						if (isExternal[neighbourAtomIdx])
						{
							distanceToExternalAtom[atomIdx] = i;
							break;
						}
				}
			}

		scores.resize(nAtoms);
		for (int i = 0; i < nAtoms; i++)
			if(isExternal[i])
				scores[i] = -1;
			else
				scores[i] = (double)distanceToExternalAtom[i];
		
	}


    void SimpleAIF_ScoreCalculator::assignScore2(
        UnitCellContent& ucContent,
        std::vector<std::vector<UnitCellContent::AtomID> >& ucConnectivity,
        const std::vector<std::pair<std::string, std::string> >& cluster,
        std::vector < double>& scores)
    {

        
        int atomIdx, nAtoms = cluster.size();
        scores.assign(nAtoms, -1);
        const int maxNeighbour = nAtoms + 1;
        int distanceToExternalAtom; 
        vector<vector<int> > connectivity;
        vector<bool> isExternal(nAtoms);
        vector<int> boundedAtomIdx;
        // nonCapAtoms
        vector<UnitCellContent::AtomID> nonCapAtoms;
        
        vector<int> internalAtoms;
        vector<int> setToSubsetIdx(nAtoms, -1);


        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            if (cluster[atomIdx].first.find('@') == string::npos)
            {
                nonCapAtoms.resize(nonCapAtoms.size() + 1);
                ucContent.findAtom(cluster[atomIdx].first, cluster[atomIdx].second, nonCapAtoms.back());
                isExternal[atomIdx] = false;
                setToSubsetIdx.push_back(internalAtoms.size());
                internalAtoms.push_back(atomIdx);
            }
        }

        structural_properties::calculateConnectivity(ucConnectivity, nonCapAtoms, connectivity);
        vector<vector<int> > neighbourShells;
        int nMissingNeighbors;
        for (int atomInSubsetIdx = 0; atomInSubsetIdx < nonCapAtoms.size(); atomInSubsetIdx++)
        {
            graph_algorithms::breadth_first_search(connectivity, atomInSubsetIdx, neighbourShells);
            distanceToExternalAtom = maxNeighbour;
            nMissingNeighbors = 0;
            for (int i = 0; i < neighbourShells.size(); i++)
            {
                if (distanceToExternalAtom < i)
                    break;
                nMissingNeighbors = 0;

                for (int neighbourAtomIdx : neighbourShells[i])
                {
                    int nNeighborsCrystal = ucConnectivity[nonCapAtoms[neighbourAtomIdx].atomIndex].size();
                    int nNeighborsCluster = connectivity[neighbourAtomIdx].size();
                    if (nNeighborsCrystal != nNeighborsCluster)
                    {
                        distanceToExternalAtom = i;
                        nMissingNeighbors += nNeighborsCrystal - nNeighborsCluster;
                    }
                }
            }
            scores[internalAtoms[atomInSubsetIdx] ] = distanceToExternalAtom - 1e-8 * nMissingNeighbors;
        }

    }
}

