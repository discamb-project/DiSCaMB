#include "discamb/StructuralProperties/RingCalculator.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/StructuralProperties/structural_properties.h"

#include <set>
#include <iostream>
#include <algorithm>

using namespace std;

namespace {

	void disconnect(
		std::vector<std::vector<int> > &connectivityMatrix,
		const std::vector<int> &atomToDisconnect)
	{
		
		int atomIdx, nAtoms = connectivityMatrix.size();
		std::vector<std::vector<int> > updatedConnectivityMatrix(nAtoms);
		vector<bool> toDisconnect(nAtoms, false);
		for (int atomIdx = 0; atomIdx < atomToDisconnect.size(); atomIdx++)
			toDisconnect[atomToDisconnect[atomIdx]] = true;

		for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
			if (!toDisconnect[atomIdx])
				for (int i = 0; i < connectivityMatrix[atomIdx].size(); i++)
					if (!toDisconnect[connectivityMatrix[atomIdx][i]])
						updatedConnectivityMatrix[atomIdx].push_back(connectivityMatrix[atomIdx][i]);
		updatedConnectivityMatrix.swap(connectivityMatrix);
	}

	void disconnectUnfitAtoms(
		    const std::vector<std::vector<int> > &connectivityMatrix,
		    const std::vector<double> &planarityEsd,
			double atomInRingPlanarityThreshold,
			int maxNeighbourCount,
		    std::vector<std::vector<int> > &updatedConnectivityMatrix)
		{
			updatedConnectivityMatrix = connectivityMatrix;

		    int atomIdx, nAtoms = planarityEsd.size();
			vector<int> atomsToDisconnect;

			//--  max number of neighbours --

			for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
				if (updatedConnectivityMatrix[atomIdx].size() > maxNeighbourCount)
					atomsToDisconnect.push_back(atomIdx);

			disconnect(updatedConnectivityMatrix, atomsToDisconnect);

			//--  planarity --

			atomsToDisconnect.clear();

			for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
				if (planarityEsd[atomIdx] > atomInRingPlanarityThreshold)
					atomsToDisconnect.push_back(atomIdx);

			disconnect(updatedConnectivityMatrix, atomsToDisconnect);

			//-- at least two neighbours

			do {
				atomsToDisconnect.clear();
				for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
					if (updatedConnectivityMatrix[atomIdx].size() == 1)
						atomsToDisconnect.push_back(atomIdx);
				disconnect(updatedConnectivityMatrix, atomsToDisconnect);
			} while (!atomsToDisconnect.empty());
			
		}


    void disconnect34RingSearchUnfitAtoms(
        const std::vector<std::vector<int> > &connectivityMatrix,
        std::vector<std::vector<int> > &updatedConnectivityMatrix)
    {
        updatedConnectivityMatrix = connectivityMatrix;

        int atomIdx, nAtoms = connectivityMatrix.size();
        vector<int> atomsToDisconnect;

        //-- at least two neighbours

        do {
            atomsToDisconnect.clear();
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (updatedConnectivityMatrix[atomIdx].size() == 1)
                    atomsToDisconnect.push_back(atomIdx);
            disconnect(updatedConnectivityMatrix, atomsToDisconnect);
        } while (!atomsToDisconnect.empty());

    }


    // checks if r1 is larger than r2 and contains r2
    bool containsSmallerRing(vector<int> &r1, vector<int> &r2)
    {
        if (r1.size() <= r2.size())
            return false;
        int i, n = r2.size();
        for (i = 0; i < n; i++)
            if (find(r1.begin(), r1.end(), r2[i]) == r1.end())
                return false;
        return true;
    }

    // removes rings which include other rings
    void removeIncludingRings(vector<vector<int> > &rings)
    {
        int nRings;
        set<int> ringsToRemove;
        vector<vector<int> > newRings;

        nRings = rings.size();
        for (int i = 0; i < nRings; i++)
            for (int j = 0; j < nRings; j++)
            {
                if (containsSmallerRing(rings[i], rings[j]))
                    ringsToRemove.insert(i);
                if (containsSmallerRing(rings[j], rings[i]))
                    ringsToRemove.insert(j);
            }

        for (int i = 0; i < nRings; i++)
            if (find(ringsToRemove.begin(), ringsToRemove.end(), i) == ringsToRemove.end())
                newRings.push_back(rings[i]);
        newRings.swap(rings);
    }

}

namespace discamb {

    RingCalculator::RingCalculator()
    {
    }


	RingCalculator::RingCalculator(
		double maxRingPlanarityEsd,
		double maxAtomPlanarityEsd,
		int maxNeighboursCount)
	{
		set(maxRingPlanarityEsd, maxAtomPlanarityEsd, maxNeighboursCount);
	}

    RingCalculator::~RingCalculator()
    {
    }
	
	void RingCalculator::RingCalculator::set(
		double maxRingPlanarityEsd,
		double maxAtomPlanarityEsd,
		int maxNeighboursCount)
	{
		mMaxRingPlanarityEsd = maxRingPlanarityEsd;
		mMaxAtomPlanarityEsd = maxAtomPlanarityEsd;
		mMaxNeighboursCount = maxNeighboursCount;
	}

    void RingCalculator::calculate34Rings(
        const std::vector<std::vector<int> > &connectivityMatrix,
        std::set<std::set<int> > &rings3,
        std::set<std::set<int> > &rings4)
    {
        vector<int> atoms(connectivityMatrix.size());
        for (int i = 0; i < atoms.size(); i++)
            atoms[i] = i;

        calculate34Rings(atoms, connectivityMatrix, rings3, rings4);

    }


    void RingCalculator::calculate34Rings(
        const std::vector<int> &atoms,
        const std::vector<std::vector<int> > &connectivityMatrix,
        std::set<std::set<int> > &rings3,
        std::set<std::set<int> > &rings4)
    {
        rings3.clear();
        rings4.clear();
        std::set<int> ring3, ring4;
        vector<vector<int> > updatedConnectivityMatrix;
        int i, j, nNeighbours, neighbour1, neighbour2, nAtoms = connectivityMatrix.size();
        disconnect34RingSearchUnfitAtoms(connectivityMatrix, updatedConnectivityMatrix);
        const vector<vector<int> > & m = updatedConnectivityMatrix;
        //for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        for(int atomIdx: atoms)
        {
            nNeighbours = m[atomIdx].size();
            if (nNeighbours > 0)
                for (i = 0; i < nNeighbours; i++)
                    for (j = 0; j < i; j++)
                    {
                        neighbour1 = m[atomIdx][i];
                        neighbour2 = m[atomIdx][j];
                        //3
                        for (auto neighbourOfNeighbour1 : m[neighbour1])
                            if (neighbourOfNeighbour1 == neighbour2)
                            {
                                //ring3.insert({ atomIdx, neighbour1, neighbour2 });
                                rings3.insert({ atomIdx, neighbour1, neighbour2 });
                            }
                        //4
                        for (auto neighbourOfNeighbour1 : m[neighbour1])
                            for (auto neighbourOfNeighbour2 : m[neighbour2])
                                if (neighbourOfNeighbour1 == neighbourOfNeighbour2 && neighbourOfNeighbour1 != atomIdx)
                                    // if atomIdx and neighbourOfNeighbour1 are bonded or 
                                    // neighbour1 and neighbour2 are bonded, then the 4 ring contains 3 ring
                                    if (find(m[atomIdx].begin(), m[atomIdx].end(), neighbourOfNeighbour1) == m[atomIdx].end())
                                        if (find(m[neighbour2].begin(), m[neighbour2].end(), neighbour1) == m[neighbour2].end())
                                        {
                                            //ring4.insert({ atomIdx, neighbour1, neighbour2, neighbourOfNeighbour1 });
                                            rings4.insert({ atomIdx, neighbour1, neighbour2, neighbourOfNeighbour1 });
                                        }

                    }

        }

    }



    void RingCalculator::calculateRings(
		const std::vector<int> &atoms,
        const std::vector<std::vector<int> > &connectivityMatrix,
		const std::vector<Vector3d> &positions,
        const std::vector<double> &atomPlanarityEsd,
        int maxRingSize,
        std::vector<std::vector<int> > &rings,
		std::vector<double> &ringPlanarityEsd)
    {
        int atom, nAtoms, maxDistance;
        vector<vector<int> > atomRings;
        vector<int> newRingsIndices;
        vector<vector<int> > updatedConnectivityMatrix;
		vector<int> nRingsPerAtom;
        
		rings.clear();
		ringPlanarityEsd.clear();

        //-------------

		disconnectUnfitAtoms(connectivityMatrix, atomPlanarityEsd, mMaxAtomPlanarityEsd, mMaxNeighboursCount, updatedConnectivityMatrix);
        //updatedConnectivityMatrix = connectivityMatrix;

        //-------------

		nAtoms = atoms.size();// connectivityMatrix.size();
		nRingsPerAtom.assign(connectivityMatrix.size(), 0);
        maxDistance = maxRingSize / 2;

        for (atom = 0; atom < nAtoms; atom++)
        {
			if (nRingsPerAtom[atom] == 1 && updatedConnectivityMatrix[atom].size() == 2)
				continue;
			//cout << atom << "/" << nAtoms << endl;
            findRings(atoms[atom], updatedConnectivityMatrix, maxRingSize, atomRings);
            checkIfRingsDefined(atomRings, rings, newRingsIndices);
			for (int i = 0; i < newRingsIndices.size(); i++)
			{
				rings.push_back(atomRings[newRingsIndices[i]]);
				for (int j = 0; j < atomRings[newRingsIndices[i]].size(); j++)
					nRingsPerAtom[atomRings[newRingsIndices[i]][j]]++;
			}
        }

		// check rings planarity

		vector<vector<int> > planarRings;
		vector<Vector3d> ringAtomsPositions;
		double d;
		for (int i = 0; i < rings.size(); i++)
		{
			ringAtomsPositions.clear();
			for (int j = 0; j < rings[i].size(); j++)
				ringAtomsPositions.push_back(positions[rings[i][j]]);
            auto planarity = structural_properties::planarity(ringAtomsPositions, this->mMaxRingPlanarityEsd, d);
			if (planarity == Tribool::True)
  	        {
			    	planarRings.push_back(rings[i]);
			    	ringPlanarityEsd.push_back(d);
			}
		}

		planarRings.swap(rings);
        removeIncludingRings(rings);
    }

    void RingCalculator::calculateRings(
        const std::vector<std::vector<int> > &connectivityMatrix,
		const std::vector<Vector3d> &positions,
		const std::vector<double> &planarityEsd,
        int maxRingSize,
        std::vector<std::vector<int> > &rings,
		std::vector<double> &ringPlanarityEsd)
    {
		vector<int> atoms(planarityEsd.size());
		generate(atoms.begin(), atoms.end(), [i=int(0)] ( ) mutable {return i++; });
		calculateRings(atoms, connectivityMatrix, positions, planarityEsd, maxRingSize, rings, ringPlanarityEsd);
    }


    void RingCalculator::findRings(
        int centralAtom,
        const std::vector<std::vector<int> > &connectivityMatrix,
        int maxRingSize,
        std::vector<int> &neighbourhood,
        std::vector<std::vector<int> > &atomRings)
    {
        vector<vector<int> > localConnectivity, localRings;
        vector<int>::iterator it;
        int nAtoms, atom, atomIndex, neighbor, nNeighbors, neighborIndex, localIndex, ringIndex, nRings;


        atomRings.clear();

        neighbourhood.push_back(centralAtom);

        nAtoms = neighbourhood.size();

        // make local connectivity matrix (a connectivity matrix only for the central atom and its neighbors

        localConnectivity.resize(nAtoms);

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            atom = neighbourhood[atomIndex];
            nNeighbors = connectivityMatrix[atom].size();
            for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
            {
                neighbor = connectivityMatrix[atom][neighborIndex];
                it = find(neighbourhood.begin(), neighbourhood.end(), neighbor);
                if (it != neighbourhood.end())
                {
                    localIndex = distance(neighbourhood.begin(), it);
                    localConnectivity[atomIndex].push_back(localIndex);
                }
            }
        }

        // find rings in local connectivity matrix

        findRings(nAtoms - 1, localConnectivity, maxRingSize, localRings);

        // translate numeration (from local to whole molecule)

        nRings = localRings.size();
        atomRings.resize(nRings);

        for (ringIndex = 0; ringIndex < nRings; ringIndex++)
        {
            nAtoms = localRings[ringIndex].size();
            atomRings[ringIndex].resize(nAtoms);

            for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
                atomRings[ringIndex][atomIndex] = neighbourhood[localRings[ringIndex][atomIndex]];

            ringCanonicalForm(atomRings[ringIndex]);
        }


    }


    void RingCalculator::findRings(
        int centralAtom,
        const vector<vector<int> > &connectivityMatrix,
        int maxRingSize,
        vector<vector<int> > &atomRings)
    {
        Path path;
        vector<int> ring;
        int maxPathSize = maxRingSize + 1;

        firstPath(centralAtom, connectivityMatrix, maxPathSize, path);

        do {
            if (pathContainsRing(path, ring))
                addRingIfNew(ring, atomRings);
        } while (nextPath(path, connectivityMatrix, maxPathSize));

    }

    void RingCalculator::checkIfRingsDefined(
        const std::vector<std::vector<int> > &newRings,
        const std::vector<std::vector<int> > &alreadyDefinedRings,
        std::vector<int> &uniqueNewRingsIndices)
    {
        int nNewRings, newRingIndex;
        uniqueNewRingsIndices.clear();
        nNewRings = newRings.size();
        for (newRingIndex = 0; newRingIndex < nNewRings; newRingIndex++)
            if (find(alreadyDefinedRings.begin(), alreadyDefinedRings.end(), newRings[newRingIndex]) == alreadyDefinedRings.end())
                uniqueNewRingsIndices.push_back(newRingIndex);
    }


    void RingCalculator::firstPath(
        int centralAtom,
        const std::vector<std::vector<int> > &connectivityMatrix,
        int maxSize,
        Path &firstPath)
    {
        if (maxSize == 0)
            return;

        firstPath.nodesAsNeighborIndex.assign(1, -1);
        firstPath.nodes.assign(1, centralAtom);

        setPathPartToFirst(firstPath, connectivityMatrix, 0, maxSize);
    }

    bool RingCalculator::pathContainsRing(
        Path &path,
        std::vector<int> &ring)
    {
        int index, pathSize;
        bool containsRing = false;

        ring.clear();

        pathSize = path.nodes.size();

        if (pathSize < 4)
            return false;

        for (index = 3; index < pathSize; index++)
            if (path.nodes[index] == path.nodes[0])
                if (path.nodes[index - 1] != path.nodes[1])
                {
                    ring.clear();
                    ring.insert(ring.end(), path.nodes.begin(), path.nodes.begin() + index);
                    ringCanonicalForm(ring);
                    return true;
                }

        return false;
    }

    void RingCalculator::addRingIfNew(
        std::vector<int> &ring,
        std::vector<std::vector<int> > &atomRings)
    {
        if (find(atomRings.begin(), atomRings.end(), ring) == atomRings.end())
            atomRings.push_back(ring);
    }

    bool RingCalculator::nextPath(
        Path &path,
        const std::vector<std::vector<int> > &connectivityMatrix,
        int maxPathSize)
    {
        int pathComponentToChange = int(path.nodes.size()) - 1;
        int parentNode, increment;

        while (pathComponentToChange > 0)
        {
            parentNode = path.nodes[pathComponentToChange - 1];
            if (path.nodesAsNeighborIndex[pathComponentToChange] + 1 < connectivityMatrix[parentNode].size())
            {
                increment = 1;
                if (pathComponentToChange > 1) // in this case it should be checked if candidate node is not its own 'grandparent'
                {
                    if (path.nodes[pathComponentToChange - 2] == connectivityMatrix[parentNode][path.nodesAsNeighborIndex[pathComponentToChange] + 1]) // its own grandparent
                    {
                        if (path.nodesAsNeighborIndex[pathComponentToChange] + 2 < connectivityMatrix[parentNode].size()) // then we take the next neighbor if possible
                            increment = 2;
                        else // there is no next neighbor to take - the component can be changed we go one level down
                        {
                            pathComponentToChange--;
                            continue;
                        }
                    }
                }

                path.nodesAsNeighborIndex[pathComponentToChange] += increment;
                path.nodes[pathComponentToChange] = connectivityMatrix[parentNode][path.nodesAsNeighborIndex[pathComponentToChange]];
                setPathPartToFirst(path, connectivityMatrix, pathComponentToChange, maxPathSize);
                return true;
            }
            else
                pathComponentToChange--;
        }

        return false;
    }


    void RingCalculator::setPathPartToFirst(
        Path &path,
        const std::vector<std::vector<int> > &connectivityMatrix,
        int lastComponentToSave,
        int maxPathSize)
    {
        int length, firstChildIndex, savedPathSize;

        if (maxPathSize < lastComponentToSave)
            on_error::throwException("bug", __FILE__, __LINE__);

        if (path.nodes.empty())
            on_error::throwException("bug", __FILE__, __LINE__);

        if (lastComponentToSave == maxPathSize - 1)
            return;

        if (path.nodes.size() <= lastComponentToSave)
            on_error::throwException("bug", __FILE__, __LINE__);

        savedPathSize = lastComponentToSave + 1;

        path.nodes.resize(savedPathSize);
        path.nodesAsNeighborIndex.resize(savedPathSize);

        if (lastComponentToSave == 0)
        {
            if (connectivityMatrix[path.nodes[0]].empty())
                return;
            path.nodesAsNeighborIndex.push_back(0);
            path.nodes.push_back(connectivityMatrix[path.nodes[0]][0]);
            length = 2;
        }
        else
            length = savedPathSize;


        if (maxPathSize == 2)
            return;



        do
        {
            int currentNodeIndex = path.nodes.back();
            int parentNode = path.nodes[length - 2];

            if (connectivityMatrix[currentNodeIndex].size() == 1)
                return;

            connectivityMatrix[currentNodeIndex][0] == parentNode ? firstChildIndex = 1 : firstChildIndex = 0;

            path.nodes.push_back(connectivityMatrix[currentNodeIndex][firstChildIndex]);
            path.nodesAsNeighborIndex.push_back(firstChildIndex);
            length++;
        } while (length < maxPathSize);

    }

    void RingCalculator::ringCanonicalForm(
        std::vector<int> &ring)
    {
        vector<int> canonical;
        vector<int>::const_iterator it;
        int ringSize = ring.size();
        int index, left, right;
        bool goRight;

        if (ringSize < 3)
            on_error::throwException("bug - trying to write not a ring in ring canonical form", __FILE__, __LINE__);

        it = min_element(ring.begin(), ring.end());

        right = *this->right(it, ring);
        left = *this->left(it, ring);

        goRight = (left > right);
        canonical.resize(ringSize);
        canonical[0] = *it;

        for (index = 1; index < ringSize; index++)
        {
            if (goRight)
                it = this->right(it, ring);
            else
                it = this->left(it, ring);
            canonical[index] = *it;
        }

        ring = canonical;
    }

}

