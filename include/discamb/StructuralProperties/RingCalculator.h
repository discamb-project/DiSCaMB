#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <cstdlib>
#include <map>
#include <set>
#include <vector>

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    class RingCalculator
    {
    public:
        RingCalculator();
		RingCalculator(
            double maxRingPlanarityEsd, 
            double maxAtomPlanarityEsd=0.25, 
            int maxNeighboursCount=3,
            const std::map<std::pair<int, int>, double >& maxBondLengthInAromaticRing = std::map<std::pair<int, int>, double >());

        ~RingCalculator();

        void set(
            double maxRingPlanarityEsd,
            double maxAtomPlanarityEsd = 0.25,
            int maxNeighboursCount = 3,
            const std::map<std::pair<int, int>, double >& maxBondLengthInAromaticRing = std::map<std::pair<int, int>, double >());

        void calculateRings(
            const std::vector<std::vector<int> >& connectivityMatrix,
            const std::vector<Vector3d>& positions,
            const std::vector<double>& planarityEsd,
            int maxRingSize, std::vector<std::vector<int> >& rings,
            std::vector<double>& ringPlanarityEsd,
            const std::vector<int>& atomicNumbers);// = std::vector<int>(),
            //const std::map<std::pair<int, int>, double >& maxInteratomicDistance);// = std::map<std::pair<int, int>, double >());
        /*
        finds planar rings in the structure defined by connectivityMatrix and positions
        if maxInteratomicDistance and atomicNumbers are given then it also checks that all interatomic distances in the ring are below the given threshold for atoms with given atomic numbers
        atoms - indices of atoms to take into account when searching for rings
        rings - output rings as vectors of atom indices, neighbouring indices corresponds to neighboring atoms
        maxInteratomicDistance[{atomicNumber1, atomicNumber2}] = maxDistance - max distance for atoms with given atomic numbers to be considered being in aromatic ring
        atomicNumbers - atomic numbers of atoms corresponding to connectivityMatrix and positions, needed when maxInteratomicDistance is provided
        */
        void calculateRings(
            const std::vector<int>& atoms,
            const std::vector<std::vector<int> >& connectivityMatrix,
            const std::vector<Vector3d>& positions,
            const std::vector<double>& planarityEsd,
            int maxRingSize,
            std::vector<std::vector<int> >& rings,
            std::vector<double>& ringPlanarityEsd,
            const std::vector<int>& atomicNumbers);// = std::vector<int>(),
            //const std::map<std::pair<int, int>, double >& maxInteratomicDistance);

        static void calculate34Rings(const std::vector<std::vector<int> > &connectivityMatrix, 
                                     std::set<std::set<int> > &rings3,
                                     std::set<std::set<int> > &rings4);

        static void calculate34Rings(const std::vector<int> &atoms, 
                                     const std::vector<std::vector<int> > &connectivityMatrix,
                                     std::set<std::set<int> > &rings3,
                                     std::set<std::set<int> > &rings4);

    private:
		double mMaxRingPlanarityEsd = 0.1;
		double mMaxAtomPlanarityEsd = 0.1;
		int mMaxNeighboursCount = 3;
        std::map<std::pair<int, int>, double > mMaxBondLengthInAromaticRing;

        void findRings(int atom, const std::vector<std::vector<int> > &connectivityMatrix, int maxRingSize, std::vector<int> &neighbourhood,
            std::vector<std::vector<int> > &atomRings);
        void findRings(int centralAtom, const std::vector<std::vector<int> > &connectivityMatrix, int maxRingSize, std::vector<std::vector<int> > &atomRings);

        void checkIfRingsDefined(const std::vector<std::vector<int> > &atomRings, const std::vector<std::vector<int> > &rings, std::vector<int> &newRingsIndices);

        void addRingIfNew(std::vector<int> &ring, std::vector<std::vector<int> > &atomRings);

        struct Path {
            std::vector<int> nodes;
            std::vector<int> nodesAsNeighborIndex;
        };

        void firstPath(int centralAtom, const std::vector<std::vector<int> > &connectivityMatrix, int path, Path &firstPath);
        bool nextPath(Path &path, const std::vector<std::vector<int> > &connectivityMatrix, int maxPathSize);
        bool pathContainsRing(Path &path, std::vector<int> &ring);
        void setPathPartToFirst(Path &path, const std::vector<std::vector<int> > &connectivityMatrix, int lastComponentToSave, int maxPathsize);
        void ringCanonicalForm(std::vector<int> &ring);

        std::vector<int>::const_iterator left(const std::vector<int>::const_iterator &it, const std::vector<int> &vec)
        {
            if (it == vec.begin())
                return vec.end() - 1;
            else
                return it - 1;
        }

        std::vector<int>::const_iterator right(const std::vector<int>::const_iterator &it, const std::vector<int> &vec)
        {
            if (it + 1 == vec.end())
                return vec.begin();
            else
                return it + 1;
        }

    };
    /**@}*/



}


