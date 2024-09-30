#pragma once

#include "discamb/MathUtilities/Vector3.h"

#include <cstdlib>
#include <vector>
#include <set>

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    class RingCalculator
    {
    public:
        RingCalculator();
		RingCalculator(double maxRingPlanarityEsd, double maxAtomPlanarityEsd=0.25, int maxNeighboursCount=3);
        ~RingCalculator();
		void set(double maxRingPlanarityEsd, double maxAtomPlanarityEsd = 0.25, int maxNeighboursCount = 3);

        void calculateRings(const std::vector<std::vector<int> > &connectivityMatrix, const std::vector<Vector3d> &positions, const std::vector<double> &planarityEsd, int maxRingSize, std::vector<std::vector<int> > &rings, std::vector<double> &ringPlanarityEsd);
        void calculateRings(const std::vector<int> &atoms, const std::vector<std::vector<int> > &connectivityMatrix, const std::vector<Vector3d> &positions, const std::vector<double> &planarityEsd, int maxRingSize, std::vector<std::vector<int> > &rings, std::vector<double> &ringPlanarityEsd);

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


