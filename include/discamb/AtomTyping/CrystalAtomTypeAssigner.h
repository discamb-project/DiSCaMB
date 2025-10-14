#pragma once

#include "MolecularAtomTypeAssigner.h"
#include "AtomType.h"
#include "StructureWithDescriptors.h"
//#include "discamb/CrystalStructure/PredefinedStructuralDescriptors.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/Scattering/AtomRepresentativeInfo.h"

#include <memory>
#include <vector>

namespace discamb {
    /**
* \addtogroup AtomTyping
* @{
*/

    class CrystalAtomTypeAssigner
    {
    public:
        CrystalAtomTypeAssigner();
        ~CrystalAtomTypeAssigner();

        void setAtomTypes(const std::vector<AtomType> &atomTypes);
        void setDescriptorsSettings(const DescriptorsSettings &settings);

        void assign(const Crystal& crystal, std::vector<int> &typeID, std::vector<LocalCoordinateSystem<AtomInCrystalID> > &lcs) const;
        
        //void assign(
        //    const Crystal& crystal, 
        //    const PredefinedStructuralDescriptors &predefinedDescriptors,
        //    std::vector<int>& typeID, 
        //    std::vector<LocalCoordinateSystem<AtomInCrystalID> >& lcs) const;
        
        void assign(
            const Crystal& crystal, 
            const std::vector < std::vector <std::pair<int, std::string> > > &fragmentAtoms,
            const std::vector< std::vector<int> >& atomsToAssign,
            std::vector< std::vector<int> >& typeID, 
            std::vector< std::vector<LocalCoordinateSystem<AtomInCrystalID> > > & lcs) const;
        
        void assign(const Crystal& crystal, std::vector<int>& typeID, std::vector<std::vector<Vector3d> > &lcs) const;
        void assign(const Crystal& crystal, std::vector<int> &typeID, std::vector<LocalCoordinateSystem<AtomInCrystalID> > &lcs,
                    StructureWithDescriptors &descriptors) const;

        void assign(
            const Crystal& crystal, 
            const std::vector<std::vector<AtomInCrystalID> >& fragments,
            const std::vector< std::vector<AtomRepresentativeInfo> > &atomRepresentatives,
            std::vector< std::vector<int> >& typeID, 
            std::vector< std::vector<LocalCoordinateSystem<AtomInCrystalID> > > & lcs) const;

        std::string typeLabel(int typeId) const;
        void printAssignment(std::ostream &out, const Crystal &crystal, const std::vector<int> &typeID,
                                const std::vector< LocalCoordinateSystem<AtomInCrystalID> > &lcs) const;
        void printAssignmentCSV(std::ostream &out, const Crystal &crystal, const std::vector<int> &typeID,
                                const std::vector< LocalCoordinateSystem<AtomInCrystalID> > &lcs) const;
    private:

        MolecularAtomTypeAssigner mAssigner;
        std::vector<std::string> mTypeLabels;// used when printing out assignemtn info
        DescriptorsSettings mDescriptorsSettings;
        int mPlanarRingsRange, mRings34range, mTotalRange, mLabeledNeighbourRange;
        std::vector<AtomType> mAtomTypes;
        void setRanges();

        //AtomTypesDescriptionTree mAtomTypeDescriptionTree;
        
        /**
         mAtomTypeMatchalgorithms[i] - pointer to MatchAlgorithm (in mAtomTypeDescriptionTree) for i-th
                                       atom type (i.e. i == AtomType.number). Can be used for accessing
                                       the types.
        */
        //std::vector<SharedPointer<MatchAlgorithm>::type > mAtomTypeMatchalgorithms;
        //void calculateDescriptors(const std::vector<int> &atomicNumbers, const std::vector<Vector3d> &positions, StructureWithDescriptors &descriptors) const;

        ///**
        //WARNING !!!
        //numeration in match has to correspond to atoms in molecule
        //*/
        //static void localCoordinatesForMatch(int atomIndex,const std::vector<std::vector<int> > &connectivity,const std::vector<std::vector<double> > &bondLengths,
        //	                                 const AtomType, const std::vector<int> &atomMatch,const RingMatch &ringMatch,LocalCoordinateSystem &lcs);

        //void onAtomTypeNotFound(int atomIndex, const std::vector<std::vector<int> > &connectivity, const std::vector<AtomInStructureDescriptors> &atomsDescriptors) const;
        //void onAtomTypeFound(int atomIndex, int atomID) const;
    };
    /**@}*/
}

