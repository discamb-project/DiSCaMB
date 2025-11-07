#pragma once


//#include "AtomTypesAssignment.h"
#include "AtomType.h"
#include "TypeMatchAlgorithm.h"
#include "StructureWithDescriptors.h"
#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <memory>

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    class MolecularAtomTypeAssigner
    {
    public:
        MolecularAtomTypeAssigner();
        ~MolecularAtomTypeAssigner();

        void setAtomTypes(const std::vector<AtomType> &atomTypes);
        void setDescriptorsSettings(const DescriptorsSettings &settings);
        void assign(const std::vector<int> &atomicNumbers, const std::vector<Vector3d> &positions, const std::vector<std::string> &labels, std::vector<int> &typeID, std::vector<LocalCoordinateSystem<int> > &lcs) const;
        void assign(const std::vector<int> &atomicNumbers, const std::vector<Vector3d> &positions, const std::vector<std::string> &labels, std::vector<int> &atomsToAssign, std::vector<int> &typeID, std::vector<LocalCoordinateSystem<int> > &lcs) const;
        void assign(const StructureWithDescriptors &descriptors, std::vector<int> &typeID, std::vector<LocalCoordinateSystem<int> > &lcs) const;
        void assign(const StructureWithDescriptors &descriptors, const std::vector<int> &atomsToAssign, std::vector<int> &typeID, std::vector<LocalCoordinateSystem<int> > &lcs) const;

        void assign_all_possible(const StructureWithDescriptors& descriptors, std::vector< std::vector<int> >& typeID) const;


    private:
        DescriptorsSettings mDescriptorsSettings;
        //AtomTypesDescriptionTree mAtomTypeDescriptionTree;
        /**
         mAtomTypeMatchalgorithms[i] - pointer to MatchAlgorithm (in mAtomTypeDescriptionTree) for i-th
                                       atom type (i.e. i == AtomType.number). Can be used for accessing
                                       the types.
        */
        std::vector<TypeMatchAlgorithm> mAtomTypeMatchalgorithms;
        
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

