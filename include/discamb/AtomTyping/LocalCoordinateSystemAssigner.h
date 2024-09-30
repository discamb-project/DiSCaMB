#pragma once

#include "MatchMap.h"
#include "AtomType.h"
#include "LocalCoordinateSystem.h"
#include "StructureWithDescriptors.h"

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    class LocalCoordinateSystemAssigner
    {
    public:

        LocalCoordinateSystemAssigner();
        ~LocalCoordinateSystemAssigner();

        void set(const AtomType &atomType);

		void assignLCS(const StructureWithDescriptors &describedStructure,
			const std::vector<MatchMap> &matchMaps, LocalCoordinateSystem<int> &lcs) const;

		int preferredLcs(const StructureWithDescriptors &describedStructure, const std::vector< LocalCoordinateSystem<int> > &lcss) const;

        void assignLCS(const StructureWithDescriptors &describedStructure,
            const MatchMap &matchMap, LocalCoordinateSystem<int> &lcs) const;
		
        static LcsDirectionType directionType(AtomTypeLCS::LCS_AxisType type);
        //static LocalCoordinateSystem<int>::DirectionType directionType(AtomTypeLCS::LCS_AxisType type);

    private:
        
        std::string mAtomTypeID;
        std::vector<int> mChirality;
        // true if i-th axis is given by AtomTypeLCS::ATOM_LIST or  AtomTypeLCS::LABELED_ATOM
        bool mLabeledAtomsAxis1, mLabeledAtomsAxis2;
        LcsDirectionType mDirection1Type, mDirection2Type;
        std::vector<std::string> mAtomTypeRingLabels;
        // mAtomLabelsMap[label] - givs index of atom in AtomType.atoms which has the label
        std::map<std::string, int> mAtomLabelsMap;
        AtomTypeLCS mLocalCoordinateSystem;
        void setRefPoint(const StructureWithDescriptors &describedStructure, const MatchMap &matchMap, const AtomTypeLCS::LCS_AxisType &axisType,
            const std::vector<int> &refPointCode, const AtomTypeLCS::LCS_AxisType &theOtherAxisType,
            const std::vector<int> &theOtherRefPoint, std::vector<int> &referencePointAtoms) const;
        // as above but without taking care of atoms being already used in definition of LCS (used whe nhe first reference point is defined)
        void setRefPoint(const StructureWithDescriptors &describedStructure, const MatchMap &matchMap, const AtomTypeLCS::LCS_AxisType &axisType,
            const std::vector<int> &refPointCode, std::vector<int> &referencePointAtoms) const;


        void identifyCandidateLscAtoms(const StructureWithDescriptors &describedStructure, const MatchMap &matchMap, const AtomTypeLCS::LCS_AxisType &axisType,
            const std::vector<int> &refPointCode, 
            std::vector<int> &candidateAtoms,
            std::vector<double> &candidateAtomsBondLengths) const;

        // candidateAtoms[i].first atom index (can appear more than once
        // candidateAtoms[i].second distance to bonded named neighbour
        int chooseBestAtom(const StructureWithDescriptors &describedStructure,
                              const MatchMap &matchMap, 
                              const std::vector<int> &candidateAtoms,
                              const std::vector<double> &candidateAtomsBondLengths) const;
        void assignChirality(const MatchMap &matchMap, std::vector<int> &chiralityDefiningAtoms) const;
    };

    /**@}*/
}