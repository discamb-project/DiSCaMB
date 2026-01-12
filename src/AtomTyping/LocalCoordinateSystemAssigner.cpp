#include "discamb/AtomTyping/LocalCoordinateSystemAssigner.h"

#include "discamb/BasicUtilities/on_error.h"

#include <algorithm>
#include <set>

using namespace std;

namespace discamb {

    LocalCoordinateSystemAssigner::LocalCoordinateSystemAssigner()
    {
        mAtomTypeID = std::string("NONE");
        mLabeledAtomsAxis1 = false;
        mLabeledAtomsAxis2 = false;
        mDirection1Type = LcsDirectionType::ANY_ORTHOGONAL;
        mDirection2Type = LcsDirectionType::ANY_ORTHOGONAL;
    }

    LocalCoordinateSystemAssigner::~LocalCoordinateSystemAssigner()
    {
    }

    LcsDirectionType LocalCoordinateSystemAssigner::directionType(
		AtomTypeLCS::LCS_AxisType type)
	{
		if (type == AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL)
			return LcsDirectionType::ANY_ORTHOGONAL;
		if (type == AtomTypeLCS::LCS_AxisType::AVERAGE_DIRECTION)
			return LcsDirectionType::AVERAGE_DIRECTION;
		if (type == AtomTypeLCS::LCS_AxisType::NOT_SET)
			return LcsDirectionType::NOT_SET;
		return LcsDirectionType::AVERAGE_POSITION;
	}

    void LocalCoordinateSystemAssigner::set(
        const AtomType &atomType)
    {
        mLocalCoordinateSystem = atomType.localCoordinateSystem;
        mAtomTypeRingLabels = atomType.ringLabels;

        mLabeledAtomsAxis1 = (mLocalCoordinateSystem.lcs_axis_type_1 == AtomTypeLCS::LCS_AxisType::ATOM_LIST || mLocalCoordinateSystem.lcs_axis_type_1 == AtomTypeLCS::LCS_AxisType::LABELED_ATOM);
        mLabeledAtomsAxis2 = (mLocalCoordinateSystem.lcs_axis_type_2 == AtomTypeLCS::LCS_AxisType::ATOM_LIST || mLocalCoordinateSystem.lcs_axis_type_2 == AtomTypeLCS::LCS_AxisType::LABELED_ATOM);

        mAtomTypeID = atomType.id;

		mDirection1Type = directionType(mLocalCoordinateSystem.lcs_axis_type_1);
		mDirection2Type = directionType(mLocalCoordinateSystem.lcs_axis_type_2);
        mChirality = atomType.chirality;
    }

	void LocalCoordinateSystemAssigner::assignLCS(
		const StructureWithDescriptors &describedStructure,
		const std::vector<MatchMap> &matchMaps,
		LocalCoordinateSystem<int> &lcs)
		const
	{
		int mapIdx, nMaps = matchMaps.size();
		vector<LocalCoordinateSystem<int> > lcss(nMaps);
		for (mapIdx = 0; mapIdx < nMaps; mapIdx++)
			assignLCS(describedStructure, matchMaps[mapIdx], lcss[mapIdx]);
		int preferred = preferredLcs(describedStructure, lcss);
		lcs = lcss[preferred];
		lcs.direction1_type = mDirection1Type;
		lcs.direction2_type = mDirection2Type;
        assignChirality(matchMaps[preferred], lcs.chirality);
	}

    void LocalCoordinateSystemAssigner::assignChirality(
        const MatchMap &matchMap,
        std::vector<int> &chiralityDefiningAtoms)
        const
    {
        int i, n = mChirality.size();
        chiralityDefiningAtoms.resize(n);
        for (i = 0; i < n; i++)
            chiralityDefiningAtoms[i] = matchMap.atomMatch[mChirality[i]];
    }

	int LocalCoordinateSystemAssigner::preferredLcs(
		const StructureWithDescriptors &describedStructure,
		const std::vector< LocalCoordinateSystem<int> > &lcss)
    const
	{
		if (lcss.empty())
			return 0;
		int centralAtom = lcss[0].centralAtom;
		int lcsIdx, nLcs = lcss.size();
		
		Vector3d r, rCentralAtom = describedStructure.atomDescriptors[centralAtom].position;
		/*
		for direction1 of type AVERAGE_POSITION or AVERAGE_DIRECTION 
		the two first components of preferrenceMeasurePerLcs[i].first
		are (1) sum of minus atomic numbers and (2) sum of bond lengths
		lowest values are preferred
        similarily for direction2 and the next two componts of 
        preferrenceMeasurePerLcs[i].first
		*/
		vector<pair<vector<double>, int> > preferrenceMeasurePerLcs(nLcs, { {0.0,0.0,0.0,0.0},0 });
		
		for (lcsIdx = 0; lcsIdx < nLcs; lcsIdx++)
		{
			preferrenceMeasurePerLcs[lcsIdx].second = lcsIdx;

			//if (lcss[lcsIdx].direction1_type == LcsDirectionType::AVERAGE_POSITION ||
			//	lcss[lcsIdx].direction1_type == LcsDirectionType::AVERAGE_DIRECTION)
            if (mDirection1Type == LcsDirectionType::AVERAGE_POSITION ||
                mDirection1Type == LcsDirectionType::AVERAGE_DIRECTION)
			{
				for (auto idx : lcss[lcsIdx].refPoint_1)
				{
					preferrenceMeasurePerLcs[lcsIdx].first[0] += describedStructure.atomDescriptors[idx].atomicNumber;

					r = describedStructure.atomDescriptors[idx].position - rCentralAtom;
					preferrenceMeasurePerLcs[lcsIdx].first[1] -= sqrt(r*r);
				}
			}

            if (mDirection2Type == LcsDirectionType::AVERAGE_POSITION ||
                mDirection2Type == LcsDirectionType::AVERAGE_DIRECTION)
            {
                for (auto idx : lcss[lcsIdx].refPoint_2)
                {
                    preferrenceMeasurePerLcs[lcsIdx].first[2] += describedStructure.atomDescriptors[idx].atomicNumber;

                    r = describedStructure.atomDescriptors[idx].position - rCentralAtom;
                    preferrenceMeasurePerLcs[lcsIdx].first[3] -= sqrt(r*r);
                }
            }

		}
		int preferred = max_element(preferrenceMeasurePerLcs.begin(), preferrenceMeasurePerLcs.end())->second;
 		return preferred;
	}

    void LocalCoordinateSystemAssigner::assignLCS(
        const StructureWithDescriptors &describedStructure,
        const MatchMap &matchMap,
        LocalCoordinateSystem<int> &lcs)
        const
    {
		/*
		enum DirectionType { ANY_ORTHOGONAL, AVERAGE_DIRECTION, AVERAGE_POSITION, NOT_SET};
        bool isR;
		DirectionType direction1_type, direction2_type;
		AtomInCrystalID centralAtom;
		std::vector<AtomInCrystalID> refPoint_1, refPoint_2;
		int coordinate_1, coordinate_2;

		*/
        lcs.isR = mLocalCoordinateSystem.isR;

        lcs.centralAtom = matchMap.atomMatch[0];

        lcs.coordinate_1 = mLocalCoordinateSystem.lcs_coordinate_1;
        lcs.coordinate_2 = mLocalCoordinateSystem.lcs_coordinate_2;

        lcs.refPoint_1.clear();
        lcs.refPoint_2.clear();

        if (mLabeledAtomsAxis1)
            for (int i = 0; i < mLocalCoordinateSystem.lcs_axis_1_definition.size(); i++)
                lcs.refPoint_1.push_back(matchMap.atomMatch[mLocalCoordinateSystem.lcs_axis_1_definition[i]]);
        else
            setRefPoint(describedStructure, matchMap, mLocalCoordinateSystem.lcs_axis_type_1, mLocalCoordinateSystem.lcs_axis_1_definition,
                lcs.refPoint_1);

        if (mLabeledAtomsAxis2)
            for (int i = 0; i < mLocalCoordinateSystem.lcs_axis_2_definition.size(); i++)
                lcs.refPoint_2.push_back(matchMap.atomMatch[mLocalCoordinateSystem.lcs_axis_2_definition[i]]);
        else
            //setRefPoint(describedStructure, matchMap, mLocalCoordinateSystem.lcs_axis_type_2, mLocalCoordinateSystem.lcs_axis_2_definition, lcs.refPoint_2);
            setRefPoint(describedStructure, matchMap, mLocalCoordinateSystem.lcs_axis_type_2, mLocalCoordinateSystem.lcs_axis_2_definition, 
                mLocalCoordinateSystem.lcs_axis_type_1, lcs.refPoint_1, lcs.refPoint_2);
    }

    void LocalCoordinateSystemAssigner::setRefPoint(
        const StructureWithDescriptors &describedStructure,
        const MatchMap &matchMap,
        const AtomTypeLCS::LCS_AxisType &axisType,
        const std::vector<int> &refPointCode,
        std::vector<int> &referencePointAtoms)
        const
    {
        AtomTypeLCS::LCS_AxisType fakeLcsAxisType = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
        std::vector<int> fakeReferencePoint;

        setRefPoint(describedStructure, matchMap, axisType, refPointCode, fakeLcsAxisType, fakeReferencePoint, referencePointAtoms);
    }


    void LocalCoordinateSystemAssigner::setRefPoint(
        const StructureWithDescriptors &describedStructure,
        const MatchMap &matchMap,
        const AtomTypeLCS::LCS_AxisType &axisType,
        const std::vector<int> &refPointCode,
        const AtomTypeLCS::LCS_AxisType &theOtherAxisType,
        const std::vector<int> &theOtherRefPoint,
        std::vector<int> &referencePointAtoms)
        const
    {
        map<string, int>::const_iterator iterator;
        string errorMessage; 

        referencePointAtoms.clear();

        if (axisType == AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL)
            return;

        if (axisType == AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM || axisType == AtomTypeLCS::LCS_AxisType::NOT_OF_ATOMIC_NUMBER || axisType == AtomTypeLCS::LCS_AxisType::ANY_ATOM)
        {
            vector<int> candidateAtoms;
            vector<double> candidateAtomsBondLengths;
            vector<int>::iterator alreadyUsedCandidate;
            identifyCandidateLscAtoms(describedStructure, matchMap, axisType, refPointCode, candidateAtoms, candidateAtomsBondLengths);
            if (candidateAtoms.empty())
                on_error::throwException("problem with establishing local coordinate system for multipolar representation of atomic electron density", __FILE__, __LINE__);
            // extract atom alredy used in reference point definition from candidateAtoms 
            if (theOtherAxisType == AtomTypeLCS::LCS_AxisType::ANY_ATOM || theOtherAxisType == AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM || theOtherAxisType == AtomTypeLCS::LCS_AxisType::NOT_OF_ATOMIC_NUMBER)
            {
                alreadyUsedCandidate = find(candidateAtoms.begin(), candidateAtoms.end(), theOtherRefPoint[0]);
                if (alreadyUsedCandidate != candidateAtoms.end())
                    candidateAtoms.erase(alreadyUsedCandidate);
            }

            referencePointAtoms.push_back(chooseBestAtom(describedStructure, matchMap, candidateAtoms, candidateAtomsBondLengths));
            return;
        }

        if (axisType == AtomTypeLCS::LCS_AxisType::LABELED_ATOM)
        {
            referencePointAtoms.push_back(matchMap.atomMatch[refPointCode[0]]);
            return;
        }


        if (axisType == AtomTypeLCS::LCS_AxisType::LABELED_RING)
        {
            iterator = matchMap.labeledRingsMap.find(mAtomTypeRingLabels[refPointCode[0]]);
            if (iterator == matchMap.labeledRingsMap.end())
                on_error::throwException("BUG - please report it.", __FILE__, __LINE__);
            int ringIndex = iterator->second;
            referencePointAtoms = describedStructure.planarRings[ringIndex];
            return;
        }

        if (axisType == AtomTypeLCS::LCS_AxisType::NON_LABELED_RING)
        {
            int nRings, ringIndex, structureRingIndex;
            bool ringFound = false;

            nRings = matchMap.unnamedRingsMap[0].size();
            // loopp over rings attached to the central atom
            for (ringIndex = 0; ringIndex < nRings; ringIndex++)
            {
                structureRingIndex = matchMap.unnamedRingsMap[0][ringIndex];
                // is the ring of right size
                if (describedStructure.planarRings[structureRingIndex].size() == refPointCode[0])
                {
                    // could be the ring previously used
                    if (theOtherAxisType == AtomTypeLCS::LCS_AxisType::NON_LABELED_RING || theOtherAxisType == AtomTypeLCS::LCS_AxisType::LABELED_RING)
                    {
                        // is not it already used ?
                        if (theOtherRefPoint != describedStructure.planarRings[structureRingIndex])
                            ringFound = true;
                    }
                    ringFound = true;

                    if (ringFound)
                    {
                        referencePointAtoms = describedStructure.planarRings[structureRingIndex];
                        return;
                    }
                }
            }

            // appropriate ring for making reference point not found - signalize error

            errorMessage = string("Appropriate ring for defining local coordinate system of aspherical atom not found - ") +
                string("it indicates invalid definition of atom type (atom type ID : '") + mAtomTypeID + string("'");
            on_error::throwException(errorMessage, __FILE__, __LINE__);
            return;
        }

        if (axisType == AtomTypeLCS::LCS_AxisType::ATOM_LIST)
        {
            for (int i = 0; i < refPointCode.size(); i++)
                referencePointAtoms.push_back(matchMap.atomMatch[refPointCode[i]]);
            return;
        }

        if(axisType == AtomTypeLCS::LCS_AxisType::MOST_ORTHOGONAL_ATOM)
        {
            if(theOtherAxisType != AtomTypeLCS::LCS_AxisType::LABELED_ATOM && theOtherAxisType != AtomTypeLCS::LCS_AxisType::ATOM_LIST)
                on_error::throwException("MOST_ORTHOGONAL_ATOM can be used only when the other axis is defined by LABELED_ATOM or ATOM_LIST", __FILE__, __LINE__);
            
            int centralAtom = matchMap.atomMatch[refPointCode[0]];
            auto& candidateAtoms = describedStructure.connectivity[centralAtom];

            if (candidateAtoms.empty())
                on_error::throwException("problem with establishing local coordinate system for multipolar representation of atomic electron density", __FILE__, __LINE__);
            if (theOtherRefPoint.size() == 1)
                if(candidateAtoms.size()==1)
                    if(candidateAtoms[0] == matchMap.atomMatch[theOtherRefPoint[0]])
                        on_error::throwException("problem with establishing local coordinate system for multipolar representation of atomic electron density", __FILE__, __LINE__);

            int bestAtom = -1;
            double maxCosine = -2.0; 
            Vector3d rCentralAtom = describedStructure.atomDescriptors[matchMap.atomMatch[0]].position;
            Vector3d rOtherRefPoint; 
            for(int atomIdx: theOtherRefPoint)
                rOtherRefPoint += describedStructure.atomDescriptors[atomIdx].position;
            rOtherRefPoint /= theOtherRefPoint.size();
            Vector3d vOther = rOtherRefPoint - rCentralAtom;
            vOther.normalize();
            Vector3d rCandidate, vCandidate;
            double cosine;
            for (int i = 0; i < candidateAtoms.size(); i++)
            {
                rCandidate = describedStructure.atomDescriptors[candidateAtoms[i]].position;
                vCandidate = rCandidate - rCentralAtom;
                vCandidate.normalize();
                cosine = fabs(vCandidate * vOther);
                if (cosine > maxCosine)
                {
                    maxCosine = cosine;
                    bestAtom = candidateAtoms[i];
                }
            }
            referencePointAtoms.push_back(bestAtom);
            return;
        }

		if (axisType == AtomTypeLCS::LCS_AxisType::AVERAGE_DIRECTION)
		{
			int nLabeledAtoms = refPointCode[0];
			for(int i=1;i<=nLabeledAtoms;i++)
				referencePointAtoms.push_back(matchMap.atomMatch[refPointCode[i]]);

			vector<int> candidateAtoms;
			vector<double> candidateAtomsBondLengths;

			int centralAtom = matchMap.atomMatch[0];
			const vector<int> &centralAtomNeighbours = describedStructure.connectivity[centralAtom];

			for (int i = nLabeledAtoms + 1; i < refPointCode.size(); i++)
			{
				int atomicNumber = refPointCode[i];
				candidateAtoms.clear();
				candidateAtomsBondLengths.clear();
				for (int j = 0; j < centralAtomNeighbours.size(); j++)
					if (describedStructure.atomDescriptors[centralAtomNeighbours[j]].atomicNumber == atomicNumber &&
						find(referencePointAtoms.begin(), referencePointAtoms.end(), centralAtomNeighbours[j]) == referencePointAtoms.end())
					{
						candidateAtoms.push_back(centralAtomNeighbours[j]);
						candidateAtomsBondLengths.push_back(describedStructure.bondLengths[centralAtom][j]);
					}
				if (candidateAtoms.empty())
					on_error::throwException("problem with establishing local coordinate system for multipolar representation of atomic electron density", __FILE__, __LINE__);
				referencePointAtoms.push_back(chooseBestAtom(describedStructure, matchMap, candidateAtoms, candidateAtomsBondLengths));
			}

		}

    }


    void LocalCoordinateSystemAssigner::identifyCandidateLscAtoms(
        const StructureWithDescriptors &describedStructure,
        const MatchMap &matchMap,
        const AtomTypeLCS::LCS_AxisType &axisType,
        const std::vector<int> &refPointCode,
        std::vector<int> &candidateAtoms,
        std::vector<double> &candidateAtomsBondLengths)
        const
    {
        int namedAtomIdx;
        std::vector<int> initialCandidateAtoms;
        std::vector<double> initialCandidateAtomsBondLengths;
			           
        std::set<int> initialCandidates;
        
        if (axisType == AtomTypeLCS::LCS_AxisType::ANY_ATOM)
            namedAtomIdx = matchMap.atomMatch[refPointCode[0]];
        else
            namedAtomIdx = matchMap.atomMatch[refPointCode[1]];

        auto &namedAtoms = matchMap.atomMatch;
        auto &connectedAtomsList = describedStructure.connectivity[namedAtomIdx];
        
        for (int i=0; i < connectedAtomsList.size(); i++)
            if (find(namedAtoms.begin(), namedAtoms.end(), connectedAtomsList[i]) == namedAtoms.end())
            {
                initialCandidateAtoms.push_back(connectedAtomsList[i]);
                initialCandidateAtomsBondLengths.push_back(describedStructure.bondLengths[namedAtomIdx][i]);
            }

        
        if (axisType == AtomTypeLCS::LCS_AxisType::ANY_ATOM)
        {
            candidateAtoms.swap(initialCandidateAtoms);
            candidateAtomsBondLengths.swap(initialCandidateAtomsBondLengths);
            return;
        }

        candidateAtoms.clear();
        candidateAtomsBondLengths.clear();
        if (axisType != AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM && axisType != AtomTypeLCS::LCS_AxisType::NOT_OF_ATOMIC_NUMBER)
            on_error::throwException("BUG - please report it.", __FILE__, __LINE__);

        for (int i = 0; i < initialCandidateAtoms.size(); i++)
        {
            bool fulfillsCandidateCriteria = false;
            if (axisType == AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM)
                fulfillsCandidateCriteria = (describedStructure.atomDescriptors[initialCandidateAtoms[i]].atomicNumber == refPointCode[0]);
            if (axisType == AtomTypeLCS::LCS_AxisType::NOT_OF_ATOMIC_NUMBER)
                fulfillsCandidateCriteria = (describedStructure.atomDescriptors[initialCandidateAtoms[i]].atomicNumber != refPointCode[0]);

            if (fulfillsCandidateCriteria)
            {
                candidateAtoms.push_back(initialCandidateAtoms[i]);
                candidateAtomsBondLengths.push_back(initialCandidateAtomsBondLengths[i]);
            }
        }

        if (candidateAtoms.empty())
            cout << "problem here " << __FILE__ << " " << __LINE__ << endl;

    }

    int LocalCoordinateSystemAssigner::chooseBestAtom(
        const StructureWithDescriptors &describedStructure,
        const MatchMap &matchMap,
        const std::vector<int> &candidateAtoms,
        const std::vector<double> &candidateAtomsBondLengths)
        const
    {
        int i, nCandidates = candidateAtoms.size();
        vector<pair<vector<double>, int> > descriptors(nCandidates);
        double bondLength;
        int atomIdx;
        
        for (i = 0; i < nCandidates; i++)
        {
            bondLength = candidateAtomsBondLengths[i];
            atomIdx = candidateAtoms[i];

            descriptors[i].first.resize(3);
            descriptors[i].first[0] = (double) describedStructure.atomDescriptors[atomIdx].atomicNumber;
            descriptors[i].first[1] = -double(describedStructure.connectivity[atomIdx].size());
            descriptors[i].first[2] = -bondLength;
            descriptors[i].second = i;
        }

        return candidateAtoms[std::max_element(descriptors.begin(), descriptors.end())->second];
    }

}
 