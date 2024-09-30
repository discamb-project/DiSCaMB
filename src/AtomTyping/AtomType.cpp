#include "discamb/AtomTyping/AtomType.h"
#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicChemistry/periodic_table.h"

#include <iostream>
#include <map>


using namespace std;
using namespace discamb;
using namespace string_utilities;

namespace {

	void signalizeErrorInRingSpec(const std::string &ringInfo,int line)
	{
		on_error::throwException(string("invalid specification of ring membership in atom type: '") + ringInfo + string("'") , __FILE__ , line );
	}

    
    // set 3 or 4 member rings, planar rings of ringInfo should be already set (for checking consistency)
    void setN_RingInfo(
        const std::string info,
        Tribool &inRing,
        int &nRings,
        int ringSize,
        const AtomRingInfo &ringInfo)
    {
        map<char, Tribool> charToRingInfo{ { '*', Tribool::Undefined},
                                           { '+', Tribool::True },
                                           { '-', Tribool::False }};
        if (isdigit(info[0]))
        {
            inRing = Tribool::True;
            nRings = stoi(info);
        }
        else
        {
            if (info.size() > 1 || charToRingInfo.find(info[0]) == charToRingInfo.end())
                on_error::throwException(string("invalid format of 3 or 4 member ring membershio specification :'") + info + string("'"), __FILE__, __LINE__);
            inRing = charToRingInfo[info[0]];
            nRings = -1;
            if (inRing == Tribool::False)
            {
                bool error = false;
                for (auto ring : ringInfo.labeledContainingRings)
                    if (ring.first == ringSize)
                        error = true;
                for (auto ring : ringInfo.nonLabeledContainingRings)
                    if (ring == ringSize)
                        error = true;
                if (error)
                {
                    string error_message = to_string(ringSize) + string("-member ring membership specification '") +
                        info + string("' contradicts planar rings specification");
                    on_error::throwException(error_message, __FILE__, __LINE__);
                }
            }

        }
        
    }
}

namespace discamb {

    void AtomRingInfo::set(
        const std::string &ringInfo, 
        const std::string &ring3info,
        const std::string &ring4info)
    {
        vector<string> words;
        int i, j, n, nChar, ringSize;
        string ringSizeAsString;
        bool notContainingRing;
        bool labelledRing;

        //#############################
        //
        // processing info on planar rings
        //
        //#############################

        labeledContainingRings.clear();
        labeledNonContainingRings.clear();
        nonLabeledContainingRings.clear();
        nonLabeledNonContainingRings.clear();

        split(ringInfo, words, ',');

        this->inRing = Tribool::Undefined;

        if (find(words.begin(), words.end(), string("+")) != words.end() ||
            find(words.begin(), words.end(), string("Ar+")) != words.end())
        {
            if (words.size() > 1)
                signalizeErrorInRingSpec(ringInfo, __LINE__);
            this->inRing = Tribool::True;
            //return;
        }
        else
        {

            if (find(words.begin(), words.end(), string("-")) != words.end()||
                find(words.begin(), words.end(), string("Ar-")) != words.end())
            {
                if (words.size() > 1)
                    signalizeErrorInRingSpec(ringInfo, __LINE__);
                this->inRing = Tribool::False;
                //return;
            }
            else
            {
                // RING + and RING - cases already processed
                // now the remaining possibilities of type *,6,!6,6A,!6A

                this->inAnyAdditionalRing = false;

                n = words.size();

                for (i = 0; i < n; i++)
                {
                    notContainingRing = false;
                    labelledRing = false;

                    if (words[i][0] == '!')
                    {
                        notContainingRing = true;
                        words[i] = words[i].substr(1);
                    }

                    if (words[i][0] == '*')
                    {
                        // verify correctness
                        // more than one star e.g. 6A,*,*  or  label contain star + something 
                        if (this->inAnyAdditionalRing || words[i].size() > 1)
                            signalizeErrorInRingSpec(ringInfo, __LINE__);

                        this->inAnyAdditionalRing = true;
                        continue;
                    }

                    // extracts ring size from ring label e.g. 6 from 6A

                    if (isdigit(words[i][0]))
                    {
                        nChar = words[i].size();
                        ringSizeAsString.clear();

                        for (j = 0; j < nChar; j++)
                            if (isdigit(words[i][j]))
                                ringSizeAsString += words[i][j];
                            else
                            {
                                labelledRing = true;
                                break;
                            }
                    }
                    else
                        signalizeErrorInRingSpec(ringInfo, __LINE__);

                    ringSize = convertFromString <int>(ringSizeAsString);

                    // save info on ring

                    if (labelledRing)
                    {
                        if (notContainingRing)
                            labeledNonContainingRings.push_back(make_pair(ringSize, words[i]));
                        else
                            labeledContainingRings.push_back(make_pair(ringSize, words[i]));
                    }
                    else
                    {
                        if (notContainingRing)
                            nonLabeledNonContainingRings.push_back(ringSize);
                        else
                            nonLabeledContainingRings.push_back(ringSize);
                    }

                }

                if (!labeledContainingRings.empty() || !nonLabeledContainingRings.empty())
                    this->inRing = Tribool::True;
            }
        }
        //#########################################
        //
        // processing info on 3 and 4 member rings
        //
        //#########################################

        setN_RingInfo(ring3info, in3Ring, n3rings, 3, *this);
        setN_RingInfo(ring4info, in4Ring, n4rings, 4, *this);
    }


    void AtomType::setLocalCoordinateSystem(
        const std::string &lcsAsString)
    {
        //e.g. X r(5A) Y C2 R

        vector<string> lcsWords;
        split(lcsAsString, lcsWords, CharacterType::WHITE_SPACE);
        if (lcsWords.size() != 5)
            on_error::throwException(string("invalid definition of local coordinate system, expected 5 'words'separated with space, instead got ") +
                string_utilities::convertToString(lcsWords.size()) + string(" words"), __FILE__, __LINE__);

        //localCoordinateSystem.lcs_axis_1_definition
        localCoordinateSystem.lcs_coordinate_1 = getCoordinate(lcsWords[0]);
        localCoordinateSystem.lcs_coordinate_2 = getCoordinate(lcsWords[2]);

        if (lcsWords[4] == string("R"))
            localCoordinateSystem.isR = true;
        else
        {
            if (lcsWords[4] == string("L"))
                localCoordinateSystem.isR = false;
            else
                on_error::throwException(string("invalid specification of chirality in local coordinate system definition: '") + lcsAsString + string("'"), __FILE__, __LINE__);
        }

        setReferencePoint(lcsWords[1], localCoordinateSystem.lcs_axis_1_definition, localCoordinateSystem.lcs_axis_type_1);
        setReferencePoint(lcsWords[3], localCoordinateSystem.lcs_axis_2_definition, localCoordinateSystem.lcs_axis_type_2);

        if (localCoordinateSystem.lcs_axis_type_1 == AtomTypeLCS::LCS_AxisType::LABELED_ATOM)
            if (localCoordinateSystem.lcs_axis_1_definition.size() == 1)
                if (localCoordinateSystem.lcs_axis_1_definition[0] == 0)
                    on_error::throwException(string("local coordinate system definition with vector from central atom to itself for '") + this->id + string("'"), __FILE__, __LINE__);
        if (localCoordinateSystem.lcs_axis_type_2 == AtomTypeLCS::LCS_AxisType::LABELED_ATOM)
            if (localCoordinateSystem.lcs_axis_2_definition.size() == 1)
                if (localCoordinateSystem.lcs_axis_2_definition[0] == 0)
                    on_error::throwException(string("local coordinate system definition with vector from central atom to itself for '") + this->id + string("'"), __FILE__, __LINE__);

    }

    int AtomType::getCoordinate(
        const std::string &coordinate)
        const
    {
        if (coordinate == string("X"))
            return 0;
        if (coordinate == string("Y"))
            return 1;
        if (coordinate == string("Z"))
            return 2;

        on_error::throwException(string("invalid coordinate string in local coordinate definition: '") + coordinate + string("'"), __FILE__, __LINE__);

        return 0;
    }

    void AtomType::setReferencePoint(
        const std::string &refPoint,
        std::vector<int> &definition,
        AtomTypeLCS::LCS_AxisType &definitionType)
        const
    {
        definition.clear();
        definitionType = AtomTypeLCS::LCS_AxisType::NOT_SET;

        bool ringOrBondPoint = false;
        int stringSize;
        optional<int> idx;

        stringSize = refPoint.size();

		bool hasBracket = (refPoint.find('(') != string::npos && refPoint.back() == ')');
		string preBracket, bracketsContent;
		
		if (hasBracket)
		{
			preBracket = refPoint.substr(0, refPoint.find('('));
			bracketsContent = refPoint.substr(refPoint.find('(')+1, refPoint.size()- refPoint.find('(') - 2);
			// ala(ma)
			// refPoint.size()- refPoint.find('(') - 2 = 7 - 3- 2 = 2

		}

        if (stringSize > 2)
            if (islower(refPoint[0]) && refPoint[1] == '(')
                ringOrBondPoint = true;

		if (preBracket == string("r"))
		{

			if (bracketsContent.empty())
				on_error::throwException("invalid specification of local coordinate system - wrong ring definition (empty bracket)", __FILE__, __LINE__);

			if (isdigit(bracketsContent.back()))
			{
				definitionType = AtomTypeLCS::LCS_AxisType::NON_LABELED_RING;
				definition.push_back(atoi(bracketsContent.c_str()));
			}
			else
			{
				definitionType = AtomTypeLCS::LCS_AxisType::LABELED_RING;
				vector<string>::const_iterator it = find(ringLabels.begin(), ringLabels.end(), bracketsContent);
				if (it == ringLabels.end())
					on_error::throwException(string("invalid specification of local coordinate system in definition of atom type: '") + id + string("' - referencing to undefined ring "), __FILE__, __LINE__);

				definition.push_back(std::distance(ringLabels.begin(), it));
			}
			return;
		}

		if (preBracket == "average_direction")
		{
			vector<string> atomsAsString;
			int atomIndex, nAtoms;
			vector<int> labeledAtoms, nonLabeledAtoms;

			string_utilities::split(bracketsContent, atomsAsString, ',');

			nAtoms = atomsAsString.size();

            

			for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
			{
                
                if (isdigit(atomsAsString[atomIndex][atomsAsString[atomIndex].size() - 1]))
                {
                    idx = getAtomIndex(atomsAsString[atomIndex]);
                    if (idx)
                        labeledAtoms.push_back(*idx);
                    else
                        on_error::throwException("problem when defining atom type", __FILE__, __LINE__);
                }
				else
					nonLabeledAtoms.push_back(periodic_table::atomicNumber(atomsAsString[atomIndex]));
			}

			definitionType = AtomTypeLCS::LCS_AxisType::AVERAGE_DIRECTION;
			definition.push_back(labeledAtoms.size());
			definition.insert(definition.end(), labeledAtoms.begin(), labeledAtoms.end());
			definition.insert(definition.end(), nonLabeledAtoms.begin(), nonLabeledAtoms.end());

			return;
		}

        // the case when reference point is an single atom 


        if (preBracket == string("any_atom"))
        {
            definitionType = AtomTypeLCS::LCS_AxisType::ANY_ATOM;
        }


        if (refPoint[0] == '!' && definitionType == AtomTypeLCS::LCS_AxisType::NOT_SET)
        {
            definition.push_back(periodic_table::atomicNumber(preBracket.substr(1)));
            definitionType = AtomTypeLCS::LCS_AxisType::NOT_OF_ATOMIC_NUMBER;
        }

        if (isdigit(refPoint.back()) && definitionType == AtomTypeLCS::LCS_AxisType::NOT_SET) // labeled atom
        {
            idx = getAtomIndex(refPoint);
            if(idx)
                definition.push_back(*idx);
            else
                on_error::throwException("problem when defining atom type", __FILE__, __LINE__);
            definitionType = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
            return;
        }

        if (refPoint == string("any_orthogonal"))
        {
            definitionType = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
            return;
        }

        if(definitionType == AtomTypeLCS::LCS_AxisType::NOT_SET)
        {
            definition.push_back( periodic_table::atomicNumber(preBracket));
            definitionType = AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM;
        }
        // in the case of ANY_ATOM, NOT_OF_ATOMIC_NUMBER, NON_LABELED_ATOM there must be specified an atom
        // to which the reference atom is bond 
        if (definitionType == AtomTypeLCS::LCS_AxisType::ANY_ATOM || definitionType == AtomTypeLCS::LCS_AxisType::NOT_OF_ATOMIC_NUMBER ||
            definitionType == AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM)
        {
			if (!hasBracket || bracketsContent.empty())
				on_error::throwException("invalid specification of local coordinate system for type " +
					this->id, __FILE__, __LINE__);
            idx = getAtomIndex(bracketsContent);
            if(idx)
                definition.push_back(*idx);
            else
                on_error::throwException("problem when defining atom type", __FILE__, __LINE__);
        }

    }


    std::optional<int> AtomType::getAtomIndex(
        const std::string &label)
        const
    {
        int i, n = atoms.size();
        for (i = 0; i < n; i++)
            if (atoms[i].label == label)
                return i;

        on_error::throwException(string("invalid atom label : '") + label + string("' used in local coordinate system definition"), __FILE__, __LINE__);

        return -1;
    }

    void AtomType::sphericalTypes(
        std::vector<AtomType> &types)
    {
        types.clear();
        AtomType sphericalType;
        sphericalType.atoms.resize(1);
        auto &centralAtom = sphericalType.atoms.back();
        //centralAtom.neighborsAtomicNumbersUniquelyDefined = false;
        centralAtom.planar = Tribool::Undefined;
        centralAtom.ringInfo.set("*", "*", "*");
        sphericalType.commentLines = { "spherical type" };      
        string date = "Thu May  2 12:30:32 2019"; 
        sphericalType.definitionModificationDate = date;
        sphericalType.entryCreationDate = date;
        sphericalType.localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
        sphericalType.localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
        sphericalType.parameterModificationDate = date;
        sphericalType.symmetry = "sph";

        string symbol;
        for (int atomicNumber = 1; atomicNumber <= 36; atomicNumber++)
        {
            sphericalType.atoms[0].atomic_number_range.insert(atomicNumber);
            symbol = periodic_table::symbol(atomicNumber);
            sphericalType.atoms[0].label = symbol;
            sphericalType.commentLines = { string("    spherical ") + symbol };
            sphericalType.id = symbol + string("000");
            types.push_back(sphericalType);
        }
    }

    void AtomType::transformToCanonicalForm()
    {
        vector<vector<int> > shells;
        graph_algorithms::breadth_first_search(connectivity, 0, shells);

        vector<int> old2newAtomOrder;
    }

    // this and atomType have to be transformed to canonical form first
    //bool AtomType::equivalent(
    //    const AtomType& atomType)
    //    const
    //{
    //    return true;
    //}

    //void AtomType::transformToCanonicalForm()
    //{

    //}

    void AtomType::unnamedNeighboursAtomicNumbers(
        int atomIdx,
        std::vector<int>& unnamedNeighborsZ)
    {
        unnamedNeighborsZ.clear();
        unnamedNeighborsZ.insert(unnamedNeighborsZ.end(),atoms[atomIdx].neighborsAtomicNumbers.begin(), atoms[atomIdx].neighborsAtomicNumbers.end());
        for (int idx : connectivity[atomIdx])
            if (atoms[idx].atomic_number_range.size() == 1)
                unnamedNeighborsZ.erase(find(unnamedNeighborsZ.begin(), unnamedNeighborsZ.end(), *atoms[idx].atomic_number_range.begin()));
    }

    bool AtomType::setDefaultLocalCoordinateSystem()
    {
        localCoordinateSystem = AtomTypeLCS();
        if (connectivity.empty())
            return false;
        if (connectivity[0].empty())
        {
            if (atoms[0].nNeighbours == 0)
            {
                localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
                localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
                return true;
            }
            if (atoms[0].nNeighbours == 1)
            {
                localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::ANY_ATOM;
                localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
                localCoordinateSystem.lcs_axis_1_definition.push_back(0);
                return true;
            }
            if (atoms[0].nNeighbours > 1)
            {
                localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::ANY_ATOM;
                localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::ANY_ATOM;
                localCoordinateSystem.lcs_axis_1_definition.push_back(0);
                localCoordinateSystem.lcs_axis_2_definition.push_back(0);
                return true;
            }
            
            return false;
        }

        if (connectivity[0].size() == 1)
        {
            localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
            localCoordinateSystem.lcs_axis_1_definition.push_back(connectivity[0][0]);
            
            vector<int> unnamedNeighboursZ;
            unnamedNeighboursAtomicNumbers(0, unnamedNeighboursZ);

            if (unnamedNeighboursZ.size() > 0)
            {
                localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::NON_LABELED_ATOM;
                localCoordinateSystem.lcs_axis_2_definition.push_back(
                    *max_element(unnamedNeighboursZ.begin(), unnamedNeighboursZ.end()));
                localCoordinateSystem.lcs_axis_2_definition.push_back(0);
            } 
            else if (atoms[0].neighborsAtomicNumberRanges.size() > 0)
            {
                // expected manual definition with reference to the first neighbour which atomic number is given as a range
                return false;
            }
            else
            {
                int neighbourIdx = connectivity[0][0];
                int neighbourNeighbourIdx;
                if (connectivity[neighbourIdx].size() > 1)
                {
                    neighbourNeighbourIdx = (connectivity[neighbourIdx][0] == 0 ? connectivity[neighbourIdx][1] : connectivity[neighbourIdx][0]);
                    localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
                    localCoordinateSystem.lcs_axis_2_definition.push_back(neighbourNeighbourIdx);
                }
                else if (atoms[neighbourIdx].nNeighbours > 1)
                {
                    localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::ANY_ATOM;
                    localCoordinateSystem.lcs_axis_2_definition.push_back(neighbourIdx);
                }
                else
                    localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::ANY_ORTHOGONAL;
                
                return true;
            }

            return true;
        }

        if (connectivity[0].size() == 2)
        {
            localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
            localCoordinateSystem.lcs_axis_1_definition.push_back(connectivity[0][0]);

            localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
            localCoordinateSystem.lcs_axis_2_definition.push_back(connectivity[0][1]);

            return true;
        }

        if (connectivity[0].size() > 2)
        {
            localCoordinateSystem.lcs_axis_type_1 = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
            

            localCoordinateSystem.lcs_axis_type_2 = AtomTypeLCS::LCS_AxisType::LABELED_ATOM;
            
            int neighbour1, neighbour2;
            bool neighboursChosen = chooseNamedNeighboursForLcs(neighbour1, neighbour2);

            localCoordinateSystem.lcs_axis_1_definition.push_back(neighbour1);
            localCoordinateSystem.lcs_axis_2_definition.push_back(neighbour2);
            
            return neighboursChosen;
        }


        return false;
    }


    bool AtomType::chooseNamedNeighboursForLcs(
        int& neighbour1,
        int& neighbour2)
    {
        if (connectivity.size() < 3)
            return false;
        
        int neighbourIdx, nNeighbours = connectivity[0].size();

        if (nNeighbours < 2)
            return false;
        //[idx], second index:
        // 0    lowest atomic number
        // 1    -n neighbours
        // 2    1 if it can have more neighbours, 0 otherwise
        // 3    sum of atomic numbers of neighbours
        // 4    neighbour index
        vector<vector<int> > neigborsScoreAndIdx(nNeighbours,vector<int>(5));
        for(int i=0; i< nNeighbours; i++)
        {
            neighbourIdx = connectivity[0][i];
            neigborsScoreAndIdx[i][0] = *atoms[neighbourIdx].atomic_number_range.begin();
            neigborsScoreAndIdx[i][1] = -atoms[neighbourIdx].nNeighbours;
            neigborsScoreAndIdx[i][2] = (atoms[neighbourIdx].fixedNumberOfNeighbors ? 0 : 1);
            neigborsScoreAndIdx[i][3] = 0;
            for (int z : atoms[neighbourIdx].neighborsAtomicNumbers)
                neigborsScoreAndIdx[i][3] += z;
            for (const auto &range: atoms[neighbourIdx].neighborsAtomicNumberRanges)
                neigborsScoreAndIdx[i][3] += *range.begin();
            neigborsScoreAndIdx[i][4] = neighbourIdx;
        }
        sort(neigborsScoreAndIdx.begin(), neigborsScoreAndIdx.end());
        neighbour1 = neigborsScoreAndIdx[nNeighbours - 1][4];
        neighbour2 = neigborsScoreAndIdx[nNeighbours - 2][4];
        return true;
    }
}

