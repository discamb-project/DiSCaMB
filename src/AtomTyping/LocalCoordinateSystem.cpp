#include "discamb/AtomTyping/LocalCoordinateSystem.h"

#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/BasicUtilities/OnError.h"

using namespace std;

namespace {
    string atomInCrystalIdAsString(
        const discamb::AtomInCrystalID &id,
        const vector<string> &labels)
    {
        discamb::Matrix3i m, identity(1,0,0,
                                      0,1,0,
                                      0,0,1);
        discamb::Vector3<discamb::CrystallographicRational> t;
        id.getSymmetryOperation().get(m, t);
        
        string result = labels[id.index()];

        if (m == identity && t == discamb::Vector3<discamb::CrystallographicRational>(0, 0, 0))
            return result;

        string symmOpAsString;
        id.getSymmetryOperation().get(symmOpAsString);
        result += string("(") + symmOpAsString + string(")");
        return result;
    }


    void splitAtomDefinition(
        const string& atomDef,
        string& atomLabel,
        string& symmOp)
    {
        if (atomDef.find(',') == string::npos)
        {
            atomLabel = atomDef;
            symmOp = string("X,Y,Z");
        }
        else
        {
            atomLabel = atomDef.substr(0, atomDef.find(','));
            symmOp = atomDef.substr(atomDef.find(',') + 1);
        }
    }


}

namespace discamb {

    void convertUbdbLcs(
        const  LocalCoordinateSystem<int> &lcsMolecule,
        const std::vector<AtomInCrystalID> &atomMap,
        LocalCoordinateSystem<AtomInCrystalID> &lcsCrystal)
    {
        lcsCrystal.centralAtom.set(lcsMolecule.centralAtom);
        lcsCrystal.coordinate_1 = lcsMolecule.coordinate_1;
        lcsCrystal.coordinate_2 = lcsMolecule.coordinate_2;
        lcsCrystal.direction1_type = lcsMolecule.direction1_type;
        lcsCrystal.direction2_type = lcsMolecule.direction2_type;
        lcsCrystal.isR = lcsMolecule.isR;

        lcsCrystal.refPoint_1.clear();
        lcsCrystal.refPoint_2.clear();
        lcsCrystal.chirality.clear();

        lcsCrystal.centralAtom = lcsMolecule.centralAtom;

        for (auto atom : lcsMolecule.refPoint_1)
            lcsCrystal.refPoint_1.push_back(atomMap[atom]);
        for (auto atom : lcsMolecule.refPoint_2)
            lcsCrystal.refPoint_2.push_back(atomMap[atom]);
        for (auto atom : lcsMolecule.chirality)
            lcsCrystal.chirality.push_back(atomMap[atom]);

    }

    std::string ubdbLcsAsString(
        const  LocalCoordinateSystem<AtomInCrystalID> &lcs,
        const std::vector<std::string> &labels)
    {
        char xyz[] = { 'X', 'Y', 'Z' };
        string result, space(" ");
        result = xyz[lcs.coordinate_1] + space + ubdbLcsDirectionAsString(lcs.direction1_type, lcs.refPoint_1, labels) + space +
                 xyz[lcs.coordinate_2] + space + ubdbLcsDirectionAsString(lcs.direction2_type, lcs.refPoint_2, labels);
        return result;
    }

    std::string ubdbLcsAsString(
        const  LocalCoordinateSystem<int> &lcs,
        const std::vector<std::string> &labels)
    {
        char xyz[] = { 'X', 'Y', 'Z' };
        string result, space(" ");
        result = xyz[lcs.coordinate_1] + space + ubdbLcsDirectionAsString(lcs.direction1_type, lcs.refPoint_1, labels) + space +
            xyz[lcs.coordinate_2] + space + ubdbLcsDirectionAsString(lcs.direction2_type, lcs.refPoint_2, labels);
        return result;
    }

    std::string ubdbLcsDirectionAsString(
        LcsDirectionType type,
        const std::vector<int> &indices,
        const std::vector<string> &labels)
    {
        string s;
        if (indices.size() == 1)
            return labels[indices[0]];
        else
        {
            if (indices.empty())
            {
                if (type == LcsDirectionType::NOT_SET)
                    return "undefined";
                else
                    return "any_orthogonal";
            }
            else
            {
                type == LcsDirectionType::AVERAGE_POSITION ? s = "average_position(" : s = "average_direction(";
                for (int i = 0; i < indices.size(); i++)
                {
                    if (i > 0)
                        s += ",";
                    s += labels[indices[i]];
                }
                s += ")";
            }
        }
        return s;

    }

    std::string ubdbLcsDirectionAsString(
        LcsDirectionType type,
        const std::vector<AtomInCrystalID> &indices,
        const std::vector<string> &labels)
    {
        string s, symmOpAsString;
        if (indices.size() == 1)
            return atomInCrystalIdAsString(indices[0], labels);
        else
        {
            if (indices.empty())
            {
                if (type == LcsDirectionType::NOT_SET)
                    return "undefined";
                else
                    return "any_orthogonal";
            }
            else
            {
                type == LcsDirectionType::AVERAGE_POSITION ? s = "average_position(" : s = "average_direction(";
                for (int i = 0; i < indices.size(); i++)
                {
                    if (i > 0)
                        s += ",";
                    s += atomInCrystalIdAsString(indices[i],labels);
                }
                s += ")";
            }
        }
        return s;

    }

    void xdTypeLcs(
        const std::string& definition,
        const Crystal& crystal,
        LocalCoordinateSystem<AtomInCrystalID>& lcs)
    {
        /*    void splitAtomDefinition(
        const string& atomDef,
        string& atomLabel,
        string& symmOp)
*/
        lcs.chirality.clear();
        lcs.refPoint_1.clear();
        lcs.refPoint_2.clear();
        lcs.direction1_type = LcsDirectionType::AVERAGE_POSITION;
        lcs.direction2_type = LcsDirectionType::AVERAGE_POSITION;

        vector<string> words;
        string_utilities::split(definition, words);
        // central_atom X atom1 Y atom2 (optional R/L)
        if (words.size() != 5 && words.size() != 6)
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("'"), __FILE__, __LINE__);
        
        map<string, int> xyz{ {"X", 0}, {"Y", 1}, {"Z", 2} };

        if(xyz.find(words[1])==xyz.end() || xyz.find(words[3]) == xyz.end())
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("'"), __FILE__, __LINE__);

        if(words[1]== words[3])
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("',") + words[1] + string(" defined twice"), __FILE__, __LINE__);
        if(words[0]== words[2] || words[0] == words[4])
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("', central atom and atom defining direction is the same"), __FILE__, __LINE__);
        if (words[4] == words[2])
            on_error::throwException(string("invalid definition of local coordinate system '") + definition + string("', both atoms defining directions are the same"), __FILE__, __LINE__);


        lcs.coordinate_1 = xyz[words[1]];
        lcs.coordinate_2 = xyz[words[3]];

        string atomStr, symmOpStr;

        splitAtomDefinition(words[0], atomStr, symmOpStr);
        lcs.centralAtom = AtomInCrystalID(atomStr, crystal, SpaceGroupOperation(symmOpStr));
        splitAtomDefinition(words[2], atomStr, symmOpStr);
        lcs.refPoint_1.push_back(AtomInCrystalID(atomStr, crystal, SpaceGroupOperation(symmOpStr)));
        splitAtomDefinition(words[4], atomStr, symmOpStr);
        lcs.refPoint_2.push_back(AtomInCrystalID(atomStr, crystal, SpaceGroupOperation(symmOpStr)));

        lcs.isR = true;
        if(words.size()==6)
            if(words[5]==string("L"))
                lcs.isR = false;
    }

	// central_atom coordinate_1 direction_1 refpoint_2 coordinate_2 chirality 
	// C(1) C(2),C(3),C(4)[-x,-y,z] X any_orthogonal Y R 
	// C(1) average_direction:C(2)[1-x,y,z],C(3) Z C(4) Y 

	/*void LocalCoordinateSystem::set(
		const std::string &s,
		const Crystal &c)
	{
		
	}*/
}

