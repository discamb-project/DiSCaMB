#include "discamb/AtomTyping/MacromolecularAtomAssignment.h"
#include "discamb/BasicUtilities/string_utilities.h"

using namespace std;

namespace discamb {

    namespace {
        void on_wrong_format(const string& s)
        {
            string message = "Wrong format of macromolecular atom assignment definition: " + s + "\n";
            on_error::throwException(message, __FILE__, __LINE__);
        }

        // returns -1 if the input string is not a valid coordinate name
        // othrwise returns the index of the coordinate (0 for x, 1 for y, 2 for z)
        
        int coordinateIdx(const string& s)
        {
            if (s.size() != 1)
                on_wrong_format(s);
            char c = s[0];
            c = tolower(c);
            if (c != 'x' && c != 'y' && c != 'z')
                return -1;
            map<char, int> char2xyzInt{ {'x', 0}, {'y', 1}, {'z', 2} };
            return char2xyzInt[c];
        }

        // splits the input string by comma, and for each atom definition, 
        // extracts the atom name and residue offset (if specified)
        // returns false if encounters problem when pro essing the input string

        bool stringToAtomList(
            const string& s,
            vector<pair<string, int> >& atoms)
        {
            vector<string> atomDefs;
            string_utilities::split(s, atomDefs, ',');
            for (const auto& atomDef : atomDefs)
            {
                int residueOffset = 0;
                if (atomDef.back() == ']')
                {
                    size_t pos = atomDef.find('[');
                    if (pos == string::npos)
                        return false;
                    string offsetStr = atomDef.substr(pos + 1, atomDef.size() - pos - 2);
                    try
                    {
                        residueOffset = stoi(offsetStr);
                    }
                    catch (const exception&)
                    {
                        return false;
                    }
                }
                atoms.push_back({ atomDef, residueOffset });
            }
            return true;
        }
        // returns type of direction specified by the input string,
        // and fills the atoms vector with the parsed atom definitions 
        // (label and residue offset) if applicable
        // return LcsDirectionType::NOT_SET if the input string does not 
        // specify a valid direction type
        LcsDirectionType get_direction(
            const string& s,
            vector<pair<string, int> >& atoms)
        {
            atoms.clear();

            if (s == "any_orthogonal")
                return LcsDirectionType::ANY_ORTHOGONAL;

            if (s.find("average_direction:") == 0)
            {
                string atomListStr = s.substr(string("average_direction:").size());
                if(stringToAtomList(atomListStr, atoms))
                    return LcsDirectionType::AVERAGE_DIRECTION;
            }
            if (s.find("average_position:") == 0)
            {
                string atomListStr = s.substr(string("average_position:").size());
                if (stringToAtomList(atomListStr, atoms))
                    return LcsDirectionType::AVERAGE_POSITION;
            }
            if (stringToAtomList(s, atoms))
                return LcsDirectionType::AVERAGE_POSITION;

            return LcsDirectionType::NOT_SET;
        }
    }

    void MacromolecularAtomAssignment::set(
        const std::string s) 
    {
        vector<string> words;
        string_utilities::split(s, words);
        /*
        "N401a Z CA X H"
        "N401a Z CA X H R"
        "C414 X N Y C N C CB"
        "C414 X N Y C R N C CB"
        */
        int nWords = words.size();

        if (nWords < 5 || nWords > 9)
            on_wrong_format(s);

        type_name = words[0];
        if(words[1].size()!=1 || words[3].size() != 1)
            on_wrong_format(s);
        char c = words[1][0];
        c = tolower(c);
        if (c != 'x' && c != 'y' && c != 'z')
            on_wrong_format(s);
        map<char, int> char2xyzInt{ {'x', 0}, {'y', 1}, {'z', 2} };
        lcs.coordinate_1 = coordinateIdx(words[1]);
        lcs.coordinate_2 = coordinateIdx(words[3]);
        if (lcs.coordinate_1 < 0 || lcs.coordinate_2 < 0)
            on_wrong_format(s);
        lcs.direction1_type = get_direction(words[2], lcs.refPoint_1);
        lcs.direction2_type = get_direction(words[4], lcs.refPoint_2);
        if (lcs.direction1_type == LcsDirectionType::NOT_SET || lcs.direction2_type == LcsDirectionType::NOT_SET)
            on_wrong_format(s);

    }
}
