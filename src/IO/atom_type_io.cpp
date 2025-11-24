#include "discamb/IO/atom_type_io.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"

#include "json.hpp"

#include <fstream>
#include <algorithm>

using namespace std;
using namespace nlohmann;

/*

["          S3-X     ",
 "          |        ",
 "          C1{6A}   ",
 "         / \       ",
 "  X-{6A}N   C2{6A} "]

         C(2)
         |
    C(3)-N{nonplanar}
         |
         C(2)

     |                |
    -C3-N1(H){planar}-C2{planar}
     |                 \

*/

namespace {
    struct JsonDuplicateCheck {
        std::set<std::string> keyNames;
        std::string lastKey;
        int types_depth = -1;
        bool parsingTypes = false;
        bool operator()(int depth, nlohmann::json::parse_event_t event, nlohmann::json& parsed)
        {
            //lastKey
            if (event == nlohmann::json::parse_event_t::key)
            {
                lastKey = parsed.get<string>();
                if (parsingTypes && types_depth == depth)
                {
                    if (keyNames.count(lastKey))
                            discamb::on_error::throwException("repeating type name '" + lastKey + "'", __FILE__, __LINE__);
                    else
                        if(lastKey!="comment")
                            keyNames.insert(lastKey);
                }
                if (lastKey == "types")
                {
                    types_depth = depth+1;
                    parsingTypes = true;
                }
            }
            if (parsingTypes)
            {
                if (event == nlohmann::json::parse_event_t::object_end)
                {
                    if (depth == types_depth + 1)
                        parsingTypes = false;
                }
            }
            return true;
        }
    };

    struct AtomFromString {
        void set(const std::string& s,int endPosition);
        int labelPosition; // index in line of the first character of atom label (e.g. C in C2{planar})
        int labelStart, labelEnd; // indices of the first and last character in atom description
                                  // e.g. for (X3)C2{planar} it is position of ( and })
        // in C1(=O2) for O2 inBracket is true
        bool inBracket = false;
        int line;
        string name; // e.g. C1
        string type; // chemical symbol
        vector<string> properties; // in {}
        vector<string> additional_bonded_atoms; // in()
    };

    

    void AtomFromString::set(
        const std::string& _s,
        int endPosition)
    {
        *this = AtomFromString();
        if (_s.empty())
            return;
        string s = _s;
        size_t start_bracket, end_bracket, start_curly_bracket, end_curly_bracket;
        start_bracket = s.find('(');
        end_bracket = s.find(')');
        start_curly_bracket = s.find('{');
        end_curly_bracket = s.find('}');

        if (start_bracket != string::npos || end_bracket != string::npos)
        {
            if (start_bracket != string::npos && end_bracket != string::npos && start_bracket< end_bracket)
            {
                string neighborsStr = s.substr(start_bracket+1, end_bracket - start_bracket - 1);
                string neighborsStrNoBonds;
                for (char c : neighborsStr)
                    if (c != '=' && c != '-' && c != '#')
                        neighborsStrNoBonds += c;
                neighborsStr = neighborsStrNoBonds;
                discamb::string_utilities::split(neighborsStr, additional_bonded_atoms, ',');
                size_t count = end_bracket - start_bracket + 1;
                s.replace(start_bracket, count, count, ' ');
            }
            else
                discamb::on_error::throwException("error when parsing atom type - no matching ( or ) in " + _s, __FILE__, __LINE__);
        }

        if (start_curly_bracket != string::npos || end_curly_bracket != string::npos)
        {
            if (start_curly_bracket != string::npos && end_curly_bracket != string::npos && start_curly_bracket < end_curly_bracket)
            {
                string props = s.substr(start_curly_bracket+1, end_curly_bracket - start_curly_bracket - 1);
                discamb::string_utilities::split(props, properties, ',');
                size_t count = end_curly_bracket - start_curly_bracket + 1;
                s.replace(start_curly_bracket, count, count, ' ');
            }
            else
                discamb::on_error::throwException("error when parsing atom type - no matching { or } in " + _s, __FILE__, __LINE__);
        }
        
        int startPosition = (endPosition + 1) - s.size();
        labelPosition = s.find_first_not_of(' ') + startPosition;
        name = discamb::string_utilities::trim(s);
        
        for (auto& c : name)
            if (!isdigit(c))
                type += c;
        
        labelStart = startPosition;
        labelEnd = endPosition;
    }

    struct BondFromString {
        int position;
        int line;
        int length=1;
        char bondCharacter;
    };
    
    // joins neighboring bonds (-,=,#, and | but not / and \) 
    // into one of multiple length
    // e.g. two bonds given by '-' (i.e. --) are transformed into
    // one of length 2, this funcitonality is needed 
    // for processing things like this:
    //   
    //   R--C5----C4-R
    //      |     |
    //    R-N2-C1-N3-R
    //         |
    //         C
    //

    

    void concatenateBonds(
        vector<BondFromString>& bonds)
    {
        int nRows = 0;
        int nColumns = 0;
        vector<BondFromString> bondsUpdated;

        for (auto& bond : bonds)
        {
            if (bond.line > nRows)
                nRows = bond.line;
            if (bond.position > nColumns)
                nColumns = bond.position;
        }

        nRows++;
        nColumns++;

        vector<vector<int> > bondIdx(nRows,vector<int>(nColumns,-1));
        int nBonds = bonds.size();
        vector<bool> toInspect(nBonds, true);

        for (int i = 0; i < nBonds; i++)
            bondIdx[bonds[i].line][bonds[i].position] = i;

        vector<vector<int> > bondGroups;

        //cout << "#################\n\n";
        //for (int row = 0; row < nRows; row++)
        //{
        //    for (int col = 0; col < nColumns; col++)
        //        if (bondIdx[row][col] >= 0)
        //            cout << bonds[bondIdx[row][col]].bondCharacter;
        //        else
        //            cout << " ";
        //    cout << "\n";
        //}

        // concatenate --- == ##
        for(int row=0; row<nRows; row++)
            for (int col = 0; col < nColumns; col++)
                if (bondIdx[row][col] >= 0)
                {
                    vector<int> bondGroup;
                    auto& bond = bonds[bondIdx[row][col]];
                    if (bond.bondCharacter == '-' || bond.bondCharacter == '=' || bond.bondCharacter == '#')
                    {
                        bondGroup.push_back(bondIdx[row][col]);
                        bool checkNext=true;
                        int col2 = col;
                        while (checkNext)
                        {
                            col2++;
                            if (col2 == nColumns)
                                break;
                            if (bondIdx[row][col2] < 0)
                                break;
                            char bondChar = bonds[bondIdx[row][col2]].bondCharacter;
                            if (bond.bondCharacter == bondChar)
                                bondGroup.push_back(bondIdx[row][col2]);
                            else
                                discamb::on_error::throwException("consecutive bond characters should be the same, e.g. --, not =-",
                                    __FILE__, __LINE__);
                        }
                        bondGroups.push_back(bondGroup);
                        bondsUpdated.push_back(bond);
                        bondsUpdated.back().length = bondGroup.size();
                        col += bondGroup.size() - 1;
                    }
                }

        // concatenate | 
        for (int col = 0; col < nColumns; col++)
            for (int row = 0; row < nRows; row++)
                if (bondIdx[row][col] >= 0)
                {
                    vector<int> bondGroup;
                    auto& bond = bonds[bondIdx[row][col]];
                    if (bond.bondCharacter == '|')
                    {
                        bondGroup.push_back(bondIdx[row][col]);
                        bool checkNext = true;
                        int row2 = row;
                        while (checkNext)
                        {
                            row2++;
                            if (row2 == nRows)
                                break;
                            if (bondIdx[row2][col] < 0)
                                break;
                            char bondChar = bonds[bondIdx[row2][col]].bondCharacter;
                            if (bond.bondCharacter == bondChar)
                                bondGroup.push_back(bondIdx[row2][col]);
                            else
                                discamb::on_error::throwException("consecutive bond characters should be the same, e.g. --, not =-",
                                    __FILE__, __LINE__);
                        }
                        bondGroups.push_back(bondGroup);
                        bondsUpdated.push_back(bond);
                        bondsUpdated.back().length = bondGroup.size();
                        row += bondGroup.size() - 1;
                    }
                }

        // concatenate  '\'
        for (int col = 0; col < nColumns; col++)
            for (int row = 0; row < nRows; row++)
                if (bondIdx[row][col] >= 0)
                {
                    vector<int> bondGroup;
                    auto& bond = bonds[bondIdx[row][col]];
                    if (bond.bondCharacter == '\\')
                    {
                        bondGroup.push_back(bondIdx[row][col]);
                        bondIdx[row][col] = -1;
                        bool checkNext = true;
                        int row2 = row;
                        int col2 = col;
                        while (checkNext)
                        {
                            row2++;
                            col2++;
                            if (row2 == nRows)
                                break;
                            if (col2 >= bondIdx[row2].size())
                                break;
                            if (bondIdx[row2][col2] < 0)
                                break;
                            char bondChar = bonds[bondIdx[row2][col2]].bondCharacter;
                            if (bond.bondCharacter == bondChar)
                                bondGroup.push_back(bondIdx[row2][col2]);
                            else
                                discamb::on_error::throwException("consecutive bond characters should be the same, e.g. --, not =-",
                                    __FILE__, __LINE__);
                            bondIdx[row2][col2] = -1;
                        }
                        bondGroups.push_back(bondGroup);
                        bondsUpdated.push_back(bond);
                        bondsUpdated.back().length = bondGroup.size();
                    }
                }

        // concatenate / 
        for (int row = 0; row < nRows; row++)
            for (int col = 0; col < nColumns; col++)
                if (bondIdx[row][col] >= 0)
                {
                    vector<int> bondGroup;
                    auto& bond = bonds[bondIdx[row][col]];
                    if (bond.bondCharacter == '/')
                    {
                        bondGroup.push_back(bondIdx[row][col]);
                        bondIdx[row][col] = -1;
                        bool checkNext = true;
                        int row2 = row;
                        int col2 = col;
                        while (checkNext)
                        {
                            row2++;
                            col2--;
                            if (row2 == nRows)
                                break;
                            if (col2 >= bondIdx[row2].size())
                                break;
                            if (col2 <0)
                                break;
                            if (bondIdx[row2][col2] < 0)
                                break;
                            char bondChar = bonds[bondIdx[row2][col2]].bondCharacter;
                            if (bond.bondCharacter == bondChar)
                                bondGroup.push_back(bondIdx[row2][col2]);
                            else
                                discamb::on_error::throwException("consecutive bond characters should be the same, e.g. --, not =-",
                                    __FILE__, __LINE__);
                            bondIdx[row2][col2] = -1;
                        }
                        bondGroups.push_back(bondGroup);
                        bondsUpdated.push_back(bond);
                        bondsUpdated.back().length = bondGroup.size();
                    }
                }


        vector<bool> inGroup(bonds.size(), false);
        for (auto& group : bondGroups)
            for (int bondIdx : group)
                inGroup[bondIdx] = true;

        for (int bondIdx = 0; bondIdx < nBonds; bondIdx++)
            if (!inGroup[bondIdx])
                bondsUpdated.push_back(bonds[bondIdx]);
                
        bonds.swap(bondsUpdated);
    }

    int firstOnTheLeft(const vector<int>& v, int pos)
    {
        if (pos == 0)
            return -1;

        if (v[pos - 1] >= 0)
            return v[pos - 1];
        return -1;
    }

    int firstOnTheRight(const vector<int>& v, int _pos, int length)
    {
        int pos = _pos + length - 1;
        if (pos + 1 > v.size() - 1)
            return -1;
        if (v[pos+1] >= 0)
                return v[pos + 1];
        return -1;
    }

    int valueIfExists(const vector<vector<int> >& v, int i, int j)
    {
        if (i < 0 || j < 0)
            return -1;
        if (v.size() > i)
            if (v[i].size() > j)
                return v[i][j];
        return -1;
    }

    void findConnectivity(
        const vector<AtomFromString> &atoms,
        const vector<BondFromString> &bonds,
        vector<vector<int> > &connectivity,
        vector<int> & nUnspecifiedNeighbours,
        int nLines)
    {
        
        int nAtoms = atoms.size();
        connectivity.clear();
        connectivity.resize(nAtoms);
        nUnspecifiedNeighbours.assign(nAtoms, 0);
        int maxX=0;

        for (const auto& atom : atoms)
            if (atom.labelEnd > maxX)
                maxX = atom.labelEnd;

        vector<vector<int> > atomIdx(nLines,vector<int>(maxX+1,-1));
        vector<vector<int> > wideAtomIdx(nLines, vector<int>(maxX + 1, -1));

        for (int i = 0; i < nAtoms; i++)
        {
            atomIdx[atoms[i].line][atoms[i].labelPosition] = i;
            for (int k = atoms[i].labelStart; k <= atoms[i].labelEnd; k++)
                wideAtomIdx[atoms[i].line][k] = i;
        }
        
        for (const auto& bond : bonds)
        {
            int atom1 = -2, atom2 = -2;
            if (bond.bondCharacter == '-' || bond.bondCharacter == '=' || bond.bondCharacter == '#')
            {
                atom1 = firstOnTheLeft(wideAtomIdx[bond.line], bond.position);
                atom2 = firstOnTheRight(wideAtomIdx[bond.line], bond.position, bond.length);
            } 
            else if (bond.bondCharacter == '|')
            {
                atom1 = valueIfExists(atomIdx, bond.line - 1, bond.position);
                atom2 = valueIfExists(atomIdx, bond.line + bond.length, bond.position);
            }
            else if (bond.bondCharacter == '\\')
            {
                atom1 = valueIfExists(atomIdx, bond.line - 1, bond.position - 1);
                atom2 = valueIfExists(atomIdx, bond.line + bond.length, bond.position + bond.length);
            }
            else if (bond.bondCharacter == '/')
            {
                atom1 = valueIfExists(atomIdx, bond.line - 1, bond.position + 1);
                atom2 = valueIfExists(atomIdx, bond.line + bond.length, bond.position - bond.length);
            }
            else if (bond.bondCharacter == '>')
            {
                atom1 = firstOnTheRight(wideAtomIdx[bond.line], bond.position, 1);
                atom2 = -1;
            }
            else if (bond.bondCharacter == '<')
            {
                atom1 = firstOnTheLeft(wideAtomIdx[bond.line], bond.position);
                atom2 = -1;
            }

            if (atom1 == -2 || atom2 == -2)
                discamb::on_error::throwException("bug", __FILE__, __LINE__);
            if (atom1 == -1 && atom2 == -1)
                discamb::on_error::throwException("error in atom type definition - bond not connected to any atom", __FILE__, __LINE__);
            if (atom1 == -1)
                nUnspecifiedNeighbours[atom2]++;
            if (atom2 == -1)
                nUnspecifiedNeighbours[atom1]++;
            if (bond.bondCharacter == '<' || bond.bondCharacter == '>')
                nUnspecifiedNeighbours[atom1]++;
            if (atom1 != -1 && atom2 != -1)
            {
                connectivity[atom1].push_back(atom2);
                connectivity[atom2].push_back(atom1);
            }

        }
    }

    

    void parseStructuralFormulaLine(
        const string &line,
        vector<AtomFromString> &atoms,
        vector<BondFromString> &bonds)
    {
        atoms.clear();
        bonds.clear();
        // bonds
        vector<char> bond_char{'-','=','#','|','\\','/','<','>'};
        vector<int> pos;
        BondFromString bond;
        bool inBracket = false;
        vector<bool> charInBracket(line.size());

        for (int i = 0; i < line.size(); i++)
        {
            char c = line[i];
            charInBracket[i] = inBracket;
            if (c == '(')
                inBracket = true;
            if (c == ')')
            {
                inBracket = false;
                charInBracket[i] = false;
            }
            
        }

        for (char c : bond_char)
        {
            discamb::string_utilities::findAllOccurences(line, c, pos);
            for (int idx : pos)
                if(!charInBracket[idx])
                {
                    bond.position = idx;
                    bond.bondCharacter = c;
                    bonds.push_back(bond);
                }
        }
        // atoms
        vector<char> splitters;
        splitters.swap(bond_char);
        splitters.push_back(' ');
        string atomString;
        bool inString = false;
        AtomFromString atomFromString;
        
        for (int i=0;i<line.size();i++)
        {
            char c = line[i];
            bool splitter = (find(splitters.begin(), splitters.end(), c) != splitters.end());

            bool bondInBracket = false;
            if (charInBracket[i])
                if (c == '-' || c == '=' || c == '#')
                    bondInBracket = true;
            if (splitter && inString && !bondInBracket)
            {
                inString = false;
                atomFromString.set(atomString, i-1);
                atoms.push_back(atomFromString);
                atomString.clear();
            }
            if (!splitter || bondInBracket)
            {
                inString = true;
                atomString += c;
            }
        }
        if (!atomString.empty())
        {
            atomFromString.set(atomString, line.size() - 1);
            atoms.push_back(atomFromString);
        }

    }

    void toAtomicNumbers(
        const string &_s,
        const std::map<std::string, std::set<int> >& element_groups, 
        int &neighborAtomicNumber,
        std::set<int> &neighborAtomicNumberRange)
    {
        neighborAtomicNumberRange.clear();
        // group, element, !group, !element
        bool notOfType = (_s[0] == '!');
        string s = _s;
        if (notOfType)
            s = s.substr(1);
        
        set<int> allElements;
        
        if (notOfType)
            for (int i = 0; i < 120; i++)
                allElements.insert(i);

        if (element_groups.find(s) != element_groups.end())
        {
            if (notOfType)
            {
                for (int z : element_groups.find(s)->second)
                    allElements.erase(z);
                neighborAtomicNumberRange = allElements;
            }
            else
                neighborAtomicNumberRange = element_groups.find(s)->second;
        }
        else
        {
            if (notOfType)
            {
                allElements.erase(discamb::periodic_table::atomicNumber(s));
                neighborAtomicNumberRange = allElements;
            }
            else
                neighborAtomicNumber = discamb::periodic_table::atomicNumber(s);
        }

    }

    string triboolToStr(const discamb::Tribool& tbl)
    {
        if (tbl == discamb::Tribool::False)
            return "-";
        if (tbl == discamb::Tribool::True)
            return "+";
        return "*";
    }

    void parseProperties(
        const vector<std::string> &properties,
        discamb::AtomDescriptors& descriptors)
    {   
        // 6, 6A, !6A, !3, 3
        string ring5plusStr;// = triboolToStr(descriptors.ringInfo.inRing);
        string ring3str = triboolToStr(descriptors.ringInfo.in3Ring);
        string ring4str = triboolToStr(descriptors.ringInfo.in4Ring);

        for (const string& p : properties)
        {
            string s;
            vector<string> words;
            discamb::string_utilities::split(p, words,'=');
            if (words.size() == 2)
                s = words[1];
            else
                s = p;

            if (s == "any")
            {
                descriptors.fixedNumberOfNeighbors = false;
                descriptors.planar = discamb::Tribool::Undefined;
                descriptors.ringInfo.in3Ring = discamb::Tribool::Undefined;
                descriptors.ringInfo.in4Ring = discamb::Tribool::Undefined;
                descriptors.ringInfo.inRing = discamb::Tribool::Undefined;
            }
            else if (s == "planar")
                descriptors.planar = discamb::Tribool::True;
            else if (s == "nonplanar")
                descriptors.planar = discamb::Tribool::False;
            else if (s == "!3")
                ring3str = "-";
            else if (s == "3")
                ring3str = "+";
            else if (s == "3*")
                ring3str = "*";
            else if (s == "!4")
                ring4str = "-";
            else if (s == "4")
                ring4str = "+";
            else if (s == "4*")
                ring4str = "*";
            else if (s == "!Ar")
            {
                if (!ring5plusStr.empty())
                    ring5plusStr += ",";
                ring5plusStr += "-";
            }
            else if (s == "Ar")
            {
                if (!ring5plusStr.empty())
                    ring5plusStr += ",";
                ring5plusStr += "+";
            }
            else if (ring5plusStr.empty())
                ring5plusStr += s;
            else
                ring5plusStr += "," + s;
        }
        if (ring5plusStr.empty())
            ring5plusStr = triboolToStr(descriptors.ringInfo.inRing);
        descriptors.ringInfo.set(ring5plusStr, ring3str, ring4str);
    }

    int labeledAtomsFromBrackets(
        std::vector<AtomFromString>& atoms,
        std::vector<std::vector<int> >& connectivity)
    {
        int nAtoms = atoms.size();
        int nNewAtoms = 0;
        for (int atomIdx=0; atomIdx<nAtoms; atomIdx++)
        {
            vector<int> additional_bonded_atoms_to_keep;
            for (int j = 0; j < atoms[atomIdx].additional_bonded_atoms.size(); j++)
            {
                string& word = atoms[atomIdx].additional_bonded_atoms[j];
                if (isdigit(word.back()) && !isdigit(word.front()))
                {
                    if (word[0] == '-' || word[0] == '=' || word[0] == '#')
                        word = word.substr(1);
                    AtomFromString atom;
                    atom.set(word, -1);
                    atom.inBracket = true;
                    atoms.push_back(atom);
                    connectivity[atomIdx].push_back(connectivity.size());
                    connectivity.resize(connectivity.size() + 1);
                    connectivity.back().push_back(atomIdx);
                    nNewAtoms++;
                }
                else
                    additional_bonded_atoms_to_keep.push_back(j);
            }
            vector<string> additional_bonded_atoms_new;
            for (int idx : additional_bonded_atoms_to_keep)
                additional_bonded_atoms_new.push_back(atoms[atomIdx].additional_bonded_atoms[idx]);
            atoms[atomIdx].additional_bonded_atoms = additional_bonded_atoms_new;



        }
        return nNewAtoms;
    }

    string jsonValueAsLowerCaseString(
        const json& data,
        const string& key,
        const string& defaultValue)
    {
        string result = defaultValue;
        if (data.find(key) != data.end())
        {
            if (data.find(key)->is_boolean())
                result = (data[key].get<bool>() ? "true" : "false");
            else
                result = data[key].get<string>();
        }
        return discamb::string_utilities::toLower(result);
    }

    discamb::Tribool jsonToTribool(const json& data,
        const string& key,
        const string& defaultValue)
    {
        string s = jsonValueAsLowerCaseString(data, key, defaultValue);
        if (s == "true" || s == "+")
            return discamb::Tribool::True;
        else if (s == "false" || s == "-")
            return discamb::Tribool::False;
        else if (s == "undefined" || s == "*")
            return discamb::Tribool::Undefined;
        else
            discamb::on_error::throwException("cannot parse atom type descriptor with key '" + key + "'", __FILE__, __LINE__);

        return discamb::Tribool::Undefined;
    }

    void checkForRepeatingTypeIds(
        const std::vector<discamb::AtomType>& atomTypes)
    {
        set<string> repeatingIds;
        vector<string> ids;
        for (const auto& type : atomTypes)
            ids.push_back(type.id);
        for (string id : ids)
            if (count(ids.begin(), ids.end(), id) > 1)
                repeatingIds.insert(id);

        if (!repeatingIds.empty())
        {
            string message = "there are type definitions with repeating ids: ";
            for (string id : repeatingIds)
                message += " " + id;
            discamb::on_error::throwException(message, __FILE__, __LINE__);
        }
    }

    // checks if it is a label beginning (uppercase letter not preceded by letter)
    bool isLabelBeginningOrEmpty(const string line, int position)
    {
        if (position >= line.size())
            return true;

        char c = line[position]; 

        if(c == ' ')
            return true;

        if (!isalpha(c))
            return false;

        if (c != toupper(c))
            return false;
        if (position > 0)
            if (isalpha(line[position - 1]))
                return false;
        return true;
    }

    void checkVerticalBondsInStructuralFormula(
        const std::vector<std::string>& lines)
    {
        int nLines = lines.size();
        bool ok = true;
        for (int lineIdx = 0; lineIdx < nLines; lineIdx++)
        {
            int lineSize = lines[lineIdx].size();
            for (int charIdx = 0; charIdx < lineSize; charIdx++)
            {
                char c = lines[lineIdx][charIdx];
                if (c == '|')
                {
                    // check above
                    if (lineIdx > 0)
                        if(!isLabelBeginningOrEmpty(lines[lineIdx - 1], charIdx))
                            ok = false;
                    // check below
                    if (lineIdx < nLines - 1)
                        if (!isLabelBeginningOrEmpty(lines[lineIdx + 1], charIdx))
                            ok = false;
                    
                }
            }
        }
        if(!ok)
            discamb::on_error::throwException("vertical bonds '|' in structural formula should end with chemical element symbol or space", __FILE__, __LINE__);
    }

}

namespace discamb {

    namespace atom_type_io {


        void readAtomTypes(
            const std::string& fileName,
            std::vector<AtomType>& atomTypes,
            DescriptorsSettings& descriptorsSettings)
        {
            if(fileName.size()>5)
                if (fileName.substr(fileName.size() - 5) == string(".json"))
                {
                    readAtomTypesJson(fileName, atomTypes, descriptorsSettings);
                    return;
                }

            MATTS_BankReader reader;
            vector<AtomTypeHC_Parameters> parameters;
            BankSettings bankSettings;
            reader.read(fileName, atomTypes, parameters, bankSettings);
            descriptorsSettings = bankSettings.descriptorsSettings;
        }

        // checks if vertical bonds '|' end with chemical element symbol or space


        void parseStructuralFormulaTxt(
            const std::vector<std::string>& lines,
            const std::map<std::string, std::set<int> >& element_groups,
            const AtomDescriptors& defaultDescriptors,
            const AtomDescriptors& defaultCentralAtomDescriptors,
            std::vector<AtomDescriptors>& atoms, 
            std::vector<std::vector<int> >& connectivityFinal)
        {
            atoms.clear();
            connectivityFinal.clear();

            checkVerticalBondsInStructuralFormula(lines);

            vector<AtomFromString> atomsFromString;
            vector<BondFromString> bondsFromString;


            for (int lineIdx=0; lineIdx<lines.size(); lineIdx++)
            {
                vector<AtomFromString> newAtoms;
                vector<BondFromString> newBonds;

                parseStructuralFormulaLine(lines[lineIdx], newAtoms, newBonds);
                for (auto& atom : newAtoms)
                {
                    atom.line = lineIdx;
                    atomsFromString.push_back(atom);
                }
                for (auto& bond : newBonds)
                {
                    bond.line = lineIdx;
                    bondsFromString.push_back(bond);
                }
            }
            concatenateBonds(bondsFromString);


            vector<vector<int> > connectivity;
            vector<int> nUnspecifiedNeighbours;
            findConnectivity(atomsFromString, bondsFromString, connectivity, nUnspecifiedNeighbours, lines.size());
            int nAtomsFromBrackets = labeledAtomsFromBrackets(atomsFromString, connectivity);
            nUnspecifiedNeighbours.resize(nUnspecifiedNeighbours.size() + nAtomsFromBrackets);
            // un/named atoms
            vector<int> namedAtoms, unnamedAtoms, namedAtomIdx;
            vector<string> namedAtomLabels;
            int atomIdx, nAtoms = atomsFromString.size();
            vector<bool> isNamed(nAtoms, false);

            for(int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (isdigit(atomsFromString[atomIdx].name.back()))
                {
                    namedAtoms.push_back(atomIdx);
                    namedAtomLabels.push_back(atomsFromString[atomIdx].name);
                    string idxStr;
                    for (char c : atomsFromString[atomIdx].name)
                        if (isdigit(c))
                            idxStr += c;
                    namedAtomIdx.push_back(stoi(idxStr));
                    isNamed[atomIdx] = true;
                }
                else
                    unnamedAtoms.push_back(atomIdx);

            //###### check for correcness

            string bondingString;
            for (const auto& line : lines)
                bondingString += line + "\n";

            // check that there is one and only one atom with idx 1
            int nOne = count(namedAtomIdx.begin(), namedAtomIdx.end(), 1);
            if (nOne != 1)
            {
                string s;
                for (const auto& line : lines)
                    s += line + "\n";
                if (nOne == 0)
                    on_error::throwException("invalid naming in atom type definition - no atom with idx 1 found:\n" + bondingString, __FILE__, __LINE__);
                if (nOne > 1)
                    on_error::throwException("invalid naming in atom type definition - more than one atom with idx 1 found:\n" + bondingString, __FILE__, __LINE__);
            }
            // check that there are no repeating atom labels
            set<string> uniqueLabels(namedAtomLabels.begin(), namedAtomLabels.end());
            if (uniqueLabels.size() != namedAtomLabels.size())
                on_error::throwException("non unique labels in atom type definition:\n" + bondingString, __FILE__, __LINE__);
            
            //
            set<int> allElements;
            for (int i = 0; i < 120; i++)
                allElements.insert(i);
            vector<int> newIdx2Old;
            map<int, int> oldIdx2New;
            for(atomIdx=0; atomIdx<nAtoms; atomIdx++)
                if (isNamed[atomIdx])
                {
                    AtomDescriptors atomDescriptors = defaultDescriptors;
                    // check if central atom and if yes assign to defaultCentralAtomDescriptors
                    string numberStr;
                    for (char c : atomsFromString[atomIdx].name)
                        if (isdigit(c))
                            numberStr += c;
                    if (!numberStr.empty())
                        if (stoi(numberStr) == 1)
                            atomDescriptors = defaultCentralAtomDescriptors;

                    // sets AtomDescriptors::anyAtomicNumber
                    atomDescriptors.anyAtomicNumber = (atomsFromString[atomIdx].type == string("X"));
                    
                    // sets AtomDescriptors::atomic_number_range 
                    if (element_groups.find(atomsFromString[atomIdx].type) != element_groups.end())
                        atomDescriptors.atomic_number_range = element_groups.find(atomsFromString[atomIdx].type)->second;
                    else if (!atomDescriptors.anyAtomicNumber)
                        atomDescriptors.atomic_number_range.insert(periodic_table::atomicNumber(atomsFromString[atomIdx].type));
                    
                    // sets AtomDescriptors::fixedNumberOfNeighbors
                    const vector<string>& v_str = atomsFromString[atomIdx].additional_bonded_atoms;
                    atomDescriptors.fixedNumberOfNeighbors = (find(v_str.begin(), v_str.end(), "*") == v_str.end());

                    // sets AtomDescriptors::label
                    atomDescriptors.label = atomsFromString[atomIdx].name;
                    
                    // process atomsFromString[atomIdx].additional_bonded_atoms
                    atomDescriptors.nNeighbours = 0;

                    for (const string& _s : atomsFromString[atomIdx].additional_bonded_atoms)
                    {
                        // group, element, !group, !element, X, *, Ngroup (e.g. 3Ch),Nelement (e.g. 3H)
                        int nCopies = 1;
                        string s;
                        if (isdigit(_s[0]))
                        {
                            nCopies = stoi(_s.substr(0, 1));
                            s = _s.substr(1);
                            if (s.empty())
                                s = "X";
                        }
                        else
                            s = _s;
                        
                        if (s[0] == '*')
                            atomDescriptors.fixedNumberOfNeighbors = false;
                        else
                            atomDescriptors.nNeighbours += nCopies;
                        
                        if (s != string("X") && s[0] != '*')
                        {
                            int neighborAtomicNumber;
                            std::set<int> neighborAtomicNumberRange;
                            //int nCopies = 1;
                            //if (isdigit(s[0]))
                            //{
                            //    nCopies = stoi(s.substr(0, 1));
                            //    s = s.substr(1);
                            //}
                            toAtomicNumbers(s, element_groups, neighborAtomicNumber, neighborAtomicNumberRange);

                            for (int i = 0; i < nCopies; i++)
                                if (neighborAtomicNumberRange.empty())
                                    //atomDescriptors.neighborsAtomicNumbers.push_back(neighborAtomicNumber);
                                    atomDescriptors.neighborsAtomicNumbers.insert(neighborAtomicNumber);
                                else
                                    atomDescriptors.neighborsAtomicNumberRanges.push_back(neighborAtomicNumberRange);
                        }

                    }
                    atoms.push_back(atomDescriptors);
                    oldIdx2New[atomIdx] = newIdx2Old.size();
                    newIdx2Old.push_back(atomIdx);
                }

            // finish setting atomDescriptors.neighborsAtomicNumbers, atomDescriptors.neighborsAtomicNumberRanges
            // and atomDescriptors.nNeighbours using info about neighbors

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (isNamed[atomIdx])
                {
                    AtomDescriptors& currentAtom = atoms[oldIdx2New[atomIdx]];
                    currentAtom.nNeighbours += nUnspecifiedNeighbours[atomIdx] + connectivity[atomIdx].size();
                    for (int neighborIdx : connectivity[atomIdx])
                        if (isNamed[neighborIdx])
                        {
                            AtomDescriptors& neighborAtom = atoms[oldIdx2New[neighborIdx]];
                            if (neighborAtom.atomic_number_range.size() == 1)
                            {
                                //currentAtom.neighborsAtomicNumbers.push_back(*neighborAtom.atomic_number_range.begin());
                                if (*neighborAtom.atomic_number_range.begin() != 0)
                                    currentAtom.neighborsAtomicNumbers.insert(*neighborAtom.atomic_number_range.begin());
                            }
                            else
                                currentAtom.neighborsAtomicNumberRanges.push_back(neighborAtom.atomic_number_range);
                        }
                        else
                        {
                            // group, element, !group, !element
                            int neighborAtomicNumber;
                            std::set<int> neighborAtomicNumberRange;
                            toAtomicNumbers(atomsFromString[neighborIdx].type, element_groups, neighborAtomicNumber, neighborAtomicNumberRange);
                            if (neighborAtomicNumberRange.empty())
                            {
                                //currentAtom.neighborsAtomicNumbers.push_back(neighborAtomicNumber);
                                if(neighborAtomicNumber!=0)
                                    currentAtom.neighborsAtomicNumbers.insert(neighborAtomicNumber);
                            }
                            else
                                currentAtom.neighborsAtomicNumberRanges.push_back(neighborAtomicNumberRange);

                        }
                }

            // parse properties
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                if (isNamed[atomIdx])
                {
                    AtomDescriptors& currentAtom = atoms[oldIdx2New[atomIdx]];
                    // this is already set by defaultDescriptors ("default atomic properties" in json bank file)
                    //currentAtom.planar = Tribool::Undefined;
                    //currentAtom.ringInfo.in3Ring = Tribool::Undefined;
                    //currentAtom.ringInfo.in4Ring = Tribool::Undefined;
                    //currentAtom.ringInfo.inRing = Tribool::Undefined;

                    parseProperties(atomsFromString[atomIdx].properties, currentAtom);
                }

            // find connectivity for named atoms only
            //vector<int> newIdx2Old;
            //map<int, int> oldIdx2New;
            int nNamedAtoms = newIdx2Old.size();
            connectivityFinal.resize(nNamedAtoms);
            
            for (int i = 0; i < nNamedAtoms; i++)
            {
                int oldIdx = newIdx2Old[i];
                for (int j = 0; j < connectivity[oldIdx].size(); j++)
                {
                    int neighbourOldIdx = connectivity[oldIdx][j];
                    if (isNamed[neighbourOldIdx])
                        connectivityFinal[i].push_back(oldIdx2New[neighbourOldIdx]);
                }
            }
            // set atom with label element1 e.g. C1 at the beginning of the list
            int idxOfCentralAtom = -1;
            for (int i = 0; i < nNamedAtoms; i++)
            {
                string numberStr;
                for (char c : atoms[i].label)
                    if (isdigit(c))
                        numberStr += c;
                if(!numberStr.empty())
                    if (stoi(numberStr) == 1)
                    {
                        if (idxOfCentralAtom != -1)
                            on_error::throwException("more than one atom with index 1 in atom type definition", __FILE__, __LINE__);
                        else
                            idxOfCentralAtom = i;
                    }
            }

            if (idxOfCentralAtom != 0)
            {
                AtomDescriptors atom = atoms[0];
                atoms[0] = atoms[idxOfCentralAtom];
                atoms[idxOfCentralAtom] = atom;

                vector<int> neighbours = connectivityFinal[0];
                connectivityFinal[0] = connectivityFinal[idxOfCentralAtom];
                connectivityFinal[idxOfCentralAtom] = neighbours;

                for(auto &v: connectivityFinal)
                    for (auto& idx : v)
                    {
                        if (idx == 0)
                            idx = idxOfCentralAtom;
                        else if (idx == idxOfCentralAtom)
                            idx = 0;
                    }
            }
        }

        void json2atomDescriptors(
            const nlohmann::json& data,
            AtomDescriptors& atomDescriptors)
        {
            atomDescriptors.planar = jsonToTribool(data,"planar", "*");
            atomDescriptors.ringInfo.inRing =  jsonToTribool(data,"in 5 or 6 ring", "-");
            atomDescriptors.ringInfo.in4Ring =   jsonToTribool(data,"in 4 ring", "-");
            atomDescriptors.ringInfo.in3Ring =   jsonToTribool(data,"in 3 ring", "-");

        }

        void json2atomType(
            const nlohmann::json& data,
            const std::map<std::string, std::set<int> >& element_groups,
            const string atomTypeId,
            AtomType& atomType,
            const AtomDescriptors& defaultAtomDescriptors,
            const AtomDescriptors& defaultCetralAtomDescriptors)
        {
            atomType = AtomType();
            atomType.id = atomTypeId;
            try {

                if (data.find("atoms and bonds") != data.end())
                {
                    nlohmann::json atomsAndBonds = data.find("atoms and bonds").value();
                    if (atomsAndBonds.begin()->is_string())
                    {
                        vector<string> atomsAndBondsLines;
                        for (auto it = atomsAndBonds.begin(); it != atomsAndBonds.end(); it++)
                            atomsAndBondsLines.push_back(it->get<string>());
                        parseStructuralFormulaTxt(atomsAndBondsLines,
                            element_groups,
                            defaultAtomDescriptors,
                            defaultCetralAtomDescriptors,
                            atomType.atoms,
                            atomType.connectivity);

                    }
                    else
                    {
                        on_error::not_implemented(__FILE__, __LINE__);
                    }
                }
            }
            catch (...)
            {
                on_error::throwException("problem when parsing \"atoms and bonds\" for atom type " + atomType.id, __FILE__, __LINE__);
            }
            /*
        std::vector<std::string> ringLabels;

            std::vector<int> ringSizes;

            */

            set<string> ringLabels;
            for (auto& atom : atomType.atoms)
            {
                for (const auto& ringInfo : atom.ringInfo.labeledContainingRings)
                    ringLabels.insert(ringInfo.second);
                for (const auto& ringInfo : atom.ringInfo.labeledNonContainingRings)
                    ringLabels.insert(ringInfo.second);
            }
            atomType.ringLabels.clear();
            atomType.ringLabels.insert(atomType.ringLabels.end(), ringLabels.begin(), ringLabels.end());
            atomType.ringSizes.clear();
            for (string& label : atomType.ringLabels)
                atomType.ringSizes.push_back(stoi(label.substr(0, 1)));
            if (data.find("lcs") != data.end())
            {
                string lcsString = data.find("lcs")->get<string>();
                atomType.setLocalCoordinateSystem(lcsString + " R");
            }
            else
                if (!atomType.setDefaultLocalCoordinateSystem())
                    on_error::throwException("cannot set default local coordinate sytem for atom type " + atomTypeId, __FILE__, __LINE__);
            

        }

        void json2elementGroups(
            const nlohmann::json& data,
            std::map<std::string, std::set<int> >& elementGroups)
        {
            elementGroups.clear();

            if (!data.is_object())
            {
                string msg = string("problem when interpreting atom types from JSON, ") +
                    "\"element groups\" is expected to be an object";
                on_error::throwException(msg, __FILE__, __LINE__);
            }
            else
                for (const auto group : data.items())
                {
                    std::set<int> atomicNumbers;
                    basic_chemistry_utilities::getElementsList(group.value().get<string>(), atomicNumbers);
                    elementGroups[group.key()] = atomicNumbers;
                }

            
        }

        void readAtomTypesJson(
            const std::string& jsonFileName,
            std::vector<AtomType>& atomTypes,
            std::vector<nlohmann::json>& typeData,
            DescriptorsSettings& descriptorsSettings)
        {
            atomTypes.clear();
            typeData.clear();

            // read json file

            ifstream in(jsonFileName);

            if (!in.good())
                on_error::throwException(string("cannot read file '") + jsonFileName + string("'"), __FILE__, __LINE__);
            JsonDuplicateCheck jsonDuplicateCheck;
            json data = json::parse(in, jsonDuplicateCheck);
            in.close();

            // set element groups

            map<string, set<int> > element_groups;
            if (data.find("element groups") != data.end())
                json2elementGroups(data["element groups"], element_groups);
            // 
            AtomDescriptors defaultAtomDescriptors, defaultCetralAtomDescriptors;
            if (data.find("default atomic properties") != data.end())
                json2atomDescriptors(data.find("default atomic properties").value(), defaultAtomDescriptors);
            if (data.find("default central atom properties") != data.end())
                json2atomDescriptors(data.find("default central atom properties").value(), defaultCetralAtomDescriptors);

            
            set<string> typeLabels;
            if (data.find("types") != data.end())
            {
                auto types = data.find("types");
                for (auto type : types->items())
                    if (type.key() != string("comment"))
                    {
                        AtomType atomType;
                        json2atomType(type.value(), element_groups, type.key(), atomType, defaultAtomDescriptors, defaultCetralAtomDescriptors);
                        //atomType.id = type.key();
                        //atomTypesAndData.push_back({ atomType, type.value()});
                        if(typeLabels.count(atomType.id)!=0)
                            on_error::throwException("repeating atom type id '" + atomType.id + "' in json file", __FILE__, __LINE__);
                        typeLabels.insert(atomType.id);
                        atomTypes.push_back(atomType);
                        typeData.push_back(type.value());
                    }
            }
            checkForRepeatingTypeIds(atomTypes);
            if (data.find("settings") != data.end())
            {
                vector<int> upgradeUnnamedAtoms = data["settings"].value("upgrade unnamed atoms", vector<int>());
                if (!upgradeUnnamedAtoms.empty())
                    for (auto& type : atomTypes)
                        for (int z : upgradeUnnamedAtoms)
                            type.transformUnnamedAtomsToNamed(z, defaultAtomDescriptors);
            }
        }

        void readAtomTypesJson(
            const std::string& jsonFileName,
            std::vector<AtomType>& atomTypes,
            DescriptorsSettings& descriptorsSettings)
        {
            vector<nlohmann::json> typesData;
            readAtomTypesJson(jsonFileName, atomTypes, typesData, descriptorsSettings);


            //// read json file

            //ifstream in(jsonFileName);

            //if (!in.good())
            //    on_error::throwException(string("cannot read file '") + jsonFileName + string("'"), __FILE__, __LINE__);

            //json data = json::parse(in);
            //in.close();

            //// set element groups

            //map<string, set<int> > element_groups;
            //if (data.find("element groups") != data.end())
            //    json2elementGroups(data["element groups"], element_groups);
            //// 
            //AtomDescriptors defaultAtomDescriptors;
            //if (data.find("default atomic properties") != data.end())
            //    json2atomDescriptors(data.find("default atomic properties").value(), defaultAtomDescriptors);

            //if (data.find("types") != data.end())
            //{
            //    auto types = data.find("types");
            //    for (auto type : types->items())
            //        if(type.key() != string("comment"))
            //        {
            //            AtomType atomType;
            //            json2atomType(type.value(), element_groups, type.key(), atomType, defaultAtomDescriptors);
            //            //atomType.id = type.key();
            //            atomTypes.push_back(atomType);
            //        }
            //}

            //checkForRepeatingTypeIds(atomTypes);
        }

        void writeAtomTypesTxt(
            const std::string& fileName,
            const std::vector<AtomType>& atomTypes,
            const std::map<std::string, std::set<int> >& element_groups)
        {
            ofstream out(fileName);

            if (!out.good())
                on_error::throwException("cannot write atom types to file '" + fileName + "'", __FILE__, __LINE__);

            /*
            ENTRY
  ID
    H101
  COMMENT 
    H101 x-cH3 
    N 2001?  DATE Fri Apr 22 15:44:52 2011  
  ATOM DESCRIPTORS
#  central atom
    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -
#  1-st neighbors
    C2    CONNECTED_TO  H1,H,H,!H     PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *

            */

            for (const auto &type : atomTypes)
            {

                out << "ENTRY\n"
                    << "  ID\n"
                    << "    " << type.id << "\n"
                    << "  ATOM DESCRIPTORS\n";
                for (int atomIdx = 0; atomIdx< type.atoms.size(); atomIdx++)
                {
                    const auto& atom = type.atoms[atomIdx];
                    out << "    " << left << setw(6) << atom.label << "CONNECTED_TO  ";
                    string neighboursStr;
                    //vector<int> neighborsAtomicNumbers, labeledNeighborsAtomicNumbers;
                    multiset<int> neighborsAtomicNumbers, labeledNeighborsAtomicNumbers;
                    vector<string> connectivityStrVec;
                    for (int neighbourIdx : type.connectivity[atomIdx])
                    {
                        connectivityStrVec.push_back(type.atoms[neighbourIdx].label);
                        if (!type.atoms[neighbourIdx].anyAtomicNumber)
                            if (type.atoms[neighbourIdx].atomic_number_range.size() == 1)
                                labeledNeighborsAtomicNumbers.insert(*type.atoms[neighbourIdx].atomic_number_range.begin());
                                //labeledNeighborsAtomicNumbers.push_back(*type.atoms[neighbourIdx].atomic_number_range.begin());

                    }
                    neighborsAtomicNumbers = atom.neighborsAtomicNumbers;
                    for (auto z : labeledNeighborsAtomicNumbers)
                        neighborsAtomicNumbers.erase(find(neighborsAtomicNumbers.begin(), neighborsAtomicNumbers.end(),z));

                    for (auto z : neighborsAtomicNumbers)
                        connectivityStrVec.push_back(periodic_table::symbol(z));

                    for (int i = connectivityStrVec.size(); i < atom.nNeighbours; i++)
                        connectivityStrVec.push_back("X");

                    if (!atom.fixedNumberOfNeighbors)
                        connectivityStrVec.push_back("*");
                    if(connectivityStrVec.empty())
                        connectivityStrVec.push_back("-");

                    string connectivityStr = connectivityStrVec[0];
                    for (int i = 1; i < connectivityStrVec.size(); i++)
                        connectivityStr += "," + connectivityStrVec[i];
                    out << setw(14) << connectivityStr;
//H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -
                     
                    // PLANARITY
                    
                    out << "PLANARITY  " << triboolToStr(atom.planar) << "  PLANAR_RING_WITH_PLANAR_ATOMS ";
                    
                    // PLANAR_RING_WITH_PLANAR_ATOMS
                    
                    vector<string> ringStrVec;
                    for (const auto& ring : atom.ringInfo.labeledContainingRings)
                        ringStrVec.push_back(ring.second);

                    for (const auto& ring : atom.ringInfo.labeledNonContainingRings)
                        ringStrVec.push_back("!"+ring.second);

                    for (const auto& ring : atom.ringInfo.nonLabeledContainingRings)
                        ringStrVec.push_back(std::to_string(ring));

                    for (const auto& ring : atom.ringInfo.nonLabeledNonContainingRings)
                        ringStrVec.push_back("!" + std::to_string(ring));
                    
                    if (ringStrVec.empty())
                        ringStrVec.push_back(triboolToStr(atom.ringInfo.inRing));

                    string s = ringStrVec[0];
                    for (int i = 1; i < ringStrVec.size(); i++)
                        s += "," + ringStrVec[i];
                    out << setw(8) << s;

                    out << "IN_3_MEMBER_RING " << triboolToStr(atom.ringInfo.in3Ring)
                        << " IN_4_MEMBER_RING " << triboolToStr(atom.ringInfo.in4Ring) << "\n";


                }
                out << "\n\n";
            }
        }
    }
}
