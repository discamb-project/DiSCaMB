#include "discamb/IO/shelx_io.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"

#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <iomanip>


using namespace std;

namespace {

    double get_precision(const string& _s)
    {
        // possible formats 12, 0.12
        string s = _s;
        vector<string> words;
        discamb::string_utilities::split(s, words, '.');
        if (words.size() == 1)
            return 1.0;
        s = string("0.") + string(words[1].size() - 1, '0') + string("1");
        return stod(s);
    }


    void processParametersLine(
        vector<string>::const_iterator const &start,
        vector<string>::const_iterator const &end,
        const vector<double> &fvars,
        vector<double> &params,
        vector<double>& precision)
    {
        double v, k, p;
        int freeVariableIdx;
        params.clear();
        precision.clear();
        for (auto it = start; it != end; it++)
        {
            v = stod(*it);
            precision.push_back(get_precision(*it));
            if (v > 5 && v <= 15)
                v -= 10;
            else
                if (fabs(v) > 15)
                { // v = 10*freeVariableIdx + p => p = v - 10*freeVariableIdx
                    // 31

                    k = floor ((fabs(v) + 5.0) / 10);
                    p = fabs(v) - 10 * k;
                    freeVariableIdx = int(k) - 2;
                    
                    if (freeVariableIdx + 1 > fvars.size())
                        discamb::on_error::throwException(string("problem when reading shelx file, index of free variable out of range, index deduced from parameter value : ") + to_string(freeVariableIdx + 1), __FILE__, __LINE__);

                    if (v > 0)
                        v = p * fvars[freeVariableIdx];
                    else
                        v = p * (1 - fvars[freeVariableIdx]);
                }
            params.push_back(v);
        }
    }
}

namespace discamb {
    namespace shelx_io {

        void save(
            const std::string &fName,
            const Crystal &crystal)
        {
            ofstream out(fName);

            if (!out.good())
                on_error::throwException(string("can not write to '") + fName + string("'"), __FILE__, __LINE__);

            string title;
            out << "TITL some_title\n"
                << "CELL 0.71073"
                << " " << setprecision(4) << fixed << crystal.unitCell.a()
                << " " << setprecision(4) << fixed << crystal.unitCell.b()
                << " " << setprecision(4) << fixed << crystal.unitCell.c()
                << " " << setprecision(4) << fixed << crystal.unitCell.alpha()
                << " " << setprecision(4) << fixed << crystal.unitCell.beta()
                << " " << setprecision(4) << fixed << crystal.unitCell.gamma()
                << "\nZERR " << crystal.spaceGroup.nSymmetryOperations() << " 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000\n";


            bool isObverse;
            char centering = crystal.spaceGroup.centering(isObverse);
            map<char, int> centeringAsIndex{ {'P', 1}, {'I', 2}, {'R', 3}, {'F', 4}, {'A', 5}, {'B', 6}, {'C', 7} };
            if (centering == 'H')
                on_error::throwException("hexagonal centering not supported for save to shelx", __FILE__, __LINE__);
            if (centering == 'R')
                if (!isObverse)
                    on_error::throwException("for rombohedral lattice only rhombohedral obverse is supported", __FILE__, __LINE__);

            int latticeIdx = centeringAsIndex[centering];

            if (!crystal.spaceGroup.isCentrosymmetric())
                latticeIdx *= -1;

            //if (latticeIdx != -1)
            out << "LATT " << latticeIdx << "\n";

            int i, n = crystal.spaceGroup.nSymmetryOperationsInSubset();
            string symmOperationAsString;
            Matrix3i m, identity;
            Vector3<CrystallographicRational> t, zero;
            identity.setToIdentity();
            for (i = 0; i < n; i++)
            {
                crystal.spaceGroup.getSpaceGroupOperation(0, 0, i).get(symmOperationAsString);
                crystal.spaceGroup.getSpaceGroupOperation(0, 0, i).get(m, t);
                if(!( m==identity && t == zero ) )
                    out << "SYMM " << symmOperationAsString << "\n";
            }
            int nAtoms = crystal.atoms.size();
            vector<int> atomicNumber;
            crystal_structure_utilities::atomicNumbers(crystal, atomicNumber);

            set<int> uniqueZ;
            for (auto z : atomicNumber)
                uniqueZ.insert(z);
            out << "SFAC";
            map<int, int> z2sfacIdx;
            int sfacIdx = 1;
            for (auto z : uniqueZ)
            {
                z2sfacIdx[z] = sfacIdx++;
                out << " " << periodic_table::symbol(z);
            }
            out << "\nUNIT";
            map<int, int> unit;

            for (i = 0; i < nAtoms; i++)
            {
                auto it = unit.find(atomicNumber[i]);
                if (it == unit.end())
                    unit[atomicNumber[i]] = crystal.atoms[i].multiplicity;
                else
                    it->second += crystal.atoms[i].multiplicity;
            }
            for (auto x : unit)
                out << " " << x.second;
            out << "\n" << std::left;

            int nSymmOperations = crystal.spaceGroup.nSymmetryOperations();

            for (i = 0; i < nAtoms; i++)
            {
                out << setw(8) << crystal.atoms[i].label << z2sfacIdx[atomicNumber[i]]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].coordinates[0]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].coordinates[1]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].coordinates[2]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].occupancy * 
                                                          crystal.atoms[i].multiplicity / static_cast<double>(nSymmOperations);
                if (crystal.atoms[i].adp.size() == 0)
                    out << " 0.000000\n";
                if (crystal.atoms[i].adp.size() == 1)
                    out << " " << setprecision(6) << fixed << crystal.atoms[i].adp[0] << "\n";
                if (crystal.atoms[i].adp.size() == 6)
                    out << " " << setprecision(6) << fixed << crystal.atoms[i].adp[0]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].adp[1] << " =\n        "
                    << " " << setprecision(6) << fixed << crystal.atoms[i].adp[2]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].adp[5]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].adp[4]
                    << " " << setprecision(6) << fixed << crystal.atoms[i].adp[3] << "\n";

            }

            out << "\nEND\n";
            out.close();

        }
        //void read(const std::string& fName, Crystal& crystal, std::map<std::string, std::string> &data);
        void read(
            const std::string& fName,
            Crystal& crystal)
        {
            map<string, string> data;
            read(fName, crystal, data);
        }


        void read(
            const std::string &fName,
            Crystal &crystal,
            std::map<std::string, std::string>& data)
        {
            ifstream in(fName);
            read(in, crystal, data);
            in.close();
        } 

        void read(
            std::istream& input,
            Crystal& crystal)
        {
            map<string, string> data;
            read(input, crystal, data);
        }

        void read(
            std::istream& input,
            Crystal& crystal,
            std::map<std::string, std::string>& data)
        {
            // list of SHELXL instructions from  http://shelx.uni-goettingen.de/shelxl_html.php
            vector<string> instructions = { "ABIN","ACTA","AFIX","ANIS","ANSC","ANSR","BASF","BIND",
                "BLOC","BOND","BUMP","CELL","CGLS","CHIV","CONF","CONN","DAMP","DANG","DEFS","DELU",
                "DFIX","DISP","EADP","END","EQIV","EXTI","EXYZ","FEND","FLAT","FMAP","FRAG","FREE",
                "FVAR","GRID","HFIX","HKLF","HTAB","ISOR","LATT","LAUE","LIST","L.S.","MERG","MORE",
                "MOVE","MPLA","NCSY","NEUT","OMIT","PART","PLAN","PRIG","REM","RESI","RIGU","RTAB",
                "SADI","SAME","SFAC","SHEL","SIMU","SIZE","SPEC","STIR","SUMP","SWAT","SYMM","TEMP",
                "TITL","TWIN","TWST","UNIT","WGHT","WIGL","WPDB","XNPD","ZERR" };

            vector<double> freeVariables;
            vector<string> lines, fVarsLine;
            vector<string> words;
            vector<vector<string> > atomParameterLines; // one per atom
            vector<string> symmetryLines{ "X,Y,Z" };
            vector<vector<SpaceGroupOperation> > atomPointGroups;
            vector<string> sfac;
            vector<double> adp, adp_prec;
            bool lastKeywordIsTitle = false;
            int lattice = 1;
            bool continuationLine, lineShouldBeContinued, readParameters, readFvars, cellFound = false;
            // -- read file lines
            string line, error_message;

            while (input.good())
            {
                getline(input, line);
                lines.push_back(line);
            }
            // --

            lineShouldBeContinued = false;
            readParameters = false;
            readFvars = false;
            for (auto& line : lines)
            {
                continuationLine = lineShouldBeContinued;
                lineShouldBeContinued = false;

                // All characters following '!' or '=' in an instruction line are ignored.

                if (line.find('!') != string::npos)
                    line = line.substr(0, line.find('!'));

                if (!line.empty())
                    if (line.back() == '=')
                    {
                        line = line.substr(0, line.find('='));
                        lineShouldBeContinued = true;
                    }

                //

                if (line.empty())
                    continue;

                // some format validity check (depending on white space presence on the beginning of a line

                if (continuationLine)
                {
                    if (line[0] != ' ')
                    {
                        if (lastKeywordIsTitle)
                            continue;
                        error_message = string("invalid line in shelix file/input, it should start with space character, the line is : '") +
                            line + string("'");
                        on_error::throwException(error_message, __FILE__, __LINE__);
                    }
                }

                string_utilities::split(line, words);

                if (!words.empty() && !continuationLine && line[0] == ' ')
                {
                    if (lastKeywordIsTitle)
                        continue;

                    error_message = string("invalid line in shelix file/input, it should start with space character, the line is : '") +
                        line + string("'");
                    on_error::throwException(error_message, __FILE__, __LINE__);
                }

                // ------------------------

                if (words.empty())
                    continue;


                if (continuationLine)
                {
                    if (readParameters)
                        atomParameterLines.back().insert(atomParameterLines.back().end(), words.begin(), words.end());
                    if (readFvars)
                        fVarsLine.insert(fVarsLine.end(), words.begin(), words.end());

                    readParameters = readParameters && lineShouldBeContinued;
                    readFvars = readFvars && lineShouldBeContinued;
                }
                else
                {
                    // shelx instruction line
                    readFvars = false;
                    string keywordCandidate = words[0];
                    std::transform(keywordCandidate.begin(), keywordCandidate.end(), keywordCandidate.begin(), ::toupper);

                    if (find(instructions.begin(), instructions.end(), keywordCandidate) != instructions.end())
                    {
                        string keyword = keywordCandidate;
                        readParameters = false;
                        data[keyword] = line.substr(keyword.size());


                        lastKeywordIsTitle = (keyword == "TITL");

                        // CELL, LATT, SYMM, FVAR, SFAC, END    
                        error_message.clear();
                        if (keyword == string("END"))
                            break;
                        if (keyword == string("CELL"))
                        {
                            cellFound = true;
                            if (words.size() != 8)
                                error_message = string("invalid shelix instruction file line: '") + line + string("'");
                            else
                                crystal.unitCell.set(stod(words[2]), stod(words[3]), stod(words[4]), stod(words[5]), stod(words[6]), stod(words[7]));

                        }
                        if (keyword == string("LATT"))
                        {
                            if (words.size() != 2)
                                error_message = string("invalid shelix instruction file line: '") + line + string("'");
                            else
                                lattice = stoi(words[1]);
                        }
                        if (keyword == string("SYMM"))
                            symmetryLines.push_back(line.substr(line.find('M') + 2));
                        if (keyword == string("FVAR"))
                        {
                            fVarsLine.insert(fVarsLine.end(), words.begin() + 1, words.end());
                            readFvars = true;
                        }
                        if (keyword == string("SFAC"))
                        {
                            bool element_list = true;
                            // check if it is 'SFAC E a1 b1 a2 b2 a3 b3 a4 b4 c f' f" mu r wt'
                            if (words.size() > 2)
                                if (words[2][0] == '-' || isdigit(words[2][0]) || words[2][0] == '.')
                                {
                                    element_list = false;
                                    sfac.push_back(words[1]);
                                }
                            if (element_list)
                                sfac.insert(sfac.end(), words.begin() + 1, words.end());
                        }

                        if (!error_message.empty())
                            on_error::throwException(error_message, __FILE__, __LINE__);

                    }
                    else // atom parameters line
                    {
                        atomParameterLines.push_back(words);
                        readParameters = true;
                    }
                } // if(continuationLine)else

            } // for (auto &line : lines)

            if (!cellFound)
                on_error::throwException("no information on unit cell was found in shelx file", __FILE__, __LINE__);

            // change chemical element symbols so e.g. ZN is converted to Zn
            for (auto& symbol : sfac)
                for (int i = 1; i < symbol.size(); i++)
                    symbol[i] = tolower(symbol[i]);

            // process space group 
            //A, B, C, I, F, R, H.
            map<int, char> centering{ {1,'P'}, {2,'I'}, {3,'R'},{4,'F'}, {5,'A'}, {6,'B'}, {7,'C'} };

            if (centering.find(abs(lattice)) == centering.end())
                on_error::throwException("invalid lattice specification in shelx file", __FILE__, __LINE__);

            vector<SpaceGroupOperation> spaceGroupOperations;
            for (auto const& symmOpAsString : symmetryLines)
                spaceGroupOperations.push_back(SpaceGroupOperation(symmOpAsString));

            crystal.spaceGroup.set(spaceGroupOperations, centering[abs(lattice)], true, (lattice >= 0), Vector3<CrystallographicRational>());

            // process atoms
            AtomInCrystal atom;
            vector<double> atomicParameters, precision;
            //vector<double> fvarsValues(fVarsLine.size()-1);
            vector<double> fvarsValues;
            if (!fVarsLine.empty())
            {
                fvarsValues.resize(fVarsLine.size() - 1);
                std::transform(fVarsLine.begin() + 1, fVarsLine.end(), fvarsValues.begin(), [](const string& s) {return stod(s); });
            }
            //atomParameterLines
            for (auto const& parametersLine : atomParameterLines)
            {
                if (parametersLine.size() < 6)
                {
                    error_message = "invalid atomic data specification in shelx file '";
                    for (auto& word : parametersLine)
                        error_message += word + " ";
                    error_message += "'";
                    on_error::throwException(error_message, __FILE__, __LINE__);
                }


                processParametersLine(parametersLine.begin() + 2, parametersLine.end(), fvarsValues, atomicParameters, precision);
                atom.label = parametersLine[0];
                atom.type = sfac[static_cast<int>(stoi(parametersLine[1])) - int(1)];

                atom.coordinates = Vector3d(atomicParameters[0], atomicParameters[1], atomicParameters[2]);
                atom.coordinates_precision = Vector3d(precision[0], precision[1], precision[2]);
                atom.occupancy = atomicParameters[3];
                atom.occupancy_precision = precision[3];
                atom.occupancy_sigma = 0.0;
                adp.clear();
                adp.insert(adp.end(), atomicParameters.begin() + 4, atomicParameters.end());
                adp_prec.clear();
                adp_prec.insert(adp_prec.end(), precision.begin() + 4, precision.end());
                atom.adp.clear();

                if (adp.size() == 1)
                {
                    atom.adp = adp;
                    atom.adp_precision = adp_prec;
                    atom.adp_sigma = { 0 };
                }
                else
                    if (adp.size() == 6)
                    {
                        atom.adp = { adp[0], adp[1], adp[2], adp[5], adp[4], adp[3] };
                        atom.adp_precision = { adp_prec[0], adp_prec[1], adp_prec[2], adp_prec[5], adp_prec[4], adp_prec[3] };
                        atom.adp_sigma = { 0, 0, 0, 0, 0, 0 };
                    }

                //atom.adp.clear();
                //atom.adp.insert(atom.adp.end(), atomicParameters.begin() + 4, atomicParameters.end());

                crystal.atoms.push_back(atom);

                // correct occupancy
                crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
                crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;
                crystal_structure_utilities::findAtomSymmetry(crystal, crystal.atoms.size() - 1, atomPointGroups, 0.05);
                crystal.atoms.back().multiplicity = crystal.spaceGroup.nSymmetryOperations() / atomPointGroups[0].size();

                crystal.atoms.back().occupancy *= atomPointGroups[0].size();

            }

        }


    } //namespace shelx_io
} // namespace discamb
