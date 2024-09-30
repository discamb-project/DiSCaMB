#include "discamb/IO/MATTS_BankReader.h"

#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/BasicUtilities/Exception.h"

#include "discamb/BasicChemistry/PeriodicTable.h"

//#include "Molecule/ChemicalElementData.h"

#include <fstream>
#include <set>
#include <cctype>
#include <cmath>
#include <iostream>


using namespace std;



namespace {


    vector<double> pvalSpherical =
    {
         1.0,                                                                                 0.0,
         1.0, 2.0,                                                   3.0, 4.0, 5.0, 6.0, 7.0, 0.0,
         1.0, 2.0,                                                   3.0, 4.0, 5.0, 6.0, 7.0, 0.0,
         1.0, 2.0, 1.0, 2.0, 3.0, 5.0, 5.0, 6.0, 7.0, 8.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 0.0
    };

    // s = number or chemical element symbol
    int getAtomicNumber(const string& s)
    {
        if (isdigit(s[0]))
            return stoi(s);
        string _s = s;
        if (_s[0] == '!')
            _s = _s.substr(1);
        return discamb::periodic_table::atomicNumber(_s);
    }

	// number e.g. -112.23(3)e-7 is divided into 3 components fist(second)third
	// returns true if uncertainty was defined
	bool numberWithUncertaintyFromString(const std::string &s,double &value,double &uncertainty)
	{
		string first,second,third,first_changed;
		int lastDigit,uncertaintyInteger;
		vector<string> words,words2;
		double valueBeforeChange,valueAfterChange;

		discamb::string_utilities::split(s,words,'(');

		if(words.empty()) // invalid string
            discamb::on_error::throwException(string("problem when attempting to convert string to number - invalid string: ")+s,__FILE__,__LINE__);
		if(words.size() == 1) // no uncertainty
		{
            discamb::string_utilities::convertFromString(s,value);
			uncertainty = 0;
			return false;
		}
		if(words.size() == 2)
		{
            discamb::string_utilities::split(words[1],words2,')');
			if( !(words2.size() == 2 || words2.size() == 1) )
                discamb::on_error::throwException(string("problem when attempting to convert string to number - invalid string: ")+s,__FILE__,__LINE__);
			if(!isdigit(words[0][words[0].size()-1]))
                discamb::on_error::throwException(string("problem when attempting to convert string to number - invalid string: ")+s,__FILE__,__LINE__);
			
			first = words[0];
            value = stod(words[0]);
			second = words2[0];
			if(words2.size()==2)
				third = words2[1];
			
			lastDigit = discamb::string_utilities::convertFromString<int>(string()+first[first.size()-1]);
			
			first_changed = first;
			if(lastDigit == 9)
				first_changed[first_changed.size()-1] = '8';
			else
				first_changed[first_changed.size()-1] = discamb::string_utilities::convertToString(lastDigit+1)[0];

			uncertaintyInteger = discamb::string_utilities::convertFromString<int>(second);
			valueBeforeChange = discamb::string_utilities::convertFromString<double>(first+third);
			valueAfterChange = discamb::string_utilities::convertFromString<double>(first_changed+third);
			uncertainty = fabs(valueAfterChange-valueBeforeChange)*uncertaintyInteger;
			return true;
		}
		return false;
	}

    void addSphericalTypeLines(std::ostream &out)
    {
        out << "\n";
        for (int atomicNumber = 1; atomicNumber <= 36; atomicNumber++)
        {
            double kappa;
            atomicNumber == 1 ? kappa = 1.16 : kappa = 1.0;
            string symbol = discamb::periodic_table::symbol(atomicNumber);
            out << "ENTRY\n"
                << "  ID\n"
                << "    " << symbol << "000\n"
                << "  COMMENT\n"
                << "    spherical " << symbol << "\n"
                << "  ATOM DESCRIPTORS\n"
                << "#  central atom\n"
                << "    " << symbol << "   CONNECTED_TO  *            PLANARITY  *    PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *\n"
                << "  LOCAL COORDINATE SYSTEM\n"
                << "    X any_orthogonal Y any_orthogonal R\n"
                << "  SYMMETRY\n"
                << "    sph\n"
                << "  PARAMETER MODIFICATION DATE\n"
                << "    Thu May  2 17:26 : 32 2019\n"
                << "  ENTRY CREATION DATE\n"
                << "    Thu May  2 17:26:32 2019\n"
                << "  DEFINITION MODIFICATION DATE\n"
                << "    Thu May  2 17:26:32 2019\n"
                << "  MULTIPOLE MODEL PARAMETERS\n"
                << "    PVAL " << fixed << pvalSpherical[atomicNumber-1] << " KAPPA " << fixed << setprecision(3) << kappa << " KPRIM 1.000\n"
                << "    PLMS 0 0 0.000\n\n";
        }
    }
}

namespace discamb {

    using namespace string_utilities;

    // initialize static members

    map<string, MATTS_BankReader::Section> MATTS_BankReader::mKeywordToSectionMap = MATTS_BankReader::createKeywordToSectionMap();
    map<MATTS_BankReader::Section, std::string> MATTS_BankReader::mSectionToKeywordMap = MATTS_BankReader::createSectionToKeywordMap();

    //


    MATTS_BankReader::MATTS_BankReader()
    {
    }


    MATTS_BankReader::~MATTS_BankReader()
    {
    }


    void MATTS_BankReader::read(
        const string &filename,
        vector<AtomType> &types,
        vector<AtomTypeHC_Parameters> &parameters,
        BankSettings &bankSettings,
        bool addSphericalTypes,
        RequiredSections requiredSections)
        const
    {
        ifstream input(filename.c_str());

        types.clear();
        parameters.clear();

        if (!input.good())
            on_error::throwException(string("Problem with accessing/openning file '") + filename + string("' for reading UBDB type bank"), __FILE__, __LINE__);

        //read(input, types, parameters, bankSettings);

        stringstream ss;

        if (addSphericalTypes)
        {
            string line;
            while (input.good())
            {
                getline(input, line);
                ss << line << "\n";
            }
            addSphericalTypeLines(ss);
            read(ss, types, parameters, bankSettings, false, requiredSections);
        }
        else
            read(input, types, parameters, bankSettings, false, requiredSections);


        input.close();

    }

    void MATTS_BankReader::read(
        const std::string &fName,
        std::vector<AtomType> &types,
        std::vector<AtomTypeHC_Parameters> &parameters,
        std::vector<std::string> &lines,
        std::vector<std::pair<int, int> > &parameter_lines,
        BankSettings &bankSettings,
        bool addSphericalTypes,
        RequiredSections requiredSections)
        const
    {
        ifstream in(fName);
        types.clear();
        parameters.clear();

        if (!in.good())
            on_error::throwException(string("Problem with accessing/openning file '") + fName + string("' for reading UBDB type bank"), __FILE__, __LINE__);

        stringstream ss;

        if (addSphericalTypes)
        {
            string line;
            while (in.good())
            {
                getline(in, line);
                ss << line << "\n";
            }
            addSphericalTypeLines(ss);
            
            read(ss, types, parameters, lines, parameter_lines, bankSettings, false, requiredSections);
        }
        else
            read(in, types, parameters, lines, parameter_lines, bankSettings, false, requiredSections);

        in.close();

    }

    std::string MATTS_BankReader::sectionToString(
        Section s)
    {
        auto m = createSectionToKeywordMap();
        return m[s];
    }

    void MATTS_BankReader::error_if_required(
        Section section, 
        RequiredSections sections,
        int lineIdx)
    {
        if (sections.undefined)
        {
            sections.undefined = false;
            sections.sections = 
                vector<Section>({ Section::ID, Section::COMMENT, Section::ATOM_DESCRIPTORS, Section::LOCAL_COORDINATE_SYSTEM, 
                    Section::SYMMETRY, Section::MULTIPOLE_MODEL_PARAMETERS });
            error_if_required(section, sections, lineIdx);
            return;
        }
        if (find(sections.sections.begin(), sections.sections.end(), section) != sections.sections.end())
            on_error::throwException("problem when processing atom type bank - expected section " + sectionToString(section) +
                " near line " + to_string(lineIdx), __FILE__, __LINE__);
    }

    void MATTS_BankReader::read(
        std::istream &_input,
        std::vector<AtomType> &types,
        std::vector<AtomTypeHC_Parameters> &parameters,
        std::vector<std::string> &lines,
        std::vector<std::pair<int, int> > &parameter_lines,
        BankSettings &bankSettings,
		bool addSphericalTypes,
        RequiredSections requiredSections)
        const
    {
		
		if (addSphericalTypes)
		{
			stringstream ss;
			string line;
			while (_input.good())
			{
				getline(_input, line);
				ss << line << "\n";
			}
			addSphericalTypeLines(ss);

			read(ss, types, parameters, lines, parameter_lines, bankSettings, false, requiredSections);
			return;
		}

		std::istream& input = _input;

		//-------
        types.clear();
        parameters.clear();
        lines.clear();
        parameter_lines.clear();
        int first, last;

        string line;
        vector<string> words;
        map<string, set<int> > element_groups;
        int lineNumber;

        lineNumber = 0;

        while (portable_getline(input, line))
        {
            lines.push_back(line);

            lineNumber++;

            line = trim_left(line);

            if (line.empty())  // skip empty line
                continue;

            if (line[0] == '#') // skip comment line
                continue;

            split(line, words, CharacterType::WHITE_SPACE);
            if (words.size() == 1)
            {
                if (words[0] == string("ELEMENTS_GROUP"))
                {
                    std::set<int> group;
                    portable_getline(input, line);
                    lines.push_back(line);
                    lineNumber++;
                    string_utilities::split(line, words);
                    if (words.size() != 2)
                        on_error::throwException("invlaid format of ELEMENTS_GROUP in atomic types specification", __FILE__, __LINE__);
                    stringToAtomicNumbersSet(words[1], group, element_groups);
                    element_groups[words[0]] = group;
                }
                if (words[0] == string("ENTRY"))
                {
                    types.resize(types.size() + 1);
                    parameters.resize(parameters.size() + 1);
                    bool read_1st_multipoles_data_line;
                    readType(input, types.back(), lineNumber, lines, read_1st_multipoles_data_line, element_groups, requiredSections);
                    first = lineNumber;
//                  types.back().number = types.size() - 1;
                    if (read_1st_multipoles_data_line)
                    {
                        if (!readMultipoleParameters(input, parameters.back(), lineNumber, lines))
                        {
                            //cout << "no multipole parameters " << __LINE__ << "\n";
                            error_if_required(Section::MULTIPOLE_MODEL_PARAMETERS, requiredSections, lineNumber);
                        }
                        //else
                        //    cout << "multipole parameters " << __LINE__ << "\n";
                    }
                    last = lineNumber - 1;
                    if (string_utilities::trim(lines.back()).empty())
                        last--;
                    parameter_lines.push_back({ first, last });
                }
                if (words[0] == string("SETTINGS"))
                    readSettings(input, bankSettings, lineNumber, lines);
            }

            
        }

    }



    void MATTS_BankReader::read(
        istream &_input,
        vector<AtomType> &types,
        vector<AtomTypeHC_Parameters> &parameters,
        BankSettings &bankSettings,
		bool addSphericalTypes,
        RequiredSections requiredSections)
        const
    {

		if (addSphericalTypes)
		{
			stringstream ss;
			string line;
			while (_input.good())
			{
				getline(_input, line);
				ss << line << "\n";
			}
			addSphericalTypeLines(ss);

			read(ss, types, parameters, bankSettings, false, requiredSections);
			return;
		}

		std::istream& input = _input;
		 
		//---------------

        types.clear();
        parameters.clear();

        string line;
        vector<string> words;
        map<string, set<int> > element_groups;
        int lineNumber;

        lineNumber = 0;

        while (portable_getline(input, line))
        {
            lineNumber++;

            line = trim_left(line);

            if (line.empty())  // skip empty line
                continue;

            if (line[0] == '#') // skip comment line
                continue;

            split(line, words, CharacterType::WHITE_SPACE);
            if (words.size() == 1)
            {
                if (words[0] == string("ELEMENTS_GROUP"))
                {
                    std::set<int> group;
                    portable_getline(input, line);
                    lineNumber++;
                    string_utilities::split(line, words);
                    if (words.size() != 2)
                        on_error::throwException("invlaid format of ELEMENTS_GROUP in atomic types specification", __FILE__, __LINE__);
                    stringToAtomicNumbersSet(words[1], group, element_groups);
                    element_groups[words[0]] = group;
                }

                if (words[0] == string("ENTRY"))
                {
                    types.resize(types.size() + 1);
                    parameters.resize(parameters.size() + 1);
                    bool read_1st_multipoles_data_line;
                    readType(input, types.back(), lineNumber, read_1st_multipoles_data_line, element_groups, requiredSections);
//                    types.back().number = types.size() - 1;
                    if(read_1st_multipoles_data_line)
                        if (!readMultipoleParameters(input, parameters.back(), lineNumber))
                            error_if_required(Section::MULTIPOLE_MODEL_PARAMETERS, requiredSections, lineNumber);

                }
                if (words[0] == string("SETTINGS"))
                {
                    vector<string> settingsLines;
                    readSettings(input, bankSettings, lineNumber, settingsLines);
                }

            }
        }
    }



    bool MATTS_BankReader::readMultipoleParameters(
        istream &in,
        AtomTypeHC_Parameters &parameters,
        int &lastReadLineNumber,
        vector<string> &lines)
    const
    {
        bool emptyLine = false;
        string line;
        vector<string> words;
        int nWords, nPlm;

        // works on the first line e.g.: PVAL 0.806 KAPPA 1.167 KPRIM 1.248 SIGPV 0.017
        portable_getline(in, line);
        lines.push_back(line);
        lastReadLineNumber++;
        split(line, words, CharacterType::WHITE_SPACE);

        if (words.empty())
            return false;

        nWords = words.size();

        if (nWords == 6 || nWords == 8)
        {
            if (nWords == 6)
                numberWithUncertaintyFromString(words[1], parameters.p_val, parameters.p_val_sigma);
            else
            {
                convertFromString(words[1], parameters.p_val);
                convertFromString(words[7], parameters.p_val_sigma);
            }

            numberWithUncertaintyFromString(words[3], parameters.kappa, parameters.kappa_sigma);
            numberWithUncertaintyFromString(words[5], parameters.kappa_prime, parameters.kappa_prime_sigma);
        }
        else
            on_error::throwException(string("invalid line in aspherical atom data bank : ") + line, __FILE__, __LINE__);

        // works on lines with Plms e.g.: PLMS 1 0 0.256 PLMS 2 0 0.111

        do
        {
            portable_getline(in, line);
            lines.push_back(line);
            lastReadLineNumber++;
            line = trim(line);
            split(line, words, CharacterType::WHITE_SPACE);
            nPlm = words.size() / 4;
            for (int i = 0; i < nPlm; i++)
            {
                parameters.p_lms.resize(parameters.p_lms.size() + 1);
                parameters.p_lms_sigma.resize(parameters.p_lms_sigma.size() + 1);
                parameters.p_lm_indices.resize(parameters.p_lm_indices.size() + 1);
                parameters.p_lm_indices.back() = make_pair(convertFromString<int>(words[4 * i + 1]), convertFromString<int>(words[4 * i + 2]));
                numberWithUncertaintyFromString(words[4 * i + 3], parameters.p_lms.back(), parameters.p_lms_sigma.back());
            }
        } while (!line.empty());
        return true;
    }


    bool MATTS_BankReader::readMultipoleParameters(
        istream &in,
        AtomTypeHC_Parameters &parameters,
        int &lastReadLineNumber)
        const
    {
        bool emptyLine = false;
        string line;
        vector<string> words;
        int nWords, nPlm;

        // works on the first line e.g.: PVAL 0.806 KAPPA 1.167 KPRIM 1.248 SIGPV 0.017
        portable_getline(in, line);
        lastReadLineNumber++;
        split(line, words, CharacterType::WHITE_SPACE);
        if (words.empty())
            return false;

        nWords = words.size();

        if (nWords == 6 || nWords == 8)
        {
            if (nWords == 6)
                numberWithUncertaintyFromString(words[1], parameters.p_val, parameters.p_val_sigma);
            else
            {
                convertFromString(words[1], parameters.p_val);
                convertFromString(words[7], parameters.p_val_sigma);
            }

            numberWithUncertaintyFromString(words[3], parameters.kappa, parameters.kappa_sigma);
            numberWithUncertaintyFromString(words[5], parameters.kappa_prime, parameters.kappa_prime_sigma);
        }
        else
            on_error::throwException(string("invalid line in aspherical atom data bank : ") + line, __FILE__, __LINE__);

        // works on lines with Plms e.g.: PLMS 1 0 0.256 PLMS 2 0 0.111

        do
        {
            portable_getline(in, line);
            lastReadLineNumber++;
            line = trim(line);
            split(line, words, CharacterType::WHITE_SPACE);
            nPlm = words.size() / 4;
            for (int i = 0; i < nPlm; i++)
            {
                parameters.p_lms.resize(parameters.p_lms.size() + 1);
                parameters.p_lms_sigma.resize(parameters.p_lms_sigma.size() + 1);
                parameters.p_lm_indices.resize(parameters.p_lm_indices.size() + 1);
                parameters.p_lm_indices.back() = make_pair(convertFromString<int>(words[4 * i + 1]), convertFromString<int>(words[4 * i + 2]));
                numberWithUncertaintyFromString(words[4 * i + 3], parameters.p_lms.back(), parameters.p_lms_sigma.back());
            }
        } while (!line.empty());
        return true;
    }


    void MATTS_BankReader::readType(
        istream &in,
        AtomType &type,
        int &lastReadLineNumber,
        bool& readMultipolesLine,
        const std::map<std::string, std::set<int> >& namedElementSets,
        RequiredSections requiredSections)
        const
    {
        readMultipolesLine = false;
        bool readMore;
        string line;
        vector<string> sectionLines, words;
        Section section, previousSection, currentSection;
        //vector<bool> sectionEncountered(4, false);
        map<MATTS_BankReader::Section, bool> requiredSectionsEncounter {
            {Section::COMMENT, false},{Section::ATOM_DESCRIPTORS,false}, {Section::LOCAL_COORDINATE_SYSTEM, false}, {Section::SYMMETRY,false} };
        if (!requiredSections.undefined)
        {
            requiredSectionsEncounter.clear();
            for (auto s : requiredSections.sections)
                if (s != Section::MULTIPOLE_MODEL_PARAMETERS)
                    requiredSectionsEncounter[s] = false;
        }

        readMore = true;
        previousSection = Section::NONE;
        currentSection = Section::NONE;
        Section endSection = Section::MULTIPOLE_MODEL_PARAMETERS;
        if (!requiredSections.undefined)
            if (find(requiredSections.sections.begin(), requiredSections.sections.end(), Section::MULTIPOLE_MODEL_PARAMETERS) == requiredSections.sections.end())
                endSection = Section::BLANK_LINE;
        

        while (currentSection != endSection && currentSection != Section::MULTIPOLE_MODEL_PARAMETERS && currentSection != Section::END)
        {
            portable_getline(in, line);

            string_utilities::split(line, words);
            if (words.empty() && endSection == Section::BLANK_LINE)
            {
                currentSection = Section::BLANK_LINE;
            }
            else
            {
                if (!in.good())
                    on_error::throwException("reached end of UBDB bank file without finding expected MULTIPOLE MODEL PARAMETERS section", __FILE__, __LINE__);

                lastReadLineNumber++;

                line = trim(line);

                if (line.empty())
                    on_error::throwException(string("unexpected empty line in UBDB bank at line ") + convertToString(lastReadLineNumber), __FILE__, __LINE__);

                if (line[0] == '#') // skip comment line
                    continue;

                if (isNewSection(line, section))
                {
                    if (section == Section::MULTIPOLE_MODEL_PARAMETERS)
                        readMultipolesLine = true;
                    previousSection = currentSection;
                    currentSection = section;
                    if (requiredSectionsEncounter.find(section) != requiredSectionsEncounter.end())
                        requiredSectionsEncounter[section] = true;
                    //if(section < 4)
                       // sectionEncountered[section] = true;

                    if (previousSection != Section::NONE)
                        setSection(sectionLines, previousSection, type, lastReadLineNumber - sectionLines.size(), namedElementSets);

                    sectionLines.clear();
                }
                else
                {
                    if (section == Section::NONE)
                        on_error::throwException(string("invalid format of UBDB bank, line ") + convertToString(lastReadLineNumber++), __FILE__, __LINE__);
                    sectionLines.push_back(line);
                }
            }
        }

        for(auto const &sectionEncounter: requiredSectionsEncounter)
            if (!sectionEncounter.second)
                on_error::throwException(string("invalid atom type specification in UBDB data bank - missing section ") + mSectionToKeywordMap[sectionEncounter.first] + string(" for entry near line ") + convertToString(lastReadLineNumber), __FILE__, __LINE__);

    }

    void MATTS_BankReader::readType(
        istream &in,
        AtomType &type,
        int &lastReadLineNumber,
        vector<string> &lines,
        bool& readMultipolesLine,
        const std::map<std::string, std::set<int> >& namedElementSets,
        RequiredSections requiredSections)
        const
    {
        readMultipolesLine = false;
        bool readMore;
        string line;
        vector<string> sectionLines;
        Section section, previousSection, currentSection;
        //vector<int> sectionEncountered(6, false);
        map<Section, bool> sectionPresence;

        for (auto &x : mSectionToKeywordMap)
            sectionPresence[x.first] = false;

        readMore = true;
        previousSection = Section::NONE;
        currentSection = Section::NONE;

        Section endSection = Section::MULTIPOLE_MODEL_PARAMETERS;
        if (!requiredSections.undefined)
            if (find(requiredSections.sections.begin(), requiredSections.sections.end(), Section::MULTIPOLE_MODEL_PARAMETERS) == requiredSections.sections.end())
                endSection = Section::BLANK_LINE;

        vector<string> words;
        while (currentSection != endSection && currentSection != Section::MULTIPOLE_MODEL_PARAMETERS)
        {
            portable_getline(in, line);
            lines.push_back(line);
            string_utilities::split(line, words);
            if (words.empty() && endSection == Section::BLANK_LINE)
            {
                currentSection = Section::BLANK_LINE;
            }
            else
            {

                if (!in.good())
                    on_error::throwException("reached end of UBDB bank file without finding expected MULTIPOLE MODEL PARAMETERS section", __FILE__, __LINE__);

                lastReadLineNumber++;

                line = trim(line);

                if (line.empty())
                    on_error::throwException(string("unexpected empty line in UBDB bank at line ") + convertToString(lastReadLineNumber), __FILE__, __LINE__);

                if (line[0] == '#') // skip comment line
                    continue;

                if (isNewSection(line, section))
                {
                    if (section == Section::MULTIPOLE_MODEL_PARAMETERS)
                        readMultipolesLine = true;

                    previousSection = currentSection;
                    currentSection = section;
                    sectionPresence[section] = true;
                    //sectionEncountered[section] = true;

                    if (previousSection != Section::NONE)
                        setSection(sectionLines, previousSection, type, lastReadLineNumber - sectionLines.size(), namedElementSets);

                    sectionLines.clear();
                }
                else
                {
                    if (section == Section::NONE)
                        on_error::throwException(string("invalid format of UBDB bank, line ") + convertToString(lastReadLineNumber++), __FILE__, __LINE__);
                    sectionLines.push_back(line);
                }
            }
        }

        vector<Section> requiredSectionsList{ Section::COMMENT, Section::ATOM_DESCRIPTORS, Section::LOCAL_COORDINATE_SYSTEM, Section::SYMMETRY };

        if (!requiredSections.undefined)
        {
            requiredSectionsList.clear();
            for (auto s : requiredSections.sections)
                if (s != Section::MULTIPOLE_MODEL_PARAMETERS)
                    requiredSectionsList.push_back(s);
        }


        vector<Section> missingRequiredSections;
        for (auto requiredSection : requiredSectionsList)
            if (!sectionPresence[requiredSection])
                missingRequiredSections.push_back(requiredSection);

        if( !missingRequiredSections.empty())
            on_error::throwException(string("invalid atom type specification in UBDB data bank - missing section ") + mSectionToKeywordMap[missingRequiredSections[0]] + string(" for entry near line ") + convertToString(lastReadLineNumber), __FILE__, __LINE__);

    }



    bool MATTS_BankReader::isNewSection(
        const std::string &line,
        Section &section)
        const
    {
        map<string, Section>::const_iterator it = mKeywordToSectionMap.find(line);

        if (it == mKeywordToSectionMap.end())
            return false;

        section = it->second;

        return true;
    }

    void MATTS_BankReader::setChirality(
        const std::string &line,
        AtomType &type)
    {
        vector<string> words;
        map<string, int> labelToIndex;
        for (int i = 0; i < type.atoms.size(); i++)
            labelToIndex[type.atoms[i].label] = i;

        split(line, words, CharacterType::WHITE_SPACE);
        if (words.size() != 3)
            on_error::throwException(string("invalid number of labels defining chirality in UBDB entry: '") + line + string("'"), __FILE__, __LINE__);
        type.chirality.clear();
        for (auto &word : words)
            if (labelToIndex.find(word) == labelToIndex.end())
                on_error::throwException(string("problem when processing chirality line in UBDB :'") + line + string("'"), __FILE__, __LINE__);
            else
                type.chirality.push_back(labelToIndex[word]);
    }

    void  MATTS_BankReader::setSection(
        const vector<string> &lines,
        Section section,
        AtomType &type,
        int firstLineNumber,
        const std::map<std::string, std::set<int> >& namedElementSets)
    {
        vector<string> words;
        string errorMessage;

        //if( section!= COMMENT)
            //if(lines.empty())
                //on_error::throwException(string("invalid section in UBDB bank, line ") + convertToString(firstLineNumber-1), __FILE__ , __LINE__ );

        switch (section)
        {
        case Section::ID:
            type.id = lines[0];
            break;
        case Section::COMMENT:
            type.commentLines = lines;
            split(lines[0], words, CharacterType::WHITE_SPACE);
            if(type.id.empty())
                type.id = words[0];
            break;

        case Section::SYMMETRY:
            type.symmetry = lines[0];
            break;
        case Section::PARAMETER_MODIFCATION_DATE:
            type.parameterModificationDate = lines[0];
            break;
        case Section::ENTRY_CREATION_DATE:
            type.entryCreationDate = lines[0];
            break;
        case Section::DEFINITION_MODIFICATION_DATE:
            type.definitionModificationDate = lines[0];
            break;
        case Section::ATOM_DESCRIPTORS:
            setAtomDescriptors(lines, type, firstLineNumber, namedElementSets);
            break;
        case Section::CHIRALITY:
            setChirality(lines[0], type);
            break;
        case Section::BLANK_LINE:
            break;
        case Section::LOCAL_COORDINATE_SYSTEM:
            try
            {
                type.setLocalCoordinateSystem(lines[0]);
            }
            catch (Exception &e)
            {
                string errorMessage = string("invalid specification of local coordinate system in UBDB bank, line : ") + convertToString(firstLineNumber) +
                    string(". Details: ") + e.errorMessage();
                on_error::throwException(errorMessage, __FILE__, __LINE__);
            }
            break;
        default:
            on_error::throwException("REPORT BUG", __FILE__, __LINE__);
        }

    }

    std::map<std::string, MATTS_BankReader::Section> MATTS_BankReader::createKeywordToSectionMap()
    {
        std::map<std::string, Section> result;

        result["COMMENT"] = Section::COMMENT;
        result["CHIRALITY"] = Section::CHIRALITY;
        result["SYMMETRY"] = Section::SYMMETRY;
        result["ATOM DESCRIPTORS"] = Section::ATOM_DESCRIPTORS;
        result["LOCAL COORDINATE SYSTEM"] = Section::LOCAL_COORDINATE_SYSTEM;
        result["MULTIPOLE MODEL PARAMETERS"] = Section::MULTIPOLE_MODEL_PARAMETERS;
        result["MULTIPOLE MODEL PARAMTERS"] = Section::MULTIPOLE_MODEL_PARAMETERS;
        //result["DATE"] = DATE;
        result["ID"] = Section::ID;
        result["PARAMETER MODIFCATION DATE"] = Section::PARAMETER_MODIFCATION_DATE;
        result["ENTRY CREATION DATE"] = Section::ENTRY_CREATION_DATE;
        result["DEFINITION MODIFICATION DATE"] = Section::DEFINITION_MODIFICATION_DATE;
        result["END"] = Section::END;
        return result;
    }

    std::map<MATTS_BankReader::Section, std::string> MATTS_BankReader::createSectionToKeywordMap()
    {
        std::map<Section, std::string> result;

        result[Section::COMMENT] = "COMMENT";
        result[Section::SYMMETRY] = "SYMMETRY";
        result[Section::ATOM_DESCRIPTORS] = "ATOM DESCRIPTORS";
        result[Section::LOCAL_COORDINATE_SYSTEM] = "LOCAL COORDINATE SYSTEM";
        result[Section::MULTIPOLE_MODEL_PARAMETERS] = "MULTIPOLE MODEL PARAMETERS";
        result[Section::CHIRALITY] = "CHIRALITY";
        result[Section::ID] = "ID";
        result[Section::PARAMETER_MODIFCATION_DATE] = "PARAMETER MODIFCATION DATE";
        result[Section::ENTRY_CREATION_DATE] = "ENTRY CREATION DATE";
        result[Section::DEFINITION_MODIFICATION_DATE] = "DEFINITION MODIFICATION DATE";
        result[Section::END] = "END";
    
        return result;
    }

    void MATTS_BankReader::setAtomDescriptors(
        const std::vector<std::string> &lines,
        AtomType &type,
        int firstLineNumber,
        const std::map<std::string, std::set<int> >& namedElementSets)
    {
        /*
      ATOM DESCRIPTORS
    #  central atom
        C1    CONNECTED_TO  O2,O3,C4      PLANARITY  +    RING  -
    #  1-st neighbors
        O2    CONNECTED_TO  C1            PLANARITY  *    RING  -
        O3    CONNECTED_TO  C1,C          PLANARITY  *    RING  -
        C4    CONNECTED_TO  C1,*          PLANARITY  *    RING  * IN_3_MEMBER_RING * IN_4_MEMBER_RING  *

        IN_3_MEMBER_RING * IN_4_MEMBER_RING  *


        */

        int atomIndex, nAtoms, neighborIndex, nNeighbors, neighborAtomIndex, ringIndex, nRings;
        vector<string> words, labels;
        vector<string>::iterator vecStrIterator;
        vector<vector<string> > linesWords;
        std::set<string> ringLabels;

        nAtoms = lines.size();

        type.atoms.resize(nAtoms);
        labels.resize(nAtoms);
        linesWords.resize(nAtoms);
        type.connectivity.resize(nAtoms);

        bool has3and4memberRingsInfo = false;

        // split describing lines into words

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            split(lines[atomIndex], linesWords[atomIndex], CharacterType::WHITE_SPACE);

            if (linesWords[atomIndex].size() != 7 && linesWords[atomIndex].size() != 11)
                on_error::throwException(string("invalid line in UBDB bank - expected 7 or 11 'words' description of atom, line ") + convertToString(firstLineNumber + atomIndex) + string(";\n") + lines[atomIndex], __FILE__, __LINE__);
            
        }

        // assigns atomic numbers and labels to atoms, fill labels vector
        // and processes planarity info

        //map<string, set<int> > atom2atomicNumber;

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            // atomic numbers, labels and labels vector
            labelToAtomicNumber(linesWords[atomIndex][0], type.atoms[atomIndex].atomic_number_range, namedElementSets, type.atoms[atomIndex].anyAtomicNumber);
            
            //type.atoms[atomIndex].atomic_number_range = labelToAtomicNumber(linesWords[atomIndex][0]);
            string label = linesWords[atomIndex][0].substr(0, linesWords[atomIndex][0].find(','));
           // if (type.atoms[atomIndex].atomic_number_range.size() == 1)
             //   atom2atomicNumber[label] = type.atoms[atomIndex].atomic_number_range;
            //else
            //    atom2atomicNumber[label] = 0;
            //type.atoms[atomIndex].label = linesWords[atomIndex][0];
            type.atoms[atomIndex].label = label;
            //labels[atomIndex] = linesWords[atomIndex][0];
            labels[atomIndex] = label;// linesWords[atomIndex][0];

            // planarity info

            switch (linesWords[atomIndex][4][0])
            {
            case '+':
                type.atoms[atomIndex].planar = Tribool::True;
                break;
            case '-':
                type.atoms[atomIndex].planar = Tribool::False;
                break;
            case '*':
                type.atoms[atomIndex].planar = Tribool::Undefined;
                break;
            default:
                on_error::throwException(string("wrong specification of planarity in UBDB file at line ") + convertToString(firstLineNumber + atomIndex), __FILE__, __LINE__);
            }

        }

        // get connectivity and formula

        set<int> atomicNumbersSet;
        bool anyAtomicNumber;
        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            split(linesWords[atomIndex][2], words, ',');
            type.atoms[atomIndex].fixedNumberOfNeighbors = true;
            type.atoms[atomIndex].nNeighbours = words.size();
            if (linesWords[atomIndex][2][0] == '-')
                type.atoms[atomIndex].nNeighbours = 0;
            
            //type.atoms[atomIndex].neighborsAtomicNumbersUniquelyDefined = true;
            //type.atoms[atomIndex].allNonX_NeighborsNamed = true;

            for (int i = 0; i < words.size(); i++)
            {
                // if the neighbor is among the labeled atoms update connectivity list

                vecStrIterator = find(labels.begin(), labels.end(), words[i]);

                if (vecStrIterator != labels.end()) // yes it is
                {
                    neighborIndex = std::distance(labels.begin(), vecStrIterator);
                    type.connectivity[atomIndex].push_back(neighborIndex);
                }
                //else // non named element or element range
                //    if (words[i][0] != '*' && words[i][0] != '-' && words[i][0] != 'X')
                //        type.atoms[atomIndex].allNonX_NeighborsNamed = false;

                // update neighbors atomic numbers list (kind of formula)
                if (words[i][0] != '*' && words[i][0] != '-' && words[i][0] != 'X')
                {
                    labelToAtomicNumber(words[i], atomicNumbersSet, namedElementSets, anyAtomicNumber);
                    if (atomicNumbersSet.size() == 1)
                        type.atoms[atomIndex].neighborsAtomicNumbers.insert(*atomicNumbersSet.begin());
                        //type.atoms[atomIndex].neighborsAtomicNumbers.push_back(*atomicNumbersSet.begin());
                    else
                        type.atoms[atomIndex].neighborsAtomicNumberRanges.push_back(atomicNumbersSet);
                    //if (atom2atomicNumber[words[i]].size()==1)
                    //    type.atoms[atomIndex].neighborsAtomicNumbers.push_back(*atom2atomicNumber[words[i]].begin());
                    //else
                    //    type.atoms[atomIndex].neighborsAtomicNumberRanges.push_back(atom2atomicNumber[words[i]]);
                }

                //if (words[i][0] == '*' || words[i][0] != 'X')
                //    type.atoms[atomIndex].neighborsAtomicNumbersUniquelyDefined = false;
                
                if (words[i][0] == '*') // since it was counted as a neighbour before
                {
                    type.atoms[atomIndex].nNeighbours--;
                    type.atoms[atomIndex].fixedNumberOfNeighbors = false;
                }
            }
            // sort 1-sth neighbors atomic numbers (make canonical form of 'formula')
            //sort(type.atoms[atomIndex].neighborsAtomicNumbers.begin(), type.atoms[atomIndex].neighborsAtomicNumbers.end());
            //type.atoms[atomIndex].neighborsAtomicNumbersUniquelyDefined = true;
            //if (!type.atoms[atomIndex].neighborsAtomicNumberRanges.empty() || 
            //    !type.atoms[atomIndex].fixedNumberOfNeighbors ||
            //    type.atoms[atomIndex].nNeighbours != type.atoms[atomIndex].neighborsAtomicNumbers.size())
            //    type.atoms[atomIndex].neighborsAtomicNumbersUniquelyDefined = false;
            
        }

        // verify if connectivity matrix is symmetric

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            nNeighbors = type.connectivity[atomIndex].size();

            for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
            {
                neighborAtomIndex = type.connectivity[atomIndex][neighborIndex];

                if (find(type.connectivity[neighborAtomIndex].begin(), type.connectivity[neighborAtomIndex].end(), atomIndex) == type.connectivity[neighborAtomIndex].end())
                    on_error::throwException(string("invalid entry in UBDB file: non-symmetric connectivity matrix for type ") + type.id, __FILE__, __LINE__);
            }
        }

        // collect information (per atom) on rings

        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            has3and4memberRingsInfo = (linesWords[atomIndex].size() == 11);

            if (has3and4memberRingsInfo)
                type.atoms[atomIndex].ringInfo.set(linesWords[atomIndex][6], linesWords[atomIndex][8], linesWords[atomIndex][10]);
            else
                type.atoms[atomIndex].ringInfo.set(linesWords[atomIndex][6], "*", "*");
        }

        // extract info on labeled rings


        string ringLabel;
        type.ringLabels.clear();
        type.ringSizes.clear();
        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
        {
            nRings = type.atoms[atomIndex].ringInfo.labeledContainingRings.size();
            for (ringIndex = 0; ringIndex < nRings; ringIndex++)
            {
                ringLabel = type.atoms[atomIndex].ringInfo.labeledContainingRings[ringIndex].second;
                if (find(type.ringLabels.begin(), type.ringLabels.end(), ringLabel) == type.ringLabels.end())
                {
                    type.ringLabels.push_back(ringLabel);
                    type.ringSizes.push_back(type.atoms[atomIndex].ringInfo.labeledContainingRings[ringIndex].first);
                }
            }
        }
    }

    void MATTS_BankReader::stringToAtomicNumbersSet(
        const std::string& s,
        std::set<int>& atomicNumbers,
        const std::map<std::string, std::set<int> >& namedSets)
    {
        atomicNumbers.clear();
        // format (1) C1
        // format (2) !H
        // format (3) Label1,atomic_numbers
        // format (1) i.e. SymbolNumber
        // can include namedSets, e.g. Hal2
        // and there should be set labeled Hal
        // the seam for format (2), it can be also !Hal
        // atomic_numbers in format(3) can be:
        // - single symbol/atomic number Atom1,C Atom1,6 Atom1,Hal
        // - list Atom1,C,H or Atom1,1,6 also e.g. Atom1,-1,4
        // - range Atom1,B-F or Atom1,4-9 also e.g. Atom1,!4-9
        // - combination of ranges e.g. Atom1,B-F,Si-Cl

        set<int> toInclude, toRemove;


        vector<string> words;
        set<int> theSet;
        string_utilities::split(s, words, ',');
        for (auto& word : words)
        {
            theSet.clear();
            bool elementsToRemove = word[0] == '!';
            if (elementsToRemove)
                word = word.substr(1);

            bool predefinedSet = namedSets.find(word) != namedSets.end();

            if (predefinedSet)
                theSet = namedSets.find(word)->second;
            else
            {
                auto position = word.find('-');
                if (position != string::npos)// range
                {
                    int start = getAtomicNumber(word.substr(0, position));
                    int stop = getAtomicNumber(word.substr(position + 1));
                    for (int i = start; i <= stop; i++)
                        theSet.insert(i);
                }
                else
                {
                    int z = getAtomicNumber(word);
                    if(z>0)
                        theSet.insert(getAtomicNumber(word));
                }
            }
            if (elementsToRemove)
                toRemove.insert(theSet.begin(), theSet.end());
            else
                toInclude.insert(theSet.begin(), theSet.end());

        }

        if (toInclude.empty())
            for (int i = 1; i <= 120; i++)
                toInclude.insert(i);

        for (auto z : toInclude)
            if (toRemove.find(z) == toRemove.end())
                atomicNumbers.insert(z);

    }


    void MATTS_BankReader::labelToAtomicNumber(
        const std::string &label,
        set<int> &atomicNumbers,
        const map<string, set<int> > &namedSets,
        bool &anyAtomicNumber)
    {
        anyAtomicNumber = false;
        atomicNumbers.clear();
        // format (1) C1
        // format (2) !H
        // format (3) Label1,atomic_numbers
        // format (1) i.e. SymbolNumber
        // can include namedSets, e.g. Hal2
        // and there should be set labeled Hal
        // the seam for format (2), it can be also !Hal
        // atomic_numbers in format(3) can be:
        // - single symbol/atomic number Atom1,C Atom1,6 Atom1,Hal
        // - list Atom1,C,H or Atom1,1,6 also e.g. Atom1,-1,4
        // - range Atom1,B-F or Atom1,4-9 also e.g. Atom1,!4-9
        // - combination of ranges e.g. Atom1,B-F,Si-Cl
              

        auto commaPosition = label.find(',');
        string definition;
        if (commaPosition == string::npos)
        {
            int i, labelSize = label.size();
            for (i = 0; i < labelSize; i++)
                if (!isdigit(label[i]))
                    definition += label[i];
        }
        else
            definition = label.substr(0, commaPosition - 1);

        if (definition == string("X"))
            anyAtomicNumber = true;
        else
            stringToAtomicNumbersSet(definition, atomicNumbers, namedSets);
        


            //string symbol;
            //
            //int result;

            //for (i = 0; i < labelSize; i++)
            //    if (!isdigit(label[i]))
            //        symbol += label[i];

            //if (symbol[0] == '!')
            //{
            //    symbol = symbol.substr(1);
            //    if
            //    result = periodic_table::atomicNumber(symbol);
            //    //result = ChemicalElementData::atomicNumber(symbol);
            //    result *= -1;
            //}
            //else
            //    result = periodic_table::atomicNumber(symbol);// ChemicalElementData::atomicNumber(symbol);

    }

    void MATTS_BankReader::readSettings(
        std::istream &in,
        BankSettings &settings,
        int &lastReadLineNumber,
        std::vector<std::string> &lines)
        const
    {
        int nSettingsLines=0;
        settings = BankSettings();
        string line;
        vector<string> words{ "SETTINGS" };
        
        while (!words.empty())
        {
            nSettingsLines++;

            // process last line

            if (words[0][0] != '#')
            {
                if (line.find("covalent bond threshold") != string::npos)
                    settings.descriptorsSettings.covalentBondThreshold = stod(words.back());

                if (line.find("atom planarity threshold") != string::npos)
                    settings.descriptorsSettings.atomPlanarityThreshold = stod(words.back());

                if (line.find("ring planarity threshold") != string::npos)
                    settings.descriptorsSettings.ringPlanarityThreshold = stod(words.back());
                    
                if (line.find("atom in ring planarity threshold") != string::npos)
                    settings.descriptorsSettings.atomInRingPlanarityThreshold = stod(words.back());

                if (line.find("atom in planar ring max number of neighbours") != string::npos)
                    settings.descriptorsSettings.atomInRingMaxNeighbourCount = stoi(words.back());

                if (line.find("min plm") != string::npos)
                    settings.min_plm = stod(words.back());

                if (line.find("n sigma") != string::npos)
                    settings.nSigma = stod(words.back());

                if (line.find("minimal number of type instances") != string::npos)
                    settings.min_n_instaces = stoi(words.back());
            }

            getline(in, line);
            lastReadLineNumber++;
            lines.push_back(line);
            string_utilities::split(line, words, CharacterType::WHITE_SPACE);
        }

        if (nSettingsLines != 9)
            on_error::throwException("invalid specification of SETTINGS entry in multipole model atoms databank - expected 9 lines", __FILE__, __LINE__);
    }

} // namespace discamb
