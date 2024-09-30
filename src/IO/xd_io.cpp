#include "discamb/IO/xd_io.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/HC_Model/ClementiRoettiData.h"
#include "discamb/HC_Model/DeformationValenceParameters.h"


#include <fstream>
#include <cstdint>

using namespace std;

namespace {

    void setLocalCoordinateSystems(
        const discamb::XdMasData &xd_mas_data, 
        const discamb::Crystal &crystal, 
        std::vector<discamb::XdLocalCoordinateSystem> &lcs);
    
    void setWfnTypes(
        const discamb::XdMasData &xd_mas_data, 
        std::vector<discamb::HC_WfnType> &wfn_parameters, 
        std::vector<int> &atom_to_wfn_type_map);

    void setBnkTypes(
        const discamb::XdMasData &xd_mas_data, 
        const discamb::XdParData &xd_par_data,
        std::vector<discamb::HC_AtomTypeParameters> &type_parameters,
        std::vector<int> &atom_to_type_map,
        bool group_types, 
        bool comparePval);

    struct BnkTypeRelatedStr
    {
        double occupancy = 1.0;
        int kappaSetIdx = 0;
        double p_val = 1.0;
        vector<vector<double> > p_lm;
        int atomIdx = 0;
    };

    bool occupancyBasedArg1GreaterThanArg2(const BnkTypeRelatedStr &t1,const BnkTypeRelatedStr &t2)
    { return t1.occupancy>t2.occupancy;}

    bool equivalentTypes(const BnkTypeRelatedStr &t1, const BnkTypeRelatedStr &t2,bool comparePval);

    /**
    reads next non-comment line from xd master or parameter file if it is possible
    ommits comments and merges multiple lines if neccessery (as indicated by \ character at the nd of line)
    returns true if non-empty line was read
    */
    bool readXdLine(std::ifstream &in, std::string &line);
    void processXdMasGeneralSection(const std::vector<std::string> &generalSection,
        discamb::UnitCell &unitCell, discamb::SpaceGroup &spaceGroup);
    void processSymmetryLine(const std::string &symmLine,
        std::vector<discamb::SpaceGroupOperation> &symmetryOperations);
    void processSymmetryInfo(const std::vector<discamb::SpaceGroupOperation> &listedOperations, bool isCentrosymmetric, char lattice,
        discamb::SpaceGroup &spaceGroup);
    void extractNonequivalentOperations(std::vector<discamb::SpaceGroupOperation> &operations);

    void processXdMasAtomTable(const std::vector<std::pair<std::string, std::vector<std::string> > > &sections,
        std::vector<std::string> &atomLabels,
        std::vector<std::vector<std::string> > &lcs,
        std::vector<int> &kappa_set,
        std::vector<int> &scattering_table,
        std::vector<discamb::Vector3d> &dummyAtomsFractionalCoordinates,
        bool &dummystartsFrom0);

    void getXdMasScatteringTable(const std::vector<std::pair<std::string, std::vector<std::string> > > &sections,
        std::vector<std::vector<std::string> > &scattering_table);
    int sectionIndex(const std::vector<std::pair<std::string, std::vector<std::string> > > &sections, const std::string &sectionIndex);
    void readXdParAtomEntry(std::ifstream &in, std::string &atomLabel, discamb::Vector3d &positionFractional, std::vector<double> &adps,
        double &p_val, std::vector<std::vector<double> > &plms, int &atomToKappaSetIndex, double &occupancy, std::string &lastWordRead);
    void openXdPar(std::ifstream &in, const std::string &fileName);
    // get next non-comment and non-empty line (if it is impossible returns false)

    bool isFixedFormatRealNumberWithDot(const std::string &s)
    {

        bool hasDot = false;
        // [ digits ] . [ digits]
        if (s.empty())
            return false;
        

        for (int i = 0; i < s.size(); i++)
        {
            if (s[i] == '.')
            {
                if (hasDot)
                    return false;
                else
                    hasDot = true;
            }
            else
                if(i==0)
                {
                    if (!isdigit(s[i]) && s[i]!='-')
                        return false;
                }
                else
                    if (!isdigit(s[i]))
                        return false;
        }

        return hasDot;
    }
    
    void getXdMasScatteringTable(
        const vector<pair<string, vector<string> > > &sections,
        vector<vector<string> > &scattering_table)
    {
        int xdlsmSectionIndex;
        int lineIndex, nSectionLines;
        vector<string> words;
        bool inScatTable;

        scattering_table.clear();

        xdlsmSectionIndex = sectionIndex(sections, "xdlsm");

        if (xdlsmSectionIndex < 0)
            discamb::on_error::throwException("lack of XDLSM section in XD master file", __FILE__, __LINE__);

        const vector<string> &section = sections[xdlsmSectionIndex].second;

        nSectionLines = section.size();

        inScatTable = false;

        for (lineIndex = 0; lineIndex < nSectionLines; lineIndex++)
        {
            discamb::string_utilities::split(section[lineIndex], words);

            if (!words.empty())
            {
                if (!inScatTable)
                {
                    if (words[0] == string("SCAT"))
                        inScatTable = true;
                }
                else
                {
                    if (words.size() == 2)
                        if (words[0] == string("END") && words[1] == string("SCAT"))
                            return;

                    if (isalpha(words[0][0]))
                        scattering_table.push_back(vector<string>(1, section[lineIndex]));
                    else
                    {
                        if (scattering_table.empty())
                            discamb::on_error::throwException(string("invalid line in SCAT table in XD master file: '") + section[lineIndex] + string("'"), __FILE__, __LINE__);
                        else
                            scattering_table.back().push_back(section[lineIndex]);
                    }
                }
            } // !words.empty()
        } // for( lineIndex = 0 ; lineIndex < nSectionLines ; lineIndex++ )

    }

    int sectionIndex(
        const std::vector<std::pair<std::string, std::vector<std::string> > > &sections,
        const std::string &sectionLabel)
    {
        int i, n = sections.size();
        string lowerCaseString, lowerCaseSectionLabel, starredSectionLabel;

        discamb::string_utilities::toLower(sectionLabel, lowerCaseSectionLabel);
        starredSectionLabel = string("*") + lowerCaseSectionLabel;

        for (i = 0; i < n; i++)
        {
            discamb::string_utilities::toLower(sections[i].first, lowerCaseString);
            if (lowerCaseString == lowerCaseSectionLabel || lowerCaseString == starredSectionLabel)
                return static_cast<int>(i);
        }

        return -1;
    }

    void processXdMasAtomTable(
        const std::vector<std::pair<std::string, std::vector<std::string> > > &sections,
        std::vector<std::string> &atomLabels,
        std::vector<std::vector<std::string> > &lcs,
        std::vector<int> &kappa_set,
        std::vector<int> &scattering_table,
        std::vector<discamb::Vector3d> &dummyAtomsFractionalCoordinates,
        bool &dummyStartsFrom0)
    {
        int xdlsmSectionIndex;
        int nLines, lineIndex;
        string lowerCaseString, space = string(" ");
        vector<string> words;
        bool processingAtomsSection, atomsSectionProcessed = false;
        bool processingDummyAtoms;

        atomLabels.clear();
        lcs.clear();
        kappa_set.clear();
        scattering_table.clear();
        dummyAtomsFractionalCoordinates.clear();

        xdlsmSectionIndex = sectionIndex(sections, "xdlsm");

        if (xdlsmSectionIndex < 0)
            discamb::on_error::throwException("lack of XDLSM section in XD master file", __FILE__, __LINE__);

        const vector<string> &xdlsmSection = sections[xdlsmSectionIndex].second;

        nLines = xdlsmSection.size();
        processingAtomsSection = false;
        processingDummyAtoms = false;

        for (lineIndex = 0; lineIndex < nLines && !atomsSectionProcessed; lineIndex++)
        {
            discamb::string_utilities::split(xdlsmSection[lineIndex], words, discamb::string_utilities::CharacterType::WHITE_SPACE);
            discamb::string_utilities::toLower(words[0], lowerCaseString);
            if (lowerCaseString == string("atom"))
            {
                processingAtomsSection = true;
            }
            else
            {
                if (processingAtomsSection)
                {

                    if (lowerCaseString == string("end"))
                        atomsSectionProcessed = true;
                    else
                    {
                        if (!processingDummyAtoms)
                        {
                            if (words[0] == string("DUM0") || words[0] == string("DUM1"))
                            {
                                dummyStartsFrom0 = (words[0] == string("DUM0"));
                                processingDummyAtoms = true;
                                dummyAtomsFractionalCoordinates.push_back(discamb::Vector3d(atof(words[1].c_str()), atof(words[2].c_str()), atof(words[3].c_str())));
                                continue;
                            }
                        }
                        else
                        {
                            dummyAtomsFractionalCoordinates.push_back(discamb::Vector3d(atof(words[1].c_str()), atof(words[2].c_str()), atof(words[3].c_str())));
                            continue;
                        }


                        if (words.size() < 10)
                            discamb::on_error::throwException(string("invalid line in atom table of XDLSM section in XD master file:'") + xdlsmSection[lineIndex] + string("'"), __FILE__, __LINE__);

                        atomLabels.push_back(words[0]);
                        lcs.push_back(vector<string>(words.begin() + 1, words.begin() + 7));
                        kappa_set.push_back(int(atoi(words[9].c_str())) - 1);
                        scattering_table.push_back(int(atoi(words[8].c_str())) - 1);
                    }
                }
            }
        }
    }


    void processXdMasGeneralSection(
        const std::vector<std::string> &generalSection,
        discamb::UnitCell &unitCell,
        discamb::SpaceGroup &spaceGroup)
    {
        int lineIndex, nLines;
        string uppercaseString;
        double cellParameters[] = { 1,1,1,90,90,90 };
        char lattice[] = { 'C','P' };
        vector<string> words;
        vector<discamb::SpaceGroupOperation> symmOperations, lineSymmetryOperations;
        int i;

        unitCell.set(1, 1, 1, 90, 90, 90);
        spaceGroup.set(vector<discamb::SpaceGroupOperation>(1));
        nLines = generalSection.size();

        for (lineIndex = 0; lineIndex < nLines; lineIndex++)
        {
            discamb::string_utilities::split(generalSection[lineIndex], words);
            discamb::string_utilities::toUpper(words[0], uppercaseString);

            if (uppercaseString == string("CELL"))
            {
                for (i = 1; i < words.size() && i <= 6; i++)
                    cellParameters[i - 1] = atof(words[i].c_str());
                unitCell.set(cellParameters[0], cellParameters[1], cellParameters[2],
                    cellParameters[3], cellParameters[4], cellParameters[5]);
            }

            if (uppercaseString == string("SYMM"))
            {
                processSymmetryLine(generalSection[lineIndex], lineSymmetryOperations);
                symmOperations.insert(symmOperations.end(), lineSymmetryOperations.begin(), lineSymmetryOperations.end());
            }
            if (uppercaseString == string("LATT"))
            {
                discamb::string_utilities::split(generalSection[lineIndex], words);
                for (i = 1; i < words.size() && i < 3; i++)
                    lattice[i - 1] = words[i][0];
            }
        }
        // chceck if identity operation is included
        discamb::Matrix3<double> rotation, identity;
        discamb::Vector3d translation, zero_vector;
        bool has_identity = false;

        identity.setToIdentity();

        for (i = 0; i < symmOperations.size(); i++)
        {
            symmOperations[i].get(rotation, translation);
            if (rotation == identity && translation == zero_vector)
                has_identity = true;
        }

        //---

        if (!has_identity)
            symmOperations.push_back(discamb::SpaceGroupOperation());

        spaceGroup.set(symmOperations, lattice[1], toupper(lattice[0]) == 'C');


    }


    void extractNonequivalentOperations(
        std::vector<discamb::SpaceGroupOperation> &operations)
    {
        int i, j, n, nUnique;
        std::vector<discamb::SpaceGroupOperation> symm;
        bool unique;
        discamb::Vector3i translation;
        n = operations.size();
        nUnique = 0;

        for (i = 0; i < n; i++)
        {
            unique = true;
            for (j = 0; j < nUnique; j++)
                if (symm[j].isLatticeTranslationEquivalent(operations[i], translation))
                {
                    unique = false;
                    break;
                }
            if (unique)
            {
                symm.push_back(operations[i]);
                ++nUnique;
            }
        }
        operations.swap(symm);
    }

    // SYMM X, Y, Z
    // SYMM X, Y, Z; -X,-Y,-Z
    // SYMM 0. 1. 0. 0.    0.  

    void processSymmetryLine(
        const std::string &symmLine,
        std::vector<discamb::SpaceGroupOperation> &symmetryOperations)
    {
        bool symmAsMatrix;
        int i, operationIndex;
        string errorMessage;
        string operationsLine; // the line with keyword SYMM extracted
        vector<string> words, operations;
        bool commaSeparatedComponents;
        discamb::Vector3d translation;
        discamb::Matrix3d rotation;

        symmetryOperations.clear();

        operationsLine = symmLine.substr(symmLine.find("SYMM") + 5);



        // check if symmetry operations are given as string (e.g. X,Y,Z) or numerical (matrix) 

        symmAsMatrix = true;

        for (i = 1; i < operationsLine.size(); i++)
            if (isalpha(operationsLine[i]))
                symmAsMatrix = false;

        // 

        discamb::string_utilities::split(operationsLine, operations, ';');

        for (operationIndex = 0; operationIndex < operations.size(); operationIndex++)
        {

            commaSeparatedComponents = (operations[operationIndex].find(',') != string::npos);
            symmetryOperations.resize(symmetryOperations.size() + 1);

            if (symmAsMatrix)
            {
                discamb::on_error::throwException("definition of symmetry operation with matrix in XD file not supported", __FILE__, __LINE__);
            }
            else
            {
                if (!commaSeparatedComponents)
                {
                    discamb::string_utilities::split(operations[operationIndex], words);
                    operations[operationIndex] = discamb::string_utilities::merge(words, ',');
                }

                symmetryOperations.back().set(operations[operationIndex]);
            }
        }

    }

    void openXdPar(
        std::ifstream &in,
        const std::string &fileName)
    {
        if (fileName.empty()) // if no file name is given tries to find file with default name
        {
            in.open("xd.par");
            if (!in.good())
            {
                in.close();
                in.open("xd.res");
                if (!in.good())
                {
                    in.close();
                    in.open("xd.inp");
                    if (in.good())
                        discamb::on_error::throwException("can not read xd parameter file with default name (tried xd.par, xd.res, xd.inp)",
                            __FILE__, __LINE__);
                }
            }
        }
        else
        {
            in.open(fileName.c_str());
            if (!in.good())
                discamb::on_error::throwException(string("can not read xd parameter file : '") + fileName + string("'"), __FILE__, __LINE__);
        }

    }


    void readXdParAtomEntry(
        std::ifstream &in,
        std::string &atomLabel,
        discamb::Vector3d &positionFractional,
        std::vector<double> &adps,
        double &p_val,
        std::vector<std::vector<double> > &plms,
        int &atomToKappaSetIndex,
        double &occupancy,
        std::string &lastWordRead)
    {
        int l, maxL, auxInt, m, abs_m, sign_m;
        int adpTensorOrder;
        double auxReal;
        static const int nAdpParameters[] = { 1, 1, 6, 16, 31 };
        static const int nAdpWords[] = { 6,6,6,16,31 };

        string word,zeroFloatAsStr;

        if (lastWordRead.empty())
            in >> atomLabel;
        else
            atomLabel = lastWordRead;


        // read lcs
        for (int i = 0; i < 5; i++)
            in >> auxInt;

        in >> adpTensorOrder >> auxInt >> atomToKappaSetIndex >> maxL >> auxInt >> auxInt;

        atomToKappaSetIndex--;

        if (adpTensorOrder > 4)
            discamb::on_error::throwException(string("problem when reading XD parameter file - invalid order if ADP tensor: '") +
                discamb::string_utilities::convertToString(adpTensorOrder) + string("'"), __FILE__, __LINE__);
        else
        {
            // xyz and occupancy

            in >> positionFractional[0] >> positionFractional[1] >> positionFractional[2] >> occupancy;

            // adps

            adps.resize(nAdpParameters[adpTensorOrder]);


            for (int i = 0; i < nAdpWords[adpTensorOrder]; i++)
            {
                if (i < nAdpParameters[adpTensorOrder])
                    in >> adps[i];
                else
                    in >> auxReal;
            }

            // plms

            in >> p_val;

            if (occupancy != 0.0)
                p_val /= occupancy;
            int nL;
            maxL < 0 ? nL = 0 : nL = int(maxL) + 1; // nL = maxL + 1
            plms.resize(nL);
            plms[0].resize(1);
            in >> plms[0][0];
            if (maxL > 0)
            {
                plms[1].resize(3);
                in >> plms[1][2] >> plms[1][0] >> plms[1][1];
            }

            for (l = 2; l <= maxL; l++)
            {
                plms[l].resize(2 * int(l) + 1);
                for (int i = 0; i < 2 * l + 1; i++)
                {
                    abs_m = (i + 1) / 2;
                    i % 2 == 1 ? sign_m = 1 : sign_m = -1;
                    m = sign_m * abs_m;
                    int idx = m + l;
                    in >> plms[l][int(idx)];
                }
            }

            if (occupancy != 0)
                for (l = 0; l <= maxL; l++)
                    for (int i = 0; i < 2 * l + 1; i++)
                        plms[l][i] /= occupancy;
            /*
            for (l = maxL + 1; l <= 4; l++)
                for (int i = 0; i < 2 * l + 1; i++)
                    in >> auxReal;
             */

             // handling the case when there is more Plms than corresponding to L_max 

            in >> word;

            zeroFloatAsStr = string("0.0000");

            if (isFixedFormatRealNumberWithDot(word))
                do
                {
                    in >> word;
                } while (isFixedFormatRealNumberWithDot(word));

                //in.seekg(in.tellg()-std::streampos(word.size()) - 1);
                lastWordRead = word;
        }
    }



    bool readXdLine(
        std::ifstream &in,
        std::string &newLine)
    {
        bool readLine = false;
        string totalLine, line, trimmedLine;

        while (!readLine && in.good())
        {
            discamb::string_utilities::portable_getline(in, line);
            if (line.empty())
                continue;
            if (line[0] == '!')
                continue;

            trimmedLine = discamb::string_utilities::trim(line);
            if (trimmedLine.empty())
                continue;

            if (line[line.size() - 1] == '\\')
                totalLine += line.substr(0, line.size() - 1);
            else
            {
                totalLine += line;
                readLine = true;
            }

        }

        if (!readLine && !totalLine.empty())
            discamb::on_error::throwException(string("problem with reading XD master file - ") +
                string("continuation of line with forwardslash") +
                string(" charater (\\) but no next line present"),
                __FILE__, __LINE__);
        newLine.swap(totalLine);
        return readLine;

    }


    void setLocalCoordinateSystems(
        const discamb::XdMasData &xd_mas_data,
        const discamb::Crystal &crystal,
        std::vector<discamb::XdLocalCoordinateSystem> &lcs)
    {
        int atom_index, nAtoms;

        nAtoms = xd_mas_data.lcs.size();

        lcs.resize(nAtoms);

        //xd_mas_data.lcs[i][0]  [1]  [2]   [3]  [4]  [5]
        //   O(1)          C(1)   X   O(1)  N(1)  Y    R   0  1   1   4  NO      

        for (atom_index = 0; atom_index<nAtoms; atom_index++)
            lcs[atom_index].set(xd_mas_data.atomLabels[atom_index], xd_mas_data.lcs[atom_index][0],
                                xd_mas_data.lcs[atom_index][2], xd_mas_data.lcs[atom_index][3],
                                xd_mas_data.lcs[atom_index][1], xd_mas_data.lcs[atom_index][4],
                                xd_mas_data.lcs[atom_index][5][0] == 'R', crystal,
                                xd_mas_data.dummyAtomIndexingFrom0, xd_mas_data.dummy_atoms);

    }


    void setWfnTypes(
        const discamb::XdMasData &xd_mas_data,
        std::vector<discamb::HC_WfnType> &wfn_parameters,
        std::vector<int> &atom_to_wfn_type_map)
    {
        int nScatTypes;
        int typeIdx;
        vector<double> orb_coeff, orb_exp;
        vector<int> orb_pow;
        string wfn_type;
        vector<string> words;
        const vector<vector<string> > &scat_tables = xd_mas_data.scat_tables;

        nScatTypes = scat_tables.size();
        wfn_parameters.resize(nScatTypes);
        
        discamb::ClementiRoettiData wfns;
        discamb::DeformationValenceParameters def_val;//(discamb::DeformationValenceParameters::STANDARD);
        bool known_type;

        for (typeIdx = 0; typeIdx<nScatTypes; typeIdx++)
        {
            discamb::string_utilities::split(scat_tables[typeIdx][0], words);
            wfn_type = words[0];            

            wfn_parameters[typeIdx] = wfns.getEntry(wfn_type, known_type);

            if(!known_type)
                discamb::on_error::throwException("a request for wavefunction data in mulitpolar model for an unknown type of atom",__FILE__,__LINE__);


            known_type = def_val.getParameters(wfn_type, wfn_parameters[typeIdx].deformation_valence_exponent,
                                               wfn_parameters[typeIdx].deformation_valence_power);

            if (!known_type)
            {
                string errorMessage  = string("a request for deformation valence data in mulitpolar model for an unknown type of atom: '")+
                                       wfn_type + string("'");
                discamb::on_error::throwException(errorMessage, __FILE__, __LINE__);
            }

            if (words[3] == string("RDSD"))
            {
                discamb::string_utilities::split(scat_tables[typeIdx][1], words);

                wfn_parameters[typeIdx].deformation_valence_power.resize(5);
                wfn_parameters[typeIdx].deformation_valence_exponent = atof(words[1].c_str());

                for (int i = 0; i<5; i++)
                    wfn_parameters[typeIdx].deformation_valence_power[i] = atoi(words[2 * i].c_str());
            }

        }

        atom_to_wfn_type_map = xd_mas_data.atomScatTableIndex;
    }

    void setBnkTypes(
        const discamb::XdMasData &xd_mas_data,
        const discamb::XdParData &xd_par_data,
        std::vector<discamb::HC_AtomTypeParameters> &type_parameters,
        std::vector<int> &atom_to_type_map,
        bool group_types,
        bool comparePval)
    {

        // exact comparison p_val, p_lm & kappa set index
        // kappa set taken from mas

        int i, j, refAtom, kappa_set, nTypes, nAtoms = xd_par_data.plms.size();
        
        if(nAtoms == 0)
            return;

        atom_to_type_map.resize(nAtoms);

        if (group_types)
        {
            vector<BnkTypeRelatedStr> type_params,atom_params(nAtoms);

            // plms, kappa, pval

            for (i = 0; i<nAtoms; i++)
            {
                atom_params[i].atomIdx = i;
                atom_params[i].kappaSetIdx = xd_mas_data.atomKappaSetIndex[i];
                atom_params[i].occupancy = xd_par_data.occupancy[i];
                atom_params[i].p_lm = xd_par_data.plms[i];
                atom_params[i].p_val = xd_par_data.p_val[i];
            }

            // sort parameters so the ones for atoms with higher occupancy are at the beginning of the table
            sort(atom_params.begin(),atom_params.end(), occupancyBasedArg1GreaterThanArg2);

            // find atom types
            type_params.push_back(atom_params[0]);
            atom_to_type_map[atom_params[0].atomIdx] = 0;


            bool found;
            for(i=1;i<nAtoms;i++)
            {
                nTypes = type_params.size();
                found = false;
                for( j=0; j<nTypes ; j++)
                    if( equivalentTypes(atom_params[i], type_params[j], comparePval ))
                    {
                        atom_to_type_map[atom_params[i].atomIdx] = j;
                        found = true;
                        break;
                    }
                
                if(!found)
                {
                    type_params.push_back(atom_params[i]);
                    atom_to_type_map[atom_params[i].atomIdx] = type_params.size()-1;
                }
            }

            nTypes = type_params.size();
            type_parameters.clear();
            type_parameters.resize(nTypes);

            for (i = 0; i<nTypes; i++)
            {
                refAtom = type_params[i].atomIdx;
                kappa_set = xd_mas_data.atomKappaSetIndex[refAtom];

                type_parameters[i].kappa_deformation_valence = xd_par_data.kappa_sets[kappa_set][1];
                type_parameters[i].kappa_spherical_valence = xd_par_data.kappa_sets[kappa_set][0];
                type_parameters[i].p_lm = type_params[i].p_lm;
                type_parameters[i].p_val = xd_par_data.p_val[refAtom];
            }
        } // if(group_types)
        else
        {
            
            type_parameters.clear();
            type_parameters.resize(nAtoms);

            for (i = 0; i<nAtoms; i++)
            {
                atom_to_type_map[i] = i;
                kappa_set = xd_mas_data.atomKappaSetIndex[i];

                type_parameters[i].kappa_deformation_valence = xd_par_data.kappa_sets[kappa_set][1];
                type_parameters[i].kappa_spherical_valence = xd_par_data.kappa_sets[kappa_set][0];
                type_parameters[i].p_lm = xd_par_data.plms[i];
                type_parameters[i].p_val = xd_par_data.p_val[i];
            }

        }
    }


    bool equivalentTypes(
        const BnkTypeRelatedStr &t1,
        const BnkTypeRelatedStr &t2,
        bool comparePval)
    {
        

        if(t1.kappaSetIdx != t2.kappaSetIdx)
            return false;

        int l, mIdx, nL, nM;

        // check if p_lm tables are of equal size
        nL = t1.p_lm.size();
        if (nL != t2.p_lm.size())
            return false;
        

        if(t1.occupancy == 0)
        {
            if( t2.p_val !=0 )
                return false;
            
            for( l=0 ; l<nL; l++ )
            {
                nM = 2*l+1;
                for ( mIdx = 0 ; mIdx< nM ; mIdx++)
                    if( t2.p_lm[l][mIdx] !=0 )
                        return false;
            }
            return true;
        }
        
        if (t2.occupancy == 0)
        {
            if (t1.p_val != 0)
                return false;

            for ( l = 0; l<nL; l++)
            {
                nM = 2*l+1;
                for ( mIdx = 0 ; mIdx<nM; mIdx++)
                    if (t1.p_lm[l][mIdx] != 0)
                        return false;
            }
            return true;
        }

        double threshold = 4.999e-5*(1.0/t1.occupancy + 1.0/t2.occupancy);
        

        for(l=0;l<nL;l++)
        {
            nM = 2*l+1;
            for(mIdx=0;mIdx<nM;mIdx++)
                if(fabs(t1.p_lm[l][mIdx]-t2.p_lm[l][mIdx])> threshold)
                    return false;
        }

        if(comparePval)
            if (fabs(t1.p_val - t2.p_val)> threshold)
                return false;

        return true;

    } // end of : bool equivalentTypes(..)

}


namespace discamb {

    namespace xd_io {

        void read(
            const std::string & masterFileName,
            const std::string & parameterFileName,
            Crystal &crystal)
        {
            XdMasData xdMasData;
            XdParData xdParData;

            readXdMas(masterFileName, xdMasData);
            readXdPar(parameterFileName, xdParData);

            getCrystal(xdMasData, xdParData, crystal);

        }


        void readXdMas(
            const std::string & masterFileName,
            XdMasData &xdMasData)
        {
            vector<string> generalSection;
            vector< pair< string, vector< string > > > sections;

            readXdMasSections(masterFileName, generalSection, sections);
            processXdMasGeneralSection(generalSection, xdMasData.unitCell, xdMasData.spaceGroup);

            processXdMasAtomTable(sections, xdMasData.atomLabels, xdMasData.lcs, xdMasData.atomKappaSetIndex,
                xdMasData.atomScatTableIndex, xdMasData.dummy_atoms, xdMasData.dummyAtomIndexingFrom0);

            getXdMasScatteringTable(sections, xdMasData.scat_tables);
        }


        void readXdPar(
            const std::string & parameterFileName,
            XdParData &xdParData)
        {
            int i, atomIndex, nAtoms, nDummyAtoms, dummyAtomIndex, nKappaSets;
            ifstream in;
            string line, word;
            vector<string> words;

            openXdPar(in, parameterFileName);

            readXdLine(in, line); // version line
            readXdLine(in, line); // compoun ID etc
            readXdLine(in, line); // limits
            readXdLine(in, line); // arrays dimensions

            string_utilities::split(line, words);

            nAtoms = atoi(words[1].c_str());
            nDummyAtoms = atoi(words[14].c_str());

            nKappaSets = atoi(words[4].c_str());

            // read statistics of the fit and parameters of the lsq weight
            // 8+6 = 14 words
            for (i = 0; i < 14; i++)
                in >> word;

            xdParData.dummyAtomsPositionsFractional.resize(nDummyAtoms);

            for (dummyAtomIndex = 0; dummyAtomIndex < nDummyAtoms; dummyAtomIndex++)
                for (i = 0; i < 3; i++)
                    in >> xdParData.dummyAtomsPositionsFractional[dummyAtomIndex][i];

            xdParData.adps.resize(nAtoms);
            xdParData.atomLabels.resize(nAtoms);
            xdParData.atomToKappaSetMap.resize(nAtoms);
            xdParData.plms.resize(nAtoms);
            xdParData.p_val.resize(nAtoms);
            xdParData.positionsFractional.resize(nAtoms);
            xdParData.occupancy.resize(nAtoms);
            string lastWordRead;
            for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
                readXdParAtomEntry(in, xdParData.atomLabels[atomIndex], xdParData.positionsFractional[atomIndex],
                    xdParData.adps[atomIndex], xdParData.p_val[atomIndex], xdParData.plms[atomIndex], xdParData.atomToKappaSetMap[atomIndex],
                    xdParData.occupancy[atomIndex], lastWordRead);

            xdParData.kappaSetToScatTableMap.resize(nKappaSets);
            xdParData.kappa_sets.resize(nKappaSets, vector<double>(6));
            for (i = 0; i < nKappaSets; i++)
            {
                if (i == 0)
                    xdParData.kappaSetToScatTableMap[i] = stoi(lastWordRead);
                else
                    in >> xdParData.kappaSetToScatTableMap[i];
                in >> xdParData.kappa_sets[i][0] >> xdParData.kappa_sets[i][1]
                    >> xdParData.kappa_sets[i][2] >> xdParData.kappa_sets[i][3] >> xdParData.kappa_sets[i][4]
                    >> xdParData.kappa_sets[i][5];
            }
            in.close();
        }

        void readXdMasSections(
            const std::string & masterFileName,
            std::vector<std::string> &generalSection,
            std::vector<std::pair<std::string, std::vector<string> > > &sections)
        {
            ifstream in(masterFileName.c_str());

            if (!in.good())
                on_error::throwException(string("can not read xd master file: '") + masterFileName + string("'"), __FILE__, __LINE__);

            vector<string> words;
            string line;
            string uppercaseString;
            bool generalSectionRead;

            generalSection.clear();
            sections.clear();
            generalSectionRead = false;

            while (readXdLine(in, line))
            {
                string_utilities::split(line, words);

                // is beginning of next module section
                string_utilities::toUpper(words[0], uppercaseString);

                if (uppercaseString == string("MODULE"))
                {
                    if (words.size() == 1)
                        on_error::throwException("problem when reading XD master file, keyword MODULE not followed by module name", __FILE__, __LINE__);

                    sections.resize(sections.size() + 1);
                    sections.back().first = words[1];
                    generalSectionRead = true;
                }

                if (generalSectionRead)
                    sections.back().second.push_back(line);
                else
                    generalSection.push_back(line);
            }

            in.close();
        }

        void read(
            const std::string & mas_file,
            const std::string & par_file,
            HC_ModelParameters &hc_parameters,
            Crystal &crystal,
            std::vector<XdLocalCoordinateSystem> &localCoordinateSystems,
            bool group_types,
            bool comparePval)
        {
            XdMasData xdMasData;
            XdParData xdParData;

            readXdMas(mas_file, xdMasData);
            readXdPar(par_file, xdParData);


            setWfnTypes(xdMasData, hc_parameters.wfn_parameters, hc_parameters.atom_to_wfn_map);


            setBnkTypes(xdMasData, xdParData, hc_parameters.type_parameters, hc_parameters.atom_to_type_map,
                group_types, comparePval);

            getCrystal(xdMasData, xdParData, crystal);

            setLocalCoordinateSystems(xdMasData, crystal, localCoordinateSystems);

        }

        void getCrystal(
            const XdMasData &xdMasData,
            const XdParData &xdParData,
            Crystal &crystal)
        {
            string scattere_type;
            vector<string> words;
            int nSpaceGroupOperations = xdMasData.spaceGroup.nSymmetryOperations();
            crystal.unitCell = xdMasData.unitCell;
            crystal.spaceGroup = xdMasData.spaceGroup;
            crystal.atoms.resize(xdMasData.atomLabels.size());

            for (int i = 0; i < xdMasData.atomLabels.size(); i++)
            {
                crystal.atoms[i].label = xdMasData.atomLabels[i];

                string_utilities::split(xdMasData.scat_tables[xdMasData.atomScatTableIndex[i]][0], words);
                scattere_type = words[0];
                crystal.atoms[i].type = scattere_type;
                crystal.atoms[i].occupancy = xdParData.occupancy[i];
                crystal.atoms[i].multiplicity = nSpaceGroupOperations;
            }

            for (int i = 0; i < xdParData.adps.size(); i++)
            {
                crystal.atoms[i].adp = xdParData.adps[i];
                crystal.atoms[i].coordinates = xdParData.positionsFractional[i];
            }
            crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
            crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;

        }


        void write(
            const std::string &xdMasName, 
            const std::string &xdParName, 
            const Crystal &c, 
            const HC_ModelParameters &params,
            const std::vector<XdLocalCoordinateSystem> &lcs)
        {
            on_error::not_implemented(__FILE__, __LINE__);
        }

        void writeMasterFile(
            const std::string &xdMasName,
            const Crystal &c,
            const HC_ModelParameters &params,
            const std::vector<XdLocalCoordinateSystem> &lcs)
        {
            on_error::not_implemented(__FILE__, __LINE__);
            ofstream out(xdMasName.c_str());
        }
#ifdef _MSC_VER
        void readXdFou(
            const std::string &fileName,
            std::vector<Vector3i> &hkl,
            std::vector<double> &f_obs,
            std::vector<double> &sigma,
            std::vector<double> &phase_obs,
            std::vector<std::complex<double> > &f_model_1,
            std::vector<std::complex<double> > &f_model_2)
        {
            hkl.clear();
            f_obs.clear();
            sigma.clear();
            phase_obs.clear();
            f_model_1.clear();
            f_model_2.clear();

            ifstream in(fileName, ios::in | ios::binary);
            char c;
            int32_t h, k, l;
            bool read = true;

            float f, g;

            while (read && in.good())
            {
                in.read(&c, 1);

                in.read(reinterpret_cast<char*>(&h), 4);
                in.read(reinterpret_cast<char*>(&k), 4);
                in.read(reinterpret_cast<char*>(&l), 4);

                if (h == 0 && k == 0 && l == 0)
                    read = false;

                hkl.push_back(Vector3i(h, k, l));

                in.read(reinterpret_cast<char*>(&f), 4);
                f_obs.push_back(f);
                
                in.read(reinterpret_cast<char*>(&f), 4);
                sigma.push_back(f);
                
                in.read(reinterpret_cast<char*>(&f), 4);
                phase_obs.push_back(f);

                in.read(reinterpret_cast<char*>(&f), 4);
                in.read(reinterpret_cast<char*>(&g), 4);
                f_model_1.push_back(complex<double>(f, g));

                in.read(reinterpret_cast<char*>(&f), 4);
                in.read(reinterpret_cast<char*>(&g), 4);
                f_model_2.push_back(complex<double>(f, g));

                in.read(&c, 1);

            }
            in.close();
        }

        void writeXdFou(const std::string &fileName,
            const std::vector<Vector3i> &hkl,
            const std::vector<double> &f_obs,
            const std::vector<double> &sigma,
            const std::vector<double> &phase_obs,
            const std::vector<std::complex<double> > &fModel1,
            const std::vector<std::complex<double> > &fModel2)
        {
            ofstream out(fileName, ios::out | ios::binary);
            //__int8 i;
            char c = 40;
            __int32 h, k, l;
            bool read = true;


            float f;
            for (int i = 0; i < hkl.size(); i++)
            {
                out.write(&c, 1);

                h = hkl[i][0];
                k = hkl[i][1];
                l = hkl[i][2];

                out.write(reinterpret_cast<char*>(&h), 4);
                out.write(reinterpret_cast<char*>(&k), 4);
                out.write(reinterpret_cast<char*>(&l), 4);

                /*
                fobs observed structure factor (normally with anomalous dispersion removed)
                sig error of fobs
                phase phase angle calculated with the final parameters according to the model the
                refinement was based on
                amod1 real part of the calculated structure factor (fmod1) based on an input dependent
                model (MODEL1)
                bmod1 imaginary part of fmod1
                amod2 real part of fmod2
                bmod2 imaginary part of fmod2
                */

                f = static_cast<float>(f_obs[i]);
                out.write(reinterpret_cast<char*>(&f), 4);
                f = static_cast<float>(sigma[i]);
                out.write(reinterpret_cast<char*>(&f), 4);
                f = static_cast<float>(phase_obs[i]);
                out.write(reinterpret_cast<char*>(&f), 4);
                f = static_cast<float>(fModel1[i].real());
                out.write(reinterpret_cast<char*>(&f), 4);
                f = static_cast<float>(fModel1[i].imag());
                out.write(reinterpret_cast<char*>(&f), 4);
                f = static_cast<float>(fModel2[i].real());
                out.write(reinterpret_cast<char*>(&f), 4);
                f = static_cast<float>(fModel2[i].imag());
                out.write(reinterpret_cast<char*>(&f), 4);
                out.write(&c, 1);
            }


            out.close();

        }

#endif

    }
}