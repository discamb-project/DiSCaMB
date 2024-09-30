
#include "discamb/IO/mol2_io.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/StructuralProperties/CovalentRadiousBondDetector.h"
#include "discamb/StructuralProperties/GenericConnectivityAlgorithm.h"
#include "discamb/BasicChemistry/PeriodicTable.h"

#include <fstream>

using namespace std;

namespace discamb
{
    namespace mol2_io {


        namespace {
            void line_read_error(
                int lineIdx,
                const string& line,
                const string& fileName)
            {
                string error_message = "a problem occured when processing mol2 file '" + fileName +
                    "' at line " + to_string(lineIdx) + " :\n" + line + "\n";
                on_error::throwException(error_message, __FILE__, __LINE__);
            }

        }
        void Mol2Data::toMoleculeData(
            MoleculeData& data)
            const
        {
            data.atomLabels = atomName;
            data.atomPositions =atomPosition;
            
            for (const string& sybylType : atomType)
            {
                string s = sybylType.substr(0, sybylType.find('.'));
                int z = periodic_table::atomicNumber(s);
                // no conversion to atomic number for atom types like Het, Hev, Hal
                if ( z == 0 && (s != "Any" && s != "Du"))
                    on_error::throwException("problem when transforming mol2 file data - cannot assign atomic number to atom type '" + sybylType,
                        __FILE__, __LINE__);
                data.atomicNumbers.push_back(z);
            }
            data.bonds = bonds;
            data.comment = moleculeComment;
        }

        void Mol2Data::split(
            std::vector<Mol2Data>& substructures)
            const
        {
            substructures.clear();
            if (substructureIdx.empty())
                return;

            int structureIdx, nSubstructures = *max_element(substructureIdx.begin(), substructureIdx.end());

            if (nSubstructures == 1)
            {
                substructures.push_back(*this);
                return;
            }

            substructures.resize(nSubstructures);
            vector<vector<int> > atomIndices(nSubstructures);
            int nAtoms = atomId.size();
            vector<int> structure2substructureIdx(nAtoms);

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                structureIdx = substructureIdx[atomIdx] - 1;
                structure2substructureIdx[atomIdx] = atomIndices[structureIdx].size();
                atomIndices[structureIdx].push_back(atomIdx);

                substructures[structureIdx].atomId.push_back(atomId[atomIdx]);
                substructures[structureIdx].atomName.push_back(atomName[atomIdx]);
                substructures[structureIdx].atomPosition.push_back(atomPosition[atomIdx]);
                substructures[structureIdx].atomType.push_back(atomType[atomIdx]);
            }

            if(!substructureName.empty())
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    structureIdx = substructureIdx[atomIdx] - 1;
                    substructures[structureIdx].substructureName.push_back(substructureName[atomIdx]);
                }


            if (!atomicCharge.empty())
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    structureIdx = substructureIdx[atomIdx] - 1;
                    substructures[structureIdx].atomicCharge.push_back(atomicCharge[atomIdx]);
                }


            for (int i = 0; i < bonds.size(); i++)
            {
                int substructureIdx1 = substructureIdx[bonds[i].first - 1] - 1;
                int substructureIdx2 = substructureIdx[bonds[i].second - 1] - 1;

                if (substructureIdx1 != substructureIdx2)
                    continue;

                int idx1 = structure2substructureIdx[bonds[i].first - 1] + 1;
                int idx2 = structure2substructureIdx[bonds[i].second - 1] + 1;

                substructures[substructureIdx1].bonds.push_back({ idx1, idx2 });
                if (!bondTypes.empty())
                    substructures[substructureIdx1].bondTypes.push_back(bondTypes[i]);
            }

            for (int i = 0; i < nSubstructures; i++)
            {
                if(!substructures[i].substructureName.empty())
                    substructures[i].moleculeName = substructures[i].substructureName[0];
                substructures[i].chargeType = chargeType;
                substructures[i].substructureIdx.assign(substructures[i].atomId.size(), 1);
            }
            
        }

        void read(
            const std::string fileName,
            MoleculeData& data)
        {
            Mol2Data mol2Data;
            read(fileName, mol2Data);
            mol2Data.toMoleculeData(data);
        }

        void read(
            const std::string fileName,
            Mol2Data& data)
        {
            // clear data

            data.atomicCharge.clear();
            data.atomId.clear();
            data.atomName.clear();
            data.atomPosition.clear();
            data.atomType.clear();
            data.bonds.clear();
            data.bondTypes.clear();
            data.chargeType = "";
            data.moleculeComment = "";
            data.substructureIdx.clear();
            data.substructureName.clear();

            //

            ifstream in(fileName);
            //cout << "mol_io::read(" << fileName << ")\n";
            if (!in.good())
                on_error::throwException("cannot open mole2 file '" + fileName + "' for reading", __FILE__, __LINE__);

            string line;
            vector<string> words, words2;
            enum class Section { ATOM, BOND, MOLECULE, NONE, OTHER };
            Section section = Section::NONE;
            bool newSectionLine;
            int lineIdx = 0;
            while (in.good())
            {
                newSectionLine = false;
                getline(in, line);
                lineIdx++;
                string_utilities::split(line, words);

                if (words.empty())
                    continue;
                if (words[0][0] == '#')
                    continue;

                if (words[0][0] == '@')
                {
                    newSectionLine = true;
                    string_utilities::split(words[0], words2, '>');

                    if (words2.size() != 2)
                        line_read_error(lineIdx, line, fileName);

                    section = Section::OTHER;
                    if (words2[1] == string("ATOM"))
                        section = Section::ATOM;
                    if (words2[1] == string("MOLECULE"))
                        section = Section::MOLECULE;
                    if (words2[1] == string("BOND"))
                        section = Section::BOND;
                }

                //                @<TRIPOS>ATOM
                //                    1 Br1 - 1.7294 - 6.2022 - 0.1102 Br        1 RES1 - 1.0000
                //                    2 F1 - 1.2479 - 8.5568 - 7.2869 F         2 RES2    0.0000
                //                    3 F2 - 0.1615 - 7.1127 - 2.9737 F         2 RES2    0.0000


                if (section == Section::ATOM && !newSectionLine)
                {
                    int nWords = words.size();
                    if (nWords < 6)
                        line_read_error(lineIdx, line, fileName);

                    data.atomId.push_back(stoi(words[0]));
                    data.atomName.push_back(words[1]);
                    data.atomPosition.push_back(Vector3d(stod(words[2]), stod(words[3]), stod(words[4])));
                    data.atomType.push_back(words[5]);
                    if (nWords > 6)
                        data.substructureIdx.push_back(stoi(words[6]));
                    if (nWords > 7)
                        data.substructureName.push_back(words[7]);
                    if (nWords > 8)
                        data.atomicCharge.push_back(stod(words[8]));
                }
                if (section == Section::BOND && !newSectionLine)
                {
                    int nWords = words.size();
                    if (nWords < 4)
                        line_read_error(lineIdx, line, fileName);
                    data.bonds.push_back({ stoi(words[1]),stoi(words[2]) });
                    //data.bondTypes.push_back(bondType(words[3]));
                    data.bondTypes.push_back(words[3]);
                }
                if (section == Section::MOLECULE && !newSectionLine)
                {
                    // 1
                    //getline(in, line);
                    //lineIdx++;
                    data.moleculeName = line;
                    // 2, 3
                    getline(in, line);
                    getline(in, line);
                    lineIdx += 2;
                    data.moleculeType = line;
                    // 4
                    getline(in, line);
                    data.chargeType = line;
                    lineIdx++;

                    // 5, 6
                    getline(in, line);
                    getline(in, line);
                    data.moleculeComment = line;
                    lineIdx += 2;
                }
            }
            in.close();
        }

        void atomicNumbers(
            const Mol2Data& molData,
            std::vector<int>& atomicNumber)
        {
            vector<string> words;
            int atomIdx, nAtoms = molData.atomType.size();
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                string_utilities::split(molData.atomType[atomIdx], words, '.');
                atomicNumber.push_back(periodic_table::atomicNumber(words[0]));
            }

        }



        void write(
            const std::string fileName,
            const std::vector<Mol2Data>& _data)
        {
            if (_data.size() > 1)
                on_error::throwException("multiple molecules not supported for writing mol2 file", __FILE__, __LINE__);

            if (_data.empty())
                on_error::throwException("no data for writing mol2 file", __FILE__, __LINE__);

            ofstream out(fileName);

            Mol2Data const& data = _data[0];

            if (!out.good())
                on_error::throwException("cannot write mol2 file '" + fileName + "'", __FILE__, __LINE__);

            // ---------------------------------------
            //           @<TRIPOS>MOLECULE 
            int nAtoms = data.atomPosition.size();
            out << "@<TRIPOS>MOLECULE\n";
            if (data.moleculeName.empty())
                out << fileName;
            else
                out << data.moleculeName;
            out<< "\n" << nAtoms;
            string chargeType = "NO_CHARGES";
            if (!data.atomicCharge.empty())
                chargeType = "USER_CHARGES";
            out << " " << data.bonds.size() << "\nSMALL\n" << chargeType << "\n\n\n";

            // ---------------------------------------
            //           @<TRIPOS>ATOM

            out << "@<TRIPOS>ATOM\n";


            int longestLabelSize = 0;
            string label;
            vector<string> labels = data.atomName;
            if (labels.empty())
                for (int i = 0; i < nAtoms; i++)
                {
                    label = data.atomType[i] + to_string(i + 1);
                    if (label.size() > longestLabelSize)
                        longestLabelSize = label.size();
                    labels.push_back(label);
                }

            string atomType;
            bool hasCharge = !data.atomicCharge.empty();
            bool hasResidueName = !data.moleculeName.empty();
            bool hasResidueIdx = !data.substructureIdx.empty();
            int nAdditionalFields = 0;
            if (hasResidueIdx)
                nAdditionalFields = 1;
            if (hasResidueName)
                nAdditionalFields = 2;
            if (hasCharge)
                nAdditionalFields = 3;

            for (int i = 0; i < nAtoms; i++)
            {
                out << setw(7) << i + 1 << " ";

                out << setw(longestLabelSize + 2) << left << labels[i] << " " << right
                    << setw(14) << setprecision(6) << fixed << data.atomPosition[i].x << " "
                    << setw(14) << setprecision(6) << fixed << data.atomPosition[i].y << " "
                    << setw(14) << setprecision(6) << fixed << data.atomPosition[i].z << " " << data.atomType[i];

                if (nAdditionalFields > 0)
                {
                    if (hasResidueIdx)
                        out << " " << data.substructureIdx[i];
                    else
                        out << " 1";
                }
                if (nAdditionalFields > 1)
                {
                    if (hasResidueName)
                        out << " " << data.moleculeName[i];
                    else
                        out << " res1";
                }
                if (nAdditionalFields > 2)
                {
                    if (hasCharge)
                        out << " " << data.atomicCharge[i];
                    else
                        out << " 0";
                }


                out << "\n";
            }

            // ---------------------------------------
            //           @<TRIPOS>BOND

            if (!data.bonds.empty())
            {
                out << "@<TRIPOS>BOND\n";
                for (int i = 0; i < data.bonds.size(); i++)
                {
                    out << setw(7) << i + 1 << " " << setw(8) << data.bonds[i].first + 1 << " " << setw(8)
                        << data.bonds[i].second + 1;
                    if (data.bondTypes.empty())
                        out << "  1\n";
                    else
                        out << " " << data.bondTypes[i] << "\n";
                }
            }
            out.close();
        }


        

        void write(
            const std::string fileName,
            const std::vector<int>& atomicNumber,
            const std::vector<Vector3d>& position,
            bool detectBonds,
            const std::vector<std::string>& _atomLabel)
        {
            Mol2Data data;


            //ofstream out(fileName);


            //if (!out.good())
            //    on_error::throwException("cannot write mol2 file '" + fileName + "'", __FILE__, __LINE__);

            // ---------------------------------------
            //           @<TRIPOS>MOLECULE 

            //out << "@<TRIPOS>MOLECULE\nsome molecule\n" << atomicNumber.size();
            vector<pair<int, int> > bonds;
            vector<vector<int> > connectivity;
            if (detectBonds)
            {
                GenericConnectivityAlgorithm<CovalentRadiousBondDetector> connectivityAlgorithm;
                connectivityAlgorithm.set("0.4");
                connectivityAlgorithm.calculateConnectivity(position, atomicNumber, connectivity);
                for (int i = 0; i < connectivity.size(); i++)
                    for (int j : connectivity[i])
                        if (i < j)
                            bonds.push_back({ i, j });

            }

            //out << " " << bonds.size() << "\nSMALL\nNO_CHARGES\n\n\n";

            // ---------------------------------------
            //           @<TRIPOS>ATOM

            //out << "@<TRIPOS>ATOM\n";

            vector<string> labels, atomTypes;
            int longestLabelSize = 0;
            string label;

            for (int i = 0; i < atomicNumber.size(); i++)
            {
                if (_atomLabel.empty())
                    label = periodic_table::symbol(atomicNumber[i]) + to_string(i + 1);
                else
                    label = _atomLabel[i];
                if (label.size() > longestLabelSize)
                    longestLabelSize = label.size();
                labels.push_back(label);
                atomTypes.push_back(periodic_table::symbol(atomicNumber[i]));
            }

            data.atomName = labels;
            data.atomPosition = position;
            data.atomType = atomTypes;
            data.bonds = bonds;
            write(fileName, vector<Mol2Data>(1, data));
            return;

            //    string atomType;
            //    for (int i = 0; i < atomicNumber.size(); i++)
            //    {
            //        out << setw(7) << i + 1 << " ";
            //        
            //        out << setw(longestLabelSize+2) << left << labels[i] << " " << right
            //            << setw(14) << setprecision(6) << fixed << position[i].x << " "
            //            << setw(14) << setprecision(6) << fixed << position[i].y << " "
            //            << setw(14) << setprecision(6) << fixed << position[i].z << " " << periodic_table::symbol(atomicNumber[i]) << "\n";
            //    }

            //    // ---------------------------------------
            //    //           @<TRIPOS>BOND
            //    
            //    if (!bonds.empty())
            //    {
            //        out << "@<TRIPOS>BOND\n";
            //        for (int i = 0; i < bonds.size(); i++)
            //            out << setw(7) << i + 1 << " " << setw(8) << bonds[i].first + 1 << " " << setw(8) << bonds[i].second + 1 << "  1\n";
            //    }
            //    out.close();
            //}

        }

        void write(
            const std::string fileName,
            const std::vector<ChemicalElement>& chemicalElement,
            const std::vector<Vector3d>& position,
            bool detectBonds,
            const std::vector<std::string>& atomLabel)
        {
            vector<int> atomicNumber;
            for (auto el : chemicalElement)
                atomicNumber.push_back(el.atomicNumber());
            
            write(fileName, atomicNumber, position, detectBonds, atomLabel);
        }

    }
}