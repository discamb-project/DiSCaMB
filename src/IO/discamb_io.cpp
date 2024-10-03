#include "discamb/IO/discamb_io.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/UnitCellContent.h"

#include <fstream>
#include <algorithm>

using namespace std;

namespace discamb {

    namespace discamb_io {

        void read_qm_systems_HirshFrag(
            const std::string& file_name,
            std::vector<int>& charge,
            std::vector<int>& spinMultiplicity,
            std::vector<std::string>& systemLabels,
            std::vector< std::vector<std::pair<std::string, std::string> > >& clusters)
        {
            
            systemLabels.clear();
            clusters.clear();
            charge.clear();
            spinMultiplicity.clear();

            ifstream in(file_name);
            if (!in.good())
                on_error::throwException("cannot read file '" + file_name + "'", __FILE__, __LINE__);

            string line;
            vector<string> words, words2;
            getline(in, line);
            string_utilities::split(line, words);
            systemLabels.push_back(words[1]);
            charge.push_back(stoi(words[2]));
            spinMultiplicity.push_back(stoi(words[3]));
            clusters.resize(clusters.size() + 1);

            while (getline(in, line))
            {
                string_utilities::split(line, words);

                if (words.empty())
                    continue;

                if (words.size() == 4)
                {
                    clusters.resize(clusters.size() + 1);
                    systemLabels.push_back(words[1]);
                    charge.push_back(stoi(words[2]));
                    spinMultiplicity.push_back(stoi(words[3]));
                }
                if (words.size() == 2)
                {
                    if (words[0].find('@') != string::npos)
                    {
                        replace(words[0].begin(), words[0].end(), ';', ',');
                        replace(words[1].begin(), words[1].end(), ';', ',');
                        clusters.back().push_back({ words[0], words[1] });
                    }
                    else
                    {
                        if (words[0].find(';')!=string::npos)
                        {
                            string_utilities::split(words[0], words2, ';');
                            clusters.back().push_back({ words2[0], words2[1] });
                        }
                        else
                            clusters.back().push_back({ words[0], "X,Y,Z"});
                    }
                }
            }
            
        }

        void read_representatives(
            const string& fileName,
            const Crystal& crystal,
            const std::vector<std::string>& clustersLabels,
            const std::vector< std::vector<std::pair<std::string, std::string> > >& clusterAtoms,
            vector<vector<AtomRepresentativeInfo> >& representatives)
        {
            ifstream in(fileName);
            string line, crystalAtomLabel, clusterAtomLabel, symmOpStr, weight, transformType, transformSetting;
            vector<string> words;
            stringstream ssLine;
            int refAtomIdx, clusterAtomIdx, clusterIdx, nClusters;
            map<pair<string, string>, int> labelToClusterAtomIdxMap;
            map<string, int> labelToCrystalAtomIdxMap;
            SpaceGroupOperation sgOp;

            nClusters = clustersLabels.size();


            for (int i = 0; i < crystal.atoms.size(); i++)
                labelToCrystalAtomIdxMap[crystal.atoms[i].label] = i;

            representatives.resize(crystal.atoms.size());

            if (!in.good())
                on_error::throwException(string("cannot read 'qm atom use file' ") + fileName + string("'"), __FILE__, __LINE__);

            bool firstLine = true;
            bool isClusterLabelLine;

            UnitCellContent unitCellContent;
            unitCellContent.set(crystal);


            while (in.good())
            {
                getline(in, line);
                string_utilities::split(line, words);

                int nWords = words.size();

                if (nWords != 0)
                {
                    auto it = find(clustersLabels.begin(), clustersLabels.end(), words[0]);
                    isClusterLabelLine = (it != clustersLabels.end());

                    if (isClusterLabelLine)
                    {
                        clusterIdx = distance(clustersLabels.begin(), it);

                        labelToClusterAtomIdxMap.clear();

                        for (int i = 0; i < clusterAtoms[clusterIdx].size(); i++)
                        {
                            sgOp.set(clusterAtoms[clusterIdx][i].second);
                            sgOp.get(symmOpStr);
                            labelToClusterAtomIdxMap[{clusterAtoms[clusterIdx][i].first, symmOpStr}] = i;
                        }
                    }
                    else
                    {
                        if (firstLine) // old format(?)
                        {
                            clusterIdx = 0;


                            for (int i = 0; i < clusterAtoms[clusterIdx].size(); i++)
                            {
                                sgOp.set(clusterAtoms[clusterIdx][i].second);
                                sgOp.get(symmOpStr);
                                labelToClusterAtomIdxMap[{clusterAtoms[clusterIdx][i].first, symmOpStr}] = i;
                            }
                        }


                        // --------

                        bool useNewImp = false;

                        // ------- NEW IMP

                        if (useNewImp)
                        {

                            if (nWords >= 5)
                            {
                                stringstream ssLine(line);
                                ssLine >> clusterAtomLabel >> symmOpStr >> crystalAtomLabel >> weight >> transformType;
                                SpaceGroupOperation sgOp(symmOpStr);
                                sgOp.get(symmOpStr);
                                //line.clear();
                                getline(ssLine, transformSetting);
                                //ssLine.str(string());
                            }
                            else
                            {
                                clusterAtomLabel = words[0];
                                if (nWords < 2)
                                    symmOpStr = string("X,Y,Z");
                                if (nWords < 3)
                                    crystalAtomLabel = clusterAtomLabel;
                                if (nWords < 4)
                                    weight = "1.0";
                                if (nWords < 5)
                                {
                                    transformType = "symmetry";
                                    transformSetting = string("invert ") + symmOpStr;
                                }
                            }

                            if (labelToCrystalAtomIdxMap.count(clusterAtomLabel) == 0)
                                on_error::throwException(string("unknown atom label '") + clusterAtomLabel +
                                    string("' in atom representatives file '") + fileName + string("'"), __FILE__, __LINE__);


                            if (labelToClusterAtomIdxMap.find({ clusterAtomLabel,symmOpStr }) == labelToClusterAtomIdxMap.end())
                                on_error::throwException(string("invalid atom specification in 'qm atom use file': '") + clusterAtomLabel + string(" ") + symmOpStr + string("'"), __FILE__, __LINE__);
                            clusterAtomIdx = labelToClusterAtomIdxMap[{clusterAtomLabel, symmOpStr}];


                            refAtomIdx = labelToCrystalAtomIdxMap[crystalAtomLabel];


                            representatives[refAtomIdx].resize(representatives[refAtomIdx].size() + 1);
                            representatives[refAtomIdx].back().idxInSubsystem = clusterAtomIdx;
                            representatives[refAtomIdx].back().fragmentIdx = clusterIdx;
                            representatives[refAtomIdx].back().atomLabel = crystalAtomLabel;
                            representatives[refAtomIdx].back().symmetryCode = symmOpStr;

                            representatives[refAtomIdx].back().isWeightFixed = isdigit(weight[0]);

                            if (representatives[refAtomIdx].back().isWeightFixed)
                                representatives[refAtomIdx].back().fixedWeightValue = stod(weight);
                            else
                                representatives[refAtomIdx].back().weightRepresentativeLabel = weight;

                            representatives[refAtomIdx].back().transformType = transformType;
                            representatives[refAtomIdx].back().transformSetting = transformSetting;

                        }
                        else
                        {
                            // ------- OLD IMP

                            if (nWords >= 5)
                            {
                                stringstream ssLine(line);
                                ssLine >> clusterAtomLabel >> symmOpStr >> crystalAtomLabel >> weight >> transformType;
                                sgOp.set(symmOpStr);
                                sgOp.get(symmOpStr);
                                getline(ssLine, transformSetting);

                            }
                            else
                            {
                                clusterAtomLabel = words[0];
                                words.resize(4);


                                symmOpStr = words[1];
                                if (nWords < 2)
                                    symmOpStr = string("X,Y,Z");

                                UnitCellContent::AtomID atomId;
                                unitCellContent.findAtom(clusterAtomLabel, symmOpStr, atomId);
                                unitCellContent.interpreteAtomID(atomId, clusterAtomLabel, symmOpStr);

                                crystalAtomLabel = words[2];
                                weight = words[3];


                                SpaceGroupOperation sgOp(symmOpStr);
                                sgOp.get(symmOpStr);

                                if (nWords < 3)
                                    crystalAtomLabel = clusterAtomLabel;
                                if (nWords < 4)
                                    weight = "1.0";
                                if (nWords < 5)
                                {
                                    transformType = "symmetry";
                                    transformSetting = string("invert ") + symmOpStr;
                                }
                            }

                            if (labelToCrystalAtomIdxMap.count(clusterAtomLabel) == 0)
                                on_error::throwException(string("unknown atom label '") + clusterAtomLabel +
                                    string("' in atom representatives file '") + fileName + string("'"), __FILE__, __LINE__);

                            if (labelToClusterAtomIdxMap.find({ clusterAtomLabel,symmOpStr }) == labelToClusterAtomIdxMap.end())
                                on_error::throwException("invalid atom specification in 'qm atom use file': '" + clusterAtomLabel + " " + symmOpStr + "'", __FILE__, __LINE__);
                            clusterAtomIdx = labelToClusterAtomIdxMap[{clusterAtomLabel, symmOpStr}];

                            if (labelToCrystalAtomIdxMap.find(crystalAtomLabel) == labelToCrystalAtomIdxMap.end())
                                on_error::throwException("invalid atom specification in 'qm atom use file': '" + crystalAtomLabel + "'", __FILE__, __LINE__);

                            refAtomIdx = labelToCrystalAtomIdxMap[crystalAtomLabel];


                            representatives[refAtomIdx].resize(representatives[refAtomIdx].size() + 1);
                            representatives[refAtomIdx].back().idxInSubsystem = clusterAtomIdx;
                            representatives[refAtomIdx].back().fragmentIdx = clusterIdx;
                            representatives[refAtomIdx].back().atomLabel = crystalAtomLabel;
                            representatives[refAtomIdx].back().symmetryCode = symmOpStr;

                            representatives[refAtomIdx].back().isWeightFixed = isdigit(weight[0]);

                            if (representatives[refAtomIdx].back().isWeightFixed)
                                representatives[refAtomIdx].back().fixedWeightValue = stod(weight);
                            else
                                representatives[refAtomIdx].back().weightRepresentativeLabel = weight;

                            representatives[refAtomIdx].back().transformType = transformType;
                            representatives[refAtomIdx].back().transformSetting = transformSetting;

                        }
                    }

                }

                firstLine = false;

            } //in.good()

            in.close();

        }

    }
}