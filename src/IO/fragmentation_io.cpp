#include "discamb/IO/fragmentation_io.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"

#include <fstream>
#include <set>

using namespace std;

namespace{


    void readClusters(
        std::istream& in,
        std::vector<discamb::FragmentConstructionData>& fragmentsConstructionData)
    {
        fragmentsConstructionData.clear();

        string line;
        vector<string> words;
        
        discamb::FragmentPartConstructionData fragmentPartConstructionData;

        while (getline(in, line))
        {
            discamb::string_utilities::split(line, words);
            int nWords = words.size();
            if (nWords == 0)
                continue;

            if (words[0][0] == '$' || nWords == 1 || nWords==2)
            {
                if (fragmentsConstructionData.empty())
                    discamb::on_error::throwException("invalid fragment definition line: '" + line + "'", __FILE__, __LINE__);
                fragmentPartConstructionData.set(words);
                fragmentsConstructionData.back().fragmentPartConstructionData.push_back(fragmentPartConstructionData);
            }
            else
            {
                if (words.size() == 4)
                {
                    fragmentsConstructionData.resize(fragmentsConstructionData.size() + 1);
                    fragmentsConstructionData.back().label = words[1];
                    fragmentsConstructionData.back().charge = stoi(words[2]);
                    fragmentsConstructionData.back().spin_multiplicity = stoi(words[3]);
                }
                else
                    discamb::on_error::throwException("invalid fragment definition line: '" + line + "'", __FILE__, __LINE__);
            }

        }
        // check if fragments/subsystems labels are unique
        
        std::set<string> uniqueLabels;
        
        for (auto& fragment : fragmentsConstructionData)
        {
            if (uniqueLabels.find(fragment.label) != uniqueLabels.end())
            {
                string errorMessage = string("non unique label '") + fragment.label +
                    string("' of subsystem/fragment for quantum mechanical calculations");
                discamb::on_error::throwException(errorMessage, __FILE__, __LINE__);
            }
            else
                uniqueLabels.insert(fragment.label);
        }

    }


    void replaceSemicolonWithComma(string& s)
    {
        while (s.find(';') != string::npos)
            s.replace(s.find(';'), 1, ",");
    }



    bool splitAtomAndSymmOp(
        const string& s,
        string& atomLabel,
        string& symmOp)
    {
        vector<string> words;
        discamb::string_utilities::split(s, words, ',');
        if (words.size() != 4)
            return false;

        atomLabel = words[0];
        symmOp = s.substr(words[0].size() + 1);
        return true;
    }

}

namespace discamb {

    namespace fragmentation_io {

        void parseCappingAtomInfo(
            const string word1,
            const string word2,
            string& bondedAtom,
            string& bondedAtomSymmOp,
            string& directingAtom,
            string& directingAtomSymmOp)
        {


            if (!splitAtomAndSymmOp(word1.substr(2), bondedAtom, bondedAtomSymmOp))
            {
                bondedAtom = word1.substr(2);
                bondedAtomSymmOp = "X,Y,Z";
            }

            if (!splitAtomAndSymmOp(word2, directingAtom, directingAtomSymmOp))
            {
                directingAtom = word2;
                directingAtomSymmOp = "X,Y,Z";
            }

        }

        void readPlainText(
            const std::string& fileName, std::vector<FragmentConstructionData>& fragmentsConstructionData)
        {
            ifstream in(fileName);
            if (!in.good())
                on_error::throwException("cannot read file definining wavefunction calculation units: '"
                    + fileName + "'", __FILE__, __LINE__);

            readClusters(in, fragmentsConstructionData);
            in.close();
        }

        void readJson(
            const std::string& fileName, 
            std::vector<FragmentConstructionData>& fragmentsConstructionData)
        {
            on_error::not_implemented(__FILE__, __LINE__);
            std::ifstream f(fileName);
            nlohmann::json data = nlohmann::json::parse(f);
            f.close();

        }



        void read(
            const std::string& fileName,
            std::vector<FragmentConstructionData>& fragments)
        {
            on_error::not_implemented(__FILE__, __LINE__);
            //fragments.clear();

            //ifstream in(fileName);
            //if (!in.good())
            //    on_error::throwException("cannot read file definining wavefunction calculation units: '"
            //        + fileName + "'", __FILE__, __LINE__);

            //string line;
            //bool hirshfrag = false;
            //getline(in, line);
            //if (line.find("Fragment") != string::npos)
            //    hirshfrag = true;
            //in.seekg(0);
            //
            //if (hirshfrag)
            //    readClustersHirshFrag(in, fragments);
            //else
            //    readClusters(in, fragments);
            //in.close();

        }

        void readHirsFrag(
            const std::string& fileName,
            std::vector<FragmentConstructionData>& fragments)
        {
            ifstream in(fileName);
            if (!in.good())
                on_error::throwException("cannot read file definining wavefunction calculation units: '"
                    + fileName + "'", __FILE__, __LINE__);

            stringstream ss;
            string line;
            vector<string> words, words2;
            while (in.good())
            {
                getline(in, line);
                discamb::string_utilities::split(line, words);
                if (words.empty())
                {
                    ss << line << "\n";
                    continue;
                }

                if (words[0] == "Fragment")
                {
                    ss << "System " << words[1] << " " << words[2] << " " << words[3] << "\n";
                    continue;
                }

                if (words[0].find('@') != string::npos)
                {
                    replaceSemicolonWithComma(words[0]);
                    replaceSemicolonWithComma(words[1]);
                    ss << words[0] << " " << words[1] << "\n";
                }
                else
                {
                    if (words[0].find(',') != string::npos)
                    {
                        discamb::string_utilities::split(words[0], words2, ',');
                        ss << words2[0] << " " << words2[1] << "," << words2[2] << "," << words2[3] << "\n";
                    }
                    else
                        ss << words[0] << "\n";
                }
            }
            readClusters(ss, fragments);
            in.close();
        }


        //void read(
        //    const std::string& fileName, 
        //    std::vector<FragmentData>& data)
        //{
        //    ifstream in(fileName);
        //    if(!in.good())
        //        on_error::throwException("cannot read file definining wavefunction calculation units: '"
        //            + fileName + "'", __FILE__, __LINE__);          

        //    data.clear();

        //    UnitCellContent ucContent;
        //    ucContent.set(crystal);
        //    vector<vector<UnitCellContent::AtomID> > connectivity;
        //    vector<vector<string> > atomsBasisSet;
        //    structural_properties::calcUnitCellConnectivity(ucContent, connectivity, 0.4);

        //    string line;
        //    vector<string> words;
        //    string currentBasisSet;
        //    bool usesClusterSpecificBasisSet = false;
        //    getline(in, line);
        //    string_utilities::split(line, words);
        //    if (words.size() != 2 && words.size() != 4)
        //        on_error::throwException(string("invalid format of qm systems file"), __FILE__, __LINE__);
        //    bool singleSystemFormat = (words.size() == 2);

        //    clusters.resize(1);
        //    atomsBasisSet.resize(1);
        //    if (singleSystemFormat)
        //    {
        //        clusters[0].push_back({ words[0],words[1] });
        //        systemLabels.push_back("1");
        //        charge.push_back(0);
        //        spinMultiplicity.push_back(1);
        //    }
        //    else
        //    {
        //        systemLabels.push_back(words[1]);
        //        charge.push_back(stoi(words[2]));
        //        spinMultiplicity.push_back(stoi(words[3]));
        //    }



        //    while (getline(in, line))
        //    {
        //        string_utilities::split(line, words);

        //        if (words.empty())
        //            continue;

        //        if (words[0][0] == '$')
        //        {
        //            if (words[0] == string("$basis_set"))
        //            {
        //                if (words.size() > 1)
        //                    currentBasisSet = words[1];
        //                else
        //                    currentBasisSet = string();
        //            }
        //            else
        //                processClusterCreationCommand2(
        //                    ucContent,
        //                    connectivity,
        //                    words,
        //                    currentBasisSet,
        //                    atomsBasisSet.back(),
        //                    clusters.back());
        //        }
        //        else
        //        {
        //            if (words.size() == 4)
        //            {
        //                clusters.resize(clusters.size() + 1);
        //                atomsBasisSet.resize(atomsBasisSet.size() + 1);
        //                systemLabels.push_back(words[1]);
        //                charge.push_back(stoi(words[2]));
        //                spinMultiplicity.push_back(stoi(words[3]));
        //                currentBasisSet.clear();

        //            }
        //            if (words.size() == 2)
        //            {
        //                clusters.back().push_back({ words[0], words[1] });
        //                atomsBasisSet.back().push_back(currentBasisSet);
        //            }
        //            if (words.size() == 1)
        //            {
        //                clusters.back().push_back({ words[0], string("X,Y,Z") });
        //                atomsBasisSet.back().push_back(currentBasisSet);
        //            }
        //        }

        //    }

        //    int clusterIdx, nClusters, atomIdx, nAtoms;
        //    nClusters = clusters.size();
        //    atomIdx2BasisSetMap.resize(nClusters);
        //    for (clusterIdx = 0; clusterIdx < nClusters; clusterIdx++)
        //    {
        //        nAtoms = clusters[clusterIdx].size();
        //        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        //            if (!atomsBasisSet[clusterIdx][atomIdx].empty())
        //                atomIdx2BasisSetMap[clusterIdx][atomIdx] = atomsBasisSet[clusterIdx][atomIdx];
        //    }

        //}
    }

}

