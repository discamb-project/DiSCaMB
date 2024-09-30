#include "discamb/IO/hkl_io.h"
//#include "discamb/IO/xd_io.h"
#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/BasicUtilities/OnError.h"

#include <fstream>

using namespace std;

namespace{
    bool intepretCompact(const string &s, discamb::Vector3i &v, int &count)
    {
        vector<string> words;
        discamb::string_utilities::split(s, words, 'x');
        string vStr;
        if (words.empty())
            return false;
        if (words.size()>1)
        {
            count = atoi(words[0].c_str());
            vStr = words[1];
        }
        else
        {
            count = 1;
            vStr = words[0];
        }

        discamb::string_utilities::split(vStr, words);

        if (words.size() == 3)
            v.set(atoi(words[0].c_str()), atoi(words[1].c_str()), atoi(words[2].c_str()));
        if (words.size() == 2)
            v.set(0, atoi(words[0].c_str()), atoi(words[1].c_str()));
        if (words.size() == 1)
            v.set(0, 0, atoi(words[0].c_str()));
        return true;

    }

    void writeVec(ofstream &out, const discamb::Vector3i &vec)
    {
        if (vec[0] == 0 && vec[1] == 0)
            out<< vec[2] << endl;
        else
        {
            if(vec[0] == 0)
                out<< vec[1] << " " << vec[2] << endl;
            else 
                out<< vec[0] << " " << vec[1] << " " << vec[2] << endl;

        }
    }

}

        // FORMAT(3I4,2F8.2,I4) 
namespace discamb {

    namespace hkl_io
    {
/*            void readReflectionFile(
        const std::string &reflection_file,
        std::vector<Vector3i> &hkl,
        std::vector<std::complex<double> > &sf)
    {
        int i, n;
        double a, b;
        string s;
        ifstream in(reflection_file.c_str());

        if (!in.good())
            on_error::throwException(string("can not read reflection file ") + reflection_file, __FILE__, __LINE__);

        in >> s >> s >> s >> n;
        hkl.resize(n);
        sf.resize(n);

        for (i = 0; i<n; i++)
        {
            in >> hkl[i][0] >> hkl[i][1] >> hkl[i][2] >> a >> b;
            sf[i] = complex<double>(a, b);
        }

        in.close();
    }

*/
void readHklIndices(
    const char *fileName,
    vector<Vector3i> &hklIndices,
    int hIndexColumn)
{
    hklIndices.clear();

    ifstream in(fileName);
    if (!in.good())
        on_error::throwException(string("can not open hkl file:") + string(fileName), __FILE__, __LINE__);

    vector<string> words;
    string line;

    while (in.good())
    {
        getline(in, line);
        string_utilities::split(line, words);
        if (words.size() >= hIndexColumn + 2)
            hklIndices.push_back(Vector3i(atoi(words[hIndexColumn - 1].c_str()),
                                          atoi(words[hIndexColumn].c_str()),
                                          atoi(words[hIndexColumn + 1].c_str())));
    }
    in.close();
}
/*
void readCompactHklFile(const string &fileName, vector<Vector3i> &hkl)
{
    ifstream in(fileName.c_str());
    hkl.clear();
    string line;
    vector<string> words;
    Vector3i v;
    int i, count;
    bool firstLine = true;

    if (!in.good())
        on_error::throwException(string("can not read hkl file '") + fileName + string("'"), __FILE__, __LINE__);

    while (in.good())
    {
        getline(in, line);
        string_utilities::split(line, words);
        if (intepretCompact(line, v, count))
        {
            if (firstLine)
            {
                hkl.push_back(v);
                for (i = 1; i<count; i++)
                    hkl.push_back(hkl.back() + v);
                firstLine = false;
            }
            else
                for (i = 0; i<count; i++)
                    hkl.push_back(hkl.back() + v);
        }
    }
    in.close();
}
*/
/*
void writeCompactHklFile(
    const std::string &fileName,
    const std::vector<Vector3i> &hkl)
{
    int counter = 0;
    int same_vec_counter = 0;
    Vector3i previous_vec, vec, vec_old;
    ofstream out(fileName.c_str());
    int hklIdx, nHkl = hkl.size();

    for (hklIdx = 0; hklIdx<nHkl; hklIdx++)
    {
        previous_vec = vec;

        vec = hkl[hklIdx];
        if (counter == 0)
        {
            vec = hkl[hklIdx];
            counter = 1;
            same_vec_counter = 1;
        }
        else
        {
            vec = vec - vec_old;
            if (vec == previous_vec)
                same_vec_counter += 1;
            else
            {
                if (same_vec_counter != 1)
                    out << same_vec_counter << 'x';
                writeVec(out, previous_vec);
                same_vec_counter = 1;
            }

        }
        vec_old = hkl[hklIdx];

    }

    if (same_vec_counter != 1)
    {
        out << same_vec_counter << 'x';
        writeVec(out, previous_vec);
    }
    else
        writeVec(out, vec);
}
*/
        
        void readShelxHkl(
            const std::string &fileName,
            std::vector<Vector3i> &hkls,
            std::vector<double> &intensities,
            std::vector<double> &sigmas,
            std::vector<int> &batchNumbers,
            bool freeFormat)
        {
            // reads shelx hkl file
            try
            {
                ifstream in(fileName);
                vector<string> words;
                string line, s;
                Vector3i hkl;
                double sigma, intensity;
                int h, k, l, batchNumber;
                bool read000 = false;
                bool hasBatchNumber;
                bool firstLine = true;

                if (!in.good())
                    on_error::throwException(string("cannot read hkl file '") + fileName + string("'"), __FILE__, __LINE__);
                    

                hkls.clear();
                intensities.clear();
                sigmas.clear();

                while (in.good() && !read000)
                {
                    getline(in, line);

                    bool LineHasAlphabeticCharacter = false;
                    for (auto c : line)
                        if (isalpha(c))
                            LineHasAlphabeticCharacter = true;
                    if (LineHasAlphabeticCharacter)
                        break;

                    if (freeFormat)
                        string_utilities::split(line, words);

                    if (firstLine)
                    {
                        if (freeFormat)
                            hasBatchNumber = (words.size() == 6);
                        else
                            hasBatchNumber = (line.size() >= 32);
                        firstLine = false;
                    }



                    if (freeFormat)
                    {
                        if (words.empty())
                            continue;
                        h = stoi(words[0]);
                        k = stoi(words[1]);
                        l = stoi(words[2]);
                        intensity = stod(words[3]);
                        sigma = stod(words[4]);
                        if (hasBatchNumber)
                            batchNumber = stoi(words[5]);
                    }
                    else
                    {
                        if (line.size() < 28)
                            continue;
                        h = stoi(line.substr(0, 4));
                        k = stoi(line.substr(4, 4));
                        l = stoi(line.substr(8, 4));
                        intensity = stod(line.substr(12, 8));
                        sigma = stod(line.substr(20, 8));
                        if (hasBatchNumber)
                            batchNumber = stoi(line.substr(28, 4));
                    }

                    if (h == 0 && k == 0 && l == 0)
                        read000 = true;



                    hkls.push_back(Vector3i(h, k, l));
                    intensities.push_back(intensity);
                    sigmas.push_back(sigma);
                    if (hasBatchNumber)
                        batchNumbers.push_back(batchNumber);
                }
                in.close();
            }
            catch (...) {
                on_error::throwException(string("problem when reading hkl file"), __FILE__, __LINE__);
            }

            

        } //void readShelxHkl

        void writeShelxHkl(
            const std::string &fileName,
            const std::vector<Vector3i> &hkl,
            const std::vector<double> &intensities,
            const std::vector<double> &sigma,
            const std::vector<int> &batchNumber,
            bool freeFormat)
        {
            ofstream out(fileName);

            if (!out.good())
                on_error::throwException(string("cannot write to shelx reflection file '") + fileName 
                                         + string("'"), __FILE__, __LINE__);

            for (int i = 0; i < hkl.size(); i++)
            {
                if (freeFormat)
                {
                    out << hkl[i][0] << " " << hkl[i][1] << " " << hkl[i][2] << " "
                        << fixed << intensities[i] << " " << fixed << sigma[i];
                    if (!batchNumber.empty())
                        out << " " << batchNumber[i];
                    out << endl;
                }
                else
                {
                    out << setw(4) << hkl[i][0]
                        << setw(4) << hkl[i][1]
                        << setw(4) << hkl[i][2]
                        << fixed << setw(8) << setprecision(2) << intensities[i]
                        << fixed << setw(8) << setprecision(2) << sigma[i];
                    if (!batchNumber.empty())
                        out << setw(4) << batchNumber[i];
                    out << endl;
                }
                
            } 
            out.close();
        }// writeShelxHkl

 /*       void readHklSf(
            const std::string &fileName,
            std::vector<Vector3i> &hkl,
            std::vector<std::complex<double> > &sf)
        {
            hkl.clear();
            sf.clear();

            ifstream in(fileName);

            if (!in.good())
                on_error::throwException(string("can not open hklSf input file '") + fileName + string("'"), __FILE__, __LINE__);

            string line;
            vector<string> words;
            int nReflections = 0;
            while (in.good())
            {

                getline(in, line);

                string_utilities::split(line, words);
                if (words.size() >= 5)
                {
                    nReflections++;
                    hkl.resize(nReflections);
                    


                    hkl.back()[0] = stod(words[0]);
                    hkl.back()[1] = stod(words[1]);
                    hkl.back()[2] = stod(words[2]);
                    sf.push_back( { stod(words[3]), stod(words[4]) });
                }

            }
            in.close();

        }
        */
/*
        void writeHklSf(
            const std::string &fileName,
            const std::vector<Vector3i> &hkl,
            const std::vector<std::complex<double> > &sf)
        {
            ofstream out(fileName);

            if (!out.good())
                on_error::throwException(string("can not open output file '")+fileName+string("'"), __FILE__, __LINE__);

            for (int i = 0; i < hkl.size(); i++)
                out << setw(6) << hkl[i][0] << setw(6) << hkl[i][1] << setw(6) << hkl[i][2] << setw(14) << sf[i].real() << setw(14)
                << sf[i].imag() << "\n";
            out.close();
        }
*/

    }
}
