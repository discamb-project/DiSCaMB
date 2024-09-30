#include "discamb/IO/rsp_io.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"

#include <fstream>

using namespace std;

namespace discamb {
    namespace rsp_io {
        void read(
            const std::string& fileName, 
            std::vector<Vector3d>& q)
        {
            q.clear();
            ifstream in(fileName);
            if (!in.good())
                on_error::throwException("cannot read file '" + fileName + "'", __FILE__, __LINE__);
            string line;
            getline(in, line);
            vector<string> words;

            if (line.find("POINTS LIST")!=string::npos)
            {
                while (getline(in, line))
                {
                    string_utilities::split(line, words);
                    if (words.size() == 3)
                    {
                        Vector3d v(stod(words[0]), stod(words[1]), stod(words[2]));
                        q.push_back(v);
                    }
                }
            }
            else if (line.find("POINTS GRID") != string::npos)
            {
                /*
start -3 0 0
size 601 1 1
step1 0.01  0.0 0.0
step2 0.0   0.0 0.0
step3 0.0   0.0 0.0
                */
                string s;
                Vector3d start, step1, step2, step3;
                int n1, n2, n3;
                in >> s >> start.x >> start.y >> start.z
                    >> s >> n1 >> n2 >> n3
                    >> s >> step1.x >> step1.y >> step1.z
                    >> s >> step2.x >> step2.y >> step2.z
                    >> s >> step3.x >> step3.y >> step3.z;
                for (int i = 0; i < n1; i++)
                    for (int j = 0; j < n2; j++)
                        for (int k = 0; k < n3; k++)
                            q.push_back(start + double(i) * step1 + double(j) * step2 + double(k) * step3);
            }
            else
                on_error::throwException("incorrect format of rsp file '" + fileName + "'", __FILE__, __LINE__);
            in.close();
        }
    }
}