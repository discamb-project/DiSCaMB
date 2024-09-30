#include "discamb/IO/proatom_db_io.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"

#include <fstream>

using namespace std;

namespace discamb
{
    namespace proatom_db_io
    {
        void read_proatom_db(
            const std::string& fileName,
            std::vector<int>& atomicNumber,
            std::vector<int>& charge,
            std::vector<std::vector<double> >& data)
        {
            ifstream in(fileName);
            if (!in.good())
                on_error::throwException(string("cannot read proatom database file '") + fileName + string("'"), __FILE__, __LINE__);
            
            string line;
            vector<string> words;
            int z;
            int q;
            vector<double> v(10000);

            while (in.good())
            {
                getline(in, line);
                string_utilities::split(line, words);
                if (words.size() == 4)
                {
                    z = stoi(words[1]);
                    q = stoi(words[3]);
                    for (int i = 0; i < 10000; i++)
                        in >> v[i];
                    atomicNumber.push_back(z);
                    charge.push_back(q);
                    data.push_back(v);
                }
            }
            in.close();
        }
    }
}
