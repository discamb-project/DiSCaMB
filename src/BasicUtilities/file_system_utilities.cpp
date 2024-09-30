#include "discamb/BasicUtilities/file_system_utilities.h"
#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb//MathUtilities/MixedRadixNumberIterator.h"

#include <filesystem>
#include <fstream>
#include <cstdlib>



using namespace std;


namespace discamb {

    namespace file_system_utilities {

        void file2string(
            const std::string& fileName,
            std::string& s)
        {
            s.clear();
            ifstream in(fileName,ios::binary);
            if (!in.good())
                on_error::throwException("cannot read file " + fileName, __FILE__, __LINE__);

            std::stringstream buffer;
            buffer << in.rdbuf();
            s = buffer.str();
            in.close();
        }

        void d2u(const std::string& fName)
        {
            d2u(fName, fName);
        }

        void d2u(
            const std::string& inputFileName,
            const std::string& outputFileName)
        {
            ifstream in(inputFileName, ios::binary);
            string s;
            char c, d;
            while (in.get(c))
            {
                if (c == '\r')
                {
                    if (in.get(d))
                    {
                        if (d == '\n')
                            s += '\n';
                        else
                        {
                            s += '\r';
                            s += d;
                        }
                    }
                    else
                        s += '\r';
                }
                else
                    s += c;
            }
            in.close();
            ofstream out(outputFileName, ios::binary);
            out << s;
            out.close();

        }

        void u2d(
            const std::string& fName)
        {
            u2d(fName, fName);
        }

        void u2d(
            const std::string& inputFileName,
            const std::string& outputFileName)
        {
            ifstream in(inputFileName, ios::binary);
            string s;
            char c;
            while (in.get(c))
            {
                if (c == '\n')
                    s += '\r';
                s += '\n';
            }
            in.close();
            ofstream out(outputFileName, ios::binary);
            out << s;
            out.close();
        }



        std::string find_first_file(
            const std::string& extension,
            bool caseSensitive)
        {
            vector<string> files;
            find_files(extension, files, caseSensitive);
           
            if (files.empty())
                return string();
            return *min_element(files.begin(), files.end());
        }

        std::string find_newest_file(
            const std::string& extension,
            bool caseSensitive)
        {
            vector<string> files;
            find_files(extension, files, caseSensitive);
            if (files.empty())
                return string();

            vector<pair<filesystem::file_time_type, string> > time_and_file;
            for (auto file : files)
                time_and_file.push_back({ filesystem::last_write_time(filesystem::path(file)), file});
            
            return max_element(time_and_file.begin(), time_and_file.end())->second;
        }

        void find_files(
            const std::string& extension,
            const std::string& folder,
            std::vector<std::string>& file,
            bool caseSensitive)
        {
            filesystem::path currentPath = filesystem::current_path();
            filesystem::current_path(folder);
            find_files(extension, file, caseSensitive);
            filesystem::current_path(currentPath);
        }

        void find_files(
            const std::string& _requiredExtension,
            std::vector<std::string>& files,
            bool caseSensitive)
        {
            string requiredExtension = _requiredExtension;
            string extension;
            files.clear();

            if (caseSensitive)
                string_utilities::toLower(requiredExtension, requiredExtension);

            for (auto it : filesystem::directory_iterator(filesystem::current_path()))
                if (filesystem::is_regular_file(it.status()))
                {
                    extension = it.path().extension().string();
                    
                    if (extension.size() > 1)
                    {
                        extension = extension.substr(1);

                        if (caseSensitive)
                            string_utilities::toLower(extension, extension);

                        if (extension == requiredExtension)
                            files.push_back(it.path().filename().string());
                    }
                }
        }

        bool files_are_identical(
            const std::string& file1,
            const std::string& file2)
        {
            int n;
            return files_are_identical(file1, file2, n);
        }

        bool files_are_identical(
            const std::string& file1,
            const std::string& file2,
            int& first_differing_line)
        {
            ifstream in1(file1, ios::in);
            if (!in1.good())
            {
                in1.close();
                return false;
            }
            ifstream in2(file2, ios::in);
            if (!in2.good())
            {
                in2.close();
                return false;
            }
            string line1, line2;
            int lineIdx = 0;
            while (getline(in1, line1))
            {
                getline(in2, line2);
                lineIdx++;
                if (line1 != line2)
                {
                    in1.close();
                    in2.close();
                    first_differing_line = lineIdx;
                    return false;
                }
                if (line1.find("PARAMETER MODIFICATION DATE") != string::npos)
                {
                    getline(in1, line1);
                    getline(in2, line2);
                    lineIdx++;
                }
            }
            in1.close();
            in2.close();
            return true;
        }

        std::string make_new_file_name(
            const std::string& pathStr,
            const std::string& extension,
            int length)
        {
            if (length == 0)
                return string();

            filesystem::path p;
            if (pathStr.empty())
                p = filesystem::current_path();
            else
                p = pathStr;
            //vector<char> c; {'a', 'b', 'c'}
            // 97-122
            vector<int> value, radix(length, 122 - 97 + 1);
            MixedRadixNumberIterator it(radix);
            bool stop = false;
            filesystem::path candidatePath;
            while (!stop)
            {
                string name;
                it.get(value);
                for (int v : value)
                    name += char(v+97);
                if (!extension.empty())
                    name += "." + extension;
                candidatePath = p / name;
                if (!filesystem::exists(candidatePath))
                    return name;
                stop = !it.increment();
            }
            return string();
        }

        /**
        searches for the file
        (1) in the current directory
        (2) in PATH specified folders
        */
        bool findFile(
            const std::string& fName,
            std::string& filePath)
        {
            filePath.clear();
            filesystem::path p = filesystem::current_path() / fName;
            if (filesystem::exists(p))
                if (filesystem::is_regular_file(p))
                {
                    filePath = p.string();
                    return true;
                }

            string env = getenv("PATH");
            vector<string> folderNames;
            string_utilities::split(env, folderNames, ';');

            for (auto& folder : folderNames)
            {
                filesystem::path folderPath(folder);
                p = folderPath / fName;
                if (filesystem::exists(p))
                    if (filesystem::is_regular_file(p))
                    {
                        filePath = p.string();
                        return true;
                    }
            }

            return false;
        }


    }
}
