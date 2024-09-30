#include <vector>
#include <string>

namespace discamb {

    /**
    * \addtogroup BasicUtilities BasicUtilities
    * @{
    */


    namespace file_system_utilities {
        void file2string(const std::string& fileName, std::string& s);
        // not including ., eg 'cif' not '.cif'
        std::string find_first_file(const std::string& extension, bool caseSensitive=true);
        std::string find_newest_file(const std::string& extension, bool caseSensitive = true);
        void find_files(const std::string &extension, std::vector<std::string> &file, bool caseSensitive = true);
        void find_files(const std::string& extension, const std::string& folder, std::vector<std::string>& file, bool caseSensitive = true);
        std::string make_new_file_name(const std::string& path = std::string(), const std::string& extension = std::string(), int length = 8);
        /**
        searches for the file 
        (1) in the current directory
        (2) in PATH specified folders
        */
        bool findFile(const std::string& fName, std::string& filePath);

        bool files_are_identical(const std::string& file1, const std::string& file2);
        bool files_are_identical(const std::string& file1, const std::string& file2, int &first_differing_line);

        void d2u(const std::string& fName);
        void d2u(const std::string& inputFileName, const std::string& outputFileName);
        void u2d(const std::string& fName);
        void u2d(const std::string& inputFileName, const std::string& outputFileName);

    }
    /**@}*/
}
