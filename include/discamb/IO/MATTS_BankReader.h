#pragma once


#include "discamb/AtomTyping/AtomType.h"
#include "discamb/Scattering/AtomTypeHC_Parameters.h"
#include "discamb/AtomTyping/StructureWithDescriptors.h"
#include "discamb/HC_Model/SlaterOrbitalWfnData.h"



#include <string>
#include <istream>
#include <map>

namespace discamb {

    /**
    * \addtogroup IO IO
    * @{
    */
    struct BankSettings {
        DescriptorsSettings descriptorsSettings;
        double min_plm = 0.002;
        double nSigma = 1.0;
        int min_n_instaces = 3;
        SlaterOrbitalWfnData::WfnDataBank wfn_databank = SlaterOrbitalWfnData::WfnDataBank::CR;
    };

    class MATTS_BankReader
    {
    public:
        MATTS_BankReader();
        ~MATTS_BankReader();

        enum class Section {
            COMMENT, NOI, ATOM_DESCRIPTORS, LOCAL_COORDINATE_SYSTEM, SYMMETRY, MULTIPOLE_MODEL_PARAMETERS, CHIRALITY,
            ID, PARAMETER_MODIFCATION_DATE, ENTRY_CREATION_DATE, DEFINITION_MODIFICATION_DATE, NONE, BLANK_LINE,
            END};

        struct RequiredSections {
            RequiredSections() {};
            bool undefined = true;
            std::vector<Section>  sections;
        };


        void read(const std::string &filename, std::vector<AtomType> &types, std::vector<AtomTypeHC_Parameters> &parameters,
                  BankSettings &bankSettings, bool addSphericalTypes = false, 
                  RequiredSections requiredSections = MATTS_BankReader::RequiredSections()) const;

        void read(const std::string &filename, std::vector<AtomType> &types, std::vector<AtomTypeHC_Parameters> &parameters,
                  std::vector<std::string> &lines, std::vector<std::pair<int, int> > &parameter_lines, BankSettings &bankSettings, 
                  bool addSphericalTypes = false, RequiredSections requiredSections = RequiredSections()) const;

        void read(std::istream &in, std::vector<AtomType> &types, std::vector<AtomTypeHC_Parameters> &parameters,
                  BankSettings &bankSettings, bool addSphericalTypes = false, 
                  RequiredSections requiredSections = RequiredSections()) const;

        void read(std::istream &in, std::vector<AtomType> &types, std::vector<AtomTypeHC_Parameters> &parameters,
                  std::vector<std::string> &lines, std::vector<std::pair<int, int> > &parameter_lines, 
                  BankSettings &bankSettings, bool addSphericalTypes = false,
                  RequiredSections requiredSections = RequiredSections()) const;


        bool readMultipoleParameters(std::istream &in, AtomTypeHC_Parameters &parameters, int &lastReadLineNumber) const;
        bool readMultipoleParameters(std::istream &in, AtomTypeHC_Parameters &parameters, int &lastReadLineNumber, std::vector<std::string> &lines) const;
        
        void readType(
            std::istream &in,
            AtomType &type, 
            int &lastReadLineNumber, 
            bool& readMultipolesLine, 
            const std::map<std::string, std::set<int> >& namedElementSets,
            RequiredSections requiredSections = RequiredSections()) const;
        
        void readType(
            std::istream &in,
            AtomType &type, 
            int &lastReadLineNumber,
            std::vector<std::string> &lines, 
            bool& readMultipolesLine, 
            const std::map<std::string, std::set<int> >& namedElementSets,
            RequiredSections requiredSections = RequiredSections()) const;
        
        void readSettings(std::istream &in, BankSettings &settings, int &lastReadLineNumber, std::vector<std::string> &lines) const;
        static std::string sectionToString(Section s);
        static void stringToAtomicNumbersSet(
            const std::string& s,
            std::set<int>& atomicNumbers,
            const std::map<std::string, std::set<int> >& namedSets);

    private:
        static std::map<Section, std::string> mSectionToKeywordMap;
        static std::map<Section, std::string> createSectionToKeywordMap();
        static std::map<std::string, Section> mKeywordToSectionMap;
        static std::map<std::string, Section> createKeywordToSectionMap();
        static void error_if_required(Section section, RequiredSections sections, int lineIdx);
        bool isNewSection(const std::string &line, Section &section) const;
        
        static void setSection(
            const std::vector<std::string> &lines, 
            Section section, 
            AtomType &type, 
            int firstLineNumber,
            const std::map<std::string, std::set<int> >& namedElementSets);
        
        static void setAtomDescriptors(
            const std::vector<std::string> &lines, 
            AtomType &type, 
            int firstLineNumber, 
            const std::map<std::string, std::set<int> >& namedElementSets);

        static void setChirality(const std::string &line, AtomType &type);
        /** negative for 'not of element' e.g. !H1 gives -1 */
        static void labelToAtomicNumber(
            const std::string &label, 
            std::set<int>& atomicNumbers, 
            const std::map<std::string, std::set<int> >& namedSets,
            bool &anyAtomicNumber);
    };


    /**@}*/

} //namespace

