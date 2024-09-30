#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/string_utilities.h"

using namespace std;

namespace discamb {
    namespace basic_chemistry_utilities {
        
        int atomicNumberFromLabel(const string &label)
        {
            //OS = O + S or Os, HF = H + F or Hf, CN = C + N or Cn 
            int nCharacters = label.size();

            if (nCharacters == 0)
                return 0;
            if (nCharacters == 1)
                return periodic_table::atomicNumber(string(1,toupper(label[0])));

            if(isalpha(label[1]))
                return periodic_table::atomicNumber(string(1, toupper(label[0]))+string(1,tolower(label[1])));

            return periodic_table::atomicNumber(string(1, toupper(label[0])));
            
        }

        // s can be atomic number or symbol
        int atomicNumberFromString(const std::string& s)
        {
            return std::isdigit(s[0]) ? std::stoi(s) : periodic_table::atomicNumber(s);
        }

        std::string formulaAsString(const std::map<int, int>& formula)
        {
            string result;
            for (const auto& item : formula)
                result += periodic_table::symbol(item.first) + to_string(item.second);
            return result;
        }

        void getElementsList(
            const string& s,
            set<int>& atomic_numbers)
        {
            atomic_numbers.clear();

            vector<string> words, words2;
            string_utilities::split(s, words, ',');

            for (auto& word : words)
            {
                string_utilities::split(word, words2, '-');
                if (words2.size() > 1)
                {
                    int first = atomicNumberFromString(words2[0]);
                    int last = atomicNumberFromString(words2[1]);
                    for (int z = first; z <= last; z++)
                        atomic_numbers.insert(z);
                }
                else
                    atomic_numbers.insert(atomicNumberFromString(word));
            }
            
        }

        void getFormula(
            const std::vector<int>& atomic_numbers,
            std::map<int, int>& formula)
        {
            formula.clear();
            for (int z : atomic_numbers)
                if (formula.find(z) != formula.end())
                    formula[z]++;
                else
                    formula[z] = 1;
        }

    } //namespace basic_chemistry_utilities
} //namespace discamb

