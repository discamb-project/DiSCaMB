#include "discamb/QuantumChemistry/WaveFunctionCalculationData.h"
#include "discamb/BasicUtilities/file_system_utilities.h"

#include "discamb/BasicChemistry/PeriodicTable.h"
#include "discamb/BasicUtilities/StringUtilities.h"

using namespace std;

namespace discamb {

    void QmSettings::set(
        const nlohmann::json& data)
    {
        //qmMethod = data.value("qm method", "B3LYP");
        //relativisticMethod = data.value("relativistic method", data.value("rel", false) ? "DKH2" : "");
        //basisSet = data.value("basis set", "cc-pVDZ");

        qmMethod = data.value("qm method", string());
        relativisticMethod = data.value("relativistic method", data.value("rel", false) ? "DKH2" : "");
        basisSet = data.value("basis set", string());


        atomicNumber2BasisSetMap.clear();
        vector<string> words, words2;
        if (data.find("basis set per element") != data.end())
        {
            string basisSetPerElement = data["basis set per element"].get<string>();
            string_utilities::split(basisSetPerElement, words);
            for (string& word : words)
            {
                string_utilities::split(word, words2, ',');
                for (int i = 0; i < words2.size() - 1; i++)
                    atomicNumber2BasisSetMap[periodic_table::atomicNumber(words2[i])] = words2.back();
            }
        }

        
        inputTemplate = data.value("qm template", string());
        if (data.find("qm template file") != data.end())
            file_system_utilities::file2string(data["qm template file"].get<string>(), inputTemplate);
        tryToReadGuess = data.value("read wfn guess", false);
    }
    
}

