#pragma once

#include "discamb/HC_Model/HC_ModelParameters.h"
#include "discamb/Scattering/AtomTypeHC_Parameters.h"

#include <iostream>
#include <complex>
#include <map>

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


class BankHcFunctions {
    public:
        BankHcFunctions();
        BankHcFunctions(const std::string &hcBankFileName, bool addSpherical = false);
        
        void set(const std::string& hcBankFileName, bool addSpherical = false);
        void set(std::istream& bankStream, bool addSpherical = false);

        int nTypes() const { return mTypeParameters.size(); }
        HC_ModelParameters const& getTypeParameters(int typeIdx) const { return mTypeParameters[typeIdx]; }
        int getTypeIdx(const std::string& typeLabel) const;
        std::string getTypeLabel(int idx) const;
        bool hasType(const std::string& typeLabel) const;
        void setDefault(bool addSpherical = false);
        double f_core(int typeIdx, double h) const;
        double f_val(int typeIdx, double h) const;
        double p_val(int typeIdx) const;
        double fourier_bessel_transform(int typeIdx, int l, double h) const;

    private:
        mutable std::vector<HC_ModelParameters> mTypeParameters;
        std::map<std::string, int> mTypeIdx;
        std::vector<std::string> mTypeLabel;
        //HC_ModelParameters mTypeParameters;

        void setHC_ModelParameters(
            const std::vector<AtomTypeHC_Parameters>& bankMultipoleParameters,
            int atomType,
            int atomicNumber,
            HC_ModelParameters& parameters);

        double calculate(HC_ModelParameters& hcParams, double h) const;
        std::complex<double> calculateComplex(HC_ModelParameters& hcParams, double h) const;
        

};
/** @}*/
}
