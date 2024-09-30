#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/MathUtilities/MathUtilities.h"
#include "discamb/BasicUtilities/OnError.h"

#include "json.hpp"

#include <iostream>
#include <fstream>

using namespace std;

namespace discamb {

namespace scattering_utilities
{

    std::unique_ptr<SfCalculator> scatteringFactorCalculatorFromJsonFile(
        const Crystal& crystal,
        const std::string& jsonFile)
    {
        ifstream jsonFileStream("aspher.json");
        nlohmann::json jsonData;

        if (jsonFileStream.good())
                jsonFileStream >> jsonData;
        else
        {
            on_error::throwException("can not read aspher.json file, expected to be present in the current directory", __FILE__, __LINE__);
            return std::unique_ptr<SfCalculator>(nullptr);
        }

    
        string format = "none";
    
        if (jsonData.find("form factor engine") != jsonData.end())
            format = "1.0";
    
        if (jsonData.find("model") != jsonData.end() || jsonData.find("template") != jsonData.end())
            format = "1.1";
    
        if (format == string("none"))
        {
            on_error::throwException("unrecognized aspher.json format, must contain \"form factor engine\" or \"model\" element", __FILE__, __LINE__);
            return std::unique_ptr<SfCalculator>(nullptr);
        }
    
        nlohmann::json formFactorEngineJSON;

        if (format == string("1.0"))
        {

            formFactorEngineJSON = *jsonData.find("form factor engine");

            if (formFactorEngineJSON.find("type") == formFactorEngineJSON.end())
                return std::unique_ptr<SfCalculator>(nullptr);
        }

        if (format == string("1.1"))
        {

            //if (jsonData.find("template") != jsonData.end())
              //  getJsonFromTemplate(jsonData["template"].get<string>(), jsonData);

            //mergeWithDefaults(jsonData, "defaults.json");

            string modelName = jsonData["model"].get<string>();
            jsonData.erase("model");
            formFactorEngineJSON["type"] = modelName;
            formFactorEngineJSON["data"] = jsonData;
        }


        return std::unique_ptr<SfCalculator>(SfCalculator::create(crystal, formFactorEngineJSON));
    }

    void centeringSfMultipliers(
        char centering, 
        const vector<Vector3i> &hkl, 
        vector<double> &multipliers,
        bool obverse)
    {
        int hklIdx, nHkl = hkl.size();
        int hMod2,kMod2,lMod2,sumMod2;
        multipliers.resize(nHkl);

        if(centering == 'P' || centering == 'p')
        {
            multipliers.assign(nHkl,1.0);
        }
        else if (centering == 'R' || centering == 'r')
        {
            if(obverse)
                for (hklIdx = 0; hklIdx<nHkl; hklIdx++)
                {
                    if(abs(-hkl[hklIdx].x + hkl[hklIdx].y + hkl[hklIdx].z) % 3 == 0)
                        multipliers[hklIdx] = 3.0;
                    else
                        multipliers[hklIdx] = 0.0;
                }
            else
                for (hklIdx = 0; hklIdx<nHkl; hklIdx++)
                {
                    if (abs(hkl[hklIdx].x - hkl[hklIdx].y + hkl[hklIdx].z) % 2 == 0)
                        multipliers[hklIdx] = 3.0;
                    else
                        multipliers[hklIdx] = 0.0;
                }
        }
        else if (centering == 'F' || centering == 'f')
        {
            for ( hklIdx = 0; hklIdx<nHkl; hklIdx++ )
            {
                hMod2 = abs(hkl[hklIdx].x % 2);
                kMod2 = abs(hkl[hklIdx].y % 2);
                lMod2 = abs(hkl[hklIdx].z % 2);

                if( hMod2 == kMod2 && hMod2 == lMod2 )
                    multipliers[hklIdx] = 4.0;
                else
                    multipliers[hklIdx] = 0.0;
                  
            }
        }
        else if (centering == 'H' || centering == 'h')
            on_error::throwException("reflection conditions for H centering not implemented",__FILE__,__LINE__);
        else 
        {
            for (hklIdx = 0; hklIdx<nHkl; hklIdx++)
            {
                if ( centering == 'I' || centering == 'i' )
                    sumMod2 = (hkl[hklIdx].x + hkl[hklIdx].y + hkl[hklIdx].z) % 2;
                else if ( centering == 'A' || centering == 'a')
                    sumMod2 = (hkl[hklIdx].y + hkl[hklIdx].z) % 2;
                else if ( centering == 'B' || centering == 'b')
                    sumMod2 = (hkl[hklIdx].x + hkl[hklIdx].z ) % 2;
                else if ( centering == 'C' || centering == 'c' )
                    sumMod2 = (hkl[hklIdx].x + hkl[hklIdx].y) % 2;

                sumMod2 == 0 ? multipliers[hklIdx] = 2.0 : multipliers[hklIdx] = 0.0;
            }
        }
        
    }

    double scaleFactorFcalc(
        const std::vector<double>& i_obs,
        const std::vector<double>& i_calc,
        const std::vector<double>& i_sigma,
        double a,
        double b)
    {
        double s;
        double oc = 0;
        double cc = 0;
        double w, p;
        int i, n = i_obs.size();

        for (i = 0; i < n; i++)
        {
            w = 1.0 / (i_sigma[i] * i_sigma[i]);
            oc += w * i_obs[i] * i_calc[i];
            cc += w * i_calc[i] * i_calc[i];
        }

        s = oc / cc;

        oc = 0;
        cc = 0;

        for (i = 0; i < n; i++)
        {
            p = (2 * s * i_calc[i] + std::max(i_obs[i], 0.0)) / 3.0;
            w = 1.0 / (i_sigma[i] * i_sigma[i] + a * a * p * p + b * p);
            oc += w * i_obs[i] * i_calc[i];
            cc += w * i_calc[i] * i_calc[i];
        }
        
        s = oc / cc;

        return s;
    }


    void combineScatteringFactorsSets(
        const vector< vector< complex< REAL > > > &scatteringFactorsSets,
        vector< complex< REAL > > &scatteringFactors)
    {
        int nHkl, nHklInSet, hklIdx, nDataSets;

        nDataSets = scatteringFactorsSets.size();

        // merge scattering factors

        nHkl = 0;
        for (int i = 0; i<nDataSets; i++)
            nHkl += scatteringFactorsSets[i].size();
        scatteringFactors.resize(nHkl);
        hklIdx=0;
        for (int i = 0; i<nDataSets; i++)
        {
            nHklInSet = scatteringFactorsSets[i].size();
            for(int j=0 ; j < nHklInSet ; j++ )
                scatteringFactors[hklIdx++] = scatteringFactorsSets[i][j];
        }
    }

     

    void merge_dTarget_dParameterSets(
        const std::vector<std::vector<TargetFunctionAtomicParamDerivatives> > &dTarget_dParamSets,
        std::vector<TargetFunctionAtomicParamDerivatives> &dTarget_dParam)
    {
        int nAdpComponents, atomIdx, nAtoms, setIdx, nSets = dTarget_dParamSets.size();
        
        if(nSets == 0)
        {
            on_error::throwException(" BUG(?) - Empty set of derivatives of target function w.r.t. structural parameters passed as parameter",
                                     __FILE__, __LINE__);
            return;
        }

        dTarget_dParam = dTarget_dParamSets[0];
        nAtoms = dTarget_dParamSets[0].size();
        
        
        for (setIdx = 1; setIdx < nSets; setIdx++)
        {
            for (atomIdx = 0; atomIdx<nAtoms; atomIdx++)
            {
                nAdpComponents = dTarget_dParam[atomIdx].adp_derivatives.size();
                for (int i = 0; i<nAdpComponents; i++)
                    dTarget_dParam[atomIdx].adp_derivatives[i] += dTarget_dParamSets[setIdx][atomIdx].adp_derivatives[i];

                dTarget_dParam[atomIdx].atomic_position_derivatives += dTarget_dParamSets[setIdx][atomIdx].atomic_position_derivatives;
                dTarget_dParam[atomIdx].occupancy_derivatives += dTarget_dParamSets[setIdx][atomIdx].occupancy_derivatives;
            }

        }

    }
    
}

} // namespace discamb