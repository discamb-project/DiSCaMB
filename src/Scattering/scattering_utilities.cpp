#include "discamb/Scattering/scattering_utilities.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/Timer.h"

#include "json.hpp"


#include <fstream>

using namespace std;

namespace discamb {

namespace scattering_utilities
{

    void generate_symmetry_equivalent_hkls(
        const SpaceGroup& spaceGroup,
        const std::vector<Vector3i>& hkl,
        std::vector<Vector3i>& symmEquivalentHkls)
    {
        int nSymm = spaceGroup.nSymmetryOperationsInSubset();
        vector<Matrix3i> rotations(nSymm);

        for (int i = 0; i < nSymm; i++)
            spaceGroup.getSpaceGroupOperation(0, 0, i).getRotation(rotations[i]);
        for (int i = 0; i < nSymm; i++)
            rotations.push_back(-1 * rotations[i]);


        symmEquivalentHkls.clear();
        set<Vector3i> uniqueHkl;

        for (auto const& rotation : rotations)
            for (auto const& h : hkl)
                uniqueHkl.insert(h * rotation);

        symmEquivalentHkls.assign(uniqueHkl.begin(), uniqueHkl.end());

    }


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
    
    //int  findPreferredHklOrderingDirection(
    //    const std::vector<Vector3i>& hkl,
    //    std::vector<std::vector<Vector3i> >& orderedHklLines,
    //    std::vector<std::vector<int> >& mapToOriginalSetIndices)
    //{
    //    WallClockTimer timer;

    //    timer.start();

    //    orderedHklLines.clear();

    //    set<vector<int> > lines[3];
    //    for (auto const& h : hkl)
    //    {
    //        lines[0].insert({ h[1], h[2] });
    //        lines[1].insert({ h[0], h[2] });
    //        lines[2].insert({ h[0], h[1] });
    //    }

    //    cout << "step 4.1.1 time = " << timer.stop() << "\n";
    //    timer.start();


    //    int preferredDirection = 0;
    //    int nLines = lines[0].size();
    //    for (int i = 1; i < 3; i++)
    //    {
    //        if (nLines > lines[i].size())
    //        {
    //            preferredDirection = i;
    //            nLines = lines[i].size();
    //        }
    //    }

    //    cout << "step 4.1.2 time = " << timer.stop() << "\n";
    //    timer.start();


    //    map<vector<int>, int> line2idx;
    //    int idx = 0;
    //    for (auto& pq : lines[preferredDirection])
    //        line2idx[pq] = idx++;

    //    cout << "step 4.1.3 time = " << timer.stop() << "\n";
    //    timer.start();


    //    orderedHklLines.resize(nLines);
    //    vector<int> otherIndices;
    //    if (preferredDirection == 0)
    //        otherIndices = { 1, 2 };
    //    else if (preferredDirection == 1)
    //        otherIndices = { 0, 2 };
    //    else
    //        otherIndices = { 0, 1 };

    //    vector<vector<pair<Vector3i, int> > > orderedHklLinesWithOrgIndices(nLines);
    //    
    //    int nHkl = hkl.size();

    //    for(int hklIdx = 0; hklIdx<nHkl; hklIdx++)
    //    {
    //        auto const h = hkl[hklIdx];
    //        vector<int> pq = { h[otherIndices[0]], h[otherIndices[1]] };
    //        orderedHklLines[line2idx[pq]].push_back(h);
    //    }

    //    //for (auto const& h : hkl)
    //    //{
    //    //    vector<int> pq = { h[otherIndices[0]], h[otherIndices[1]] };
    //    //    orderedHklLines[line2idx[pq]].push_back(h);
    //    //}

    //    cout << "step 4.1.2 time = " << timer.stop() << "\n";
    //    timer.start();


    //    for (auto& line : orderedHklLines)
    //        std::sort(line.begin(), line.end());

    //    cout << "step 4.1.3 time = " << timer.stop() << "\n";
    //    timer.start();


    //    mapToOriginalSetIndices.clear();
    //    mapToOriginalSetIndices.resize(nLines);
    //    for (int i = 0; i < nLines; i++)
    //    {
    //        for (auto const& h : orderedHklLines[i])
    //        {
    //            auto it = find(hkl.begin(), hkl.end(), h);
    //            int idx = distance(hkl.begin(), it);
    //            mapToOriginalSetIndices[i].push_back(idx);
    //        }
    //    }

    //    cout << "step 4.1.4 time = " << timer.stop() << "\n";


    //    return preferredDirection;
    //}

    int  findPreferredHklOrderingDirection(
        const std::vector<Vector3i>& hkl,
        std::vector<std::vector<Vector3i> >& orderedHklLines,
        std::vector<std::vector<int> >& mapToOriginalSetIndices)
    {
        orderedHklLines.clear();

        set<pair<int, int> > lines[3];
        for (auto const& h : hkl)
        {
            lines[0].insert({ h[1], h[2] });
            lines[1].insert({ h[0], h[2] });
            lines[2].insert({ h[0], h[1] });
        }

        int preferredDirection = 0;
        int nLines = lines[0].size();
        for (int i = 1; i < 3; i++)
        {
            if (nLines > lines[i].size())
            {
                preferredDirection = i;
                nLines = lines[i].size();
            }
        }

        map<pair<int, int>, int> line2idx;
        int idx = 0;
        for (auto& pq : lines[preferredDirection])
            line2idx[pq] = idx++;

        orderedHklLines.resize(nLines);
        vector<int> otherIndices;
        if (preferredDirection == 0)
            otherIndices = { 1, 2 };
        else if (preferredDirection == 1)
            otherIndices = { 0, 2 };
        else
            otherIndices = { 0, 1 };

        vector<vector<pair<Vector3i, int> > > orderedHklLinesWithOrgIndices(nLines);

        int nHkl = hkl.size();

        for (int hklIdx = 0; hklIdx < nHkl; hklIdx++)
        {
            auto const h = hkl[hklIdx];
            pair<int, int> pq = { h[otherIndices[0]], h[otherIndices[1]] };
            int lineIdx = line2idx[pq];
            orderedHklLines[lineIdx].push_back(h);
            orderedHklLinesWithOrgIndices[lineIdx].push_back({ h, hklIdx });
        }


        for (auto& line : orderedHklLines)
            std::sort(line.begin(), line.end());

        mapToOriginalSetIndices.clear();
        mapToOriginalSetIndices.resize(nLines);
        for (int i = 0; i < nLines; i++)
            for (auto const& h_idx : orderedHklLinesWithOrgIndices[i])
                mapToOriginalSetIndices[i].push_back(h_idx.second);

        return preferredDirection;
    }


    void splitHklLines(
        int nSets,
        const std::vector<std::vector<Vector3i> >& orderedHklLines,
        std::vector<std::vector<int> >& lineGroups)
    {
        lineGroups.clear();
        lineGroups.resize(nSets);
        vector<pair<int, int> > lineSizeIdx;
        vector<int> setSizes(nSets, 0);

        for (int lineIdx = 0; lineIdx < orderedHklLines.size(); lineIdx++)
        {
            int lineSize = orderedHklLines[lineIdx].size();
            lineSizeIdx.push_back(make_pair(lineSize, lineIdx));
        }

        sort(lineSizeIdx.begin(), lineSizeIdx.end());
        
        for (int i = lineSizeIdx.size() - 1; i >= 0; i--)
        {
            int smallestSetIdx = distance(setSizes.begin(), min(setSizes.begin(), setSizes.end()));
            lineGroups[smallestSetIdx].push_back(i);
        }
    }

    void splitHklLines(
        int nSubstes,
        const std::vector<std::vector<Vector3i> >& orderedHklLines,
        const std::vector<std::vector<int> > &hklIdxInOriginalSet,
        std::vector < std::vector<std::vector<Vector3i> > > & newHklLines,
        std::vector < std::vector <std::vector<std::pair<int, int> > > > & hklMap,
        std::vector < std::vector<std::vector<int> > >& subsetDataHklIdxInOryginalSet)
    {
        newHklLines.clear();
        hklMap.clear();
        newHklLines.resize(nSubstes);
        hklMap.resize(nSubstes);
        subsetDataHklIdxInOryginalSet.resize(nSubstes);

        int lineIdx, nLines = orderedHklLines.size();
        int nHkl = 0;
        for (auto const& line : orderedHklLines)
            nHkl += line.size();
        int nHklPerLine = nHkl / nSubstes;
        int nExtrendedLines = nHkl - nHklPerLine * nSubstes;
        
        //merge lines
        vector<Vector3i> mergedHkl;
        vector<pair<int, int> > mergedIndices;
        for (lineIdx = 0; lineIdx < nLines; lineIdx++)
        {
            mergedHkl.insert(mergedHkl.end(), orderedHklLines[lineIdx].begin(), orderedHklLines[lineIdx].end());
            for (int hklIdx = 0; hklIdx < orderedHklLines[lineIdx].size(); hklIdx++)
                mergedIndices.push_back(make_pair(lineIdx, hklIdx));
        }

        int idx0 = 0;
        
        for (int setIdx = 0; setIdx < nSubstes; setIdx++)
        {
            int n = (setIdx < nExtrendedLines ? nHklPerLine + 1 : nHklPerLine);
            int lastInputLineIdx = -1;
            int idx = idx0;
            for ( ; idx < idx0 + n; idx++)
            {
                if (lastInputLineIdx != mergedIndices[idx].first)
                {
                    int currentSize = newHklLines[setIdx].size();
                    newHklLines[setIdx].resize(currentSize + 1);
                    hklMap[setIdx].resize(currentSize + 1);
                    subsetDataHklIdxInOryginalSet[setIdx].resize(currentSize + 1);
                }
                lastInputLineIdx = mergedIndices[idx].first;
                newHklLines[setIdx].back().push_back(mergedHkl[idx]);
                hklMap[setIdx].back().push_back(mergedIndices[idx]);
                subsetDataHklIdxInOryginalSet[setIdx].back().push_back(hklIdxInOriginalSet[mergedIndices[idx].first][mergedIndices[idx].second]);
            }
            idx0 = idx;
        }
    }


    void init_line_multipliers(
        const Vector3d& lineStepCart,
        const std::vector< std::vector<Vector3d> >& r_atom_symm,
        const std::vector<std::vector<double> >& adps,
        const std::vector<Matrix3d>& symmOpRotationCart,
        std::vector<std::vector<std::complex<double> > >& phase_factor_multiplier,
        std::vector<std::vector<double> >& temp_factor_multiplier_multiplier)
    {
        phase_factor_multiplier.clear();
        temp_factor_multiplier_multiplier.clear();

        int nAtoms = r_atom_symm.size();
        if (nAtoms == 0)
            return;
        int nSymmOps = r_atom_symm[0].size();
        double two_pi = 2 * M_PI;

        phase_factor_multiplier.resize(nAtoms, vector<complex<double> >(nSymmOps));
        temp_factor_multiplier_multiplier.resize(nAtoms, vector<double>(nSymmOps));

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
            {
                double phase_angle = two_pi * r_atom_symm[atomIdx][symOpIdx] * lineStepCart;
                phase_factor_multiplier[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
            }

        //temp_factor_multiplier_multiplier
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            vector<double> const & adp = adps[atomIdx];
            int nADPComponents = adp.size();
            if (nADPComponents == 1)
                temp_factor_multiplier_multiplier[atomIdx][0] = exp(-2.0 * lineStepCart * lineStepCart * adp[0]);
            else if (nADPComponents == 6)
            {

                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                {

                    Vector3d step_rot = lineStepCart * symmOpRotationCart[symOpIdx];
                    double exponent = step_rot[0] * step_rot[0] * adp[0] +
                        step_rot[1] * step_rot[1] * adp[1] +
                        step_rot[2] * step_rot[2] * adp[2] +
                        2.0 * (step_rot[0] * step_rot[1] * adp[3] +
                            step_rot[0] * step_rot[2] * adp[4] +
                            step_rot[1] * step_rot[2] * adp[5]);
                    temp_factor_multiplier_multiplier[atomIdx][symOpIdx] = exp(-2.0 * exponent);
                }
            }
        }

    }


    void calculate_line_phase_factors(
        //in:
        const std::vector<std::vector<Vector3d> >& r_at,
        const std::vector<Vector3d>& hkl_cart,
        int hkl_idx,
        const std::vector<Vector3i>& hkl_line,
        int idx_in_line,
        int line_direction,
        const std::vector<Vector3<double> >& rotated_h,
        const std::vector< std::vector<std::complex<double > > >& phase_factor_multiplier,
        //std::vector<Vector3d> &rotated_h_2pi,
        //std::vector<double> &translation_factor_2pi,

        //out:
        std::vector< std::vector<std::complex<double > > >& line_phase_factor)
    {
        int nAtoms = r_at.size();
        if (nAtoms == 0)
            return;
        int nSymmOps = r_at[0].size();
        double two_pi = 2 * M_PI;
        if (idx_in_line == 0)
        {
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                {
                    double phase_angle = two_pi * r_at[atomIdx][symOpIdx] * hkl_cart[hkl_idx];
                    line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                }
            }
        }
        else
        {


            int lineBreak = hkl_line[idx_in_line][line_direction] - hkl_line[idx_in_line - 1][line_direction] - 1;
            if (lineBreak == 0)
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                        line_phase_factor[atomIdx][symOpIdx] *= phase_factor_multiplier[atomIdx][symOpIdx];
            else
                for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        double phase_angle = two_pi * r_at[atomIdx][symOpIdx] * hkl_cart[hkl_idx];
                        //double phase_angle = rotated_h_2pi[symOpIdx] * mR_at[atomIdx] + translation_factor_2pi[symOpIdx];
                        line_phase_factor[atomIdx][symOpIdx] = { cos(phase_angle), sin(phase_angle) };
                    }

        }

    }

    void calculate_line_temperature_factors(
        //in:
        const std::vector<Vector3d> &coordinates,
        const std::vector< std::vector<double> >& atoms_adps,
        const std::vector<Vector3d>& hkl_cart,
        int hkl_idx,
        const std::vector<Vector3i>& hkl_line,
        int idx_in_line,
        int line_direction,
        const std::vector<Vector3<double> >& rotated_h,
        std::vector< std::vector<double> >& temp_factor_multiplier,
        const std::vector< std::vector<double> >& temp_factor_multiplier_multiplier,
        std::vector<std::vector<double> >& adpMultipliers,
        //out:
        std::vector< std::vector<double> >& line_temperature_factors)
    {
        int nAtoms = coordinates.size();
        int nSymmOps = rotated_h.size();
        double two_pi = 2 * M_PI;

        bool calc_tf_from_scratch; // first point after break/start or second point after break/start, but no continuation
        bool calc_tf_and_multiplier_from_scratch; // second point after break/start and there will be continuation
        bool calc_tf_from_multiplier; // third or further points afger start/break

        int points_after_break = 0;
        int points_before_break = 0;
        vector<int> line;
        for (int i = 0; i < hkl_line.size(); i++)
            line.push_back(hkl_line[i][line_direction]);

        for (int i = idx_in_line + 1; i < line.size(); i++)
        {
            if (line[i] == line[i - 1] + 1)
                points_before_break++;
            else
                break;
        }
        for (int i = idx_in_line - 1; i >= 0; i--)
        {
            if (line[i] == line[i + 1] - 1)
                points_after_break++;
            else
                break;
        }
        double const* adps;
        double hVectorLength2 = hkl_cart[hkl_idx] * hkl_cart[hkl_idx];

        // first point after break/start or second point after break/start, but no continuation
        // just calculate temperature factors from scratch
        if (points_after_break == 0 || (points_after_break == 1 && points_before_break == 0))
        {
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {

                int n_adp_components = atoms_adps[atomIdx].size();// atomic_displacement_parameters[atomIdx].size();

                if (n_adp_components > 0)
                    adps = &atoms_adps[atomIdx][0];

                if (n_adp_components == 6)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        double* multipliers = &adpMultipliers[symOpIdx][0];
                        line_temperature_factors[atomIdx][symOpIdx] = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                            - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                            - multipliers[4] * adps[4] - multipliers[5] * adps[5]);
                    }
                else if (n_adp_components == 1)
                    line_temperature_factors[atomIdx][0] = exp(-hVectorLength2 * (*adps));

            }
        }
        else if (points_after_break == 1 && points_before_break > 0)
        {
            // second point after break/start and there will be continuation

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {

                int n_adp_components = atoms_adps[atomIdx].size();

                if (n_adp_components > 0)
                    adps = &atoms_adps[atomIdx][0];

                if (n_adp_components == 6)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        double* multipliers = &adpMultipliers[symOpIdx][0];
                        double previous_temp_factor = line_temperature_factors[atomIdx][symOpIdx];
                        line_temperature_factors[atomIdx][symOpIdx] = exp(-multipliers[0] * adps[0] - multipliers[1] * adps[1]
                            - multipliers[2] * adps[2] - multipliers[3] * adps[3]
                            - multipliers[4] * adps[4] - multipliers[5] * adps[5]);
                        temp_factor_multiplier[atomIdx][symOpIdx] =
                            line_temperature_factors[atomIdx][symOpIdx] / previous_temp_factor;
                    }
                else if (n_adp_components == 1)
                {
                    double previous_temp_factor = line_temperature_factors[atomIdx][0];
                    line_temperature_factors[atomIdx][0] = exp(-hVectorLength2 * (*adps));
                    temp_factor_multiplier[atomIdx][0] =
                        line_temperature_factors[atomIdx][0] / previous_temp_factor;
                }

            }
        }
        else
        {
            // third or further points after start/break
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                int n_adp_components = atoms_adps[atomIdx].size();

                if (n_adp_components == 6)
                    for (int symOpIdx = 0; symOpIdx < nSymmOps; symOpIdx++)
                    {
                        temp_factor_multiplier[atomIdx][symOpIdx] *= temp_factor_multiplier_multiplier[atomIdx][symOpIdx];
                        line_temperature_factors[atomIdx][symOpIdx] *= temp_factor_multiplier[atomIdx][symOpIdx];
                    }
                else if (n_adp_components == 1)
                {
                    temp_factor_multiplier[atomIdx][0] *= temp_factor_multiplier_multiplier[atomIdx][0];
                    line_temperature_factors[atomIdx][0] *= temp_factor_multiplier[atomIdx][0];
                }

            }
        }

    }

}

} // namespace discamb