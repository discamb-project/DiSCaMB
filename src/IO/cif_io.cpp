#include "discamb/IO/cif_io.h"

#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/MathUtilities/math_utilities.h"

#include "builder.h"
#include "parser.h"

#include <fstream>
#include <iomanip>
#include <list>

using namespace std;

namespace {

    double get_precision(const string& _s)
    {
        // possible formats 12(3), 12, 0.12, 0.12(3), 12.12, 122.12(3)
        string s = _s;
        vector<string> words;
        discamb::string_utilities::split(s, words, '(');
        s = words[0];
        discamb::string_utilities::split(s, words, '.');
        if (words.size() == 1)
            return 1.0;

        s = string("0.") + string(words[1].size() - 1, '0') + string("1");
        return stod(s);
    }

    struct CifArrayWrapper : ucif::array_wrapper_base
    {
        std::vector<std::string> array;

        CifArrayWrapper()
            : array()
        {}

        virtual void push_back(std::string const& value)
        {
            array.push_back(value);
        }



        virtual std::string operator[](unsigned const& i) const
        {
            return array[i];
        }

        virtual unsigned size() const
        {
            return static_cast<unsigned int>(array.size());
        }
    };

    struct CifBuilder : ucif::builder_base
    {
        std::vector<discamb::cif_io::DataSet> dataSets;
        virtual void start_save_frame(std::string const& save_frame_heading) {}
        virtual void end_save_frame() {}

        virtual void add_data_item(std::string const& tag, std::string const& value)
        {
            dataSets.back().key_value_items[tag] = value;
        }

        //virtual void add_loop(ucif::array_wrapper_base const& loop_headers,
        //    ucif::array_wrapper_base const& values) {}

        virtual void add_loop(
            ucif::array_wrapper_base const& loop_headers,
            std::vector<ucif::array_wrapper_base*> const& values)
        {
            int nItems, nValues, itemIdx, valueIdx;
            dataSets.back().loops.resize(dataSets.back().loops.size() + 1);
            discamb::cif_io::Loop &loop = dataSets.back().loops.back();

            if (values.size() == 0)
                return;

            nItems = loop_headers.size();
            nValues = values[0]->size();
            loop.loop_headers.resize(nItems);
            loop.values.resize(nItems, std::vector<std::string>(nValues));

            for (itemIdx = 0; itemIdx < nItems; itemIdx++)
            {
                loop.loop_headers[itemIdx] = loop_headers[static_cast<unsigned int>(itemIdx)];
                for (valueIdx = 0; valueIdx < nValues; valueIdx++)
                    loop.values[itemIdx][valueIdx] = values[static_cast<unsigned int>(itemIdx)]->operator[](static_cast<unsigned int>(valueIdx));
            }

        }


        virtual void add_data_block(std::string const& data_block_heading)
        {
            dataSets.resize(dataSets.size() + 1);
            dataSets.back().dataSetName = data_block_heading;
        }

        virtual ucif::array_wrapper_base* new_array()
        {
            return new CifArrayWrapper();
        }
    };


    bool getAnisoData(
        const discamb::cif_io::DataSet &data,
        vector<int> &atomLoopToAnisoLoopEntryMap, 
        vector<vector<double> > &adp,
        vector<vector<double> > &adp_sigma,
        vector<vector<double> > &adp_precision)
    {
        int atomSiteLoopIdx, anisoLoopIdx, anisoLabelCol, atomLabelCol;

        atomLoopToAnisoLoopEntryMap.clear();
        adp.clear();
        adp_sigma.clear();
        
        if (!findLoopIndex(data, "_atom_site_label", atomSiteLoopIdx, atomLabelCol))
            return false;
        
        if (!findLoopIndex(data, "_atom_site_aniso_label", anisoLoopIdx, anisoLabelCol))
            return false;

        const discamb::cif_io::Loop &atomLoop = data.loops[atomSiteLoopIdx];
        const discamb::cif_io::Loop &adpLoop = data.loops[anisoLoopIdx];

        int nAtoms, nAnisoAtoms, atomIdx, anisoAtomIdx;

        nAtoms = atomLoop.values[atomLabelCol].size();
        nAnisoAtoms = adpLoop.values[anisoLabelCol].size();

        // has all U_ij or B_ij


        int adp11Col, adp22Col, adp33Col, adp12Col, adp13Col, adp23Col;
        bool hasAdp11Col, hasAdp22Col, hasAdp33Col, hasAdp12Col, hasAdp13Col, hasAdp23Col;
        bool hasUani, hasBani;
        hasUani = findLoopTagIndex(adpLoop, "_atom_site_aniso_U_11", adp11Col);
        hasBani = findLoopTagIndex(adpLoop, "_atom_site_aniso_B_11", adp11Col);

        if (!(hasUani || hasBani))
            return false;

        hasAdp11Col = true;

        // check if all U_ij components are given
        if (hasUani)
        {
            hasAdp22Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_U_22", adp22Col);
            hasAdp33Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_U_33", adp33Col);
            hasAdp12Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_U_12", adp12Col);
            hasAdp13Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_U_13", adp13Col);
            hasAdp23Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_U_23", adp23Col);
        }
        else        
        {
            hasAdp22Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_B_22", adp22Col);
            hasAdp33Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_B_33", adp33Col);
            hasAdp12Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_B_12", adp12Col);
            hasAdp13Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_B_13", adp13Col);
            hasAdp23Col = findLoopTagIndex(adpLoop, "_atom_site_aniso_B_23", adp23Col);
        }

        if (!(hasAdp11Col && hasAdp22Col && hasAdp33Col && hasAdp12Col && hasAdp13Col && hasAdp23Col))
            return false;

        //-------

        adp.resize(nAnisoAtoms,vector<double>(6));
        adp_sigma.resize(nAnisoAtoms, vector<double>(6));
        adp_precision.resize(nAnisoAtoms, vector<double>(6));
        for (anisoAtomIdx = 0; anisoAtomIdx < nAnisoAtoms; anisoAtomIdx++)
        {

            discamb::string_utilities::string_to_number_and_uncertainty(adpLoop.values[adp11Col][anisoAtomIdx], adp[anisoAtomIdx][0], adp_sigma[anisoAtomIdx][0]);
            discamb::string_utilities::string_to_number_and_uncertainty(adpLoop.values[adp22Col][anisoAtomIdx], adp[anisoAtomIdx][1], adp_sigma[anisoAtomIdx][1]);
            discamb::string_utilities::string_to_number_and_uncertainty(adpLoop.values[adp33Col][anisoAtomIdx], adp[anisoAtomIdx][2], adp_sigma[anisoAtomIdx][2]);
            discamb::string_utilities::string_to_number_and_uncertainty(adpLoop.values[adp12Col][anisoAtomIdx], adp[anisoAtomIdx][3], adp_sigma[anisoAtomIdx][3]);
            discamb::string_utilities::string_to_number_and_uncertainty(adpLoop.values[adp13Col][anisoAtomIdx], adp[anisoAtomIdx][4], adp_sigma[anisoAtomIdx][4]);
            discamb::string_utilities::string_to_number_and_uncertainty(adpLoop.values[adp23Col][anisoAtomIdx], adp[anisoAtomIdx][5], adp_sigma[anisoAtomIdx][5]);

            adp_precision[anisoAtomIdx][0] = get_precision(adpLoop.values[adp11Col][anisoAtomIdx]);
            adp_precision[anisoAtomIdx][1] = get_precision(adpLoop.values[adp22Col][anisoAtomIdx]);
            adp_precision[anisoAtomIdx][2] = get_precision(adpLoop.values[adp33Col][anisoAtomIdx]);
            adp_precision[anisoAtomIdx][3] = get_precision(adpLoop.values[adp12Col][anisoAtomIdx]);
            adp_precision[anisoAtomIdx][4] = get_precision(adpLoop.values[adp13Col][anisoAtomIdx]);
            adp_precision[anisoAtomIdx][5] = get_precision(adpLoop.values[adp23Col][anisoAtomIdx]);
                                     
            if (hasBani)
                for (int i = 0; i < 6; i++)
                {
                    adp[anisoAtomIdx][i] /= 8 * M_PI * M_PI;
                    adp_sigma[anisoAtomIdx][i] /= 8 * M_PI * M_PI;
                    adp_precision[anisoAtomIdx][i] /= 8 * M_PI * M_PI;
                }
        }

        //-------

        vector<string>::const_iterator it;
        
        atomLoopToAnisoLoopEntryMap.resize(nAtoms, -1);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            it = find(adpLoop.values[anisoLabelCol].begin(), adpLoop.values[anisoLabelCol].end(), atomLoop.values[atomLabelCol][atomIdx]);
            if (it != adpLoop.values[anisoLabelCol].end())
                atomLoopToAnisoLoopEntryMap[atomIdx] = static_cast<int>(distance(adpLoop.values[anisoLabelCol].begin(), it));
        }
        return true;
    }
}


namespace discamb {



    namespace cif_io {

        string Loop::to_string(
            int nItemsPerLine)
            const
        {
            if (values.empty())
                return string();

            string s = "loop_\n";
            for (const string& header: loop_headers)
                s += header + "\n";
            int nEntries = values[0].size();
            for (int entryIdx = 0; entryIdx < nEntries; entryIdx++)
            {
                bool isMultiline;
                for (int headerIdx = 0; headerIdx < values.size(); headerIdx++)
                {
                    stringstream ss;
                    ss << values[headerIdx][entryIdx];

                    string firstLine;
                    getline(ss, firstLine);
                    isMultiline = (firstLine != values[headerIdx][entryIdx]);

                    if (isMultiline)
                        s += "\n;\n" + values[headerIdx][entryIdx] + "\n;\n";
                    else
                    {
                        if (headerIdx > 0)
                            s += " ";
                        if (values[headerIdx][entryIdx].find(' ') != string::npos)
                            s += "'" + values[headerIdx][entryIdx] + "'";
                        else
                            s += values[headerIdx][entryIdx];
                    }

                    if (nItemsPerLine > 0)
                        if ((headerIdx + 1) % nItemsPerLine == 0)
                            if (headerIdx + 1 != values.size())
                                s += "\n";

                }
                if (!isMultiline)
                    s += "\n";
            }
            return s;
        }

        void readCif(
            const std::string &fName, 
            std::vector<DataSet> &data) 
        {
            std::string input_string;
            std::ifstream cifFile(fName.c_str(), std::ifstream::in);
            
            if (!cifFile.is_open())
                on_error::throwException(string("can not open CIF file '") + fName + string("'"), __FILE__, __LINE__);

            std::string tmp;
            while (getline(cifFile, tmp)) {
                input_string += tmp;
                input_string += "\n";
            }
            cifFile.close();


            //

            CifBuilder builder;
            //ucif::parser parsed(&builder, input_string, fName, /*strict=*/true);
            ucif::parser parsed(&builder, input_string, fName, false);
            data = builder.dataSets;

/*
            for (int dataIdx = 0; dataIdx < builder.dataSets.size(); dataIdx++)
            {
                DataSet &data = builder.dataSets[dataIdx];

                std::cout << data.dataSetName << std::endl;


                for (auto &item : data.key_value_items)
                    std::cout << item.first << "   " << item.second << std::endl;


                for (int loopIdx = 0; loopIdx < data.loops.size(); loopIdx++)
                {
                    for (int itemIdx = 0; itemIdx < data.loops[loopIdx].loop_headers.size(); itemIdx++)
                        std::cout << data.loops[loopIdx].loop_headers[itemIdx] << std::endl;

                    for (int valueIdx = 0; valueIdx < data.loops[loopIdx].values[0].size(); valueIdx++)
                    {
                        for (int itemIdx = 0; itemIdx < data.loops[loopIdx].loop_headers.size(); itemIdx++)
                            std::cout << data.loops[loopIdx].values[itemIdx][valueIdx] << "  ";
                        std::cout << std::endl;
                    }


                }
            }
            */
        }

        bool findLoopIndex(
            const DataSet &data,
            const std::string &tag,
            int &index,
            int &headerIdx)
        {
            for (int i = 0; i < data.loops.size(); i++)
                for (int j = 0; j < data.loops[i].loop_headers.size();j++)
                    if (tag == data.loops[i].loop_headers[j])
                    {
                        index = i;
                        headerIdx = j;
                        return true;
                    }
            return false;
        }

        bool findLoopTagIndex(
            const Loop &loop,
            const std::string &tag,
            int &headerIdx)
        {
            vector<string>::const_iterator it = find(loop.loop_headers.begin(), loop.loop_headers.end(), tag);
            if (it != loop.loop_headers.end())
            {
                headerIdx = distance(loop.loop_headers.begin(), it);
                return true;
            }
            return false;
        }
        

        bool findLoopIndex(
            const DataSet &data,
            const std::vector<std::string> &tags,
            int &loopIdx,
            std::vector<int> &tagsIndices)
        {
            vector<int> loops(tags.size());
            int previousLoopIndex;

            tagsIndices.resize(tags.size());

            for (int i=0; i<tags.size(); i++)
            {
                if (!findLoopIndex(data, tags[i], loopIdx, tagsIndices[i]))
                    return false;
                if (i !=0 )
                    if (previousLoopIndex != loopIdx)
                        return false;
                
                previousLoopIndex = loopIdx;
            }
            return true;
        }


        /*  
        uani in separate loop 
        has loop with _atom_site_type_symbol or _atom_site_label or both (not separated)
        if uani loop then has all _atom_site_aniso_B_ij or all _atom_site_aniso_U_ij otherwise do process info on B_ij/U_ij
        */

        void extractAtomSiteData(
            const DataSet &data,
            std::vector<AtomInCrystal> &atoms)
        {
            int loopIdx, tagIdx, nAtoms, atomIdx;
            int uIsoOrEqCol, adpTypeCol, xFracCol, yFracCol, zFracCol, atomLabelCol, siteTypeCol, occupancyCol, muliplicityCol;
            bool hasIsoOrEq, hasAdpType, hasX_Frac, hasY_Frac, hasZ_Frac, hasAtomLabel, hasSiteType, hasOccupancy, hasMuliplicity;
            bool hasAnisoData;
            int adpLoopIdx, atomAnisoLabelCol;
			double xFrac, yFrac, zFrac, xFracSigma, yFracSigma, zFracSigma, occupancy, occupancySigma;


            atoms.clear();

            hasSiteType = findLoopIndex(data, "_atom_site_type_symbol", loopIdx, siteTypeCol);

            if (!hasSiteType)
                hasAtomLabel = findLoopIndex(data, "_atom_site_label", loopIdx, atomLabelCol);
            else
                hasAtomLabel = findLoopTagIndex(data.loops[loopIdx], "_atom_site_label", atomLabelCol);

            if ( !hasAtomLabel && !hasSiteType)
                return;
            
            const Loop &atomSiteLoop = data.loops[loopIdx];

            // check ADP info location
            vector<int> atomLoopToAnisoLoopEntryMap;
            vector<vector<double> > adp, adpSigma, adpPrecision;


            hasAnisoData = getAnisoData(data, atomLoopToAnisoLoopEntryMap, adp, adpSigma, adpPrecision);
            
            //-----

            hasIsoOrEq = findLoopTagIndex(atomSiteLoop, "_atom_site_U_iso_or_equiv", uIsoOrEqCol);
            hasX_Frac = findLoopTagIndex(atomSiteLoop, "_atom_site_fract_x", xFracCol);
            hasY_Frac = findLoopTagIndex(atomSiteLoop, "_atom_site_fract_y", yFracCol);
            hasZ_Frac = findLoopTagIndex(atomSiteLoop, "_atom_site_fract_z", zFracCol);
            hasOccupancy = findLoopTagIndex(atomSiteLoop, "_atom_site_occupancy", occupancyCol);
            hasMuliplicity = findLoopTagIndex(atomSiteLoop, "_atom_site_symmetry_multiplicity", muliplicityCol);
            hasAdpType = findLoopTagIndex(atomSiteLoop, "_atom_site_adp_type", adpTypeCol);
            hasIsoOrEq = findLoopTagIndex(atomSiteLoop, "_atom_site_U_iso_or_equiv", uIsoOrEqCol);


            if (hasAtomLabel)
                nAtoms = atomSiteLoop.values[atomLabelCol].size();
            else
                nAtoms = atomSiteLoop.values[siteTypeCol].size();

            atoms.resize(nAtoms);

            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                if (hasAtomLabel)
                    atoms[atomIdx].label = atomSiteLoop.values[atomLabelCol][atomIdx];

                if (hasSiteType)
                    atoms[atomIdx].type = atomSiteLoop.values[siteTypeCol][atomIdx];
                else
                    if (hasAtomLabel)
                        atoms[atomIdx].type = periodic_table::symbol(basic_chemistry_utilities::atomicNumberFromLabel(atoms[atomIdx].label));

                    

                if (hasAdpType)
                {
                    if (atomSiteLoop.values[adpTypeCol][atomIdx] == string("Uiso"))
                        if (hasIsoOrEq)
                        {
                            atoms[atomIdx].adp.resize(1);
                            atoms[atomIdx].adp_sigma.resize(1);
                            atoms[atomIdx].adp_precision.resize(1);
                            discamb::string_utilities::string_to_number_and_uncertainty(
                                atomSiteLoop.values[uIsoOrEqCol][atomIdx],
                                atoms[atomIdx].adp[0],
                                atoms[atomIdx].adp_sigma[0]);

                            atoms[atomIdx].adp_precision[0] = get_precision(atomSiteLoop.values[uIsoOrEqCol][atomIdx]);
                        }

                    if (atomSiteLoop.values[adpTypeCol][atomIdx] == string("Uani"))
                        if (hasAnisoData)
							if (atomLoopToAnisoLoopEntryMap[atomIdx] >= 0)
							{
								atoms[atomIdx].adp = adp[atomLoopToAnisoLoopEntryMap[atomIdx]];
								atoms[atomIdx].adp_sigma = adpSigma[atomLoopToAnisoLoopEntryMap[atomIdx]];
                                atoms[atomIdx].adp_precision = adpPrecision[atomLoopToAnisoLoopEntryMap[atomIdx]];
							}
                }
                else 
                    if (hasAnisoData)
						if (atomLoopToAnisoLoopEntryMap[atomIdx] >= 0)
						{
							atoms[atomIdx].adp = adp[atomLoopToAnisoLoopEntryMap[atomIdx]];
							atoms[atomIdx].adp_sigma = adpSigma[atomLoopToAnisoLoopEntryMap[atomIdx]];
						}

				if (hasX_Frac && hasY_Frac && hasZ_Frac)
				{
					string_utilities::string_to_number_and_uncertainty(atomSiteLoop.values[xFracCol][atomIdx], xFrac, xFracSigma);
					string_utilities::string_to_number_and_uncertainty(atomSiteLoop.values[yFracCol][atomIdx], yFrac, yFracSigma);
					string_utilities::string_to_number_and_uncertainty(atomSiteLoop.values[zFracCol][atomIdx], zFrac, zFracSigma);
					atoms[atomIdx].coordinates = Vector3d(xFrac, yFrac, zFrac);
					atoms[atomIdx].coordinates_sigma = Vector3d(xFracSigma, yFracSigma, zFracSigma);
                    atoms[atomIdx].coordinates_precision = Vector3d(
                        get_precision(atomSiteLoop.values[xFracCol][atomIdx]),
                        get_precision(atomSiteLoop.values[yFracCol][atomIdx]),
                        get_precision(atomSiteLoop.values[zFracCol][atomIdx]));

					//atoms[atomIdx].coordinates = Vector3d(atof(atomSiteLoop.values[xFracCol][atomIdx].c_str()),
					//	atof(atomSiteLoop.values[yFracCol][atomIdx].c_str()),
					//	atof(atomSiteLoop.values[zFracCol][atomIdx].c_str()));
				}

				if (hasOccupancy)
				{
					string_utilities::string_to_number_and_uncertainty(atomSiteLoop.values[occupancyCol][atomIdx], 
						                                               atoms[atomIdx].occupancy, atoms[atomIdx].occupancy_sigma);
                    atoms[atomIdx].occupancy_precision = get_precision(atomSiteLoop.values[occupancyCol][atomIdx]);
					//atoms[atomIdx].occupancy = atof(atomSiteLoop.values[occupancyCol][atomIdx].c_str());
				}
				else
				{
					atoms[atomIdx].occupancy = 1.0;
					atoms[atomIdx].occupancy_sigma = 0.0;
				}

                if (hasMuliplicity)
                    atoms[atomIdx].multiplicity = atoi(atomSiteLoop.values[muliplicityCol][atomIdx].c_str());
                else
                    atoms[atomIdx].multiplicity = 0;

            }

        }

        void extractSpaceGroupData(
            const DataSet &data,
            SpaceGroup &spaceGroup)
        {
            bool hasSymmOps;
            int symmOpsColumn, symmLoopIdx, symmOpIdx, nSymmOps;
            spaceGroup = SpaceGroup();

            hasSymmOps = findLoopIndex(data, "_space_group_symop_operation_xyz", symmLoopIdx, symmOpsColumn);
            if(!hasSymmOps)
                hasSymmOps = findLoopIndex(data, "_symmetry_equiv_pos_as_xyz", symmLoopIdx, symmOpsColumn);

            if (!hasSymmOps)
                return;

            nSymmOps = data.loops[symmLoopIdx].values[symmOpsColumn].size();
            vector<SpaceGroupOperation> symmetryOperations;
            for (symmOpIdx = 0; symmOpIdx < nSymmOps; symmOpIdx++)
                symmetryOperations.push_back(SpaceGroupOperation(data.loops[symmLoopIdx].values[symmOpsColumn][symmOpIdx]));
                
            spaceGroup.set(symmetryOperations);
        }

        void extractUnitCellData(
            const DataSet &data,
            UnitCell &unitCell)
        {
            unitCell = UnitCell();

            map<string, string>::const_iterator it;
            double a, b, c, alpha, beta, gamma;
            
            it = data.key_value_items.find("_cell_length_a");
            if (it != data.key_value_items.end())
                a = atof(it->second.c_str());
            else
                return;

            it = data.key_value_items.find("_cell_length_b");
            if (it != data.key_value_items.end())
                b = atof(it->second.c_str());
            else
                return;

            it = data.key_value_items.find("_cell_length_c");
            if (it != data.key_value_items.end())
                c = atof(it->second.c_str());
            else
                return;

            it = data.key_value_items.find("_cell_angle_alpha");
            if (it != data.key_value_items.end())
                alpha = atof(it->second.c_str());
            else
                return;

            it = data.key_value_items.find("_cell_angle_beta");
            if (it != data.key_value_items.end())
                beta = atof(it->second.c_str());
            else
                return;

            it = data.key_value_items.find("_cell_angle_gamma");
            if (it != data.key_value_items.end())
                gamma = atof(it->second.c_str());
            else
                return;

            unitCell.set(a, b, c, alpha, beta, gamma);
        }

        /*
        required:
        loop with _refln_index_h, _refln_index_k, _refln_index_l (in the same loop)
        */

        void extractReflectionData(
            const DataSet &data,
            ReflectionsData &reflectionsData,
            std::map<string, double> &reflns_scale_meas_intensity)
        {
            reflectionsData.clear();
            reflns_scale_meas_intensity.clear();

            int loopIdx, loopIdx2, loopIdx3, hIdx, kIdx, lIdx;
            bool hasH, hasK, hasL;
            int columnIdx, intensity_calculated_idx, intensity_measured_idx, intensity_sigma_idx, scale_group_code_idx, refln_observed_status_idx;
            int scale_meas_intensity_idx;
            vector<int> h, k, l;

            // check if there are data

            hasH = findLoopIndex(data, "_refln_index_h", loopIdx, hIdx);
            hasK = findLoopIndex(data, "_refln_index_k", loopIdx2, kIdx);
            hasL = findLoopIndex(data, "_refln_index_l", loopIdx3, lIdx);
                

            if (! (hasH && hasK && hasL)) // no data
                return;

            if ( (loopIdx != loopIdx2) || (loopIdx != loopIdx3) ) // refln indices in different loops
                return;

            // 
            //int intensity_measured_idx, intensity_observed_idx, intensity_sigma_idx, scale_group_code_idx, refln_observed_status_idx;

            const Loop &loop = data.loops[loopIdx];

            convertVectorOfStrings(loop.values[hIdx], h);
            convertVectorOfStrings(loop.values[kIdx], k);
            convertVectorOfStrings(loop.values[lIdx], l);

            for (int i = 0; i < h.size(); i++)
                reflectionsData.hkl.push_back(Vector3i(h[i], k[i], l[i]));

            if (findLoopTagIndex(loop, "_refln_intensity_calc", columnIdx))
                convertVectorOfStrings(loop.values[columnIdx], reflectionsData.intensity_calculated);

            if (findLoopTagIndex(loop, "_refln_intensity_meas", columnIdx))
                convertVectorOfStrings(loop.values[columnIdx], reflectionsData.intensity_measured);

            if (findLoopTagIndex(loop, "_refln_intensity_sigma", columnIdx))
                convertVectorOfStrings(loop.values[columnIdx], reflectionsData.intensity_sigma);

            if (findLoopTagIndex(loop, "_refln_scale_group_code", columnIdx))
                convertVectorOfStrings(loop.values[columnIdx], reflectionsData.scale_group_code);

            if (findLoopTagIndex(loop, "_refln_observed_status", columnIdx))
                convertVectorOfStrings(loop.values[columnIdx], reflectionsData.observed_status);

            //reflns_scale_meas_intensity
            /*
             _reflns_scale_group_code
             _reflns_scale_meas_intensity

             scale_meas_intensity_idx
             scale_group_code_idx
            */

            if (findLoopIndex(data, "_reflns_scale_meas_intensity", loopIdx, scale_meas_intensity_idx))
            {
                if (findLoopTagIndex(data.loops[loopIdx], "_reflns_scale_group_code", scale_group_code_idx))
                {
                    int i, n = data.loops[loopIdx].values[0].size();
                    for (i = 0; i < n; i++)
                        reflns_scale_meas_intensity[data.loops[loopIdx].values[scale_group_code_idx][i]] =
                            atof(data.loops[loopIdx].values[scale_meas_intensity_idx][i].c_str());
                }
            }
            


        }


        void cifDataToCrystal(
            const DataSet &data, 
            Crystal &crystal) 
        {
            crystal.atoms.clear();
            crystal.spaceGroup = SpaceGroup();
            crystal.unitCell = UnitCell();

            
            extractSpaceGroupData(data, crystal.spaceGroup);
            extractUnitCellData(data, crystal.unitCell);
            extractAtomSiteData(data, crystal.atoms);

            crystal.adpConvention = structural_parameters_convention::AdpConvention::U_cif;
            crystal.xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;

            vector<vector<SpaceGroupOperation> > atomPointGroups;
            for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
            {
                crystal_structure_utilities::findAtomSymmetry(crystal, atomIdx, atomPointGroups, 0.05);
                crystal.atoms[atomIdx].multiplicity = crystal.spaceGroup.nSymmetryOperations() / atomPointGroups[0].size();
            }


            // _atom_site_type_symbol.
            // Examples: ‘C’, ‘Cu2+’, ‘H(SDS)’, ‘dummy’, ‘FeNi’

        }

        void saveCif(
            const std::string& fName,
            const Crystal& _crystal,
            const HC_ModelParameters& parameters,
            const std::vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs)
        {
            bool hasMultipoleModel = !parameters.atom_to_type_map.empty();
            
            //for(parameters.atom_to_type_map

            Crystal crystal = _crystal;
            DataSet multipoleData;
            std::vector<Vector3d> dummyAtomFractionalPosition;

            if (hasMultipoleModel)
            {
                saveMultipoleModel(crystal, multipoleData,
                    dummyAtomFractionalPosition,
                    parameters, lcs);
                for (int i = 0; i < dummyAtomFractionalPosition.size(); i++)
                {
                    crystal.atoms.resize(crystal.atoms.size() + 1);
                    crystal.atoms.back().adp = { 0.01 };
                    crystal.atoms.back().coordinates = dummyAtomFractionalPosition[i];
                    crystal.atoms.back().label = "DUM" + to_string(i + 1);
                    crystal.atoms.back().type = ".";
                    crystal.atoms.back().occupancy = 0.0;
                    //crystal_structure_utilities::findAtomSymmetry()
                }
            }
                    


            ofstream out(fName);
            out << setprecision(4) << fixed;
            out << "data_discamb\n\n"
                << "_cell_angle_alpha  " << setw(8) << crystal.unitCell.alpha() << "\n"
                << "_cell_angle_beta   " << setw(8) << crystal.unitCell.beta() << "\n"
                << "_cell_angle_gamma  " << setw(8) << crystal.unitCell.gamma() << "\n"
                << "_cell_length_a     " << setw(8) << crystal.unitCell.a() << "\n"
                << "_cell_length_b     " << setw(8) << crystal.unitCell.b() << "\n"
                << "_cell_length_c     " << setw(8) << crystal.unitCell.c() << "\n"
                << "\nloop_\n"
                << "_space_group_symop_id\n"
                << "_space_group_symop_operation_xyz\n";

            string symmOp;

            for (int i = 0; i < crystal.spaceGroup.nSymmetryOperations(); i++)
            {
                crystal.spaceGroup.getSpaceGroupOperation(i).get(symmOp);
                out << setw(2) << i+1 << " " << string_utilities::toLower(symmOp) << "\n";
            }

            out << "\nloop_\n"
                << "_atom_site_label\n"
                << "_atom_site_type_symbol\n"
                << "_atom_site_fract_x\n"
                << "_atom_site_fract_y\n"
                << "_atom_site_fract_z\n"
                << "_atom_site_U_iso_or_equiv\n"
                << "_atom_site_adp_type\n"
                << "_atom_site_occupancy\n"
                << "_atom_site_symmetry_multiplicity\n";
            double u_iso;
            out << setprecision(5);
            for (int i = 0; i < crystal.atoms.size(); i++)
            {
                out << setw(8) << crystal.atoms[i].label
                    << setw(4) << crystal.atoms[i].type
                    << setw(9) << crystal.atoms[i].coordinates[0]
                    << setw(9) << crystal.atoms[i].coordinates[1]
                    << setw(9) << crystal.atoms[i].coordinates[2];
                if (crystal.atoms[i].adp.empty())
                    u_iso = 0.0;
                else
                    u_iso = crystal.atoms[i].adp[0];
                out << setw(9) << u_iso;
                if (crystal.atoms[i].adp.size() == 6)
                    out << " Uani ";
                else
                    out << " Uiso ";
                out << setw(9) << crystal.atoms[i].occupancy << setw(6) << crystal.atoms[i].multiplicity << "\n";
            }

            out << "\nloop_\n"
                << "_atom_site_aniso_label\n"
                << "_atom_site_aniso_U_11\n"
                << "_atom_site_aniso_U_22\n"
                << "_atom_site_aniso_U_33\n"
                << "_atom_site_aniso_U_12\n"
                << "_atom_site_aniso_U_13\n"
                << "_atom_site_aniso_U_23\n";

            for (int i = 0; i < crystal.atoms.size(); i++)
                if (crystal.atoms[i].adp.size() == 6)
                {
                    out << setw(9) << crystal.atoms[i].label << " ";
                    for (int j = 0; j < 6; j++)
                        out << setw(9) << crystal.atoms[i].adp[j];
                    out << endl;
                }

            if (hasMultipoleModel)
                for (auto& loop : multipoleData.loops)
                    out << "\n" << loop.to_string(10);

        }

        void saveCif(
            const std::string &fName, 
            const Crystal &crystal)
        {
            ofstream out(fName);
            out << setprecision(4) << fixed;
            out << "data_discamb\n"
                << "_cell_angle_alpha  " << setw(8) << crystal.unitCell.alpha() << "\n"
                << "_cell_angle_beta   " << setw(8) << crystal.unitCell.beta() << "\n"
                << "_cell_angle_gamma  " << setw(8) << crystal.unitCell.gamma() << "\n"
                << "_cell_length_a     " << setw(8) << crystal.unitCell.a() << "\n"
                << "_cell_length_b     " << setw(8) << crystal.unitCell.b() << "\n"
                << "_cell_length_c     " << setw(8) << crystal.unitCell.c() << "\n"
                << "loop_\n"
                << "_symmetry_equiv_pos_as_xyz\n";
            
            string symmOp;

            for (int i = 0; i < crystal.spaceGroup.nSymmetryOperations(); i++)
            {
                crystal.spaceGroup.getSpaceGroupOperation(i).get(symmOp);
                out << symmOp << "\n";
            }
            
            out << "loop_\n"
                << "_atom_site_type_symbol\n"
                << "_atom_site_label\n"
                << "_atom_site_fract_x\n"
                << "_atom_site_fract_y\n"
                << "_atom_site_fract_z\n"
                << "_atom_site_U_iso_or_equiv\n"
                << "_atom_site_adp_type\n"
                << "_atom_site_occupancy\n";
            double u_iso;
            out << setprecision(5);
            for (int i = 0; i < crystal.atoms.size(); i++)
            {
                out << setw(4) << crystal.atoms[i].type
                    << setw(8) << crystal.atoms[i].label
                    << setw(9) << crystal.atoms[i].coordinates[0]
                    << setw(9) << crystal.atoms[i].coordinates[1]
                    << setw(9) << crystal.atoms[i].coordinates[2];
                if (crystal.atoms[i].adp.empty())
                    u_iso = 0.0;
                else
                    u_iso = crystal.atoms[i].adp[0];
                out << setw(9) << u_iso;
                if (crystal.atoms[i].adp.size() == 6)
                    out << " Uani ";
                else
                    out << " Uiso ";
                out << setw(9) << crystal.atoms[i].occupancy << "\n";
            }
            
            out << "loop_\n"
                << "_atom_site_aniso_label\n"
                << "_atom_site_aniso_U_11\n"
                << "_atom_site_aniso_U_22\n"
                << "_atom_site_aniso_U_33\n"
                << "_atom_site_aniso_U_12\n"
                << "_atom_site_aniso_U_13\n"
                << "_atom_site_aniso_U_23\n";

            for (int i = 0; i < crystal.atoms.size(); i++)
                if (crystal.atoms[i].adp.size() == 6)
                {
                    out << setw(9) << crystal.atoms[i].label << " ";
                    for (int j = 0; j < 6; j++)
                        out << setw(9) << crystal.atoms[i].adp[j];
                    out << endl;
                }



        }

        template<>
        void convertVectorOfStrings<char>(const std::vector<std::string> stringData, std::vector<char> &convertedData)
        {
            int i, n = stringData.size();
            convertedData.resize(n);
            for (i = 0; i < n; i++)
                if (stringData[i].empty())
                    convertedData[i] = (char)0;
                else
                    convertedData[i] = stringData[i][0];
        }

        void readFcf(
            const std::string &fName,
            SpaceGroup &spaceGroup,
            UnitCell &unitCell,
            std::vector<Vector3i> &hkl,
            std::vector<double> &_refln_F_meas,
            std::vector<double> &_refln_F_sigma,
            std::vector<double> &_refln_A_calc,
            std::vector<double> &_refln_B_calc)
        {
            vector<DataSet> data;
            readCif(fName, data);
            if (data.empty())
                on_error::throwException(string("no data_ in fcf file '") + fName + string("'"), __FILE__, __LINE__);
            extractSpaceGroupData(data[0], spaceGroup);
            extractUnitCellData(data[0], unitCell);
            vector<string> tags{ "_refln_index_h", "_refln_index_k", "_refln_index_l",
                                  "_refln_F_meas", "_refln_F_sigma", "_refln_A_calc", "_refln_B_calc" };
                
            int loopIdx;
            vector<int> tagsIndices;

            if (!cif_io::findLoopIndex(data[0], tags, loopIdx, tagsIndices))
                on_error::throwException("problem when reading fcf file, cannot find loop with tags:\
                                         _refln_index_h, _refln_index_k, _refln_index_l, _refln_F_meas,\
                                         _refln_F_sigma, _refln_A_calc, _refln_B_calc", __FILE__, __LINE__);
            cif_io::Loop &loop = data[0].loops[loopIdx];
            int nHkl = loop.values[0].size();
            hkl.resize(nHkl);
            for (int i = 0; i < nHkl; i++)
                hkl[i] = Vector3i(stoi(loop.values[tagsIndices[0]][i]),
                                  stoi(loop.values[tagsIndices[1]][i]),
                                  stoi(loop.values[tagsIndices[2]][i]));

            convertVectorOfStrings(loop.values[tagsIndices[3]], _refln_F_meas);
            convertVectorOfStrings(loop.values[tagsIndices[4]], _refln_F_sigma);
            convertVectorOfStrings(loop.values[tagsIndices[5]], _refln_A_calc);
            convertVectorOfStrings(loop.values[tagsIndices[6]], _refln_B_calc);
        }

        void writeFcf(
            const std::string& fName,
            const SpaceGroup& sg,
            const UnitCell& uc,
            const std::vector<Vector3i>& hkl,
            const std::vector<double>& _refln_F_squared_meas,
            const std::vector<double>& _refln_F_squared_sigma,
            const std::vector<double>& _refln_F_squared_calc,
            bool scale,
            double a,
            double b)
        {
            if (scale)
                on_error::not_implemented(__FILE__,__LINE__);

            ofstream out(fName);
            out << "data_test\n\n";

            out << "loop_\n"
                << " _space_group_symop_operation_xyz\n";
            
            vector<string> symmOpsStr;
            int i, n = sg.nSymmetryOperations();

            for (i = 0; i < n; i++)            
                out << " " << sg.getSpaceGroupOperation(i).string() << "\n";
            
            out << "\n";

            out << "_cell_length_a    " << uc.a() << "\n"
                << "_cell_length_b    " << uc.b() << "\n"
                << "_cell_length_c    " << uc.c() << "\n"
                << "_cell_angle_alpha " << uc.alpha() << "\n"
                << "_cell_angle_beta  " << uc.beta() << "\n"
                << "_cell_angle_gamma " << uc.gamma() << "\n";
            
            out << "\n";

            out << "loop_\n"
                << " _refln_index_h\n"
                << " _refln_index_k\n"
                << " _refln_index_l\n"
                << " _refln_F_squared_calc\n"
                << " _refln_F_squared_meas\n"
                << " _refln_F_squared_sigma\n";
            
            n = hkl.size();
            out << setprecision(4) << fixed;
            for (int i = 0; i < n; i++)
                out << setw(4) << hkl[i][0]
                    << setw(4) << hkl[i][1]
                    << setw(4) << hkl[i][2]
                    << setw(16) << _refln_F_squared_calc[i]
                    << setw(16) << _refln_F_squared_meas[i]
                    << setw(16) << _refln_F_squared_sigma[i] << "\n";
                
            out << "\n";
            out.close();
        }

        void extractMultipoleModel(
            const Crystal& crystal,
            const DataSet& data,
            HC_ModelParameters& parameters,
            vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs)
        {
            parameters = HC_ModelParameters();
            lcs.clear();

            // look for dummy atoms
            vector<string> dummyAtomLabels;
            vector<int> dummyAtomIndices;
            for (int i=0;i<crystal.atoms.size();i++)
            {
                const auto& atom = crystal.atoms[i];
                string lowerCaseAtomLabel = string_utilities::toLower(atom.label);
                if (lowerCaseAtomLabel.find("dummy") != string::npos)
                {
                    ;
                }
            }

            // 

            int loopIdx;
            int headerIdx;
            if (!findLoopIndex(data, "_atom_rho_multipole_atom_label", loopIdx, headerIdx))
                return;
            
            

            const Loop &loop = data.loops[loopIdx];

            // check if atoms are in the same order
            // does not support the case when they are not
            vector<string> atomLabels;
            for (const auto& atom : crystal.atoms)
            {
                string lowerCaseAtomLabel = string_utilities::toLower(atom.label);
                //if(lowerCaseAtomLabel.find("dummy" != )
                //atomLabels.push_back(atom.label);
            }
                

            /*
 _atom_rho_multipole_coeff_Pc
 _atom_rho_multipole_coeff_Pv
 _atom_rho_multipole_coeff_P00
 _atom_rho_multipole_coeff_P10
 _atom_rho_multipole_coeff_P11
 _atom_rho_multipole_coeff_P1-1
 _atom_rho_multipole_coeff_P20
 _atom_rho_multipole_kappa
 _atom_rho_multipole_kappa_prime0
 _atom_rho_multipole_kappa_prime1
 _atom_rho_multipole_kappa_prime2
 _atom_rho_multipole_kappa_prime3
 _atom_rho_multipole_kappa_prime4
 _atom_rho_multipole_radial_slater_n0
 _atom_rho_multipole_radial_slater_zeta0
 _atom_rho_multipole_radial_slater_n1
 _atom_rho_multipole_radial_slater_zeta1
 _atom_rho_multipole_radial_slater_n2
 _atom_rho_multipole_radial_slater_zeta2
 _atom_rho_multipole_radial_slater_n3
 _atom_rho_multipole_radial_slater_zeta3
 _atom_rho_multipole_radial_slater_n4
 _atom_rho_multipole_radial_slater_zeta4
 _atom_rho_multipole_configuration
 _atom_rho_multipole_scat_core
 _atom_rho_multipole_scat_valence

            */
            
        }


        void saveMultipoleModel(
            const Crystal& crystal,
            DataSet& data,
            std::vector<Vector3d>& dummyAtomFractionalPosition,
            const HC_ModelParameters& parameters,
            const vector < LocalCoordinateSystem<AtomInCrystalID> >& lcs)
        {
            data = DataSet();
            data.loops.push_back(Loop());
            dummyAtomFractionalPosition.clear();

            auto& loop = data.loops[0];
            int atomIdx, nAtoms = parameters.atom_to_type_map.size();

            loop.loop_headers.push_back("_atom_rho_multipole_atom_label");
            loop.loop_headers.push_back("_atom_rho_multipole_coeff_Pv");

            map<pair<int, int>, int> lm2idx;
            int counter = 2;
            for(int l=0; l<=4; l++)
                for (int m = -l; m <= l; m++)
                {
                    loop.loop_headers.push_back("_atom_rho_multipole_coeff_P" + to_string(l) + to_string(m));
                    lm2idx[{l, m}] = counter++;
                }
            loop.loop_headers.push_back("_atom_rho_multipole_kappa");
            for (int l = 0; l <= 4; l++)
                loop.loop_headers.push_back("_atom_rho_multipole_kappa_prime" + to_string(l));

            for (int l = 0; l <= 4; l++)
            {
                loop.loop_headers.push_back("_atom_rho_multipole_radial_slater_n" + to_string(l));
                loop.loop_headers.push_back("_atom_rho_multipole_radial_slater_zeta" + to_string(l));
            }

                
            loop.values.resize(loop.loop_headers.size());

            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                int typeIdx = parameters.atom_to_type_map[atomIdx];
                int wfnIdx = parameters.atom_to_wfn_map[atomIdx];
                loop.values[0].push_back(crystal.atoms[atomIdx].label);
                loop.values[1].push_back(string_utilities::realToString(parameters.type_parameters[typeIdx].p_val, 4, true));

                for (int i = 2; i < 27; i++)
                    loop.values[i].push_back("0.0000");

                int maxL = parameters.type_parameters[typeIdx].p_lm.size();
                maxL--;
                for (int l = 0; l <= maxL; l++)
                    for (int m = -l; m <= l; m++)
                        loop.values[lm2idx[{l, m}]][atomIdx] = 
                            string_utilities::realToString(parameters.type_parameters[typeIdx].p_lm[l][l + m], 4, true);
                loop.values[27].push_back(string_utilities::realToString(parameters.type_parameters[typeIdx].kappa_spherical_valence, 4, true));

                for (int l = 0; l <= 4; l++)
                    if (l <= maxL)
                        loop.values[28 + l].push_back(string_utilities::realToString(parameters.type_parameters[typeIdx].kappa_deformation_valence, 4, true));
                    else
                        loop.values[28 + l].push_back("1.0000");

                for (int l = 0; l <= 4; l++)
                {
                    double value = l; 
                    if (l <= maxL)
                        value = parameters.wfn_parameters[wfnIdx].deformation_valence_power[l];
                    loop.values[33 + 2*l].push_back(string_utilities::realToString(value, 4, true));
                    value = 1.0;
                    if (l <= maxL)
                        value = parameters.wfn_parameters[wfnIdx].deformation_valence_exponent;
                    loop.values[33 + 2 * l + 1].push_back(string_utilities::realToString(value, 4, true));

                }

            }

            // lcs
            data.loops.push_back(Loop());
            auto &lcsLoop = data.loops[1];

            lcsLoop.loop_headers = {
                "_atom_local_axes_atom_label",
                "_atom_local_axes_atom0",
                "_atom_local_axes_atom1",
                "_atom_local_axes_atom2",
                "_atom_local_axes_ax1",
                "_atom_local_axes_ax2" };
            lcsLoop.values.resize(6);

            vector<string> xyzStr{ "X","Y","Z" };
            for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                string atom1Label, atom2Label;

                LocalCoordinateSystemCalculator lcsCalc(lcs[atomIdx], crystal);
                Matrix3d lcsMatrix;
                lcsCalc.calculate(lcsMatrix, crystal);
                Vector3d xyz[3];
                lcsCalc.calculate(xyz[0], xyz[1], xyz[2], crystal);
                Vector3d centralAtomFrac, centralAtomCart, dummyCart, dummyFrac;

                bool atomInAsymmetricUnitAsRefPoint = false;
                if (lcs[atomIdx].direction1_type == LcsDirectionType::AVERAGE_POSITION && lcs[atomIdx].refPoint_1.size() == 1)
                    if (lcs[atomIdx].refPoint_1[0].getSymmetryOperation().isIdentity())
                        atomInAsymmetricUnitAsRefPoint = true;

                if (atomInAsymmetricUnitAsRefPoint)
                    atom1Label = crystal.atoms[lcs[atomIdx].refPoint_1[0].index()].label;
                else // use dummy
                {
                    centralAtomFrac = crystal.atoms[atomIdx].coordinates;
                    crystal.unitCell.fractionalToCartesian(centralAtomFrac, centralAtomCart);
                    dummyCart = centralAtomCart + xyz[lcs[atomIdx].coordinate_1];
                    crystal.unitCell.cartesianToFractional(dummyCart, dummyFrac);
                        
                    dummyAtomFractionalPosition.push_back(dummyFrac);
                    atom1Label = "DUM" + to_string(dummyAtomFractionalPosition.size());
                }
                
                atomInAsymmetricUnitAsRefPoint = false;
                if (lcs[atomIdx].direction2_type == LcsDirectionType::AVERAGE_POSITION && lcs[atomIdx].refPoint_2.size() == 1)
                    if (lcs[atomIdx].refPoint_2[0].getSymmetryOperation().isIdentity())
                        atomInAsymmetricUnitAsRefPoint = true;

                if (atomInAsymmetricUnitAsRefPoint)
                    atom2Label = crystal.atoms[lcs[atomIdx].refPoint_2[0].index()].label;
                else // use dummy
                {
                    centralAtomFrac = crystal.atoms[atomIdx].coordinates;
                    crystal.unitCell.fractionalToCartesian(centralAtomFrac, centralAtomCart);
                    dummyCart = centralAtomCart + xyz[lcs[atomIdx].coordinate_2];
                    crystal.unitCell.cartesianToFractional(dummyCart, dummyFrac);

                    dummyAtomFractionalPosition.push_back(dummyFrac);
                    atom2Label = "DUM" + to_string(dummyAtomFractionalPosition.size());
                }

                lcsLoop.values[0].push_back(crystal.atoms[atomIdx].label);
                lcsLoop.values[1].push_back(atom1Label);
                lcsLoop.values[2].push_back(crystal.atoms[atomIdx].label);
                lcsLoop.values[3].push_back(atom2Label);
                lcsLoop.values[4].push_back(xyzStr[lcs[atomIdx].coordinate_1]);
                lcsLoop.values[5].push_back(xyzStr[lcs[atomIdx].coordinate_2]);
            }
        }


    } // namespace cif_io


}