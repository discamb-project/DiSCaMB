#include "discamb/IO/olex2_io.h"
#include "npy.hpp"

#include "discamb/BasicUtilities/StringUtilities.h"
#include "discamb/BasicUtilities/OnError.h"

#include <fstream>
#include <sstream>

using namespace std;

namespace {

    void rawVcovToCrystalVcov(
        vector<vector<double> > const& data,
        vector<string> const& labels,
        discamb::CrystalVarianceCovarianceMatrix& vcov)
    {
        vector<string> words;

        vector< pair<int, discamb::CrystalVarianceCovarianceMatrix::variable_type > > variableList;
        int atomIdx = -1;
        string atomLabel;
        int i, n = data.size();

        //int atom
        string varLabel, currentVarLabel;
        bool newVariable;
        map<string, discamb::CrystalVarianceCovarianceMatrix::variable_type> str2varType{ {"xyz", discamb::CrystalVarianceCovarianceMatrix::variable_type::xyz},
            {"u_iso", discamb::CrystalVarianceCovarianceMatrix::variable_type::u_iso}, {"u_aniso", discamb::CrystalVarianceCovarianceMatrix::variable_type::u_aniso},
            {"occ", discamb::CrystalVarianceCovarianceMatrix::variable_type::occupancy} };


        for (i = 0; i < n; i++)
        {
            newVariable = false;

            discamb::string_utilities::split(labels[i], words, '.');

            if (words[1].size() == 1)
                currentVarLabel = string("xyz");
            else
            {
                if (words[1].substr(0, 2) == string("ui"))
                    currentVarLabel = string("u_iso");
                else
                {
                    if (words[1][0] == 'u')
                        currentVarLabel = string("u_aniso");
                    else
                        currentVarLabel = string("occ");
                }
            }

            if (currentVarLabel != varLabel)
            {
                newVariable = true;
                varLabel = currentVarLabel;
            }

            if (atomLabel != words[0])
            {
                atomLabel = words[0];
                newVariable = true;
                atomIdx++;
            }

            if (newVariable)
                variableList.push_back({ atomIdx, str2varType[varLabel] });
        }

        vcov.set(data, variableList);
    }


}



namespace discamb
{
    namespace olex2_io {

        void read_vcov(
            const std::string& fName,
            std::vector<std::vector<double> >& vcov,
            std::vector<std::string>& dataLabels)
        {
            ifstream in(fName);

            if (!in.good())
                on_error::throwException(string("cannot read olex2 variance-covariance matrix file  '") + fName + string("'"), __FILE__, __LINE__);

            string line;
            vector<string> words;
            getline(in, line);
            getline(in, line);
            string_utilities::split(line, dataLabels);
            getline(in, line);
            string_utilities::split(line, words);
            int i, j, k, nData = dataLabels.size();

            vcov.assign(nData, vector<double>(nData, 0.0));
            k = 0;
            for(i=0; i<nData; i++)
                for (j = i; j <nData; j++)
                {
                    vcov[i][j] = stod(words[k++]);
                    vcov[j][i] = vcov[i][j];
                }
            in.close();
        }

        void read_vcov(
            const std::string& fName,
            CrystalVarianceCovarianceMatrix& vcov)
        {
            vector<vector<double> > data;
            //vector<string> labels, words;
            vector<string> labels;
            olex2_io::read_vcov(fName, data, labels);

            rawVcovToCrystalVcov(data, labels, vcov);

            //vector< pair<int, CrystalVarianceCovarianceMatrix::variable_type > > variableList;
            //int atomIdx = -1;
            //string atomLabel;
            //int i, n = data.size();

            ////int atom
            //string varLabel, currentVarLabel;
            //bool newVariable;
            //map<string, CrystalVarianceCovarianceMatrix::variable_type> str2varType{ {"xyz", CrystalVarianceCovarianceMatrix::xyz},
            //    {"u_iso", CrystalVarianceCovarianceMatrix::u_iso}, {"u_aniso", CrystalVarianceCovarianceMatrix::u_aniso},
            //    {"occ", CrystalVarianceCovarianceMatrix::occupancy} };


            //for (i = 0; i < n; i++)
            //{
            //    newVariable = false;

            //    string_utilities::split(labels[i], words, '.');

            //    if (words[1].size() == 1)
            //        currentVarLabel = string("xyz");
            //    else
            //    {
            //        if (words[1].substr(0, 2) == string("ui"))
            //            currentVarLabel = string("u_iso");
            //        else
            //        {
            //            if (words[1][0] == 'u')
            //                currentVarLabel = string("u_aniso");
            //            else
            //                currentVarLabel = string("occ");
            //        }
            //    }

            //    if (currentVarLabel != varLabel)
            //    {
            //        newVariable = true;
            //        varLabel = currentVarLabel;
            //    }

            //    if (atomLabel != words[0])
            //    {
            //        atomLabel = words[0];
            //        newVariable = true;
            //        atomIdx++;
            //    }

            //    if (newVariable)
            //        variableList.push_back({ atomIdx, str2varType[varLabel] });
            //}

            //vcov.set(data, variableList);
        }

        void read_vcov_npy(
            const std::string& fName,
            CrystalVarianceCovarianceMatrix& vcov)
        {
            vector<vector<double> > data;
            vector<string> labels;
            olex2_io::read_vcov_npy(fName, data, labels);

            rawVcovToCrystalVcov(data, labels, vcov);
        }

        void read_vcov_npy(
            const std::string& fName,
            std::vector<std::vector<double> >& vcov,
            std::vector<std::string>& dataLabels)
        {
            ifstream in(fName, ifstream::in | ifstream::binary);
            ofstream out("discamb_tmp_yuturymauysuiauy5_dusiau.npy", ofstream::out | ofstream::binary);
            stringstream sstream(stringstream::binary);
            string labelsLine;

            if (!in.good())
                on_error::throwException(string("cannot read file '") + fName + string("'"), __FILE__, __LINE__);
            if (!out.good())
                on_error::throwException(string("cannot write to file out.npy"), __FILE__, __LINE__);
            char c;
            int nBreakLines = 0;
            while (in.get(c))
            {
                if (nBreakLines == 1 && c != '\n')
                    labelsLine += c;
                if (nBreakLines > 1)
                    out.put(c);
                if (c == '\n')
                    nBreakLines++;
            }
            in.close();
            out.close();
            string_utilities::split(labelsLine, dataLabels);
            //out.close();
            vector<unsigned long> shape;
            bool fortran_order;
            vector<double> data1d;
            npy::LoadArrayFromNumpy("discamb_tmp_yuturymauysuiauy5_dusiau.npy", shape, fortran_order, data1d);

            int i, j, k, nData = dataLabels.size();

            vcov.assign(nData, vector<double>(nData, 0.0));
            k = 0;
            for (i = 0; i < nData; i++)
                for (j = i; j < nData; j++)
                {
                    vcov[i][j] = data1d[k++];
                    vcov[j][i] = vcov[i][j];
                }
            
        } 
    }
}
