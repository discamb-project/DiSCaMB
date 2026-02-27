#include "discamb/IO/wfn_io.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/QuantumChemistry/ecp_electron_density_tables.h"


#include <fstream>
#include <sstream>
#include <set>


using namespace std;
/*
struct WfnFileData {
std::string label;
bool isGaussian;
std::vector<std::string> center_label;
std::vector<Vector3d> center_position;
std::vector<double> center_charge;
std::vector<int> primitive_to_center;
std::vector<int> primitive_type;
std::vector<double> primitive_exponents;
std::vector<double> molecular_orbital_energy;
std::vector<double> molecular_orbital_occupancy;
std::vector<std::vector<double> > molecular_orbitals;
double energy;
double virial;
};

*/

namespace {

    bool isWfx(const string& fileName)
    {

        /*
                string last4 = fileName.substr(fileName.size() - 4);
                
                if(!isWfx)
                    if(last4 != string(".wfn"))
                        on_error::throwException("cannot deduce wavefunction file type from file extension: "+ last4, __FILE__, __LINE__);

        */
        if (fileName.size() >= 4)
        {
            string last4 = fileName.substr(fileName.size() - 4);
            bool isWfx = (last4 == string(".wfx"));
            if (!isWfx)
            {
                if (last4 == string(".wfn"))
                    return false;
            }
            else
                return true;
        }


        ifstream in(fileName);
        if (!in.good())
            discamb::on_error::throwException("cannot read file '" + fileName + "'", __FILE__, __LINE__);
        string line;
        int lineIdx = 0;
        while (getline(in, line))
        {
            if (lineIdx++ > 100)
            {
                in.close();
                return false;
            }
            if (line.find("<") != string::npos && line.find(">") != string::npos)
            {
                in.close();
                return true;
            }
        }
        
        in.close();
        return false;
    }

    double fortranStod(string& s)
    {
        auto e(s.find_first_of("Dd"));
        if (e != std::string::npos)
            s[e] = 'E';
        return stod(s);
    }

    // 5.27044619 10.54089238-14.11639406
    discamb::Vector3d getPosition(string str)
    {
        string s;
        for (auto c : str)
        {
            if (c == '-')
                s += ' ';
            s += c;
        }
        vector<string> words;
        discamb::string_utilities::split(s, words);
        if (words.size() != 3)
            discamb::on_error::throwException(
                string("error when processing wfn, cannot process line contatining:'") + str + string("'"),
                __FILE__, __LINE__);
        return discamb::Vector3d(stod(words[0]), stod(words[1]), stod(words[2]) );
    }

    bool keywordLine(const string& line, bool& start, string& name)
    {
        if (line.empty())
            return false;
        if (line[0] != '<')
            return false;

        start = true;
        name = line.substr(1, line.find('>') - 1);
        if (name[0] == '/')
        {
            start = false;
            name = name.substr(1);
        }
        return true;
    }

    template<typename T>
    void linesToItems(
        const vector<string>& lines,
        vector<T>& v)
    {
        v.clear();
        vector<string> words;
        for (auto const& line : lines)
        {
            discamb::string_utilities::split(line, words);
            for(auto word: words)
                v.push_back(discamb::string_utilities::convertFromString<T>(word));
        }
    }

    void getRequiredNodeLines(
        const string& key,
        map<string, vector<string> > data,
        vector<string>& lines)
    {
        auto it = data.find(key);
        if (it == data.end())
            discamb::on_error::throwException(string("cannot find entry '") + key + string("' in wfx file"), __FILE__, __LINE__);
        lines.swap(it->second);
    }

    void processWfxOrbitalLines(
        const vector<string> &lines,
        vector<vector<double> > &coefficients)
    {
        coefficients.clear();
        vector<double> coeff;
        vector<string> moLines;
        int i, n = lines.size();
        for (i=0; i<n; i++)
        {
            if (lines[i].find("<MO Number>") != string::npos)
            {
                i += 2;
                if (!moLines.empty())
                {
                    linesToItems(moLines, coeff);
                    coefficients.push_back(coeff);
                }
                moLines.clear();
            }
            else
                moLines.push_back(lines[i]);
        }

        if (!moLines.empty())
        {
            linesToItems(moLines, coeff);
            coefficients.push_back(coeff);
        }

    }
    /*
    read wfx line which have fields in brackets in separate lines, 
    e.g. not like this
    1.819000000000e-001  1.819000000000e-001 </Primitive Exponents>
    but like this:
    1.819000000000e-001  1.819000000000e-001 
    </Primitive Exponents>
    and fields in brackets are not nested, ie not like e.g.:
    <Additional Electron Density Function (EDF)>
    <Number of EDF Primitives>
    <> is followed by </> not by another <>
    */

    void processFormatedWfxLines(
        const vector<string> &lines,
        std::map<std::string, std::vector<std::string> >& data)
    {
        vector<string> nodeLines;
        string line, currentNodeName, nodeName;
        bool inNode, start;
        


        inNode = false;
        int lineIdx = 1;
        //while (getline(input, line))
        for (string line : lines)
        {

            if (keywordLine(line, start, nodeName))
            {
                if (inNode)
                {
                    if (start)
                        nodeLines.push_back(line);
                    else
                    {
                        if (currentNodeName == nodeName) //end of node
                        {
                            data[currentNodeName] = nodeLines;
                            currentNodeName.clear();
                            nodeLines.clear();
                            inNode = false;
                        }
                        else
                        {
                            nodeLines.push_back(line);
                        }
                    }
                }
                else
                {
                    if (start)
                    {
                        currentNodeName = nodeName;
                        inNode = true;
                    }
                    else
                        discamb::on_error::throwException("error when processing wfx stram", __FILE__, __LINE__);
                }
            }
            else
                if (inNode)
                    nodeLines.push_back(line);

            lineIdx++;
        }

    }

    // EDFs lines are removed from lines
    void extractEDFs(
        vector<string>& lines,
        vector < map<string, vector<string> > >& edfsData)
    {
        edfsData.clear();

        bool inEdf = false;
        vector<string> noEdfLines, edfLines;

        for (auto& line : lines)
        {
            if (line.find("<Additional Electron Density Function (EDF)>") != string::npos)
            {
                inEdf = true;
                edfLines.clear();
                edfsData.resize(edfsData.size());
                continue;
            }

            if (line.find("</Additional Electron Density Function (EDF)>") != string::npos)
            {
                inEdf = false;
                edfsData.resize(edfsData.size() + 1);
                processFormatedWfxLines(edfLines, edfsData.back());
                continue;
            }

            if (inEdf)
                edfLines.push_back(line);
            else
                noEdfLines.push_back(line);
        }
        lines = noEdfLines;
    }

/*------------------------------------------------------------------
  red wfx file lines and splits orca wfx lines of this type:
 
  1.819000000000e-001  1.819000000000e-001 </Primitive Exponents>

  in such a way that brackets <> are in sepparate lines, like this:

  1.819000000000e-001  1.819000000000e-001
  </Primitive Exponents>

*/

    void readAndFormatWfxLines(
        std::istream& input,
        vector<string>& lines)
    {
        string line;
        while (getline(input, line))
        {
            auto startPosition = line.find('<');
            if (startPosition != string::npos)
            {
                string subLine;
                auto endPosition = line.find('>');
                if (startPosition != 0)
                {
                    subLine = line.substr(0, startPosition);
                    discamb::string_utilities::trim(subLine);
                    if (!subLine.empty())
                        lines.push_back(subLine);
                }

                lines.push_back(line.substr(startPosition, endPosition - startPosition + 1));

                if (endPosition != line.size() - 1)
                {
                    subLine = line.substr(endPosition + 1);
                    discamb::string_utilities::trim(subLine);
                    if (!subLine.empty())
                        lines.push_back(subLine);
                }

            }
            else
                lines.push_back(line);
        }

    }

    void processEdfNode(
        const std::map<std::string, std::vector<std::string> > &nodes,
        discamb::wfn_io::AdditionalElectronDensity &edf)
    {
        vector<string> nodeLines;
        // EDF Primitive Centers
        getRequiredNodeLines(string("EDF Primitive Centers"), nodes, nodeLines);
        linesToItems(nodeLines, edf.primitive_to_center);

        getRequiredNodeLines(string("EDF Primitive Types"), nodes, nodeLines);
        linesToItems(nodeLines, edf.primitive_type);

        getRequiredNodeLines(string("EDF Primitive Exponents"), nodes, nodeLines);
        linesToItems(nodeLines, edf.primitive_exponents);

        getRequiredNodeLines(string("EDF Primitive Coefficients"), nodes, nodeLines);
        linesToItems(nodeLines, edf.primitive_coefficients);
    }
}

namespace discamb {
    namespace wfn_io {

        void read_wfx(
            const std::string& fileName,
            std::map<std::string, std::vector<std::string> >& data,
            std::vector< std::map<std::string, std::vector<std::string> > >& edfs)
        {
            ifstream in(fileName);
            read_wfx(in, data, edfs);
            in.close();
        }

        void read_wfx(
            std::istream& input,
            std::map<std::string, std::vector<std::string> >& data,
            std::vector< std::map<std::string, std::vector<std::string> > >& edfs)
        {
            data.clear();
            vector<string> lines;

            if (!input.good())
                on_error::throwException("cannot read wfx file", __FILE__, __LINE__);

            readAndFormatWfxLines(input, lines);
            extractEDFs(lines, edfs);
            processFormatedWfxLines(lines, data);
        }

        void read_wfx(
            const std::string& fileName,
            WfnFileData& data)
        {
            ifstream in(fileName);
            if (!in.good())
                on_error::throwException("cannot read wfx file '" + fileName + "'", __FILE__, __LINE__);

            read_wfx(in, data);

            in.close();
        }

        void read_wfx(
            std::istream& in,
            WfnFileData& data)
        {
            data = WfnFileData();

            map<string, vector<string> > nodes;
            std::vector< std::map<std::string, std::vector<std::string> > > edfNodes;
            read_wfx(in, nodes, edfNodes);

            for (auto& node : edfNodes)
            {
                data.edfs.resize(data.edfs.size() + 1);
                processEdfNode(node, data.edfs.back());
            }

            //string s, s2, xStr, yStr, zStr, line;
            int nPrimitives, nMo, nNuclei, nucleiIdx;
            vector<string> words, nodeLines;

            
            getRequiredNodeLines(string("Title"), nodes, nodeLines);
            data.label = nodeLines[0];

            getRequiredNodeLines(string("Keywords"), nodes, nodeLines);
            string_utilities::split(nodeLines[0], words);
            data.isGaussian = (find(words.begin(), words.end(), string("GTO")) != words.end());

            getRequiredNodeLines(string("Number of Nuclei"), nodes, nodeLines);
            nNuclei = stoi(nodeLines[0]);
            getRequiredNodeLines(string("Number of Primitives"), nodes, nodeLines);
            nPrimitives = stoi(nodeLines[0]);
            getRequiredNodeLines(string("Number of Occupied Molecular Orbitals"), nodes, nodeLines);
            nMo = stoi(nodeLines[0]);

            getRequiredNodeLines(string("Nuclear Names"), nodes, nodeLines);
            linesToItems(nodeLines, data.center_label);

            getRequiredNodeLines(string("Atomic Numbers"), nodes, nodeLines);
            linesToItems(nodeLines, data.atomic_numbers);

            getRequiredNodeLines(string("Nuclear Charges"), nodes, nodeLines);
            linesToItems(nodeLines, data.center_charge);


            getRequiredNodeLines(string("Nuclear Cartesian Coordinates"), nodes, nodeLines);
            vector<double> xyz;
            linesToItems(nodeLines, xyz);
            data.center_position.resize(nNuclei);
            for (nucleiIdx = 0; nucleiIdx < nNuclei; nucleiIdx++)
                data.center_position[nucleiIdx] = Vector3d(xyz[3 * nucleiIdx], xyz[3 * nucleiIdx + 1], xyz[3 * nucleiIdx + 2]);

            getRequiredNodeLines(string("Primitive Centers"), nodes, nodeLines);
            linesToItems(nodeLines, data.primitive_to_center);

            getRequiredNodeLines(string("Primitive Types"), nodes, nodeLines);
            linesToItems(nodeLines, data.primitive_type);

                                        
            getRequiredNodeLines(string("Primitive Exponents"), nodes, nodeLines);
            linesToItems(nodeLines, data.primitive_exponents);


            getRequiredNodeLines(string("Molecular Orbital Occupation Numbers"), nodes, nodeLines);
            linesToItems(nodeLines, data.molecular_orbital_occupancy);

            getRequiredNodeLines(string("Molecular Orbital Energies"), nodes, nodeLines);
            linesToItems(nodeLines, data.molecular_orbital_energy);

            getRequiredNodeLines(string("Molecular Orbital Primitive Coefficients"), nodes, nodeLines);
            processWfxOrbitalLines(nodeLines, data.molecular_orbitals);

        }

        void read_wavefunction(
            const std::string& fileName,
            WfnFileData& data)
        {
            try 
            {
                bool isWfxFile = isWfx(fileName);
                ifstream in(fileName);
                
                if (!in.good())
                    on_error::throwException("cannot read wave function file '" + fileName + "'", __FILE__, __LINE__);

                if (isWfxFile)
                    read_wfx(in, data);
                else
                    read_wfn(in, data);
                in.close();
            }
            catch (exception& e)
            {
                on_error::throwException(string("prome when trying to read file '") + fileName +
                    string("' ") + e.what(), __FILE__, __LINE__);
            }

        }

        void add_edf_from_library(
            WfnFileData& data)
        {
            if (data.atomic_numbers.size() != data.center_charge.size())
                return;

            if (!data.edfs.empty())
                return;

            AdditionalElectronDensity edf;
            for (int atomIdx = 0; atomIdx < data.atomic_numbers.size(); atomIdx++)
            {
                if (fabs(data.atomic_numbers[atomIdx] - data.center_charge[atomIdx]) > 0.001)
                {
                    int nCore = int(fabs(data.atomic_numbers[atomIdx] - data.center_charge[atomIdx] + 0.01));
                    vector<double> exponents, coefficients;

                    ecp_electron_density_tables::ecp_electron_density(
                        data.atomic_numbers[atomIdx],
                        nCore, 
                        exponents, 
                        coefficients);

                    edf.primitive_coefficients.insert(
                        edf.primitive_coefficients.end(),
                        coefficients.begin(),
                        coefficients.end());

                    edf.primitive_exponents.insert(
                        edf.primitive_exponents.end(),
                        exponents.begin(),
                        exponents.end());

                    int nPrimitives = exponents.size();

                    for (int primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                    {
                        edf.primitive_to_center.push_back(atomIdx + 1);
                        edf.primitive_type.push_back(1);
                    }
                }
            }
            data.edfs.push_back(edf);
        }

        void read_wfn(
            const std::string &fileName,
            WfnFileData &data,
            bool skipNonAtomicCenters)
        {
            try {
                ifstream in(fileName);
                read_wfn(in, data, skipNonAtomicCenters);

                in.close();
            }
            catch (...)
            {
                on_error::throwException(string("prome when trying to read file '") + fileName +
                    string("'"), __FILE__, __LINE__);
            }
        }

        void read_wfn(
            std::istream &in,
            WfnFileData &data,
            bool skipNonAtomicCenters)
        {
            string s, s2, xStr, yStr, zStr, line;
            int nPrimitives, nMo, moIdx, primitiveIdx, nNuclei, nucleiIdx;
            vector<string> words;


            if (!in.good())
            {
                on_error::throwException("cannot read wfn stream", __FILE__, __LINE__);
                return;
            }

            getline(in, data.label);

            in >> s;
            data.isGaussian = (s == string("GAUSSIAN"));

            // we assume that it is always GAUSSIAN
            if (!data.isGaussian)
                on_error::throwException("incorrect format of wfn stream", __FILE__, __LINE__);

            //     5    MOL ORBITALS     30   PRIMITIVES  3    NUCLEI
            in >> nMo >> s >> s >> nPrimitives >> s >> nNuclei >> s;
            /*data.center_charge.resize(nNuclei);
            data.center_label.resize(nNuclei);
            data.center_position.resize(nNuclei);
            data.atomic_numbers.resize(nNuclei);
            */
            //int nCentersToCount = 0;
            string centerLabel;
            
            for (nucleiIdx = 0; nucleiIdx < nNuclei; nucleiIdx++)
            {
                in >> centerLabel;
                if (centerLabel[0] != 'Q' || !skipNonAtomicCenters)
                {
                    data.center_label.push_back(centerLabel);
                    getline(in, line);
                    string_utilities::split(line, words);
                    data.center_charge.push_back(stod(words.back()));
                    string_utilities::split(line, words,')');
                    //  Q 1120    (CENTRE1120)   5.27044619 10.54089238-14.11639406  CHARGE =  3.3
                    string xyzStr = words[1].substr(0, words[1].find('C'));

                    data.center_position.push_back(getPosition(xyzStr));

                   // in >> data.center_label[nucleiIdx] >> s >> s;
                   // if (s.back() != ')')
                   //     in >> s;
                   // in >> xStr >> yStr >> zStr >> s >> s >> s;
                   // data.center_position[nucleiIdx].set(fortranStod(xStr), fortranStod(yStr), fortranStod(zStr));
                   // data.atomic_numbers[nucleiIdx] = periodic_table::atomicNumber(data.center_label[nucleiIdx]);
                    data.atomic_numbers.push_back(periodic_table::atomicNumber(centerLabel));
                }
                else
                    getline(in, line);
            }

            data.primitive_to_center.clear();
            data.primitive_exponents.clear();
            data.primitive_type.clear();

            in >> s;
            while (in.good() && s != string("MO"))
            {

                if (string("CENT") == s.substr(0,4))
                {
                    in >> s;
                    getline(in, s);
                    string_utilities::split(s, words);
                    for (auto word : words)
                        data.primitive_to_center.push_back(stoi(word));
                }
                if (string("TYPE") == s)
                {
                    //TYPE ASSIGNMENTS      1  1  1  1  1  1  1  1  1  2  2  2  3  3  3  4  4  4  1  2
                    in >> s;
                    getline(in, s);
                    string_utilities::split(s, words);
                    for (auto word : words)
                        data.primitive_type.push_back(stoi(word));

                }
                if (string("EXPONENTS") == s)
                {
                    getline(in, s);
                    string_utilities::split(s, words);
                    for (auto word : words)
                        data.primitive_exponents.push_back(fortranStod(word));

                }
                in >> s;
            }

            if (!in.good())
                on_error::throwException("error when processing data in wfn format", __FILE__, __LINE__);

            data.molecular_orbitals.assign(nMo, vector<double>(nPrimitives, 0.0));
            data.molecular_orbital_energy.resize(nMo);
            data.molecular_orbital_occupancy.resize(nMo);

            string occStr, energyStr, orbitalCoeffStr, virialStr;
            moIdx = 0;
            while (s.substr(0,2) == string("MO") && in.good())
            {
                //MO   1   OCC  NO    =   2.00000000 ORB. ENERGY = -26.09206896
                in >> s;
                //moIdx = stoi(s) - 1;
                std::getline(in, line);
                
                occStr = line.substr(line.find("OCC NO =") + 8, 10);
                energyStr = line.substr(line.find("ORB. ENERGY =") + 13, 10);

                data.molecular_orbital_occupancy[moIdx] = fortranStod(occStr);
                data.molecular_orbital_energy[moIdx] = fortranStod(energyStr);

                for (primitiveIdx = 0; primitiveIdx < nPrimitives; primitiveIdx++)
                {
                    in >> orbitalCoeffStr;
                    data.molecular_orbitals[moIdx][primitiveIdx] = fortranStod(orbitalCoeffStr);
                }

                in >> s;
                moIdx++;
            }
// THE  HF ENERGY =   -378.454977754666 THE VIRIAL(-V/T)=   2.00644967
// THE SCF ENERGY =   -225.239160350067 THE VIRIAL(-V/T)=   2.00541908
//RHF      ENERGY =     -224.0027785184   VIRIAL(-V/T)  =   2.00098495
        }

        void read_atomic_wfn_database(
            const std::string &fileName, 
            std::vector<WfnFileData> &data) 
        {
            stringstream wfnStringStream;
            string wfnString, line;
            data.clear();
            ifstream in(fileName);
            vector<string> words;
            WfnFileData wfnEntry;
            
            while (in.good())
            {
                getline(in, line);
                string_utilities::split(line, words);
                if (words.size() == 0)
                {
                    wfnStringStream.str(wfnString);
                    wfnStringStream.clear();
                    wfnString = "";
                    read_wfn(wfnStringStream, wfnEntry);
                    data.push_back(wfnEntry);
                }
                else
                    wfnString += line + string("\n");
            }
            in.close();
        }

        void mergeEdfs(
            const std::vector<AdditionalElectronDensity>& edfs,
            AdditionalElectronDensity &edf)
        {
            edf = AdditionalElectronDensity();

            for (auto const& x : edfs)
            {
                edf.primitive_coefficients.insert(
                    edf.primitive_coefficients.end(), 
                    x.primitive_coefficients.begin(),
                    x.primitive_coefficients.end());

                edf.primitive_exponents.insert(
                    edf.primitive_exponents.end(),
                    x.primitive_exponents.begin(),
                    x.primitive_exponents.end());

                edf.primitive_to_center.insert(
                    edf.primitive_to_center.end(),
                    x.primitive_to_center.begin(),
                    x.primitive_to_center.end());

                edf.primitive_type.insert(
                    edf.primitive_type.end(),
                    x.primitive_type.begin(),
                    x.primitive_type.end());
            }
        }

        //void removeCentersWithNoBasisFunctions(
        //    WfnFileData& data)
        //{
        //    std::set<int> centersWithFunctions;
        //    centersWithFunctions.insert(data.primitive_to_center.begin(), data.primitive_to_center.end());
        //    vector<int> centersToRemove;
        //    int centerIdx, nCenters = data.atomic_numbers.size();
        //    for (centerIdx = 0; centerIdx < nCenters; centerIdx++)
        //        if (centersWithFunctions.find(centerIdx + 1) == centersWithFunctions.end())
        //            centersToRemove.push_back(centerIdx);

        //    // in numeration starting from 0
        //    vector<int> newIdices(data.atomic_numbers.size(),0);

        //    for (auto c : centersWithFunctions)
        //    {
        //        int centerIdx = c - 1;
        //        int shift = 0;
        //        for (auto center2remove : centersToRemove)
        //            if (centerIdx > center2remove)
        //                shift++;
        //        newIdices[centerIdx] = centerIdx - shift;
        //    }

        //    nCenters = centersWithFunctions.size();
        //    vector<int> atomic_numbers(nCenters);
        //    vector<double> charge(nCenters);
        //    vector<string> label(nCenters);
        //    vector<Vector3d> position(nCenters);
        //    
        //    for (auto idx : centersWithFunctions)
        //    {
        //        int oldIdx = idx - 1;
        //        int newIdx = newIdices[oldIdx];
        //        atomic_numbers[newIdx] = data.atomic_numbers[oldIdx];
        //        charge[newIdx] = data.center_charge[oldIdx];
        //        label[newIdx] = data.center_label[oldIdx];
        //        position[newIdx] = data.center_position[oldIdx];
        //    }


        //    

        //    //for(auto &x: data.primitive_to_center)
        //    //for()
        //    //data.atomic_numbers
        //    //data.center_charge
        //    //data.center_label
        //    //data.center_position
        //    //data.primitive_to_center
        //    
        //}
    }
}

/*
GAUSSIAN              5 MOL ORBITALS     30 PRIMITIVES        3 NUCLEI
O    1    (CENTRE  1)  -0.18157451  0.03321874  0.12439407  CHARGE =  8.0
H    2    (CENTRE  2)   1.11644029 -1.20653392  0.35861804  CHARGE =  1.0
H    3    (CENTRE  3)   0.33615580  0.94078402 -1.35377058  CHARGE =  1.0
CENTRE ASSIGNMENTS    1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
CENTRE ASSIGNMENTS    1  1  2  2  2  2  3  3  3  3
TYPE ASSIGNMENTS      1  1  1  1  1  1  1  1  1  2  2  2  3  3  3  4  4  4  1  2
TYPE ASSIGNMENTS      3  4  1  1  1  1  1  1  1  1
EXPONENTS  5.4846717E+03 8.2523495E+02 1.8804696E+02 5.2964500E+01 1.6897570E+01
EXPONENTS  5.7996353E+00 1.5539616E+01 3.5999336E+00 1.0137618E+00 1.5539616E+01
EXPONENTS  3.5999336E+00 1.0137618E+00 1.5539616E+01 3.5999336E+00 1.0137618E+00
EXPONENTS  1.5539616E+01 3.5999336E+00 1.0137618E+00 2.7000580E-01 2.7000580E-01
EXPONENTS  2.7000580E-01 2.7000580E-01 1.8731137E+01 2.8253937E+00 6.4012170E-01
EXPONENTS  1.6127780E-01 1.8731137E+01 2.8253937E+00 6.4012170E-01 1.6127780E-01
MO  1                     OCC NO =   2.00000000 ORB. ENERGY = -26.09206896
8.29792992E-01  1.52723284E+00  2.47136461E+00  3.24867638E+00  2.78637184E+00
9.52707785E-01 -6.63409449E-03 -2.96012215E-03  8.74127832E-03 -1.51275850E-03
-1.16549825E-03 -5.11716610E-04 -5.49670031E-03 -4.23490898E-03 -1.85935352E-03
*/