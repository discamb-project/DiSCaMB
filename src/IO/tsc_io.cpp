#include "discamb/IO/tsc_io.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/BasicUtilities/StringUtilities.h"

#include <fstream>
//#include <iostream>


using namespace std;

namespace discamb {

    namespace {
        int maxNChar(const std::vector<Vector3i>& hkl)
        {
            int nChar = 0;
            for (auto h : hkl)
            {
                int n = max({ to_string(h[0]).size(), to_string(h[1]).size(), to_string(h[2]).size() });
                if (n > nChar)
                    nChar = n;
            }
            cout << "nChar " << nChar << endl;
            return nChar;
        }

        string formFactorAsString(const std::complex<double>& ff)
        {
            stringstream ss;
            ss << " " << setprecision(6) << fixed << ff.real() << ","
                << setprecision(6) << fixed << ff.imag();
            string s;
            ss >> s;
            return s;

        }


        int maxNChar(const std::vector<vector<complex<double> > > & hkl)
        {
            int nChar = 0;

            for (auto h_vec : hkl)
                for (auto h: h_vec)
                {
                    int n = formFactorAsString(h).size();


                    if (n > nChar)
                        nChar = n;
                }
            cout << "nChar " << nChar << endl;
            return nChar;
        }


    }

    namespace tsc_io {
        void read_tsc(
            const std::string& fileName,
            std::vector<std::string>& atomLabels,
            std::vector<Vector3i>& hklSet,
            std::vector<std::vector<std::complex<double> > >& atomicFormFactors)
        {
            atomLabels.clear();
            hklSet.clear();
            atomicFormFactors.clear();

            ifstream in(fileName);

            if (!in.good())
                on_error::throwException("can not read olex tsc file", __FILE__, __LINE__);

            bool dataSection = false;
            string line;
            
            /*
            TITLE: form factors for urea generated with discamb
AD: FALSE
SYMM: expanded
SCATTERERS: C O N H1 H2
DATA:
            */

            
            vector<string> words, ff_words;
            while (in.good() && !dataSection)
            {
                getline(in, line);
                string_utilities::split(line, words);
                if (words[0] == string("DATA:"))
                    dataSection = true;
                if (words[0] == string("SCATTERERS:"))
                    atomLabels.insert(atomLabels.end(), words.begin() + 1, words.end());
            }
            
            int nAtoms = atomLabels.size();
            Vector3i hkl;
            vector<complex<double> > ff(nAtoms);

            if (dataSection)
            {
                //getline(in, line);
                while (in.good())
                {
                    getline(in, line);
                    string_utilities::split(line, words);
                    if (words.size() == nAtoms + 3)
                    {
                        hkl = { stoi(words[0]), stoi(words[1]), stoi(words[2]) };
                        //cout << hkl[0] << " " << hkl[1] << " " << hkl[2] << "\n";
                        hklSet.push_back(hkl);
                        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                        {
                            string_utilities::split(words[atomIdx + 3], ff_words, ',');
                            ff[atomIdx] = { stod(ff_words[0]), stod(ff_words[1]) };
                        }
                        atomicFormFactors.push_back(ff);
                    }
                }
            }

            in.close();
        }

        void read_tsc2(
            const std::string& fileName,
            std::vector<std::string>& atomLabels,
            std::vector<Vector3i>& hklSet,
            //[hklIdx][atomIdx]
            std::vector<std::vector<std::complex<double> > >& atomicFormFactors)
        {
            atomLabels.clear();
            hklSet.clear();
            atomicFormFactors.clear();

            ifstream in(fileName);

            if (!in.good())
                on_error::throwException("can not read olex tsc file", __FILE__, __LINE__);

            bool dataSection = false;
            string line;

            /*
            TITLE: form factors for urea generated with discamb
AD: FALSE
SYMM: expanded
SCATTERERS: C O N H1 H2
DATA:
            */


            vector<string> words, ff_words;
            while (in.good() && !dataSection)
            {
                getline(in, line);
                string_utilities::split(line, words);
                if (words[0] == string("DATA:"))
                    dataSection = true;
                if (words[0] == string("SCATTERERS:"))
                    atomLabels.insert(atomLabels.end(), words.begin() + 1, words.end());
            }

            int atomIdx, nAtoms = atomLabels.size();
            Vector3i hkl;
            vector<complex<double> > ff(nAtoms);
            string s;
            int k, l;

            bool run = true;
            if (dataSection)
                do
                {
                    in >> s;
                    if (in.good())
                    {
                        if (s.size() != 0)
                        {
                            if (!isdigit(s[0]) && s[0]!='-')
                                run = false;
                        }
                        else
                            run = false;
                    }
                    else
                        run = false;
                    if (run)
                    {
                        in >> k >> l;
                        hklSet.push_back(Vector3i(stoi(s),k,l));
                        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                        {
                            in >> s;
                            auto pos = s.find(',');
                            ff[atomIdx].real(stod(s.substr(0, pos)));
                            ff[atomIdx].imag(stod(s.substr(pos+1)));
                        }
                        atomicFormFactors.push_back(ff);
                    }

                } while (run);

            in.close();
        }


        void write_tsc(
            const std::string& fileName,
            const std::vector<std::string>& atomLabels,
            const std::vector<Vector3i>& hkl,
            const std::vector<std::vector<std::complex<double> > >& atomicFormFactors,
            const std::string& additionalText)
        {
            ofstream out(fileName);

            if (!out.good())
                on_error::throwException(string("can not write olex tsc file '")+fileName+string("'"), __FILE__, __LINE__);


            out << "TITLE: " << "form factors generated with discamb\n";

            out << "AD: FALSE\nSYMM: expanded\nSCATTERERS:";
            for (auto const& label : atomLabels)
                out << " " << label;
            out << "\nDATA:\n";
            out << additionalText;
            vector<complex<double> > formFactors;
            int hklIdx, nHkl, atomIdx, nAtoms = atomLabels.size();
            vector<bool> includeAtom(nAtoms, true);
            int maxNChars = maxNChar(hkl);
            int maxSizeFfStr = maxNChar(atomicFormFactors);
            nHkl = hkl.size();

            for(hklIdx=0; hklIdx<nHkl; hklIdx++)
            {
                //out << setw(maxNChars+1) << hkl[hklIdx][0] << setw(maxNChars + 1) << hkl[hklIdx][1] << setw(maxNChars + 1) << hkl[hklIdx][2];
                out << hkl[hklIdx][0] << " " << hkl[hklIdx][1] << " " << hkl[hklIdx][2];
                for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    string s = formFactorAsString(atomicFormFactors[hklIdx][atomIdx]);
                    //out << setw(maxSizeFfStr + 1) << s;
                    out << " " << s;
                }

                out << "\n";
            }

            out.close();

        }

        void write_tscd(
            const std::string& fileName,
            const std::vector<std::string>& atomLabels,
            const std::vector<Vector3d>& hkl,
            const std::vector<std::vector<std::complex<double> > >& atomicFormFactors,
            const std::string& additionalText)
        {
            ofstream out(fileName);

            if (!out.good())
                on_error::throwException(string("can not write olex tsc file '") + fileName + string("'"), __FILE__, __LINE__);


            out << "TITLE: " << "form factors generated with discamb\n";

            out << "AD: FALSE\nSYMM: expanded\nSCATTERERS:";
            for (auto const& label : atomLabels)
                out << " " << label;
            out << "\nDATA:\n";
            out << additionalText;
            vector<complex<double> > formFactors;
            int hklIdx, nHkl, atomIdx, nAtoms = atomLabels.size();
            vector<bool> includeAtom(nAtoms, true);
            
            int maxSizeFfStr = maxNChar(atomicFormFactors);
            nHkl = hkl.size();

            for (hklIdx = 0; hklIdx < nHkl; hklIdx++)
            {
                out << hkl[hklIdx][0] << " " << hkl[hklIdx][1] << " " << hkl[hklIdx][2] << " ";
                for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                {
                    string s = formFactorAsString(atomicFormFactors[hklIdx][atomIdx]);
                    out << setw(maxSizeFfStr + 3) << s;
                }

                out << "\n";
            }

            out.close();

        }

    }
}
