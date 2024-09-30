#include "discamb/IO/xyz_io.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/BasicUtilities/on_error.h"

#include <fstream>
#include <iomanip>

using namespace std;

namespace {
    void atomicNumbersToSymbols(
        const vector<int> &atomicNumbers,
        vector<std::string> &symbols)
    {
        symbols.clear();
        for (int z : atomicNumbers)
            symbols.push_back(discamb::periodic_table::symbol(z));
    }

    void symbolsToAtomicNumbers(
        const vector<std::string> &symbols,
        vector<int> &atomicNumbers)
    {
        atomicNumbers.clear();
        for (auto & symbol : symbols)
            atomicNumbers.push_back(discamb::periodic_table::atomicNumber(symbol));
    }

}

namespace discamb {
    namespace xyz_io {

        void readXyz(
            const std::string &fileName,
            std::vector<std::string> &symbol,
            std::vector<Vector3d> &position)
        {
            string line;
            ifstream in(fileName);

            if (!in.good())
                on_error::throwException(string("can not read *.xyz file : '") + fileName + string("'"), __FILE__, __LINE__);
            int n;
            
            in >> n;
            symbol.resize(n);
            position.resize(n);

            getline(in, line);
            getline(in, line);

            for (int i = 0; i < n; i++)
                in >> symbol[i] >> position[i][0] >> position[i][1] >> position[i][2];
            in.close();
        }

        void readXyz(
            const std::string &fileName,
            std::vector<int> &atomicNumber,
            std::vector<Vector3d> &position)
        {
            vector<string> symbols;
            readXyz(fileName, symbols, position);
            symbolsToAtomicNumbers(symbols, atomicNumber);
        }

        void readXyz(
            const std::string& fileName,
            std::vector<ChemicalElement>& element,
            std::vector<Vector3d>& position)
        {
            vector<string> symbol;
            readXyz(fileName, symbol, position);
            element.clear();
            for (auto& s : symbol)
                element.push_back(ChemicalElement(s));
        }

        void readXyz(
            const std::string& fileName,
            MoleculeData& data)
        {
            data = MoleculeData();
            
            string line, s;
            ifstream in(fileName);

            if (!in.good())
                on_error::throwException(string("can not read *.xyz file : '") + fileName + string("'"), __FILE__, __LINE__);
            int n;

            in >> n;
            data.atomicNumbers.resize(n);
            data.atomPositions.resize(n);

            getline(in, line);
            getline(in, data.comment);

            for (int i = 0; i < n; i++)
            {
                in >> s >> data.atomPositions[i][0] >> data.atomPositions[i][1] >> data.atomPositions[i][2];
                data.atomicNumbers[i] = periodic_table::atomicNumber(s);
            }
            in.close();
        }

        void writeXyz(
            const std::string &fileName,
            const  std::vector<std::string> &symbol,
            const std::vector<Vector3d> &position)
        {
            ofstream out(fileName);
            if (!out.good())
                on_error::throwException(string("can not writ to *.xyz file : '") + fileName + string("'"), __FILE__, __LINE__);
            int i, n = symbol.size();
            out << n << endl << "comment line" << endl;
            for (i = 0; i < n; i++)
                out << symbol[i]
                    << " " << setw(16) << setprecision(7) << fixed << position[i][0]
                    << " " << setw(16) << setprecision(7) << fixed << position[i][1]
                    << " " << setw(16) << setprecision(7) << fixed << position[i][2] << endl;
            out.close();
        }

        void writeXyz(
            const std::string &fileName,
            const std::vector<int> &atomicNumber,
            const std::vector<Vector3d> &position)
        {
            vector<string> symbol;
            atomicNumbersToSymbols(atomicNumber, symbol);
            writeXyz(fileName, symbol, position);
        }

        void writeXyz(
            const std::string& fileName,
            const std::vector<ChemicalElement>& element,
            const std::vector<Vector3d>& position)
        {
            vector<string> symbol;
            for (auto e : element)
                symbol.push_back(e.symbol());
            writeXyz(fileName, symbol, position);
        }

        void writeXyz(
            const std::string& fileName,
            const MoleculeData &data)
        {
            ofstream out(fileName);
            if (!out.good())
                on_error::throwException(string("can not writ to *.xyz file : '") + fileName + string("'"), __FILE__, __LINE__);
            int i, n = data.atomicNumbers.size();
            out << n << endl << "comment line" << endl;
            for (i = 0; i < n; i++)
                out << periodic_table::symbol(data.atomicNumbers[i])
                << " " << setw(16) << setprecision(7) << fixed << data.atomPositions[i][0]
                << " " << setw(16) << setprecision(7) << fixed << data.atomPositions[i][1]
                << " " << setw(16) << setprecision(7) << fixed << data.atomPositions[i][2] << endl;
            out.close();
        }

    }
}
