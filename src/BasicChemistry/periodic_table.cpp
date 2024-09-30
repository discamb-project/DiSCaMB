#include "discamb/BasicChemistry/periodic_table.h"

#include "discamb/BasicUtilities/string_utilities.h"

using namespace std;

namespace {

    discamb::periodic_table::Block char2block(char c)
    {
        if (c == 'S')
            return discamb::periodic_table::Block::S;
        if (c == 'P')
            return discamb::periodic_table::Block::P;
        if (c == 'D')
            return discamb::periodic_table::Block::D;
        if (c == 'F')
            return discamb::periodic_table::Block::F;

        return discamb::periodic_table::Block::NONE;

    }

    const char* symbols[119] = {
"X","H","He","Li","Be","B","C","N","O","F",
"Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K",
"Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu",
"Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y",
"Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In",
"Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr",
"Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm",
"Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au",
"Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac",
"Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es",
"Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt",
"Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og" };

    char block_[119] = {
    'X','S','S','S','S','P','P','P','P','P',
    'P','S','S','P','P','P','P','P','P','S',
    'S','D','D','D','D','D','D','D','D','D',
    'D','P','P','P','P','P','P','S','S','D',
    'D','D','D','D','D','D','D','D','D','P',
    'P','P','P','P','P','S','S','F','F','F',
    'F','F','F','F','F','F','F','F','F','F',
    'F','F','D','D','D','D','D','D','D','D',
    'D','P','P','P','P','P','P','S','S','F',
    'F','F','F','F','F','F','F','F','F','F',
    'F','F','F','F','D','D','D','D','D','D',
    'D','D','D','P','P','P','P','P','P' };

    int period_[119] = {
    0, 1, 1, 2, 2, 2, 2, 2, 2, 2,
    2, 3, 3, 3, 3, 3, 3, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 5, 5, 5, 6, 6, 5, 5, 5,
    5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
    5, 5, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 6, 6, 6, 7, 7, 6,
    6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 7, 7, 7, 7, 7, 7,
    7, 7, 7, 7, 7, 7, 7, 7, 7 };

    int group_[119] = {
    0, 1, 18, 1, 2, 13, 14, 15, 16, 17,
    18, 1, 2, 13, 14, 15, 16, 17, 18, 1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 1, 2, 3,
    4, 5, 6, 7, 8, 9, 10, 11, 12, 13,
    14, 15, 16, 17, 18, 1, 2, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 4, 5, 6, 7, 8, 9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 1, 2, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 4, 5, 6, 7, 8, 9,
    10, 11, 12, 13, 14, 15, 16, 17, 18 };




    struct Element {
        std::string symbol;
        char block = 'S';
        int period = 1;
        int group = 1;
        int atomic_number = 1;
    };

    struct Data {
        std::map<std::string, int> symbolToAtomicNumber;
        std::vector<Element> elements;
        Data();
    };


    Data::Data()
    {
        string pt = 
            "   H   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -  He  " // 1
            "  Li  Be   -   -   -   -   -   -   -   -   -   -   B   C   N   O   F  Ne  " // 2
            "  Na  Mg   -   -   -   -   -   -   -   -   -   -  Al  Si   P   S  Cl  Ar  " // 3
            "   K  Ca  Sc  Ti   V  Cr  Mn  Fe  Co  Ni  Cu  Zn  Ga  Ge  As  Se  Br  Kr  " // 4
            "  Rb  Sr   Y  Zr  Nb  Mo  Tc  Ru  Rh  Pd  Ag  Cd  In  Sn  Sb  Te   I  Xe  " // 5
            "  Cs  Ba   -  Hf  Ta   W  Re  Os  Ir  Pt  Au  Hg  Tl  Pb  Bi  Po  At  Rn  " // 6
            "  Fr  Ra   -  Rf  Db  Sg  Bh  Hs  Mt  Ds  Rg  Cn  Nh  Fl  Mc  Lv  Ts  Og  " // 7
            "   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -  "
            "   -   -   -  La  Ce  Pr  Nd  Pm  Sm  Eu  Gd  Tb  Dy  Ho  Er  Tm  Yb  Lu  "
            "   -   -   -  Ac  Th  Pa   U  Np  Pu  Am  Cm  Bk  Cf  Es  Fm  Md  No  Lr  ";

        vector<string> symbols;
        discamb::string_utilities::split(pt, symbols);
        Element element;
        int atomic_number, index;
        elements.resize(119);

        elements[0].symbol = string("X");
        elements[0].atomic_number = 0;
        elements[0].block = 'X';
        elements[0].period = 0;
        elements[0].group = 0;


        atomic_number = 1;
        index = 0;
        for (int period = 1; period < 6; period++)
            for (int group = 1; group < 19; group++)
            {
                if (symbols[index] != string("-"))
                {
                    element.symbol = symbols[index];
                    element.atomic_number = atomic_number;
                    element.block = (group < 3 ? 'S' : (group > 12 ? 'P' : 'D'));
                    if (element.symbol == string("He"))
                        element.block = 'S';
                    element.group = group;
                    element.period = period;
                    elements[atomic_number] = element;
                    symbolToAtomicNumber[symbols[index]] = atomic_number;
                    atomic_number++;
                }
                index++;
            }

        for (int period = 6; period < 8; period++)
            for (int group = 1; group < 19; group++)
            {
                if (index == 92)
                    atomic_number = 72;
                if (index == 110)
                    atomic_number = 104;

                if (symbols[index] != string("-"))
                {
                    element.symbol = symbols[index];
                    element.atomic_number = atomic_number++;
                    element.block = (group < 3 ? 'S' : (group > 12 ? 'P' : 'D'));
                    element.group = group;
                    element.period = period;
                    //elements.push_back(element);
                    elements[element.atomic_number] = element;
                    symbolToAtomicNumber[symbols[index]] = element.atomic_number;
                }
                index++;
            }

        index += 18;
        atomic_number = 57;


        for (int period = 5; period < 7; period++)
            for (int group = 1; group < 19; group++)
            {
                if (atomic_number == 72)
                    atomic_number = 89;

                if (symbols[index] != string("-"))
                {
                    element.symbol = symbols[index];
                    element.atomic_number = atomic_number++;
                    element.block = 'F';
                    element.group = 0;
                    element.period = period;
                    //elements.push_back(element);
                    elements[element.atomic_number] = element;
                    symbolToAtomicNumber[symbols[index]] = element.atomic_number;
                }
                index++;
            }

    }

    Data periodicTableData;


}


namespace discamb {

    namespace periodic_table
    {

        std::string symbol(
            int atomicNumber)
        {
            if (atomicNumber < 119)
                return symbols[atomicNumber];
                //return periodicTableData.elements[atomicNumber].symbol;
            return "X";
        }

        int atomicNumber(
            const std::string &_symbol, 
            bool caseSensitive)
        {
            string symbol = _symbol;
            if (symbol.empty())
                return 0;

            if (!caseSensitive)
                symbol[0] = toupper(symbol[0]);
            for(int i = 1; i< symbol.size(); i++)
                symbol[i] = tolower(symbol[i]);
                
            for (int i = 0; i < 119; i++)
                if (string(symbols[i]) == symbol)
                    return i;
            return 0;
            if (periodicTableData.symbolToAtomicNumber.find(symbol) != periodicTableData.symbolToAtomicNumber.end())
                return periodicTableData.symbolToAtomicNumber.find(symbol)->second;
            if (symbol == string("D"))
                return 1;
            return 0;
        }

        Block block(
            int atomicNumber)
        {
            if (atomicNumber < 119)
                return char2block(block_[atomicNumber]);
                //return periodicTableData.elements[atomicNumber].block;
            return Block::NONE;
        }

        char block2str(Block b)
        {
            if (b == Block::S)
                return 'S';
            if (b == Block::P)
                return 'P';
            if (b == Block::D)
                return 'D';
            if (b == Block::F)
                return 'F';

            return 'X';
        }

        int period(
            int atomicNumber)
        {
            if (atomicNumber < 119)
                return period_[atomicNumber];
                //return periodicTableData.elements[atomicNumber].period;
            return 0;
        }

        int group(
            int atomicNumber)
        {
            if (atomicNumber < 119)
                return group_[atomicNumber];
                //return periodicTableData.elements[atomicNumber].group;
            return 0;

        }

    }
}