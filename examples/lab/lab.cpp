#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/IO/fragmentation_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"

using namespace std;
using namespace discamb;

void test_taam_disorder(int argc, char *argv[])
{
    if (argc != 5)
    {
        puts("expected 4 arguments:\n"
            " (1) cif/res file\n"
            " (2) bank path\n"
            " (3) fragment definitions\n"
            " (4) output file\n");
        exit(0);
    }
    string structureFile = argv[1];
    string bankPath = argv[2];
    string fragmentDefinitionsFile = argv[3];
    string outputFile = argv[4];

    // structure

    Crystal crystal;
    structure_io::read_structure(structureFile, crystal);

    // atom types

    MATTS_BankReader bankReader;
    vector<AtomType> atomTypes;
    vector<AtomTypeHC_Parameters> typeParameters;
    BankSettings bankSettings;

    bankReader.read(bankPath, atomTypes, typeParameters, bankSettings);

    // read fragments
    //vector<FragmentData> _fragments;
    //fragmentation_io::readPlainText(fragmentDefinitionsFile, _fragments);
    //UnitCellContent unitCellContent(crystal);

    //// convert fragmetns
    //vector< vector<AtomInCrystalID> > fragments(_fragments.size());
    //for (int i=0;i<fragments.size();i++)
    //    for(int j=0;j< fragments[i].size();j++)
    //    {
    //        fragments[i].push_back(AtomInCrystalID(_fragments[i].))
    //    }
    //// assign types

    //CrystalAtomTypeAssigner assigner;
    //assigner.setAtomTypes(atomTypes);
    //assigner.setDescriptorsSettings(bankSettings.descriptorsSettings);
    //assigner.assign(crystal, )
}
 
int main(int argc, char* argv[])
{
    try { 
     

    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
}

