#include "discamb/Scattering/MATTS_Default.h"

#include <vector>

using namespace std;



namespace discamb{
	void default_ubdb_bank_string(std::string &s)
	{

		//s = 
		vector<string> v = {
"",
"SETTINGS",
"  covalent bond threshold 0.23",
"  atom planarity threshold 0.1",
"  ring planarity threshold 0.1",
"  atom in ring planarity threshold 0.1",
"  atom in planar ring max number of neighbours 3",
"  min plm 0.002",
"  n sigma 1.0",
"  minimal number of type instances 2",
"",
"",
"ENTRY",
"  ID",
"    H101",
"  COMMENT",
"    in UBDB2018: H101 KJ  ",
"  NOI",
"    12234",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2                PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,H,H,!H         PLANARITY *  PLANAR_RING_WITH_PLANAR_ATOMS -        IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 X !H(C2) R",
"  SYMMETRY",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Sun Nov 12 12:51:24 2023",
"  MULTIPOLE MODEL PARAMETERS",
"    PVAL          1.092(38) KAPPA         1.091(13) KPRIM         1.144(18)",
"    PLMS   1  0    0.1866(94) PLMS   2  0    0.0857(76)",
"",
""
        };

	s.clear();
	for (auto const& line : v)
		s += line + string("\n");
	}
}
