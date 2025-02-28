#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicChemistry/chemical_element_data.h"

#include "discamb/BasicUtilities/constants.h"
#include "discamb/BasicUtilities/on_error.h"

#include "discamb/CrystalStructure/UnitCellContent.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"

#include "discamb/IO/xyz_io.h"

#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/MathUtilities/geometry3d.h"

#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/StructuralProperties/CovalentRadiousBondDetector.h"
#include "discamb/StructuralProperties/GenericConnectivityAlgorithm.h"

#include "discamb/MathUtilities/graph_algorithms.h"
#include "discamb/IO/MATTS_BankReader.h"
#include "discamb/AtomTyping/MolecularAtomTypeAssigner.h"

#include "argedit.h"
#include "argraph.h"
#include "argedit.h"
#include "vf2_state.h"
#include "match.h"

#include <set>

//--------
#include <fstream>
#include "discamb/BasicChemistry/periodic_table.h"
//--------

using namespace std;


namespace discamb {


namespace 
{
    const vector<string> hydrogenTypesBank = {
"",
"SETTINGS",
"  covalent bond threshold 0.4",
"  atom planarity threshold 0.1",
"  ring planarity threshold 0.1",
"  atom in ring planarity threshold 0.1",
"  atom in planar ring max number of neighbours 3",
"  min plm 0.002",
"  n sigma 1.0",
"  minimal number of type instances 3",
"",
"",
"ENTRY",
"  ID",
"    C(cyclopropyl)-H",
"  COMMENT ",
"    C(cyclopropyl)-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,C3,C4,X    PLANARITY  -  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING + IN_4_MEMBER_RING -",
"    C3    CONNECTED_TO  C2,C4,X,*     PLANARITY  -  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING + IN_4_MEMBER_RING -",
"    C4    CONNECTED_TO  C2,C3,X,*     PLANARITY  -  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING + IN_4_MEMBER_RING -",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    C3-Csp3-H",
"  COMMENT ",
"    C3-Csp3-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,C,C,C      PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    C2-Csp3-H2",
"  COMMENT ",
"    C2-Csp3-H2 ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,C,C,H      PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    C-Csp3-H3",
"  COMMENT ",
"    C-Csp3-H3 ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,C,H,H      PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    Z-Csp3-H3",
"  COMMENT ",
"    C-Csp3-H3 ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,!H,H,H      PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    C=Csp2-H",
"  COMMENT ",
"    C=Csp2-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,C3,X       PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  C2,X,X        PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    C#Csp1-H",
"  COMMENT ",
"    C#Csp1-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,C3         PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  C2,X          PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",

"",
"ENTRY",
"  ID",
"    C(ar)-H",
"  COMMENT ",
"    C(ar)-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,!H,!H      PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS +       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    Z2-Csp3-H2",
"  COMMENT ",
"    C2-Csp3-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,!H,!H,H      PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.06419(66) KAPPA       1.09836(23) KPRIM       1.14173(41)",
"    PLMS   1  0     0.179(14) PLMS   2  0     0.085(10)",
"",
"ENTRY",
"  ID",
"    Z3-Csp3-H",
"  COMMENT ",
"    Z3-Csp3-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,!H,!H,!H     PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y !H(C2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    H-O-H",
"  COMMENT ",
"    H-O-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,H3           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    H3    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y H3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Csp3-O-H",
"  COMMENT ",
"    Csp3-O-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,C3           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  O2,X,X,X        PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    C(ar)-O-H",
"  COMMENT ",
"    C(ar)-O-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,C3           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  O2,X,X          PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS +       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    O=Csp2-O-H",
"  COMMENT ",
"    O=Csp2-O-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,C3           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  O2,O3,X         PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    O3    CONNECTED_TO  C3              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    C(any)-O-H",
"  COMMENT ",
"    C(any)-O-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,C3           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  O2,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Z-O-H",
"  COMMENT ",
"    Z-O-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,!H           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y !H(O2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Z-N-H",
"  COMMENT ",
"    Z-N-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,!H           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y !H(N2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    N+-H",
"  COMMENT ",
"    N+-H ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,X3,X,X        PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    X3    CONNECTED_TO  N2,*           PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y X3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Csp3-N-H2",
"  COMMENT ",
"    Csp3-N-H2 ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,H,C3         PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  N2,X,X,X        PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Csp2-N-H2",
"  COMMENT ",
"    Csp2-N-H2 ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,H,C3         PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  N2,X,X          PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Csp2-N(pl)-H2",
"  COMMENT ",
"    Csp2-N(pl)-H2 ",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,H,C3         PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  N2,X,X          PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Csp2-N-H2(amido)",
"  COMMENT ",
"    Csp2-N-H2(amido)",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,H,C3         PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  N2,O4,X         PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    O4    CONNECTED_TO  C3              PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    C(ar)-N-H2",
"  COMMENT ",
"    Csp2-N-H2(amido)",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,H,C3         PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  N2,X,X          PLANARITY  +  PLANAR_RING_WITH_PLANAR_ATOMS +       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Z2-N-H",
"  COMMENT ",
"    Z2-N-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    N2    CONNECTED_TO  H1,!H,!H        PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y !H(N2) R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    (Csp2)2-N-H",
"  COMMENT ",
"    (Csp2)2-N-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,C3,C4        PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C3    CONNECTED_TO  N2,X,X          PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS +       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"    C4    CONNECTED_TO  N2,X,X          PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS +       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y C3 R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    B-H",
"  COMMENT ",
"    B-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  B2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    B2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z B2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    C-H",
"  COMMENT ",
"    C-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  C2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    C2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z C2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    N-H",
"  COMMENT ",
"    N-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  N2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    N2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z N2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    O-H",
"  COMMENT ",
"    O-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  O2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    O2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z O2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    P-H",
"  COMMENT ",
"    P-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  P2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    P2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z P2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    S-H",
"  COMMENT ",
"    S-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  S2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    S2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z S2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
"",
"ENTRY",
"  ID",
"    Si-H",
"  COMMENT ",
"    Si-H",
"  ATOM DESCRIPTORS",
"    H1    CONNECTED_TO  Si2              PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS -       IN_3_MEMBER_RING - IN_4_MEMBER_RING -",
"    Si2    CONNECTED_TO  H1,*            PLANARITY  *  PLANAR_RING_WITH_PLANAR_ATOMS *       IN_3_MEMBER_RING * IN_4_MEMBER_RING *",
"  LOCAL COORDINATE SYSTEM",
"    Z Si2 Y any_orthogonal R",
"  SYMMETRY  ",
"    cyl",
"  PARAMETER MODIFICATION DATE",
"    Fri Mar  6 15:30:55 2020",
"  ENTRY CREATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  DEFINITION MODIFICATION DATE",
"    Tue Apr 30 17:40:04 2019",
"  MULTIPOLE MODEL PARAMTERS",
"    PVAL        1.08802(47) KAPPA       1.09322(19) KPRIM       1.14868(29)",
"    PLMS   1  0     0.184(17) PLMS   2  0     0.085(11)",
""
    };

    const map<string, double> hydrogenType2BondLength = {
    {"C(cyclopropyl)-H",1.08}, // Allen 2010
    {"C3-Csp3-H",1.099},       // Allen 2010
    {"C2-Csp3-H2",1.092},      // Allen 2010
    {"C-Csp3-H3",1.077},       // Allen 2010
    {"Z-Csp3-H3",1.077},       // Allen 2010
    {"C=Csp2-H",1.082},        // Allen 2010
    {"C#Csp1-H",1.055},        // Allen 2010
    {"C(ar)-H",1.083},         // Allen 2010
    {"Z2-Csp3-H2",1.091},      // Allen 2010
    {"Z3-Csp3-H",1.098},       // Allen 2010
    {"H-O-H",0.959},           // Woiñska et al. (2016, Sci. Adv. 2, e1600192) 
    {"Csp3-O-H",0.97},         // Allen 2010
    {"C(ar)-O-H",0.992},       // Allen 2010
    {"O=Csp2-O-H",1.018},      // Allen 2010
    {"C(any)-O-H",0.98},       // Allen 2010
    {"Z-O-H",0.983},           // Allen 2010
    {"Z-N-H",1.015},           // Allen 2010
    {"N+-H",1.036},            // Allen 2010
    {"Csp3-N-H2",1.002},       // Allen 2010
    {"Csp2-N-H2",1.013},       // Allen 2010
    {"Csp2-N(pl)-H2",1.012},   // Allen 2010
    {"Csp2-N-H2(amido)",1.01}, // Allen 2010
    {"C(ar)-N-H2",1.011},      // Allen 2010
    {"Z2-N-H",1.027},          // Allen 2010
    {"(Csp2)2-N-H",1.030},     // Allen 2010
    {"B-H",1.185},             // Allen 2010
    {"C-H",1.089},             // Mercury 2020.2.0
    {"N-H",1.015},             // Mercury 2020.2.0
    {"O-H",0.993},             // Mercury 2020.2.0
    {"P-H",1.389},             // average from CSD and ICSD from 6 neutron structures with all atoms U_aniso (13.09.2023)
    {"S-H",1.338},             // Allen 2010
    {"Si-H",1.506}             // Allen 2010
    };



	void makeCellsShell(
		int shell,
		set<Vector3i>& cells)
	{
		cells.clear();
		Vector3i x(1, 0, 0), y(0, 1, 0), z(0, 0, 1);

		for (int i = -shell; i <= shell; i++)
			for (int j = -shell; j <= shell; j++)
			{
				cells.insert(shell * x + i * y + j * z);
				cells.insert(-shell * x + i * y + j * z);

				cells.insert(shell * y + i * x + j * z);
				cells.insert(-shell * y + i * x + j * z);

				cells.insert(shell * z + i * y + j * x);
				cells.insert(-shell * z + i * y + j * x);
			}

	}

	void calculateShellMolecules(
		int shell,
		const vector<vector<Vector3i> >& moleculesIntersectingUnitCell,
		const set<pair<int, Vector3i> >& previousShellMolecules,
		set<pair<int, Vector3i> >& shellMolecules)
	{
		shellMolecules.clear();

		// get molecules intersecting unit cell in useful format

		vector< pair<int, Vector3i> > cellIntersectingMolecules;

		for (int i = 0; i < moleculesIntersectingUnitCell.size(); i++)
			for (auto& cell : moleculesIntersectingUnitCell[i])
				cellIntersectingMolecules.push_back({ i,cell });

		// find cells forming shell

		set<Vector3i> shellCells;
		makeCellsShell(shell, shellCells);

		// find molecules intersecting shell

		for (auto& cell : shellCells)
			for (auto& molecule : cellIntersectingMolecules)
				shellMolecules.insert({ molecule.first, molecule.second + cell });

		// remove molecules from previous (sub)shell

		for (auto& molecule : previousShellMolecules)
		{
			auto it = shellMolecules.find(molecule);
			if (it != shellMolecules.end())
				shellMolecules.erase(it);
		}

	}

	void findClusterExtendingMolecules(
		const vector<Vector3d>& centralMolAtomicCoords,
        const vector<int> &centralMolZ, // central molecule atomic numbers
		const vector<vector<Vector3d> >& moleculesAtomicCoords,
        const vector < vector <int> > &moleculesZ, // atomic numbers
		const vector<Vector3d>& latticeVec,
		const set<pair<int, Vector3i> >& shellMolecules,
		set<pair<int, Vector3i> >& shellMoleculesclusterExtension,
		double threshold,
        bool useVdW)
	{
		Vector3d rShellAtom, cellShift, diff;
		int moleculeIdx;
		bool moleculeAddedToClusterExtension;
		shellMoleculesclusterExtension.clear();
		for (auto const& mol : shellMolecules)
		{
			moleculeAddedToClusterExtension = false;
			moleculeIdx = mol.first;
			cellShift = double(mol.second[0]) * latticeVec[0] + double(mol.second[1]) * latticeVec[1] + double(mol.second[2]) * latticeVec[2];
			double threshold2 = threshold * threshold;
            for(int i=0; i< moleculesAtomicCoords[moleculeIdx].size(); i++)
			//for (auto const& r : moleculesAtomicCoords[moleculeIdx])
			{
                auto const& r = moleculesAtomicCoords[moleculeIdx][i];
				rShellAtom = r + cellShift;
                for(int j=0;j< centralMolAtomicCoords.size(); j++)
				//for (auto const& rCentralMolAtom : centralMolAtomicCoords)
				{
                    auto const& rCentralMolAtom = centralMolAtomicCoords[j];
					diff = rShellAtom - rCentralMolAtom;
                    bool include = false;
                    if (useVdW)
                    {
                        double vanDerWaalsRadiiSum = chemical_element_data::vdwRadius(moleculesZ[moleculeIdx][i]) +
                                                     chemical_element_data::vdwRadius(centralMolZ[j]);
                        include = threshold * vanDerWaalsRadiiSum >= sqrt(diff * diff);
                    }
                    else
                        include = threshold2 >= diff * diff;

					if (include)
					{
						shellMoleculesclusterExtension.insert(mol);
						moleculeAddedToClusterExtension = true;
						break;
					}

				}
				if (moleculeAddedToClusterExtension)
					break;
			}
		}
	}

}

namespace structural_properties {

    /*
Vector3d StockholderAtomFormFactorCalcManager::capAtomPosition(
        const Crystal& crystal,
        int bondedAtom,
        const SpaceGroupOperation& bondedAtomSymmOp,
        int directingAtom,
        const SpaceGroupOperation& directingAtomSymmOp)
    {
        //vector<double> rH{ 1.089, 1.015, 0.993 };

        map<string, double> rH;

        for (int i = 1; i < 113; i++)
            rH[periodic_table::symbol(i)] = chemical_element_data::covalentRadius(1) + chemical_element_data::covalentRadius(i);

        map<string, double> rH_standrized { {"B", 1.185}, { "C", 1.089 }, { "N",1.015 }, { "O", 0.993 }, { "Si",1.506 }, { "P",1.42 }, { "S",1.338 } };

        for (auto x : rH_standrized)
            rH[x.first] = x.second;

        double bond_to_H_length;
        Vector3d a1, a2, a12, norm_a12;

        a1 = atomPosition(bondedAtom, bondedAtomSymmOp, crystal);

        a2 = atomPosition(directingAtom, directingAtomSymmOp, crystal);

        a12 = a2 - a1;

        norm_a12 = a12 / sqrt(a12 * a12);

        string symbol = crystal.atoms[bondedAtom].type;

        bond_to_H_length = 1.0;

        if (rH.find(symbol) != rH.end())
            bond_to_H_length = rH[symbol];
        else
            on_error::throwException("standarized bond length for capped hydrogen atom not available for '" + symbol + "'", __FILE__, __LINE__);

        //int atomicNumber = periodic_table::atomicNumber(crystal.atoms[bondedAtom].type);

        //if (atomicNumber < 6 || atomicNumber > 8)
        //	on_error::throwException(string("standarized bond length for capped hydrogen atom not available for atomic number ")
        //		+ to_string(atomicNumber), __FILE__, __LINE__);

        return a1 + norm_a12 * bond_to_H_length;

    }
    */

    Vector3d capping_atom_position(
        const Vector3d& bonded_atom_r,
        const Vector3d& directing_atom_r,
        int bonding_atom_atomic_number)
    {
        double bond_to_H_length;
        Vector3d a12, norm_a12, bond_direction_normalized;
        map<int, double> rH_standrized{ {5, 1.185}, { 6, 1.089 }, { 7,1.015 }, { 8, 0.993 }, { 14, 1.506 }, { 15, 1.42 }, { 16, 1.338 } };

        if (rH_standrized.find(bonding_atom_atomic_number) != rH_standrized.end())
            bond_to_H_length = rH_standrized[bonding_atom_atomic_number];
        else
            bond_to_H_length = chemical_element_data::covalentRadius(1) + chemical_element_data::covalentRadius(bonding_atom_atomic_number);
  
        bond_direction_normalized  = directing_atom_r - bonded_atom_r;
        bond_direction_normalized = bond_direction_normalized / sqrt(bond_direction_normalized * bond_direction_normalized);

        return bonded_atom_r + bond_direction_normalized * bond_to_H_length;
    }

    /*void asymmetricUnitConnectivity(
        const Crystal &c,
        std::vector<std::vector<std::pair<int, std::string> > > &connectivity)
    {

    }*/

    /**
    Algorith for finding plane with minimal rmsd of points' distance from the plane as described
    in A.G.URZHUMTSEV 'How to Calculate Planarity Restraints' Acta Cryst. (1991). A47,723-727
    */

    void find_plane(
        const std::vector<Vector3d>& positions, 
        Vector3d& average, 
        Vector3d& normal)
    {
        int i, j, pointIdx, nPoints = positions.size();
        if (nPoints < 3)
            return;

        

        Matrix3<double> matrix;
        Vector3d position;
        Vector3d eigenVectors[3], planeNormal;
        vector<pair<double, int> > eigenValues(3);

        // clculate average position

        for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
            average += positions[pointIdx];
        average /= (double) nPoints;

        // calculate 'displacement matrix'

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
                    matrix(i, j) += (positions[pointIdx](i) - average(i)) * (positions[pointIdx](j) - average(j));

        // find smallest eigenvector of the matrix = normal to the 'plane'

        for (i = 0; i < 3; i++)
            eigenValues[i].second = i;
        algebra3d::eigensystemRealSymm(matrix, eigenVectors[0], eigenVectors[1], eigenVectors[2], eigenValues[0].first, eigenValues[1].first, eigenValues[2].first);

        planeNormal = eigenVectors[std::min_element(eigenValues.begin(), eigenValues.end())->second];

        planeNormal /= sqrt(planeNormal * planeNormal); // normalize



    }

    /**
    Calculates estimated standard deviation of point distance from plane.
    Algorith for finding the plane is described in A.G.URZHUMTSEV 'How to Calculate Planarity Restraints' Acta Cryst. (1991). A47,723-727
    */
    double planarity_esd(const std::vector<Vector3d>& positions)
    {

        Vector3d average, normal;

        find_plane(positions, average, normal);

        // calculate sigma of atom distance to plane

        double averageOntoPlaneNormalProjection = average * normal;
        double distance_esd = 0;

        for (auto const &r: positions)
        {
            double distance = normal * r - averageOntoPlaneNormalProjection;
            distance_esd += distance * distance;
        }

        distance_esd = sqrt(distance_esd / double(positions.size() - 3));

        return distance_esd;

    }




    void calculateConnectivity(
        const std::vector<Vector3d>& positions,
        const std::vector<int>& atomicNumbers,
        std::vector<std::vector<int> >& connectivity,
        double threshold)
    {
        GenericConnectivityAlgorithm<CovalentRadiousBondDetector> connectivityAlgorithm;
        connectivityAlgorithm.set(to_string(threshold));
        
        //MolecularDisorder disorder;

        connectivityAlgorithm.calculateConnectivity(positions, atomicNumbers, connectivity);
    }

    void normalizeXhBondLengths(
        std::vector<Vector3d>& positions,
        const std::vector<int>& atomicNumbers)
    {
        vector<string> labels;
        for (int i = 0; i < atomicNumbers.size(); i++)
            labels.push_back(periodic_table::symbol(atomicNumbers[i]) + to_string(i + 1));


        stringstream bankStream;
        MATTS_BankReader bankReader;
        vector<AtomType> atomTypes;
        vector<AtomTypeHC_Parameters> hcParameters;
        BankSettings bankSettings;
        //map<pair<string, string>, double> bondLengths;

        for (auto& line : hydrogenTypesBank)
            bankStream << line << "\n";
        bankReader.read(bankStream, atomTypes, hcParameters, bankSettings);

        MolecularAtomTypeAssigner assigner;
        assigner.setAtomTypes(atomTypes);
        assigner.setDescriptorsSettings(DescriptorsSettings());
        vector < LocalCoordinateSystem<int> > lcs;
        vector<int> types;
        assigner.assign(atomicNumbers, positions, labels, types, lcs);

        for (int i = 0; i < types.size(); i++)
        {

            //cout << periodic_table::symbol(atomicNumbers[i]) << i + 1 << " ";
            //if (types[i] >= 0)
            //    cout << atomTypes[types[i]].id << "\n";
            //else
            //    cout << "----" << "\n";
            if (atomicNumbers[i] == 1 && types[i] < 0)
                cout << i + 1 << "\n";
            
        }
                


        // bond lengths from types

        vector<double> newLengths(positions.size(), -1.0); // -1.0 = do not change
        for (int i = 0; i < types.size(); i++)
            if (types[i] >= 0)
                newLengths[i] = hydrogenType2BondLength.find(atomTypes[types[i]].id)->second;

        // change H positions

        GenericConnectivityAlgorithm< CovalentRadiousBondDetector> alg;
        vector<vector<int> > connectivity;
        alg.calculateConnectivity(positions, atomicNumbers, connectivity);


        Vector3d rX, rH, rXH, new_rH, new_rH_frac, frac, unitCellShift;
        double bondLength;

        for (int i = 0; i < newLengths.size(); i++)
            if (newLengths[i] > 0)
            {
                if (connectivity[i].size() != 1)
                    on_error::throwException("invalid connectivity for " + labels[i], __FILE__, __LINE__);

                rXH = positions[i] - positions[connectivity[i][0]];
                bondLength = sqrt(rXH * rXH);
                positions[i] = positions[connectivity[i][0]] + rXH * (newLengths[i] / bondLength);
            }

    }


    void asymmetricUnitConnectivity(
        const Crystal& c,
        std::vector<std::vector<std::pair<int, std::string> > >& connectivity,
        double threshold)
    {
        connectivity.clear();

        UnitCellContent unitCellContent;
        vector<vector<UnitCellContent::AtomID> > ucConnectivity;
        unitCellContent.set(c);
        calcUnitCellConnectivity(unitCellContent, ucConnectivity, threshold);
        int atomIdx, nAtoms = c.atoms.size();
        int neighbourIdx, nNeighbours;
        SpaceGroupOperation spaceGroupOperation, translationSymmOp, neighbourTranslationSymmOp;
        Vector3i latticeTranslation;
        Vector3<CrystallographicRational> translation;
        /*
            spaceGroupOperation = unitCellContent.getGeneratingOperation(atomIdx, 0);
            spaceGroupOperation.getTranslation(translation);
            latticeTranslation.set(-translation[0].numerator(), -translation[1].numerator(), -translation[2].numerator());
            asymmetricUnit.push_back(UnitCellContent::AtomID(atomIdx, latticeTranslation));
        */
        map<string, int> label2index;
        int index = 0;
        for (auto const& a : c.atoms)
            label2index[a.label] = index++;

        connectivity.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            nNeighbours = ucConnectivity[atomIdx].size();
            spaceGroupOperation = unitCellContent.getGeneratingOperation(atomIdx, 0);
            spaceGroupOperation.getTranslation(translation);
            translationSymmOp.setToTranslation(-translation);

            connectivity[atomIdx].resize(nNeighbours);

            for (neighbourIdx = 0; neighbourIdx < nNeighbours; neighbourIdx++)
            {
                string label, symmOpAsString;
                unitCellContent.interpreteAtomID(ucConnectivity[atomIdx][neighbourIdx].atomIndex, label, symmOpAsString);
                spaceGroupOperation.set(symmOpAsString);
                Vector3i &v = ucConnectivity[atomIdx][neighbourIdx].unitCellPosition;
                neighbourTranslationSymmOp.setToTranslation(v[0], 1, v[1], 1, v[2], 1);
                spaceGroupOperation = translationSymmOp * neighbourTranslationSymmOp * spaceGroupOperation;
                spaceGroupOperation.get(symmOpAsString);
                connectivity[atomIdx][neighbourIdx] = { label2index[label], symmOpAsString };
            }
        }
    }

    void assymetricUnitWithNeighbours(
        const Crystal &crystal,
        std::vector< std::pair<int, std::string> > &asuWithNeighbours,
        int neighbourRange,
        double threshold)
    {
        vector<int> shellSizes;
        assymetricUnitWithNeighbours(crystal, asuWithNeighbours, neighbourRange, threshold, shellSizes);

        //asuWithNeighbours.clear();

        //UnitCellContent unitCellContent;
        //vector<UnitCellContent::AtomID> asymmetricUnit, graph;
        //unitCellContent.set(crystal);
        //SpaceGroupOperation spaceGroupOperation;
        //int atomIdx, nAtoms = crystal.atoms.size();
        //Vector3<CrystallographicRational> translation;
        //Matrix3i rotation;
        //Vector3i latticeTranslation;

        //for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        //{
        //    spaceGroupOperation = unitCellContent.getGeneratingOperation(atomIdx, 0);
        //    spaceGroupOperation.getTranslation(translation);
        //    latticeTranslation.set(-translation[0].numerator(), -translation[1].numerator(), -translation[2].numerator());
        //    asymmetricUnit.push_back(UnitCellContent::AtomID(atomIdx, latticeTranslation));
        //}

        //structural_properties::graphToNthNeighbour(unitCellContent, asymmetricUnit, graph, 8, threshold);
        //
        //AtomInCrystal atomInCrystal;
        //string label, symmetryOperationAsString;
        //for (auto atom : graph)
        //{
        //    atomIdx = unitCellContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex);
        //    unitCellContent.interpreteAtomID(atom, label, symmetryOperationAsString);
        //    asuWithNeighbours.push_back({ atomIdx, symmetryOperationAsString });
        //}

    }

    void assymetricUnitWithNeighbours(const Crystal &crystal,
        std::vector< std::pair<int, std::string> > &asuWithNeighbours,
        int neighbourRange,
        double threshold,
        std::vector<int> &shellSizes)
    {
        asuWithNeighbours.clear();

        UnitCellContent unitCellContent;
        vector<UnitCellContent::AtomID> asymmetricUnit, graph;
        unitCellContent.set(crystal);
        SpaceGroupOperation spaceGroupOperation;
        int atomIdx, nAtoms = crystal.atoms.size();
        Vector3<CrystallographicRational> translation;
        Matrix3i rotation;
        Vector3i latticeTranslation;

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            spaceGroupOperation = unitCellContent.getGeneratingOperation(atomIdx, 0);
            spaceGroupOperation.getTranslation(translation);
            latticeTranslation.set(-translation[0].numerator(), -translation[1].numerator(), -translation[2].numerator());
            asymmetricUnit.push_back(UnitCellContent::AtomID(atomIdx, latticeTranslation));
        }

        //structural_properties::graphToNthNeighbour(unitCellContent, asymmetricUnit, graph, 8, threshold, shellSizes);
        structural_properties::graphToNthNeighbour(unitCellContent, asymmetricUnit, graph, neighbourRange, threshold, shellSizes);

        AtomInCrystal atomInCrystal;
        string label, symmetryOperationAsString;
        for (auto atom : graph)
        {
            atomIdx = unitCellContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex);
            unitCellContent.interpreteAtomID(atom, label, symmetryOperationAsString);
            asuWithNeighbours.push_back({ atomIdx, symmetryOperationAsString });
        }

    }


    void assymetricUnitWithNeighbours(
        const Crystal &crystal,
        std::vector<int> &atomicNumbers,
        std::vector<Vector3d> &positions,
        std::vector<std::string> &labels,
        int neighbourRange,
        double threshold)
    {
        vector<int> shellSizes;
        assymetricUnitWithNeighbours(crystal, atomicNumbers, positions, labels, neighbourRange, threshold, shellSizes);

        //int atomIdx, nAtoms;
        //Vector3d cartesianXyz, fractionalXyz, fractionalTransformedXyz;
    
        //atomicNumbers.clear();
        //positions.clear();
        //labels.clear();
        //SpaceGroupOperation symmetryOperation;
    
        //vector< pair<int, string> > asuWithNeighbours;
        //vector<int> atomsInAsymmetricUnitAtomicNumbers;
        //
        //structural_properties::assymetricUnitWithNeighbours(crystal, asuWithNeighbours, neighbourRange, threshold);
        //crystal_structure_utilities::atomicNumbers(crystal, atomsInAsymmetricUnitAtomicNumbers);
    
        //for (int i = 0; i < asuWithNeighbours.size(); i++)
        //{
        //    atomicNumbers.push_back(atomsInAsymmetricUnitAtomicNumbers[asuWithNeighbours[i].first]);
        //    symmetryOperation.set(asuWithNeighbours[i].second);
        //    //cartesianXyz, fractionalXyz
        //    fractionalXyz = crystal.atoms[asuWithNeighbours[i].first].coordinates;
        //    symmetryOperation.apply(fractionalXyz, fractionalTransformedXyz);
        //    crystal.unitCell.fractionalToCartesian(fractionalTransformedXyz, cartesianXyz);
        //    positions.push_back(cartesianXyz);
        //    labels.push_back(crystal.atoms[asuWithNeighbours[i].first].label);
        //    if (asuWithNeighbours[i].second != string("X,Y,Z"))
        //        labels.back() += string("(") + asuWithNeighbours[i].second + string(")");
        //}
    
    
    }
    
    void assymetricUnitWithNeighbours(
        const Crystal &crystal,
        std::vector<int> &atomicNumbers,
        std::vector<Vector3d> &positions,
        std::vector<std::string> &labels,
        int neighbourRange,
        double threshold,
        std::vector<int> &shellSizes)
    {

        Vector3d cartesianXyz, fractionalXyz, fractionalTransformedXyz;

        atomicNumbers.clear();
        positions.clear();
        labels.clear();
        SpaceGroupOperation symmetryOperation;

        vector< pair<int, string> > asuWithNeighbours;
        vector<int> atomsInAsymmetricUnitAtomicNumbers;

        structural_properties::assymetricUnitWithNeighbours(crystal, asuWithNeighbours, neighbourRange, threshold, shellSizes);
        crystal_structure_utilities::atomicNumbers(crystal, atomsInAsymmetricUnitAtomicNumbers);

        for (int i = 0; i < asuWithNeighbours.size(); i++)
        {
            atomicNumbers.push_back(atomsInAsymmetricUnitAtomicNumbers[asuWithNeighbours[i].first]);
            symmetryOperation.set(asuWithNeighbours[i].second);
            //cartesianXyz, fractionalXyz
            fractionalXyz = crystal.atoms[asuWithNeighbours[i].first].coordinates;
            symmetryOperation.apply(fractionalXyz, fractionalTransformedXyz);
            crystal.unitCell.fractionalToCartesian(fractionalTransformedXyz, cartesianXyz);
            positions.push_back(cartesianXyz);
            labels.push_back(crystal.atoms[asuWithNeighbours[i].first].label);
            if (asuWithNeighbours[i].second != string("X,Y,Z"))
                labels.back() += string("(") + asuWithNeighbours[i].second + string(")");
        }

    }

    void assymetricUnitWithNeighbours(
        const Crystal &crystal,
        std::vector< std::pair<int, std::string> > &asuWithNeighbours,
        std::vector<int> &atomicNumbers,
        std::vector<Vector3d> &positions,
        std::vector<std::string> &labels,
        int neighbourRange,
        double threshold,
        std::vector<int> &shellSizes)
    {

        Vector3d cartesianXyz, fractionalXyz, fractionalTransformedXyz;

        atomicNumbers.clear();
        positions.clear();
        labels.clear();
        SpaceGroupOperation symmetryOperation;

        asuWithNeighbours.clear();
        vector<int> atomsInAsymmetricUnitAtomicNumbers;

        structural_properties::assymetricUnitWithNeighbours(crystal, asuWithNeighbours, neighbourRange, threshold, shellSizes);
        crystal_structure_utilities::atomicNumbers(crystal, atomsInAsymmetricUnitAtomicNumbers);

        for (int i = 0; i < asuWithNeighbours.size(); i++)
        {
            atomicNumbers.push_back(atomsInAsymmetricUnitAtomicNumbers[asuWithNeighbours[i].first]);
            symmetryOperation.set(asuWithNeighbours[i].second);
            //cartesianXyz, fractionalXyz
            fractionalXyz = crystal.atoms[asuWithNeighbours[i].first].coordinates;
            symmetryOperation.apply(fractionalXyz, fractionalTransformedXyz);
            crystal.unitCell.fractionalToCartesian(fractionalTransformedXyz, cartesianXyz);
            positions.push_back(cartesianXyz);
            labels.push_back(crystal.atoms[asuWithNeighbours[i].first].label);
            if (asuWithNeighbours[i].second != string("X,Y,Z"))
                labels.back() += string("(") + asuWithNeighbours[i].second + string(")");
        }

    }


    void calcUnitCellConnectivity(
        const UnitCellContent &uc,
        std::vector<std::vector<UnitCellContent::AtomID> > &connectivity,
        double treshold)
    {
        int atomIndex, nAtomsInUnitCell, i;
        Vector3d aCartesian, bCartesian, cCartesian, latticeVectorCartesian, fractional, position, diff;
        vector<Vector3d> cartesianCoordinates;
        vector<int> atomicNumbers;
        const UnitCell &unitCell = uc.getCrystal().unitCell;
        int a, b, c;
        CovalentRadiousBondDetector bondDetector;
        double distance;
        bondDetector.setThreshold(treshold);
        connectivity.clear();

        nAtomsInUnitCell = uc.nAtoms();

        atomicNumbers.resize(nAtomsInUnitCell);

        vector<int> z_asymm_unit;
        crystal_structure_utilities::atomicNumbers(uc.getCrystal(), z_asymm_unit);

        for (int i = 0; i < nAtomsInUnitCell; i++)
            //atomicNumbers[i]  = basic_chemistry_utilities::atomicNumberFromLabel(uc.getAtom(i).label);
            atomicNumbers[i] = z_asymm_unit[uc.indexOfSymmetryEquivalentAtomInCrystal(i)];

        connectivity.resize(nAtomsInUnitCell);
        cartesianCoordinates.resize(nAtomsInUnitCell);

        for (atomIndex = 0; atomIndex < nAtomsInUnitCell; atomIndex++)
            unitCell.fractionalToCartesian(uc.getAtom(atomIndex).coordinates, cartesianCoordinates[atomIndex]);

        //------
        //ofstream out("unit_cell.xyz");

        //out << cartesianCoordinates.size() << endl;
        //for (int i = 0; i < cartesianCoordinates.size(); i++)
        //    out << periodic_table::symbol(atomicNumbers[i]) << " " << cartesianCoordinates[i][0]
        //        << " " << cartesianCoordinates[i][1] << " " << cartesianCoordinates[i][2] << endl;
 
        //out.close();

        //------

        unitCell.fractionalToCartesian(Vector3d(1, 0, 0), aCartesian);
        unitCell.fractionalToCartesian(Vector3d(0, 1, 0), bCartesian);
        unitCell.fractionalToCartesian(Vector3d(0, 0, 1), cCartesian);

        connectivity.resize(nAtomsInUnitCell);

        for (atomIndex = 0; atomIndex < nAtomsInUnitCell; atomIndex++)
        {
            // iterate over atoms in 3x3x3 unit cells box

            for (a = -1; a<2; a++)
                for (b = -1; b<2; b++)
                    for (c = -1; c<2; c++)
                        for (i = 0; i<nAtomsInUnitCell; i++)
                        {
                            if (atomIndex == i)
                                if (a == 0 && b == 0 && c == 0)
                                    continue;

                            position = cartesianCoordinates[i] + double(a)*aCartesian + double(b)*bCartesian + double(c)*cCartesian;
                            diff = cartesianCoordinates[atomIndex] - position;
                            distance = sqrt(diff*diff);
                            if (bondDetector.areBonded(atomicNumbers[atomIndex], atomicNumbers[i], distance))
                                connectivity[atomIndex].push_back(UnitCellContent::AtomID(i, Vector3i(a, b, c)));
                        }
        }

    }

    void calcUnitCellConnectivity2(
        const UnitCellContent &uc,
        std::vector<std::vector<UnitCellContent::AtomID> > &_connectivity,
        double treshold)
    {
        int atomIndex, nAtomsInUnitCell;
        Vector3d aCartesian, bCartesian, cCartesian, latticeVectorCartesian, fractional, position, diff;
        vector<Vector3d> cartesianCoordinates;
        vector<int> atomicNumbers;
        vector<vector<int> > connectivity;
        const UnitCell &unitCell = uc.getCrystal().unitCell;
        int a, b, c;
        CovalentRadiousBondDetector bondDetector;

        bondDetector.setThreshold(treshold);
        _connectivity.clear();

        nAtomsInUnitCell = uc.nAtoms();

        atomicNumbers.resize(nAtomsInUnitCell);

        //for (int i = 0; i < nAtomsInUnitCell; i++)
          //  atomicNumbers[i] = basic_chemistry_utilities::atomicNumberFromLabel(uc.getAtom(i).label);


        //--
        vector<int> z_asymm_unit;
        crystal_structure_utilities::atomicNumbers(uc.getCrystal(), z_asymm_unit);

        for (int i = 0; i < nAtomsInUnitCell; i++)
            //atomicNumbers[i]  = basic_chemistry_utilities::atomicNumberFromLabel(uc.getAtom(i).label);
            atomicNumbers[i] = z_asymm_unit[uc.indexOfSymmetryEquivalentAtomInCrystal(i)];

        //--

        
        cartesianCoordinates.resize(nAtomsInUnitCell);

        connectivity.resize(nAtomsInUnitCell);

        for (atomIndex = 0; atomIndex < nAtomsInUnitCell; atomIndex++)
            unitCell.fractionalToCartesian(uc.getAtom(atomIndex).coordinates, cartesianCoordinates[atomIndex]);

        unitCell.fractionalToCartesian(Vector3d(1, 0, 0), aCartesian);
        unitCell.fractionalToCartesian(Vector3d(0, 1, 0), bCartesian);
        unitCell.fractionalToCartesian(Vector3d(0, 0, 1), cCartesian);

        

        // 
        vector<pair<int, Vector3i> > ids(nAtomsInUnitCell);

        // find bounding box
        Vector3d startCorner, endCorner;
        double boxSize = 6; // in Angstroms
        startCorner = endCorner = cartesianCoordinates[0];
        for (int atomIdx = 0; atomIdx < nAtomsInUnitCell; atomIdx++)
        {
            for (int i = 0; i < 3; i++)
            {
                startCorner[i] = min(startCorner[i], cartesianCoordinates[atomIdx][i]);
                endCorner[i] = max(endCorner[i], cartesianCoordinates[atomIdx][i]);
            }
            ids[atomIdx] = { atomIdx, Vector3i{0,0,0} };
        }

        // adjsut startCorner, endCorner, set nBoxes

        startCorner -= Vector3d(0.1, 0.1, 0.1);
        endCorner += Vector3d(0.1, 0.1, 0.1);

        Vector3i nBoxes;
        for (int i = 0; i < 3; i++)
            nBoxes[i] = int((endCorner[i] - startCorner[i]) / boxSize) + 3;

        startCorner -= Vector3d(boxSize, boxSize, boxSize);
        endCorner = startCorner + boxSize * Vector3d(nBoxes);

        //

        vector<vector<vector<vector<int> > > > nonCoreAtomsInBox(nBoxes[0], vector<vector<vector<int> > >(nBoxes[1], vector<vector<int> >(nBoxes[2])));
        vector<vector<vector<vector<int> > > > coreAtomsInBox(nBoxes[0], vector<vector<vector<int> > >(nBoxes[1], vector<vector<int> >(nBoxes[2])));


        // find non core atoms and assign to cells
        Vector3d unitCellPosition;
        
        int i, j, k;
        for (a = -1; a < 2; a++)
            for (b = -1; b < 2; b++)
                for (c = -1; c < 2; c++)
                {
                    if (a == 0 && b == 0 && c == 0)
                        continue;

                    unitCellPosition = double(a)*aCartesian + double(b)*bCartesian + double(c)*cCartesian;
                    for (atomIndex = 0; atomIndex < nAtomsInUnitCell; atomIndex++)
                    {
                        position = cartesianCoordinates[atomIndex] + unitCellPosition;
                        i = int((position[0] - startCorner[0]) / boxSize);
                        j = int((position[1] - startCorner[1]) / boxSize);
                        k = int((position[2] - startCorner[2]) / boxSize);
                        if (i >= 0 && i < nBoxes[0] && j >= 0 && j < nBoxes[1] && k >= 0 && k < nBoxes[2])
                        {
                            ids.push_back({ atomIndex,{a,b,c} });
                            nonCoreAtomsInBox[i][j][k].push_back(cartesianCoordinates.size());
                            cartesianCoordinates.push_back(position);
                        }
                    }
                }
        // find cells for core atoms

        for (atomIndex = 0; atomIndex < nAtomsInUnitCell; atomIndex++)
        {
            i = int((cartesianCoordinates[atomIndex][0] - startCorner[0]) / boxSize);
            j = int((cartesianCoordinates[atomIndex][1] - startCorner[1]) / boxSize);
            k = int((cartesianCoordinates[atomIndex][2] - startCorner[2]) / boxSize);
            coreAtomsInBox[i][j][k].push_back(atomIndex);
        }

        // 

        // core core connectivity
        double d;
        vector<Vector3i> neighbourCell; // neighbours
        for (i = -1; i < 2; i++)
            for (j = -1; j < 2; j++)
                for (k = -1; k < 2; k++)
                    if (!(i == 0 && j == 0 && k == 0))
                        neighbourCell.push_back({ i,j,k });

        for (i = 1; i < nBoxes[0] - 1; i++)
            for (j = 1; j < nBoxes[1] - 1; j++)
                for (k = 1; k < nBoxes[2] - 1; k++)
                {
                    vector<int> &idx = coreAtomsInBox[i][j][k];

                    // inter cell connectivity
                    for (int idx1 : idx)
                        for (int idx2 : idx)
                        {
                            diff = cartesianCoordinates[idx1] - cartesianCoordinates[idx2];
                            d = sqrt(diff*diff);
                            if (bondDetector.areBonded(atomicNumbers[idx1], atomicNumbers[idx2], d))
                                if (idx1 != idx2)
                                    connectivity[idx1].push_back(idx2);
                        }
                    // cell - neighbour cell connectivity
                    for (auto &nb : neighbourCell)
                    {
                        vector<int> &idxNeighbour = coreAtomsInBox[i + nb[0]][j + nb[1]][k + nb[2]];
                        for (int idx1 : idx)
                            for (int idx2 : idxNeighbour)
                            {
                                diff = cartesianCoordinates[idx1] - cartesianCoordinates[idx2];
                                d = sqrt(diff*diff);
                                if (bondDetector.areBonded(atomicNumbers[idx1], atomicNumbers[idx2], d))
                                    if (idx1 != idx2)
                                        connectivity[idx1].push_back(idx2);
                            }
                    }
                }

    // core non-core connectivity
        neighbourCell.push_back({ 0,0,0 });

        for (i = 1; i < nBoxes[0] - 1; i++)
            for (j = 1; j < nBoxes[1] - 1; j++)
                for (k = 1; k < nBoxes[2] - 1; k++)
                {
                    vector<int> &idx = coreAtomsInBox[i][j][k];

                    // cell - neighbour cell connectivity
                    for (auto &nb : neighbourCell)
                    {
                        vector<int> &idxNeighbour = nonCoreAtomsInBox[i + nb[0]][j + nb[1]][k + nb[2]];

                        for (int idx1 : idx)
                            for (int idx2 : idxNeighbour)
                            {
                                diff = cartesianCoordinates[idx1] - cartesianCoordinates[idx2];
                                d = sqrt(diff*diff);
                                if (bondDetector.areBonded(atomicNumbers[idx1], atomicNumbers[idx2], d))
                                        connectivity[idx1].push_back(idx2);
                            }
                    }
                }

        //
        _connectivity.resize(nAtomsInUnitCell);
        for(atomIndex=0; atomIndex<nAtomsInUnitCell; atomIndex++)
            for (auto &neighbour : connectivity[atomIndex])
                _connectivity[atomIndex].push_back(UnitCellContent::AtomID(ids[neighbour].first, ids[neighbour].second));

        //for (atomIndex = 0; atomIndex < nAtomsInUnitCell; atomIndex++)
        //{
        //    // iterate over atoms in 3x3x3 unit cells box

        //    for (a = -1; a < 2; a++)
        //        for (b = -1; b < 2; b++)
        //            for (c = -1; c < 2; c++)
        //                for (i = 0; i < nAtomsInUnitCell; i++)
        //                {
        //                    if (atomIndex == i)
        //                        if (a == 0 && b == 0 && c == 0)
        //                            continue;

        //                    position = cartesianCoordinates[i] + double(a)*aCartesian + double(b)*bCartesian + double(c)*cCartesian;
        //                    diff = cartesianCoordinates[atomIndex] - position;
        //                    distance = sqrt(diff*diff);
        //                    if (bondDetector.areBonded(atomicNumbers[atomIndex], atomicNumbers[i], distance))
        //                        connectivity[atomIndex].push_back(UnitCellContent::AtomID(i, Vector3i(a, b, c)));
        //                }
        //}

    }

    void groupSymmetryRelatedMolecules(
        const UnitCellContent& uc,
        std::vector<std::vector<UnitCellContent::AtomID> >& molecules,
        //[group][molecule in group] .first-index .second-symmetry operation
        std::vector<std::vector<std::pair<int, SpaceGroupOperation> > >& moleculeGroups)
    {
        moleculeGroups.clear();

        int moleculeIdx, groupMolIdx, nMolecules = molecules.size();     
        string atom1Label, atom2Label, atom1SymmOp, atom2SymmOp;
        
        for (moleculeIdx = 0; moleculeIdx < nMolecules; moleculeIdx++)
        {
            bool groupFound = false;
            
            for (auto& group : moleculeGroups)
            {
                // if an atom with the same label as the first atom in molecule can be found in group[0]
                // i.e. if it contains symmetry equivalent atom
                // then the molecule should be added to the group
                
                groupMolIdx = group[0].first;
                for (auto atomInGroupRepresentative : molecules[groupMolIdx])
                {
                    uc.interpreteAtomID(molecules[moleculeIdx][0].atomIndex, atom1Label, atom1SymmOp);
                    uc.interpreteAtomID(atomInGroupRepresentative.atomIndex, atom2Label, atom2SymmOp);
                    

                    if (atom1Label == atom2Label)
                    {
                        groupFound = true;
                        auto transformingOperation = uc.getTransformingOperation(atomInGroupRepresentative, molecules[moleculeIdx][0]);
                        group.push_back({ moleculeIdx, transformingOperation });
                    }
                }
            }
            if (!groupFound)
            {
                moleculeGroups.resize(moleculeGroups.size() + 1);
                moleculeGroups.back().push_back({ moleculeIdx, SpaceGroupOperation(string("X,Y,Z")) });
            }

        }
        
    }


    void findIncludingMolecule(
        const UnitCellContent::AtomID& atomId,
        const std::vector<UnitCellContent::AtomID>& atomsToDisconnect,
        const std::vector<std::vector<UnitCellContent::AtomID> >& connectivity,
        std::vector<UnitCellContent::AtomID>& molecule,
        std::vector<std::pair<UnitCellContent::AtomID, UnitCellContent::AtomID> >& disconnectedBonds,
        const UnitCellContent& ucContent,
        int maxSize)
    {
        int nAtoms = connectivity.size();
        int atomIdx = atomId.atomIndex;
        vector<bool> alreadyChosen(nAtoms, false);
        vector<Vector3i> chosenAtomUnitCell(nAtoms);
        Vector3i neighbourUnitCell;
        vector<UnitCellContent::AtomID> lastAdded, nextLastAdded;
        // first - atom from molecule, second - disconnected atom
        vector<pair<UnitCellContent::AtomID, UnitCellContent::AtomID> > recreateBond;
        UnitCellContent::AtomID neighbour;
        int nNeighbours, neighbourCounter, currentAtom;

        molecule.clear();
        

        lastAdded.push_back(UnitCellContent::AtomID(atomIdx, atomId.unitCellPosition));
        molecule.push_back(lastAdded.back());
        //alreadyChosen[atomIdx] = true;
        //int shell = 1;
        while (!lastAdded.empty())
        {


            for (int i = 0; i < lastAdded.size(); i++)
            {
                currentAtom = lastAdded[i].atomIndex;
                nNeighbours = connectivity[currentAtom].size();

                for (neighbourCounter = 0; neighbourCounter < nNeighbours; neighbourCounter++)
                {
                    neighbour.atomIndex = connectivity[currentAtom][neighbourCounter].atomIndex;

                    neighbour.unitCellPosition = lastAdded[i].unitCellPosition + connectivity[currentAtom][neighbourCounter].unitCellPosition;


                    if (find(molecule.begin(), molecule.end(), neighbour) == molecule.end())
                        if (find(nextLastAdded.begin(), nextLastAdded.end(), neighbour) == nextLastAdded.end())
                        {
                            if (find(atomsToDisconnect.begin(), atomsToDisconnect.end(), neighbour) == atomsToDisconnect.end())
                                nextLastAdded.push_back(neighbour);
                            else
                                disconnectedBonds.push_back({lastAdded[i], neighbour});
                        }
                }
            }

            lastAdded = nextLastAdded;
            nextLastAdded.clear();

            if (maxSize <= molecule.size() + nextLastAdded.size())
            {
                on_error::throwException("too many atoms = " + to_string(molecule.size() + nextLastAdded.size()) +
                    " in molecule definition, unable to locate atom ", __FILE__, __LINE__);
            }
            else
            {
                molecule.insert(molecule.end(), lastAdded.begin(), lastAdded.end());
                //vector<Vector3d> positions;
                //vector<ChemicalElement> elements;
                //crystal_structure_utilities::convertToXyzAndElementList(ucContent.getCrystal(), ucContent, molecule, elements, positions);
                //xyz_io::writeXyz("mol_" + to_string(shell) + ".xyz", elements, positions);
                //ofstream out("mol_" + to_string(shell));
                //
                //for (auto atom : molecule)
                //{
                //    string label, symmOp;
                //    ucContent.interpreteAtomID(atom, label, symmOp);
                //    out << setw(12) << left << label << " " << symmOp << "\n";
                //}
                //out.close();


                //shell++;
            }
        }

    }


    void findIncludingMolecule(
        int atomIdx,
        const vector< vector< UnitCellContent::AtomID > > &connectivity,
        vector< UnitCellContent::AtomID > &molecule,
        vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > &networkBonds)
    {
        int nAtoms = connectivity.size();
        vector<bool> alreadyChosen(nAtoms, false);
        vector<Vector3i> chosenAtomUnitCell(nAtoms);
        Vector3i neighbourUnitCell;
        vector<UnitCellContent::AtomID> lastAdded, nextLastAdded;
        int nNeighbours, neighbourCounter, neighbourAtom, currentAtom;

        molecule.clear();
        networkBonds.clear();

        lastAdded.push_back(UnitCellContent::AtomID(atomIdx, Vector3i()));
        molecule.push_back(UnitCellContent::AtomID(atomIdx, Vector3i()));
        alreadyChosen[atomIdx] = true;

        while (!lastAdded.empty())
        {
            

            for (int i = 0; i < lastAdded.size(); i++)
            {
                currentAtom = lastAdded[i].atomIndex;
                nNeighbours = connectivity[currentAtom].size();
                
                for (neighbourCounter = 0; neighbourCounter < nNeighbours; neighbourCounter++)
                {
                    neighbourAtom = connectivity[currentAtom][neighbourCounter].atomIndex;

                    neighbourUnitCell = lastAdded[i].unitCellPosition + connectivity[currentAtom][neighbourCounter].unitCellPosition;

                    if (alreadyChosen[neighbourAtom])
                    {
                        if (neighbourUnitCell != chosenAtomUnitCell[neighbourAtom])
                            networkBonds.push_back(make_pair(lastAdded[i], UnitCellContent::AtomID(neighbourAtom, neighbourUnitCell)));
                    }
                    else
                    {
                        alreadyChosen[neighbourAtom] = true;
                        chosenAtomUnitCell[neighbourAtom] = neighbourUnitCell;
                        nextLastAdded.push_back(UnitCellContent::AtomID(neighbourAtom, neighbourUnitCell));
                    }

                }
            }

            lastAdded = nextLastAdded;
            nextLastAdded.clear();
            molecule.insert(molecule.end(), lastAdded.begin(), lastAdded.end());
            
        }

    }

    void findIncludingMolecule(
        const std::vector<int> &atoms,
        const UnitCellContent &uc,
        std::vector<std::vector<UnitCellContent::AtomID> > &molecules,
        std::vector< std::vector< std::pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > &networkBonds,
        double threshold)
    {
        vector<vector<UnitCellContent::AtomID> > connectivity;
        std::vector<UnitCellContent::AtomID > molecule;
        set<int> notInMolecule;
        vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > moleculeNetworkBonds;

        molecules.clear();
        networkBonds.clear();

        for (int i = 0, n = atoms.size(); i < n; i++)
            notInMolecule.insert(atoms[i]);

        calcUnitCellConnectivity(uc, connectivity, threshold);

        while (!notInMolecule.empty())
        {
            findIncludingMolecule(*notInMolecule.begin(), connectivity, molecule, moleculeNetworkBonds);

            for (int i = 0, n = molecule.size(); i < n; i++)
                notInMolecule.erase(molecule[i].atomIndex);

            molecules.push_back(molecule);
            networkBonds.push_back(moleculeNetworkBonds);

            molecule.clear();
            moleculeNetworkBonds.clear();
        }

    }

    void splitUnitCellIntoMolecules(
        const UnitCellContent &uc,
        std::vector<std::vector<UnitCellContent::AtomID > > &molecules,
        std::vector< std::vector< std::pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > &networkBonds,
        double threshold,
		bool eachAtomPresentOnlyOnce)
    {
		vector<int> atoms(uc.nAtoms());
		 
        for (int i = 0, n = uc.nAtoms(); i < n; i++)
            atoms[i] = i;

        findIncludingMolecule(atoms, uc, molecules, networkBonds, threshold);

		if (eachAtomPresentOnlyOnce)
			return;

		vector<vector<UnitCellContent::AtomID > > moreMolecules;
		set<Vector3i> uniqueUnitCells;
		Vector3i vectorZero(0, 0, 0);
		for (auto& molecule : molecules)
		{
			uniqueUnitCells.clear();
			for (auto& atom : molecule)
				uniqueUnitCells.insert(atom.unitCellPosition);
			for(auto &cell: uniqueUnitCells)
				if (cell != vectorZero)
				{
					moreMolecules.resize(moreMolecules.size()+1);
					for (auto& atom : molecule)
						moreMolecules.back().push_back({ atom.atomIndex, atom.unitCellPosition - cell });
				}
		}
		molecules.insert(molecules.end(), moreMolecules.begin(), moreMolecules.end());
    }

    void splitIntoMolecules(
        const std::vector<int> &atomicNumbers,
        const std::vector<Vector3d> &positionsInAngstroms,
        std::vector<std::vector<int> > &molecules,
        double threshold)
    {
        GenericConnectivityAlgorithm<CovalentRadiousBondDetector> connectivityAlgorithm;
        connectivityAlgorithm.set(to_string(threshold));
        vector < vector<int> > connectivity;
        //MolecularDisorder disorder;
        
        connectivityAlgorithm.calculateConnectivity(positionsInAngstroms, atomicNumbers, connectivity);
        graph_algorithms::split(connectivity, molecules);
        //set<int> notProcessedAtoms;
        //vector<vector<int> > neighbors;
        //molecules.clear();
        //int atomIdx, nAtoms = atomicNumbers.size();
        //for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        //    notProcessedAtoms.insert(atomIdx);


        //while (!notProcessedAtoms.empty())
        //{
        //    graph_algorithms::breadth_first_search(connectivity, *notProcessedAtoms.begin(), neighbors);
        //    molecules.resize(molecules.size() + 1);
        //    for (auto &shell : neighbors)
        //    {
        //        molecules.back().insert(molecules.back().end(), shell.begin(), shell.end());
        //        for (auto &atom : shell)
        //            notProcessedAtoms.erase(atom);
        //    }
        //}
    }



    void splitAsymmetricUnitIntoMolecules(
        const UnitCellContent &uc,
        std::vector<std::vector<UnitCellContent::AtomID> > &molecules,
        std::vector< std::vector< std::pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > &networkBonds,
        double threshold)
    {
        int nAtomsInAsymmetricUnit = uc.getCrystal().atoms.size();
        vector<int> atoms(nAtomsInAsymmetricUnit);

        for (int i = 0, n = nAtomsInAsymmetricUnit; i < n; i++)
            atoms[i] = i;

        findIncludingMolecule(atoms, uc, molecules, networkBonds, threshold);
    }


    void graphToNthNeighbour(
        const UnitCellContent &unitCellContent,
        const vector<UnitCellContent::AtomID> &startingSet,
        vector<UnitCellContent::AtomID> &graph,
        int n,
        double threshold)
    {
        vector<int> shellSizes;
        graphToNthNeighbour(unitCellContent, startingSet, graph, n, threshold, shellSizes);



        /*
        std::vector<std::vector<UnitCellContent::AtomID> > unitCellConnectivity;
        calcUnitCellConnectivity(unitCellContent, unitCellConnectivity, threshold);
                        
        //ofstream out("uc_connect");
        //for (auto &atomConn : unitCellConnectivity)
        //{
        //    
        //    for (auto &neighb : atomConn)
        //    {
        //        out << neighb.atomIndex << ",(" << neighb.unitCellPosition[0] << ", " << neighb.unitCellPosition[1]
        //            << ", " << neighb.unitCellPosition[2] << ")  ";
        //    }
        //    out << endl;
        //}

        //out.flush();
        //out.close();

        vector<UnitCellContent::AtomID> newShell;
        vector<vector<UnitCellContent::AtomID> > graphShells;
                
        int i, j; 
        
        graphShells.resize(1);
        for (i = 0; i < startingSet.size(); i++)
            graphShells[0].push_back(startingSet[i]);

        UnitCellContent::AtomID neighbourId;

        for (int step = 1; step <= n; step++)
        {
            newShell.clear();

            auto &outerShell = graphShells.back();

            for (auto node: outerShell)
            {

                for (auto neighbour : unitCellConnectivity[node.atomIndex])
                {
                    neighbourId = neighbour;
                    neighbourId.unitCellPosition += node.unitCellPosition;
                    // check if neighbbour already selected, if not add to new shell
                    bool alreadySelected = false;
                    for (auto shell : graphShells)
                        if (find(shell.begin(), shell.end(), neighbourId) != shell.end())
                        {
                            alreadySelected = true;
                            break;
                        }

                    if (!alreadySelected)
                        if(find(newShell.begin(), newShell.end(), neighbourId) == newShell.end())
                            newShell.push_back(neighbourId);
                }
            }

            if (newShell.empty())
                break;
            else
                graphShells.push_back(newShell);
        }

        // generate graph

        graph.clear();
        for (auto &shell : graphShells)
            graph.insert(graph.end(), shell.begin(), shell.end());
            */
    }

    void graphToNthNeighbour(
        const UnitCellContent &unitCellContent,
        const std::vector<UnitCellContent::AtomID> &startingSet,
        std::vector<UnitCellContent::AtomID> &graph,
        int n,
        double threshold,
        std::vector<int> &shellSizes)
    {
        std::vector<std::vector<UnitCellContent::AtomID> > unitCellConnectivity;
        calcUnitCellConnectivity(unitCellContent, unitCellConnectivity, threshold);

        //ofstream out("uc_connect");
        //for (auto &atomConn : unitCellConnectivity)
        //{
        //    
        //    for (auto &neighb : atomConn)
        //    {
        //        out << neighb.atomIndex << ",(" << neighb.unitCellPosition[0] << ", " << neighb.unitCellPosition[1]
        //            << ", " << neighb.unitCellPosition[2] << ")  ";
        //    }
        //    out << endl;
        //}

        //out.flush();
        //out.close();

        vector<UnitCellContent::AtomID> newShell;
        vector<vector<UnitCellContent::AtomID> > graphShells;

        int i;

        graphShells.resize(1);
        for (i = 0; i < startingSet.size(); i++)
            graphShells[0].push_back(startingSet[i]);

        UnitCellContent::AtomID neighbourId;

        for (int step = 1; step <= n; step++)
        {
            newShell.clear();

            auto &outerShell = graphShells.back();

            for (auto node : outerShell)
            {

                for (auto neighbour : unitCellConnectivity[node.atomIndex])
                {
                    neighbourId = neighbour;
                    neighbourId.unitCellPosition += node.unitCellPosition;
                    // check if neighbbour already selected, if not add to new shell
                    bool alreadySelected = false;
                    for (auto shell : graphShells)
                        if (find(shell.begin(), shell.end(), neighbourId) != shell.end())
                        {
                            alreadySelected = true;
                            break;
                        }

                    if (!alreadySelected)
                        if (find(newShell.begin(), newShell.end(), neighbourId) == newShell.end())
                            newShell.push_back(neighbourId);
                }
            }

            if (newShell.empty())
                break;
            else
                graphShells.push_back(newShell);
        }

        // generate graph

        graph.clear();
        shellSizes.clear();
        for (auto &shell : graphShells)
        {
            graph.insert(graph.end(), shell.begin(), shell.end());
            shellSizes.push_back(shell.size());
        }

    }


    void calculateConnectivity(
        const std::vector<std::vector<UnitCellContent::AtomID> > &unitCellConnectivity,
        const std::vector<UnitCellContent::AtomID> &atoms,
        std::vector<std::vector<int> > &connectivity)
    {
        int i, j, n = atoms.size();
        connectivity.clear();
        connectivity.resize(n);

        for (i = 0; i < n; i++)
            for ( j = 0; j < i; j++)
                for (auto neighbour : unitCellConnectivity[atoms[i].atomIndex])
                    if (atoms[j].atomIndex == neighbour.atomIndex)
                        if (atoms[j].unitCellPosition - atoms[i].unitCellPosition == neighbour.unitCellPosition)
                        {
                            connectivity[i].push_back(j);
                            connectivity[j].push_back(i);
                        }

    }

    void makeCluster(
        const UnitCellContent& ucContent,
        const std::vector<UnitCellContent::AtomID>& centralMolecule,
        const std::vector<std::vector<UnitCellContent::AtomID> >& molecules,
        std::vector<UnitCellContent::AtomID>& clusterAtoms,
        double threshold, bool vdwThreshold)
    {
        // [molecule ndex][atom in molecule index]
        vector<vector<Vector3d> > coordinates;
        vector<Vector3d> centralMoleculeAtomicCoordinates;
        vector<Vector3d> latticeVec(3);
        Vector3d fractionalCoordinates, cartesianCoordinates;;

        auto& crystal = ucContent.getCrystal();
        crystal.unitCell.fractionalToCartesian({ 1,0,0 }, latticeVec[0]);
        crystal.unitCell.fractionalToCartesian({ 0,1,0 }, latticeVec[1]);
        crystal.unitCell.fractionalToCartesian({ 0,0,1 }, latticeVec[2]);

        //------- central molecule --------------
        // set coordinates and atomic numbers

        vector<int> centralMolZ;
        vector<int> atomicNumbers;
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);
        for (auto const& atom : centralMolecule)
        {
            fractionalCoordinates = ucContent.getAtom(atom.atomIndex).coordinates + Vector3d(atom.unitCellPosition);
            crystal.unitCell.fractionalToCartesian(fractionalCoordinates, cartesianCoordinates);
            centralMoleculeAtomicCoordinates.push_back(cartesianCoordinates);
            centralMolZ.push_back(atomicNumbers[ucContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex)]);
        }

        //------- unit cell molecules --------
        // set coordinates and atomic numbers
        vector<vector<int> > moleculesZ;
        int nMolecules = molecules.size();
        coordinates.resize(nMolecules);
        moleculesZ.resize(nMolecules);
        for (int i = 0; i < nMolecules; i++)
            for (auto const& atom : molecules[i])
            {
                fractionalCoordinates = ucContent.getAtom(atom.atomIndex).coordinates + Vector3d(atom.unitCellPosition);
                crystal.unitCell.fractionalToCartesian(fractionalCoordinates, cartesianCoordinates);
                coordinates[i].push_back(cartesianCoordinates);
                moleculesZ[i].push_back(atomicNumbers[ucContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex)]);
            }

        //-----------------


        vector<vector<Vector3i> > moleculesIntersectingUnitCell(molecules.size());


        set<Vector3i> uniqueUnitCells;
        Vector3i vectorZero(0, 0, 0);
        for (int molIdx = 0; molIdx < molecules.size(); molIdx++)
        {
            auto& molecule = molecules[molIdx];
            uniqueUnitCells.clear();
            for (auto& atom : molecule)
                uniqueUnitCells.insert(atom.unitCellPosition);
            for (auto& cell : uniqueUnitCells)
                moleculesIntersectingUnitCell[molIdx].push_back(-cell);

        }

        set<pair<int, Vector3i> > shellMolecules, cluster, clusterExtension;

        bool clusterWasExtended = true;
        int shell = 0;

        while (clusterWasExtended)
        {
            calculateShellMolecules(shell, moleculesIntersectingUnitCell, set<pair<int, Vector3i> >(), shellMolecules);
            findClusterExtendingMolecules(centralMoleculeAtomicCoordinates, centralMolZ, coordinates, moleculesZ, latticeVec,
                shellMolecules, clusterExtension, threshold, vdwThreshold);
            clusterWasExtended = !clusterExtension.empty();
            shell++;
            cluster.insert(clusterExtension.begin(), clusterExtension.end());
        }

        for (auto& molecule : cluster)
        {
            auto const& mol = molecules[molecule.first];

            for (int i = 0; i < mol.size(); i++)
                clusterAtoms.push_back({ mol[i].atomIndex, mol[i].unitCellPosition + molecule.second });
        }

    }



	void makeCluster(
		const UnitCellContent& ucContent,
		const std::vector<UnitCellContent::AtomID>& centralMolecule,
		std::vector<UnitCellContent::AtomID>& clusterAtoms,
		double threshold, 
        bool vdwThreshold)
	{
		vector<vector<UnitCellContent::AtomID> > molecules;
		vector<UnitCellContent::AtomID> uc;
		vector< vector< pair< UnitCellContent::AtomID, UnitCellContent::AtomID > > > networkBonds;

		discamb::structural_properties::splitUnitCellIntoMolecules(ucContent, molecules, networkBonds, 0.4);

		//
		// [molecule ndex][atom in molecule index]
		vector<vector<Vector3d> > coordinates;
		vector<Vector3d> centralMoleculeAtomicCoordinates;
		vector<Vector3d> latticeVec(3);
		Vector3d fractionalCoordinates, cartesianCoordinates;;

		auto &crystal = ucContent.getCrystal();
		crystal.unitCell.fractionalToCartesian({ 1,0,0 }, latticeVec[0]);
		crystal.unitCell.fractionalToCartesian({ 0,1,0 }, latticeVec[1]);
		crystal.unitCell.fractionalToCartesian({ 0,0,1 }, latticeVec[2]);

        //------- central molecule --------------
        // set coordinates and atomic numbers

        vector<int> centralMolZ;
        vector<int> atomicNumbers;
        crystal_structure_utilities::atomicNumbers(crystal, atomicNumbers);
		for (auto const& atom : centralMolecule)
		{
			fractionalCoordinates = ucContent.getAtom(atom.atomIndex).coordinates + Vector3d(atom.unitCellPosition);
			crystal.unitCell.fractionalToCartesian(fractionalCoordinates, cartesianCoordinates);
			centralMoleculeAtomicCoordinates.push_back(cartesianCoordinates);
            centralMolZ.push_back(atomicNumbers[ucContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex)]);
		}

        //------- unit cell molecules --------
        // set coordinates and atomic numbers
        vector<vector<int> > moleculesZ;
		int nMolecules = molecules.size();
		coordinates.resize(nMolecules);
        moleculesZ.resize(nMolecules);
		for (int i = 0; i < nMolecules; i++)
			for (auto const& atom : molecules[i])
			{
				fractionalCoordinates = ucContent.getAtom(atom.atomIndex).coordinates + Vector3d(atom.unitCellPosition);
				crystal.unitCell.fractionalToCartesian(fractionalCoordinates, cartesianCoordinates);
				coordinates[i].push_back(cartesianCoordinates);
                moleculesZ[i].push_back(atomicNumbers[ucContent.indexOfSymmetryEquivalentAtomInCrystal(atom.atomIndex)]);
			}

		//-----------------


		vector<vector<Vector3i> > moleculesIntersectingUnitCell(molecules.size());


		set<Vector3i> uniqueUnitCells;
		Vector3i vectorZero(0, 0, 0);
		for (int molIdx = 0; molIdx < molecules.size(); molIdx++)
		{
			auto& molecule = molecules[molIdx];
			uniqueUnitCells.clear();
			for (auto& atom : molecule)
				uniqueUnitCells.insert(atom.unitCellPosition);
			for (auto& cell : uniqueUnitCells)
				moleculesIntersectingUnitCell[molIdx].push_back(-cell);

		}

		set<pair<int, Vector3i> > shellMolecules, cluster, clusterExtension;

		bool clusterWasExtended = true;
		int shell = 0;
        /*
                const vector<Vector3d>& centralMolAtomicCoords,
        const vector<int> &centralMolZ, // central molecule atomic numbers
        const vector<vector<Vector3d> >& moleculesAtomicCoords,
        const vector < vector <int> > &moleculesZ, // atomic numbers
        const vector<Vector3d>& latticeVec,
        const set<pair<int, Vector3i> >& shellMolecules,
        set<pair<int, Vector3i> >& shellMoleculesclusterExtension,
        double threshold,
        bool useVdW)
        */
        //cout << "calculate cluster atoms\n";
		while (clusterWasExtended)
		{
            //cout << "shell " << shell << "\n";
            //cout << "calculating shell molecules\n";
			calculateShellMolecules(shell, moleculesIntersectingUnitCell, set<pair<int, Vector3i> >(), shellMolecules);
            //cout << "finding extending molecules\n";
			findClusterExtendingMolecules(centralMoleculeAtomicCoordinates, centralMolZ, coordinates, moleculesZ, latticeVec,
				                          shellMolecules, clusterExtension, threshold, vdwThreshold);
            //cout << "done\n";
			clusterWasExtended = !clusterExtension.empty();
			shell++;
			cluster.insert(clusterExtension.begin(), clusterExtension.end());
		}

		for (auto& molecule : cluster)
		{
			auto const& mol = molecules[molecule.first];

			for (int i = 0; i < mol.size(); i++)
				clusterAtoms.push_back({ mol[i].atomIndex, mol[i].unitCellPosition + molecule.second });
		}

	}


    double interatomicDistance(
        const Crystal& crystal,
        int atom1,
        int atom2,
        const SpaceGroupOperation& s1,
        const SpaceGroupOperation& s2)
    {
        Vector3d r1_fractional, r2_fractional, r1_cartesian, r2_cartesian, r12;

        s1.apply(crystal.atoms[atom1].coordinates, r1_fractional);
        s2.apply(crystal.atoms[atom2].coordinates, r2_fractional);

        crystal.unitCell.fractionalToCartesian(r1_fractional, r1_cartesian);
        crystal.unitCell.fractionalToCartesian(r2_fractional, r2_cartesian);

        r12 = r1_cartesian - r2_cartesian;

        return sqrt(r12 * r12);
    }

    double interatomicDistance(
        const Crystal& crystal,
        int atom1,
        int atom2,
        const CrystalVarianceCovarianceMatrix& vcov,
        double& stadardDeviation,
        const SpaceGroupOperation& s1,
        const SpaceGroupOperation& s2)
    {
        Matrix3d m, m11, m12, m21, m22, p, rot1_frac, rot2_frac, rot1_cart, rot1_cart_t, rot2_cart, rot2_cart_t, frac2cart, cart2frac;
        vector<vector<double> > c11, c12, c21, c22;
        vector<vector<double> > mVcov(6,vector<double>(6));
        int i, j;

        bool standardDeviationAvailable = true;
        structural_parameters_convention::XyzCoordinateSystem  xyzConvention;
        structural_parameters_convention::AdpConvention adpConvention;

        vcov.getConvention(xyzConvention, adpConvention);
        if (xyzConvention == structural_parameters_convention::XyzCoordinateSystem::fractional)
        {
            CrystalVarianceCovarianceMatrix cart_vcov(vcov);
            cart_vcov.changeXyzConvention(structural_parameters_convention::XyzCoordinateSystem::cartesian, crystal.unitCell);
            if (!cart_vcov.get(atom1, atom1, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c11))
                standardDeviationAvailable = false;
            if (!cart_vcov.get(atom1, atom2, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c12))
                standardDeviationAvailable = false;
            if (!cart_vcov.get(atom2, atom1, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c21))
                standardDeviationAvailable = false;
            if (!cart_vcov.get(atom2, atom2, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c22))
                standardDeviationAvailable = false;
        }
        else
        {
            if (!vcov.get(atom1, atom1, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c11))
                standardDeviationAvailable = false;
            if (!vcov.get(atom1, atom2, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c12))
                standardDeviationAvailable = false;
            if (!vcov.get(atom2, atom1, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c21))
                standardDeviationAvailable = false;
            if (!vcov.get(atom2, atom2, CrystalVarianceCovarianceMatrix::variable_type::xyz, CrystalVarianceCovarianceMatrix::variable_type::xyz, c22))
                standardDeviationAvailable = false;
        }

        if (standardDeviationAvailable)
        {
            //&&&&&&&&&&&&&
            s1.getRotation(rot1_frac);
            s2.getRotation(rot2_frac);
            //&&&&&&&&&&&&&
            frac2cart = crystal.unitCell.getFractionalToCartesianMatrix();
            cart2frac = crystal.unitCell.getCartesianToFractionalMatrix();

            rot1_cart = frac2cart * rot1_frac * cart2frac;
            rot1_cart_t = rot1_cart;
            rot1_cart_t.transpose();

            rot2_cart = frac2cart * rot1_frac * cart2frac;
            rot2_cart_t = rot2_cart;
            rot2_cart_t.transpose();

            //&&&&&&&&&&&&&

            m.set(c11[0][0], c11[0][1], c11[0][2],
                c11[1][0], c11[1][1], c11[1][2],
                c11[2][0], c11[2][1], c11[2][2]);

            m11 = rot1_cart * m * rot1_cart_t;
            m.set(c21[0][0], c21[0][1], c21[0][2],
                c21[1][0], c21[1][1], c21[1][2],
                c21[2][0], c21[2][1], c21[2][2]);

            m21 = rot2_cart * m * rot1_cart_t;
            m.set(c12[0][0], c12[0][1], c12[0][2],
                c12[1][0], c12[1][1], c12[1][2],
                c12[2][0], c12[2][1], c12[2][2]);

            m12 = rot1_cart * m * rot2_cart_t;
            m.set(c22[0][0], c22[0][1], c22[0][2],
                c22[1][0], c22[1][1], c22[1][2],
                c22[2][0], c22[2][1], c22[2][2]);

            m22 = rot2_cart * m * rot2_cart_t;



            

            for (i = 0; i < 3; i++)
                for (j = 0; j < 3; j++)
                {
                    mVcov[i][j] = m11(i, j);
                    mVcov[i + 3][j] = m21(i, j);
                    mVcov[i][j + 3] = m12(i, j);
                    mVcov[i + 3][j + 3] = m22(i, j);
                }
        }
        ////////
         
        double r12_length;
        Vector3d r1_fractional, r2_fractional, r1_cartesian, r2_cartesian, r12;

        s1.apply(crystal.atoms[atom1].coordinates, r1_fractional);
        s2.apply(crystal.atoms[atom2].coordinates, r2_fractional);

        crystal.unitCell.fractionalToCartesian(r1_fractional, r1_cartesian);
        crystal.unitCell.fractionalToCartesian(r2_fractional, r2_cartesian);

        r12 = r1_cartesian - r2_cartesian;

        r12_length = sqrt(r12 * r12);

        ////
        if (standardDeviationAvailable)
        {
            vector<double> jacobian(6);
            for (i = 0; i < 3; i++)
            {
                jacobian[i] = (r1_cartesian[i] - r2_cartesian[i]) / r12_length;
                jacobian[i + 3] = -jacobian[i];
            }

            double variance = 0;

            for (i = 0; i < 6; i++)
                for (j = 0; j < 6; j++)
                    variance += mVcov[i][j] * jacobian[i] * jacobian[j];
            stadardDeviation = sqrt(variance);
        }
        else
            stadardDeviation = 0;
        return r12_length;
    }

    //------ PLANARITY

    void planarity(
        const std::vector<int> atomicNumbers,
        const std::vector<Vector3d>& positions,
        const std::vector<std::vector<int> >& connectivityMatrix,
        vector<Tribool>& _planarity,
        double threshold)
    {
        vector<double> distance_esd;
        planarity(atomicNumbers, positions, connectivityMatrix, _planarity, distance_esd, threshold);
    }

    void planarity(
        const std::vector<int> atomicNumbers,
        const std::vector<Vector3d>& positions,
        const std::vector<std::vector<int> >& connectivityMatrix,
        vector<Tribool>& _planarity,
        std::vector<double>& distance_esd,
        double threshold)
    {
        int atomIndex, nAtoms = atomicNumbers.size();
        _planarity.resize(nAtoms);
        distance_esd.resize(nAtoms);
        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
            _planarity[atomIndex] = atomPlanarity(atomIndex, atomicNumbers, positions, connectivityMatrix[atomIndex], distance_esd[atomIndex], threshold);
    }


    void planarity(
        const std::vector<int> atomicNumbers, const std::vector<Vector3d>& positions,
        const std::vector<int>& atoms,
        const std::vector<std::vector<int> >& connectivity_matrix,
        std::vector<Tribool>& _planarity,
        double threshold)
    {
        vector<double> distance_esd;
        planarity(atomicNumbers, positions, atoms, connectivity_matrix, _planarity, distance_esd, threshold);
    }

    void planarity(
        const vector<int> atomicNumbers,
        const vector<Vector3d>& positions,
        const vector<int>& atoms,
        const vector<vector<int> >& connectivity_matrix,
        vector<Tribool>& _planarity,
        vector<double>& distance_esd,
        double threshold)
    {
        int atomIndex, nAtoms = atoms.size();
        _planarity.resize(nAtoms);
        distance_esd.resize(nAtoms);
        for (atomIndex = 0; atomIndex < nAtoms; atomIndex++)
            _planarity[atomIndex] = atomPlanarity(atoms[atomIndex], atomicNumbers, positions, connectivity_matrix[atoms[atomIndex]], distance_esd[atomIndex], threshold);
    }


    Tribool atomPlanarity(
        int atomIndex,
        const vector<int> atomicNumbers,
        const vector<Vector3d>& _positions,
        const vector<int>& neighbors,
        double& distance_esd,
        double threshold,
        double valence_threshold)
    {
        distance_esd = 0.0;
        int nNeighbors = neighbors.size();
        if (nNeighbors < 3)
            return Tribool::Undefined;

        int neighborIndex, nAtoms;
        vector<Vector3d> positions;
        Matrix3<double> matrix;
        Vector3d average, position;
        Vector3d eigenVectors[3], planeNormal;
        vector<pair<double, int> > eigenValues(3);


        nAtoms = nNeighbors + 1;

        // collect positions of atoms from the 'plane'

        positions.resize(nAtoms);
        positions[0] = _positions[atomIndex];
        for (neighborIndex = 0; neighborIndex < nNeighbors; neighborIndex++)
            positions[neighborIndex + 1] = _positions[neighbors[neighborIndex]];


        //vector<double> valence_angles;
        //for (int i = 0; i < nNeighbors; i++)
        //    for (int j = i; j < nNeighbors; j++)
        //        valence_angles.push_back(geometry3d::angle(positions[i + 1], positions[0], positions[j + 1]);

        return planarity(positions, threshold, distance_esd);
    }

        Tribool  planarity(
        const std::vector<Vector3d>& positions,
        double threshold,
        double& distance_esd)
    {
        int i, j, atom, nAtoms = positions.size();
        if (nAtoms < 3)
            return Tribool::Undefined;
        if (nAtoms == 3)
            return Tribool::True;

        double distance, averageOntoPlaneNormalProjection;

        Matrix3<double> matrix;
        Vector3d average, position;
        Vector3d eigenVectors[3], planeNormal;
        vector<pair<double, int> > eigenValues(3);

        // clculate average position

        for (atom = 0; atom < nAtoms; atom++)
            average += positions[atom];
        average /= (double) nAtoms;

        // calculate 'displacement matrix'

        for (i = 0; i < 3; i++)
            for (j = 0; j < 3; j++)
                for (atom = 0; atom < nAtoms; atom++)
                    matrix(i, j) += (positions[atom](i) - average(i)) * (positions[atom](j) - average(j));

        // find smallest eigenvector of the matrix = normal to the 'plane'

        for (i = 0; i < 3; i++)
            eigenValues[i].second = i;
        algebra3d::eigensystemRealSymm(matrix, eigenVectors[0], eigenVectors[1], eigenVectors[2], eigenValues[0].first, eigenValues[1].first, eigenValues[2].first);

        planeNormal = eigenVectors[std::min_element(eigenValues.begin(), eigenValues.end())->second];

        // normalize

        planeNormal /= sqrt(planeNormal * planeNormal);

        // calculate sigma of atom distance to plane

        averageOntoPlaneNormalProjection = average * planeNormal;
        distance_esd = 0;

        for (atom = 0; atom < nAtoms; atom++)
        {
            distance = planeNormal * positions[atom] - averageOntoPlaneNormalProjection;
            distance_esd += distance * distance;
        }

        distance_esd = sqrt(distance_esd / double(nAtoms - 3));

        if (threshold < distance_esd)
            return Tribool::False;
        else
            return Tribool::True;

    }

    bool molecularGraphIsomorphism(
        const std::vector<int>& z1,
        const std::vector<std::vector<int> >& connectivity1,
        const std::vector<int>& z2,
        const std::vector<std::vector<int> >& connectivity2)
    {
        if (z1.size() != z2.size())
            return false;

        map<int, int> formula1, formula2;
        basic_chemistry_utilities::getFormula(z1, formula1);
        basic_chemistry_utilities::getFormula(z2, formula2);

        if (formula1 != formula2)
            return false;

        ARGEdit ed1, ed2;
        vector<int> _z1 = z1;
        vector<int> _z2 = z2;

        for (int i = 0; i < z1.size(); i++)
        {
            ed1.InsertNode(&_z1[i]);
            ed2.InsertNode(&_z2[i]);
        }

        for (int i = 0; i < z1.size(); i++)
        {
            for (int j = 0; j < connectivity1[i].size(); j++)
                ed1.InsertEdge(i, connectivity1[i][j], NULL);
            for (int j = 0; j < connectivity2[i].size(); j++)
                ed2.InsertEdge(i, connectivity2[i][j], NULL);

        }

        class NodeComparator : public AttrComparator {
        public:
            NodeComparator() {}
            ~NodeComparator() {}
            virtual bool compatible(void* a1, void* a2)
            {
                int z1 = *(int*)(a1);
                int z2 = *(int*)(a2);
                return z1 == z2;
            }
        };

        Graph g1(&ed1), g2(&ed2);
        g1.SetNodeComparator(new NodeComparator());

        VF2State search_state(&g1, &g2);
        vector<node_id> n1(z1.size()), n2(z1.size());
        int nMatchingNodes;
        bool result = match(&search_state, &nMatchingNodes, &n1[0], &n2[0]);
        return result;
    }


//------- END PLANARITY 

}

}
