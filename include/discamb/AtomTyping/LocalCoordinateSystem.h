#pragma once

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/AtomInCrystalID.h"

#include <vector>
#include <string>

namespace discamb {
    /**
    * \addtogroup AtomTyping
    * @{
    */

    enum class LcsDirectionType { ANY_ORTHOGONAL, AVERAGE_DIRECTION, AVERAGE_POSITION, NOT_SET };

	template<typename T>
    struct LocalCoordinateSystem
    {
		
        bool isR = true;
		/* reference points  */
        LcsDirectionType direction1_type = LcsDirectionType::NOT_SET;
        LcsDirectionType direction2_type = LcsDirectionType::NOT_SET;
		T centralAtom;
		/*
		each refPoint is defined with list of atoms
		first coordinate is defined as direction from central atom to 
		refPoint_1 - defines first coordinate as direction from central atom to point defined by atom indices in refPoint_1 (transformed with refPoint_1_symmOps)
		
		if refPoint_3 is empty then the second direction  (being an averega of the atoms positions)
		*/
        // list of atoms defining reference point r, if refPoint_3 is empty then the second direction  (being an averega of the atoms positions)
		// 
        //std::vector<int> refPoint_1, refPoint_2a, refPoint_2b;
		// 
		// AtomInCrystalID for crystal and int for molecule
		std::vector<T> refPoint_1, refPoint_2;

        
        // X - 0 , Y - 1 , Z - 2
        int coordinate_1 = 0;
        int coordinate_2 = 1;

        std::vector<T> chirality;
		// central_atom coordinate_1 direction_1 refpoint_2 coordinate_2 chirality chirality_defining_atoms
		// C(1) C(2),C(3),C(4)[-x,-y,z] X any_orthogonal Y R C(2),C(3),C(6)
		// C(1) average_direction:C(2)[1-x,y,z],C(3) Z C(4) Y L C(2),C(3),C(6)

		//void set(const std::string &s, const Crystal &c);
        
    };

    void convertUbdbLcs(const  LocalCoordinateSystem<int> &lcsMolecule, const std::vector<AtomInCrystalID> &atomMap, LocalCoordinateSystem<AtomInCrystalID> &lcsCrystal);
    std::string ubdbLcsAsString(const  LocalCoordinateSystem<AtomInCrystalID> &lcs, const std::vector<std::string> &labels);
    std::string ubdbLcsAsString(const  LocalCoordinateSystem<int> &lcs, const std::vector<std::string> &labels);
    std::string ubdbLcsDirectionAsString(LcsDirectionType type, const std::vector<int> &indices, const std::vector<std::string> &labels);
    std::string ubdbLcsDirectionAsString(LcsDirectionType type, const std::vector<AtomInCrystalID> &indices, const std::vector<std::string> &labels);
    // central_atom X atom1 Y atom2 
	void xdTypeLcs(const std::string& definition, const Crystal& c, LocalCoordinateSystem<AtomInCrystalID>& lcs);

    /**@}*/
}

