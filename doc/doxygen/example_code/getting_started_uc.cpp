#include "discamb/CrystalStructure/UnitCell.h"
#include <iostream>

int main()
{
	discamb::UnitCell unitCell(90.0, 111.781, 90.0, 5.0999, 11.9516, 5.4594);
	Vector3d positionFractional(0.30478, 0.09443, 0.23515);
	Vector3d positionCartesian;
	unitCell.fractionalToCartesian(positionFractional, positionCartesian);
	std::cout<< positionCartesian[0] << " " 
	         << positionCartesian[1] << " " 
			 << positionCartesian[2] << std::endl;
}
