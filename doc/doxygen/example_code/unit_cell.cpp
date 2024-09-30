
    discamb::UnitCell unitCell(5.565, 5.565, 4.684, 90, 90, 90);
    discamb::Vector3d cartesian;
    unitCell.fractionalToCartesian(c1.coordinates, cartesian);
    std::cout << c1.type << " " << cartesian[0] << " " << cartesian[1] << " " << cartesian[2] << std::endl;
