
    // lets define some hkl indices
    vector<Vector3i> hkl{ {1, 0, 0}, {1, -1, 2}, {1, 2, 3} };

    // and local coordinate systems
    vector< discamb::XdLocalCoordinateSystem > lcs(5);

    // set coordinate system for the first atom in crystal.atoms

    lcs[0].set("C1", "O1", "C1", "DUM1", "Z", "X", true, crystal, false, { { -0.144680, 0.355320, 0.179010} });

    // and get local coordinates for current geometry

    vector< Matrix3d > lcsMatrices(5);

    lcs[0].calculate(lcsMatrices[0], crystal); // .. and so on for other atoms
    