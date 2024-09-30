
    // set data specific for atom type

    discamb::HC_AtomTypeParameters c_type, n_type, h_type, o_type;
    
    c_type.kappa_spherical_valence = 0.99;
    c_type.kappa_deformation_valence = 0.835;
    c_type.p_val = 1.0488;
    
    c_type.p_lm = { { 0.0000 },
                    { 0.0000,  0.0203,  0.0000 }, // P_{1,-1}, P_{1,0}, P_{1,1}
                    { 0.0000,  0.0000,  0.0578,  0.0000,  0.0532 }, // P_{2,-2}, P_{2,-1}, P_{2,0}, P_{2,1}, P_{2,2}
                    { 0.0000,  0.0000,  0.0000,  0.1060,  0.0000, -0.0742,  0.0000 },
                    { 0.0000,  0.0000,  0.0000,  0.0000, -0.0120,  0.0000,  0.0275,  0.0000,  0.0043 } };
                    