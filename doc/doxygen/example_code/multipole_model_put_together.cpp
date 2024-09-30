
    // define multipole model parameters

    discamb::HC_ModelParameters multipoleModelParameters;

    multipoleModelParameters.wfn_parameters = { c_wfn, o_wfn, n_wfn, h_wfn };
    multipoleModelParameters.atom_to_wfn_map = { 0, 1, 2, 3, 3 };
    multipoleModelParameters.type_parameters = { c_type, o_type, n_type, h_type };
    multipoleModelParameters.atom_to_type_map = { 0, 1, 2, 3, 3 };

    // alternatively it can be read from XD files
    
    vector<XdLocalCoordinateSystem> lcs;
    discamb::xd_io::read("xd.mas", "xd.par", multipoleModelParameters, crystal, lcs, false, true);
    