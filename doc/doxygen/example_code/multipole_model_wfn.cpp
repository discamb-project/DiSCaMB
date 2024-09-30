    
    // set wave function related data

    discamb::DeformationValenceParameters deformationValenceParameters;
    discamb::HC_WfnType h_wfn, c_wfn, o_wfn, n_wfn;
    discamb::ClementiRoettiData clementiRoettiData;
        
    c_wfn = clementiRoettiData.getEntry("C");
    deformationValenceParameters.getParameters("C", c_wfn.deformation_valence_exponent,
                                                    c_wfn.deformation_valence_power);
