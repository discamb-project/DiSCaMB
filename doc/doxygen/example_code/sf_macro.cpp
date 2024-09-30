
    // the derivatives of target function with respect to structure factor components
    vector<complex<double> > dT_dF{ { 1.2, -0.34 }, { 0.13, 3.1 }, { 1, 2 } };

    // the deravatives to be computed - derivatives of target function
    // with respect to structural parameters
    std::vector< TargetFunctionAtomicParamDerivatives > derivatives;
    
    // container for structure factors
    vector<complex<double> > structureFactors;

    // set up the calculator
    discamb::HansenCoppensStructureFactorCalculator calculator(crystal, multipoleModelParameters);

    // calculate
    calculator.calculateStructureFactorsAndDerivatives( crystal.atoms, lcsMatrices, 
                                                        hkl, structureFactors, derivatives,
                                                        dT_dF, vector<bool>(5, true));
    