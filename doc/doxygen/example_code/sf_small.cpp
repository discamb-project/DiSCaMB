    
    // the type implementing funtction onPerHklCalculation
    // which is used for data collection
    struct DataCollector
    {
        std::vector<std::complex<double> > structureFactorsData;
        std::vector<discamb::SfDerivativesAtHkl> derivativesData;
        
        // this function is called by discamb::HansenCoppensStructureFactorCalculator
        // for each hkl index in order to pass results for this hkl
        void onPerHklCalculation(size_t hkl_idx, 
                                std::complex<double> &structureFactor, 
                                discamb::SfDerivativesAtHkl &derivatives)
                                {
                                       structureFactorsData[hkl_idx] = structureFactor;
                                       derivativesData[hkl_idx] = derivatives;
                                }
    };
    

    DataCollector collector;
    
    // get number of hkl vectors
    size_t nHkl = hkl.size();
    
    // and allocate memory for collector structures 
    collector.structureFactorsData.resize(nHkl);
    collector.derivativesData.resize(nHkl);
    
    // set up the calculator
    discamb::HansenCoppensStructureFactorCalculator calculator(crystal, multipoleModelParameters);

    // calculate
    
    calculator.calculateStructureFactorsAndDerivatives( crystal.atoms, lcsMatrices, 
                                                        hkl, collector); 

    // print out structure factors
    
    std::cout << "h k l and structure factor value" << std::endl;
    
    for( size_t i=0; i<nHkl; i++)
        std::cout << hkl[i][0] << " " << hkl[i][1] << " " << hkl[i][2]
                  << " " << collector.structureFactorsData[i] << std::endl;
                                                        
