
    // default constructor creates identity operation
    discamb::SpaceGroupOperation identity; 
    
    // space group operations can be set using either matrix - vector notation
    discamb::SpaceGroupOperation fourBar( discamb::Matrix3i{ 0,  1,  0,
                                                            -1,  0,  0,
                                                             0,  0, -1 });

    // or string x,y,z notation
    discamb::SpaceGroupOperation mirrorPlane( "1/2-Y, 1/2-X, Z" );
    // or alternatively
    mirrorPlane.set( { 0, -1,  0,
                      -1,  0,  0,
                       0,  0,  1}, {1/2, 1/2, 0});                       

    // lets generata all the operations necessary to define space group P -4 21 m in DiSCaMB
    vector< discamb::SpaceGroupOperation > operations = 
                           { identity                    , mirrorPlane * identity                   ,
                             fourBar                     , mirrorPlane * fourBar                    ,
                             fourBar * fourBar           , mirrorPlane * fourBar * fourBar          , 
                             fourBar * fourBar * fourBar , mirrorPlane * fourBar * fourBar * fourBar };
   
   discamb::SpaceGroup p4bar21m;
   p4bar21m.set(operations);
        

