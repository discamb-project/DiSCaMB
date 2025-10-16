add_example_exec(harBankMaker
                 harBankMaker.cpp
                 FilterStructures.h
                 FilterStructures.cpp
                 AssignAtomTypes.h
                 AssignAtomTypes.cpp
                 CalcDensities.cpp
                 CalcDensities.h 
                 ChooseStructures.h
                 ChooseStructures.cpp
                 RunAll.h
                 RunAll.cpp
                 Program.h
                 Program.cpp
                 WfnCalc.cpp
                 WfnCalc.h)


target_link_libraries(harBankMaker OpenMP::OpenMP_CXX)
