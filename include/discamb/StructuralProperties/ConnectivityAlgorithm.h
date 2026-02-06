#pragma once

#include "discamb/MathUtilities/Vector3.h"
//#include "discamb/StructuralProperties/MolecularDisorder.h"
#include <string>
#include <vector>

namespace discamb {

    /**
    * \addtogroup StructuralProperties
    * @{
    */


    class ConnectivityAlgorithm
    {
    public:
        virtual ~ConnectivityAlgorithm() {};

        virtual void set(const std::string &settings) = 0;
        
        virtual void calculateConnectivity(
            const std::vector<Vector3d> &positions,
            const std::vector<int> &atomicNumbers,
            std::vector<std::vector<int> > &connectivity) const = 0 ;


        //virtual void calculateConnectivity(
        //                 const std::vector<Vector3d> &positions,
        //                 const std::vector<int> &atomicNumbers,
        //                 const MolecularDisorder &molecularDisorder,
        //                 std::vector<std::vector<int> > &connectivity) const = 0;
        //
        //virtual void calculateConnectivity(
        //                 const std::vector<Vector3d> &positions,
        //                 const std::vector<int> &atomicNumbers,
        //                 const MolecularDisorder &molecularDisorder, 
        //                 const std::vector<std::vector<int> > &disjonedAtomGroups, 
        //                 std::vector<std::vector<int> > &connectivity) const = 0;
    };
    /** @}*/
}
