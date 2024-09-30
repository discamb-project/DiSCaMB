#pragma once
#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"
#include "discamb/AtomTyping/LocalCoordinateSystem.h"

namespace discamb {

    /**
    * \addtogroup Scattering Scattering
    * @{
    */


    class AtomicDensityTransform 
    {
        Matrix3d mIdentityMatrix;
    public:
        AtomicDensityTransform(); // identity
        ~AtomicDensityTransform();
        
        virtual Matrix3d getTransformMatrix(const Crystal& crystal) const;
        static AtomicDensityTransform *create(
            const std::string &type, 
            const Crystal &crystal, 
            const std::string& settings,
            const std::pair<std::string, std::string>& clusterAtom = std::pair<std::string,std::string>(std::string(), std::string()),
            const std::string crystalAtom = std::string());
    };

    class SymmetryBasedAtomicDensityTransform: public AtomicDensityTransform 
    {
        Matrix3d mMatrix;
    public:
        SymmetryBasedAtomicDensityTransform(const Crystal &, const std::string &symmOp);
        ~SymmetryBasedAtomicDensityTransform();
        virtual Matrix3d getTransformMatrix(const Crystal &crystal) const;

    };

    class AutomaticLcsDensityTransform : public AtomicDensityTransform
    {
        Matrix3d mMatrix;
    public:
        /*
        C2 X,Y,Z C2' 1.0 auto_lcs toulene_2
        pair<string,string> clusterAtom, string crystalAtom, string clusterLabel, clusters, clustersLabels
        */
        AutomaticLcsDensityTransform(
            const Crystal&, 
            const std::string& definition,
            const std::pair<std::string, std::string> &clusterAtom,
            const std::string crystalAtom, 
            const std::string &clusterLabel,
            const std::vector< std::vector<std::pair<std::string, std::string> > >& clusters, 
            const std::vector < std::string> & clustersLabels);
        ~AutomaticLcsDensityTransform();
        virtual Matrix3d getTransformMatrix(const Crystal& crystal) const;

    };

    class UserDefinedLcsDensityTransform : public AtomicDensityTransform
    {
        Matrix3d mMatrix;
    public:
        /*
        C2 X,Y,Z C2' 1.0 auto_lcs toulene_2
        pair<string,string> clusterAtom, string crystalAtom, string clusterLabel, clusters, clustersLabels
        */
        UserDefinedLcsDensityTransform(
            const Crystal&,
            const std::string& definition,
            const std::pair<std::string, std::string>& clusterAtom,
            const std::string crystalAtom);
        ~UserDefinedLcsDensityTransform();
        virtual Matrix3d getTransformMatrix(const Crystal& crystal) const;
    private:
 
        LocalCoordinateSystemCalculator mLcsCluster, mLcsTarget;
 
        /** @}*/
    };


}