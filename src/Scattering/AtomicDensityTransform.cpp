#include "discamb/Scattering/AtomicDensityTransform.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/BasicUtilities/StringUtilities.h"

using namespace std;



namespace discamb {



    AtomicDensityTransform::AtomicDensityTransform() 
    {
        mIdentityMatrix.set(1.0, 0.0, 0.0,
                            0.0, 1.0, 0.0,
                            0.0, 0.0, 1.0);
    }
        
    AtomicDensityTransform::~AtomicDensityTransform() {}

    Matrix3d AtomicDensityTransform::getTransformMatrix(
        const Crystal &crystal) 
        const
    {
        return mIdentityMatrix;
    };
        
    AtomicDensityTransform *AtomicDensityTransform::create(
        const std::string &type, 
        const Crystal &crystal, 
        const std::string& settings,
        const std::pair<std::string, std::string>& clusterAtom,
        const std::string crystalAtom)
    {
        if (type == string("identity"))
            return new AtomicDensityTransform;
        if (type == string("symmetry"))
            return new SymmetryBasedAtomicDensityTransform(crystal, settings);
        if (type == string("user_lcs"))
            return new UserDefinedLcsDensityTransform(crystal, settings, clusterAtom, crystalAtom);
        on_error::throwException(string("invalid type of AtomicDensityTransform: '") + type + string("'"), __FILE__, __LINE__);
        return NULL;
    };

    Matrix3d SymmetryBasedAtomicDensityTransform::getTransformMatrix(
        const Crystal &crystal)
        const
    {
        return mMatrix;
    }

    SymmetryBasedAtomicDensityTransform::SymmetryBasedAtomicDensityTransform(
        const Crystal &crystal,
        const std::string &symmOp)
    {
        vector<string> words;
        string_utilities::split(symmOp, words);
        SpaceGroupOperation symmetryOperation(words.back());

        if (words.size() == 2)
        {
            if (words[0] == string("invert"))
                symmetryOperation.invert();
            else
                on_error::throwException(string("invalid specification of atomic transformation for form factor calculation :'") + symmOp + string("'"), __FILE__, __LINE__);
        }

        Vector3d v;
        crystal_structure_utilities::symmetryOperationCartesian(symmetryOperation, crystal.unitCell, mMatrix, v);
    }

    SymmetryBasedAtomicDensityTransform::~SymmetryBasedAtomicDensityTransform() {}

    
    
    //#################################################################################################################


    AutomaticLcsDensityTransform::AutomaticLcsDensityTransform(
        const Crystal&,
        const std::string& symm_and_cluster,
        const std::pair<std::string, std::string>& clusterAtom,
        const std::string crystalAtom,
        const std::string& clusterLabel,
        const std::vector< std::vector<std::pair<std::string, std::string> > >& clusters,
        const std::vector < std::string>& clustersLabels)
    {
        /*
        toulene_1
        .
        .
        C2 X,Y,Z C2' 1.0 auto_lcs -X,-Y,Z toulene_2

        find lcs matrix M1 for C2 in toluene_1
        find lcs matrix M2 for C2',-X,-Y,Z in toluene_2
        M3 is -X,-Y,Z matrix

        local coordinates x,y,z vectors are columns of M1 and M2

        we look for transformation which transforms C2 density into C2' density


        */

        /* jak liczyæ M1 maj¹c polo¿enia atomów w krysztale
        liczymy geometrie klastra 
        wyznaczamy wi¹zania
        wyznaczamy lcs
        atomy typy 'H@C2 C5' zamieniamy na C5
        tlumaczymy na lcs dla krysztalu
        */

    }

    AutomaticLcsDensityTransform::~AutomaticLcsDensityTransform()
    {}

    Matrix3d AutomaticLcsDensityTransform::getTransformMatrix(
        const Crystal& crystal)
        const
    {
        return Matrix3d();
    }

    //#################################################################################################################

    UserDefinedLcsDensityTransform::UserDefinedLcsDensityTransform(
        const Crystal& crystal,
        const std::string& definition,
        const std::pair<std::string, std::string>& clusterAtom,
        const std::string crystalAtom)
    {
        // [X C3 Y C4] [X C3' Y C4']
        string lcs1str, lcs2str;
        lcs1str = definition.substr(definition.find('[')+1, definition.find(']') - definition.find('[') - 1 );
        lcs2str = definition.substr(definition.rfind('[')+1, definition.rfind(']') - definition.rfind('[') - 1 );

        string clusterAtomLabel = clusterAtom.first+string(",")+ clusterAtom.second;

        // central_atom X atom1 Y atom2 
        LocalCoordinateSystem<AtomInCrystalID> lcsCluster, lcsTarget;
        //LocalCoordinateSystemCalculator mLcsCluster, mLcsTarget;
        //LocalCoordinateSystemCalculator(const LocalCoordinateSystem<AtomInCrystalID> & lcs, const Crystal & c);

        xdTypeLcs(clusterAtomLabel + string(" ") + lcs1str, crystal, lcsCluster);
        xdTypeLcs(crystalAtom + string(" ") + lcs2str, crystal, lcsTarget);

        mLcsCluster.set(lcsCluster, crystal);
        mLcsTarget.set(lcsTarget, crystal);


        Matrix3d lcsMatrixCluster, lcsMatrixClusterTransposed, lcsMatrixTarget, lcsMatrixTargetTransposed;
        mLcsCluster.calculate(lcsMatrixCluster, crystal);
        mLcsTarget.calculate(lcsMatrixTarget, crystal);

        lcsMatrixClusterTransposed = lcsMatrixCluster;
        lcsMatrixClusterTransposed.transpose();

        //lcsMatrixTargetTransposed = lcsMatrixTarget;
        //lcsMatrixTargetTransposed.transpose();
        //result = lcsMatrixTargetTransposed * lcsMatrixCluster;
        //result = lcsMatrixClusterTransposed * lcsMatrixTarget;

        mMatrix = lcsMatrixTarget * lcsMatrixClusterTransposed;

    }

    UserDefinedLcsDensityTransform::~UserDefinedLcsDensityTransform()
    {
    }

    Matrix3d UserDefinedLcsDensityTransform::getTransformMatrix(
        const Crystal& crystal)
        const
    {
        Matrix3d lcsMatrixCluster, lcsMatrixClusterTransposed, lcsMatrixTarget, lcsMatrixTargetTransposed, result;
        mLcsCluster.calculate(lcsMatrixCluster, crystal);
        mLcsTarget.calculate(lcsMatrixTarget, crystal);

        lcsMatrixClusterTransposed = lcsMatrixCluster;
        lcsMatrixClusterTransposed.transpose();

        //lcsMatrixTargetTransposed = lcsMatrixTarget;
        //lcsMatrixTargetTransposed.transpose();
        //result = lcsMatrixTargetTransposed * lcsMatrixCluster;
        //result = lcsMatrixClusterTransposed * lcsMatrixTarget;

        //mMatrix = lcsMatrixTarget * lcsMatrixClusterTransposed;
        result = lcsMatrixTarget * lcsMatrixClusterTransposed;
        return result;
    }

}