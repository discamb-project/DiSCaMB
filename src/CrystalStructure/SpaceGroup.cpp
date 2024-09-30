#include "discamb/CrystalStructure/SpaceGroup.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/MathUtilities/algebra3d.h"

#include <algorithm>

#include <fstream>

using namespace std;

namespace {
    void add_to_log(int n, const string &s = string())
    {
       // ofstream out("f:\\tmp\\1\\log", std::ofstream::out | std::ofstream::app);
       // out << n << " " << s << endl;
       // out.close();
    }
}

namespace discamb {

SpaceGroup::SpaceGroup()
{
    set(std::vector<SpaceGroupOperation>(1, SpaceGroupOperation()));
}

SpaceGroup::~SpaceGroup()
{
}

void SpaceGroup::set(
const std::vector<SpaceGroupOperation> &_operations)
{
    Matrix3<int> rotation,identity_matrix,inversion_matrix;

    add_to_log(__LINE__, "SpaceGroup::set (list of operations): ");
    
    for (int i = 0; i < _operations.size(); i++)
    {
        string xxx;
        _operations[i].get(xxx);
        add_to_log(__LINE__, xxx);
    }
    
    std::vector< Matrix3<int> > rotations, rotations_unique_rot, rotations_subset;

    Vector3<CrystallographicRational> translation,translation_0(0, 0, 0);
    std::vector<Vector3<CrystallographicRational> > translations,centering_translations,
                                                 translation_unique_rot, translations_subset,
                                                 inversion_center_translations;

    std::vector<SpaceGroupOperation> operations_subset;

    
    mIsCentrosymmetric = false;
    
    int n_symm = _operations.size();
    int n_centering_translations;
    int n_symm_P;

    // set translations components to [0,1) interval

    std::vector<SpaceGroupOperation> operations;
    
    for(int i=0; i<n_symm ; i++)
    {
        _operations[i].get(rotation, translation);
        moveTo01(translation);
        operations.push_back(SpaceGroupOperation(rotation,translation));
    }

    // end of set translations to 0,1


    // check for centering and inversion center, collect centering translations

    inversion_matrix.set(-1, 0, 0,
                          0,-1, 0,
                          0, 0,-1);

    identity_matrix.setToIdentity();

    rotations.resize(n_symm);
    translations.resize(n_symm);

    for(int i=0;i<n_symm;++i)
    {
        operations[i].get(rotations[i],translations[i]);
        
        if( inversion_matrix == rotations[i])
        {
            mIsCentrosymmetric = true;
            inversion_center_translations.push_back(translations[i]);
        }
        if( identity_matrix == rotations[i])
            centering_translations.push_back(translations[i]);
    }

    mCentering = latticeCentering(centering_translations,mObverse);
    n_centering_translations = centering_translations.size();

    
    // subset of operations with unique rotations is chosen
    // stored in rotations_unique_rot and translation_unique_rot

    std::vector<bool> is_processed(n_symm, false);

    for (int i = 0; i<n_symm; ++i)
        if(!is_processed[i])
        {
            is_processed[i]=true;
            rotations_unique_rot.push_back(rotations[i]);
            translation_unique_rot.push_back(translations[i]);
            for (int j = i+1; j<n_symm; ++j)
                if(rotations[i]==rotations[j])
                    is_processed[j] = true;
        }


    // divide operations such that one set contains operations inversion*operation_from_the_other_set
    
    n_symm_P = rotations_unique_rot.size();

    if(mIsCentrosymmetric)
    {
        // set inversion center translation
        for (int i=0;i<inversion_center_translations.size(); ++i)
            moveTo01(inversion_center_translations[i]);
        sort(inversion_center_translations.begin(), inversion_center_translations.end());
        mInversionCenterTranslation = inversion_center_translations[0];

        /*
        Group G containig inversion center can be divided into two sets with equal number of elements - 
        one with operations of 1-st kind (det(matrix)=1), ant the other one with operations of the second kind
        (det(matrix)=-1). The first set H is a subgroup of G, and the second set is inversion_centre x H.
        */

        int determinant;
        
        for(int i=0;i<n_symm_P;i++)
        {
            determinant = algebra3d::det3d(rotations_unique_rot[i]);

            if (determinant != -1 && determinant != 1)
                on_error::throwException("problem when processing space group operations, absolute value of rotation determinant not equal 1", __FILE__, __LINE__);

            if(determinant==1)
            {
                rotations_subset.push_back(rotations_unique_rot[i]);
                translations_subset.push_back(translation_unique_rot[i]);
            }

        }

    }
    else
    {
        rotations_subset.swap(rotations_unique_rot);
        translations_subset.swap(translation_unique_rot);
    }

    // replace symmetry operations of the form S = {R|t} by S' = {R|0} if S' belongs to the space group

    removeCenteringTranslations(translations_subset, mCentering, mObverse);

    // make symmetry operations from translations and rotations

    n_symm = rotations_subset.size();
    for(int i=0;i<n_symm;i++)
        operations_subset.push_back(SpaceGroupOperation(rotations_subset[i],translations_subset[i]));
    

    set(operations_subset, mCentering, mObverse, mIsCentrosymmetric, mInversionCenterTranslation);
}

void SpaceGroup::removeCenteringTranslations(
    std::vector<Vector3<CrystallographicRational> > &t,
    char centering, 
    bool obverse)
{
    std::vector<Vector3<CrystallographicRational> > translations;
    centeringTranslations(centering, translations, obverse);

    if(translations.size()==1)
        return;

    Vector3<CrystallographicRational> zeroTranslation;
    
    // remove translation (0,0,0) from set of centering translations
    for(int i=0;i<translations.size();i++)
        if(translations[i] == zeroTranslation)
        {
            translations.erase(translations.begin()+i);
            break;
        }

    // find out if translation in each symmetry operation is one of 
    // centering translations and if yes change it to (0,0,0)

    Vector3<CrystallographicRational> translation;
    Matrix3i rotation;
    std::vector<Vector3<CrystallographicRational> >::iterator it;
    int nSymmOps = t.size();

    for(int i=0;i<nSymmOps;i++)
    {
        it = find(translations.begin(), translations.end(),t[i]);
        if (it != translations.end())
            t[i] = zeroTranslation;
    }

}

void SpaceGroup::setPlainSpaceGroupOperations()
{
    mOperationsPlain.clear();
    int i,j,k,n_centro,n_translations,n;

    mIsCentrosymmetric ? n_centro=2 : n_centro=1;
    n_translations = mOperations.size();
    n_centro = mOperations[0].size();
    n = mOperations[0][0].size();

    for(i=0;i<n_translations;++i)
        for(j=0;j<n_centro;j++)
            for(k=0;k<n;k++)
                mOperationsPlain.push_back(mOperations[i][j][k]);

}

bool SpaceGroup::obverseSetting()
const
{
    return mObverse;
}

void SpaceGroup::centeringTranslations(
char centering, 
std::vector<Vector3<CrystallographicRational> > &centeringTranslations,
bool obverse)
{
    centeringTranslations.clear();
    centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(0), CrystallographicRational(0), CrystallographicRational(0)));

    if (centering == 'A')
        //centeringTranslations.push_back({ { 0 },{ 1,2 },{ 1,2 } }); C++11
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(0), CrystallographicRational(1,2), CrystallographicRational(1,2)));
    if (centering == 'B')
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1,2), CrystallographicRational(0), CrystallographicRational(1,2)));
    if (centering == 'C')
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1,2), CrystallographicRational(1,2), CrystallographicRational(0)));
    if (centering == 'I')
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1,2), CrystallographicRational(1,2), CrystallographicRational(1,2)));
    if (centering == 'F')
    {
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(0), CrystallographicRational(1,2), CrystallographicRational(1,2)));
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1,2), CrystallographicRational(0), CrystallographicRational(1,2)));
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1,2), CrystallographicRational(1,2), CrystallographicRational(0)));
    }
    if(centering == 'R')
    {
        if(obverse)
        {
            centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1, 3), CrystallographicRational(2, 3), CrystallographicRational(2, 3)));
            centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(2, 3), CrystallographicRational(1, 3), CrystallographicRational(1, 3)));
        }
        else
        {
            centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1, 3), CrystallographicRational(2, 3), CrystallographicRational(1, 3)));
            centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(2, 3), CrystallographicRational(1, 3), CrystallographicRational(2, 3)));
        }
    }
    if (centering == 'H')
    {
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(1, 3), CrystallographicRational(2, 3), CrystallographicRational(0, 1)));
        centeringTranslations.push_back(Vector3<CrystallographicRational>(CrystallographicRational(2, 3), CrystallographicRational(1, 3), CrystallographicRational(0, 1)));
    }

}

char SpaceGroup::latticeCentering(
const std::vector<Vector3<CrystallographicRational> > &centering_translations,
bool &obverse)
{
    // identify centering translation

    for (int i = 0; i < centering_translations.size(); i++)
    {
        string xxx;
        CrystallographicRational x, y, z;
        centering_translations[i].get(x, y, z);
        xxx = x.asString() + string(", ") + y.asString() + string(", ") + z.asString();
        add_to_log(__LINE__, xxx);
    }


    Vector3<CrystallographicRational> translation;
    std::vector<Vector3<CrystallographicRational> > translationsCanonical;
    std::vector<Vector3<CrystallographicRational> > P,A,B,C,I,F,H,R_obverse,R_reverse;

    centeringTranslations('A', A, false);
    centeringTranslations('B', B, false);
    centeringTranslations('C', C, false);
    centeringTranslations('F', F, false);
    centeringTranslations('H', H, false);
    centeringTranslations('I', I, false);
    centeringTranslations('P', P, false);
    centeringTranslations('R', R_obverse, true);
    centeringTranslations('R', R_reverse, false);

    int n_translations = centering_translations.size();
    translationsCanonical = centering_translations;

    for (int i = 0; i<n_translations ; i++)
        moveTo01(translationsCanonical[i]);

    sort(translationsCanonical.begin(), translationsCanonical.end());

    if (translationsCanonical == P)
        return 'P';
    if (translationsCanonical == A)
        return 'A';
    if (translationsCanonical == B)
        return 'B';
    if (translationsCanonical == C)
        return 'C';
    if (translationsCanonical == F)
        return 'F';
    if (translationsCanonical == H)
        return 'H';
    if (translationsCanonical == I)
        return 'I';
    if (translationsCanonical == R_obverse)
    {
        obverse = true;
        return 'R';
    }
    if (translationsCanonical == R_reverse)
    {
        obverse = false;
        return 'R';
    }
    
    on_error::throwException("problem when defining space group, symmetry operations with unknown lattice centering ",__FILE__,__LINE__);

    return 'X';
}

void SpaceGroup::moveTo01(
Vector3<CrystallographicRational> &translation)
{
    int denominator,numerator;

    for (int j = 0; j<3; j++)
    {
        denominator = translation[j].denominator();
        numerator = translation[j].numerator();
        if (denominator<0)
        {
            denominator *= -1;
            numerator *= -1;
        }

        numerator = numerator % denominator;
        if (numerator<0)
            numerator += denominator;

        translation[j].assign(numerator, denominator);
    }
    
}



int SpaceGroup::nSymmetryOperations() 
const
{
    return mOperationsPlain.size();
}

const SpaceGroupOperation &SpaceGroup::getSpaceGroupOperation(
int index)
const
{
    return mOperationsPlain[index];
}

int SpaceGroup::nSymmetryOperationsInSubset() 
const
{
    return mOperations[0][0].size();
}

bool SpaceGroup::isCentrosymmetric() 
const
{
    return mIsCentrosymmetric;
}


void SpaceGroup::inversionCenterTranslation(
Vector3<CrystallographicRational> &translation)
const
{
    translation = mInversionCenterTranslation;
}

char SpaceGroup::centering(
    bool &obverse)
const
{
    obverse = mObverse;
    return mCentering;
}

int SpaceGroup::nCenteringNodes() 
const
{
    if( mCentering == 'P' )
        return 1;
    if( mCentering == 'A' || mCentering == 'B' || mCentering == 'C' || mCentering == 'I' )
        return 2;
    if( mCentering == 'R')
        return 3;
    if( mCentering == 'F' )
        return 4;

    return 0;
}

const SpaceGroupOperation &SpaceGroup::getSpaceGroupOperation(
int centering_node, 
int inversion, 
int index) 
const
{
    return mOperations[centering_node][inversion][index];
}

void SpaceGroup::show_info(
std::ostream &out) 
const
{
    out<< "centering : " << mCentering << std::endl;
    out<< "has inversion center : ";
    if(mIsCentrosymmetric)
        out<< "yes";
    else
        out<< "no";
    out<< std::endl << "inversion center translation : ";
    if (mIsCentrosymmetric)
        out<< mInversionCenterTranslation[0].asString() << " "
           << mInversionCenterTranslation[1].asString() << " "
           << mInversionCenterTranslation[2].asString();
    else
        out<< "-";

    out<< std::endl << "number of operations (subset): " << mOperations[0][0].size() << std::endl;

    std::string symmOperationsAsString;
    for(int i=0;i<mOperations[0][0].size();++i)
    {
        mOperations[0][0][i].get(symmOperationsAsString);
        out << symmOperationsAsString << std::endl;
    }
    
    out << std::endl << "number of operations (all): " << mOperationsPlain.size() << std::endl;

    for (int i = 0; i<mOperationsPlain.size(); ++i)
    {
        mOperationsPlain[i].get(symmOperationsAsString);
        out << symmOperationsAsString << std::endl;
    }

}

void SpaceGroup::set(
const std::vector<SpaceGroupOperation> &operations,
char centering,
bool inversion_center)
{
    Vector3<CrystallographicRational> inv_center_translation(CrystallographicRational(0,1), CrystallographicRational(0,1), CrystallographicRational(0,1));
    set(operations, centering, false, inversion_center, inv_center_translation);
}


void SpaceGroup::set(
const std::vector<SpaceGroupOperation> &operations, 
char centering,
bool obverse,
bool inversion_center, 
const Vector3<CrystallographicRational> &inv_center_translation)
{
    mOperations.clear();
    mCentering = centering;
    mIsCentrosymmetric = inversion_center;
    mInversionCenterTranslation = inv_center_translation;
    mObverse = obverse;

    int n_centering, n_translations, n_symm, translation_idx, symm_op_idx;
    Vector3<CrystallographicRational> translation;
    std::vector<Vector3<CrystallographicRational> > centering_translations;
    Matrix3i rotation;
    
    centeringTranslations(centering, centering_translations, obverse);
    n_translations = centering_translations.size();
    n_symm = operations.size();
    inversion_center ? n_centering=2 : n_centering=1;
    mOperations.resize(n_translations);

    for(translation_idx=0; translation_idx<n_translations; translation_idx++)
    {
        mOperations[translation_idx].resize(n_centering);
        for(symm_op_idx=0;symm_op_idx<n_symm;symm_op_idx++)
        {
            operations[symm_op_idx].get(rotation,translation);
            translation += centering_translations[translation_idx];
            moveTo01(translation);
            mOperations[translation_idx][0].push_back(SpaceGroupOperation(rotation,translation));
            if(inversion_center)
            {
                translation = inv_center_translation - translation;
                moveTo01(translation);
                mOperations[translation_idx][1].push_back(SpaceGroupOperation(-rotation,translation));
            }
        }
    }
    setPlainSpaceGroupOperations();
}

}//namespace discamb{
