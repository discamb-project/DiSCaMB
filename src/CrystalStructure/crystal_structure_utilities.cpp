#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/BasicChemistry/basic_chemistry_utilities.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"
#include "discamb/MathUtilities/algebra3d.h"
#include "discamb/MathUtilities/lebedev_laikov.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/StructuralProperties/structural_properties.h"
#include "discamb/CrystalStructure/StructuralParametersConverter.h"

#include <limits>
#include <cmath>
#include <cfenv>
#include <random>

using namespace std;

namespace{

// returns such integer n, thaty x+n belongs to interval [0, 1) 

int integerShiftTo01(
    double x)
{
    if (x < 0)
        return abs(int(x)) + 1;
    return -abs(int(x));
}

void integerShiftTo01(
    const discamb::Vector3d &r,
    discamb::Vector3i &t)
{
    for (int i = 0; i < 3; i++)
        t[i] = integerShiftTo01(r[i]);
}

// calculates distance btween two points in cristal given in fractional coordinates

double calcDistance(
    const discamb::Vector3d &r1,
    const discamb::Vector3d &r2,
    const discamb::Crystal &crystal)
{
    discamb::Vector3d r1cart, r2cart, diff;

    crystal.unitCell.fractionalToCartesian(r1, r1cart);
    crystal.unitCell.fractionalToCartesian(r2, r2cart);

    diff = r1cart - r2cart;

    return std::sqrt(diff*diff);
}

// finds the shortest distance between r1 and r2 + t, where t is lattice translation
// vectors are in fractional coordinates

double findShortestDistance(
    const discamb::Vector3d &r1,
    const discamb::Vector3d &r2,
    const discamb::Crystal &crystal,
    //const discamb::Vector3d &a,
    //const discamb::Vector3d &b,
    //const discamb::Vector3d &c,
    discamb::Vector3i &translation)
{
    double minDistance, distance;
    int i, j, k, l;
    discamb::Vector3i t1, t2, t0, t0_min;
    discamb::Vector3d r1_0, r2_0,r2_t;

    integerShiftTo01(r1, t1);
    integerShiftTo01(r2, t2);

    for (i = 0; i < 3; i++)
    {
        r1_0[i] = r1[i] + t1[i];
        r2_0[i] = r2[i] + t2[i];
    }

    minDistance = std::numeric_limits<double>::max();

    for(i = -1; i<2; i++)
        for (j = -1; j<2; j++)
            for (k = -1; k < 2; k++)
            {
                t0[0] = i;
                t0[1] = j;
                t0[2] = k;

                for (l = 0; l < 3; l++)
                    r2_t[l] = r2_0[l] + t0[l];
                distance = calcDistance(r1_0, r2_t, crystal);
                if (distance < minDistance)
                {
                    minDistance = distance;
                    t0_min = t0;
                }
            }
    translation = t0_min + t2 - t1;
    return minDistance;
}

} // namespace

namespace discamb {

    namespace crystal_structure_utilities {


        void distort_structure(
            const Crystal& crystal,
            Crystal& distortedCrystal,
            const StructureDistortionParameters& distortionParameters)
        {
            std::default_random_engine generator;
            std::uniform_real_distribution<double> distribution(-1.0, 1.0);

            double number = distribution(generator);

            distortedCrystal = crystal;
            for (auto& atom : distortedCrystal.atoms)
            {
                if (distortionParameters.distort_positions)
                {
                    if (distortionParameters.useCartesianShift)
                    {
                        Vector3d shiftCart;
                        for (int i = 0; i < 3; i++)
                            shiftCart[i] = distribution(generator) * distortionParameters.maxCartesianShift;

                        Vector3d fracShift;
                        distortedCrystal.unitCell.cartesianToFractional(shiftCart, fracShift);
                        atom.coordinates += fracShift;
                    }
                    else
                    {
                        for (int i = 0; i < 3; i++)
                            atom.coordinates[i] += number * distortionParameters.maxFractionalShift;
                    }
                }
                if (distortionParameters.distort_adps)
                {
                    if (atom.adp.size() == 1)
                    {
                        double delta_u = distribution(generator) * distortionParameters.maxAdpChange;
                        double newUiso = atom.adp[0] + delta_u;
                        if (newUiso < 0)
                            newUiso = 0.0;
                        atom.adp[0] = newUiso;
                    }

                    if (atom.adp.size() == 6)
                    {
                        auto& u = atom.adp;
                        Matrix3d uMatrix(u[0], u[3], u[4],
                            u[3], u[1], u[5],
                            u[4], u[5], u[2]);
                        Vector3d eigVec[3];
                        double eigVal[3];
                        algebra3d::eigensystemRealSymm(uMatrix, eigVec[0], eigVec[1], eigVec[2], eigVal[0], eigVal[1], eigVal[2]);
                        for (int i = 0; i < 3; i++)
                        {
                            double delta_u = distribution(generator) * distortionParameters.maxAdpChange;
                            eigVal[i] += delta_u;
                            if (eigVal[i] < 0.0)
                                eigVal[i] = 0.0;
                        }
                        // reconstruct U matrix
                        Matrix3d uMatrixNew = eigVal[0] * algebra3d::outer(eigVec[0], eigVec[0]) +
                            eigVal[1] * algebra3d::outer(eigVec[1], eigVec[1]) +
                            eigVal[2] * algebra3d::outer(eigVec[2], eigVec[2]);
                        // set back to atom.adp
                        atom.adp[0] = uMatrixNew(0, 0);
                        atom.adp[1] = uMatrixNew(1, 1);
                        atom.adp[2] = uMatrixNew(2, 2);
                        atom.adp[3] = uMatrixNew(0, 1);
                        atom.adp[4] = uMatrixNew(0, 2);
                        atom.adp[5] = uMatrixNew(1, 2);
                    }
                }
                if (distortionParameters.distort_occupancies)
                {
                    double delta_occupancy = distribution(generator) * distortionParameters.maxOccupancyChange;
                    double newOccupancy = atom.occupancy + delta_occupancy;
                    if (newOccupancy < 0.0)
                        newOccupancy = 0.0;
                    if (newOccupancy > 1.0)
                        newOccupancy = 1.0;
                    atom.occupancy = newOccupancy;
                }
            }

        }


       /* void translateAtomTo01(
            const Crystal &crystal,
            int atomIdx,
            SpaceGroupOperation &operation,
            SpaceGroupOperation &operationAfterTranslation,
            Vector3i &translation)
        {


        }*/

        double u_eq(
            const UnitCell& unit_cell,
            const std::vector<double>& uij_cif)
        {
            if (uij_cif.size() == 0)
                return 0.0;

            if (uij_cif.size() == 1)
                return uij_cif[0];

            StructuralParametersConverter converter;
            converter.set(unit_cell);
            vector<double> u_cart(6);
            converter.convertADP(uij_cif, u_cart, structural_parameters_convention::AdpConvention::U_cif,
                structural_parameters_convention::AdpConvention::U_cart);
            return (u_cart[0] + u_cart[1] + u_cart[2]) / 3.0;
        }

        void u_eq(
            const Crystal& crystal,
            std::vector<double>& u_eqivalent)
        {
            u_eqivalent.clear();
            StructuralParametersConverter converter;
            converter.set(crystal.unitCell);
            vector<double> u_cart(6);

            for (const auto& atom : crystal.atoms)
            {
                if (atom.adp.size() == 0)
                    u_eqivalent.push_back(0.0);

                if (atom.adp.size() == 1)
                    u_eqivalent.push_back(atom.adp[0]);
                    
                if (atom.adp.size() > 1)
                {
                    converter.convertADP(atom.adp, u_cart, structural_parameters_convention::AdpConvention::U_cif,
                        structural_parameters_convention::AdpConvention::U_cart);
                    u_eqivalent.push_back((u_cart[0] + u_cart[1] + u_cart[2]) / 3.0);
                }
            }
        }

        void stringToAtomList(
            const std::string& str,
            std::vector<std::pair<std::string, std::string> >& atomList,
            char separator)
        {
            atomList.clear();
            vector<string> words, words2;
            string_utilities::split(str, words, separator);

            for (auto& word : words)
            {
                string_utilities::split(word, words2, ',');
                if (words2.size() == 4)
                    atomList.push_back({ words2[0], word.substr(words2[0].size() + 1) });
                else
                {
                    if (words2.size() == 1)
                        atomList.push_back({ words2[0], string("X,Y,Z") });
                    else
                        on_error::throwException(string("invalid entry for 'remove' key in aspher.json file: '") + str + string("'"), __FILE__, __LINE__);

                }
            }

        }

        /*
        Cartesian position of atom in crystal
        */
        Vector3d atomPosition(
            int atomIdx,
            const SpaceGroupOperation& symmOp,
            const Crystal crystal)
        {
            Vector3d result, transformedFractional;

            symmOp.apply(crystal.atoms[atomIdx].coordinates, transformedFractional);
            crystal.unitCell.fractionalToCartesian(transformedFractional, result);

            return result;
        }

        /*
        Cartesian position of atom in crystal
        */
        Vector3d atomPosition(
            const std::string& atomLabel,
            const SpaceGroupOperation& symmOp,
            const Crystal crystal)
        {
            int idx = crystal.atomIdx(atomLabel);
            return atomPosition(idx, symmOp, crystal);
        }

        void splitIntoAtomAndSymmOp(
            const std::string& s,
            std::string& atomLabel,
            std::string& symmOp,
            bool throwException)
        {
            if (!splitIntoAtomAndSymmOp(s, atomLabel, symmOp) && throwException)
                on_error::throwException("encountered invalid definition of atom in crystal: '" + s + "'", __FILE__, __LINE__);
        }

        void splitIntoAtomAndSymmOp(
            const std::vector<std::string>& s,
            std::vector < std::pair<std::string, std::string> >& atoms,
            bool throwException)
        {
            atoms.clear();
            string symmOperationStr, atomLabel;
            for (auto const& word : s)
            {
                splitIntoAtomAndSymmOp(word, atomLabel, symmOperationStr, throwException);
                atoms.push_back({ atomLabel, symmOperationStr });
            }
        }


        bool splitIntoAtomAndSymmOp(
            const string& s,
            string& atomLabel,
            string& symmOp)
        {
            vector<string> words;
            discamb::string_utilities::split(s, words, ',');
            if (words.size() != 4 && words.size()!=1)
                return false;


            atomLabel = words[0];

            if (words.size() == 1)
            {
                symmOp = "X,Y,Z";
                return true;
            }

            symmOp = s.substr(words[0].size() + 1);
            return true;
        }


        // r = rInUnitCell + unitCellOrigin, 0 <= rInUnitCell[i] < 1 

        void extractLatticeVector(
            const Vector3d &r,
            Vector3d &rInUnitCell,
            Vector3i &unitCellOrigin)
        {
            std::fesetround(FE_DOWNWARD);

            for (int i = 0; i < 3; i++)
            {
                unitCellOrigin[i] = lrint(r[i]);
                rInUnitCell[i] = r[i] - unitCellOrigin[i];
            }

        }

        void getLatticeVectors(
            const UnitCell& unitCell, 
            Vector3d& a, 
            Vector3d& b, 
            Vector3d& c)
        {
            unitCell.fractionalToCartesian(Vector3d(1, 0, 0), a);
            unitCell.fractionalToCartesian(Vector3d(0, 1, 0), b);
            unitCell.fractionalToCartesian(Vector3d(0, 0, 1), c);
        }

        void convertToXyzAndElementList(
            const Crystal& crystal,
            const std::vector<std::pair<std::string, std::string> >& labelAndSymmOperation,
            std::vector<ChemicalElement>& symbols,
            std::vector<Vector3d>& positions)
        {
            map<string, int> label2idx;

            for (int atomIdx = 0; atomIdx < crystal.atoms.size(); atomIdx++)
                label2idx[crystal.atoms[atomIdx].label] = atomIdx;

            Vector3d frac, cart;
            SpaceGroupOperation spaceGroupOperation;

            positions.clear();
            symbols.clear();

            for (auto const& atom : labelAndSymmOperation)
            {
              
                 if (label2idx.find(atom.first) != label2idx.end())
                 {
                     int atomIdx = label2idx[atom.first];
                     if (!atom.second.empty())
                     {
                         spaceGroupOperation.set(atom.second);
                         spaceGroupOperation.apply(crystal.atoms[atomIdx].coordinates, frac);
                     }
                     else
                         frac = crystal.atoms[atomIdx].coordinates;
                     crystal.unitCell.fractionalToCartesian(frac, cart);
                     positions.push_back(cart);
                     symbols.push_back(crystal.atoms[atomIdx].type);
                 }
                 else
                     on_error::throwException(string("invalid atom label '") + atom.first + string("'"), __FILE__, __LINE__);
                
            }

        }

        void convertToXyzAndElementList(
            const Crystal& crystal,
            const UnitCellContent& ucContent,
            std::vector<UnitCellContent::AtomID>& atoms,
            std::vector<ChemicalElement>& element,
            std::vector<Vector3d>& position)
        {
            vector<pair<string, string> > labelAndSymmOperation;
            string label, symmetryOperation;
            
            for (auto atom : atoms)
            {
                ucContent.interpreteAtomID(atom, label, symmetryOperation);
                labelAndSymmOperation.push_back({ label, symmetryOperation });
            }
            convertToXyzAndElementList(crystal, labelAndSymmOperation, element, position);
            
        }

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector<UnitCellContent::AtomID>& atomList,
            std::vector < std::pair < std::string, std::string > >& atomListConverted)
        {
            atomListConverted.clear();
            string label, symmOp;
            for (auto const& atom : atomList)
            {
                ucContent.interpreteAtomID(atom, label, symmOp);
                atomListConverted.push_back({ label, symmOp });
            }
        }

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector < std::pair < std::string, std::string > >& atomList,
            std::vector<UnitCellContent::AtomID>& atomsListConverted)
        {
            atomsListConverted.clear();
            UnitCellContent::AtomID atomId;
            for (auto const& atom : atomList)
            {
                ucContent.findAtom(atom.first, atom.second, atomId);
                atomsListConverted.push_back(atomId);
            }
        }

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector<UnitCellContent::AtomID>& atomList,
            std::vector < std::pair < int, std::string > >& atomListConverted)
        {
            vector < pair < string, string > > atomListConverted0;
            convertAtomList(ucContent, atomList, atomListConverted0);
            convertAtomList(ucContent.getCrystal(), atomListConverted0, atomListConverted);
        }

        void convertAtomList(
            const UnitCellContent& ucContent,
            const std::vector < std::pair < int, std::string > >& atomList,
            std::vector<UnitCellContent::AtomID>& atomListConverted)
        {
            vector < pair < string, string > > atomListConverted0;
            convertAtomList(ucContent.getCrystal(), atomList, atomListConverted0);
            convertAtomList(ucContent, atomListConverted0, atomListConverted);

        }

        void convertAtomList(
            const Crystal& crystal,
            const vector < pair < int, string > >& atomList,
            vector < pair < string, string > >& atomListConverted)
        {
            atomListConverted.clear();
            for (auto const& atom : atomList)
                atomListConverted.push_back({ crystal.atoms[atom.first].label, atom.second });
        }

        void convertAtomList(
            const Crystal& crystal,
            const std::vector < std::pair < std::string, std::string > >& atomList,
            std::vector < std::pair < int, std::string > >& atomListConverted)
        {
            atomListConverted.clear();

            map<string, int> label2idx;
            for (int idx = 0; idx<int(crystal.atoms.size()); idx++)
                label2idx[crystal.atoms[idx].label] = idx;

            for (auto const& atom : atomList)
                atomListConverted.push_back({ label2idx[atom.first], atom.second });
        }

        void convertAtomList(
            const Crystal& crystal,
            const std::vector<std::pair<std::string, std::string> >& atomList,
            std::vector<AtomInCrystalID>& atomListConverted)
        {
            vector<pair<int, string > > atomListIntStr;
            convertAtomList(crystal, atomList, atomListIntStr);
            convertAtomList(atomListIntStr, atomListConverted);
        }

        void convertAtomList(
            const vector<pair<int, string> >& atomList,
            vector<AtomInCrystalID>& atomListConverted)
        {
            atomListConverted.clear();
            for (auto const& atom : atomList)
                atomListConverted.push_back(AtomInCrystalID(atom.first, SpaceGroupOperation(atom.second)));
        }



        double interatomicDistance(
            const Crystal& crystal,
            const std::string& atomLabel1,
            const std::string& symmOp1,
            const std::string& atomLabel2,
            const std::string& symmOp2)
        {
            int atomIdx1, atomIdx2;
            atomIdx1 = atomIdx2 = 0;
            bool found1, found2;
            found1 = false;
            found2 = false;

            for (int i = 0; i < crystal.atoms.size(); i++)
            {
                if (crystal.atoms[i].label == atomLabel1)
                {
                    found1 = true;
                    atomIdx1 = i;
                }

                if (crystal.atoms[i].label == atomLabel2)
                {
                    found2 = true;
                    atomIdx2 = i;
                }

            }

            if (!found1)
                on_error::throwException("cannot find an atom with label '" + atomLabel1, __FILE__, __LINE__);
            if (!found2)
                on_error::throwException("cannot find an atom with label '" + atomLabel2, __FILE__, __LINE__);

            Vector3d r1, r2, frac, fracTransformed;
            SpaceGroupOperation symmOp;
            symmOp.set(symmOp1);
            frac = crystal.atoms[atomIdx1].coordinates;
            symmOp.apply(frac, fracTransformed);
            crystal.unitCell.fractionalToCartesian(fracTransformed, r1);

            symmOp.set(symmOp2);
            frac = crystal.atoms[atomIdx2].coordinates;
            symmOp.apply(frac, fracTransformed);
            crystal.unitCell.fractionalToCartesian(fracTransformed, r2);

            Vector3d diff = r1 - r2;
            return sqrt(diff * diff);
        }

        void atomicNumbers(
            const Crystal &crystal, 
            std::vector<int> &atomicNumbers)
        {
            atomicNumbers.clear();
            int atomicNumber;
            for (auto const &atom : crystal.atoms)
            {
                atomicNumber = basic_chemistry_utilities::atomicNumberFromString(atom.type);
                if (atomicNumber == 0)
                    atomicNumber = basic_chemistry_utilities::atomicNumberFromString(atom.label);

                atomicNumbers.push_back(atomicNumber);
            }

        }

        bool findAtomSymmetry(
            const Crystal &crystal,
            int atomIdx,
            std::vector<std::vector<SpaceGroupOperation> > &pointGroups,
            double threshold)
        {
            int operationIdx, nOperations, groupIdx, nGroups, memberIdx, nMembers;
            //std::vector< std::vector<Vector3d> > positionsCart;
            std::vector< std::vector<Vector3d> > positionFrac;
            //std::vector< std::vector<SpaceGroupOperation> > symmOps;
            Vector3d rCart, rFrac, rTransformedFractional, rTransformedCartesian;
            bool addedToGroup;
            Vector3i translation, translation01;
            double distance;
            Matrix3i identity(1, 0, 0, 0, 1, 0, 0, 0, 1);
            // 
            int firstCloseGroup;
            int firstCloseMember;

            pointGroups.clear();
            
            if (crystal.xyzCoordinateSystem == structural_parameters_convention::XyzCoordinateSystem::fractional)
                rFrac = crystal.atoms[atomIdx].coordinates;
            else
                crystal.unitCell.cartesianToFractional(crystal.atoms[atomIdx].coordinates, rFrac);

            nOperations = crystal.spaceGroup.nSymmetryOperations();

            for (operationIdx = 0; operationIdx < nOperations; operationIdx++)
            {
                firstCloseGroup = -1;
                firstCloseMember = -1;

                //crystal.spaceGroup.getSpaceGroupOperation(operationIdx).apply(rFrac, rTransformedFractional);
                auto sgOp = crystal.spaceGroup.getSpaceGroupOperation(operationIdx);
                sgOp.apply(rFrac, rTransformedFractional);

                integerShiftTo01(rTransformedFractional, translation01);
                rTransformedFractional += translation01;
                addedToGroup = false;

                //groupIdx, nGroups, memberIdx, nMembers

                nGroups = positionFrac.size();

                for (groupIdx = 0; groupIdx < nGroups; groupIdx++)
                {
                    nMembers = positionFrac[groupIdx].size();

                    for (memberIdx = 0; memberIdx < nMembers; memberIdx++)
                    {
                        distance = findShortestDistance(positionFrac[groupIdx][memberIdx], rTransformedFractional, crystal, translation);

                        if (distance < threshold)
                        {
                            // already added to another group
                            if (addedToGroup && firstCloseGroup != groupIdx)
                                    return false;
                            else
                                firstCloseGroup = int(groupIdx);

                            // previous member of this grouop was not close enough
                            if (memberIdx != 0 && !addedToGroup) 
                                return false;
                            

                            if (!addedToGroup)
                            {
                                addedToGroup = true;

                                positionFrac[groupIdx].push_back(rTransformedFractional);

                                pointGroups[groupIdx].push_back(SpaceGroupOperation(identity, translation + translation01)*
                                    crystal.spaceGroup.getSpaceGroupOperation(operationIdx));
                            }
                        } // if (distance < threshold)
                        else
                        {
                            // not close enough to member of its own group
                            if (firstCloseGroup == groupIdx)
                                return false;
                        }

                    } // for (memberIdx = 0;...
                } // for (groupIdx = 0;...

                if (!addedToGroup)
                {
                    pointGroups.resize(pointGroups.size() + 1);
                    pointGroups.back().push_back(SpaceGroupOperation(identity, translation01) *
                                                 crystal.spaceGroup.getSpaceGroupOperation(operationIdx));
                    positionFrac.resize(positionFrac.size() + 1);
                    positionFrac.back().push_back(rTransformedFractional);
                }

            } // for (operationIdx = 0;...

            // validate - all groups have the same size
            //
            bool error=false;

            error = pointGroups.empty();
            if (!error)
            {
                nGroups = pointGroups.size();
                nMembers = pointGroups[0].size();
                for (groupIdx = 0; groupIdx < nGroups; groupIdx++)
                    if (pointGroups[groupIdx].size() != nMembers)
                        error = true;
                if (nGroups * nMembers != crystal.spaceGroup.nSymmetryOperations())
                    error = true;
            }

            if (error)
                on_error::throwException("problem when checking atom position symmetry, either bug or invalid space group specification", __FILE__, __LINE__);


            return true;
        } // bool findAtomSymmetry

        void symmetryOperationCartesian(
            const SpaceGroupOperation &symmOp,
            const UnitCell & unitCell,
            Matrix3d &rotationCartesian,
            Vector3d &translationCartesian)
        {
            Matrix3d mf,mc2f;
            Vector3d vf;
            symmOp.get(mf, vf);
            rotationCartesian = unitCell.getFractionalToCartesianMatrix() * mf *
                                unitCell.getCartesianToFractionalMatrix();
            unitCell.fractionalToCartesian(vf, translationCartesian);
        }


        bool equivalentPositions(
            const Vector3d& v1,
            const Vector3d& v2,
            const UnitCell& unitCell,
            Vector3i& shift)
        {
            int i;
            double threshold = 0.007;
            Vector3d diff = v1 - v2;
            Vector3d diffAfterShiftFract, diffAfterShiftCart;

            for (i = 0; i < 3; i++)
                shift(i) = math_utilities::roundInt(diff(i));

            for (i = 0; i < 3; i++)
                diffAfterShiftFract(i) = diff(i) - shift(i);

            //mCrystal.unitCell.fractionalToCartesian(diffAfterShiftFract, diffAfterShiftCart);
            unitCell.fractionalToCartesian(diffAfterShiftFract, diffAfterShiftCart);

            if (sqrt(diffAfterShiftCart * diffAfterShiftCart) > threshold)
                return false;

            return true;
        }


        double s12AdpSimilarityIndex(
            const UnitCell& uc,
            const std::vector<double>& cifAdp1,
            const std::vector<double>& cifAdp2)
        {
            vector<double> cartAdp1(6), cartAdp2(6);
            StructuralParametersConverter converter(uc);
            converter.convertADP(
                cifAdp1, cartAdp1,
                structural_parameters_convention::AdpConvention::U_cif,
                structural_parameters_convention::AdpConvention::U_cart);

            converter.convertADP(
                cifAdp2, cartAdp2,
                structural_parameters_convention::AdpConvention::U_cif,
                structural_parameters_convention::AdpConvention::U_cart);
            
            return s12AdpSimilarityIndex(cartAdp1, cartAdp2);
        }

        double s12AdpSimilarityIndex(
            const std::vector<double>& cartAdp1,
            const std::vector<double>& cartAdp2)
        {
            Matrix3d u1, u2, u1_inv, u2_inv;

            u1.set(cartAdp1[0], cartAdp1[3], cartAdp1[4],
                   cartAdp1[3], cartAdp1[1], cartAdp1[5],
                   cartAdp1[4], cartAdp1[5], cartAdp1[2]);

            u2.set(cartAdp2[0], cartAdp2[3], cartAdp2[4],
                   cartAdp2[3], cartAdp2[1], cartAdp2[5],
                   cartAdp2[4], cartAdp2[5], cartAdp2[2]);

            algebra3d::invert3d(u1, u1_inv);
            algebra3d::invert3d(u2, u2_inv);
            double det1 = algebra3d::det3d(u1_inv * u2_inv);
            double det2 = algebra3d::det3d(u1_inv + u2_inv);
            double r12 = sqrt(8.0 * sqrt(det1) / det2);
            return 100 * (1.0 - r12);
        }



        double msdCorrelation_2(
            const std::vector<double>& u1,
            const std::vector<double>& u2)
        {
            Matrix3d m1(u1[0], u1[3], u1[4],
                        u1[3], u1[1], u1[5],
                        u1[4], u1[5], u1[2]);

            Matrix3d m2(u2[0], u2[3], u2[4],
                        u2[3], u2[1], u2[5],
                        u2[4], u2[5], u2[2]);

            Matrix3d m12 = m1 * m2;
            Matrix3d m11 = m1 * m1;
            Matrix3d m22 = m2 * m2;

            double cov = algebra3d::trace(m12) - (1.0 / 3.0) * algebra3d::trace(m1) * algebra3d::trace(m2);
            double var1 = algebra3d::trace(m11) - (1.0 / 3.0) * algebra3d::trace(m1) * algebra3d::trace(m1);
            double var2 = algebra3d::trace(m22) - (1.0 / 3.0) * algebra3d::trace(m2) * algebra3d::trace(m2);

            double result = cov / sqrt(var1 * var2);

            return result;
        }


        double msdCorrelation(
            const std::vector<double>& u1,
            const std::vector<double>& u2)
        {
            Vector3d d1(u1[0], u1[1], u1[2]),
                     d2(u2[0], u2[1], u2[2]),
                     o1(u1[3], u1[4], u1[5]),
                     o2(u2[3], u2[4], u2[5]);
            double u_eq_1 = (u1[0] + u1[1] + u1[2]) / 3.0;
            double u_eq_2 = (u2[0] + u2[1] + u2[2]) / 3.0;

            double cov = (2.0 / 15.0) * (d1 * d2 + 2.0 * o1 * o2 - 3 * u_eq_1 * u_eq_2);
            double var1 = (2.0 / 15.0) * (d1 * d1 + 2.0 * o1 * o1 - 3 * u_eq_1 * u_eq_1);
            double var2 = (2.0 / 15.0) * (d2 * d2 + 2.0 * o2 * o2 - 3 * u_eq_2 * u_eq_2);

            return cov / sqrt(var1 * var2);
        }
        
        // old implementation
        
        //double msdCorrelation(
        //    const std::vector<double>& u1, 
        //    const std::vector<double>& u2)
        //{

        //    Matrix3d u1_m(u1[0], u1[3], u1[4],
        //        u1[3], u1[1], u1[5],
        //        u1[4], u1[5], u1[2]);

        //    Matrix3d u2_m(u2[0], u2[3], u2[4],
        //        u2[3], u2[1], u2[5],
        //        u2[4], u2[5], u2[2]);

        //    vector<Vector3d> r;
        //    vector<double> w;
        //    lebedev_laikov::get_grid(5810, r, w);

        //    double msd_var_1 = 0;
        //    double msd_var_2 = 0;
        //    double msd_cov = 0;
        //    double msd_av_1 = (u1[0] + u1[1] + u1[2]) / 3;
        //    double msd_av_2 = (u2[0] + u2[1] + u2[2]) / 3;

        //    double e_msd2 = 0;

        //    for (int i = 0; i < r.size(); i++)
        //    {
        //        double msd1 = r[i] * u1_m * r[i];
        //        double msd2 = r[i] * u2_m * r[i];

        //        e_msd2 += w[i] * msd1 * msd1;

        //        double dMsd = msd1 - msd2;

        //        //double absDiffMsd = fabs(dMsd);
        //        //double rootAbsDiffMsd = sqrt(absDiffMsd);

        //        msd_var_1 += w[i] * (msd1 - msd_av_1) * (msd1 - msd_av_1);
        //        msd_var_2 += w[i] * (msd2 - msd_av_2) * (msd2 - msd_av_2);
        //        msd_cov += w[i] * (msd1 - msd_av_1) * (msd2 - msd_av_2);

        //    }

        //    double result = msd_cov / sqrt(msd_var_1 * msd_var_2);
        //    
        //    return result;
        //}


        double rmsdCorrelation(
            const std::vector<double>& u1,
            const std::vector<double>& u2)
        {
            

            Matrix3d u1_m(u1[0], u1[3], u1[4],
                u1[3], u1[1], u1[5],
                u1[4], u1[5], u1[2]);

            Matrix3d u2_m(u2[0], u2[3], u2[4],
                u2[3], u2[1], u2[5],
                u2[4], u2[5], u2[2]);

            vector<Vector3d> r;
            vector<double> w;
            const int integrationGridSize = 5810;
            lebedev_laikov::get_grid(integrationGridSize, r, w);

            double rmsd_var_1 = 0;
            double rmsd_var_2 = 0;
            double rmsd_cov = 0;
            double rmsd_av_1 = 0;
            double rmsd_av_2 = 0;

            for (int i = 0; i < integrationGridSize; i++)
            {
                double rmsd1 = sqrt(r[i] * u1_m * r[i]);
                double rmsd2 = sqrt(r[i] * u2_m * r[i]);
                rmsd_av_1 += w[i] * rmsd1;
                rmsd_av_2 += w[i] * rmsd2;
            }


            for (int i = 0; i < integrationGridSize; i++)
            {
                double rmsd1 = sqrt(r[i] * u1_m * r[i]);
                double rmsd2 = sqrt(r[i] * u2_m * r[i]);

                double dRmsd = rmsd1 - rmsd2;

                rmsd_var_1 += w[i] * (rmsd1 - rmsd_av_1) * (rmsd1 - rmsd_av_1);
                rmsd_var_2 += w[i] * (rmsd2 - rmsd_av_2) * (rmsd2 - rmsd_av_2);
                rmsd_cov += w[i] * (rmsd1 - rmsd_av_1) * (rmsd2 - rmsd_av_2);

            }

            double result = rmsd_cov / sqrt(rmsd_var_1 * rmsd_var_2);

            return result;

        }


        double overlappingAdpSimilarityIndex(
            const vector<double>& u1,
            const vector<double>& u2)
        {
            //double norm1, norm2;
            Matrix3d u1_m(u1[0], u1[3], u1[4],
                u1[3], u1[1], u1[5],
                u1[4], u1[5], u1[2]);
            Matrix3d u2_m(u2[0], u2[3], u2[4],
                u2[3], u2[1], u2[5],
                u2[4], u2[5], u2[2]);
            double det1 = algebra3d::det3d(u1_m);
            double det2 = algebra3d::det3d(u2_m);
            double n1 = sqrt(1.0 / (algebra3d::det3d(u1_m) * 8 * M_PI * M_PI * M_PI));
            double n2 = sqrt(1.0 / (algebra3d::det3d(u2_m) * 8 * M_PI * M_PI * M_PI));

            Matrix3d m1, m2;
            algebra3d::invert3d(u1_m, m1);
            algebra3d::invert3d(u2_m, m2);

            vector<Vector3d> r;
            vector<double> w;
            lebedev_laikov::get_grid(5810, r, w);
            double d = 0;
            double e = 0;
            for (int i = 0; i < r.size(); i++)
            {
                double exp1, exp2;
                exp1 = r[i] * m1 * r[i] / 2.0;
                exp2 = r[i] * m2 * r[i] / 2.0;
                double radialIntegral;
                if ((fabs(1.0 - exp1 / exp2) < 1e-10)) // roughly exp1 = exp2 
                    radialIntegral = min(n1, n2) * sqrt(M_PI) / (2 * exp1 * sqrt(exp1));
                else
                {
                    // cp2 - square of cross point
                    double cp2 = log(n1 / n2) / (exp1 - exp2);
                    if (cp2 < 0) // there is no cross point, one function is always larger
                    {
                        if (n1 > n2)
                            radialIntegral = n2 * sqrt(M_PI) / (2 * exp2 * sqrt(exp2));
                        else
                            radialIntegral = n1 * sqrt(M_PI) / (2 * exp1 * sqrt(exp1));
                    }
                    else
                    {
                        if (n1 == n2)
                        {
                            double e = max(exp1, exp2);
                            radialIntegral = n1 * sqrt(M_PI) / (2 * e * sqrt(e));
                        }
                        else
                        {
                            // cp - cross point
                            double cp = sqrt(cp2);
                            // multiplier and exponent for the function which has lower value between -cp and cp
                            // i.e. the 'interior' part of the integral
                            double exp_int, n_int;
                            // the same for 'exterior'
                            double exp_ext, n_ext;
                            if (n1 > n2) // n1*exp(-
                            {
                                exp_int = exp2;
                                n_int = n2;
                                exp_ext = exp1;
                                n_ext = n1;
                            }
                            else
                            {
                                exp_int = exp1;
                                n_int = n1;
                                exp_ext = exp2;
                                n_ext = n2;
                            }

                            double exp_pow_3by2 = exp_int * sqrt(exp_int);
                            double d = sqrt(exp_int) * cp;
                            radialIntegral = -n_int * (2 * exp(-d * d) * d - sqrt(M_PI) * erf(d)) / (2 * exp_pow_3by2);
                            //radialIntegral = -(2 * exp(-a * cp * cp) * sqrt(a) * cp - sqrt(M_PI) * erf(sqrt(a) * cp)) / (2 * a * sqrt(a));
                            d = sqrt(exp_ext) * cp;
                            exp_pow_3by2 = exp_ext * sqrt(exp_ext);
                            radialIntegral += 2 * n_ext * (-sqrt(M_PI) * erf(d) / (4 * exp_pow_3by2) + cp * exp(-d * d) / (2 * exp_ext) + sqrt(M_PI) / (4 * exp_pow_3by2));
                        }
                    }
                }
                e += w[i] * radialIntegral;
            }
            double eta = e * 2 * M_PI;
            return eta;// 100 * (1.0 - eta);
        }


        bool isSuperLatticeNode(
            const Vector3i& node,
            const Vector3i& superLatticeNode,
            const std::vector<Vector3i>& superLatticeBase)
        {
            if (superLatticeBase.empty())
                return node == superLatticeNode;
            if (node == superLatticeNode)
                return true;
            Vector3i diff = node - superLatticeNode;

            if (superLatticeBase.size() == 1)
            {
                double aLength = sqrt(superLatticeBase[0] * superLatticeBase[0]);
                double wLength = sqrt(diff * diff);

                // is w collinear with l?

                double cos_wl = diff * superLatticeBase[0] / (aLength * wLength);
                if (fabs(1.0 - fabs(cos_wl)) > 1e-10)
                    return false;

                // is wLength multiple of lLength
                double abs_x = wLength / aLength;

                return (fabs(abs_x - round(abs_x)) < 1e-10);
            }

            // get the third lattice vector perpendicular to a and b
            Vector3d a = superLatticeBase[0];
            Vector3d b = superLatticeBase[1];
            Vector3d c;

            superLatticeBase.size() == 3 ? c = superLatticeBase[2] : c = cross_product(a, b);

            // localToGlobal*x transforms x to from local to global coordinates
            Matrix3d localToGlobal(a[0], b[0], c[0],
                a[1], b[1], c[1],
                a[2], b[2], c[2]);

            Matrix3d globalToLocal;
            algebra3d::invert3d(localToGlobal, globalToLocal);

            Vector3d wLocal = globalToLocal * Vector3d(diff);

            // in 2D case is it on the plane span by the lattice vectors?
            if (superLatticeBase.size() == 2)
                if (fabs(wLocal[2]) > 1e-10)
                    return false;

            // are new coordinates integers?
            double outOfInt = fabs(wLocal[0] - round(wLocal[0])) + fabs(wLocal[1] - round(wLocal[1])) +
                fabs(wLocal[2] - round(wLocal[2]));
            return (outOfInt < 1e-10);


        }

        void generate_hkl(
            const UnitCell& uc,
            double resolution,
            std::vector<Vector3i>& hkl)
        {
            hkl.clear();
            ReciprocalLatticeUnitCell rluCell(uc);
            Vector3i minVector, maxVector;
            // find bounding box in hkl space
            for(int i=0;i<8;i++)
            {
                Vector3d cubeVertex, cubeVertexFractional;
                cubeVertex[0] = (i & 1) ? 1.0 : -1.0;
                cubeVertex[1] = (i & 2) ? 1.0 : -1.0;
                cubeVertex[2] = (i & 4) ? 1.0 : -1.0;
                rluCell.cartesianToFractional(cubeVertex, cubeVertexFractional);
                cubeVertexFractional *= resolution;

                for(int j=0;j<3;j++)
                {
                    int intUp = ceil(cubeVertexFractional[j]);
                    int intDown = floor(cubeVertexFractional[j]);
                    if (intDown < minVector[j])
                        minVector[j] = intDown;
                    if (intUp > maxVector[j])
                        maxVector[j] = intUp;
                }
            }
            // generate hkl points in the bounding box and filter them by resolution
            Vector3d hklVectorCartesian;
            for (int h = minVector[0]; h <= maxVector[0]; h++)
                for (int k = minVector[1]; k <= maxVector[1]; k++)
                    for (int l = minVector[2]; l <= maxVector[2]; l++)
                    {
                        Vector3i hklVector(h, k, l);
                        rluCell.fractionalToCartesian(Vector3d(h, k, l) / resolution, hklVectorCartesian);
                        double dSpacing = 1.0 / sqrt(hklVectorCartesian * hklVectorCartesian);
                        if (dSpacing >= resolution)
                            hkl.push_back(hklVector);
                    }

        }

        void filter_hkl(
            const UnitCell& uc,
            double maxResolution,
            const std::vector<Vector3i>& hkl,
            std::vector<Vector3i>& filteredHkl)
        {
            vector<int> indicesOfOriginals;
            filter_hkl(uc, maxResolution, hkl, filteredHkl, indicesOfOriginals);
        }

        void filter_hkl(
            const UnitCell& uc,
            double maxResolution,
            const std::vector<Vector3i>& hkl,
            std::vector<Vector3i>& filteredHkl,
            std::vector<int>& indicesOfOriginals)
        {
            filteredHkl.clear();
            indicesOfOriginals.clear();
            ReciprocalLatticeUnitCell rluCell(uc);
            Vector3d hklVectorCartesian;
            
            for (int i=0; i<hkl.size(); i++)
            {
                rluCell.fractionalToCartesian(Vector3d(hkl[i]), hklVectorCartesian);
                double dSpacing = 1.0 / sqrt(hklVectorCartesian * hklVectorCartesian);
                if (dSpacing >= maxResolution)
                {
                    filteredHkl.push_back(hkl[i]);
                    indicesOfOriginals.push_back(i);
                }
            }

        }

    } //namespace crystal_structure_utilities

} // namespace discamb


