#ifndef _DISCAMB_CRYSTALSTRUCTURE_CRYSTAL_HPP_
#define _DISCAMB_CRYSTALSTRUCTURE_CRYSTAL_HPP_

#include "SpaceGroup.h"
#include "UnitCell.h"
#include "StructuralParametersConverter.h"

#include "discamb/MathUtilities/Vector3.h"

#include <vector>
#include <map>
#include <string>

namespace discamb {

    /**
    * \defgroup CrystalStructure CrystalStructure
    \brief Information on crystal structure.
    * @{
    */


/** \brief Container for data describing an atom in crystal. 

Some members do not have strictly defined meaning - e.g. structural parameters (AtomInCrystal::coordinates, AtomInCrystal::adp) 
do not have assigned units and AtomInCrystal::type does not have fixed meaning.
*/

struct AtomInCrystal
{
    /** \brief Position of the atom.
    (No fixed units assumed).*/
    Vector3d coordinates;
	Vector3d coordinates_sigma;
    Vector3d coordinates_precision;
    /** \brief Type of atom. 
    Does not have fixed meaning, can be e.g. type of scatterer.*/
    std::string type;
    /** \brief Number of symmetry equivalent atoms per unit cell (0 if undefined) */
    int multiplicity = 1;
    /**
    \brief Atomic displacement parameters.
    In the case of anisotropic ADPs the parameters are
    in the following order U11 U22 U33 U12 U13 U23
    */
    std::vector<double> adp;
	std::vector<double> adp_sigma;
    std::vector<double> adp_precision;
    
    /** \brief Atom label. */ 
    std::string label;
    /** \brief Occupancy.*/
    double occupancy = 1.0;
	double occupancy_sigma = 0.0;
    double occupancy_precision = 0.0;
    /** \brief Symmetry operations forming point group of the atom site.*/
    std::vector<SpaceGroupOperation> siteSymetry;
};

/** \brief Container for data describing crystal structure. */

struct Crystal
{
    /** \brief Unit cell. */
    UnitCell unitCell;
    /** \brief Space group. */
    SpaceGroup spaceGroup;
    /** \brief Atoms in asymmetric unit.*/
    std::vector<AtomInCrystal> atoms;
    /** \brief Specifies if the coordinates are in fractional or Cartesian coordinates. In the latter case Angstroms are assumed as units.*/
    structural_parameters_convention::XyzCoordinateSystem xyzCoordinateSystem = structural_parameters_convention::XyzCoordinateSystem::fractional;
    /** \brief Specifies in which way ADPs are parameterized (\f$U_{cif}\f$, \f$U_{cart}\f$ or \f$U^*\f$)*/
    structural_parameters_convention::AdpConvention adpConvention = structural_parameters_convention::AdpConvention::U_cif;
    int atomIdx(const std::string& atomLabel) const;
    bool atomIdx(const std::string& atomLabel,int &idx) const;
};

/**@}*/

}//namespace discamb{





#endif /*_DISCAMB_CRYSTALSTRUCTURE_CRYSTAL_HPP_*/
