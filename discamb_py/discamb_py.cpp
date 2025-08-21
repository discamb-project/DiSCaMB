#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "CrystalStructure.h"

#include "har_utilities.h"
#include "taam_utilities.h"

#include "discamb/AtomTyping/CrystalAtomTypeAssigner.h"
#include "discamb/AtomTyping/LocalCoordinateSystemCalculator.h"

#include "discamb/BasicChemistry/periodic_table.h"

#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicUtilities/string_utilities.h"

#include "discamb/CrystalStructure/Crystal.h"
#include "discamb/CrystalStructure/crystal_structure_utilities.h"
#include "discamb/CrystalStructure/LocalCoordinateSystemInCrystal.h"

#include "discamb/IO/hkl_io.h"
#include "discamb/IO/structure_io.h"
#include "discamb/IO/MATTS_BankReader.h"

#include "discamb/Scattering/AnyIamCalculator.h"
#include "discamb/Scattering/AnyScattererStructureFactorCalculator.h"
#include "discamb/Scattering/ElectronFromXrayFormFactorCalculationsManager.h"
#include "discamb/Scattering/HcFormFactorCalculationsManager.h"


#include <vector>
#include <string>
//#include <iostream>
#include <iomanip>


void set_crystal(
    discamb::Crystal& crystal,
    const std::vector<std::string>& symm_ops,
    const std::vector<double>& unit_cell,
    const std::vector<std::string>& symbols,
    const std::vector<double>& positions_fractional,
    const std::vector<double>& occupancy,
    const std::vector<int>& adp_size,
    const std::vector<double>& adps,
    const std::vector<int>& multiplicity)
{
    
    std::vector<discamb::SpaceGroupOperation> space_group_operations;
    for (auto const& symm_op : symm_ops)
        space_group_operations.push_back(discamb::SpaceGroupOperation(symm_op));
    crystal.spaceGroup.set(space_group_operations);
    
    crystal.unitCell.set(unit_cell[0], unit_cell[1], unit_cell[2], unit_cell[3], unit_cell[4], unit_cell[5]);
    int nAtoms = symbols.size();
    crystal.atoms.resize(nAtoms);
    int adpVecPosition = 0;
    for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    {
        crystal.atoms[atomIdx].adp.resize(adp_size[atomIdx]);
        for (int i = 0; i < adp_size[atomIdx]; i++)
            crystal.atoms[atomIdx].adp[i] = adps[adpVecPosition + i];
        adpVecPosition += adp_size[atomIdx];
        
        crystal.atoms[atomIdx].label = symbols[atomIdx] + std::to_string(atomIdx + 1);
        crystal.atoms[atomIdx].coordinates.set(
            positions_fractional[3 * atomIdx],
            positions_fractional[3 * atomIdx + 1],
            positions_fractional[3 * atomIdx + 2]);
        
        crystal.atoms[atomIdx].multiplicity = multiplicity[atomIdx];
        crystal.atoms[atomIdx].occupancy = occupancy[atomIdx];
    }
    
}

void set_hkl(
    std::vector<discamb::Vector3i>& hkl,
    const std::vector<int>& miller_indices)
{
    //deserialize Miller indices
    int n_hkl = miller_indices.size() / 3;
    hkl.resize(n_hkl);
    for (int i = 0; i < n_hkl; i++)
        hkl[i].set(miller_indices[3 * i], miller_indices[3 * i + 1], miller_indices[3 * i + 2]);
}

void calculate_sf_iam(
    const discamb::Crystal& crystal,
    const std::vector<discamb::Vector3i>& hkl,
    std::vector< std::complex<double> > &structureFactors,
    bool electron_scattering)
{
    bool electronScattering = false;
    //table can be "Waasmeier-Kirfel", "IT92", "electron-IT"
    std::string table = "Waasmeier-Kirfel";
    if (electron_scattering)
        table = "electron-IT";
    discamb::AnyIamCalculator iamCalculator(crystal, electronScattering, table);
    std::vector<bool> countAtomContribution(crystal.atoms.size(), true);
    iamCalculator.calculateStructureFactors(crystal.atoms, hkl, structureFactors, countAtomContribution);
}



void pack_sf(
    std::vector<double>& sf_packed,
    std::vector< std::complex<double> >& structureFactors)
{
    int n_sf = structureFactors.size();
    sf_packed.resize(2*n_sf);
    
    for (int i = 0; i < n_sf; i++)
    {
        sf_packed[2 * i] = structureFactors[i].real();
        sf_packed[2 * i + 1] = structureFactors[i].imag();
    }
}


std::vector<double> calc_sf(
    const std::vector<std::string>& symm_ops,
    const std::vector<double>& unit_cell,
    const std::vector<std::string>& symbols,
    const std::vector<double>& positions_fractional,
    const std::vector<double>& occupancy,
    const std::vector<int>& adp_size,
    const std::vector<double>& adps,
    const std::vector<int> &multiplicity,
    const std::vector<int> &miller_indices,
    bool taam,
    bool electron_scattering)
{
    discamb::Crystal crystal;
    set_crystal(crystal, symm_ops, unit_cell,
        symbols, positions_fractional, occupancy,
        adp_size, adps, multiplicity);
    
    std::vector<discamb::Vector3i> hkl;
    set_hkl(hkl, miller_indices);

    // calculate structure factors
    std::vector< std::complex<double> > structureFactors;

    calculate_sf_iam(crystal, hkl, structureFactors, electron_scattering);


    // pack structure factors
    std::vector<double> result;
    pack_sf(result, structureFactors);

    return result;
}

double discamb_add(double a, double b)
{
    return a + b;
}


PYBIND11_MODULE(discamb_py, m) {
    m.doc() = "DiSCaMB bindings"; // optional module docstring

    pybind11::class_<har_utilities::Representative>(m, "Representatives")
        .def_readwrite("substructureIdx", &har_utilities::Representative::substructureIdx)
        .def_readwrite("atomIdx", &har_utilities::Representative::atomIdx)
        .def_readwrite("weight", &har_utilities::Representative::weight);

    m.def("calc_sf", &calc_sf, "calc sf");
    m.def("discamb_add", &discamb_add, "discamb_add");
    m.def("find_default_representatives", &har_utilities::find_default_representatives);
    m.def("find_atom_types",&taam_utilities::find_atom_types);
    m.def("find_atom_types_for_fragments", &taam_utilities::find_atom_types_for_fragments);
    m.def("get_atom_type_names", &taam_utilities::get_atom_type_names);
    pybind11::enum_<PackInlcudeMode>(m, "PackInlcudeMode")
        .value("ATOM", PackInlcudeMode::ATOM)
        .value("MOLECULE_ATOM", PackInlcudeMode::MOLECULE_ATOM)
        .value("MOLECULE_CENTER", PackInlcudeMode::MOLECULE_CENTER)
        .export_values();

    pybind11::class_<CrystalStructure>(m, "CrystalStructure")
        .def(pybind11::init<const std::string&>())
        .def("getAtomNames", &CrystalStructure::getAtomNames)
        .def("getPositions", &CrystalStructure::getPositions)
        .def("getOccupancies", &CrystalStructure::getOccupancies)
        .def("getAdps", &CrystalStructure::getAdps)
        .def("getAtomSymbols", &CrystalStructure::getAtomSymbols)
        .def("getBondedAtoms", &CrystalStructure::getBondedAtoms)
        .def("getBonds", pybind11::overload_cast<void> (&CrystalStructure::getBonds, pybind11::const_))
        .def("getBonds", pybind11::overload_cast<const std::vector<std::vector<int> > &> (&CrystalStructure::getBonds, pybind11::const_))
        //.def("getBonds2", &CrystalStructure::getBonds2)
        .def("getGrouppedAsymetricUnitBondedAtoms", &CrystalStructure::getGrouppedAsymetricUnitBondedAtoms)
        .def("getIncludedAtoms", &CrystalStructure::getIncludedAtoms)
        .def("atomicNumber", &CrystalStructure::atomicNumber)
        .def("getSymmetryOperationStr", &CrystalStructure::getSymmetryOperationStr)
        .def("getAtomLabel", &CrystalStructure::getAtomLabel)
        .def("getAsymetricUnitAtoms", &CrystalStructure::getAsymetricUnitAtoms)
        .def("getAtomPositionCart", &CrystalStructure::getAtomPositionCart)
        .def("getThermalEllipsoid", &CrystalStructure::getThermalEllipsoid)
        .def("getUnitCellVectors", &CrystalStructure::getUnitCellVectors)
        .def("getNeighboringAtoms",&CrystalStructure::getNeighboringAtoms)
        .def("splitIntoChemicalUnits", &CrystalStructure::splitIntoChemicalUnits)
        .def("getDistance",&CrystalStructure::getDistance)
        .def("groupIntoChemicalUnits",&CrystalStructure::groupIntoChemicalUnits)
        .def("numberOfAtoms", &CrystalStructure::numberOfAtoms);

}

