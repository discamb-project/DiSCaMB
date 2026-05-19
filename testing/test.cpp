#include <discamb/Scattering/HansenCoppens_SF_Engine4.h>

#include <cereal/archives/binary.hpp>
#include <cereal/types/complex.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/vector.hpp>
#include <fstream>

#include "cereal/details/helpers.hpp"
#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/MathUtilities/Vector3.h"
#include "discamb/Scattering/Real.h"
#include "discamb/Scattering/SF_CalcDataTypes.h"
#include "discamb/Scattering/SF_Engine_DataTypes.h"
template <typename T>
void load_from_file(const std::string &filename, T &value) {
    std::ifstream is(filename, std::ios::binary);
    if (!is) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    cereal::BinaryInputArchive ar(is);
    ar(value);
}

int main() {
    try {
        // 1. Load cell
        discamb::UnitCell unit_cell;
        load_from_file("unit_cell.bin", unit_cell);

        std::vector<discamb::sf_engine_data_types::HC_WfnParam> wfn_parameters;
        load_from_file("wfn_parameters.bin", wfn_parameters);

        std::vector<discamb::sf_engine_data_types::HC_TypeParam>
            type_parameters;
        load_from_file("type_parameters.bin", type_parameters);

        std::vector<int> atom_to_wfn_map;
        load_from_file("atom_to_wfn_map.bin", atom_to_wfn_map);

        std::vector<int> atom_to_type_map;
        load_from_file("atom_to_type_map.bin", atom_to_type_map);

        std::vector<discamb::Vector3<discamb::REAL>> atomic_positions;
        load_from_file("atomic_positions.bin", atomic_positions);

        std::vector<std::vector<discamb::REAL>> atomic_displacement_parameters;
        load_from_file("atomic_displacement_parameters.bin",
                       atomic_displacement_parameters);

        std::vector<discamb::REAL> atomic_occupancy;
        load_from_file("atomic_occupancy.bin", atomic_occupancy);

        std::vector<std::complex<discamb::REAL>> anomalous_dispersion;
        load_from_file("anomalous_dispersion.bin", anomalous_dispersion);

        std::vector<discamb::REAL> atomic_multiplicity_factor;
        load_from_file("atomic_multiplicity_factor.bin",
                       atomic_multiplicity_factor);

        std::vector<discamb::Matrix3<discamb::REAL>> local_coordinate_systems;
        load_from_file("local_coordinate_systems.bin",
                       local_coordinate_systems);

        std::vector<discamb::sf_engine_data_types::SymmetryOperation>
            symmetry_operations;
        load_from_file("symmetry_operations.bin", symmetry_operations);

        bool centrosymmetric;
        load_from_file("centrosymmetric.bin", centrosymmetric);

        discamb::Vector3<discamb::REAL> inversion_translation;
        load_from_file("inversion_translation.bin", inversion_translation);

        std::vector<discamb::Vector3<discamb::REAL>> h_vectors;
        load_from_file("h_vectors.bin", h_vectors);

        std::vector<discamb::Vector3i> hkl_indices;
        load_from_file("hkl_indices.bin", hkl_indices);

        std::vector<std::complex<discamb::REAL>> f;
        load_from_file("f.bin", f);

        std::vector<discamb::TargetFunctionAtomicParamDerivatives>
            d_target_dparam;
        load_from_file("dtarget_dparam.bin", d_target_dparam);

        std::vector<std::complex<discamb::REAL>> d_target_df;
        load_from_file("dtarget_df.bin", d_target_df);

        std::vector<bool> include_atom_contribution;
        load_from_file("include_atom_contribution.bin",
                       include_atom_contribution);

        int n_threads;
        load_from_file("n_threads.bin", n_threads);

        discamb::DerivativesSelector derivatives_switch;
        load_from_file("derivatives_switch.bin", derivatives_switch);

        bool electron;
        load_from_file("electron.bin", electron);

        std::vector<int> atomic_number;
        load_from_file("atomic_number.bin", atomic_number);

        discamb::HansenCoppens_SF_Engine4 engine;
        engine.calculateSF(unit_cell,
                           wfn_parameters,
                           type_parameters,
                           atom_to_wfn_map,
                           atom_to_type_map,
                           atomic_positions,
                           atomic_displacement_parameters,
                           atomic_occupancy,
                           anomalous_dispersion,
                           atomic_multiplicity_factor,
                           local_coordinate_systems,
                           symmetry_operations,
                           centrosymmetric,
                           inversion_translation,
                           h_vectors,
                           hkl_indices,
                           f,
                           d_target_dparam,
                           d_target_df,
                           include_atom_contribution,
                           n_threads,
                           derivatives_switch,
                           electron,
                           atomic_number);

    } catch (const std::exception &e) {
        std::cerr << "Error during deserialization: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
