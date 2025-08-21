from HirshfeldAtomModelSettings import HirshfeldAtomModelSettings
from CrystalStructurePresenter import CrystalStructurePresenter
from CappingHydrogenAcceptDialg import CappingHydrogenAcceptDialg
from HarSettingsBox import HarSettingsBox
from Settings import Settings

from PySide6 import QtCore
from PySide6 import QtWidgets
import json
import copy
import os
import subprocess
import asyncio

class HirshfeldAtomModelPresenter:
    def __init__(self, har_settings_box: HarSettingsBox, structure_presenter: CrystalStructurePresenter):
        self.settings = HirshfeldAtomModelSettings()
        self.structure = None
        self.structure_presenter = structure_presenter
        self.n_threads = 1
        self.memory = 2
        self.memory_unit = "GB"
        self.basis_set = "def2-SVP"
        self.qm_method = "PBE"
        self.use_distributed_multipoles = False
        self.distributed_multipoles_range = 8.0
        self.qm_program = "orca"
        self.har_settings_box = har_settings_box
        self.last_selected_chemical_unit = -1
        self.use_xray = True
        self.har_settings_box.har_x_ray.connect(self.onRadiationTypeChanged)
        self.har_settings_box.memory_units_changed.connect(self.set_memory_unit)
        self.har_settings_box.qm_method_changed.connect(self.set_qm_method)
        self.har_settings_box.basis_set_changed.connect(self.set_basis_set)
        self.har_settings_box.n_threads_changed.connect(self.set_n_threads)
        self.har_settings_box.qm_software_changed.connect(self.set_qm_program)
        self.har_settings_box.memory_changed.connect(self.set_memory)
        self.har_settings_box.distributed_multipole_range.connect(self.set_distributed_multipoles_range)
        self.har_settings_box.use_distributed_multipoles.connect(self.set_use_distributed_multipoles)

        self.har_settings_box.chemicalUnitBox.selectedChemicalUnitIndexChanged.connect(self.showChemicalUnitDetails)
        self.har_settings_box.chemicalUnitBox.colorByChemicalUnits.connect(self.colorByChemicalUnitsChanged)
        self.har_settings_box.chemicalUnitBox.mergeChemicalUnits.connect(self.mergeChemicalUnits)
        self.har_settings_box.chemicalUnitBox.deleteChemicalUnits.connect(self.deleteChemicalUnits)
        self.har_settings_box.chemicalUnitBox.addChemicalUnit.connect(self.addChemicalUnit)
        self.har_settings_box.chemicalUnitBox.selectTableSelectedSubsystemAtoms.connect(self.select_subsystem_atoms)
        self.har_settings_box.chemicalUnitBox.splitChemicalUnit.connect(self.splitChemicalUnit)
        self.har_settings_box.chemicalUnitBox.showOnlyChemicalUnits.connect(self.showOnlyChemicalUnits)
        self.har_settings_box.chemicalUnitBox.showAtomWeights.connect(self.show_atom_weights)
        self.har_settings_box.chemicalUnitBox.setAtomWeights.connect(self.set_atom_weights)
        self.har_settings_box.chemicalUnitBox.setMultipoleSites.connect(self.set_mutlipole_sites)
        self.har_settings_box.chemicalUnitBox.showMultipoleSites.connect(self.show_distributed_multipoles_sites)

        self.colour_by_chemical_units_checked = (self.har_settings_box.chemicalUnitBox.colorByChemicalUnitCheckBox.
                                                 isChecked())
        self.har_settings_box.chemicalUnitBox.chemicalUnitChargeChanged.connect(self.on_chemical_unit_charge_changed)
        self.har_settings_box.chemicalUnitBox.chemicalUnitSpinMultiplicityChanged.connect(
            self.on_chemical_unit_spin_multiplicity_changed)
        self.har_settings_box.chemicalUnitBox.atomWeightChanged.connect(self.onAtomWeightChanged)

    def onRadiationTypeChanged(self, x_ray):
        self.use_xray = x_ray

    def set_mutlipole_sites(self, subsystem_idx):
        subsystem = self.settings.subsystemInfo.subsystems[subsystem_idx]
        subsystem.distributed_multipole_sites = self.structure_presenter.getSelectedAtoms()
        subsystem.has_user_defined_distributed_multipole_sites = True
    def set_use_distributed_multipoles(self, use_multipoles):
        self.use_distributed_multipoles = use_multipoles

    def set_distributed_multipoles_range(self, range):
        self.distributed_multipoles_range = range

    def onAtomWeightChanged(self, atom_idx, weight):
        self.settings.subsystemInfo.subsystems[self.last_selected_chemical_unit].weights[atom_idx] = weight
        self.settings.subsystemInfo.subsystems[self.last_selected_chemical_unit].user_defined_weight[atom_idx] = True

    def show_atom_weights(self, subsystem_idx):
        if subsystem_idx<0:
            self.structure_presenter.molCanvas.unsetLabels()
            return
        atom_id_list = self.settings.subsystemInfo.subsystems[subsystem_idx].atoms
        labels = []
        for weight in self.settings.subsystemInfo.subsystems[subsystem_idx].weights:
            labels.append(f'{weight:.3f}')
        self.structure_presenter.show_listed_atom_labels(atom_id_list, labels)

    def show_distributed_multipoles_sites(self, subsystem_idx):
        if subsystem_idx<0:
            return
        self.structure_presenter.deselect_all()
        subsystem = self.settings.subsystemInfo.subsystems[subsystem_idx]
        atoms_to_show = copy.deepcopy(subsystem.atoms)
        atoms_to_select = []
        if self.use_distributed_multipoles:
            if subsystem.has_user_defined_distributed_multipole_sites:
                atoms_to_show += subsystem.distributed_multipole_sites
                atoms_to_select = subsystem.distributed_multipole_sites
            else:
                atoms_to_select= self.structure_presenter.structureData.getNeighboringAtoms(
                    subsystem.atoms, self.distributed_multipoles_range, True)
                atoms_to_show += atoms_to_select
        self.structure_presenter.set_visible_atoms(atoms_to_show, save_selection=False)
        self.structure_presenter.select_atoms(atoms_to_select)

    def set_atom_weights(self, subsystem_idx, weight):
        subsystem = self.settings.subsystemInfo.subsystems[subsystem_idx]
        selected_atoms = self.structure_presenter.getSelectedAtoms()
        if selected_atoms:
            for atom in selected_atoms:
                if atom in subsystem.atoms:
                    idx = subsystem.atoms.index(atom)
                    subsystem.weights[idx] = weight
                    subsystem.user_defined_weight[idx] = True
        else:
            subsystem.weights = [weight]*len(subsystem.atoms)
            subsystem.user_defined_weight = [True]*len(subsystem.atoms)
        if self.har_settings_box.chemicalUnitBox.chemicalUnitAtomsView.isVisible():
            self.showChemicalUnitDetails(subsystem_idx)

    def select_subsystem_atoms(self, subsystem_idx, atom_indices):
        atom_list = []
        for idx in atom_indices:
            atom_list.append(self.settings.subsystemInfo.subsystems[subsystem_idx].atoms[idx])
        self.structure_presenter.select_atoms(atom_list)

    def set_memory_unit(self, memory_unit: str):
        self.memory_unit = memory_unit

    def set_basis_set(self, basis_set: str):
        self.basis_set = basis_set

    def set_qm_method(self, qm_method: str):
        self.qm_method = qm_method

    def set_qm_program(self, qm_program: str):
        self.qm_program = qm_program

    def set_n_threads(self, n):
        self.n_threads = n

    def set_memory(self, n):
        self.memory = n
    def on_chemical_unit_charge_changed(self, idx: int, charge: int):
        self.settings.subsystemInfo.subsystems[idx].charge = charge

    def on_chemical_unit_spin_multiplicity_changed(self, idx: int, spin: int):
        self.settings.subsystemInfo.subsystems[idx].spin_multiplicity = spin

    def run(self, directory, structure_file, hkl_file):
        if not self.save(directory):
            return
        subprocess.run(["./bin/discamb_ff.exe", structure_file, hkl_file], cwd=directory)
        msgBox = QtWidgets.QMessageBox()
        msgBox.setText("done")
        msgBox.exec()

    def save(self, directory):
        settings = Settings()

        if self.qm_program == "orca":
            orca_path = settings.get('paths/orca')
            if not orca_path:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setText("Set orca folder in Edit/Paths")
                msgBox.exec()
                return False
        if self.qm_program == "gaussian":
            gaussian_path = settings.get('paths/gaussian')
            if not gaussian_path:
                msgBox = QtWidgets.QMessageBox()
                msgBox.setText("Set gaussian folder in Edit/Paths")
                return False
        model = {"model": "HAR"}
        model["electron scattering"] = not self.use_xray
        partition_name = self.har_settings_box.partition_box.partition_combo_box.currentText()
        if partition_name == "exponential Hirshfeld":
            model["power"] = self.har_settings_box.partition_box.expHarParameters_doubleSpinBox.value()
        '''
        "qm method": "B3LYP",
        "model": "har",
        "basis set": "cc-pVTZ",
        "qm program": "orca",
        "qm folder": "c:\\Orca5",
        "n cores": 12,
        "memory": "40GB",
        '''
        model["qm program"] = self.qm_program
        model["qm method"] = self.qm_method
        model["basis set"] = self.basis_set
        model["n cores"] = self.n_threads
        model["qm folder"] = settings.get("paths/"+self.qm_program)
        model["memory"] = str(self.memory) + self.memory_unit
        if self.use_distributed_multipoles:
            model["multipoles"] = {
                "threshold": self.distributed_multipoles_range
            }
        '''
        susbsystems_with_user_defined_multipoles = []
        for subsystem in self.settings.subsystemInfo.subsystems:
            if subsystem.distributed_multipole_sites:
                susbsystems_with_user_defined_multipoles.append(subsystem)
        if susbsystems_with_user_defined_multipoles:
            #"qm settings per subsystem"
            #"qm settings per subsystem"
            susbystem_multipoles_dict = {}
            for subsystem in susbsystems_with_user_defined_multipoles:
                sites_str = ""
                for site in subsystem.distributed_multipole_sites:
                    site_label = self.structure.getAtomLabel(site[0])
                    sites_str += " " + site_label
                    site_symm_op = self.structure.getSymmetryOperationStr(site)
                    if site_symm_op.lower() != "x,y,z":
                        sites_str += "," + site_symm_op
                susbystem_multipoles_dict[subsystem.label] = {
                    "multipoles": {"sites": sites_str}
                }
            model["qm settings per subsystem"] = susbystem_multipoles_dict
        '''
        cwd = os.getcwd()
        model["partition settings"] = {"atoms file": os.path.join(cwd, "data", "atomic_densities.txt")}
        user_defined_representatives = False
        for subsystem in self.settings.subsystemInfo.subsystems:
            if True in subsystem.user_defined_weight:
                user_defined_representatives = True
        qm_structure = []
        atom_representatives = []
        for subsystem in self.settings.subsystemInfo.subsystems:
            atoms = ""
            for atom in subsystem.atoms:
                atom_label = self.structure.getAtomLabel(atom[0])
                atoms += " " + atom_label
                atom_symm_op = self.structure.getSymmetryOperationStr(atom)
                if atom_symm_op.lower() != "x,y,z":
                    atoms += ","+atom_symm_op
            # user defined distributed multipole sites
            sites_str = ""
            if subsystem.distributed_multipole_sites:
                for site in subsystem.distributed_multipole_sites:
                    site_label = self.structure.getAtomLabel(site[0])
                    sites_str += " " + site_label
                    site_symm_op = self.structure.getSymmetryOperationStr(site)
                    if site_symm_op.lower() != "x,y,z":
                        sites_str += "," + site_symm_op

            capping_hydrogen_atoms = []
            for capping_h in subsystem.capping_hydrogen_atoms:
                label_bonded = self.structure.getAtomLabel(capping_h[0][0])
                symm_op_bonded = self.structure.getSymmetryOperationStr(capping_h[0])
                label_directing = self.structure.getAtomLabel(capping_h[1][0])
                symm_op_directing = self.structure.getSymmetryOperationStr(capping_h[1])
                cap_h_dict = {
                    "bonded": label_bonded,
                    "bonded symm": symm_op_bonded,
                    "directing": label_directing,
                    "directing symm": symm_op_directing
                }
                capping_hydrogen_atoms.append(cap_h_dict)
            subsystem_dict = {
                "label": subsystem.label,
                "charge": subsystem.charge,
                "spin multiplicity": subsystem.spin_multiplicity,
                "atoms": atoms}
            if subsystem.distributed_multipole_sites:
                subsystem_dict["multipole sites"] = sites_str
            if capping_hydrogen_atoms:
                subsystem_dict["capping hydrogen atoms"] = capping_hydrogen_atoms
            representatives_dict = {
                "label": subsystem.label,
                "atom weights": subsystem.weights
                }
            qm_structure.append(subsystem_dict)
            atom_representatives.append(representatives_dict)
        model["qm structure"] = qm_structure
        if user_defined_representatives:
            model["atom representatives"] = atom_representatives
        with open(os.path.join(directory, "aspher.json"),'w') as out:
            out.write(json.dumps(model, indent=4))
        return True

    def showOnlyChemicalUnits(self, chemicalUnits):
        if not chemicalUnits:
            chemicalUnits = list(range(0, len(self.settings.subsystemInfo.subsystems)))
        #if not chemicalUnits:
        #    return
        atoms = []
        for idx in chemicalUnits:
            atoms += self.settings.subsystemInfo.subsystems[idx].atoms
        self.structure_presenter.set_visible_atoms(atoms)
    
    def splitChemicalUnit(self, idx):
        self.settings.subsystemInfo.splitChemicalUnit(idx)
        self.har_settings_box.setChemicalUnits(self.settings.subsystemInfo.subsystems)
        if self.colour_by_chemical_units_checked:
            self.colorByChemicalUnit()

    def find_broken_bonds(self, atoms):
        """Finds atoms which do not have all bonded neighbours
        included in the atoms set and the missing neighbours

        :param atoms:
        :return:
        """
        includingMolecule = self.structure.getNeighboringAtoms(atoms, 0, True)
        missing_neighbours = {}
        for i in range(0, len(atoms)):
            atom = atoms[i]
            neighbours = self.structure.getBondedAtoms(atom)
            has_missing_neighbours = False
            for neighbour in neighbours:
                if neighbour not in atoms:
                    if not has_missing_neighbours:
                        has_missing_neighbours = True
                        missing_neighbours[i]=[neighbour]
                    else:
                        missing_neighbours[i].append(neighbour)

        return missing_neighbours

    def addChemicalUnit(self):
        chemicalUnitAtoms = self.structure_presenter.getSelectedAtoms()
        missing_neighbours = self.find_broken_bonds(chemicalUnitAtoms)
        capping_hydrogen_indices = []
        if missing_neighbours:
            capping_h_dialog = CappingHydrogenAcceptDialg(self.structure, chemicalUnitAtoms, missing_neighbours)
            capping_h_dialog.exec()
            capping_hydrogen_indices = capping_h_dialog.get_values()
        self.settings.subsystemInfo.addSubsystem(chemicalUnitAtoms, capping_hydrogens=capping_hydrogen_indices)
        self.har_settings_box.setChemicalUnits(self.settings.subsystemInfo.subsystems)
        if self.colour_by_chemical_units_checked:
            self.colorByChemicalUnit()

    def mergeChemicalUnits(self, chemicalUnitsIndices):
        selectedAtoms = self.structure_presenter.getSelectedAtoms()
        self.settings.subsystemInfo.mergeChemicalUnits(chemicalUnitsIndices, selectedAtoms)
        self.har_settings_box.setChemicalUnits(self.settings.subsystemInfo.subsystems)
        if self.colour_by_chemical_units_checked:
            self.colorByChemicalUnit()

    def deleteChemicalUnits(self, chemicalUnitsIndices):
        #selectedAtoms = self.structure_presenter.getSelectedAtoms()
        #self.settings.subsystemInfo.deleteChemicalUnits(chemicalUnitsIndices, selectedAtoms)
        self.settings.subsystemInfo.deleteChemicalUnits(chemicalUnitsIndices)
        self.har_settings_box.setChemicalUnits(self.settings.subsystemInfo.subsystems)
        if self.colour_by_chemical_units_checked:
            self.colorByChemicalUnit()

        
    def colorByChemicalElement(self):
        self.structure_presenter.defaultAtomColors()
    
    def colorByChemicalUnitsChanged(self, selected):
        if not selected:
            self.colour_by_chemical_units_checked = False
            self.colorByChemicalElement()
        else:
            self.colour_by_chemical_units_checked = True
            self.colorByChemicalUnit()
#            atoms = []
#            colors = []
#            nChemicalUnits = len(self.settings.subsystemInfo.subsystems)
#            if nChemicalUnits == 0:
#                return
#            for chemicalUnitIdx in range(0, nChemicalUnits):
#                atoms += self.settings.subsystemInfo.subsystems[chemicalUnitIdx].atoms
#                nAtoms = len(self.settings.subsystemInfo.subsystems[chemicalUnitIdx].atoms)
#                qtColor = self.harSettingsBox.chemicalUnitBox.chemicalUnitsModel.item(chemicalUnitIdx, 0).background().color()
#                color = [qtColor.red()/255.0, qtColor.green()/255.0, qtColor.blue()/255.0]
#                colors += nAtoms*[color]
#            self.structurePresenter.colorAtoms(atoms, colors)

    def colorByChemicalUnit(self):
        atoms = []
        colors = []
        nChemicalUnits = len(self.settings.subsystemInfo.subsystems)
        if nChemicalUnits == 0:
            return
        for chemicalUnitIdx in range(0, nChemicalUnits):
            atoms += self.settings.subsystemInfo.subsystems[chemicalUnitIdx].atoms
            nAtoms = len(self.settings.subsystemInfo.subsystems[chemicalUnitIdx].atoms)
            qtColor = self.har_settings_box.chemicalUnitBox.chemicalUnitsModel.item(chemicalUnitIdx, 0).background().color()
            color = [qtColor.red()/255.0, qtColor.green()/255.0, qtColor.blue()/255.0]
            colors += nAtoms*[color]
        self.structure_presenter.colorAtoms(atoms, colors)
    
    
    def setCrystalStructure(self, structure):
        self.settings.setStructure(structure)
        self.structure = structure
        self.har_settings_box.setChemicalUnits(self.settings.subsystemInfo.subsystems)

    def unsetCrystalStructure(self):
        self.settings.subsystemInfo.subsystems = []
        self.har_settings_box.setChemicalUnits([])
        
    
    def showChemicalUnitDetails(self, idx):
        self.last_selected_chemical_unit = idx
        chemicalUnit = self.settings.subsystemInfo.subsystems[idx]
        atomLabels = []
        atomSymmOps = []
        for atom in chemicalUnit.atoms:
            atomSymmOps.append(self.structure.getSymmetryOperationStr(atom))
            atomLabels.append(self.structure.getAtomLabel(atom[0]))
        self.har_settings_box.chemicalUnitBox.showChemicalUnitDetails(idx, atomLabels, atomSymmOps, chemicalUnit.weights)
        
