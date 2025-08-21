import os
import json
import subprocess
from Settings import Settings
from PySide6 import QtWidgets
from PySide6 import QtCore
from CrystalAsSubsystems import CrystalAsSubsystems
from CrystalAsSubsystems import Subsystem
import discamb_py
class TaamPresenter():
    def __init__(self, taamSettingsBox, structurePresenter):
        self.view = taamSettingsBox
        self.structure_presenter = structurePresenter
        self.radiation = "X-ray"
        self.subsystems = None
        self.structure = None
        self.atom_type_labels = None
        self.unassigned = None
        self.last_selected_chemical_unit = -1
        self.view.show_atom_types.connect(self.onShowAtomTypes)
        self.view.hide_atom_types.connect(self.onHideAtomTypes)
        self.view.select_unassigned.connect(self.onSelectUnassigned)
        self.view.setChemicalUnitsCheckBox.checkStateChanged.connect(self.onSetChemicalUnitCheckBoxStateChanged)
        self.view.chemicalUnitBox.selectedChemicalUnitIndexChanged.connect(self.showChemicalUnitDetails)
        self.view.chemicalUnitBox.colorByChemicalUnits.connect(self.colorByChemicalUnitsChanged)
        self.view.chemicalUnitBox.deleteChemicalUnits.connect(self.deleteChemicalUnits)
        self.view.chemicalUnitBox.addChemicalUnit.connect(self.addChemicalUnit)
        self.view.chemicalUnitBox.selectTableSelectedSubsystemAtoms.connect(self.select_subsystem_atoms)
        self.view.chemicalUnitBox.showOnlyChemicalUnits.connect(self.showOnlyChemicalUnits)
        self.view.chemicalUnitBox.showAtomWeights.connect(self.show_atom_weights)
        self.view.chemicalUnitBox.setAtomWeights.connect(self.set_atom_weights)
        self.colour_by_chemical_units_checked = (self.view.chemicalUnitBox.colorByChemicalUnitCheckBox.
                                                 isChecked())
        self.view.chemicalUnitBox.atomWeightChanged.connect(self.onAtomWeightChanged)

    def getCurrentChemicalUnit(self):
        if self.subsystems is not None:
            if len(self.subsystems.subsystems) == 1:
                return 0
            else:
                selected_chemical_units = self.view.chemicalUnitBox.selectedChemicalUnitsIndices()
                if len(selected_chemical_units) == 1:
                    return selected_chemical_units[0]
                else:
                    return -1
        else:
            return -1

    def onSelectUnassigned(self):
        self.onAssignAtomTypes()
        if self.subsystems is not None:
            cu_index = self.getCurrentChemicalUnit()
            if cu_index<0:
                return
            if self.unassigned:
                atom_list = []
                for idx in self.unassigned[cu_index]:
                    atom_list.append(self.subsystems.subsystems[cu_index].atoms[idx])
                self.structure_presenter.select_atoms(atom_list)
        else:
            if self.unassigned:
                atom_list = []
                for idx in self.unassigned:
                    atom_list.append(self.structure_presenter.atomsShown[idx])
                self.structure_presenter.select_atoms(atom_list)

    def onHideAtomTypes(self):
        #self.structure_presenter.show_listed_atom_labels([], [])
        self.structure_presenter.molCanvas.unsetLabels()
    def getAtomTypeNames(self, typeIds, typeNames):
        atomTypeNames = []
        for id in typeIds:
            if id<0:
                atomTypeNames.append('')
            else:
                atomTypeNames.append(typeNames[id])
        return atomTypeNames

    def getUnassigned(self, atoms, typeIds):
        i = 0
        unassigned = []
        for id, atom in zip (typeIds, atoms):
            atom_symm_op = self.structure.getSymmetryOperationStr(atom)
            identity = (atom_symm_op.lower() == "x,y,z")
            n_digits = 0
            for c in id:
                if c.isdigit():
                    n_digits += 1
            if len(id) > 3:
                last3 = id[-3:]
                if last3 == '000' and n_digits == 3 and identity:
                    unassigned.append(i)
            else:
                if len(id) == 0:
                    unassigned.append(i)
            i += 1
        return unassigned

    def onAssignAtomTypes(self):
        self.atom_type_labels = []
        self.unassigned = []
        settings = Settings()
        matts_path = settings.get('paths/MATTS')
        if not matts_path:
            self.no_MATTS_message()
            return
        atom_type_names = discamb_py.get_atom_type_names(matts_path)
        atom_types = []
        if self.subsystems is not None:
            fragments = []
            for subsystem in self.subsystems.subsystems:
                fragments.append(subsystem.atoms)
            atom_type_ids = discamb_py.find_atom_types_for_fragments(self.structure, matts_path, fragments)
            for fragment, ids in zip(fragments, atom_type_ids):
                type_names = self.getAtomTypeNames(ids, atom_type_names)
                self.unassigned.append(self.getUnassigned(fragment, type_names))
                self.atom_type_labels.append(type_names)
        else:
            atom_type_ids = discamb_py.find_atom_types(self.structure, matts_path)
            type_names = self.getAtomTypeNames(atom_type_ids, atom_type_names)
            self.unassigned = self.getUnassigned(self.structure.getAsymetricUnitAtoms(), type_names)
            self.atom_type_labels = type_names

    def onShowAtomTypes(self):
        if self.structure is None:
            return
        self.onAssignAtomTypes()
        if self.atom_type_labels is None:
            self.onAssignAtomTypes()
        if not self.atom_type_labels:
            self.onAssignAtomTypes()
        cu_index = 0
        if self.subsystems is not None:
            cu_index = self.getCurrentChemicalUnit()
            if cu_index<0:
                return
            else:
                self.structure_presenter.show_listed_atom_labels(self.subsystems.subsystems[cu_index].atoms, self.atom_type_labels[cu_index])
        else:
            self.structure_presenter.show_listed_atom_labels(self.structure_presenter.atomsShown, self.atom_type_labels)
        msgBox = QtWidgets.QMessageBox()
        n_unassigned = 0
        if self.subsystems is None:
            n_unassigned = str(len(self.unassigned))
        else:
            n_unassigned = str(len(self.unassigned[cu_index]))
        msgBox.setText("Number of atoms with unassigned atom type: " + n_unassigned)
        msgBox.exec()

    def unsetCrystalStructure(self):
        self.subsystems = None
        self.view.setChemicalUnits([])
    def onSetChemicalUnitCheckBoxStateChanged(self, state):
        show = (state != QtCore.Qt.Unchecked)
        if show:
            self.view.chemicalUnitBox.show()
        else:
            self.view.chemicalUnitBox.hide()

    def setCrystalStructure(self, crystalStructure):
        self.structure = crystalStructure
        self.view.setChemicalUnits([])
        occupancies = crystalStructure.getOccupancies()
        ordered = True
        for occupancy in occupancies:
            if occupancy != 1.0:
                ordered = False
        if ordered:
            self.subsystems = None
            self.view.setChemicalUnitsCheckBox.hide()
            self.view.setChemicalUnitsCheckBox.setCheckState(QtCore.Qt.Unchecked)
        else:
            self.subsystems = CrystalAsSubsystems()
            self.subsystems.setDefault(crystalStructure, taam=True)
            self.view.setChemicalUnitsCheckBox.show()
            self.view.setChemicalUnitsCheckBox.setCheckState(QtCore.Qt.Unchecked)
            self.view.setChemicalUnits(self.subsystems.subsystems)

    def no_MATTS_message(self):
        msgBox = QtWidgets.QMessageBox()
        msgBox.setText("Set MATTS databank path\nEdit->Settings->Set Paths")
        msgBox.exec()

    def save(self, directory):
        settings = Settings()
        matts_path = settings.get('paths/MATTS')
        if not matts_path:
            self.no_MATTS_message()
            return False
        else:
            if not os.path.exists(matts_path):
                msgBox = QtWidgets.QMessageBox()
                msgBox.setText("MATTS databank file not found at the specified path " + matts_path)
                msgBox.exec()
                return False
        model = {
            "model": "TAAM",
            "electron scattering": self.radiation == "electron",
            "assignment info": "atom_types.log",
            "bank path": matts_path
        }
        if self.subsystems is not None:
            subsystems = []
            for subsystem in self.subsystems.subsystems:
                atom_list = []
                for atom, weight in zip(subsystem.atoms, subsystem.weights):
                    atom_label = self.structure.getAtomLabel(atom[0])
                    atom_symm_op = self.structure.getSymmetryOperationStr(atom)
                    is_identity = (atom_symm_op.lower() == "x,y,z")
                    if not is_identity and weight != 0.0:
                        weight = 0.0
                    weight_and_atom = [weight, atom_label]
                    if not is_identity:
                        weight_and_atom.append(atom_symm_op)
                    atom_list.append(weight_and_atom)
                subsystems.append(atom_list)
            model["ordered fragments"] = subsystems

        with open(os.path.join(directory, "aspher.json"),'w') as out:
            out.write(json.dumps(model, indent=4))
        return True

    def run(self, directory, structure_file, hkl_file):
        if not self.save(directory):
            return
        subprocess.run(["./bin/discamb_ff.exe", structure_file, hkl_file], cwd=directory)
        msgBox = QtWidgets.QMessageBox()
        msgBox.setText("done")
        msgBox.exec()


    def onAtomWeightChanged(self, atom_idx, weight):
        self.subsystems.subsystems[self.last_selected_chemical_unit].weights[atom_idx] = weight
        self.subsystems.subsystems[self.last_selected_chemical_unit].user_defined_weight[atom_idx] = True

    def show_atom_weights(self, subsystem_idx):
        if subsystem_idx<0:
            self.structure_presenter.molCanvas.unsetLabels()
            return
        atom_id_list = self.subsystems.subsystems[subsystem_idx].atoms
        labels = []
        for weight in self.subsystems.subsystems[subsystem_idx].weights:
            labels.append(f'{weight:.3f}')
        self.structure_presenter.show_listed_atom_labels(atom_id_list, labels)

    def set_atom_weights(self, subsystem_idx, weight):
        subsystem = self.subsystems.subsystems[subsystem_idx]
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
        if self.view.chemicalUnitBox.chemicalUnitAtomsView.isVisible():
            self.showChemicalUnitDetails(subsystem_idx)

    def select_subsystem_atoms(self, subsystem_idx, atom_indices):
        atom_list = []
        for idx in atom_indices:
            atom_list.append(self.subsystems.subsystems[subsystem_idx].atoms[idx])
        self.structure_presenter.select_atoms(atom_list)

    def showOnlyChemicalUnits(self, chemicalUnits):
        if not chemicalUnits:
            chemicalUnits = list(range(0, len(self.subsystems.subsystems)))
        # if not chemicalUnits:
        #    return
        atoms = []
        for idx in chemicalUnits:
            atoms += self.subsystems.subsystems[idx].atoms
        self.structure_presenter.set_visible_atoms(atoms, save_labels=False)

    def addChemicalUnit(self):
        subsystem = Subsystem()
        subsystem.atoms = self.structure_presenter.getSelectedAtoms()
        subsystem.label = "subsystem_" + str(len(self.subsystems.subsystems)+1)
        subsystem.weights = [1.0]*len(subsystem.atoms)
        self.subsystems.subsystems.append(subsystem)
        self.view.setChemicalUnits(self.subsystems.subsystems)
        if self.colour_by_chemical_units_checked:
            self.colorByChemicalUnit()

    def deleteChemicalUnits(self, chemicalUnitsIndices):
        for index in sorted(chemicalUnitsIndices, reverse=True):
            del self.subsystems.subsystems[index]
        self.view.setChemicalUnits(self.subsystems.subsystems)
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

    def colorByChemicalUnit(self):
        atoms = []
        colors = []
        nChemicalUnits = len(self.subsystems.subsystems)
        if nChemicalUnits == 0:
            return
        for chemicalUnitIdx in range(0, nChemicalUnits):
            atoms += self.subsystems.subsystems[chemicalUnitIdx].atoms
            nAtoms = len(self.subsystems.subsystems[chemicalUnitIdx].atoms)
            qtColor = self.view.chemicalUnitBox.chemicalUnitsModel.item(chemicalUnitIdx,
                                                                                    0).background().color()
            color = [qtColor.red() / 255.0, qtColor.green() / 255.0, qtColor.blue() / 255.0]
            colors += nAtoms * [color]
        self.structure_presenter.colorAtoms(atoms, colors)

    def showChemicalUnitDetails(self, idx):
        self.last_selected_chemical_unit = idx
        chemicalUnit = self.subsystems.subsystems[idx]
        atomLabels = []
        atomSymmOps = []
        for atom in chemicalUnit.atoms:
            atomSymmOps.append(self.structure.getSymmetryOperationStr(atom))
            atomLabels.append(self.structure.getAtomLabel(atom[0]))
        self.view.chemicalUnitBox.showChemicalUnitDetails(idx, atomLabels, atomSymmOps, chemicalUnit.weights)
