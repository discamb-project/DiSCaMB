
import discamb_py
from PySide6 import QtCore
from PySide6.QtCore import Signal, QObject
from PySide6 import QtWidgets
from ElementData import ElementData
from AtomDisplayStyle import AtomDisplayStyle
import cpk_colors
import copy




class DisplayStyleSettings:
    def __init__(self):
        self.wireframe_line_width = 1
        self.bond_radius = 0.05
        self.ball_radius = 0.2
        self.small_ball_radius = 0.05


class CrystalStructurePresenter(QObject):
    structureSetSignal = Signal(object, str)

    def __init__(self, mol_canvas):
        super(CrystalStructurePresenter, self).__init__()
        self.molCanvas = mol_canvas
        self.molCanvas.atomSelectionCallback = self.onAtomsPick
        self.structureData = None
        self.molecules = []
        # 4-int (atomIdx, a, b, c) ids of atoms shown
        self.atomsShown = []
        self.bonds = []
        self.includingMolecule = []
        self.atomColors = []
        self.atomStyle = []
        self.atomStyleAll = AtomDisplayStyle.default()
        self.exclusiveSelection = True
        self.selectByMolecule = False
        self.displayStyleSettings = DisplayStyleSettings()
        # atom indices in self.molCanvas
        self.labeled_atoms = []
        #atom indices in self.molCanvas
        self.selected_atoms = []
        self.measure_distance = False
        self.measure_angle = False

    '''
    def changeAtomStyle(self, style):
        selected = self.molCanvas.selectedAtoms
        if selected:
            for idx in selected:
                self.atomStyle[idx] = style
        else:
            self.atomStyle = [style]*len(self.atomStyle)
            self.atomStyleAll = style
        #self.show_atoms(self.atomsShown, )
        #self.molCanvas.setAtomStyle(atomStyle=style)  # , atomIndices = self.molCanvas.selected)
    '''

    def getSelectedAtoms(self):
        return [self.atomsShown[idx] for idx in self.molCanvas.selectedAtoms]

    def setExclusiveSelection(self, exclusiveSelection):
        self.exclusiveSelection = exclusiveSelection

    def setSelectByMolecule(self, selectByMolecule):
        self.selectByMolecule = selectByMolecule

    def select_atoms(self, atom_list):
        atom_indices = []
        atoms_shown = self.atomsShown
        for atom in atom_list:
            if atom in atoms_shown:
                atom_indices.append(atoms_shown.index(atom))
            self.molCanvas.selectedAtoms = atom_indices
            self.molCanvas.updateAtomMesh()
            self.molCanvas.updateBondMesh()

    def selectAll(self):
        if self.molCanvas is not None:
            self.molCanvas.selectedAtoms = list(range(0, len(self.atomsShown)))
            self.molCanvas.updateAtomMesh()
            self.molCanvas.updateBondMesh()

    def deselect_all(self):
        if self.molCanvas is not None:
            self.molCanvas.selectedAtoms = []
            self.molCanvas.updateAtomMesh()
            self.molCanvas.updateBondMesh()

    def invert_selection(self):
        if self.molCanvas is not None:
            inverted_selection = []
            n_atoms = len(self.molCanvas.atomPositionsScaled)
            for i in range(0, n_atoms):
                if i not in self.molCanvas.selectedAtoms:
                    inverted_selection.append(i)
            self.molCanvas.selectedAtoms = inverted_selection
            self.molCanvas.updateAtomMesh()
            self.molCanvas.updateBondMesh()

    def addContainingMoleculeAtoms(self, atoms):
        molIndices = []
        for atom in atoms:
            includingMolecule = self.includingMolecule[atom]
            if includingMolecule not in molIndices:
                molIndices.append(includingMolecule)
        containingMoleculeAtoms = []
        for moleculeIdx in molIndices:
            containingMoleculeAtoms += self.molecules[moleculeIdx]
        return containingMoleculeAtoms

    def onAtomsPick(self, molCanvas3d, atoms):
        if self.measure_distance:
            if atoms:
                if len(atoms)==1:
                    if len(molCanvas3d.selectedAtoms)==1:
                        molCanvas3d.selectedAtoms.append(atoms[0])
                        molCanvas3d.updateAtomMesh()
                        molCanvas3d.updateBondMesh()
                        dist = self.structureData.getDistance(self.atomsShown[atoms[0]],
                                                              self.atomsShown[molCanvas3d.selectedAtoms[0]])
                        msg_box = QtWidgets.QMessageBox()
                        message = "interatomic distance" + '\u212b'
                        msg_box.setWindowTitle("interatomic distance")
                        #msg_box.setWindowIcon(QtWidgets.QMessageBox.Icon.Information)
                        msg_box.setText('{:.5f}'.format(dist) + ' ' + '\u212b')
                        msg_box.exec()
                        self.deselect_all()
                        return


        if self.selectByMolecule:
            atoms = self.addContainingMoleculeAtoms(atoms)

        if atoms:
            if self.exclusiveSelection:
                molCanvas3d.selectedAtoms = atoms
            else:
                molCanvas3d.selectedAtoms += atoms
            # remove duplicates if any
            molCanvas3d.selectedAtoms = list(dict.fromkeys(molCanvas3d.selectedAtoms))
            molCanvas3d.updateAtomMesh()
            molCanvas3d.updateBondMesh()
        else:
            if molCanvas3d.selectedAtoms and self.exclusiveSelection:
                molCanvas3d.selectedAtoms = []
                molCanvas3d.updateAtomMesh()
                molCanvas3d.updateBondMesh()

    def defaultAtomColors(self):
        self.atomColors = []
        for atom_id in self.atomsShown:
            atomic_number = self.structureData.atomicNumber(atom_id[0])
            color = cpk_colors.cpk_colors_rgb[ElementData.symbol[atomic_number]]
            self.atomColors.append([color[0]/255.0, color[1]/255.0, color[2]/255.0])
        self.update_atoms_view()

    def set_atom_style(self, style):
        selected = self.molCanvas.selectedAtoms
        if not selected:
            self.atomStyle = [style] * len(self.atomStyle)
            self.atomStyleAll = style
        else:
            for idx in selected:
                self.atomStyle[idx] = style
            self.molCanvas.selectedAtoms = []
        self.update_atoms_view()

    def update_atoms_view(self):
        atom_positions = []
        ellipsoids = []
        for atom_id, style in zip(self.atomsShown, self.atomStyle):
            atom_positions.append(
                self.structureData.getAtomPositionCart(atom_id[0], atom_id[1], atom_id[2], atom_id[3]))
            ellipsoids.append(self.get_ellipsoid(atom_id, style))

        wire_bonds = []
        idx = 0
        line_bond_styles = [AtomDisplayStyle.line, AtomDisplayStyle.ball_and_line, AtomDisplayStyle.ellipsoid_and_line]
        for bond in self.bonds:
            if self.atomStyle[bond[0]] in line_bond_styles and self.atomStyle[bond[1]] in line_bond_styles:
                wire_bonds.append(idx)
            idx += 1

        self.molCanvas.setAtomsAndBonds(atom_positions, self.atomColors, bonds=self.bonds, ellipsoids=ellipsoids,
                                        wire_bonds=wire_bonds)

    def set_visible_atoms(self, atom_ids, colors=None, save_labels=True, save_selection=True):
        previous_labeled_atoms = copy.deepcopy(self.labeled_atoms)
        previous_selection = copy.deepcopy(self.molCanvas.selectedAtoms)
        previous_atom_list = copy.deepcopy(self.atomsShown)
        self.deselect_all()
        #self.molCanvas.showLabels()
        existing_atoms = []
        new_atoms = []
        # index in self.molCanvas
        preserved_atom_new_idx = [-1]*len(self.atomsShown)
        for atom_id in atom_ids:
            if atom_id in self.atomsShown:
                existing_atoms.append(atom_id)
                preserved_atom_new_idx[self.atomsShown.index(atom_id)] = len(existing_atoms)-1
            else:
                new_atoms.append(atom_id)
        atom_colors = []
        atom_style = []
        for atom_id in existing_atoms:
            idx = self.atomsShown.index(atom_id)
            atom_colors.append(self.atomColors[idx])
            atom_style.append(self.atomStyle[idx])
        for atom_id in new_atoms:
            atomic_number = self.structureData.atomicNumber(atom_id[0])
            color = cpk_colors.cpk_colors_rgb[ElementData.symbol[atomic_number]]
            atom_colors.append([color[0] / 255.0, color[1] / 255.0, color[2] / 255.0])
            atom_style.append(self.atomStyleAll)
        self.atomsShown = existing_atoms + new_atoms
        self.__setMolecules()
        if colors:
            atom_colors = colors
        self.atomColors = atom_colors
        self.atomStyle = atom_style
        self.bonds = self.structureData.getBonds(self.atomsShown)
        self.update_atoms_view()
        if save_selection:
            selected_atom_indices = []
            select_by_molecule = self.selectByMolecule
            self.selectByMolecule = False
            for previous_idx in previous_selection:
                current_idx = preserved_atom_new_idx[previous_idx]
                if current_idx >= 0:
                    selected_atom_indices.append(current_idx)
            self.molCanvas.onAtomPicking(selected_atom_indices)
            self.selectByMolecule = select_by_molecule
        if save_labels:
            labels=[]
            labeled_atom_indices = []
            for previous_idx in previous_labeled_atoms:
                current_idx = preserved_atom_new_idx[previous_idx]
                if current_idx >= 0:
                    labeled_atom_indices.append(current_idx)
                    labels.append(self.atomsShown[current_idx][0])
            self.molCanvas.showLabels(labels, labeled_atom_indices)
            self.labeled_atoms = labeled_atom_indices

    def showAtomsInRange(self, min_max, mode):
        if self.structureData is not None:
            atomsToShow = self.structureData.getIncludedAtoms(
                min_max[0][0],
                min_max[0][1],
                min_max[1][0],
                min_max[1][1],
                min_max[2][0],
                min_max[2][1],
                mode)
            if atomsToShow is not None:
                self.set_visible_atoms(atomsToShow)

    def __setMolecules(self):
        self.molecules = self.structureData.groupIntoChemicalUnits(self.atomsShown)
        self.includingMolecule = [-1] * len(self.atomsShown)
        for molIdx in range(0, len(self.molecules)):
            for atomIdx in self.molecules[molIdx]:
                self.includingMolecule[atomIdx] = molIdx

    def get_ellipsoid(self, atom_id, style):
        if style in {AtomDisplayStyle.ellipsoid_and_line, AtomDisplayStyle.ellipsoid_and_stick}:
            return self.structureData.getThermalEllipsoid(atom_id[0])
        if style == AtomDisplayStyle.ball_and_stick:
            return [self.displayStyleSettings.ball_radius]
        if style == AtomDisplayStyle.ball_and_line:
            return [self.displayStyleSettings.small_ball_radius]
        if style == AtomDisplayStyle.cpk:
            return [ElementData.vdvRadious[self.structureData.atomicNumber(atom_id[0])]]
        if style == AtomDisplayStyle.line:
            return []
        if style == AtomDisplayStyle.stick:
            return [self.displayStyleSettings.bond_radius]


    def showAsymmetricUnitAtoms(self):
        if self.structureData is None:
            return
        self.set_visible_atoms(self.structureData.getAsymetricUnitAtoms())

    def showAsymmetricUnitBonded(self):
        if self.structureData is None:
            return
        molecule_atoms = self.structureData.getGrouppedAsymetricUnitBondedAtoms()
        atom_ids = []
        for molecule in molecule_atoms:
            atom_ids += molecule
        self.set_visible_atoms(atom_ids)

    @QtCore.Slot()
    def hideSelected(self):
        if self.structureData is None:
            return
        selectedAtoms = []
        for selectedAtomIdx in self.molCanvas.selectedAtoms:
            selectedAtoms.append(self.atomsShown[selectedAtomIdx])
        if not selectedAtoms:
            return
        notSelected = []
        for atom in self.atomsShown:
            if atom not in selectedAtoms:
                notSelected.append(atom)
        self.molCanvas.selectedAtoms = []
        self.set_visible_atoms(notSelected)

    def showNeighboringAtoms(self, rangeAngstrom, select_whole_molecules):
        if self.structureData is None:
            return
        selectedAtoms = []
        for selectedAtomIdx in self.molCanvas.selectedAtoms:
            selectedAtoms.append(self.atomsShown[selectedAtomIdx])
        if not selectedAtoms:
            return
        neighbours = self.structureData.getNeighboringAtoms(selectedAtoms, rangeAngstrom, select_whole_molecules)
        atomsToAdd = []
        for atom in neighbours:
            if atom not in self.atomsShown:
                atomsToAdd.append(atom)
        self.set_visible_atoms(self.atomsShown + atomsToAdd)

    def unset_structure(self):
        self.molCanvas.clear()
        self.structureData = None
        self.molecules = []
        self.atomsShown = []
        self.bonds = []
        self.includingMolecule = []
        self.atomColors = []
        self.atomStyle = []

    def setStructure(self, fileName):
        self.structureData = discamb_py.CrystalStructure(fileName)
        self.showAsymmetricUnitBonded()
        # molecule_atoms = self.structureData.getGrouppedAsymetricUnitBondedAtoms()
        # atom_ids = []
        # for molecule in molecule_atoms:
        #     atom_ids += molecule
        # self.set_visible_atoms(atom_ids)
        # self.showAsymmetricUnitBonded()
        self.structureSetSignal.emit(self.structureData, fileName)

    def showUnitCell(self, show):
        if self.structureData is not None:
            unitCellVectors = self.structureData.getUnitCellVectors()
            self.molCanvas.showUnitCell(show, unitCellVectors)

    def colorAtoms(self, atoms, colors, addAtomIfAbsent=True):
        atomsToAdd = []
        colors_all = []
        for atom in atoms:
            if atom not in self.atomsShown:
                atomsToAdd.append(atom)
        #atomsToShow = self.atomsShown
        atomsToShow = copy.deepcopy(self.atomsShown)
        if atomsToAdd:
            atomsToShow += atomsToAdd
        for atom in atomsToShow:
            if atom in atoms:
                colors_all.append(colors[atoms.index(atom)])
            else:
                colors_all.append(self.atomColors[self.atomsShown.index(atom)])
        self.set_visible_atoms(atomsToShow, colors_all)

    def show_all_labels(self):
        self.showLabels(all_atoms=True)

    def on_measure_distance_changed(self, checked):
        self.measure_distance = checked

    def on_measure_angle_changed(self, checked):
        self.measure_angle = checked

    def showLabels(self, all_atoms=False):
        labels = []
        selected = self.getSelectedAtoms()

        if selected and not all_atoms:
            self.labeled_atoms = self.molCanvas.selectedAtoms
            for atomId in selected:
                labels.append(self.structureData.getAtomLabel(atomId[0]))
            self.molCanvas.showLabels(labels, self.molCanvas.selectedAtoms)
        else:
            for atomId in self.atomsShown:
                labels.append(self.structureData.getAtomLabel(atomId[0]))
            self.labeled_atoms = [i for i in range(0, len(self.atomsShown)) ]
            self.molCanvas.showLabels(labels)

    def show_listed_atom_labels(self, atom_id_list, labels):
        if not labels:
            for atom in atom_id_list:
                labels.append(self.structureData.getAtomLabel(atom[0]))
        idx_in_mol_canvas = []
        labels_to_use = []
        for i in range(0, len(atom_id_list)):
            if atom_id_list[i] in self.atomsShown:
                idx_in_mol_canvas.append(self.atomsShown.index(atom_id_list[i]))
                labels_to_use.append(labels[i])
        self.molCanvas.showLabels(labels_to_use, idx_in_mol_canvas)

        '''
        self.structureData
        if self.selectedAtoms:
            labels = len(self.molecule.xyz)*[""]
            for idx in self.selectedAtoms:
                labels[idx] = self.molecule.atom_symbols[idx]
        else:
            labels = self.molecule.atom_symbols        
        self.labelsVisual = vispy.scene.visuals.Text(labels, pos=pos, face = 'Arial Black', font_size=16, color=(0,0.0,0), parent=self.scene)
        '''
