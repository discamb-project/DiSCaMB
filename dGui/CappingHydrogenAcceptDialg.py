from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6.QtCore import Signal


class CappingHydrogenAcceptDialg(QtWidgets.QDialog):

    def __init__(self, structure, atoms, capping_hydrogen_atoms):
        super().__init__()
        self.capping_h_accept_checkboxes = []
        layout = QtWidgets.QGridLayout()
        check_all = QtWidgets.QCheckBox("select all")
        check_all.checkStateChanged.connect(self.on_check_all_changed)
        check_all.setChecked(True)
        layout.addWidget(check_all, 0, 0)
        layout_row = 1
        self.atom_neighbour_pairs = []
        for idx, neighbours in capping_hydrogen_atoms.items():
            atom_1_str = (structure.getAtomLabel(atoms[idx][0]) + " " +
                          structure.getSymmetryOperationStr(list(atoms[idx])))
            for neighbour in neighbours:
                atom_2_str = (structure.getAtomLabel(neighbour[0]) + " " +
                              structure.getSymmetryOperationStr(neighbour))
                cb = QtWidgets.QCheckBox("H at " + atom_1_str + " pointing to " + atom_2_str)
                cb.setChecked(True)
                layout.addWidget(cb, layout_row, 0)
                self.capping_h_accept_checkboxes.append(cb)
                self.atom_neighbour_pairs.append([list(atoms[idx]), neighbour])
                layout_row += 1
        ok = QtWidgets.QPushButton("OK")
        cancel = QtWidgets.QPushButton("cancel")
        ok.clicked.connect(self.on_ok)
        cancel.clicked.connect(self.on_cancel)
        ok_cancel_layout = QtWidgets.QHBoxLayout()
        ok_cancel_layout.addWidget(ok)
        ok_cancel_layout.addWidget(cancel)
        layout.addLayout(ok_cancel_layout, layout_row, 0)
        self.setLayout(layout)
        self.setWindowTitle("Add capping hydrogens")
        self.chosen_pairs = []

    def on_check_all_changed(self, state):
        if state == QtCore.Qt.Checked:
            for cb in self.capping_h_accept_checkboxes:
                cb.setChecked(True)
        else:
            for cb in self.capping_h_accept_checkboxes:
                cb.setChecked(False)

    def on_ok(self):
        idx = 0
        for cb in self.capping_h_accept_checkboxes:
            if cb.isChecked():
                self.chosen_pairs.append(self.atom_neighbour_pairs[idx])
            idx += 1
            self.close()

    def get_values(self):
        return self.chosen_pairs

    def on_cancel(self):
        self.close()
