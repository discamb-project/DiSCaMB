from PySide6 import QtWidgets
from ChemicalUnitsBox import ChemicalUnitsBox
from PySide6 import QtCore

class TaamSettingsBox(QtWidgets.QGroupBox):
    set_chemical_units_checkbox_state_changed = QtCore.Signal(bool)
    check_atom_types = QtCore.Signal()
    show_atom_types = QtCore.Signal()
    hide_atom_types = QtCore.Signal()
    select_unassigned = QtCore.Signal()
        
    def __init__(self, parent=None):
        super().__init__("TAAM settings", parent)
        radio_button_x_ray = QtWidgets.QRadioButton("X-ray", self)
        radio_button_electrons = QtWidgets.QRadioButton("electron", self)
        radio_button_x_ray.setChecked(True)
        self.setChemicalUnitsCheckBox = QtWidgets.QCheckBox("set chemical units")
        self.setChemicalUnitsCheckBox.hide()
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(radio_button_x_ray)
        layout.addWidget(radio_button_electrons)

        atom_types_layout = QtWidgets.QHBoxLayout()
        #atomTypesLabel = QtWidgets.QLabel("atom types")
        selectUnassignedAtomsButton = QtWidgets.QPushButton("select unassigned")
        selectUnassignedAtomsButton.clicked.connect(self.select_unassigned)
        self.showAtomTypesButton = QtWidgets.QPushButton("show atom types")
        self.showAtomTypesButton.clicked.connect(self.onShowAtomTypes)
        #atom_types_layout.addWidget(atomTypesLabel)
        atom_types_layout.addWidget(self.showAtomTypesButton)
        atom_types_layout.addWidget(selectUnassignedAtomsButton)
        layout.addLayout(atom_types_layout)
        layout.addWidget(self.setChemicalUnitsCheckBox)
        self.chemicalUnitBox = ChemicalUnitsBox(self, taam=True)
        self.chemicalUnitBox.hide()
        layout.addWidget(self.chemicalUnitBox)
        self.setLayout(layout)
        layout.addStretch(10)
        self.setFlat(True)

    def onShowAtomTypes(self):
        if self.showAtomTypesButton.text() == "show atom types":
            self.showAtomTypesButton.setText("hide atom types")
            self.show_atom_types.emit()
        else:
            self.showAtomTypesButton.setText("show atom types")
            self.hide_atom_types.emit()

    def setChemicalUnits(self, chemicalUnits):
        self.chemicalUnitBox.setChemicalUnitsList(chemicalUnits)

    def unsetChemicalUnits(self, chemicalUnits):
        self.chemicalUnitBox.hide()
        self.chemicalUnitBox.setChemicalUnitsList([])

    def checkAtomTypes(self):
        pass

