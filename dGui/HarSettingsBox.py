from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6.QtCore import Signal
from PySide6 import QtGui
from ChemicalUnitsBox import ChemicalUnitsBox
from HarSettings import HarSettings
from PartitionSettingsBox import PartitionSettingsBox
from Settings import Settings

class HarSettingsBox(QtWidgets.QGroupBox):
    basis_set_changed = Signal(str)
    qm_method_changed = Signal(str)
    qm_program_changed = Signal(str)
    n_threads_changed = Signal(int)
    memory_changed = Signal(int)
    memory_units_changed = Signal(str)
    qm_software_changed = Signal(str)
    use_distributed_multipoles = Signal(bool)
    distributed_multipole_range = Signal(float)
    har_x_ray = Signal(bool)

    def __init__(self, parent=None):
        super().__init__("HAR settings", parent)
        self.molView = None
        self.harSettings = HarSettings()

        ############################
        # partition
        ############################

        self.partition_box = PartitionSettingsBox()

        ############################
        # wave function calculation
        ############################

        self.qm_group_box = QtWidgets.QGroupBox("wave function calculation")
        self._set_qm_group_box()

        #########################################
        # more settings/ chemical units switches
        #########################################

        self.show_more_switches_layout = QtWidgets.QHBoxLayout()
        self.expand_more_settings_tool_button = None
        self.expand_chemical_units_box_tool_button = None
        self._set_show_more_switches_layout()

        #########################################
        # more settings
        #########################################

        self.more_settings_group = QtWidgets.QGroupBox()
        self.qm_software_comboBox = QtWidgets.QComboBox(self)
        self.memory_spin_box = QtWidgets.QSpinBox(self)
        self.memory_unit_combo_box = QtWidgets.QComboBox()
        self.threads_spin_box = QtWidgets.QSpinBox(self)
        self._set_more_settings_group()

        #########################################

        separator = QtWidgets.QFrame()
        separator.setFrameShape(QtWidgets.QFrame.HLine)
        separator.setFrameShadow(QtWidgets.QFrame.Raised)

        self.chemicalUnitBox = ChemicalUnitsBox(self)
        self.chemicalUnitBox.hide()

        layout = QtWidgets.QGridLayout()
        layout.addWidget(self.partition_box, 0, 0, 1, 3)
        layout.addWidget(self.qm_group_box, 1, 0, 1, 3)
        layout.addLayout(self.show_more_switches_layout, 2, 0, 1, 3)

        layout.addWidget(self.more_settings_group, 3, 0, 1, 3)
        layout.addWidget(separator, 4, 0, 1, 3)
        layout.addWidget(self.chemicalUnitBox, 5, 0, 1, 3)
        layout.setRowStretch(layout.rowCount(), 10)
        self.setLayout(layout)
        self.crystalStructure = None

    def set_hardware_from_settings(self):
        settings = Settings()
        n_cpu = int(settings.get("hardware/nCPU"))
        memory = settings.get("hardware/memory")
        memory_number = int(memory[0:-2])
        memory_units = memory[-2:]

        self.memory_spin_box.setValue(memory_number)
        self.memory_unit_combo_box.setCurrentText(memory_units)
        self.threads_spin_box.setValue(n_cpu)

    def _set_more_settings_group(self):
        more_settings_group_layout = QtWidgets.QGridLayout()
        qm_software_label = QtWidgets.QLabel("QM software")
        self.qm_software_comboBox.currentTextChanged.connect(self.qm_software_changed)
        self.qm_software_comboBox.addItems(["ORCA", "Gaussian"])

        memory_label = QtWidgets.QLabel("memory")
        self.memory_spin_box.valueChanged.connect(self.memory_changed)
        self.memory_unit_combo_box = QtWidgets.QComboBox()
        self.memory_unit_combo_box.currentTextChanged.connect(self.memory_units_changed)
        self.memory_unit_combo_box.addItems(["MB", "GB"])

        n_threads_label = QtWidgets.QLabel("#CPUs")
        self.threads_spin_box.valueChanged.connect(self.on_n_threads_changed) #self.n_threads_changed)
        self.set_hardware_from_settings()

        radiation_type_label = QtWidgets.QLabel("radiation type")
        radiation_combo_box = QtWidgets.QComboBox()
        radiation_combo_box.addItems(["X-ray", "electron"])
        radiation_combo_box.currentTextChanged.connect(self.on_radiation_changed)

        more_settings_group_layout.addWidget(qm_software_label, 0, 0)
        more_settings_group_layout.addWidget(self.qm_software_comboBox, 0, 1, 1, 2)
        hardware_layout = QtWidgets.QHBoxLayout()
        hardware_layout.addWidget(memory_label)
        hardware_layout.addWidget(self.memory_spin_box)
        hardware_layout.addWidget(self.memory_unit_combo_box)
        hardware_layout.addWidget(n_threads_label)
        hardware_layout.addWidget(self.threads_spin_box)
        #radiation_type_layout = QtWidgets.QHBoxLayout()
        #radiation_type_layout.addWidget(radiation_type_label)
        #radiation_type_layout.addWidget(radiation_combo_box)
        #more_settings_group_layout.addLayout(radiation_type_layout, 2, 0, 1, 3)
        more_settings_group_layout.addWidget(radiation_type_label, 1, 0)
        more_settings_group_layout.addWidget(radiation_combo_box, 1, 1, 1, 2)
        more_settings_group_layout.addLayout(hardware_layout, 2, 0, 1, 3)
        self.more_settings_group.setLayout(more_settings_group_layout)
        self.more_settings_group.hide()

    def on_radiation_changed(self, radiation_type):
        self.har_x_ray.emit(radiation_type == "X-ray")

    def _set_qm_group_box(self):
        qm_group_box_layout = QtWidgets.QGridLayout()
        basis_set_label = QtWidgets.QLabel("basis set")
        basis_set_comboBox = QtWidgets.QComboBox(self)
        basis_set_comboBox.currentTextChanged.connect(self.basis_set_changed)
        basis_set_comboBox.addItems(["def2-SVP", "def2-TZVP", "def2-TZVPP", "def2-QZVP", "cc-pVDZ", "cc-pVTZ",
                                     "cc-pVQZ"])
        basis_set_comboBox.setEditable(True)
        method_label = QtWidgets.QLabel("method")
        method_comboBox = QtWidgets.QComboBox(self)
        method_comboBox.addItems(["PBE", "BLYP", "B3LYP", "Hartree-Fock", "MP2"])
        method_comboBox.setEditable(True)
        method_comboBox.currentTextChanged.connect(self.qm_method_changed)
        method_comboBox.setEditable(True)
        method_comboBox.lineEdit().setPlaceholderText("choose or type")
        distributed_multipole_checkbox = QtWidgets.QCheckBox("distributed multipoles")
        distributed_multipole_checkbox.clicked.connect(self.on_distributed_multipole_checkbox_state_changed)
        range_label = QtWidgets.QLabel("range")
        range_spin_box = QtWidgets.QDoubleSpinBox()
        range_spin_box.valueChanged.connect(self.distributed_multipole_range)
        range_spin_box.setValue(8.0)
        qm_group_box_layout.addWidget(basis_set_label, 0, 0)
        qm_group_box_layout.addWidget(basis_set_comboBox, 0, 1, 1, 2)
        qm_group_box_layout.addWidget(method_label, 1, 0)
        qm_group_box_layout.addWidget(method_comboBox, 1, 1, 1, 2)
        qm_group_box_layout.addWidget(distributed_multipole_checkbox,2,0)
        qm_group_box_layout.addWidget(range_label, 2, 1)
        qm_group_box_layout.addWidget(range_spin_box, 2, 2)
        self.qm_group_box.setLayout(qm_group_box_layout)

    def on_distributed_multipole_checkbox_state_changed(self, state):
        self.use_distributed_multipoles.emit(state != QtCore.Qt.Unchecked)

    def _set_show_more_switches_layout(self):
        self.expand_more_settings_tool_button = QtWidgets.QToolButton()
        self.expand_more_settings_tool_button.setArrowType(QtCore.Qt.RightArrow)
        self.expand_more_settings_Action = QtGui.QAction("show more settings")
        self.expand_more_settings_tool_button.setDefaultAction(self.expand_more_settings_Action)
        self.expand_more_settings_tool_button.clicked.connect(self.on_more_settings_clicked)

        self.more_settings_push_button = QtWidgets.QPushButton("more settings")
        self.more_settings_push_button.clicked.connect(self.on_more_settings_clicked)

        self.expand_chemical_units_box_tool_button = QtWidgets.QToolButton()
        self.expand_chemical_units_box_tool_button.setArrowType(QtCore.Qt.RightArrow)
        self.expandChemicalUnitsAction = QtGui.QAction("show chemical units settings")
        self.expand_chemical_units_box_tool_button.setDefaultAction(self.expandChemicalUnitsAction)
        self.expand_chemical_units_box_tool_button.clicked.connect(self.onExpandChemicalUnitsBoxClicked)
        self.show_more_switches_layout.addWidget(self.expand_more_settings_tool_button)
        self.show_more_switches_layout.addWidget(QtWidgets.QLabel("more settings"))
        self.show_more_switches_layout.addWidget(self.expand_chemical_units_box_tool_button)
        self.show_more_switches_layout.addWidget(QtWidgets.QLabel("Chemical Units"))

    def on_more_settings_clicked(self):
        if self.more_settings_group.isHidden():
            self.more_settings_group.show()
            self.expand_more_settings_tool_button.setArrowType(QtCore.Qt.DownArrow)
        else:
            self.more_settings_group.hide()
            self.expand_more_settings_tool_button.setArrowType(QtCore.Qt.RightArrow)

    def on_n_threads_changed(self, n):
        self.n_threads_changed.emit(n)

    def onExpandChemicalUnitsBoxClicked(self):
        if self.chemicalUnitBox.isVisible():
            self.chemicalUnitBox.hide()
            self.expand_chemical_units_box_tool_button.setArrowType(QtCore.Qt.RightArrow)
        else:
            self.chemicalUnitBox.show()
            self.expand_chemical_units_box_tool_button.setArrowType(QtCore.Qt.DownArrow)
        
    def setMolView(self, molView):
        self.molView = molView
        self.chemicalUnitBox.setMolView(molView)

    def setChemicalUnits(self, chemicalUnits):
        self.chemicalUnitBox.setChemicalUnitsList(chemicalUnits)
        #self.chemicalUnitBox.show()
    
    def setCrystalStructure(self,crystalStructure):
        #self.harSettings.setCrystalStructure(crystalStructure)
        self.chemicalUnitBox.setChemicalUnits(crystalStructure)
        #self.chemicalUnitBox.show()
     
         
     
        

        