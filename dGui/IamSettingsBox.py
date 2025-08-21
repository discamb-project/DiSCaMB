from PySide6 import QtWidgets
from PySide6.QtCore import Signal

class IamSettingsBox(QtWidgets.QGroupBox):
    radiation_type_changed = Signal(str)
    form_factors_table_changed = Signal(str)
        
    def __init__(self, parent=None):
        super().__init__("IAM settings", parent)
        radio_button_x_ray = QtWidgets.QRadioButton("X-ray", self)
        radio_button_electrons = QtWidgets.QRadioButton("electron", self)
        radio_button_x_ray.setChecked(True)
        self.radiation = "X-ray"
        radio_button_x_ray.clicked.connect(self.on_x_ray)
        radio_button_electrons.clicked.connect(self.on_electrons)
        self.combo_box_x_ray_table = QtWidgets.QComboBox()
        self.combo_box_x_ray_table.addItems(["Waasmeier-Kirfel", "IT92"])
        self.combo_box_x_ray_table.currentTextChanged.connect(self.form_factors_table_changed)
        self.combo_box_electron = QtWidgets.QComboBox()
        self.combo_box_electron.addItems(["International Tables 6A", "International Tables 2A"])
        self.combo_box_electron.currentTextChanged.connect(self.form_factors_table_changed)
        label_scattering_table = QtWidgets.QLabel("scattering table")
        separator = QtWidgets.QFrame()
        separator.setFrameShape(QtWidgets.QFrame.HLine)
        separator.setFrameShadow(QtWidgets.QFrame.Raised)
        layout = QtWidgets.QVBoxLayout()
        layout.addWidget(radio_button_x_ray)
        layout.addWidget(radio_button_electrons)
        layout.addWidget(separator)
        layout.addWidget(label_scattering_table)
        layout.addWidget(self.combo_box_x_ray_table)
        layout.addWidget(self.combo_box_electron)
        layout.addStretch()
        self.setLayout(layout)
        self.on_x_ray()
        self.setSizePolicy(QtWidgets.QSizePolicy.Minimum,QtWidgets.QSizePolicy.Minimum)
    
    def on_x_ray(self):
        self.combo_box_x_ray_table.setVisible(True)
        self.combo_box_electron.setVisible(False)
        self.radiation = "X-ray"
    
    def on_electrons(self):
        self.combo_box_x_ray_table.setVisible(False)
        self.combo_box_electron.setVisible(True)
        self.radiation = "electron"

    def current_radiation(self):
        return self.radiation

    def current_table(self):
        if self.current_radiation == "electron":
            return self.combo_box_electron.currentText()
        else:
            return self.combo_box_x_ray_table.currentText()

