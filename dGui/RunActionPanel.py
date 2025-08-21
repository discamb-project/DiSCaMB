from PySide6 import QtWidgets
from enum import IntEnum

class ActionId(IntEnum):
    MAKE_CONFIG = 1
    MAKE_TSC = 2



class RunActionPanel(QtWidgets.QWidget):

    def __init__(self, form_factor_model_presenter, parent=None):
        super().__init__(parent)
        self.hkl_file = ""
        self.form_factor_model_presenter = form_factor_model_presenter
        self.layout = QtWidgets.QVBoxLayout(self)
        self.make_tsc_button = QtWidgets.QPushButton("Make tsc file")
        self.actionGroupBox = QtWidgets.QGroupBox("Action")
        self.checkBoxGroup = QtWidgets.QButtonGroup()
        self.actionGroupBoxLayout = QtWidgets.QVBoxLayout(self.actionGroupBox)

        self.make_configuration_file_check_box = QtWidgets.QCheckBox("make configuration file")
        self.checkBoxGroup.addButton(self.make_configuration_file_check_box, ActionId.MAKE_CONFIG)
        self.actionGroupBoxLayout.addWidget(self.make_configuration_file_check_box)

        self.make_tsc_file_check_box = QtWidgets.QCheckBox("make tsc file")
        self.checkBoxGroup.addButton(self.make_tsc_file_check_box, ActionId.MAKE_TSC)
        self.checkBoxGroup.setExclusive(True)

        self.make_tsc_file_check_box.setChecked(True)
        self.actionGroupBoxLayout.addWidget(self.make_tsc_file_check_box)
        self.layout.addWidget(self.actionGroupBox)

        self.checkBoxGroup.idClicked.connect(self.on_action_chosen)

        self.make_config_group = QtWidgets.QGroupBox()
        self.make_config_group_layout = QtWidgets.QHBoxLayout(self.make_config_group)
        self.make_config_run_button = QtWidgets.QPushButton("Run")
        self.make_config_run_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.make_config_run_button.clicked.connect(self.make_config)
        self.make_config_group_layout.addWidget(self.make_config_run_button)
        self.make_config_group_layout.addStretch(10)
        self.make_config_group.hide()
        self.layout.addWidget(self.make_config_group)

        self.make_tsc_group = QtWidgets.QGroupBox()
        self.make_tsc_group_layout = QtWidgets.QGridLayout(self.make_tsc_group)

        hkl_file_label = QtWidgets.QLabel("Hkl indices file")
        self.hkl_file_name_edit = QtWidgets.QLineEdit()
        self.hkl_file_name_edit.setReadOnly(True)
        self.set_hkl_path_button = QtWidgets.QPushButton("Set")
        self.set_hkl_path_button.clicked.connect(self.set_hkl_file)

        self.make_tsc_run_button = QtWidgets.QPushButton("Run")
        self.make_tsc_run_button.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.make_tsc_run_button.clicked.connect(self.make_tsc)
        self.make_tsc_run_button.setEnabled(False)
        self.make_tsc_group_layout.addWidget(hkl_file_label, 0, 0)
        self.make_tsc_group_layout.addWidget(self.hkl_file_name_edit, 1, 0, 1, 4)
        self.make_tsc_group_layout.addWidget(self.set_hkl_path_button, 1, 4)
        self.make_tsc_group_layout.addWidget(self.make_tsc_run_button, 2, 0)

        #self.make_tsc_group_layout.addStretch(10)
        self.layout.addWidget(self.make_tsc_group)


        self.layout.addStretch(10)

    def set_hkl_file(self):
        file = QtWidgets.QFileDialog.getOpenFileName(self, "hkl indices file",
                                                     filter="hkl indices files (*.hkl *.tsc *.tscb)")
        if file:
            self.hkl_file_name_edit.setText(file[0])
            self.make_tsc_run_button.setEnabled(True)
            self.hkl_file = file[0]
        else:
            self.make_tsc_run_button.setEnabled(False)


    def make_config(self):
        self.form_factor_model_presenter.save()

    def make_tsc(self):
        self.form_factor_model_presenter.run(self.hkl_file)

    def on_action_chosen(self, action_id):
        if action_id == ActionId.MAKE_TSC:
            self.make_config_group.hide()
            self.make_tsc_group.show()
        if action_id == ActionId.MAKE_CONFIG:
            self.make_tsc_group.hide()
            self.make_config_group.show()

