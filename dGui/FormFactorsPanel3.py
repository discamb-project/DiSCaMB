from PySide6 import QtWidgets
from PySide6.QtCore import Signal
from IamSettingsBox import IamSettingsBox
from TaamSettingsBox import TaamSettingsBox
from HarSettingsBox import HarSettingsBox
from HarSettings import HarSettings
from enum import IntEnum

class FormFactorModelEnum(IntEnum):
    IAM = 0
    TAAM = 1
    HAR = 2

class FormFactorsPanel3(QtWidgets.QWidget):
    model_changed = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        '''
        ffp = FormFactorsPanel
        ffModel = FFModel(ffp)
        harModel
        taamModel
        iamModel
        harView
        taamView
        iamView

        '''
        self.layout = QtWidgets.QVBoxLayout(self)
        self.molView = None
        self.harSettings = HarSettings()

        ff_model_label = QtWidgets.QLabel("Form factors model")
        ff_model_comboBox = QtWidgets.QComboBox()

        ff_model_comboBox.addItems(["IAM", "TAAM", "HAR"])
        ff_model_layout = QtWidgets.QHBoxLayout()
        ff_model_layout.addWidget(ff_model_label)
        ff_model_layout.addWidget(ff_model_comboBox)
        self.layout.addLayout(ff_model_layout)

        '''
        self.ff_model_box = QtWidgets.QGroupBox("Form factors model")
        self.ff_model_box.setMinimumWidth(200)
        self.radio_iam = QtWidgets.QRadioButton("IAM", self.ff_model_box)
        self.radio_taam = QtWidgets.QRadioButton("TAAM", self.ff_model_box)
        self.radio_har = QtWidgets.QRadioButton("HAR", self.ff_model_box)
        self.radio_iam.setChecked(True)
        self.model_changed.emit(FormFactorModelEnum.IAM)
        qv = QtWidgets.QVBoxLayout()
        qv.addWidget(self.radio_iam)
        qv.addWidget(self.radio_taam)
        qv.addWidget(self.radio_har)
        self.radio_iam.clicked.connect(self.onModelChange)
        self.radio_taam.clicked.connect(self.onModelChange)
        self.radio_har.clicked.connect(self.onModelChange)
        self.ff_model_box.setLayout(qv)
        self.ff_model_box.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        #self.addWidget(self.ff_model_box)
        self.layout.addWidget(self.ff_model_box)
        '''
        self.settings_layout = QtWidgets.QStackedLayout()
        self.iamSettingsBox = IamSettingsBox()
        self.settings_layout.addWidget(self.iamSettingsBox)
        self.taamSettingsBox = TaamSettingsBox()
        self.settings_layout.addWidget(self.taamSettingsBox)
        self.harSettingsBox = HarSettingsBox()
        self.settings_layout.addWidget(self.harSettingsBox)
        self.settings_layout.setSizeConstraint(QtWidgets.QLayout.SizeConstraint.SetMinimumSize)
        #self.addLayout(self.settings_layout)
        #self.addStretch(10)
        self.layout.addLayout(self.settings_layout)
        #spacer = QtWidgets.QSpacerItem()
        #self.layout.addWidget(spacer)
        self.layout.addStretch(20)
        ff_model_comboBox.currentTextChanged.connect(self.on_model_changed)


    def setMolView(self, molView):
        self.molView = molView
        self.harSettingsBox.setMolView(self.molView)

    def on_model_changed(self, model_name):
        idx = 0
        if model_name == "IAM":
            idx = 0
        if model_name == "TAAM":
            idx = 1
        if model_name == "HAR":
            idx = 2
        self.settings_layout.setCurrentIndex(idx)
        self.model_changed.emit(FormFactorModelEnum(idx))

    def onModelChange(self):
        idx = 0
        if self.radio_iam.isChecked():
            idx = 0
        if self.radio_taam.isChecked():
            idx = 1
        if self.radio_har.isChecked():
            idx = 2
        self.settings_layout.setCurrentIndex(idx)
        self.model_changed.emit(FormFactorModelEnum(idx))

