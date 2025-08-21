from PySide6 import QtWidgets
from IamSettingsBox import IamSettingsBox
from HarSettingsBox import HarSettingsBox
#from HarSettings import HarSettings 

class FormFactorsPanel(QtWidgets.QVBoxLayout): 
        
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
        self.molView = None
        #self.harSettings = HarSettings()
        self.ff_model_box = QtWidgets.QGroupBox("Form factors model")
        self.ff_model_box.setMinimumWidth(200)
        self.radio_iam = QtWidgets.QRadioButton("IAM", self.ff_model_box)
        self.radio_taam = QtWidgets.QRadioButton("TAAM", self.ff_model_box)
        self.radio_har = QtWidgets.QRadioButton("HAR", self.ff_model_box)
        self.radio_iam.setChecked(True)
        qv = QtWidgets.QVBoxLayout()
        qv.addWidget(self.radio_iam)
        qv.addWidget(self.radio_taam)
        qv.addWidget(self.radio_har)
        self.radio_iam.clicked.connect(self.onModelChange)
        self.radio_taam.clicked.connect(self.onModelChange)
        self.radio_har.clicked.connect(self.onModelChange)
        self.ff_model_box.setLayout(qv)
        self.ff_model_box.setSizePolicy(QtWidgets.QSizePolicy.Minimum,QtWidgets.QSizePolicy.Minimum)
        self.addWidget(self.ff_model_box)
        
        self.settings_layout = QtWidgets.QStackedLayout()
        self.iam_settings = IamSettingsBox()
        self.settings_layout.addWidget(self.iam_settings)
        self.harSettingsBox = HarSettingsBox()
        self.settings_layout.addWidget(self.harSettingsBox)           
        self.addLayout(self.settings_layout)
        #self.settings_layout.addStretch()
        '''
        self.settings_layout = QtWidgets.QStackedLayout()
        self.iam_settings = IamSettingsBox()
        self.addWidget(self.iam_settings)
        self.har_settings = HarSettingsBox()
        self.addWidget(self.har_settings)           
        '''
        self.addStretch(10)
    
    def setMolView(self, molView):
        self.molView = molView
        self.harSettingsBox.setMolView(self.molView)
        
    def onModelChange(self):
        if self.radio_iam.isChecked():
            self.settings_layout.setCurrentIndex(0)
        if self.radio_har.isChecked():
            self.settings_layout.setCurrentIndex(1)
        
