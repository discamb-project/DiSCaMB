from PySide6 import QtWidgets
from PySide6 import QtCore
from PySide6 import QtGui
from PySide6.QtCore import Signal

from RangeSpinBoxes import RangeSpinBoxes
from AtomDisplayStyle import AtomDisplayStyle
from DisplayStyleComboBox import DisplayStyleComboBox

class DisplayControls(QtWidgets.QGroupBox):
    packRangeChanged = Signal(object)
    packToggled = Signal(bool)
    showUnitCellChanged = Signal(bool)
    showNeigboringMolecules = Signal(float)
    showNeigboringAtoms = Signal(float)
    hide = Signal()
    exclusiveSelection = Signal(bool)
    selectByMoleculeClicked = Signal()
    selectByAtomClicked = Signal()
    displayStyleChanged = Signal(object)
    selectAll = Signal()
    showAsymmetricUnitAtoms = Signal()
    showAsymmetricUnitMolecules = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self.mainLayout = QtWidgets.QVBoxLayout()
        self.setLayout(self.mainLayout)

        ###################
        # Show
        ###################

        show_box_layout = QtWidgets.QGridLayout()
        self.show_box = QtWidgets.QGroupBox("Show")
        self.show_box.setMinimumWidth(200)

        self.showAsymmericUnitButton = QtWidgets.QPushButton("asymmetric unit atoms")
        self.showAsymmericUnitButton.clicked.connect(self.onShowAsymmetricUnitAtoms)
        show_box_layout.addWidget(self.showAsymmericUnitButton, 0, 0, 1, 3)

        self.showAsymmericUnitMoleculesButton = QtWidgets.QPushButton("asymmetric unit molecules")
        self.showAsymmericUnitMoleculesButton.clicked.connect(self.onShowAsymmetricUnitMolecules)
        show_box_layout.addWidget(self.showAsymmericUnitMoleculesButton, 1, 0, 1, 3)

        separator = QtWidgets.QFrame()
        separator.setFrameShape(QtWidgets.QFrame.HLine)
        separator.setFrameShadow(QtWidgets.QFrame.Raised)
        show_box_layout.addWidget(separator, 2, 0, 1, 3)

        neighboringLabel = QtWidgets.QLabel("neighboring")
        neighboringLabel.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)
        show_box_layout.addWidget(neighboringLabel, 3, 0)

        self.showAtomsPushButton = QtWidgets.QPushButton("Atoms")
        self.showAtomsPushButton.clicked.connect(self.onShowNeighbouringAtoms)
        self.showMoleculesPushButton = QtWidgets.QPushButton("Molecules")
        self.showMoleculesPushButton.clicked.connect(self.onShowNeighbouringMolecules)
        show_box_layout.addWidget(self.showAtomsPushButton, 4, 0)
        show_box_layout.addWidget(self.showMoleculesPushButton, 4, 1)

        withinLabel = QtWidgets.QLabel("within")
        withinLabel.setAlignment(QtCore.Qt.AlignHCenter | QtCore.Qt.AlignVCenter)

        self.atomShowRangeSpinBox = QtWidgets.QDoubleSpinBox()
        self.atomShowRangeSpinBox.setValue(2.1)
        self.atomShowRangeSpinBox.setDecimals(3)
        self.atomShowRangeSpinBox.setSingleStep(0.1)
        self.atomShowRangeSpinBox.setDisabled(False)

        angstromLabel = QtWidgets.QLabel("&#8491;")
        angstromLabel.setTextFormat(QtCore.Qt.RichText)

        fromSelectedAtomsLabel = QtWidgets.QLabel("from selected Atoms")


        self.showUnitCellCheckbox = QtWidgets.QCheckBox("unit cell")
        self.showUnitCellCheckbox.checkStateChanged.connect(self.onShowUnitCellCheckboxStateChanged)
        show_box_layout.addWidget(withinLabel, 5, 0)
        show_box_layout.addWidget(self.atomShowRangeSpinBox, 5, 1)
        show_box_layout.addWidget(angstromLabel, 5, 2)
        show_box_layout.addWidget(fromSelectedAtomsLabel, 6, 0, 1, 3)

        separator2 = QtWidgets.QFrame()
        separator2.setFrameShape(QtWidgets.QFrame.HLine)
        separator2.setFrameShadow(QtWidgets.QFrame.Raised)

        show_box_layout.addWidget(separator2, 7, 0, 1, 3)
        show_box_layout.addWidget(self.showUnitCellCheckbox, 8, 0)
        # show_unit_cell = QtWidgets.
        # line = QtWidgets.QFrame()
        # line.setFrameShape(QtWidgets.QFrame.HLine)
        # show_box_layout.addWidget(line, 3, 0, 1, 3)

        ###################
        # Hide
        ###################
        '''
        self.hidePushButton = QtWidgets.QPushButton("Hide")
        self.hidePushButton.clicked.connect(self.hide)
        self.mainLayout.addWidget(self.hidePushButton)
        '''
        ###################
        # Pack 
        ###################

        self.pack_box = QtWidgets.QGroupBox("Pack")
        self.pack_box.setCheckable(True)
        self.pack_box.setChecked(False)
        self.pack_box.clicked.connect(self.onPackChanged)
        # pack = QtWidgets.QCheckBox("Pack")
        self.pack_box.clicked.connect(self.showPackRange)
        pack_box_layout = QtWidgets.QGridLayout()

        # pack_box_layout.addWidget(pack, 4, 0)

        self.range_a = RangeSpinBoxes()
        self.range_b = RangeSpinBoxes()
        self.range_c = RangeSpinBoxes()
        self.range_a.changed.connect(self.onPackRangeAChanged)
        self.range_b.changed.connect(self.onPackRangeBChanged)
        self.range_c.changed.connect(self.onPackRangeCChanged)
        pack_box_layout.addLayout(self.range_a, 0, 0, 1, 2)
        pack_box_layout.addLayout(self.range_b, 1, 0, 1, 2)
        pack_box_layout.addLayout(self.range_c, 2, 0, 1, 2)

        self.label_a = QtWidgets.QLabel("a")
        self.label_b = QtWidgets.QLabel("b")
        self.label_c = QtWidgets.QLabel("c")
        pack_box_layout.addWidget(self.label_a, 0, 2)
        pack_box_layout.addWidget(self.label_b, 1, 2)
        pack_box_layout.addWidget(self.label_c, 2, 2)
        self.showPackRange(False)
        self.pack_box.setLayout(pack_box_layout)

        ###################
        # Select
        ###################

        select_box_layout = QtWidgets.QGridLayout()
        self.select_box = QtWidgets.QGroupBox("Selection")
        self.select_box.setMinimumWidth(200)

        self.selectAllPushButton = QtWidgets.QPushButton("Select All")
        select_box_layout.addWidget(self.selectAllPushButton, 0, 0)
        self.selectAllPushButton.clicked.connect(self.selectAll)
        self.hidePushButton = QtWidgets.QPushButton("Hide")
        self.hidePushButton.clicked.connect(self.hide)
        select_box_layout.addWidget(self.hidePushButton, 0, 1)


        self.select_atoms_radioButton = QtWidgets.QRadioButton("Atoms")
        self.select_atoms_radioButton.setChecked(True)
        self.select_atoms_radioButton.clicked.connect(self.selectByAtomClicked)
        self.select_molecules_radioButton = QtWidgets.QRadioButton("Molecules")
        self.select_molecules_radioButton.clicked.connect(self.selectByMoleculeClicked)
        self.exclusiveSelectionCheckBox = QtWidgets.QCheckBox("exclusive selection")
        self.exclusiveSelectionCheckBox.checkStateChanged.connect(self.onExclusiveSelectionCheckBoxStateChanged)
        self.exclusiveSelectionCheckBox.setCheckState(QtCore.Qt.Checked)
        select_box_layout.addWidget(self.select_atoms_radioButton, 1, 0)
        select_box_layout.addWidget(self.select_molecules_radioButton, 1, 1)
        select_box_layout.addWidget(self.exclusiveSelectionCheckBox, 2, 0, 1, 2)
        self.select_box.setLayout(select_box_layout)

        ###################
        #  Display style
        ###################
        self.displayStyleComboBox = DisplayStyleComboBox(AtomDisplayStyle.default())
        self.displayStyleComboBox.currentIndexChanged.connect(self.onDisplayStyleChanged)
        displayStyleLayout = QtWidgets.QHBoxLayout()
        displayStyleLayout.addWidget(QtWidgets.QLabel("Display style"))
        displayStyleLayout.addWidget(self.displayStyleComboBox)

        ###################

        self.mainLayout.addLayout(displayStyleLayout)
        self.mainLayout.addWidget(self.select_box)
        self.show_box.setLayout(show_box_layout)
        self.mainLayout.addWidget(self.show_box)
        self.show_box.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        self.mainLayout.addWidget(self.show_box)
        self.mainLayout.addWidget(self.pack_box)
        # self.diplayStyleComboBox.setSizePolicy(QtWidgets.QSizePolicy.Minimum,QtWidgets.QSizePolicy.Minimum)
        self.mainLayout.addStretch(10)

    def currentAtomStyle(self):
        optionStr = self.displayStyleComboBox.currentText().replace(" ", "_")
        return AtomDisplayStyle[optionStr]
    def onDisplayStyleChanged(self):
        #optionStr = self.diplayStyleComboBox.currentText().replace(" ", "_")
        self.displayStyleChanged.emit(self.currentAtomStyle())#AtomDisplayStyle[optionStr])

    def onExclusiveSelectionCheckBoxStateChanged(self, state):
        self.exclusiveSelection.emit(state == QtCore.Qt.Checked)

    def onShowAsymmetricUnitAtoms(self):
        self.showAsymmetricUnitAtoms.emit()


    def onShowAsymmetricUnitMolecules(self):
        self.showAsymmetricUnitMolecules.emit()

    def onShowNeighbouringMolecules(self):
        self.showNeigboringMolecules.emit(self.atomShowRangeSpinBox.value())

    def onShowNeighbouringAtoms(self):
        self.showNeigboringAtoms.emit(self.atomShowRangeSpinBox.value())

    def onShowUnitCellCheckboxStateChanged(self, state):
        self.showUnitCellChanged.emit(state == QtCore.Qt.Checked)

    def onPackRangeAChanged(self, minimum, maximum):
        self.onPackRangeChanged(minimum, maximum, 0)

    def onPackRangeBChanged(self, minimum, maximum):
        self.onPackRangeChanged(minimum, maximum, 1)

    def onPackRangeCChanged(self, minimum, maximum):
        self.onPackRangeChanged(minimum, maximum, 2)

    def onPackRangeChanged(self, minimum, maximum, component):
        self.packRangeChanged.emit([minimum, maximum, component])

    def onPackChanged(self, checked):
        self.showPackRange(checked)
        self.packToggled.emit(checked)

    def showPackRange(self, show):
        if show:
            self.range_a.show(True)
            self.range_b.show(True)
            self.range_c.show(True)
            self.label_a.show()
            self.label_b.show()
            self.label_c.show()
        else:
            self.label_a.hide()
            self.label_b.hide()
            self.label_c.hide()
            self.range_a.show(False)
            self.range_b.show(False)
            self.range_c.show(False)
