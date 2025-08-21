from PySide6 import QtWidgets
from PySide6 import QtGui
from PySide6 import QtCore
from PySide6.QtCore import Signal

from ColorItem import ColorItem



class SpinBoxDelegate(QtWidgets.QStyledItemDelegate):

    def __init__(self, parent=None, minimum=-100, maximum=100, init_value=0):
        super().__init__(parent)
        self.minimum = minimum
        self.maximum = maximum
        self.init_value = init_value

    def createEditor(self, parent, option, index):
        editor = QtWidgets.QSpinBox(parent)
        editor.setFrame(False)
        editor.setMinimum(-100)
        editor.setMaximum(100)
        editor.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        return editor

    def setEditorData(self, editor, index):
        value = int(index.model().data(index, QtCore.Qt.ItemDataRole.EditRole))
        if index.column() == 3 and value < 1:
            value = 1
        editor.setValue(value)

    def setModelData(self, editor, model, index):
        editor.interpretText()
        value = editor.value()
        if index.column() == 3:
            if int(value) < 1:
                value = "1"
            if int(value) > 1:
                brush = QtGui.QBrush(QtCore.Qt.yellow)
                if QtGui.QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark:
                    brush = QtGui.QBrush(QtCore.Qt.darkBlue)
                model.setData(index, brush, QtCore.Qt.BackgroundRole)
            if int(value) == 1:
                brush = QtGui.QBrush(QtCore.Qt.white)
                model.setData(index, brush, QtCore.Qt.BackgroundRole)

        model.setData(index, value, QtCore.Qt.ItemDataRole.EditRole)

    def updateEditorGeometry(self, editor, option, index):
        editor.setGeometry(option.rect)


class ChemicalUnitsBox(QtWidgets.QGroupBox):
    selectedChemicalUnitIndexChanged = Signal(int)
    colorByChemicalUnits = Signal(bool)
    mergeChemicalUnits = Signal(object)
    deleteChemicalUnits = Signal(object)
    showOnlyChemicalUnits = Signal(object)
    addChemicalUnit = Signal()
    splitChemicalUnit = Signal(int)
    chemicalUnitChargeChanged = Signal(int, int)
    chemicalUnitSpinMultiplicityChanged = Signal(int, int)
    atomWeightChanged = Signal(int, float)
    showAtomWeights = Signal(int)
    setAtomWeights = Signal(int, float)
    selectTableSelectedSubsystemAtoms = Signal(int, object)
    showMultipoleSites = Signal(int)
    hideMultipoleSites = Signal()
    setMultipoleSites = Signal(int)

    def __init__(self, parent=None, taam=False):
        # super().__init__("Chemical Units", parent)
        super().__init__(parent)
        self.chemicalUnits = []
        self.setWindowTitle('simple list view')
        self.settings = {}
        self.molView = None
        self.weights_shown = False
        self.taam = taam
        layout = QtWidgets.QGridLayout()

        self.colorByChemicalUnitCheckBox = QtWidgets.QCheckBox("color by chemical unit")
        self.colorByChemicalUnitCheckBox.checkStateChanged.connect(self.onColorAtomsByChemicalUnitChanged)
        # colorByChemicalUnitCheckBox.setCheckable(True)
        layout.addWidget(self.colorByChemicalUnitCheckBox, 0, 0, 1, 2)
        nItems = 4
        if self.taam:
            nItems = 2
        self.chemicalUnitsModel = QtGui.QStandardItemModel(0, nItems)
        self.chemicalUnitsModel.itemChanged.connect(self.on_chemical_unit_model_changed)
        self.chemicalUnitsView = QtWidgets.QTableView()
        delegateCharge = SpinBoxDelegate()

        #self.chemicalUnitsView.setItemDelegateForColumn(3, delegateCharge)
        delegateSpin = SpinBoxDelegate()
        #self.chemicalUnitsView.setItemDelegateForColumn(2, delegateCharge)
        self.chemicalUnitsView.doubleClicked.connect(self.cuViewClicked)
        self.chemicalUnitsView.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        self.chemicalUnitsView.setColumnWidth(0, 10)
        if self.taam:
            self.chemicalUnitsModel.setHorizontalHeaderLabels(["", "name"])
        else:
            self.chemicalUnitsModel.setHorizontalHeaderLabels(["", "name", "charge", "spin mult."])
        self.chemicalUnitsView.setModel(self.chemicalUnitsModel)
        self.chemicalUnitsView.resizeColumnToContents(0)
        if not self.taam:
            self.chemicalUnitsView.resizeColumnToContents(2)
            self.chemicalUnitsView.resizeColumnToContents(3)
        self.chemicalUnitsView.resizeRowsToContents()
        self.chemicalUnitsView.selectionModel().selectionChanged.connect(self.selectionChanged)
        self.chemicalUnitsView.setSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)

        layout.addWidget(self.chemicalUnitsView, 1, 0, 1, 3)
        self.chemicalUnitsModel.itemChanged.connect(self.itemChanged)


        #self.chemicalUnitsActionsComboBox = QtWidgets.QComboBox()
        #if self.taam:
        #    self.chemicalUnitsActionsComboBox.addItems(["chemical units", "merge", "delete", "add"])
        #else:
        #    self.chemicalUnitsActionsComboBox.addItems(["chemical units", "merge", "delete", "add", "split"])
        #self.chemicalUnitsActionsComboBox.currentIndexChanged.connect(self.onChemicalUnitsActionComboBoxIndexChanged)
        #layout.addWidget(self.chemicalUnitsActionsComboBox, 2, 0)

        #self.actionMergeCus = QtGui.QAction("merge", self)
        #self.actionMergeCus.triggered.connect(self.onChemicalUnitsMerge)
        # self.actionOpen.setObjectName(u"actionOpen")
        #self.deleteCus = QtGui.QAction("delete", self)
        #self.deleteCus.triggered.connect(self.onDeleteChemicalUnits)
        chemical_units_menu = QtWidgets.QMenu(self)
        #chemical_units_menu.addActions([self.actionMergeCus, self.deleteCus])
        if not self.taam:
            chemical_units_menu.addAction("merge", self.onChemicalUnitsMerge)
        chemical_units_menu.addAction("delete", self.onDeleteChemicalUnits)
        chemical_units_menu.addAction("add", self.onAddChemicalUnit)
        if not self.taam:
            chemical_units_menu.addAction("split", self.onSplitChemicalUnit)
        self.chemicalUnitsPushButton = QtWidgets.QPushButton("chemical units")
        self.chemicalUnitsPushButton.setMenu(chemical_units_menu)
        self.chemicalUnitsPushButton.clicked.connect(self.chemicalUnitsPushButton.showMenu)
        layout.addWidget(self.chemicalUnitsPushButton, 2, 0)

        #self.mergeCusPushButton = QtWidgets.QPushButton("merge")
        #self.deleteCusPushButton = QtWidgets.QPushButton("delete")
        #self.addCuPushButton = QtWidgets.QPushButton("add")
        #self.splitCuPushButton = QtWidgets.QPushButton("split")
        self.showOnlyPushButton = QtWidgets.QPushButton("show only")

        #self.mergeCusPushButton.clicked.connect(self.onChemicalUnitsMerge)
        #self.deleteCusPushButton.clicked.connect(self.onDeleteChemicalUnits)
        #self.addCuPushButton.clicked.connect(self.addChemicalUnit)
        #self.splitCuPushButton.clicked.connect(self.onSplitChemicalUnit)
        self.showOnlyPushButton.clicked.connect(self.onShowOnlyChemicalUnits)

        #layout.addWidget(self.mergeCusPushButton, 2, 0)
        #layout.addWidget(self.deleteCusPushButton, 2, 1)
        #layout.addWidget(self.addCuPushButton, 2, 2)
        #layout.addWidget(self.splitCuPushButton, 3, 0)
        layout.addWidget(self.showOnlyPushButton, 2, 1)
        ######################################
        # atom list
        ######################################

        self.showAtomsButton = QtWidgets.QPushButton("atoms list")
        self.showAtomsButton.setCheckable(True)
        self.showAtomsButton.clicked.connect(self.showAtomsButtonClicked)
        layout.addWidget(self.showAtomsButton, 2, 2)

        atom_weights_menu = QtWidgets.QMenu(self)
        atom_weights_menu.addAction("show", self.onShowAtomsWeights)
        atom_weights_menu.addAction("hide", self.onHideAtomsWeights)
        atom_weights_menu.addAction("set", self.onSetAtomsWeights)
        self.atomWeightsPushButton = QtWidgets.QPushButton("weights")
        self.atomWeightsPushButton.setMenu(atom_weights_menu)
        self.atomWeightsPushButton.clicked.connect(self.atomWeightsPushButton.showMenu)
        layout.addWidget(self.atomWeightsPushButton, 3, 0)


        #self.showAtomWeightsButton = QtWidgets.QPushButton("show weights")
        #self.showAtomWeightsButton.clicked.connect(self.onShowAtomsWeights)
        #layout.addWidget(self.showAtomWeightsButton, 3, 0)

        #self.setAtomWeightsButton = QtWidgets.QPushButton("set weights")
        #self.setAtomWeightsButton.clicked.connect(self.onSetAtomsWeights)
        #layout.addWidget(self.setAtomWeightsButton, 3, 1)

        self.selectChosenAtomsButton = QtWidgets.QPushButton("select chosen")
        self.selectChosenAtomsButton.clicked.connect(self.onSelectChosenAtoms)
        layout.addWidget(self.selectChosenAtomsButton, 3, 1)

        multipole_sites_menu = QtWidgets.QMenu(self)
        multipole_sites_menu.addAction("show", self.onShowMultipoleSites)
        multipole_sites_menu.addAction("hide", self.onHideMultipoleSites)
        multipole_sites_menu.addAction("set", self.onSetMultipoleSites)
        self.atomMultipolePushButton = QtWidgets.QPushButton("multipole sites")
        self.atomMultipolePushButton.setMenu(multipole_sites_menu)
        self.atomMultipolePushButton.clicked.connect(self.atomMultipolePushButton.showMenu)
        if not self.taam:
            layout.addWidget(self.atomMultipolePushButton, 3, 2)


        #self.multipoleComboBox = QtWidgets.QComboBox()
        #self.multipoleComboBox.addItems(["multipole sites", "show", "set"])
        #self.multipoleComboBox.currentIndexChanged.connect(self.onMultipoleComboBoxIndexChanged)
        #layout.addWidget(self.multipoleComboBox, 3, 2)
        '''
        multipole_buttons_layout = QtWidgets.QHBoxLayout()
        self.showMultipolesButton = QtWidgets.QPushButton("show multipole sites")
        self.setMultipolesButton = QtWidgets.QPushButton("set multipole sites")
        multipole_buttons_layout.addWidget(self.showMultipolesButton)
        multipole_buttons_layout.addWidget(self.setMultipolesButton)
        layout.addLayout(multipole_buttons_layout, 6, 0, 1, 3)
        '''
        self.chemicalUnitAtomsModel = QtGui.QStandardItemModel(0, 3, self)
        self.chemicalUnitAtomsModel.itemChanged.connect(self.onAtomListElementChanged)
        self.chemicalUnitAtomsView = QtWidgets.QTableView()
        self.chemicalUnitAtomsView.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectRows)
        #self.chemicalUnitAtomsModel.setHorizontalHeaderLabels(["label", "symmetry operation", "weight"])
        self.chemicalUnitAtomsModel.setHorizontalHeaderLabels(["label", "symm.", "wght."])
        self.chemicalUnitAtomsView.setModel(self.chemicalUnitAtomsModel)
        self.chemicalUnitAtomsView.resizeColumnToContents(1)

        # self.chemicalUnitAtomsView.verticalHeader().setSectionResizeMode(QtWidgets.QHeaderView.Stretch)

        self.chemicalUnitAtomsView.hide()
        layout.addWidget(self.chemicalUnitAtomsView, 4, 0, 1, 3)
        # pb = QtWidgets.QLabel("H<sup>2</sup>O")
        # pb.setTextFormat(QtCore.Qt.MarkdownText)
        # layout.addWidget(pb,4,0)

        # -----------------
        layout.setRowStretch(layout.rowCount(), 10)
        self.setLayout(layout)
        self.crystalStructure = None
        self.chemicalUnits = None

    def onChemicalUnitsActionComboBoxIndexChanged(self, idx):
        if idx == 0:
            return
        if idx == 1:
            self.onChemicalUnitsMerge()
        if idx == 2:
            self.onDeleteChemicalUnits()
        if idx == 3:
            self.addChemicalUnit()
        if idx == 4:
            self.onSplitChemicalUnit()
        self.chemicalUnitsActionsComboBox.setCurrentIndex(0)

    def onAddChemicalUnit(self):
        self.addChemicalUnit.emit()

    def onShowMultipoleSites(self):
        cu_indices = self.selectedChemicalUnitsIndices()
        if len(cu_indices)==1:
            self.showMultipoleSites.emit((cu_indices[0]))
    def onHideMultipoleSites(self):
        self.hideMultipoleSites.emit()

    def onSetMultipoleSites(self):
        cu_indices = self.selectedChemicalUnitsIndices()
        if len(cu_indices)==1:
            self.setMultipoleSites.emit((cu_indices[0]))

    def onMultipoleComboBoxIndexChanged(self, index):
        #"multipole sites", "show", "set"    showMultipoleSites = Signal(int)    setMultipolesites = Signal(int)
        if index == 0:
            return
        cu_indices = self.selectedChemicalUnitsIndices()
        if len(cu_indices) == 1:
            cu_idx = cu_indices[0]
            if index == 1:
                if self.multipoleComboBox.itemText(1) == "show":
                    self.showMultipoleSites.emit(cu_idx)
                    self.multipoleComboBox.setItemText(1, "hide")
                else:
                    self.hideMultipoleSites.emit(cu_idx)
                    self.multipoleComboBox.setItemText(1, "show")
            if index == 2:
                self.setMultipoleSites.emit(cu_idx)
        self.multipoleComboBox.setCurrentIndex(0)

    def onShowAtomsWeights(self):
        if self.weights_shown:
            self.showAtomWeights.emit(-1)
            self.weights_shown = False
        else:
            cu_indices = self.selectedChemicalUnitsIndices()
            if len(cu_indices) == 1:
                self.showAtomWeights.emit(cu_indices[0])
                self.weights_shown = True

    def onHideAtomsWeights(self):
        self.showAtomWeights.emit(-1)

    def onSetAtomsWeights(self):
        cu_indices = self.selectedChemicalUnitsIndices()
        if len(cu_indices) == 1:
            weight, ok = QtWidgets.QInputDialog.getDouble(self, "Set weights", "weigh values")
            if ok:
                self.setAtomWeights.emit(cu_indices[0], weight)

    def onSelectChosenAtoms(self):
        cu_indices = self.selectedChemicalUnitsIndices()
        atom_indices = []
        if len(cu_indices) == 1:
            if self.chemicalUnitAtomsView.isVisible():
                selectionModel = self.chemicalUnitAtomsView.selectionModel()
                if selectionModel.hasSelection():
                    rows = selectionModel.selectedRows()
                    for row in rows:
                        atom_indices.append(row.row())
                self.selectTableSelectedSubsystemAtoms.emit(cu_indices[0], atom_indices)


    def on_chemical_unit_model_changed(self, item: QtGui.QStandardItem):
        if item.column() == 2:
            self.chemicalUnitChargeChanged.emit(item.row(), int(item.text()))
        if item.column() == 3:
            self.chemicalUnitSpinMultiplicityChanged.emit(item.row(), int(item.text()))

    def onAtomListElementChanged(self, item):
        if item.column() == 2:
            self.atomWeightChanged.emit(item.row(), float(item.text()))

    def onShowOnlyChemicalUnits(self):
        chemical_unit_indices = [rowIndex.row() for rowIndex in self.chemicalUnitsView.selectionModel().selectedRows()]
        self.showOnlyChemicalUnits.emit(chemical_unit_indices)

    def onSplitChemicalUnit(self):
        indicesList = self.selectedChemicalUnitsIndices()
        if len(indicesList) == 1:
            self.splitChemicalUnit.emit(indicesList[0])

    def selectedChemicalUnitsIndices(self):
        selectionModel = self.chemicalUnitsView.selectionModel()
        if not selectionModel.hasSelection():
            return []
        chemicalUnitIndices = selectionModel.selectedRows()
        return [idx.row() for idx in chemicalUnitIndices]

    def onDeleteChemicalUnits(self):
        indicesList = self.selectedChemicalUnitsIndices()
        if indicesList:
            self.deleteChemicalUnits.emit(indicesList)

    def onChemicalUnitsMerge(self):
        indicesList = self.selectedChemicalUnitsIndices()
        if indicesList:
            self.mergeChemicalUnits.emit(indicesList)

    def showAtomsButtonClicked(self, checked):
        if checked:
            self.chemicalUnitAtomsView.show()
            self.showAtomsButton.setText("hide atom list")
        else:
            self.chemicalUnitAtomsView.hide()
            self.showAtomsButton.setText("show atom list")

    def cuViewClicked(self, index):
        if index.column() == 0:
            cuIdx = index.row()
            item = self.chemicalUnitsModel.item(cuIdx, index.column())
            background = item.background()
            color = QtWidgets.QColorDialog.getColor(initial=background.color())
            if not color.isValid():
                return
            background.setColor(color)
            item.setBackground(background)

    def setMolView(self, molView):
        self.molView = molView

    def onColorAtomsByChemicalUnitChanged(self, state):
        selected = (state == QtCore.Qt.Checked)
        self.colorByChemicalUnits.emit(selected)

    def onChemicalUnitColorClicked(self):
        currentColor = self.chemicalUnitColorButton.palette().color(QtGui.QPalette.Button)
        color = QtWidgets.QColorDialog.getColor(initial=currentColor)
        if not color.isValid():
            return
        palette = QtGui.QPalette()
        palette.setColor(QtGui.QPalette.Button, color)
        self.chemicalUnitColorButton.setPalette(palette)
        cuIdx = self.chemicalUnitsView.currentIndex().row()
        colorItem = ColorItem(color)
        self.chemicalUnitsModel.setItem(cuIdx, 0, colorItem)

    def addChemicalUnitToList(self, color, name, charge, spin):
        colorItem = ColorItem(color)
        spinItem = QtGui.QStandardItem(str(spin))
        bcg = spinItem.background()
        spinItem.setTextAlignment(QtGui.Qt.AlignHCenter | QtGui.Qt.AlignVCenter)
        '''
        if spin > 1:
            brush = QtGui.QBrush(QtCore.Qt.yellow)
            if QtGui.QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark:
                # brush = QtGui.QBrush(QtCore.Qt.darkBlue)
                brush = QtGui.QBrush(QtGui.QColor(128, 64, 0))
            spinItem.setBackground(brush)
        else:
            if QtGui.QGuiApplication.styleHints().colorScheme() == QtCore.Qt.ColorScheme.Dark:
                brush = QtGui.QBrush(QtGui.QColor(32, 32, 32))
                spinItem.setBackground(bcg)
        '''
        if self.taam:
            self.chemicalUnitsModel.appendRow(
                [colorItem,
                 QtGui.QStandardItem(name)])
        else:
            self.chemicalUnitsModel.appendRow(
            [colorItem,
             QtGui.QStandardItem(name),
             QtGui.QStandardItem(str(charge)),
             spinItem])  # QtGui.QStandardItem(str(spin))])
        nRows = self.chemicalUnitsModel.rowCount()
        rowHight = self.chemicalUnitsView.rowHeight(0)
        headerHight = self.chemicalUnitsView.horizontalHeader().rect().height()
        self.chemicalUnitsView.setMaximumHeight(
            nRows * rowHight + headerHight + self.chemicalUnitsView.verticalScrollBar().rect().height())

    def itemChanged(self, item):
        if item.column() == 3:
            if int(item.text()) < 1:
                item.setValue(1)
        self.selectedChemicalUnitIndexChanged.emit(item.row())

    def showChemicalUnitDetails(self, index):
        chemicalUnitSettings = self.settings["chemical units"][index]
        representatives_cu = []
        self.chemicalUnitNameEdit.setText(chemicalUnitSettings["label"])
        self.spinMultiplicitySpinBox.setValue(chemicalUnitSettings["spin multiplicity"])
        self.chargeSpinBox.setValue(chemicalUnitSettings["charge"])
        self.chemicalUnitAtomsModel.setRowCount(0)

        colorItem = self.chemicalUnitsModel.item(index, 0)
        palette = QtGui.QPalette()
        palette.setColor(QtGui.QPalette.Button, colorItem.background().color())
        self.chemicalUnitColorButton.setPalette(palette)

        for atom in chemicalUnitSettings["atoms"]:
            label = self.crystalStructure.getAtomLabel(atom[0])
            symmOp = self.crystalStructure.getSymmetryOperationStr(atom)
            self.chemicalUnitAtomsModel.appendRow(
                [QtGui.QStandardItem(label), QtGui.QStandardItem(symmOp)])

    def showChemicalUnitDetails(self, index, atomLabels, atomSymmOps, weights):
        self.chemicalUnitAtomsModel.setRowCount(0)
        for label, symmOp, weight in zip(atomLabels, atomSymmOps, weights):
            self.chemicalUnitAtomsModel.appendRow(
                [QtGui.QStandardItem(label), QtGui.QStandardItem(symmOp), QtGui.QStandardItem(str(weight))])

    def selectionChanged(self, selected, deselected):
        indexes = selected.indexes()
        if indexes:
            self.selectedChemicalUnitIndexChanged.emit(indexes[0].row())

    def setChemicalUnitsList(self, chemicalUnits):
        self.chemicalUnitAtomsView.hide()
        self.showAtomsButton.setText("show atom list")
        nRows = self.chemicalUnitAtomsModel.rowCount()
        self.chemicalUnitAtomsModel.removeRows(0, nRows)
        nRows = self.chemicalUnitsModel.rowCount()
        if nRows > 0:
            self.chemicalUnitsModel.removeRows(0, nRows)
        nChemicalUnits = len(chemicalUnits)
        if nChemicalUnits == 0:
            return
        color_step = 360.0 / nChemicalUnits
        index = 1
        for chemicalUnit in chemicalUnits:
            color = QtGui.QColor()
            color.setHsv(30 + (index - 1) * color_step, 255, 255, 255)
            index += 1
            self.addChemicalUnitToList(color, chemicalUnit.label, chemicalUnit.charge, chemicalUnit.spin_multiplicity)
        # self.showChemicalUnitDetails(0)

    def populateChemicalUnitsModel(self):
        nChemicalUnits = len(self.settings["chemical units"])
        if nChemicalUnits == 0:
            return
        color_step = 360.0 / nChemicalUnits
        index = 1
        for chemicalUnit in self.settings["chemical units"]:
            color = QtGui.QColor()
            color.setHsv(30 + (index - 1) * color_step, 255, 255, 255)
            index += 1
            self.addChemicalUnitToList(color, chemicalUnit["label"], chemicalUnit["charge"],
                                       chemicalUnit["spin multiplicity"])
        self.showChemicalUnitDetails(0)

    def __even_number_of_electrons(self, chemicalUnit):
        nElectrons = 0
        for atom in chemicalUnit:
            nElectrons += self.crystalStructure.atomicNumber(atom[0])
        return (nElectrons % 2 == 0)

    def setAnonymousChemicalUnits(self, chemicalUnits):
        self.settings["chemical units"] = []
        idx = 1
        for chemicalUnit in chemicalUnits:
            self.settings["chemical units"].append(
                dict(label="subsystem_" + str(idx),
                     atoms=chemicalUnit))
