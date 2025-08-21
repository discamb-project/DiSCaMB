import copy
from CrystalStructurePresenter import CrystalStructurePresenter
from DisplayControls import DisplayControls
import discamb_py


class DisplayPresenter: 
    
    def __init__(self, display_controls: DisplayControls, crystal_presenter: CrystalStructurePresenter):
        self.displayControls = display_controls
        self.crystalPresenter = crystal_presenter
        self.displayControls.packRangeChanged.connect(self.onPackRangeChanged)
        self.displayControls.packToggled.connect(self.onPackToggled)
        self.displayControls.showNeigboringMolecules.connect(self.onShowNeighboringMolecules)
        self.displayControls.showNeigboringAtoms.connect(self.onShowNeighboringAtoms)
        self.displayControls.showUnitCellChanged.connect(self.onShowUnitCellChanged)
        self.displayControls.selectAll.connect(self.crystalPresenter.selectAll)
        self.displayControls.hide.connect(crystal_presenter.hideSelected)
        self.displayControls.showAsymmetricUnitAtoms.connect(crystal_presenter.showAsymmetricUnitAtoms)
        self.displayControls.showAsymmetricUnitMolecules.connect(crystal_presenter.showAsymmetricUnitBonded)
        self.displayControls.exclusiveSelection.connect(crystal_presenter.setExclusiveSelection)
        self.displayControls.selectByAtomClicked.connect(self.onSelectByAtom)
        self.displayControls.selectByMoleculeClicked.connect(self.onSelectByMolecule)
        #optionStr = self.displayControls.diplayStyleComboBox.currentText().replace(" ", "_")
        self.crystalPresenter.set_atom_style(self.displayControls.currentAtomStyle())#AtomDisplayStyle[optionStr])
        self.displayControls.displayStyleChanged.connect(self.crystalPresenter.set_atom_style)
        self.packRange = [[0.0,1.0],[0.0,1.0],[0.0,1.0]]
    
    #def onAtomStyleChanged(self, style):
    #    self.crystalPresenter.
    
    def onSelectByAtom(self):
        self.crystalPresenter.setSelectByMolecule(False)
        
    def onSelectByMolecule(self):
        self.crystalPresenter.setSelectByMolecule(True)
    
    def onShowNeighboringMolecules(self, maxDistance):
        self.crystalPresenter.showNeighboringAtoms(maxDistance, True)

    def onShowNeighboringAtoms(self, maxDistance):
        self.crystalPresenter.showNeighboringAtoms(maxDistance, False)

    def onShowUnitCellChanged(self, show):
        self.crystalPresenter.showUnitCell(show)
    
    def onPackToggled(self, selected):
        if selected:
            self.crystalPresenter.showAtomsInRange(self.packRange, discamb_py.PackInlcudeMode.MOLECULE_CENTER)
        else:
            self.crystalPresenter.showAsymmetricUnitBonded()
    
    def onPackRangeChanged(self, rangeChange):
        self.packRange[rangeChange[2]][0]=rangeChange[0]
        self.packRange[rangeChange[2]][1]=rangeChange[1]
        self.crystalPresenter.showAtomsInRange(self.packRange, discamb_py.PackInlcudeMode.MOLECULE_CENTER)
