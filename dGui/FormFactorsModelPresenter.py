#from FormFactorsPanel2 import FormFactorsPanel2
from HirshfeldAtomModelPresenter import HirshfeldAtomModelPresenter
from TaamPresenter import TaamPresenter
from IamPresenter import IamPresenter
from FormFactorsPanel3 import FormFactorModelEnum, FormFactorsPanel3
import os
class FormFactorsModelPresenter:
    def __init__(self, formFactorsPanel: FormFactorsPanel3, structurePresenter):
        self.project_dir = os.getcwd()
        self.structure_file = None
        self.hkl_file = None
        self.currentModel = FormFactorModelEnum.IAM
        self.formFactorsPanel = formFactorsPanel
        self.harPresenter = HirshfeldAtomModelPresenter(self.formFactorsPanel.harSettingsBox, structurePresenter)
        self.taamPresenter = TaamPresenter(self.formFactorsPanel.taamSettingsBox, structurePresenter)
        self.iamPresenter = IamPresenter(self.formFactorsPanel.iamSettingsBox)
        self.structurePresenter = structurePresenter
        self.structurePresenter.structureSetSignal.connect(self.setCrystalStructure)
        self.formFactorsPanel.model_changed.connect(self.set_current_model)

    def set_current_model(self, model):
        self.currentModel = model

    def setCrystalStructure(self, crystal_structure, structure_file):
        self.harPresenter.setCrystalStructure(crystal_structure)
        self.taamPresenter.setCrystalStructure(crystal_structure)
        self.structure_file = structure_file
        
    def unsetCrystalStructure(self):
        self.harPresenter.unsetCrystalStructure()
        self.taamPresenter.unsetCrystalStructure()
        self.project_dir = ""
        self.structure_file = ""

    def save(self):
        if self.currentModel == FormFactorModelEnum.IAM:
            self.iamPresenter.save(self.project_dir)
        if self.currentModel == FormFactorModelEnum.TAAM:
            self.taamPresenter.save(self.project_dir)
        if self.currentModel == FormFactorModelEnum.HAR:
            self.harPresenter.save(self.project_dir)

    def run(self, hkl_file):
        if self.currentModel == FormFactorModelEnum.IAM:
            self.iamPresenter.run(self.project_dir, self.structure_file, hkl_file)
        if self.currentModel == FormFactorModelEnum.TAAM:
            self.taamPresenter.run(self.project_dir, self.structure_file, hkl_file)
        if self.currentModel == FormFactorModelEnum.HAR:
            self.harPresenter.run(self.project_dir, self.structure_file, hkl_file)
