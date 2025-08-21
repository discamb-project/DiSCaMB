from CrystalAsSubsystems import CrystalAsSubsystems

class HirshfeldAtomModelSettings:

    def __init__(self):
        self.subsystemInfo = CrystalAsSubsystems()

    def setStructure(self, structure):
        self.subsystemInfo.setDefault(structure)

    def set_default_representatives(self):
        pass

    def unsetStructure(self):
        self.subsystemInfo = []

    '''
    def setAnonymousChemicalUnits(self, chemicalUnits):
        self.settings["chemical units"] = []
        idx = 0
        for chemicalUnit in chemicalUnits:
            idx += 1
            self.settings["chemical units"].append(
                dict(label = "subsystem_"+str(idx),
                     atoms = chemicalUnit)) 
    '''
        
