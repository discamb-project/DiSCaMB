class HarSettings:
    def __init__(self):
        self.settings = {
        "chemical units": [],
        "basis set": "cc-pVDZ",
        "qm method": "PBE",
        "qm program": "orca",
        "power": 1,
        "n cores": 1,
        "memory": "2GB",
        }
        self.crystalStructure = None
    
    def setCrystalStructure(self, crystalStructure):
        self.crystalStructure = crystalStructure
        chemicalUnits = crystalStructure.getGrouppedAsymetricUnitBondedAtoms()
        self.setAnonymousChemicalUnits(chemicalUnits)
        
    def setAnonymousChemicalUnits(self, chemicalUnits):
        self.settings["chemical units"] = []
        idx = 0
        for chemicalUnit in chemicalUnits:
            idx += 1
            self.settings["chemical units"].append(
                dict(label = "subsystem_"+str(idx),
                     atoms = chemicalUnit))
        
        

        