import discamb_py


class Subsystem:
    def __init__(self):
        self.atoms = []
        self.weights = []
        self.user_defined_weight = []
        self.capping_hydrogen_atoms = []
        self.spin_multiplicity = 1
        self.charge = 0
        self.label = "subsystem"
        self.has_user_defined_distributed_multipole_sites = False
        self.distributed_multipole_sites = []


class AtomRepresentative:
    def __init__(self):
        self.subsystemIdx = 0
        self.atomIdx = 0
        self.weight = 0.0


class CrystalAsSubsystems:
    subsystems: list[Subsystem]

    def __init__(self):
        self.subsystems = []
        self.representatives = []
        self.structure = None

    def setDefault(self, crystalStructure, taam=False):
        self.structure = crystalStructure
        self.subsystems = []
        chemicalUnits = []
        atom_groups = crystalStructure.getGrouppedAsymetricUnitBondedAtoms()
        if taam:
            atoms = []
            for group in atom_groups:
                atoms += group
            chemicalUnits.append(atoms)
        else:
            chemicalUnits = atom_groups
        idx = 0
        for chemicalUnit in chemicalUnits:
            idx += 1
            subsystem = Subsystem()
            subsystem.atoms = chemicalUnit
            subsystem.charge = 0
            if self.__even_number_of_electrons(subsystem.atoms, len(subsystem.capping_hydrogen_atoms)):
                subsystem.spin_multiplicity = 1
            else:
                subsystem.spin_multiplicity = 2
            subsystem.label = "subsystem_" + str(idx)
            self.subsystems.append(subsystem)
        # set weights
        for subsystem in self.subsystems:
            if taam:
                subsystem.weights = []
                for atom in subsystem.atoms:
                    atom_symm_op = self.structure.getSymmetryOperationStr(atom)
                    is_identity = (atom_symm_op.lower() == "x,y,z")
                    if is_identity:
                        subsystem.weights.append(1.0)
                    else:
                        subsystem.weights.append(0.0)
            else:
                subsystem.weights = [0] * len(subsystem.atoms)
            subsystem.user_defined_weight = [False] * len(subsystem.atoms)
        if not taam:
            self.setRepresentatives()

    def setRepresentatives(self):
        if not self.subsystems:
            return
        chemical_units = []
        for subsystem in self.subsystems:
            chemical_units.append(subsystem.atoms)
        self.representatives = discamb_py.find_default_representatives(self.structure, chemical_units)
        for subsystem in self.subsystems:
            subsystem.user_defined_weight = [False] * len(subsystem.atoms)
        for rep in self.representatives:
            self.subsystems[rep.substructureIdx].weights[rep.atomIdx] = rep.weight


    def splitChemicalUnit(self, chemicalUnitIndex):
        chemicalUnits = self.structure.splitIntoChemicalUnits(self.subsystems[chemicalUnitIndex].atoms)
        if len(chemicalUnits) == 1:
            return
        self.deleteChemicalUnits([chemicalUnitIndex])
        for chemicalUnit in chemicalUnits:
            self.addSubsystem(chemicalUnit)
        self.setRepresentatives()

    def addSubsystem(self, atoms, capping_hydrogens = None, name=""):
        if not atoms:
            return
        subsystem = Subsystem()
        idx = len(self.subsystems) + 1
        subsystem.label = "subsystem_" + str(idx)
        otherLabels = [subsystem.label for subsystem in self.subsystems]
        while subsystem.label in otherLabels:
            idx += 1
            subsystem.label = "subsystem_" + str(idx)
        subsystem.atoms = atoms
        self.charge = 0
        subsystem.label = "subsystem_" + str(idx)
        subsystem.weights = [0.0]*len(subsystem.atoms)
        if capping_hydrogens is not None:
            subsystem.capping_hydrogen_atoms = capping_hydrogens
        if self.__even_number_of_electrons(subsystem.atoms, len(subsystem.capping_hydrogen_atoms)):
            subsystem.spin_multiplicity = 1
        else:
            subsystem.spin_multiplicity = 2
        self.subsystems.append(subsystem)
        self.setRepresentatives()

    #def deleteChemicalUnits(self, chemicalUnitIndices, atomsToRemove=[]):
    def deleteChemicalUnits(self, chemicalUnitIndices):
        '''
        nChemicalUnits = len(chemicalUnitIndices)
        if nChemicalUnits > 1 and atomsToRemove:
            return
        if not chemicalUnitIndices:
            return
        if atomsToRemove:
            chemicalUnit = self.subsystems[chemicalUnitIndices[0]]
            newAtomList = [atom for atom in chemicalUnit.atoms if atom not in atomsToRemove]
            chemicalUnit.atoms = newAtomList
            if self.__even_number_of_electrons(chemicalUnit.atoms, len(chemicalUnit.capping_hydrogen_atoms)):
                chemicalUnit.spin_multiplicity = 1
            else:
                chemicalUnit.spin_multiplicity = 2
        else:
            nRemoved = 0
            chemicalUnitIndices.sort()
            for idx in chemicalUnitIndices:
                del self.subsystems[idx - nRemoved]
                nRemoved += 1
        self.setRepresentatives()
        '''
        nRemoved = 0
        chemicalUnitIndices.sort()
        for idx in chemicalUnitIndices:
            del self.subsystems[idx - nRemoved]
            nRemoved += 1
        self.setRepresentatives()

    def mergeChemicalUnits(self, chemicalUnitIndices, atomsToAdd):
        if not chemicalUnitIndices:
            return
        # if len(chemicalUnitIndices)==1:
        #    merged = self.subsystems[chemicalUnitIndices[0]]
        #    merged.atoms += [atom for atom in atomsToAdd if atom not in atoms]
        #    if self.__even_number_of_electrons(merged.atoms):
        #        merged.spin_multiplicity = 1
        #    else:
        #        merged.spin_multiplicity = 2
        #    return

        atoms = []
        charge = 0
        for idx in chemicalUnitIndices:
            atoms += self.subsystems[idx].atoms
            charge += self.subsystems[idx].charge
        merged = Subsystem()
        merged.atoms = atoms + [atom for atom in atomsToAdd if atom not in atoms]
        merged.charge = charge
        if self.__even_number_of_electrons(merged.atoms, len(merged.capping_hydrogen_atoms)):
            merged.spin_multiplicity = 1
        else:
            merged.spin_multiplicity = 2
        merged.weights = [0.0]*len(merged.atoms)
        chemicalUnitIndices.sort()
        self.subsystems[chemicalUnitIndices[0]] = merged
        nRemoved = 0
        for idx in chemicalUnitIndices[1:]:
            del self.subsystems[idx - nRemoved]
            nRemoved += 1
        self.setRepresentatives()

    def __even_number_of_electrons(self, chemicalUnit, nCappingHydrogens):
        nElectrons = 0
        for atom in chemicalUnit:
            nElectrons += self.structure.atomicNumber(atom[0])
        return ((nElectrons + nCappingHydrogens) % 2 == 0)
