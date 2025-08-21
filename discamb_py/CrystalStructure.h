#pragma once

#include "discamb/CrystalStructure/UnitCellContent.h"
#include <string>
#include <vector>

struct ThermalEllispoid {
    std::vector<discamb::Vector3d> axesDirections;
    std::vector<double> adpEigenvalues;
    double u_iso = 0.0;
    bool isotropic = false;
    bool defined = false;
};

enum class PackInlcudeMode {
    ATOM = 0,
    MOLECULE_ATOM = 1,
    MOLECULE_CENTER = 2
};


class CrystalStructure
{
public:
    CrystalStructure();
    CrystalStructure(const std::string &fileName);
    
    int numberOfAtoms() const;
    std::vector<std::string> getAtomNames() const;
    std::vector<std::string> getAtomSymbols() const;
    std::vector<std::vector<double> > getPositions() const;
    std::vector<std::vector<double> > getAdps() const;
    std::vector<double> getOccupancies() const;
    std::vector<int> getBonds() const;
    std::vector<std::pair<int,int> > getBonds(const std::vector<std::vector<int> > &atomIds) const; 
    std::vector< std::vector<std::vector<int> > > splitIntoChemicalUnits(const std::vector<std::vector<int> >& atomIds) const;
    std::vector<std::vector<int> > groupIntoChemicalUnits(const std::vector<std::vector<int> >& atomIds) const;
    std::vector<std::vector<int> > getAsymetricUnitBondedAtoms() const;
    std::vector< std::vector<std::vector<int> > > getGrouppedAsymetricUnitBondedAtoms() const;
    std::vector<std::vector<int> > getAsymetricUnitAtoms() const;
    std::vector<std::vector<int> > getIncludedAtoms(
        double a_min, 
        double a_max, 
        double b_min, 
        double b_max, 
        double c_min, 
        double c_max,
        PackInlcudeMode includeMode) const;
    std::vector< std::vector<std::vector<int> > > getGrouppedAsymetricUnitAtoms() const;
    std::string getAtomLabel(int indexInUnitCell) const;
    std::vector<double> getThermalEllipsoid(int indexInUnitCell) const;
    std::vector<double> getAtomPositionCart(int indexInUnitCell, int nA, int nB, int nC) const;
    //std::string getSymmetryOperationStr(int indexInUnitCell) const;
    std::string getSymmetryOperationStr(const std::vector<int> &atomId) const;
    std::vector<std::vector<double> > getUnitCellVectors() const;
    int atomicNumber(int indexInUnitCell);
    std::vector<std::vector<int> > getNeighboringAtoms(
        std::vector<std::vector<int> > const& atoms, 
        double range, 
        bool includeAllAtomsInIncludingMolecule = true);

    std::vector<std::vector<int> > getBondedAtoms(std::vector<int> const& atom) const;

    double getDistance(const std::vector<int>& atom1, const std::vector<int>& atom2) const;
    const discamb::Crystal& getCrystal() const;
    const discamb::UnitCellContent getUnitCellContent() const;
private:
    discamb::UnitCellContent mUnitCellContent;
    std::vector<discamb::UnitCellContent::AtomID> mAtoms;
    std::vector< std::vector<discamb::UnitCellContent::AtomID> > mAsymmUnitMolecules;
    std::vector< std::vector<discamb::UnitCellContent::AtomID> > mMolecules;
    std::vector<discamb::Vector3d> mMoleculeCenters;
    std::vector<double> mMoleculeRadious;
    std::vector<int> mAtomicNumbers;
    ThermalEllispoid getThermalEllispoid(discamb::UnitCellContent::AtomID const& atomId) const;
    ThermalEllispoid getThermalEllispoid(int idx) const;
    std::vector<ThermalEllispoid> mThermalEllispoids;
    std::vector<int> mBonds;
    std::vector<std::vector<discamb::UnitCellContent::AtomID> > mUnitCellConnectivity;
    void setMolecules();
};

