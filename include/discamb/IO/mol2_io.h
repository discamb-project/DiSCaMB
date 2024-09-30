#pragma once

#include "discamb/MathUtilities/Vector3.h"
#include "discamb/BasicChemistry/ChemicalElement.h"
#include "discamb/BasicChemistry/MoleculeData.h"

#include <vector>
#include <string>

namespace discamb
{
    /**
    * \addtogroup IO IO
    * @{
    */


    namespace mol2_io {

        struct Mol2Data {
            void toMoleculeData(MoleculeData& data) const;
            void split(std::vector<Mol2Data> &substructures) const;
            std::string moleculeName;
            //enum class MoleculeType{SMALL, BIOPOLYMER, PROTEIN, NUCLEIC_ACID, SACCHARIDE};
            
            //MoleculeType moleculeTypeFromString(const std::string& s);
            //std::string moleculeTypeToString(MoleculeType s);

            //MoleculeType moleculeType;
            std::string moleculeType;
            std::string moleculeComment;
//            enum class ChargeType {
//                NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL,
//                PULLMAN, GAUSS80_CHARGES, AMPAC_CHARGES,
//                MULLIKEN_CHARGES, DICT_CHARGES, MMFF94_CHARGES,
//                USER_CHARGES};

//            ChargeType chargeTypeFromString(const std::string& s);
//            std::string chargeTypeToString(ChargeType s);
//
//            MoleculeType typeFromString(const std::string& s);
//            std::string typeToString(MoleculeType s);


            //ChargeType chargeType;
            std::string chargeType;
            std::vector<std::string> atomName;
            std::vector <Vector3d> atomPosition;
            std::vector <std::string> atomType;
            std::vector<std::pair<int, int> > bonds;
            std::vector <int> substructureIdx;
            std::vector <int> atomId;
            std::vector <std::string> substructureName;
            std::vector<double> atomicCharge;
            //enum class BondType { SINGLE, DOUBLE, TRIPLE, AMIDE, AROMATIC, DUMMY, UNKNOWN, NOT_CONNECTED };
            //std::vector< BondType> bondTypes;
            std::vector< std::string> bondTypes;
            
            //BondType bondTypeFromString(const std::string& s);
            //std::string bondTypeToString(const BondType& s);
            //std::vector<std::str
        };
        // not implemented
        void write(const std::string fileName, const std::vector<Mol2Data>& data);
        


        void read(const std::string fileName, Mol2Data& data);
        void read(const std::string fileName, MoleculeData& data);
        void atomicNumbers(const Mol2Data& data, std::vector<int>& atomicNumber);
                  

        void write(
            const std::string fileName, 
            const std::vector<int>& atomicNumber, 
            const std::vector<Vector3d>& position, 
            bool detectBonds=false,
            const std::vector<std::string> &atomLabel = std::vector<std::string>(0));

        void write(
            const std::string fileName,
            const std::vector<ChemicalElement>& chemicalElement,
            const std::vector<Vector3d>& position,
            bool detectBonds = false,
            const std::vector<std::string>& atomLabel = std::vector<std::string>(0));


    }
    /**@}*/
}