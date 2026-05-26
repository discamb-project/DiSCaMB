#include "discamb/Scattering/FragHarMacromol.h"
#include "discamb/Scattering/disordered_structure_fragments.h"

#include <set>

using namespace std;

namespace discamb {

    FragHarMacromol::FragHarMacromol(
        const Crystal& crystal,
        const HirshfeldAtomModelSettings& settings,
        const MacromolecularStructuralInformation& macromolInfo,
        bool electronScattering,
        const std::string& jobname)
    {
        set(crystal, settings, macromolInfo, electronScattering, jobname);
    }

    FragHarMacromol::FragHarMacromol(
        const Crystal& crystal,
        const nlohmann::json& data)
    {
        HirshfeldAtomModelSettings settings;
        settings.set(data, crystal);
        string jobName = data.value("job name", "noname");
        bool electronScattering = data.value("electron_scattering", false);
        electronScattering = data.value("electron scattering", electronScattering);
        MacromolecularStructuralInformation mmInfo;
        if (data.find("macromolecular structural information") != data.end())
            mmInfo.set(data["macromolecular structural information"]);
        else
            on_error::throwException("macromolecular information is missing in the json data", __FILE__, __LINE__);

        set(crystal, settings, mmInfo, electronScattering);
    }

    FragHarMacromol::~FragHarMacromol()
    {

    }

    void FragHarMacromol::getModelInformation(
        std::vector<std::pair<std::string, std::string> >& modelInfo)
        const
    {
        //mModelInfo.push_back({ "SCATTERING MODEL", "HAR" });
        
        mSfCalculator->getModelInformation(modelInfo);
        for (auto& p : modelInfo)
            if (p.first == "SCATTERING MODEL")
                p.second = "fragHAR_MX";
    }

    void FragHarMacromol::setAnomalous(
        const std::vector<std::complex<double> >& anomalous)
    {
        mSfCalculator->setAnomalous(anomalous);
    }


    // ignores count atom contrib if iam atoms present
    void FragHarMacromol::calculateStructureFactorsAndDerivatives(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        std::vector<TargetFunctionAtomicParamDerivatives>& dTarget_dparam,
        const std::vector<std::complex<double> >& dTarget_df,
        const std::vector<bool>& countAtomContribution)
    {
        mSfCalculator->calculateStructureFactorsAndDerivatives(
            atoms,
            hkl,
            f,
            dTarget_dparam,
            dTarget_df,
            countAtomContribution);
    }

    // ignores count atom contrib if iam atoms present
    void FragHarMacromol::calculateStructureFactors(
        const std::vector<AtomInCrystal>& atoms,
        const std::vector<Vector3i>& hkl,
        std::vector<std::complex<double> >& f,
        const std::vector<bool>& countAtomContribution)
    {
        mSfCalculator->calculateStructureFactors(
            atoms,
            hkl,
            f,
            countAtomContribution);
    }



    void FragHarMacromol::update(
        const std::vector<AtomInCrystal>& atoms)
    {
        mSfCalculator->update(atoms);
    }
        

    void FragHarMacromol::calculateStructureFactorsAndDerivatives(
        const Vector3i& hkl,
        std::complex<double>& scatteringFactor,
        discamb::SfDerivativesAtHkl& derivatives,
        const std::vector<bool>& countAtomContribution)
    {
        mSfCalculator->calculateStructureFactorsAndDerivatives(
            hkl,
            scatteringFactor,
            derivatives,
            countAtomContribution);
    }

    void FragHarMacromol::calculateFormFactors(
        const Vector3i& hkl,
        std::vector<std::complex<double> >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mSfCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    void FragHarMacromol::calculateFormFactors(
        const std::vector<Vector3i>& hkl,
        std::vector< std::vector<std::complex<double> > >& formFactors,
        const std::vector<bool>& includeAtom)
        const
    {
        mSfCalculator->calculateFormFactors(hkl, formFactors, includeAtom);
    }

    //std::unique_ptr<StockholderAtomSfCalculator> mSfCalculator;
    void FragHarMacromol::set(
        const Crystal& crystal,
        const HirshfeldAtomModelSettings& _settings,
        const MacromolecularStructuralInformation& macromolInfo,
        bool electronScattering,
        const std::string& jobname)
    {
        HirshfeldAtomModelSettings settings = _settings;
        
        make_fragments(crystal, macromolInfo, settings.crystalFragments);
        
        //settings.
        /*
        vector< vector<pair<string, double> > > ordered_parts;
        vector< vector<pair<string, double> > > ordered_parts_atom_indices;
        vector<StructureWithDescriptors> structureDescriptors;

        disordered_structure_fragments::split_and_describe_with_macromol_info_asymm(crystal, macromolInfo, ordered_parts, structureDescriptors);
        
        // alternatively, without structureDescriptors which includes connectivity: 
        // disordered_structure_fragments::split_with_macromol_info(crystal, macromolInfo, ordered_parts);
        
        bool hasConnectivity = false;
        for (auto& descriptors : structureDescriptors)
            if (!descriptors.connectivity.empty())
                hasConnectivity = true;

        map<string, int> label2Idx;
        for (int i = 0; i < crystal.atoms.size(); i++)
            label2Idx[crystal.atoms[i].label] = i;
        //ordered_parts_atom_indices

        if (!hasConnectivity)
        {
            // calculate connectity here
        }
        int nOrderedParts = ordered_parts.size();
        for (int partIdx = 0; partIdx < nOrderedParts; partIdx++)
        {
            MacromolecularStructuralInformation macromolInfoPart;

        }
        */
        //const Crystal &crystal, const HirshfeldAtomModelSettings& settings, bool electronScattering=false, const std::string& jobname = "ham"
        mSfCalculator = std::make_unique<StockholderAtomSfCalculator>(crystal, settings, electronScattering, jobname);
    }

    void FragHarMacromol::make_fragments(
        const Crystal& crystal,
        const MacromolecularStructuralInformation& macromolInfo,
        std::vector<QmFragmentInCrystal>& fragments)
    {
        fragments.clear();
        vector<vector<int> > residueAtoms;
        vector<bool> isAminoAcid;
        int nResidues = 0;
        int nAtoms = macromolInfo.atomNames.size();
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            if (macromolInfo.residueSequenceNumbers[atomIdx] > nResidues)
                nResidues = macromolInfo.residueSequenceNumbers[atomIdx];

        residueAtoms.resize(nResidues);
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            residueAtoms[macromolInfo.residueSequenceNumbers[atomIdx] - 1].push_back(atomIdx);

        for (int resIdx = 0; resIdx < nResidues; resIdx++)
        {
            QmFragmentInCrystal fragment;
            fragment.label = "res_" + to_string(resIdx);

            // add atoms from this residue

            for (int i = 0; i < residueAtoms[resIdx].size(); i++)
            {
                int atomIdx = residueAtoms[resIdx][i];
                fragment.atoms.atomList.push_back({ crystal.atoms[atomIdx].label, "X,Y,Z" });
            }

            // add atoms from previous residue

            if (resIdx > 0)
            {
                int previuos_res = resIdx - 1;
                // C O  CA HA (for gly HA2 HA3) H@CA->CB (if there is CB) H@CA->N
                std::set<string> atomTags;
                bool is_glycine = macromolInfo.residueNames[residueAtoms[previuos_res][0]] == "GLY";
                if(is_glycine)
                    atomTags = { "C","O","CA","HA2", "HA3" };
                else
                    atomTags = { "C", "O", "CA", "HA" };
                
                vector<int> atomIndices;
                int idx_CA = -1;
                int idx_CB = -1;
                int idx_N = -1;

                for (int i = 0; i < residueAtoms[previuos_res].size(); i++)
                {
                    int atomIdx = residueAtoms[previuos_res][i];

                    if (macromolInfo.atomNames[atomIdx] == "CA")
                        idx_CA = atomIdx;

                    if (macromolInfo.atomNames[atomIdx] == "CB")
                        idx_CB = atomIdx;

                    if (macromolInfo.atomNames[atomIdx] == "N")
                        idx_N = atomIdx;

                    if (atomTags.count(macromolInfo.atomNames[atomIdx]) == 1)
                        fragment.atoms.atomList.push_back({ crystal.atoms[atomIdx].label, "X,Y,Z" });
                }

                if (idx_CA == -1)
                    on_error::throwException("missing CA atom", __FILE__, __LINE__);
                if (idx_N == -1)
                    on_error::throwException("missing N atom", __FILE__, __LINE__);
                CappingHydrogen cappingHydrogen;
                cappingHydrogen.bondedAtom = crystal.atoms[idx_CA].label;
                cappingHydrogen.bondedAtomSymmOp = "X,Y,Z";
                cappingHydrogen.directingAtom = crystal.atoms[idx_N].label;
                cappingHydrogen.directingAtomSymmOp = "X,Y,Z";
                fragment.atoms.cappingHydrogens.push_back(cappingHydrogen);
                if(idx_CB == -1 &&  !is_glycine)
                    on_error::throwException("missing CB atom", __FILE__, __LINE__);
                if (idx_CB != -1)
                {
                    cappingHydrogen.directingAtom = crystal.atoms[idx_CB].label;
                    cappingHydrogen.directingAtomSymmOp = "X,Y,Z";
                    fragment.atoms.cappingHydrogens.push_back(cappingHydrogen);
                }
            }

            // add atoms from next residue

            if (resIdx < nResidues - 1)
            {
                int next_res = resIdx + 1;
                // N H CA HA (for gly HA2 HA3) H@CA->CB (if there is CB) H@CA->C
                std::set<string> atomTags;
                bool is_glycine = macromolInfo.residueNames[residueAtoms[next_res][0]] == "GLY";
                if (is_glycine)
                    atomTags = { "N","H","CA","HA2", "HA3" };
                else
                    atomTags = { "N", "H", "CA", "HA" };

                vector<int> atomIndices;
                int idx_CA = -1;
                int idx_CB = -1;
                int idx_C = -1;

                for (int i = 0; i < residueAtoms[next_res].size(); i++)
                {
                    int atomIdx = residueAtoms[next_res][i];

                    if (macromolInfo.atomNames[atomIdx] == "CA")
                        idx_CA = atomIdx;

                    if (macromolInfo.atomNames[atomIdx] == "CB")
                        idx_CB = atomIdx;

                    if (macromolInfo.atomNames[atomIdx] == "C")
                        idx_C = atomIdx;

                    if (atomTags.count(macromolInfo.atomNames[atomIdx]) == 1)
                        fragment.atoms.atomList.push_back({ crystal.atoms[atomIdx].label, "X,Y,Z" });
                }

                if (idx_CA == -1)
                    on_error::throwException("missing CA atom", __FILE__, __LINE__);
                if (idx_C == -1)
                    on_error::throwException("missing C atom", __FILE__, __LINE__);
                CappingHydrogen cappingHydrogen;
                cappingHydrogen.bondedAtom = crystal.atoms[idx_CA].label;
                cappingHydrogen.bondedAtomSymmOp = "X,Y,Z";
                cappingHydrogen.directingAtom = crystal.atoms[idx_C].label;
                cappingHydrogen.directingAtomSymmOp = "X,Y,Z";
                fragment.atoms.cappingHydrogens.push_back(cappingHydrogen);
                if (idx_CB == -1 && !is_glycine)
                    on_error::throwException("missing CB atom", __FILE__, __LINE__);
                if (idx_CB != -1)
                {
                    cappingHydrogen.directingAtom = crystal.atoms[idx_CB].label;
                    cappingHydrogen.directingAtomSymmOp = "X,Y,Z";
                    fragment.atoms.cappingHydrogens.push_back(cappingHydrogen);
                }
            }

            fragments.push_back(fragment);
        }
    }
}