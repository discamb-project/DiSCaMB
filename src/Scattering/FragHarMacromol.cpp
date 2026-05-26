#include "discamb/Scattering/FragHarMacromol.h"
#include "discamb/Scattering/disordered_structure_fragments.h"

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
        const HirshfeldAtomModelSettings& settings,
        const MacromolecularStructuralInformation& macromolInfo,
        bool electronScattering,
        const std::string& jobname)
    {
        vector<QmFragmentInCrystal> fragments;

        make_fragments(macromolInfo, fragments);

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

        //mSfCalculator = std::make_unique<StockholderAtomSfCalculator>();
    }

    void FragHarMacromol::make_fragments(
        const MacromolecularStructuralInformation& macromolInfo,
        std::vector<QmFragmentInCrystal>& fragments)
    {

    }
}