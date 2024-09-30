#pragma once

#include "discamb/QuantumChemistry/ElectronDensityPartition.h"
#include "discamb/QuantumChemistry/ElectronDensityCalculator.h"
#include "discamb/QuantumChemistry/PromoleculeElectronDensity.h"



#include <memory>


namespace discamb {

    /**
    * \addtogroup QuantumChemistry
    * @{
    */


    class MbisPartition: public ElectronDensityPartition
    {
    public:
        /*
        if (data.find("MBIS convergence") != data.end())
            mConvergenceThreshold = data.find("MBIS convergence").value().get<double>();
        if (data.find("MBIS n iterations") != data.end())
            mMaxN_Iterations = data.find("MBIS n iterations").value().get<int>();
        if (data.find("n cores") != data.end())
            mN_Cores = data.find("n cores").value().get<int>();

        */
        
        MbisPartition();

        MbisPartition(
            std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
            const std::vector<int> &atomicNumbers,
            const std::vector<Vector3d>& positions,
            const nlohmann::json& data = nlohmann::json());

        MbisPartition(const std::string& wfnFile, const nlohmann::json& data = nlohmann::json());

        MbisPartition(
            std::shared_ptr<ElectronDensityCalculator>& ed_calculator,
            const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& positions,
            double convergenceThreshold = 1e-6,
            int maxN_Steps = 100,
            int nCores = 1,
            bool parallelAfterInitialization = false);

        MbisPartition(
            const std::string& wfnFile,
            double convergenceThreshold,
            int maxN_Steps = 100,
            int nCores = 1,
            bool parallelAfterInitialization = false);


        virtual ~MbisPartition() override;

        virtual void set(
            const nlohmann::json& data, 
            const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& atomPositions,
            std::shared_ptr<ElectronDensityCalculator>& edCalculator);

        void set(std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
            const std::vector<int>& atomicNumbers,
            const std::vector<Vector3d>& positions);
        

        void setFromFile(const std::string& wfn);
        virtual void applySettings(const nlohmann::json& data);
        void applySettings(double convergenceThreshold = 1e-6,
            int maxN_Steps = 100,
            int nCores = 1,
            bool parallelAfterInitialization = false);

        virtual double calculate(int atom, const Vector3d &r) const;
        virtual double calculate(int atom, const Vector3d& r,int molecularOrbital) const override;
        double promoleculeDensity(const Vector3d& r) const;
        double promoleculeDensity(int atomIdx, const Vector3d& r) const;
        double promoleculeDensity(int n, int atomIdx, const Vector3d& r) const;
        virtual void setIncludeRange(double range) override;
        void unsetIncludeRange();
        void getParameters(std::vector<std::vector<double> >& populations, std::vector<std::vector<double> >& exponentsInverted) const;
        void setParameters(const std::vector<std::vector<double> >& populations, const std::vector<std::vector<double> >& exponentsInverted);
        void getAtomicCharges(std::vector<double>& charges) const;
        virtual void setAtom(int atomIdx) override;
        virtual double calculate(const Vector3d& r) const { return calculate(mAtomIdx, r); };
        MbisPartition* clone() const override;
        void setPower(double p);
        void initialize() const;
    private:
        void readParametersFile(const std::string& fileName) const;
        void saveParametersFile(const std::string& fileName) const;
        std::string mParametersInputFile;
        std::string mParametersOutputFile;
        mutable bool mReadParametersFromFile = false;
        mutable bool mSaveParametersFile = false;
        mutable bool mPrintError = false;
        bool mParallelAfterInitialization = false;
        mutable int mN_Cores = 1;
        int mMaxN_Iterations=100;
        double mConvergenceThreshold = 1e-6;
        int mN_Atoms = 0;
        int mAtomIdx = 0;
        double mPower = 1.0;
        double mIncludeRange = -1;
        mutable bool mIsInitialized = false;
        std::vector<std::vector<int> > mIncludeAtoms;
        
        
        std::shared_ptr<ElectronDensityCalculator> mEdCalculator;
        mutable std::vector<std::vector<double> > mExponentsIverted;
        mutable std::vector<std::vector<double> > mPopulations;
        // total electron density multiplied by interation weights
        // on atomic integration grids
        mutable std::vector<std::vector<double> > mElectronDensityGrids, mAtomicDensityGrids, mNextAtomicDensityGrids;
        mutable std::vector<std::vector<double> > mPromoleculeDensity;
        mutable std::vector<std::vector<std::vector<double> > > mShellWeights;
        void makeElectronDensityGrids() const;
        void setShellWeigths() const;
        void allocateShellWeights() const;
        void allocatePromoleculeDensity() const;
        void makeAtomicDensityGrids() const;
        void calculatePromoleculeDensities() const;
        void calculatePromoleculeDensities2() const;
        void calculateAtomicDensities(std::vector<std::vector<double> > &atomicDensities)  const;
        void initialGuess() const;
        void setIntegrationGrid() const;
        void findParameters() const;
        void nextParameters(
            int atomIdx,
            int n,
            double& newExponentInverted,
            double& newPopulation) const;

        void nextParameters(
            std::vector<std::vector<double> >& newExponentsIverted,
            std::vector<std::vector<double> >& newPopulations)  const;

        std::vector<int> mAtomicNumbers;
        std::vector<Vector3d> mPositions;


        mutable std::vector < std::vector<Vector3d> > mElementGridPoints;
        mutable std::vector < std::vector<double> > mElementIntegrationWeights;

        //double compareWithOldParameters();
            //std::vector<std::vector<double> >& exponentsIverted,
            //std::vector<std::vector<double> >& populations);

        double compareWithOldParameters(
            int atomIdx) const;/*,
            std::vector<std::vector<double> >& exponentsIverted,
            std::vector<std::vector<double> >& populations);*/


    };
    /**@}*/
}