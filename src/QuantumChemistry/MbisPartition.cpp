#include "discamb/QuantumChemistry/MbisPartition.h"
#include "discamb/BasicUtilities/on_error.h"
#include "discamb/BasicChemistry/periodic_table.h"
#include "discamb/IO/wfn_io.h"
#include "discamb/MathUtilities/math_utilities.h"
#include "discamb/MathUtilities/lebedev_laikov.h"
#include "discamb/MathUtilities/radial_grid.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <omp.h>

using namespace std;

namespace discamb {

    MbisPartition::MbisPartition()
    {}

    MbisPartition::MbisPartition(
        std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& positions,
        const nlohmann::json& data)
    {
        applySettings(data);
        set(ed_calculator, atomicNumbers, positions);
    }

    MbisPartition::MbisPartition(
        std::shared_ptr<ElectronDensityCalculator>& ed_calculator,
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& positions,
        double convergenceThreshold,
        int maxN_Steps,
        int nCores,
        bool parallelAfterInitialization)
    {
        applySettings(convergenceThreshold,
            maxN_Steps,
            nCores,
            parallelAfterInitialization);
        set(ed_calculator, atomicNumbers, positions);
    }


    MbisPartition::MbisPartition(
        const std::string& wfnFile,
        const nlohmann::json& data)
    {
        applySettings(data);
        setFromFile(wfnFile);
    }

    MbisPartition::MbisPartition(const std::string& wfnFile,
        double convergenceThreshold,
        int maxN_Steps,
        int nCores,
        bool parallelAfterInitialization)
    {
        applySettings(convergenceThreshold,
            maxN_Steps,
            nCores,
            parallelAfterInitialization);

        setFromFile(wfnFile);
    }

    void MbisPartition::set(
        const nlohmann::json& data,
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& atomPositions,
        std::shared_ptr<ElectronDensityCalculator>& edCalculator)
    {
        set(edCalculator, atomicNumbers, atomPositions);
        applySettings(data);
    }

    void MbisPartition::applySettings(
        double convergenceThreshold,
        int maxN_Steps,
        int nCores,
        bool parallelAfterInitialization)
    {
        mParallelAfterInitialization = parallelAfterInitialization;
        mN_Cores = nCores;
        mMaxN_Iterations = maxN_Steps;
        mConvergenceThreshold = convergenceThreshold;
    }


    MbisPartition::~MbisPartition()
    {
    }

    void MbisPartition::setParameters(
        const std::vector<std::vector<double> >& populations,
        const std::vector<std::vector<double> >& exponentsInverted)
    {
        mPopulations = populations;
        mExponentsIverted = exponentsInverted;
        
    }

    void MbisPartition::getParameters(
        std::vector<std::vector<double> >& populations,
        std::vector<std::vector<double> >& exponentsInverted)
        const
    {
        populations = mPopulations;
        exponentsInverted = mExponentsIverted;
    }

    //std::vector<std::vector<double> > mElectronDensityGrids;
    void MbisPartition::makeElectronDensityGrids()
        const
    {
        int nAtoms = static_cast<int>(mAtomicNumbers.size());
        mElectronDensityGrids.resize(nAtoms);
        cout << "n threads = " << mN_Cores << "\n";
        omp_set_num_threads(static_cast<int>(mN_Cores));
        
        vector<int> threadIdx(nAtoms);
        
        vector<shared_ptr<ElectronDensityCalculator> > ed_calcs(mN_Cores);
        ed_calcs[0] = mEdCalculator;
        for (int i = 1; i < mN_Cores; i++)
            ed_calcs[i] = shared_ptr<ElectronDensityCalculator>(mEdCalculator->clone());


#pragma omp parallel for
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            //cout << omp_get_thread_num() << "\n";
            int threadId = omp_get_thread_num();
            threadIdx[atomIdx] = threadId;
            int z = mAtomicNumbers[atomIdx];
            int pointIdx, nPoints = mElementGridPoints[z].size();
            mElectronDensityGrids[atomIdx].resize(nPoints);
            vector<Vector3d> &rGrid = mElementGridPoints[z];
            vector<double> &w = mElementIntegrationWeights[z];
            Vector3d const & r0 = mPositions[atomIdx];
            for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
            {
                Vector3d r = r0 + rGrid[pointIdx];
                mElectronDensityGrids[atomIdx][pointIdx] = ed_calcs[threadId]->calculate2(r.x, r.y, r.z); //mEdCalculator->calculate2(r.x, r.y, r.z);
            }
        }

    }

    void MbisPartition::applySettings(
        const nlohmann::json& data)
    {
        
        mConvergenceThreshold = data.value("MBIS convergence", 1e-6);
        mMaxN_Iterations = data.value("MBIS n iterations", 1000);
        mPrintError = data.value("print error", false);

        if (data.find("n cores") != data.end())
            mN_Cores = data.find("n cores").value().get<int>();
        else
        {
            cout << "no n cores" << "\n";
            exit(0);
        }
        
        if (data.find("read parameters file") != data.end())
        {
            mReadParametersFromFile = true;
            mParametersInputFile = data["read parameters file"].get<string>();
        }

        if (data.find("write parameters file") != data.end())
        {
            mSaveParametersFile = true;
            mParametersOutputFile = data["write parameters file"].get<string>();
        }

    }

    void MbisPartition::calculateAtomicDensities(
        std::vector<std::vector<double> >& atomicDensities)
        const 
    {
        int atomIdx, nAtoms = mAtomicDensityGrids.size();

        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int z = mAtomicNumbers[atomIdx];
            const vector<double>& rhoMol = mElectronDensityGrids[atomIdx];
            const vector<Vector3d>& r = mElementGridPoints[z];
            Vector3d r0 = mPositions[atomIdx];
            int i, nPoints = r.size();
            auto& rhoAtom = mAtomicDensityGrids[atomIdx];
            for (i = 0; i < nPoints; i++)
            {
                Vector3d rGlobal = r[i] + r0;
                rhoAtom[i] = rhoMol[i] * promoleculeDensity(atomIdx, rGlobal) / mPromoleculeDensity[atomIdx][i];// promoleculeDensity(rGlobal);
            }

        }

    }

    void MbisPartition::getAtomicCharges(
        std::vector<double>& charges)
        const
    {

        int atomIdx, nAtoms = mAtomicNumbers.size();
        
        charges.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int z = mAtomicNumbers[atomIdx];
            charges[atomIdx] = double(mAtomicNumbers[atomIdx]);
            const vector<double>& w = mElementIntegrationWeights[z];
            const vector<double>& rhoMol = mElectronDensityGrids[atomIdx];
            const vector<Vector3d>& r = mElementGridPoints[z];
            Vector3d r0 = mPositions[atomIdx];
            int i, nPoints = r.size();
            double rhoAtom;
            for (i = 0; i < nPoints; i++)
            {
                Vector3d rGlobal = r[i] + r0;
                //rhoAtom = mEdCalculator->calculate2(rGlobal.x, rGlobal.y, rGlobal.z) * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
                rhoAtom = rhoMol[i] * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
                charges[atomIdx] -= w[i] * rhoAtom;
            }

        }
    }


    void MbisPartition::setIncludeRange(double range)
    {
        mIncludeRange = range;
        

        int i, j, nAtoms = mPositions.size();
        mIncludeAtoms.clear();
        mIncludeAtoms.resize(nAtoms);

        if (range < 0)
        {
            for (i = 0; i < nAtoms; i++)
                mIncludeAtoms[0].push_back(i);
            for (i = 1; i < nAtoms; i++)
                mIncludeAtoms[i] = mIncludeAtoms[0];
            return;
        }

        double length2, range2 = range*range;
        Vector3d r;
        for (i = 0; i < nAtoms; i++)
        {
            mIncludeAtoms[i].push_back(i);
            for (j = 0; j < i; j++)
            {
                r = mPositions[i] - mPositions[j];
                length2 = r * r;
                if (length2 <= range2)
                {
                    mIncludeAtoms[i].push_back(j);
                    mIncludeAtoms[j].push_back(i);
                }
            }
        }

        mEdCalculator->setContributingCenters(mIncludeAtoms[mAtomIdx]);
    }
    
    void MbisPartition::unsetIncludeRange()
    {
        mIncludeRange = -1;
    }


    void MbisPartition::setFromFile(
        const std::string& wfnFile)
    {
        wfn_io::WfnFileData wfn;
        wfn_io::read_wavefunction(wfnFile, wfn);
        mN_Atoms = wfn.center_position.size();
        auto edCalculator = std::shared_ptr<ElectronDensityCalculator>(new ElectronDensityCalculator);
        
        edCalculator->set_1rdm(
            wfn.center_position,
            wfn.primitive_to_center,
            wfn.primitive_type,
            wfn.primitive_exponents,
            wfn.molecular_orbital_occupancy,
            wfn.molecular_orbitals);

        if (!wfn.edfs.empty())
            edCalculator->setAdditionalDensity(
                wfn.edfs[0].primitive_to_center,
                wfn.edfs[0].primitive_type,
                wfn.edfs[0].primitive_exponents,
                wfn.edfs[0].primitive_coefficients);

        set(edCalculator, wfn.atomic_numbers, wfn.center_position);        
    }


    void MbisPartition::setPower(double p)
    {
        mPower = p;
    }

    MbisPartition* MbisPartition::clone()
        const
    {
        shared_ptr<ElectronDensityCalculator> ed_calculator(mEdCalculator->clone());
        MbisPartition* mbisPartition = new MbisPartition();
        //mbisPartition->mEdCalculator = ed_calculator;
        mbisPartition->mEdCalculator = shared_ptr<ElectronDensityCalculator>( ed_calculator->clone());
        mbisPartition->mAtomicNumbers = mAtomicNumbers;
        mbisPartition->mAtomIdx = mAtomIdx;
        mbisPartition->mExponentsIverted = mExponentsIverted;
        mbisPartition->mIncludeAtoms = mIncludeAtoms;
        mbisPartition->mIncludeRange = mIncludeRange;
        mbisPartition->mPopulations = mPopulations;
        mbisPartition->mPositions = mPositions;
        mbisPartition->mMaxN_Iterations = mMaxN_Iterations;
        mbisPartition->mConvergenceThreshold = mConvergenceThreshold;
        mbisPartition->mIsInitialized = mIsInitialized;

        return mbisPartition;
    }

    double MbisPartition::promoleculeDensity(
        const Vector3d& r)
        const
    {
        double rho0 = 0;
        for (int atomIdx = 0; atomIdx < mAtomicNumbers.size(); atomIdx++)
            rho0 += promoleculeDensity(atomIdx, r);
        return rho0;
    }

    double MbisPartition::promoleculeDensity(
        int atomIdx,
        const Vector3d& r)
        const
    {
        double rho0 = 0;
        for (int n = 1; n <= mExponentsIverted[atomIdx].size(); n++)
            rho0 += promoleculeDensity(n, atomIdx, r);
        return rho0;
    }

    double MbisPartition::promoleculeDensity(int n, int atomIdx, const Vector3d& r)
        const
    {
        double a = mExponentsIverted[atomIdx][n-1];
        Vector3d v = r - mPositions[atomIdx];
        double rLength = sqrt(v * v);
        return mPopulations[atomIdx][n-1] / (8.0 * M_PI * a * a * a) * exp(-rLength / a);
    }

    double MbisPartition::compareWithOldParameters(
        int atomIdx//,
        //std::vector<std::vector<double> >& exponentsIverted,
        //std::vector<std::vector<double> >& populations
    )
        const
    {
        int z = mAtomicNumbers[atomIdx];
        vector<double> diff = mElementIntegrationWeights[z];

        double result = 0;
        
        vector<double>& w = mElementIntegrationWeights[z];
        //vector<Vector3d>& r = mElementGridPoints[z];
        //auto& rhoMol = mElectronDensityGrids[atomIdx];
        
        //Vector3d r0 = mPositions[atomIdx];
        int i, nPoints = mAtomicDensityGrids[atomIdx].size();
        //vector<double> rhoOld(nPoints);
        
        //for (i = 0; i < nPoints; i++)
        //{
        //    Vector3d rGlobal = r[i] + r0;
        //    
        //    //rhoOld[i] = mEdCalculator->calculate2(rGlobal.x, rGlobal.y, rGlobal.z) * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
        //    rhoOld[i] = rhoMol[i] * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
        //}

        //exponentsIverted.swap(mExponentsIverted);
        //populations.swap(mPopulations);

        //for (i = 0; i < nPoints; i++)
        //{
        //    Vector3d rGlobal = r[i] + r0;
        //    //double rho = mEdCalculator->calculate2(rGlobal.x, rGlobal.y, rGlobal.z) * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
        //    double rho = rhoMol[i] * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
        //    result += w[i] * (rho - rhoOld[i]) * (rho - rhoOld[i]);
        //}

        for (i = 0; i < nPoints; i++)
        {
            //Vector3d rGlobal = r[i] + r0;
            //double rho = mEdCalculator->calculate2(rGlobal.x, rGlobal.y, rGlobal.z) * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
            //double rho = rhoMol[i] * promoleculeDensity(atomIdx, rGlobal) / promoleculeDensity(rGlobal);
            //result += w[i] * (rho - rhoOld[i]) * (rho - rhoOld[i]);
            double d = mAtomicDensityGrids[atomIdx][i] - mNextAtomicDensityGrids[atomIdx][i];
            result += d * d * w[i];
        }


        //exponentsIverted.swap(mExponentsIverted);
        //populations.swap(mPopulations);

        return sqrt(result);
    }


    //double MbisPartition::compareWithOldParameters(
    //    std::vector<std::vector<double> >& exponentsIverted,
    //    std::vector<std::vector<double> >& populations)
    //{
    //    int atomIdx, nAtoms = mPositions.size();
    //    vector<double> diff(nAtoms);

    //    for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
    //        diff[atomIdx] = compareWithOldParameters(atomIdx, exponentsIverted, populations);
    //    return *max_element(diff.begin(), diff.end());
    //}


    void MbisPartition::allocateShellWeights()
        const
    {
        int atomIdx, nPoints, nAtoms = mAtomicNumbers.size();
        mShellWeights.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            nPoints = mElementGridPoints[mAtomicNumbers[atomIdx]].size();
            mShellWeights[atomIdx].resize(mPopulations[atomIdx].size(),vector<double>(nPoints));
        }
    }

  /*  void MbisPartition::setShellWeigths()
    {
        int atomIdx, shellIdx, pointIdx, nPoints, nShells, nAtoms = mAtomicNumbers.size();
        
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int z = mAtomicNumbers[atomIdx];
            nShells = mPopulations[atomIdx].size();
            nPoints = mElementGridPoints[z].size();
            for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
            for(shellIdx=0; shellIdx<nShells; shellIdx++)
                
                    mShellWeights[atomIdx][shellIdx][.resize(mPopulations[atomIdx].size(), vector<double>(nPoints));
        }
    }*/

    void MbisPartition::findParameters()
        const
    {

        if (mReadParametersFromFile)
        {
            mReadParametersFromFile = false;
            readParametersFile(mParametersInputFile);
            return;
        }

        /*
        int mMaxN_Iterations=100;
        double mConvergenceThreshold = 1e-6;

        */
        //int max_n_iterations = 20;
        //double threshold = 0.001;
        double diff = 2 * mConvergenceThreshold;

        int nAtoms = mAtomicNumbers.size();

        vector<vector<double> > newExponentsIverted, newPopulations;
        int iteration = 1;
        // for memory allocation only
        mAtomicDensityGrids = mElectronDensityGrids;
        mNextAtomicDensityGrids = mElectronDensityGrids;
        mPromoleculeDensity = mElementIntegrationWeights;
        // ----
        allocatePromoleculeDensity();
        calculatePromoleculeDensities();
        calculateAtomicDensities(mAtomicDensityGrids);
        vector<double> atomicCharges;
        cout<< "MBIS convergence indicator, threshold set to " << setw(14) << setprecision(8) << fixed << mConvergenceThreshold << "\n";
        double minDiff = 1e+15;
        
        ofstream out;
        if (mPrintError)
            out.open("error.mbis");
        while (diff > mConvergenceThreshold && iteration < mMaxN_Iterations)
        {
            nextParameters(newExponentsIverted, newPopulations);
            mPopulations = newPopulations;
            mExponentsIverted = newExponentsIverted;
            calculatePromoleculeDensities();
            calculateAtomicDensities(mNextAtomicDensityGrids);

            vector<double> atomicDifferences(nAtoms);
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
                atomicDifferences[atomIdx] = compareWithOldParameters(atomIdx);// , newExponentsIverted, newPopulations);

            mNextAtomicDensityGrids.swap(mAtomicDensityGrids);

            diff = *max_element(atomicDifferences.begin(), atomicDifferences.end());
            getAtomicCharges(atomicCharges);
            //diff = compareWithOldParameters(newExponentsIverted, newPopulations);
            cout << setw(14) << setprecision(8) << fixed << "\r" << diff;
            if (mPrintError)
            {
                minDiff = min(minDiff, diff);
                out << setw(14) << setprecision(8) << fixed << diff
                    << setw(14) << setprecision(8) << fixed << minDiff << endl;
            }
            //for (auto q : atomicCharges)
            //    cout << setw(14) << setprecision(6) << fixed << q;
            //cout << endl;
            mExponentsIverted.swap(newExponentsIverted);
            mPopulations.swap(newPopulations);
            iteration++;
        }
        cout << endl;
        if (mPrintError)
            out.close();
        if (mSaveParametersFile)
            saveParametersFile(mParametersOutputFile);

    }

    void MbisPartition::saveParametersFile(
        const std::string& fileName) const
    {
        ofstream out(fileName);

        if (!out.good())
            on_error::throwException("cannot write MBIS partition parameters to file '" + fileName + "'", __FILE__, __LINE__);

        int nAtoms = mAtomicNumbers.size();

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int n = mExponentsIverted[atomIdx].size();
            for (int i = 0; i < n; i++)
                out << setprecision(12) << mPopulations[atomIdx][i] << " ";
            out << "\n";
            for (int i = 0; i < n; i++)
                out << setprecision(12) << mExponentsIverted[atomIdx][i] << " ";
            out << "\n";
        }
        out.close();
    }

    void MbisPartition::readParametersFile(
        const std::string& fileName)
        const
    {
        ifstream in(fileName);

        if (!in.good())
            on_error::throwException("cannot read MBIS partition parameters from file '" + fileName + "'", __FILE__, __LINE__);

        int nAtoms = mAtomicNumbers.size();

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int n = mExponentsIverted[atomIdx].size();
            for (int i = 0; i < n; i++)
                in >> mPopulations[atomIdx][i];
            for (int i = 0; i < n; i++)
                in>> mExponentsIverted[atomIdx][i];
        }
        in.close();

    }

    void MbisPartition::nextParameters(
        std::vector<std::vector<double> > & newExponentsIverted,
        std::vector<std::vector<double> > & newPopulations)
        const
    {
        int nAtoms = static_cast<int>(mAtomicNumbers.size());

        newExponentsIverted.resize(nAtoms);
        newPopulations.resize(nAtoms);
        omp_set_num_threads(int(mN_Cores));
#pragma omp parallel for
        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int n = periodic_table::period(mAtomicNumbers[atomIdx]);
            
            newExponentsIverted[atomIdx].resize(n);
            newPopulations[atomIdx].resize(n);

            for (int i = 1; i <= n; i++)
                nextParameters(atomIdx, i, newExponentsIverted[atomIdx][i-1], newPopulations[atomIdx][i-1]);
        }

    }


    void MbisPartition::initialGuess()
        const
    {
        
        int nAtoms = mAtomicNumbers.size();

        vector<int> nElementsPerShell{ 2, 8, 8, 18, 18, 32, 32 };

        mExponentsIverted.resize(nAtoms);
        mPopulations.resize(nAtoms);

        for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            double z = double(mAtomicNumbers[atomIdx]);
            int n = periodic_table::period(mAtomicNumbers[atomIdx]);
            mExponentsIverted[atomIdx].resize(n);
            mPopulations[atomIdx].resize(n);

            // population
            int nElectrons = 0;
            for (int i = 0; i < n - 1; i++)
            {
                mPopulations[atomIdx][i] = nElementsPerShell[i];
                nElectrons += nElementsPerShell[i];
            }
            mPopulations[atomIdx].back() = z - nElectrons;

            // exponents
            mExponentsIverted[atomIdx][0]   = 0.5 / z;
            mExponentsIverted[atomIdx][n-1] = 0.5;

            for (int shellIdx = 2; shellIdx < n; shellIdx++)
                mExponentsIverted[atomIdx][shellIdx - 1] = 0.5 / pow(z, 1 - double(shellIdx - 1) / double(n - 1));

        }
    }

    void MbisPartition::set(
        std::shared_ptr<ElectronDensityCalculator> &ed_calculator,
        const std::vector<int>& atomicNumbers,
        const std::vector<Vector3d>& positions)
    {

        mIsInitialized = false;
        mEdCalculator = ed_calculator;
        mAtomicNumbers = atomicNumbers;
        mPositions = positions;
        setIncludeRange(-1);
        initialize();
    }

    void MbisPartition::initialize()
        const
    {
        mIsInitialized = true;
        initialGuess();
        setIntegrationGrid();
        makeElectronDensityGrids();
        findParameters();
        if (!mParallelAfterInitialization)
            mN_Cores = 1;
    }

    void MbisPartition::calculatePromoleculeDensities2()
        const
    {

        int atomIdx, nAtoms = mAtomicNumbers.size();
        
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int pointIdx, nPoints = mPromoleculeDensity[atomIdx].size();
            int z = mAtomicNumbers[atomIdx];
            auto& rGrid = mElementGridPoints[z];
            Vector3d r0 = mPositions[atomIdx];
            for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
            {
                Vector3d r = r0 + rGrid[pointIdx];
                mPromoleculeDensity[atomIdx][pointIdx] = 0;
                for(int atom=0;atom<nAtoms; atom++)
                    mPromoleculeDensity[atomIdx][pointIdx] += promoleculeDensity(atom, r);
            }
        }

        
    }

    void MbisPartition::calculatePromoleculeDensities()
       const
    {
        int atomIdx, nAtoms = int(mAtomicNumbers.size());
        omp_set_num_threads(int(mN_Cores));
        //cout << mN_Cores << "\n";
        //cout << omp_get_max_threads() << endl;
#pragma omp for
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            //cout << omp_get_max_threads() << endl;
            int pointIdx, nPoints = mPromoleculeDensity[atomIdx].size();
            int z = mAtomicNumbers[atomIdx];
            auto& rGrid = mElementGridPoints[z];
            Vector3d r0 = mPositions[atomIdx];
            for (pointIdx = 0; pointIdx < nPoints; pointIdx++)
            {
                Vector3d r = r0 + rGrid[pointIdx];
                mPromoleculeDensity[atomIdx][pointIdx] = 0;
                for (int atom = 0; atom < nAtoms; atom++)
                    mPromoleculeDensity[atomIdx][pointIdx] += promoleculeDensity(atom, r);
            }
        }


    }


    void MbisPartition::allocatePromoleculeDensity()
        const
    {
        int atomIdx, nAtoms = mAtomicNumbers.size();
        mPromoleculeDensity.resize(nAtoms);
        for (atomIdx = 0; atomIdx < nAtoms; atomIdx++)
        {
            int z = mAtomicNumbers[atomIdx];
            int nPoints = mElementGridPoints[z].size();
            mPromoleculeDensity[atomIdx].resize(nPoints);
        }
    }

    void MbisPartition::nextParameters(
        int atomIdx,
        int n,
        double &newExponentInverted,
        double &newPopulation)
        const
    {
        
        int z = mAtomicNumbers[atomIdx];
        vector<double>& w = mElementIntegrationWeights[z];
        vector<Vector3d>& r = mElementGridPoints[z];
        auto& rhoMol = mElectronDensityGrids[atomIdx];
        Vector3d r0 = mPositions[atomIdx];
        int i, nPoints = r.size();
        double f;
        newPopulation = 0;
        newExponentInverted = 0;
        for (i = 0; i < nPoints; i++)
        {
            Vector3d rGlobal = r[i] + r0;
            //f = mEdCalculator->calculate2(rGlobal.x, rGlobal.y, rGlobal.z) * promoleculeDensity(n, atomIdx, rGlobal) / promoleculeDensity(rGlobal);
            f = rhoMol[i] * promoleculeDensity(n, atomIdx, rGlobal) / mPromoleculeDensity[atomIdx][i];// promoleculeDensity(rGlobal);
            newPopulation += w[i] * f;
            newExponentInverted += w[i] * f * sqrt(r[i]*r[i]);
        }
        newExponentInverted *= 1.0 / (3.0 * mPopulations[atomIdx][n - 1]);
    }

    

    void MbisPartition::setIntegrationGrid()
        const
    {
        int nAngular = 590;
        int nRadial = 75;

        mElementGridPoints.resize(113);
        mElementIntegrationWeights.resize(113); 


        std::set<int> uniqueZ;

        uniqueZ.insert(mAtomicNumbers.begin(), mAtomicNumbers.end());

        vector<Vector3d> angularGridPoints, gridPoints;
        vector<double> radialPoints, radialWeights, angularWeights, weights;
        Vector3d r, atomPosition;

        lebedev_laikov::get_grid(nAngular, angularGridPoints, angularWeights);

        gridPoints.resize(nAngular * nRadial);
        weights.resize(nAngular * nRadial);

        for (auto z : uniqueZ)
        {
            radial_grid::treutler_ahlrichs(z, nRadial, radialPoints, radialWeights);
            int pointIdx = 0;
            for (int angularIdx = 0; angularIdx < nAngular; angularIdx++)
                for (int radialIdx = 0; radialIdx < nRadial; radialIdx++)
                {
                    weights[pointIdx] = radialWeights[radialIdx] * angularWeights[angularIdx] *
                        radialPoints[radialIdx] * radialPoints[radialIdx] * M_PI * 4.0;
                    gridPoints[pointIdx] = radialPoints[radialIdx] * angularGridPoints[angularIdx];
                    pointIdx++;
                }
            mElementGridPoints[z] = gridPoints;
            mElementIntegrationWeights[z] = weights;
        }

    }



    void MbisPartition::setAtom(int atomIdx)
    { 
        mAtomIdx = atomIdx; 
        mEdCalculator->setContributingCenters(mIncludeAtoms[atomIdx]);
    };


    double MbisPartition::calculate(
        int atom,
        const Vector3d& r,
        int molecularOrbitalIdx)
        const
    {

        double promol = promoleculeDensity(r);
        
        if (promol > 0)
        {
            double rho = mEdCalculator->calculate2(r.x, r.y, r.z, molecularOrbitalIdx);
            double promol_atom = promoleculeDensity(atom, r);
            return rho * promol_atom / promol;
        }

        return 0.0;
    }

    double MbisPartition::calculate(
        int atom,
        const Vector3d &r)
        const 
    {

        //mEdCalculator->setContributingCenters(mIncludeAtoms[atom]);

        if (mPower == 1.0)
        {
            double promol = promoleculeDensity(r);
            if (promol > 0)
            {
                double rho = mEdCalculator->calculate2(r.x, r.y, r.z);
                double promol_atom = promoleculeDensity(atom, r);
                return rho * promol_atom / promol;
            }
        }
        else
        {
            int nAtoms = mAtomicNumbers.size();
            double weight_denominator = 0;
            double atomicDensity;
            for (int atomIdx = 0; atomIdx < nAtoms; atomIdx++)
            {
                atomicDensity = promoleculeDensity(atomIdx,r);
                if(atomicDensity>0)
                    weight_denominator += pow(atomicDensity , mPower);
            }
        
            if (weight_denominator > 0)
            {
                double rho = mEdCalculator->calculate2(r.x, r.y, r.z);
                atomicDensity = promoleculeDensity(atom, r);
                if (atomicDensity > 0)
                {
                    double weight_numerator = pow(promoleculeDensity(atom, r), mPower);
                    return rho * weight_numerator / weight_denominator;
                }
            }
        }
            return 0.0;
    }


    
}
