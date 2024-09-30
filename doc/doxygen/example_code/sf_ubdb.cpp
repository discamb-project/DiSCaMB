#include "discamb/IO/xd_io.h"
#include "discamb/BasicUtilities/Exception.h"
#include "discamb/Scattering/HansenCoppensStructureFactorCalculator.h"

#include <iostream>

using namespace std;
using namespace discamb;

void run(const char *xdMasterFileName,const char *xdParameterFileName)
{
    HC_ModelParameters hcModelParameters;
    Crystal crystal;
    vector<Vector3i> hkl_indices;
    vector<complex<double> > scattering_factors;
    vector<XdLocalCoordinateSystem> localCoordinateSystems;
    vector<Matrix3d> localCoordinateSystemMatrices;

    // defines reflections for which the scattering factors will be calculated
    hkl_indices.push_back(Vector3i(0, 1, 0));
    hkl_indices.push_back(Vector3i(1, 1, 0));
    hkl_indices.push_back(Vector3i(2, 2, 2));

    // reads informatoin on model from XD files
    xd_io::read(xdMasterFileName, xdParameterFileName, hcModelParameters, crystal, localCoordinateSystems, true);

    // set local coordinate system matrices
    size_t nAtoms = localCoordinateSystems.size();
    localCoordinateSystemMatrices.resize(nAtoms);
    for(size_t atomIdx=0; atomIdx<nAtoms; atomIdx++)
        localCoordinateSystems[atomIdx].calculate(localCoordinateSystemMatrices[atomIdx],crystal);


    // 
    HansenCoppensStructureFactorCalculator calculator(crystal, hcModelParameters);
    

    // runs the calculations
    calculator.calculateStructureFactors(crystal.atoms, localCoordinateSystemMatrices, hkl_indices, scattering_factors);

    // prints out the results (scattering factors)
    for (size_t i = 0, n = hkl_indices.size(); i<n; i++)
        cout << "h , k , l : " << hkl_indices[i][0] << " " << hkl_indices[i][1] << " " << hkl_indices[i][2]
             << "  scattering factor : " << scattering_factors[i] << endl;

}

int main(int argc,char *argv[])
{
    
    try
    {
        if (argc != 3)
        {
            cout << "expected 2 arguments: name of xd master file and xd parameter file" << endl;
            exit(0);
        }

        run(argv[1],argv[2]);
    }
    catch (const Exception & e)
    {
        cout<< e.what();
    }
    catch (const std::exception& e)
    {
        cout<< "standad exception " << e.what() << endl;
    }
    
    return 0;
}
