#include "discamb/MathUtilities/RealSphericalHarmonics.h"
#include "discamb/BasicUtilities/OnError.h"
#include "discamb/MathUtilities/MathUtilities.h"

#include <fstream>
#include <map>

using namespace std;

namespace {
    vector<vector<double> > initializeConversionFactors()
    {
        vector<vector<double> > m;

        m.resize(8);
        for (int l = 0; l <= 7; l++)
            m[l].resize(l + 1);

        m[0][0] = 0.28209479177388;
        m[1][0] = 0.65147001587056; m[1][1] = 0.65147001587056;
        m[2][0] = 0.65552905835525; m[2][1] = 0.68646842464783; m[2][2] = 0.68646842464783;
        m[3][0] = 0.65613421114746; m[3][1] = 0.70087743932709; m[3][2] = 0.69189513695864; m[3][3] = 0.71929123343438;
        m[4][0] = 0.65620991154888; m[4][1] = 0.70847465461627; m[4][2] = 0.69879556866550; m[4][3] = 0.70616251710906; m[4][4] = 0.74899845670137;
        m[5][0] = 0.65617176926179; m[5][1] = 0.71306675774595; m[5][2] = 0.70407303654785; m[5][3] = 0.70547522628220; m[5][4] = 0.72266090165344; m[5][5] = 0.7759136537823;
        m[6][0] = 0.65611016832343; m[6][1] = 0.71609595329005; m[6][2] = 0.70800876706179; m[6][3] = 0.70703232886769; m[6][4] = 0.71554818711837; m[6][5] = 0.73945146220726; m[6][6] = 0.8004796889213;
        m[7][0] = 0.65605000193133; m[7][1] = 0.71822036312964; m[7][2] = 0.71099470247414; m[7][3] = 0.70896558772175; m[7][4] = 0.71344649338342; m[7][5] = 0.72696198663051; m[7][6] = 0.75586913625466; m[7][7] = 0.82308112861385;

        return m;
    }

    map<string, vector<pair<int,int> > > symmetry_invariants = 
    {
        {"1",{{0,0},{1,-1},{1,0},{1,1},{2,-2},{2,-1},{2,0},{2,1},{2,2},{3,-3},{3,-2},{3,-1},{3,0},{3,1},{3,2},{3,3},{4,-4},{4,-3},{4,-2},{4,-1},{4,0},{4,1},{4,2},{4,3},{4,4}}},
        {"no",{{0,0},{1,-1},{1,0},{1,1},{2,-2},{2,-1},{2,0},{2,1},{2,2},{3,-3},{3,-2},{3,-1},{3,0},{3,1},{3,2},{3,3},{4,-4},{4,-3},{4,-2},{4,-1},{4,0},{4,1},{4,2},{4,3},{4,4}}},
        {"cyl",{{0,0},{1,0},{2,0},{3,0},{4,0}}},
        {"sph",{{0,0}}},
        {"-1",{{0,0},{2,-2},{2,-1},{2,0},{2,1},{2,2},{4,-4},{4,-3},{4,-2},{4,-1},{4,0},{4,1},{4,2},{4,3},{4,4}}},
        {"2",{{0,0},{1,0},{2,-2},{2,0},{2,2},{3,-2},{3,0},{3,2},{4,-4},{4,-2},{4,0},{4,2},{4,4}}},
        {"m",{{0,0},{1,-1},{1,1},{2,-2},{2,0},{2,2},{3,-3},{3,-1},{3,1},{3,3},{4,-4},{4,-2},{4,0},{4,2},{4,4}}},
        {"2/m",{{0,0},{2,-2},{2,0},{2,2},{4,-4},{4,-2},{4,0},{4,2},{4,4}}},
        {"222",{{0,0},{2,0},{2,2},{3,-2},{4,0},{4,2},{4,4}}},
        {"mm2",{{0,0},{1,0},{2,0},{2,2},{3,0},{3,2},{4,0},{4,2},{4,4}}},
        {"mmm",{{0,0},{2,0},{2,2},{4,0},{4,2},{4,4}}},
        {"4",{{0,0},{1,0},{2,0},{3,0},{4,-4},{4,0},{4,4}}},
        {"-4",{{0,0},{2,0},{3,-2},{3,2},{4,-4},{4,0},{4,4}}},
        {"4/m",{{0,0},{2,0},{4,-4},{4,0},{4,4}}},
        {"422",{{0,0},{2,0},{4,0},{4,4}}},
        {"4mm",{{0,0},{1,0},{2,0},{3,0},{4,0},{4,4}}},
        {"-42m",{{0,0},{2,0},{3,-2},{4,0},{4,4}}},
        {"4/mmm",{{0,0},{2,0},{4,0},{4,4}}},
        {"3",{{0,0},{1,0},{2,0},{3,-3},{3,0},{3,3},{4,-3},{4,0},{4,3}}},
        {"-3",{{0,0},{2,0},{4,-3},{4,0},{4,3}}},
        {"32",{{0,0},{2,0},{3,-3},{4,0},{4,3}}},
        {"3m",{{0,0},{1,0},{2,0},{3,0},{3,3},{4,0},{4,3}}},
        {"-3m",{{0,0},{2,0},{4,0},{4,3}}},
        {"6",{{0,0},{1,0},{2,0},{3,0},{4,0}}},
        {"-6",{{0,0},{2,0},{3,-3},{3,3},{4,0}}},
        {"6/m",{{0,0},{2,0},{4,0}}},
        {"622",{{0,0},{2,0},{4,0}}},
        {"6mm",{{0,0},{1,0},{2,0},{3,0},{4,0}}},
        {"-6m2",{{0,0},{2,0},{3,3},{4,0}}},
        {"6/mmm",{{0,0},{2,0},{4,0}}}
    };
}

namespace discamb {

namespace real_spherical_harmonics{


    std::vector<std::vector<double> > getDensityToWfnNormalization()
    {
        return initializeConversionFactors();
    }

    bool symmetryInvariant(
        const std::string &pointGroupLabel,
        std::vector<std::pair<int, int> > &plms)
    {
        
        auto it = symmetry_invariants.find(pointGroupLabel); 
        if (it == symmetry_invariants.end())
            return false;
        plms = it->second;
        return true;
    }

 void getDensityToWfnMultipliers(int maxL,std::vector<std::vector<double> > &values)
 {
     static const vector<vector<double> > densityToWfnNormalization = initializeConversionFactors();
     if(maxL>7)
         on_error::throwException("density normalized spherical harmonics with L>7 are not supported",__FILE__,__LINE__); 
     values.clear();
     values.insert(values.end(), densityToWfnNormalization.begin(), densityToWfnNormalization.begin()+maxL+1);

 }

 double densityNormalized(
     const Vector3d &normalizedVector3D, 
     int l, 
     int m)
 {
    return densityNormalizationMultipliers[l][l+m]*polynomial(normalizedVector3D,l,m);
 }

 double wfnNormalized(
     const Vector3d &normalizedVector3D, 
     int l, 
     int m)
 {
     return wfnNormalizationMultipliers[l][l + m] * polynomial(normalizedVector3D, l, m);
 }



 double polynomial(
     const Vector3d &v, // noralized 3D vector, 
     int l, 
     int m)
 {
     
     switch(l)
     {
     case 0:
         return 1.0;
     case 1:
         if(m==-1)
             return v.y;
         else
         {
             if(m==0)
                 return v.z;
             else
                 return v.x;
         }
         return 0;
     case 2:
         switch(m)
         {
         case -2:
             return v.x*v.y;
         case -1:
             return v.y*v.z;
         case 0:
             return 3*v.z*v.z-1;
         case 1:
             return v.x*v.z;
         case 2:
             return v.x*v.x-v.y*v.y;
         default:
             return 0;
         }
     case 3:
         switch(m)
         {
         case -3:
             return (3*v.x*v.x - v.y*v.y)*v.y;
         case -2:
             return v.x*v.y*v.z;
         case -1:
             return v.y * (5*v.z*v.z - 1);
         case 0:
             return v.z*(5*v.z*v.z - 3);
         case 1:
             return v.x * (5 *v.z*v.z - 1);
         case 2:
             return (v.x*v.x - v.y*v.y)*v.z;
         case 3:
             return (v.x*v.x - 3*v.y*v.y) * v.x;
         default:
             return 0;
         }
     case 4:
         switch (m)
         {
         case -4:
             return v.x*v.y * (v.x*v.x - v.y*v.y);
         case -3:
             return (3*v.x*v.x - v.y*v.y) * v.y*v.z;
         case -2:
             return v.x*v.y * (7*v.z*v.z - 1);
         case -1:
             return v.y*v.z * (7*v.z*v.z - 3);
         case 0:
             return v.z*v.z*(35*v.z*v.z - 30) + 3;
         case 1:
             return v.x*v.z * (7*v.z*v.z - 3);
         case 2:
             return (v.x*v.x - v.y*v.y)*(7*v.z*v.z - 1);
         case 3:
             return (v.x*v.x - 3*v.y*v.y)*v.x*v.z;
         case 4:
             return v.x*v.x * (v.x*v.x - 3* v.y*v.y)- v.y*v.y * (3*v.x*v.x - v.y*v.y);
         default:
             return 0;
         }
     default:
         return 0;
     } // maxL>0
     return 0.0;
 }


}

} // namespace discamb


