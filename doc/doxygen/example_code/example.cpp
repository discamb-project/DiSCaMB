#include "discamb/CrystalStructure/UnitCell.h"
#include "discamb/CrystalStructure/SpaceGroup.h"
#include <iostream>

using namespace discamb;
using namespace std;

int main()
{
    Vector3d r_f, r_cart;
    // defines unit cell
    UnitCell uc(8.006, 5.482, 19.444, 90.00, 110.46, 90.00);
    r_f.set(0.2432, 0.3076, -0.04171);

    // transforms position from fractional to cartesian and print out

    uc.fractionalToCartesian(r_f, r_cart);

    cout<< r_cart[0] << " " << r_cart[1] << " " << r_cart[2] << endl;

    // apply symmetry operation

    SpaceGroupOperation operation("-x,y+1/2,-z+1/2");
    Vector3d transformed_r_f;

    operation.apply(r_f, transformed_r_f);

    // transform to cartesian and print out

    uc.fractionalToCartesian(transformed_r_f, r_cart);

    cout<< r_cart[0] << " " << r_cart[1] << " " << r_cart[2] << endl;
}

