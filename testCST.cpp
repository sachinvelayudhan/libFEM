#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <list>

#include "C_Matrix_Dense.h"
#include "C_FEM_BasisFunction_2D_TRIA.h"
#include "C_FEM_GaussPoint_2D_TRIA.h"

int main()
{
    std::vector<double> x;
    std::vector<double> y;
    // C_Matrix_Dense sp;
    // sp.reshape(2,2);
    // x = {0, 1, 0};
    // y = {0, 0, 1};

    C_GaussPoint_2D_TRIA gpData(3);
    C_TriangleBasis feLag(1,gpData);

    return 0;
}