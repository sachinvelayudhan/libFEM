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
    C_GaussPoint_2D_TRIA gpData(1);
    C_TriangleBasis feLag(1,gpData);
    C_Matrix_Dense ke(3,3);
    for (int itGp=0; itGp<gpData.num_GP; itGp++)
    {
        ke=(ke+1.0*(feLag.dsp[itGp].T()*feLag.dsp[itGp]));
    }
    std::cout<< ke;
    return 0;
}