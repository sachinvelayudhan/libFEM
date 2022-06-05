#include <eigen-3.4.0/Eigen/Dense> 
#include <eigen-3.4.0/Eigen/Sparse>

#include "C_GaussData.h" 
#include "C_FEM_BasisFunction_2D_QUAD.h"
#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"
#include "C_Mesh.h"

#include "f_StiffnessLaplace.h"
#include "f_SolverInterfaces.h"


int main()
{
    C_GaussData2DQuad gpData(2);
    C_Mesh2D mesh;
    mesh.meshRectangle({0,1},{0,1},2,2);
    C_Matrix_Sparse kGlob;
    C_QuadrilateralBasis feL(1, gpData);
    int itEl=0;
 
    for (const auto& elCon: mesh.elements)
    {
        std::cout<< itEl;
        std::vector<std::vector<double>> elNodes=mesh.getElNodes(itEl);
        stiffnessLaplace(elNodes, elCon, kGlob, feL, gpData);
        itEl++;
    }    

    return 0;
}