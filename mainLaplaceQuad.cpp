#include "./eigen-3.4.0/Eigen/Dense" 
#include "./eigen-3.4.0/Eigen/Sparse"

#include "./inc/C_GaussData.h" 
#include "./inc/C_FEM_BasisFunction_2D_QUAD.h"
#include "./inc/C_Matrix_Sparse.h"
#include "./inc/C_Matrix_Dense.h"
#include "./inc/C_Mesh.h"

#include "./inc/f_StiffnessLaplace.h"
#include "./Solver_Interfaces/f_SolverInterfaces.h"

int main()
{
    C_GaussData2DQuad gpData(2);
    C_Mesh2D mesh;

    int noX = 4, noY = 2, IEL = 1, NNM, NDF = 1, totDof;
    mesh.meshRectangle({0,0.18},{0,0.1},noX,noY);
    C_Matrix_Sparse kGlob;
    C_QuadrilateralBasis feL(1, gpData);
    int itEl = 0;

    NNM = (noX*IEL + 1)*(noY*IEL + 1);
    totDof = NNM*NDF;

    // Material coefficient matrix
    double A11 = 0.2, A22 = 0.2;
    C_Matrix_Dense A(2,2);

    A(0,0) = A11; A(1,1) = A22;
 
    for (const auto& elCon: mesh.elements)
    {
        std::cout << (itEl+1);
        std::vector<std::vector<double>> elNodes = mesh.getElNodes(itEl);
        stiffnessLaplace(elNodes, elCon, kGlob, feL, gpData, A);
        itEl++;
    }  

    //std::cout << kGlob;
    Eigen::VectorXd fGlobal(totDof);

    // i. Neumann Boundary Conditions
    // int fDof = totDof - 1;
    // fGlobal(fDof) = 10.0;

    // ii. Dirichlet Boundary Conditions
    int bc, NSPV = 6;
    std::vector<int> ISPV;
    std::vector<double> VSPV;

    ISPV = {0,4,5,9,10,14};
    VSPV = {500,300,500,300,500,300};

    for(int ii = 0; ii < NSPV; ii++)
    {
        bc = ISPV[ii];
        kGlob.col_NonSparseAssign(0.0, bc);
        kGlob.row_NonSparseAssign(0.0, bc);
        kGlob(bc,bc) = 1.0;
        fGlobal(bc) = VSPV[ii]; 
    }

    // Solve, using Cholesky Factorization of kGlob
    Eigen::SparseMatrix<double> kG_eigen(totDof, totDof);
    convert_to_Eigen(kGlob, kG_eigen);

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > chol;
    chol.compute(kG_eigen);  
    Eigen::VectorXd sol = chol.solve(fGlobal);
    std::cout << "\n";
    std::cout << sol;
    std::cout << "\n";

    return 0;

}