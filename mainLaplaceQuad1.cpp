#include "./eigen-3.4.0/Eigen/Dense" 
#include "./eigen-3.4.0/Eigen/Sparse"

#include "./inc/C_GaussData.h" 
#include "./inc/C_FEM_BasisFunction_2D_QUAD.h"
#include "./inc/C_Matrix_Sparse.h"
#include "./inc/C_Matrix_Dense.h"
#include "./inc/C_Mesh.h"

#include "./inc/f_StiffnessLaplace.h"
#include "./Solver_Interfaces/f_SolverInterfaces.h"

#include<fstream>

int main()
{
    C_GaussData2DQuad gpData(2);
    C_Mesh2D mesh;

    std::ofstream out;

    int noX = 2, noY = 2, IEL = 1, NNM, NDF = 1, totDof;
    double x0 = 0, xa = 1.0, y0 = 0, yb = 1.0;
    mesh.meshRectangle({x0,xa},{y0,yb},noX,noY);
    C_Matrix_Sparse kGlob(mesh.num_Nd,mesh.num_Nd), fGlob(mesh.num_Nd,1);
    C_QuadrilateralBasis feL(1, gpData);
    int itEl = 0;

    NNM = (noX*IEL + 1)*(noY*IEL + 1);
    totDof = NNM*NDF;

    // Material coefficient matrix
    double A11 = 1.0, A22 = 1.0;
    C_Matrix_Dense A(2,2);

    A(0,0) = A11; A(1,1) = A22;
 
    for (const auto& elCon: mesh.elements)
    {
        // std::cout << (itEl+1);
        std::vector<std::vector<double>> elNodes = mesh.getElNodes(itEl);
        stiffnessLaplace(elNodes, elCon, kGlob, fGlob, feL, gpData, A);
        itEl++;
    }  

    Eigen::VectorXd fGlobal(totDof);

     // ii. Dirichlet Boundary Conditions
    int bc, NSPV = 5;
    std::vector<int> ISPV;
    std::vector<double> VSPV;

    // 2 x 2 Mesh
    ISPV = {2,5,6,7,8};
    VSPV = {0,0,0,0,0};

    // 4 x 4 Mesh
    // ISPV = {4,9,14,19,20,21,22,23,24};
    // VSPV = {0,0,0,0,0,0,0,0,0};

    for(int ii = 0; ii < NSPV; ii++)
    {
        bc = ISPV[ii];
        kGlob.col_NonSparseAssign(0.0, bc);
        kGlob.row_NonSparseAssign(0.0, bc);
        kGlob(bc,bc) = 1.0;
        fGlob(bc,0) = VSPV[ii]; 
    }

    // Converting force vector to eigen::force vector
    for(int i = 0; i < totDof; i++)
    {
        fGlobal(i) = fGlob(i,0);
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

    out.open("solution.txt");

    //Writing solution vector to a text file
    out << sol;
    out.close();

    return 0;

}