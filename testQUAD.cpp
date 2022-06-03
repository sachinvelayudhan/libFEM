#include <eigen-3.4.0/Eigen/Dense>  // Linear Solver Libraries
#include <eigen-3.4.0/Eigen/Sparse>

#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <list>

#include "C_Matrix_Dense.h"
#include "C_Matrix_Sparse.h"
#include "C_Mesh_2D_QUAD.h"
#include "C_FEM_BasisFunction_2D_QUAD.h"
#include "C_FEM_GaussPoint_1D.h" 
#include "C_Stiffness_Laplace_2D_QUAD.h"
#include "f_SolverInterfaces.h"

int main()
{
    int NDF = 1, NPE = 4, NX = 1, NY = 1, IEL = 1; 
    double X0 = 0, Y0 = 0, a = 1, b = 1;
    int NNM, NEL, totDof, NEM;

    NNM = (NX*IEL + 1)*(NY*IEL + 1);
    NEL = NPE*NDF;
    totDof = NNM*NDF;
    NEM = NX*NY;

    double A11 = 1.0, A22 = 1.0; // Material coefficients
    C_Matrix_Sparse kGlobal;
    Eigen::VectorXd fGlobal(totDof);
    
    /* X0 - x-coordinate of first node
    Y0 - y-coordinate of first node
    NX - Number of division in x-direction
    NY - Number of division in y-direction
    IEL - =1(First order), =2(Second order)
    a - Length in x - direction
    b - Length in y - direction
    NDF - Number of degrees of freedom
    NPE - Number of nodes per element
    NEM - Number of elements in the mesh
    NNM - Total number of nodes in the mesh
    NEL - Size of element stiffness matrix
    totDof - Size of global stiffness matrix
    GLXY - Global coordinate matrix
    ELXY - Element coordinate matrix
    conmat - Connectivity matrix
    A - Material coefficient matrix
    */

    C_Matrix_Dense A(2,2);
    C_Mesh_2D_QUAD mesh2D(NX,NY,IEL,X0,Y0,a,b,NPE,NEM);
    C_GaussPoint_1D gpData(2);
    C_QuadrilateralBasis feLag(1,gpData);
    double jacob, cnst;

    // Material properties
    A(0,0) = A11; A(1,1) = A22;

    // This function returns the global stiffness matrix
    for (int ielem = 0; ielem < NEM; ielem++)
    {
        Stiffness_Laplace_2D_QUAD (NEM, NPE, NDF, NEL, ielem, kGlobal, gpData, feLag, A, mesh2D);
    }

    // i. Neumann Boundary Conditions
    int fDof = totDof - 1;
    fGlobal(fDof) = 10.0;

    // ii. Dirichlet Boundary Conditions
    int bcDof = 0;
    kGlobal.col_NonSparseAssign(0.0, bcDof);
    kGlobal.row_NonSparseAssign(0.0, bcDof);
    kGlobal(bcDof,bcDof) = 1;
    fGlobal(bcDof) = 0.0; 

    // Solve, using Cholesky Factorization of kGlob
    Eigen::SparseMatrix<double> kG_eigen(totDof, totDof);
    convert_to_Eigen(kGlobal, kG_eigen);

    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > chol;
    chol.compute(kG_eigen);  
    Eigen::VectorXd sol = chol.solve(fGlobal);
    std::cout << sol;
    std::cout << "\n";

    return 0;
}