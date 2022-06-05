#ifndef C_STIFFNESS_LAPLACE_2D_QUAD_H
#define C_STIFFNESS_LAPLACE_2D_QUAD_H

#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <list>

#include "C_Matrix_Dense.h"
#include "C_Matrix_Sparse.h"
#include "C_FEM_BasisFunction_2D_QUAD.h"
#include "C_FEM_GaussPoint_1D.h"
#include "C_Mesh_2D_QUAD.h"

void Stiffness_Laplace_2D_QUAD (int NEM, int NPE, int NDF, int NEL, int ielem, C_Matrix_Sparse& kGlobal, C_GaussPoint_1D& gpData, C_QuadrilateralBasis& feLag, C_Matrix_Dense& A, C_Mesh_2D_QUAD& mesh2D)
{
    C_Matrix_Dense J(2,2), Jinv(2,2), B(2,NEL), Ke(NEL,NEL), ELXY(NPE,2);
    double jacob, cnst;
    std::vector<int> NR;
    std::vector<int> NC;
    int start = 0;

    for (int j = 0; j < NPE; j++)
    {
        ELXY(j, 0) = mesh2D.GLXY(mesh2D.conmat(ielem,j)-1, 0);
        ELXY(j, 1) = mesh2D.GLXY(mesh2D.conmat(ielem,j)-1, 1);

        NR.push_back(mesh2D.conmat(ielem,j)-1);
        NC.push_back(mesh2D.conmat(ielem,j)-1);
    }

    for (int i = 0; i < gpData.num_GP; i++)
    {
        for (int j = 0; j < gpData.num_GP; j++)
        {
            J = feLag.dsp[start]*ELXY;
            jacob = J(0,0)*J(1,1) - J(0,1)*J(1,0);
            Jinv(0,0) = J(1,1);
            Jinv(0,1) = -1*J(0,1);
            Jinv(1,0) = -1*J(1,0);
            Jinv(1,1) = J(0,0);

            B = (1/jacob)*Jinv*feLag.dsp[start];
            cnst = jacob*gpData.wt[i]*gpData.wt[j];

            Ke =  (Ke + cnst*(B.T()*A*B));
            start = start + 1;
        }
    }

    kGlobal.add_matr(Ke, NR, NC); 
};

#endif