#include<vector>
#include "C_Matrix_Dense.h"
#include "C_Matrix_Sparse.h"

void StiffnessGraFEA2D (const std::vector<std::vector<double>>& elNodes, const std::vector<int>& elCon, C_Matrix_Sparse& kGlob, int NPLANE, double H, int NDF)
{
    double E1 = 1, E2 = 1, nu12 = 0.3, G12 = 1;
    double nu21, temp, temp2, temp1;
    double A;

    std::vector<double> beta;
    std::vector<double> gama;

    std::vector<int> con;

    C_Matrix_Dense C(3,3);
    C_Matrix_Dense B(3,6), ke(6,6);

    nu21 = nu12 * (E2/E1);
    temp = 1 - nu12*nu21;
    temp1 = (1 - nu12 -2*nu12*nu21);
    temp2 = (1 + nu12) * temp1;

    if (NPLANE == 1)
    {
        C(0,0) = E1 / temp;
        C(1,1) = C(0,0) * (E2/E1);
        C(0,1) = nu12 * C(1,1);
        C(1,0) = C(0,1);
        C(2,2) = G12;
    }
    else
    {
       C(0,0) = E1 * (temp/temp2);
       C(1,1) = C(0,0) * (E2/E1);
       C(0,1) = nu12 * E2 / temp1;
       C(1,0) = C(0,1);
       C(2,2) = G12; 
    }

    beta = {(elNodes[1][1] - elNodes[2][1]), (elNodes[2][1] - elNodes[0][1]), (elNodes[0][1] - elNodes[1][1])};
    gama = {-(elNodes[1][0] - elNodes[2][0]), -(elNodes[2][0] - elNodes[0][0]), -(elNodes[0][0] - elNodes[1][0])};

    A = 0.5 * (((elNodes[2][1] - elNodes[0][1])*(elNodes[1][0] - elNodes[0][0])) - ((elNodes[1][1] - elNodes[0][1])*(elNodes[2][0] - elNodes[0][0])));

    B(0,0) = (1/(2*A))*beta[0];
    B(0,2) = (1/(2*A))*beta[1];
    B(0,4) = (1/(2*A))*beta[2];

    B(1,1) = (1/(2*A))*gama[0];
    B(1,3) = (1/(2*A))*gama[1];
    B(1,5) = (1/(2*A))*gama[2];

    B(2,0) = (1/(2*A))*gama[0];
    B(2,1) = (1/(2*A))*beta[0];
    B(2,2) = (1/(2*A))*gama[1];
    B(2,3) = (1/(2*A))*beta[1];
    B(2,4) = (1/(2*A))*gama[2];
    B(2,5) = (1/(2*A))*beta[2];

    ke = A*H*(B.T()*C*B);

    int nr;
    for (int i = 0; i < 3; i++)
    {
        nr = elCon[i]*NDF - NDF;
        
        for(int j = 0; j < NDF; j++)
        {
            con.push_back(nr+j);
        }
    }

    kGlob.add_matr(ke, con, con); 
}