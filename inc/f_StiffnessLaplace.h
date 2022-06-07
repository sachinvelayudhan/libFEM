#include<vector>

void stiffnessLaplace(const std::vector<std::vector<double>>& elNodes, const std::vector<int>& elCon, C_Matrix_Sparse& kGlob, C_QuadrilateralBasis& feL, const C_GaussData2DQuad& gpData, C_Matrix_Dense& A)
{
    C_Matrix_Dense ke(4,4);
    double jac;
    double JxW;
    for(int itGp = 0; itGp < gpData.num_GP; itGp++)
    {
        C_Matrix_Dense& dsp = feL.getDerShapeFunction(itGp, elNodes, jac);
        JxW = jac*gpData.wt[itGp];
        ke += JxW*(dsp.T()*A*dsp);
    }
    //std::cout<< ke;

    kGlob.add_matr(ke, elCon, elCon); 
}