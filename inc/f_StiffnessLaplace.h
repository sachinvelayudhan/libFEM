#include<vector>

void stiffnessLaplace(const std::vector<std::vector<double>>& elNodes, const std::vector<int>& elCon, C_Matrix_Sparse& kGlob, C_Matrix_Sparse& fGlob, C_QuadrilateralBasis& feL, const C_GaussData2DQuad& gpData, C_Matrix_Dense& A)
{
    C_Matrix_Dense ke(4,4), fe(4,1);
    const std::vector<int> one_vec = {0}; // Vector of value 1
    double jac;
    double JxW;

    for(int itGp = 0; itGp < gpData.num_GP; itGp++)
    {
        C_Matrix_Dense& dsp = feL.getDerShapeFunction(itGp, elNodes, jac);
        C_Matrix_Dense& spR = feL.getShapeFunction(itGp);
        JxW = jac*gpData.wt[itGp];
        ke += JxW*(dsp.T()*A*dsp);
        fe += JxW*spR.T();
    }

    // std::cout<< ke;
    // std::cout<< fe;

    kGlob.add_matr(ke, elCon, elCon); 
    fGlob.add_matr(fe, elCon, one_vec); 
}