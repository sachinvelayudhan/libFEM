#ifndef C_FEM_BASISFUNCTION_2D_QUAD_H
#define C_FEM_BASISFUNCTION_2D_QUAD_H

#include <vector>
#include <math.h>
#include <iostream>
#include <tuple>

#include "C_GaussData.h" 
#include "C_Matrix_Dense.h"

//! This class produces Lagrange basis functions at Gauss Points for a Quadrilateral element. It takes the polynomial order and a Quadrilateral Gauss point object as input.
class C_QuadrilateralBasis
{
    
    private:
        int num_sp=0;
        C_Matrix_Dense* spR;
        C_Matrix_Dense* dspR;
        C_Matrix_Dense* dsp;

        //! Class Constructor
        /*!
        \param order the order of the requested interpolation functions. 
        \param GP_Data object containing Gauss Point information, such as coordinate data.
        */
    public:
    C_QuadrilateralBasis (int order, const C_GaussData2DQuad& gpData)
    {
        int itSp = 0;
        if (order == 1)
        {
            // Linear
            num_sp = 4;
            int noGp =gpData.num_GP;

            (this -> spR) = new C_Matrix_Dense [noGp];
            (this -> dspR) = new C_Matrix_Dense [noGp];
            dsp = new C_Matrix_Dense;
            (*dsp).reshape(2,num_sp);
            for (int itGp = 0; itGp < noGp; itGp++)
            {
                double s = gpData.ptx[itGp];
                double t = gpData.pty[itGp];
                (this -> spR)[itGp].reshape(1,num_sp);
                (this -> dspR)[itGp].reshape(2,num_sp);

                (this -> spR)[itGp](0,0) = 0.25*(1-s)*(1-t);
                (this -> spR)[itGp](0,1) = 0.25*(1+s)*(1-t);
                (this -> spR)[itGp](0,2) = 0.25*(1+s)*(1+t);
                (this -> spR)[itGp](0,3) = 0.25*(1-s)*(1+t);

                (this -> dspR)[itGp](0,0) = -0.25*(1-t);
                (this -> dspR)[itGp](0,1) = 0.25*(1-t);
                (this -> dspR)[itGp](0,2) = 0.25*(1+t);
                (this -> dspR)[itGp](0,3) = -0.25*(1+t);

                (this -> dspR)[itGp](1,0) = -0.25*(1-s);
                (this -> dspR)[itGp](1,1) = -0.25*(1+s);
                (this -> dspR)[itGp](1,2) = 0.25*(1+s);
                (this -> dspR)[itGp](1,3) = 0.25*(1-s);
            }
        }
        else
        {
           // Quadratic
            num_sp = 9; 
        }
    }
    C_Matrix_Dense& getShapeFunction(int itGp)
    {   
        return (this->spR)[itGp];
    }
    C_Matrix_Dense& getDerShapeFunction(int itGp, const std::vector<std::vector<double>>& elNodes, double& j)
    {
        C_Matrix_Dense elNodesMat(4,2);
        C_Matrix_Dense jac(2,2), invJac(2,2);
        int i=0;
        for (auto el:elNodes){
            elNodesMat(i,0)=el[0];
            elNodesMat(i,1)=el[1];
            i++;
        }
        
        jac=dspR[itGp]*elNodesMat;

        j = jac(0,0)*jac(1,1) - jac(0,1)*jac(1,0);
        invJac(0,0) = jac(1,1);
        invJac(0,1) = -1*jac(0,1);
        invJac(1,0) = -1*jac(1,0);
        invJac(1,1) = jac(0,0);
        (*dsp)=(1/j)*invJac*dspR[itGp];
        return (*dsp);
    }
};
    
#endif