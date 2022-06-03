#ifndef C_FEM_BASISFUNCTION_2D_QUAD_H
#define C_FEM_BASISFUNCTION_2D_QUAD_H

#include <vector>
#include <math.h>
#include <iostream>
#include <tuple>

#include "C_FEM_GaussPoint_1D.h" 
#include "C_Matrix_Dense.h"

//! This class produces Lagrange basis functions at Gauss Points for a Quadrilateral element. It takes the polynomial order and a Quadrilateral Gauss point object as input.
class C_QuadrilateralBasis
{
    public:
        int num_sp;
        C_Matrix_Dense* sp;
        C_Matrix_Dense* dsp;

        //! Class Constructor
        /*!
        \param order the order of the requested interpolation functions. 
        \param GP_Data object containing Gauss Point information, such as coordinate data.
        */

    C_QuadrilateralBasis (int order, C_GaussPoint_1D GP_Data)
    {
        int itSp = 0;
        if (order == 1)
        {
            // Linear
            num_sp = 4;
            int noGp =GP_Data.num_GP*GP_Data.num_GP;

            (this -> sp) = new C_Matrix_Dense [noGp];
            (this -> dsp) = new C_Matrix_Dense [noGp];
            for (int ii = 0; ii < GP_Data.num_GP; ii++)
            {
                for (int jj = 0; jj < GP_Data.num_GP; jj++)
                {
                    double s = GP_Data.pt[ii];
                    double t = GP_Data.pt[jj];
                    (this -> sp)[itSp].reshape(1,num_sp);
                    (this -> dsp)[itSp].reshape(2,num_sp);

                    (this -> sp)[itSp](0,0) = 0.25*(1-s)*(1-t);
                    (this -> sp)[itSp](0,1) = 0.25*(1+s)*(1-t);
                    (this -> sp)[itSp](0,2) = 0.25*(1+s)*(1+t);
                    (this -> sp)[itSp](0,3) = 0.25*(1-s)*(1+t);

                    (this -> dsp)[itSp](0,0) = -0.25*(1-t);
                    (this -> dsp)[itSp](0,1) = 0.25*(1-t);
                    (this -> dsp)[itSp](0,2) = 0.25*(1+t);
                    (this -> dsp)[itSp](0,3) = -0.25*(1+t);

                    (this -> dsp)[itSp](1,0) = -0.25*(1-s);
                    (this -> dsp)[itSp](1,1) = -0.25*(1+s);
                    (this -> dsp)[itSp](1,2) = 0.25*(1+s);
                    (this -> dsp)[itSp](1,3) = 0.25*(1-s);

                    itSp = itSp + 1;
                }
            }
        }
        else
        {
           // Quadratic
            num_sp = 9; 
        }
    }
};
    
#endif