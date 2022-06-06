#ifndef C_FEM_GAUSSPOINT_2D_TRIA_H
#define C_FEM_GAUSSPOINT_2D_TRIA_H

#include <vector>
#include <math.h>
#include <iostream>
#include "C_FEM_GaussPoint_1D.h"

class C_GaussPoint_2D_QUAD
{
    public:
        int num_GP; //!< Number of Gauss Points
        std::vector<double> ptx; //!< Vector of Gauss Points (x-coordinates)
        std::vector<double> pty; //!< Vector of Gauss Points (y-coordinates)
        std::vector<double> wt; //!< Corresponding Gauss Weights
    
    C_GaussPoint_2D_QUAD (int n_in)
    {
        C_GaussPoint_1D gaussDat(n_in);
        num_GP = n_in*n_in;
        for (int i=0; i<n_in; i++)
        {
            double xi=gaussDat.pt[i];
            double w1=gaussDat.wt[i];
            for (int j=0; j<n_in; j++)
            {
                double eta=gaussDat.pt[j];
                double w2=gaussDat.wt[j];
                ptx.push_back(xi);
                pty.push_back(eta);
                wt.push_back(w1*w2);
            }
        }
    }
};

#endif