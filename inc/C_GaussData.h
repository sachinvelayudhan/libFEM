#ifndef C_GAUSSDATA_H
#define C_GAUSSDATA_H

#include <vector>
#include <math.h>
#include <iostream>

class C_GaussData{
    public:
        int dim;
        int num_GP; //!< Number of Gauss Points
        std::vector<double> ptx; //!< Vector of Gauss Points
        std::vector<double> pty;
        std::vector<double> ptz;
        std::vector<double> wt; //!< Corresponding Gauss Weights
};

class C_GaussData1D: public C_GaussData{
    public:
    C_GaussData1D (int n_in){
        dim=1;
        num_GP = n_in;

        if (num_GP == 1){
            ptx = {0.0};
            wt = {2.0};
        }
        else if (num_GP == 2){
            ptx = {0.577350269189626, -0.577350269189626};
            wt = {1.0,                1.0};
        }
        else if (num_GP == 3){
            ptx = {-0.774596669241483, 0,   0.774596669241483};
            wt = { 0.55555555555,               0.88888888888, 0.55555555555};
        }
        else if (num_GP == 4){
            ptx = {-0.86113631, -0.33998104, 0.33998104, 0.86113631};
            wt = { 0.34785485,  0.65214515, 0.65214515, 0.34785485};
        }
        else {
            ptx = {-0.906180, -0.538469, 0.0,      0.538469, 0.906180};
            wt = { 0.236927,  0.478629, 0.568889, 0.478629, 0.236927};
            num_GP = 5;
        }
    }
};

class C_GaussData2DQuad: public C_GaussData{
    public:
        C_GaussData2DQuad(int order)
        {
            C_GaussData1D gaussDat(order);
            num_GP = order*order;
            for (int i=0; i<order; i++)
            {
                double xi=gaussDat.ptx[i];
                double w1=gaussDat.wt[i];
                for (int j=0; j<order; j++)
                {
                    double eta=gaussDat.ptx[j];
                    double w2=gaussDat.wt[j];
                    ptx.push_back(xi);
                    pty.push_back(eta);
                    wt.push_back(w1*w2);
                }
            }
        }
};
#endif