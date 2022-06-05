#ifndef C_MESH_2D_QUAD_H
#define C_MESH_2D_QUAD_H

#include <math.h>
#include <vector>

#include "f_MiscellaneousFunctions.h"
#include "C_Matrix_Dense.h"

//! This class is used to store data describing a 2D quadrilateral mesh.
class C_Mesh_2D_QUAD
{
    public:
    C_Matrix_Dense GLXY, nodemat, conmat;

    C_Mesh_2D_QUAD (int NX, int NY, int IEL, double X0, double Y0, double a, double b, int NPE, int NEM)
    {
        double ndx, ndy, x, y;
        int start, starts, ronode, culnode;
        int nodex = NX*IEL + 1, nodey = NY*IEL + 1, num_node = nodex*nodey;

        this -> GLXY = C_Matrix_Dense(num_node,2);
        this -> nodemat = C_Matrix_Dense(nodey, nodex);
        this -> conmat = C_Matrix_Dense(NEM, NPE);

        ndx = a / (nodex - 1);
        ndy = b / (nodey - 1);

        y = Y0;
        start = 0;

        for (int i = 0; i < nodey; i++)
        {
            x = X0;

            for (int j = 0; j < nodex; j++)
            {
                GLXY(start, 0) = x;
                GLXY(start, 1) = y;
                start = start + 1;
                nodemat(i,j) = start;
                x = x + ndx;
            }

            y = y + ndy;

        }

        starts = 0;

        for (int k = 0; k < NY; k++)
        {
            ronode = k;
            for (int l = 0; l < NX; l++)
            {
                culnode = l;
                conmat(starts,0) = nodemat(ronode,culnode);
                conmat(starts,1) = nodemat(ronode,(culnode+1));
                conmat(starts,2) = nodemat((ronode+1),(culnode+1));
                conmat(starts,3) = nodemat((ronode+1),culnode);
                starts = starts + 1;
            }
        }
    }       
};

#endif