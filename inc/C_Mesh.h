#ifndef C_MESH_H
#define C_MESH_H

#include <math.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>

#include "f_MiscellaneousFunctions.h"

/*Base mesh class for all 1D, 2D and 3D mesh*/
class C_Mesh{
    public:
        std::vector< std::vector<double>> nodes;
        std::vector< std::vector<int>> elements;
        int num_NPE; //! Nodes per element
        int num_Nd; //! Total number of nodes
        int num_El; //! Total number of elements
        int dim;     //! Dimension of mesh

    C_Mesh() {}
};

class C_Mesh2D: public C_Mesh
{
    public:
        void meshRectangle(std::vector<double>,std::vector<double>, int, int);
        std::vector<std::vector<double>> getElNodes(int);
};

void C_Mesh2D::meshRectangle(std::vector<double> xlim, std::vector<double> ylim, int noX, int noY)
{
    double x,y;
    double dx, dy;
    int nodex, nodey, IEL = 1; // IEL = 1 for linear element, IEL = 2 for quadratic elements
    
    nodex = noX*IEL + 1; // Number of nodes in x-direction
    nodey = noY*IEL + 1; // Number of nodes in y-direction
    dx = (xlim[1]-xlim[0])/(nodex-1);
    dy = (ylim[1]-ylim[0])/(nodey-1);
    y = ylim[0];
    num_Nd = nodex*nodey; // Total number of nodes

    // Create the nodes in the rectangular mesh
    for (int i = 0; i < nodey; i++)
    {
        x = xlim[0];
        for (int j = 0; j < nodex; j++)
        {
            nodes.push_back({x,y});
            x += dx;
        }
        y += dy;
    }

    // Create the connectivity for the quadrilateral
    num_El = noX*noY;
    int row = 0;
    std::vector<int> tmp;
    for (int j = 1; j <= noY; j++)
    {
        for (int i = 0; i < noX; i++)
        {
            tmp = {row+i, row+i+1, noX+2+i+row, noX+1+i+row};
            elements.push_back(tmp);
        }
        row = row + noX + 1;
    }   
}

std::vector<std::vector<double>> C_Mesh2D::getElNodes(int itEl)
{
    std::vector<std::vector<double>> elNodes;
    for (auto node:elements[itEl])
    {
        elNodes.push_back(nodes[node]);
    }
    return elNodes;
}

#endif