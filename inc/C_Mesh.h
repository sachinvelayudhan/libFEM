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
        std::vector< std::vector<int>>    elements;
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
    dx=(xlim[1]-xlim[0])/noX;
    dy=(ylim[1]-ylim[0])/noY;
    x=xlim[0];
    num_Nd=0;

    // Create the nodes in the rectangular mesh
    for (int i=0; i<=noX; i++)
    {
        y=ylim[0];
        for (int j=0; j<=noY; j++)
        {
            nodes.push_back({x,y});
            y+=dy;
            num_Nd++;
        }
        x+=dx;
    }

    // Create the connectivity for the quadrilateral
    num_El=noX*noY;
    int row=0;
    std::vector<int> tmp;
    for (int j=1; j<noY; j++)
    {
        for (int i=0; i<noX; i++)
        {

            tmp={row+i, row+i+1, (noY*j)+i+1,  (noY*j)+i};
            elements.push_back(tmp);
        }
        row+=noX;
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