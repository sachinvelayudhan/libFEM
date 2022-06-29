#include "./inc/C_Matrix_Sparse.h"
#include "./inc/C_Matrix_Dense.h"

#include "./inc/C_Loadtxt.h"
#include "./inc/f_StiffnessGraFEA2D.h"

/*
NEM - Number of elements
NPE - Number of nodes per element
NDF - Number of degrees of freedom per node
ELXY - Element coordinate matrix
CONMAT - Element node numbers
NPLANE - (=1: Plane stress condition), (=2: Plane strain condition)
*/

int main()
{
    int NEM = 1, NPE = 3, NI, NDF = 2, NPLANE = 1;
    std::vector < std::vector<double> > elNodes; 
    std::vector < int > elCon;  
    C_Matrix_Sparse kGlob (6,6);
    double H = 1.0;

    auto nodes = loadtxt <double> ( "nodes.txt" );
    auto elements = loadtxt <int> ( "elements.txt" );
    
    for (int i = 0; i < NEM; i++)
    {
        for (int j = 0; j < NPE; j++)
        {
            NI = elements[i][j] - 1;
            elNodes.push_back(nodes[NI]);
            elCon.push_back(elements[i][NI]);
        }

        StiffnessGraFEA2D(elNodes,elCon,kGlob,NPLANE,H,NDF);
        cout << kGlob;
    }
    
    return 0;
}