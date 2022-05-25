#include <math.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <list>

#include "C_Matrix_Sparse.h"
#include "C_Matrix_Dense.h"

int main()
{
    // Test 1: Dense Matrix Operations
    C_Matrix_Dense A1(2,2), B1(4,4), C1(4,6);
    C_Matrix_Dense A2(4,3), B2(3,4), C2(3,3);

    C_Matrix_Dense A3(2,2);
    A3 = {1, 2, 3, 4};
    A3+=A3;
    A1=A3;
    std::cout << A3;
    std::cout << A1-A1;

    return 0;
}