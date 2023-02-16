#include <omp.h>
#include <iostream>
int main()
{
    int A[1] = {-1};
    #pragma omp target
    {
        A[0] = omp_is_initial_device();
    }
    if(!A[0]) std::cout << "Able to use offloading\n";
}