
#include <cstdio>
#include <omp.h>

int main() {
    int A[1] = {-1};
    #pragma omp target
    {
    A[0] = omp_is_initial_device();
    }

    if (!A[0]) {
    printf("Able to use offloading!\n");
    }
}