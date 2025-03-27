#include <stdio.h>
#include <omp.h>


// random omp code
int main(){
    int i;
    #pragma omp parallel for
    for(i = 0; i < 10; i++){
        printf("Hello from thread %d\n", i);
    }
    return 0;
}