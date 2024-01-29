#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <random>
#include <ctime>

using namespace std;


int main(int argc, char *argv[])
{
    unsigned int seed = time(0);

    unsigned long long count = 0, n = 0;
    int numThreads;
    const double PI25DT = 3.141592653589793238462643;

    cout << "Input number of tosses: ";
    cin >> n;

    cout << "Input number of threads: ";
	cin >> numThreads;

    printf("---------------\n");

    if (numThreads <= 0) {
#pragma omp parallel
		numThreads = omp_get_num_threads();
	}

    double seqPI;
    double start = omp_get_wtime();
    for(unsigned long long i = 0; i < n; i++) 
    {
        double x = (double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0;
        double y = (double)rand_r(&seed)/RAND_MAX * 2.0 - 1.0;
        double z = (x*x) + (y*y);
        if(z <= 1.0) 
            count++;
    }
    seqPI = (((double) count/ (double) n) * 4);
    double stop = omp_get_wtime();
    double seqTime = stop - start;
    printf("Sequential pi is approximately %.16f, Error is %.16f\n", seqPI ,fabs(seqPI - PI25DT));
    printf("Time (sec.) %.16f\n", seqTime);
    printf("---------------\n");

    count = 0;
    double parPI;
    start = omp_get_wtime();
#pragma omp parallel for num_threads(numThreads) private(seed) reduction(+:count)
	for(unsigned long long i = 0; i < n; i++) 
    {   
        double x = (double)rand_r(&seed)/RAND_MAX;
        double y = (double)rand_r(&seed)/RAND_MAX;
        double z = (x*x) + (y*y);
        if(z <= 1.0) 
            count++;
    }  
    parPI = (((double) count/ (double) n) * 4);
    stop = omp_get_wtime();

    double parTime = stop - start;
    printf("Parallel pi is approximately %.16f, Error is %.16f\n", parPI, fabs(parPI - PI25DT));
	printf("Time (sec.) %.16f\n", parTime);
    printf("---------------\n");
    printf("Speedup is %.16f, Efficiency is %.16f \n", (seqTime/parTime), ((seqTime/parTime)/numThreads));
    return 0;
}