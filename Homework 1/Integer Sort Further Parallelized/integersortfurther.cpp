#include <iostream>
#include <sstream>
#include <omp.h>
#include <time.h>
#include <algorithm>

using namespace std;

void intSort(int *A, int sz, int max)
{
    int *counts = new int[max + 1]();
    int *temp = new int[sz];

    for(int i = 0; i < sz; i++)
        counts[A[i]]++;

    for(int i = 0; i < max + 1; i++)
        counts[i] += counts[i - 1];

    for(int i = 0; i < sz; i++)
        temp[--counts[A[i]]] = A[i];

    for(int i = 0; i < sz; i++)
        A[i] = temp[i];
    delete[] temp;
    delete[] counts;
}

void parallelIntSort(int *A, int sz, int max, int numThreads)
{
    int *counts = new int[max + 1]();
    int *temp = new int[sz];

#pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < sz; i++)
    {
        counts[A[i]]++;
    }

    for(int i = 0; i < max + 1; i++)
    {
        counts[i] += counts[i - 1];
    }
#pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < sz; i++)
        temp[--counts[A[i]]] = A[i];

#pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < sz; i++)
        A[i] = temp[i];

    delete[] temp;
    delete[] counts;
}

void newParallelIntSort(int *A, int sz, int max, int numThreads)
{
    int *counts = new int[max + 1]();
    int *temp = new int[sz];

    for(int i = 0; i < sz; i++)
    {
        counts[A[i]]++;
    }
    for(int i = 0; i < max + 1; i++)
    {
        counts[i] += counts[i - 1];
    }
#pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < max + 1; i++)
    {
        for(int j = 0; j < (counts[i] - counts[i - 1]); i++)
        {
            counts[i]--;
            temp[counts[i]] = i;
        }
    }
#pragma omp parallel for num_threads(numThreads)
    for(int i = 0; i < sz; i++)
        A[i] = temp[i];  

    delete[] temp;
    delete[] counts;
}

bool test(int *A, int n)
{
    for(int i = 0; i < n - 1; i++)
    {
        if(A[i] > A[i+1])
            return false;
    }
    return true;
}

int maxInteger(int *A, int n)
{
    int max = A[0];
    for(int i = 0; i < n; i++)
    {
        if(A[i] > max)
            max = A[i];
    }
    return max;
}

int main(int argc, char *argv[])
{
    int size, numThreads;
    cout << "Enter size of array: ";
    cin >> size;

    cout << "Input number of threads: ";
	cin >> numThreads;

    if (numThreads <= 0) {
#pragma omp parallel
		numThreads = omp_get_num_threads();
	}

    int *A = new int[size];
    int *seqA = new int[size];
    int *parA = new int[size];
    int *sortA = new int[size];
    int *newParA = new int[size];

    for(int i = 0; i < size; i++)
        A[i] = rand() % 100000;

    for(int i = 0; i < size; i++)
    {
        seqA[i] = A[i];
        parA[i] = A[i];
        sortA[i] = A[i];
        newParA[i] = A[i];
    }

    int max = maxInteger(seqA, size);

    double start = omp_get_wtime();
    intSort(seqA, size, max);
    double stop = omp_get_wtime();
    double seqTime = stop - start;
    cout << "-----Sequential-----" << endl;
    printf("Time (sec.) %.16f\n", seqTime);
    bool seqTest = test(seqA, size);
    cout << seqTest << endl;

    start = omp_get_wtime();
    parallelIntSort(parA, size, max, numThreads);
    stop = omp_get_wtime();
    double parTime = stop - start;
    cout << "-----Parallel-----" << endl;
    printf("Time (sec.) %.16f\n", parTime);
    printf("---------------\n");
    printf("Speedup is %0.16f, Efficiency is %.16f \n", (seqTime/parTime), ((seqTime/parTime)/numThreads));
    bool parTest = test(parA, size);
    if(parTest == 1)
        printf("Parallel is Sorted. \n");
    if(parTest == 0)
        printf("Parallel is not Sorted. \n");
   
    start = omp_get_wtime();
    std::sort(sortA, sortA + size);
    stop = omp_get_wtime();
    double sortTime = stop - start;
    cout << "-----Sort-----" << endl;
    printf("Time (sec.) %.16f\n", sortTime);
    printf("---------------\n");

    start = omp_get_wtime();
    newParallelIntSort(newParA, size, max, numThreads);
    stop = omp_get_wtime();
    double newParTime = stop - start;
    cout << "-----New Parallel-----" << endl;
    printf("Time (sec.) %.16f\n", newParTime);
    printf("---------------\n");
    printf("Speedup is %0.16f, Efficiency is %.16f \n", (seqTime/newParTime), ((seqTime/newParTime)/numThreads));
    bool newParTest = test(newParA, size);
    if(newParTest == 1)
        printf("New Parallel is Sorted. \n");
    if(newParTest == 0)
        printf("New Parallel is not Sorted. \n");

    delete[] seqA;
    delete[] parA;
    delete[] sortA;
    delete[] newParA;
    return 0;
}       