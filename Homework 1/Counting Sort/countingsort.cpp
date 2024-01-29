#include <iostream>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <cstring>
#include <algorithm>

using namespace std;

void countSort(int *A, int n)
{
    int i, j, count;
    int *temp = new int[n];

    for(i = 0; i < n; i++)
    {
        count = 0;
        for(j = 0; j < n; j++)
        {
            if(A[j] < A[i])
                count++;
            else if(A[j] == A[i] && j < i)
                count++;
        }
        temp[count] = A[i];
    }
    memcpy(A, temp, n*sizeof(int));
    delete[] temp;
}

void parallelCountSort(int *A, int n, int numThreads)
{
    int i, j, count;
    int *temp = new int[n];
#pragma omp parallel for num_threads(numThreads) private(i, j, count)
    for(i = 0; i < n; i++)
    {
        count = 0;
        for(j = 0; j < n; j++)
        {
            if(A[j] < A[i])
                count++;
            else if(A[j] == A[i] && j < i)
                count++;
        }
        temp[count] = A[i];
    }
#pragma omp parallel for
    for(int i = 0; i < n; i++)
        A[i] = temp[i];
    delete[] temp;
}

bool test(int *A, int n)
{
    for(int i = 0; i < n - 1; i++)
    {
        if(A[i] > A[i + 1])
            return false;
    }
    return true;
}

void insertionSort(int *A, int n)
{
    int i, key, j;
    for (i = 1; i < n; i++) {
        key = A[i];
        j = i - 1;
        while (j >= 0 && A[j] > key) {
            A[j + 1] = A[j];
            j = j - 1;
        }
        A[j + 1] = key;
    }
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
    int *insA = new int[size];
    int *sortA = new int[size];

    for(int i = 0; i < size; i++)
        A[i] = rand();

    for(int i = 0; i < size; i++)
    {
        seqA[i] = A[i];
        parA[i] = A[i];
        insA[i] = A[i];
        sortA[i] = A[i];
    }

    double start = omp_get_wtime();
    countSort(seqA, size);
    double stop = omp_get_wtime();
    double seqTime = stop - start;
    cout << "-----Sequential-----" << endl;
    printf("Time (sec.) %.16f\n", seqTime);
    bool seqTest = test(seqA, size);
    if(seqTest == 1)
        printf("Sequential is Sorted. \n");
    if(seqTest == 0)
        printf("Sequential is not Sorted. \n");

    cout << "--------------" << endl;

    start = omp_get_wtime();
    parallelCountSort(parA, size, numThreads);
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
    insertionSort(insA, size);
    stop = omp_get_wtime();
    double insTime = stop - start;
    cout << "-----Insertion Sort-----" << endl;
    printf("Time (sec.) %.16f\n", insTime);
    printf("---------------\n");

    start = omp_get_wtime();
    std::sort(sortA, sortA + size);
    stop = omp_get_wtime();
    double sortTime = stop - start;
    cout << "-----Sort-----" << endl;
    printf("Time (sec.) %.16f\n", sortTime);
    printf("---------------\n");

    delete[] A;
    delete[] seqA;
    delete[] parA;
    delete[] insA;
    delete[] sortA;
    return 0;
}