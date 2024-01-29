#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <random>
#include <cstring>
#include <algorithm>

using namespace std;

template<class T>
void merge(T A[], T Temp[], int startA, int startB, int end) 
{
    int aptr = startA;
    int bptr = startB;
    int i = startA;
    while (aptr < startB && bptr <= end)
        if (A[aptr] < A[bptr])
            Temp[i++] = A[aptr++];
        else
            Temp[i++] = A[bptr++];
    while (aptr < startB)
        Temp[i++] = A[aptr++];
    while (bptr <= end)
        Temp[i++] = A[bptr++];
    for (i = startA; i <= end; i++)
        A[i] = Temp[i];
}

template<class T>
void mergeSort(T A[], T Temp[], int start, int end) 
{
    if (start < end) 
    {
        int mid = (start + end) / 2;
        mergeSort(A, Temp, start, mid);
        mergeSort(A, Temp, mid + 1, end);
        merge(A, Temp, start, mid + 1, end);
    }
}

template<class T>
void mergeSort(T A[], int size) 
{
    T *Temp = new T[size];
    mergeSort(A, Temp, 0, size - 1);
    delete [] Temp;
}

template<class T>
void parallelMerge(T A[], T Temp[], int startA, int startB, int end) 
{
    int aptr = startA;
    int bptr = startB;
    int i = startA;
    while (aptr < startB && bptr <= end)
        if (A[aptr] < A[bptr])
            Temp[i++] = A[aptr++];
        else
            Temp[i++] = A[bptr++];
    while (aptr < startB)
        Temp[i++] = A[aptr++];
    while (bptr <= end)
        Temp[i++] = A[bptr++];
    for (i = startA; i <= end; i++)
        A[i] = Temp[i];
}

template<class T>
void parallelMergeSort(T A[], T Temp[], int start, int end) 
{
    if (start < end) 
    {
        int mid = (start + end) / 2;
        #pragma omp task 
        parallelMergeSort(A, Temp, start, mid);
        #pragma omp task
        parallelMergeSort(A, Temp, mid + 1, end);
        #pragma omp taskwait
        parallelMerge(A, Temp, start, mid + 1, end);
    }
}

template<class T>
void parallelMergeSort(T A[], int size, int numThreads) 
{
    T *Temp = new T[size];
    #pragma omp parallel num_threads(numThreads)
    {
        #pragma omp single
        parallelMergeSort(A, Temp, 0, size - 1);
    }
    delete [] Temp;
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

template<class T>
bool sorted(T A[], int s)
{
    bool sorted = true;
    for(int i = 0; i < s - 1; i++)
    {
        if(A[i + 1] < A[i])
        {
            sorted = false;
            break;
        }
    }
    return sorted;
}

int main()
{
    int numThreads;
    int s;
	double start, end;
	srand(time(0));

	cout << "Input the size of the array. ";
	cin >> s;

    int *A = new int[s];
    int *seqA = new int[s];
    int *parA = new int[s];
    int *insA = new int[s];
    int *sortA = new int[s];

	cout << "Input number of threads: ";
	cin >> numThreads;

	if (numThreads <= 0) {
#pragma omp parallel
		numThreads = omp_get_num_threads();
	}
    for(long int i = 0; i < s; i++)
    {
        A[i] = rand();
        seqA[i] = A[i];
        parA[i] = A[i];
        insA[i] = A[i];
        sortA[i] = A[i];
    }
    start = omp_get_wtime();
    mergeSort(seqA, s);
    end = omp_get_wtime();
    double seqTime = end - start;
    cout << "-----Sequential Merge Sort-----" << endl;
    printf("Time (sec.) %.16f\n", seqTime);
    
    start = omp_get_wtime();
    parallelMergeSort(parA, s, numThreads);
    end = omp_get_wtime();
    double parTime = end - start;
    cout << "-----Parallel Merge Sort----" << endl;
    printf("Time (sec.) %.16f\n", parTime);
    printf("---------------\n");
    printf("Speedup is %0.16f, Efficiency is %.16f \n", (seqTime/parTime), ((seqTime/parTime)/numThreads));

    start = omp_get_wtime();
    insertionSort(insA, s);
    end = omp_get_wtime();
    double insTime = end - start;
    cout << "-----Insertion Sort-----" << endl;
    printf("Time (sec.) %.16f\n", insTime);

    start = omp_get_wtime();
    std::sort(sortA, sortA + s);
    end = omp_get_wtime();
    double sortTime = end - start;
    cout << "-----Sort Algorithm-----" << endl;
    printf("Time (sec.) %.16f\n", sortTime);
    
    bool seqSorted = sorted(seqA, s);
    bool parSorted = sorted(parA, s);
    bool insSorted = sorted(insA, s);
    bool sortSorted = sorted(sortA, s);
    
    if(seqSorted)
        cout << "Sequential is sorted." << endl;
    else
        cout << "Sequential is not sorted." << endl;
    if(parSorted)
        cout << "Parallel is sorted." << endl;
    else
        cout << "Parallel is not sorted." << endl;
    if(insSorted)
        cout << "Insertion is sorted." << endl;
    else
        cout << "Insertion is not sorted." << endl;
    if(sortSorted)
        cout << "Sort Algorithm is sorted." << endl;
    else
        cout << "Sort Algorithm is not sorted." << endl;

    return 0;
}
