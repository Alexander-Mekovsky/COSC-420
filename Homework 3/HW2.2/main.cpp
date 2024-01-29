//Shell Sort Parallelization 2

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ctime>
#include <cmath>
#include <random>
#include <cstring>
#include <iostream>
#include <iomanip>

using namespace std;

bool order(int *arr, int s) 
{
    for (int i = 0; i < s - 1; i++) 
    {
        if (arr[i] > arr[i + 1])
            return false;
    }
    return true;
}

bool checkArrs(int *seq, int *para, int s) 
{
    for (int i = 0; i < s - 1; i++) 
    {
        if (seq[i] != para[i])
        {
            return false;
        }
    }
    return true;
}

void seqShellSort(int data[], int n) 
{
    int i, j, hCnt, h;
    int increments[20], k;
    for (h = 1, i = 0; h < n; i++) 
    {
        increments[i] = h;
        h = 3 * h + 1;
    }
    for (i--; i >= 0; i--) 
    {
        h = increments[i];
        for (hCnt = h; hCnt < 2 * h; hCnt++) 
        {
            for (j = hCnt; j < n; ) 
            { 
                double tmp = data[j];
                k = j;
                while (k - h >= 0 && tmp < data[k - h]) 
                {
                    data[k] = data[k - h];
                    k -= h;
                }
                data[k] = tmp;
                j += h;
            }
        }
    }
}

void parShellSort(int data[], int n, int commsz, int rank) 
{
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int local_n = (n / commsz); 
    int local_data[local_n];  
    MPI_Scatter(data, local_n, MPI_INT, local_data, local_n, MPI_INT, 0, MPI_COMM_WORLD);

    int i, j, hCnt, h;
    int increments[6], k;
    
    for(i = 0, h = commsz; h > 2; i++)
    {
        increments[i] = h;
        h = (1/3) * h;
    }
    for (i--; i >= 0; i--) 
    {
        h = increments[i];
            for (hCnt = h; hCnt < 2 * h; hCnt++) 
            { 
                for (j = hCnt; j < local_n; ) 
                {
                    double tmp = local_data[j];
                    k = j;
                    while (k - h >= 0 && tmp < local_data[k - h]) 
                    {
                        local_data[k] = local_data[k - h];
                        k -= h;
                    }
                    local_data[k] = tmp;
                    j += h;
                }
        }
    }
    
    MPI_Gather(local_data, local_n, MPI_INT, data, local_n, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
         int i, key, j;
        for (i = 1; i < n; i++) 
        {
            key = data[i];
            j = i - 1;
            while (j >= 0 && data[j] > key) 
            {
                data[j + 1] = data[j];
                j = j - 1;
            }
            data[j + 1] = key;
        }
    }
}

int main(void) 
{
    int commsz, rank, n;
    double starttime, endtime;
    
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &commsz);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    srand(time(0) * time(0));
    
    int *A = NULL;
    int *seqA = NULL;
    int *parA = NULL;
    
    double seqTime, paraTime;
    
    if (rank == 0) 
    {
        cout << "Input size of array: ";
        cin >> n;

        A = new int[n]; 
        seqA = new int[n];
        parA = new int[n];

        for (int i = 0; i < n; i++) 
        {
            A[i] = rand() % 100000000;
        }
        for (int i = 0; i < n; i++) 
        {
            seqA[i] = A[i];
            parA[i] = A[i];
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if (rank == 0) 
    {
        starttime = MPI_Wtime();

        seqShellSort(seqA, n);
        bool check = order(seqA, n);

        endtime = MPI_Wtime();
        seqTime = endtime - starttime;
        cout << "\n-----Sequential-----\n" << endl;
        printf("Shell Sort sequentially done in %.8f sec.\n", seqTime);
        if (check == true) 
        {
            printf("Array is sorted\n");
        }
        else 
        {
            printf("Array is not sorted\n");
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    if(rank ==0)
    {
        starttime = MPI_Wtime();
    }
        
    parShellSort(parA, n, commsz, rank);

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == 0) 
    {
        endtime = MPI_Wtime();
        bool check = order(parA, n);
        paraTime = endtime - starttime;
        cout << "\n-----Parallel-----\n" << endl;
        printf("Shell Sort parallel done in %.8f sec.\n", paraTime);
        if (check == true) 
        {
            printf("The array is sorted\n");
        } 
        else 
        {
            printf("The array is not sorted\n");
        }
        check = checkArrs(seqA, parA, n);
        if(check == true)
        {
            printf("Arrays are equal to eachother\n");
        }
        else
        {
            printf("Arrays are not equal to eachother\n");
        }
        printf("Speedup =  %.8f\n", seqTime / paraTime);
        printf("Efficiency = %.8f \n", (seqTime / paraTime) / commsz);
        
    }
    
    delete[] A;
    delete[] seqA;
    delete[] parA;
        
    MPI_Finalize();

    return 0;
}

