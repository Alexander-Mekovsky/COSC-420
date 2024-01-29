//Alexander Mekovsky and Gavin Beauchamp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ctime>
#include <cmath>

bool isPrime(int n)
{
    for(int i = 2; i < sqrt(n); i++)
        if(n % i == 0)
            return false;
    return true;
}

using namespace std;

int main(void) {
	int commsz; 
	int rank;
    int n;
    int numPrimes = 0;
    bool prime;
    int total = 0;

	/* Start up MPI */
	MPI_Init(NULL, NULL);

	/* Get the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);

	/* Get my rank among all the processes */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    srand(time(0) + rank);

    if (rank == 0) {
		cout << "Input number of numbers: ";
		cin >> n;
	}
	// Broadcast array size to all processes.
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int *A = new int[n];
    int *AVals = new int[n];
    for(int i = 0; i < n; i++)
        A[i] = (rand() % 999998) + 2;
    //MPI_Bcast(&A, n, MPI_INT, 0, MPI_COMM_WORLD);
    for(int i = 1; i < commsz; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(isPrime(A[i])){
                    AVals[i] = 1;
                }
                else{
                    AVals[i] = 0;
                }   
            }
        }
    
    MPI_Reduce(A, AVals, n, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0) 
    {
        for(int i = 1; i < commsz; i++)
        {
            for(int j = 0; j < n; j++)
            {
                if(AVals[i] == 1)
                    numPrimes++;
            }            
        }
        printf("Number of prime numbers: %f\n", (double)numPrimes/n);
	} 
//	else 
//    {
//        for(int i = 0; i < n; i++)
//        {
//            numPrimes = 0;
//            prime = isPrime(A[i]);
//            if(prime == true)
//                numPrimes++;
//            MPI_Send(&numPrimes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//       }
//	}
	
	delete[] A;
    delete[] AVals;
	MPI_Finalize();

	return 0;
}
