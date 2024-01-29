 //Alexander Mekovsky and Gavin Beauchamp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ctime>

using namespace std;

int main(void) {
	int number;
	int commsz; 
	int rank;
    int temp;
	int parSum = 0;
    int seqSum = 0;

	/* Start up MPI */
	MPI_Init(NULL, NULL);

	/* Get the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);

	/* Get my rank among all the processes */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    srand(time(0) + rank);
    
    int n = rand() % 10000;
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) 
    {
        for(int i = 1; i < commsz; i++)
        {
            //number = rand() % 10000;
            if(i % 2 == 0)
            {
                seqSum += (number * 2);
            }
            else
            {
                seqSum += (number * -1);
            }
            MPI_Send(&number, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }
        
       for(int i = 1; i < commsz; i++)
    {
           MPI_Recv(&number, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//           parSum += number;        
        
    }

        MPI_Reduce(&number, &parSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        printf("Sequential sum = %d\n", seqSum);
        printf("Parallel sum = %d\n", parSum);
	} 
	else 
    {
        MPI_Recv(&temp, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if(rank % 2 == 0)
        {
            temp *= 2;
        }
        else
        {
            temp *= -1;
        }
        MPI_Send(&temp, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
	}
	MPI_Finalize();

	return 0;
}
