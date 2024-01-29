//Alexander Mekovsky and Gavin Beauchamp

#include <stdio.h>
#include <string.h>
#include <mpi.h>

int main(void) {
	int commsz; 
	int my_rank; 
	int n;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Send(&my_rank, 1, MPI_INT, (my_rank+1) % commsz, 0, MPI_COMM_WORLD);
    //MPI_Recv(&n, 1, MPI_INT, (rank==0) ? commsz - 1 : rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(my_rank == 0){
        MPI_Recv(&n, 1, MPI_INT, commsz -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    else{
        MPI_Recv(&n, 1, MPI_INT, my_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);  
    }
    printf("Process %d received message %d.\n", my_rank, n);
    
	MPI_Finalize();

	return 0;
} 
