#include <cmath>
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

bool equalMatrices(double A[], double B[], double error, double c, double r)
{
    for(int i = 0; i < (r * c); i++)
    {
        if(A[i] - B[i] > error)
           return false;
    }
    return true;
}

int main(int argc, char *argv[]) 
{
	int commsz, rank, n, g, t;
    double d;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Set the seed of the RNG.  Use different seeds per process.
	srand(time(0) + rank);

	if (rank == 0) 
    {
		cout << "Input number of nodes: " << endl;
        cin >> n;
        cout << "Do you want a directed graph (1) or an undirected graph (2): " << endl;
        cin >> g;
        cout << "Node Density (Between 0 and 1): " << endl;
        cin >> d;
        cout << "Number of trial simulations: " << endl;
        cin >> t;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&g, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);

	
	double *A = new double[n * n];
	double *B = new double[n * n];
    double *C = new double[n * n];
	double *D = new double[n * n];
    
    if(rank == 0)
    {
        for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
            {
                A[i * n + j] = 0;
                B[i * n + j] = 0;
                C[i * n + j] = 0;
                D[i * n + j] = 0;
            }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    delete[] A;
    delete[] B;
    delete[] C;
    delete[] D;

	MPI_Finalize();

	return 0;
}
