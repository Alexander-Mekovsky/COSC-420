#include <math.h>
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

void printMatrix(double* A, int r, int c) 
{
    for (int i = 0; i < r; i++) 
    {
        for (int j = 0; j < c; j++) 
        {
            cout << A[i * c + j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

int main(int argc, char *argv[]) {
	int commsz, rank, r, c, n;
    double starttime, endtime;
    double seqtime, partime;
    double *A = NULL;
    double *seqA = NULL;
    double *parA = NULL;

	/* Start up MPI */
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Set the seed of the RNG.  Use different seeds per process.
	srand(time(0) + rank);
    
    if(rank == 0)
    {
        cout << "Input rows of array: ";
        cin >> r;
        cout << "Input columns of array: ";
        cin >> c;
        
        A = new double[r * c];
        seqA = new double[r * c];
        parA = new double[r * c];
        
        n = (r * c);
    
        for(int i = 0; i < n; i++)
        {
            A[i] = rand() % 10;
            seqA[i] = A[i];
            parA[i] = A[i];
        }
        
        //cout << "Original Sequential Matrix: " << endl;
        //printMatrix(seqA, r, c);
        
        starttime = MPI_Wtime();
        
        for (int i = 0; i < r; i++) 
        {
            int maxR = i;
            for (int j = i + 1; j < r; j++) 
            {
                if (abs(seqA[j * c + i]) > abs(seqA[maxR * c + i])) 
                {
                    maxR = j;
                }
            }
        
            for (int k = i; k <= c; k++) 
            {
                double temp = seqA[i * c + k];
                seqA[i * c + k] = seqA[maxR * c + k];
                seqA[maxR * c + k] = temp;
            }

            double pivot = seqA[i * c + i];
            for (int k = i + 1; k < c; k++) 
            {
                seqA[i * c + k] /= pivot;
            }
            seqA[i * c + i] = 1.0;

            for (int j = 0; j < r; j++) 
            {
                if (j != i) 
                {
                    double factor = seqA[j * c + i];
                    for (int k = i; k < c; k++) 
                    {
                        seqA[j * c + k] -= factor * seqA[i * c + k];
                    }
                }
            }
        }
        
        //cout << "Row-Reduced Sequential Matrix: " << endl;
        //printMatrix(seqA, r, c);
        
        endtime = MPI_Wtime();
        seqtime = endtime-starttime;
        cout << "\n-----Sequential-----\n" << endl;
        printf("Gauss-Jordan done in %.8f sec.\n", seqtime);
    }
    
    seqtime = endtime-starttime;
    
    MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int process = r / commsz;
    double *buffer = new double[c];
    
    starttime = MPI_Wtime();

    MPI_Bcast(&r, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Scatter(parA, process * c, MPI_DOUBLE, buffer, process * c, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    for (int i = 0; i < r; i++) 
        {
            int maxR = i;
            for (int j = i + 1; j < r; j++) 
            {
                if (abs(parA[j * c + i]) > abs(parA[maxR * c + i])) 
                {
                    maxR = j;
                }
            }
            
            MPI_Gather(&parA[maxR * c], c, MPI_DOUBLE, buffer, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);

            // Broadcast the maxR-th row to all processes
            MPI_Bcast(buffer, c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
            for (int k = i; k <= c; k++) 
            {
                double temp = seqA[i * c + k];
                parA[i * c + k] = parA[maxR * c + k];
                parA[maxR * c + k] = temp;
            }

            double pivot = parA[i * c + i];
            for (int k = i + 1; k < c; k++) 
            {
                parA[i * c + k] /= pivot;
            }
            parA[i * c + i] = 1.0;

            for (int j = 0; j < r; j++) 
            {
                if (j != i) 
                {
                    double factor = parA[j * c + i];
                    for (int k = i; k < c; k++) 
                    {
                        parA[j * c + k] -= factor * parA[i * c + k];
                    }
                }
            }
        }
        
    MPI_Gather(buffer, process * c, MPI_DOUBLE, parA, process * c, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    endtime = MPI_Wtime();
    partime = endtime-starttime;
    cout << "\n-----Parallel-----\n" << endl;
    printf("Gauss Jordan done in %.8f sec.\n",  partime);
    printf("Speedup = %.8f\n", seqtime/partime);
    printf("Efficiency = %.8f \n", (seqtime/partime)/commsz);

    
    delete[] A;
    delete[] seqA;
    delete[] parA;
    delete[] buffer;

	MPI_Finalize();
    

	return 0;
}
