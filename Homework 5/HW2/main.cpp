#include <math.h>
#include <mpi.h>
#include <iostream>
#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char *argv[]) {
    int commsz, rank, n, maxIter = 1000;
    double tol = 0.00001;
    double starttime, endtime;
    double seqtime, partime;
    double *A = NULL;
    double *seqA = NULL;
    double *parA = NULL;
    double *seqE = NULL;
    double *parE = NULL;
    double *localA = NULL;
    double *lyp = NULL;
    
    MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Set the seed of the RNG.  Use different seeds per process.
	srand(time(0) + rank);

    if(rank == 0)
    {
        
        cout << "Enter size of the square matrix: ";
        cin >> n;
        
        A = new double[n * n];
        seqA = new double[n * n];
        parA = new double[n * n];
        seqE = new double[n];
        parE = new double[n];
    
        for(int i = 0; i < (n * n); i++)
        {
            A[i] = rand() % 10;
            seqA[i] = A[i];
            parA[i] = A[i];
        }
        for(int i = 0; i < n; i++)
        {
            seqE[i] = 1.0;
            parE[i] = 1.0;
        }
        
        starttime = MPI_Wtime();

        double *temp = new double[n];

        for (int iter = 0; iter < maxIter; iter++) 
        {
            double *product = new double[n];

            for (int i = 0; i < n; i++) 
            {
                for (int j = 0; j < n; j++) 
                {
                    product[i] += seqA[i * n + j] * seqE[j];
                }
            }
            double maxValue = product[0];
            for (int i = 1; i < n; i++) 
            {
                if (abs(product[i]) > abs(maxValue)) 
                {
                    maxValue = product[i];
                }
            }

            for (int i = 0; i < n; i++) 
            {
                seqE[i] = product[i] / maxValue;
            }

            double error = 0.0;
            for (int i = 0; i < n; i++) 
            {
                temp[i] = 0.0;
                for (int j = 0; j < n; j++) 
                {
                    temp[i] += seqA[i * n + j] * seqE[j];
                }
                error += temp[i] * seqE[i];
            }

            if (error < tol) 
            {
                cout << "Converged after " << iter + 1 << " iterations." << endl;
                break;
            }

            delete[] product;
        }
        
        double eigenValue = 0.0;
        for (int i = 0; i < n; i++) 
        {
            temp[i] = 0.0;
            for (int j = 0; j < n; j++) 
            {
                temp[i] += seqA[i * n + j] * seqE[j];
            }
            eigenValue += temp[i] * seqE[i];
        }

        cout << "Dominant Eigenvalue: " << eigenValue << endl;
        cout << "Eigenvector: ";
        for (int i = 0; i < n; i++) 
        {
            cout << seqE[i] << " ";
        }
        cout << endl;
        
        delete[] temp;
        
        endtime = MPI_Wtime();
        seqtime = endtime-starttime;
        cout << "\n-----Sequential-----\n" << endl;
        printf("Gauss-Jordan done in %.8f sec.\n", seqtime);
    }
    
    starttime = MPI_Wtime();
    
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&maxIter, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    localA = new double[n * n / commsz];
	lyp = new double[n / commsz];
	// Scatter the rows of A to respective process.
	MPI_Scatter(A, n * n / commsz, MPI_DOUBLE, localA, n * n / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    double *temp = new double[n];

        for (int iter = 0; iter < maxIter; iter++) 
        {
            double *product = new double[n];

            for (int i = 0; i < n; i++) 
            {
                for (int j = 0; j < n; j++) 
                {
                    product[i] += seqA[i * n + j] * seqE[j];
                }
            }
            double maxValue = product[0];
            for (int i = 1; i < n; i++) 
            {
                if (abs(product[i]) > abs(maxValue)) 
                {
                    maxValue = product[i];
                }
            }

            for (int i = 0; i < n; i++) 
            {
                seqE[i] = product[i] / maxValue;
            }

            double error = 0.0;
            for (int i = 0; i < n; i++) 
            {
                temp[i] = 0.0;
                for (int j = 0; j < n; j++) 
                {
                    temp[i] += seqA[i * n + j] * seqE[j];
                }
                error += temp[i] * seqE[i];
            }

            if (error < tol) 
            {
                cout << "Converged after " << iter + 1 << " iterations." << endl;
                break;
            }

            delete[] product;
        }
        
        double eigenValue = 0.0;
        for (int i = 0; i < n; i++) 
        {
            temp[i] = 0.0;
            for (int j = 0; j < n; j++) 
            {
                temp[i] += seqA[i * n + j] * seqE[j];
            }
            eigenValue += temp[i] * seqE[i];
        }
        
        delete[] temp;
    
    
    MPI_Gather(lyp, n / commsz, MPI_DOUBLE, parE, n / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if(rank == 0)
    {
        cout << "Dominant Eigenvalue: " << eigenValue << endl;
        cout << "Eigenvector: ";
        for (int i = 0; i < n; i++) 
        {
            cout << seqE[i] << " ";
        }
        cout << endl;
    }
    
    endtime = MPI_Wtime();
    partime = endtime-starttime;
    cout << "\n-----Parallel-----\n" << endl;
    printf("Gauss Jordan done in %.8f sec.\n",  partime);
    printf("Speedup = %.8f\n", seqtime/partime);
    printf("Efficiency = %.8f \n", (seqtime/partime)/commsz);

    delete[] A;
    delete[] seqA;
    delete[] parA;
    delete[] seqE;
    delete[] parE;
    
    MPI_Finalize();

    return 0;
}
