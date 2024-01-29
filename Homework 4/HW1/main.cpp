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

int main(int argc, char *argv[]) {
	int commsz, rank, r1, c1, r2, c2;
    bool add = false, sub = true, mult = false, trans = false;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Set the seed of the RNG.  Use different seeds per process.
	srand(time(0) + rank);

	if (rank == 0) 
    {
		cout << "Input r1, number of matrix rows will be m*commsz: ";
		cin >> r1;
        r1 *= commsz;
		cout << "Input c1, number of matrix columns will be m*commsz: ";
		cin >> c1;
        c1 *= commsz;
        cout << "Input r2, number of matrix rows will be m*commsz: ";
		cin >> r2;
        r2 *= commsz;
        cout << "Input c2, number of matrix columns will be m*commsz: ";
		cin >> c2;
        c2 *= commsz;
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    
    MPI_Bcast(&r1, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&c1, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&c2, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	double *A = new double[r1 * c1];
	double *B = new double[r2 * c2];
    double *seqAdd = new double[r1 * c1];
    double *parAdd = new double[r1 * c1];
    double *seqSub = new double[r1 * c1];
    double *parSub = new double[r1 * c1];
    double *seqMult = new double[r2 * c1];
    double *parMult = new double[r2 * c1];
    //double *seqTrans = new double[r1 * c1];
    //double *parTrans = new double[r1 * c1];
    
    if(rank == 0)
    {
        for(int i = 0; i < r1; i++)
			for(int j = 0; j < c1; j++)
				A[i * c1 + j] = rand() % 1000;
        for(int i = 0; i < r2; i++)
			for(int j = 0; j < c2; j++)
				B[i * c2 + j] = rand() % 1000;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Addition function
    
    if(add == true)
    {
        if(rank == 0)
        {
            if(r1 != r2 || c1 != c2)
            {
                printf("Matrices are not the same size.");
                return -1;
            }
            for(int i = 0; i < r2; i++)
            {
                seqAdd[i] = 0;
                for(int j = 0; j < c2; j++)
                    seqAdd[i * c2 + j] = A[i * c2 + j] + B[i * c2 + j];

            }
        }
    
        MPI_Bcast(A, c1 * r1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(B, c2 * r2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        double *localA = new double[r1 * c1];
        double *localB = new double[r2 * c2];
        double *localPar = new double[r1 * c1];
    
        MPI_Scatter(A, c1 * r1 / commsz, MPI_DOUBLE, localA, c1 * r1 / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(B, c2 * r2 / commsz, MPI_DOUBLE, localA, c2 * r2 / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        for (int i = 0; i < r1 / commsz; i++) 
        {
            parAdd[i] = 0;
            for(int j = 0; j < c1; j++)
                localPar[i] += localA[i * c1 + j] + localB[i * c1 + j];
        }
        
        MPI_Gather(parAdd, r1 / commsz, MPI_DOUBLE, localPar, r1 / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        if(rank == 0)
        {
            bool test = equalMatrices(seqAdd, parAdd, 0.00001, c1, r1);
            if(test == false)
                cout << "The matrices are not equal." << endl;
            if(test == true)
                cout << "The matrices are equal." << endl;
        }
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Subtraction function
    if(sub == true)
    {
        if(rank == 0)
        {
            if(r1 != r2 || c1 != c2)
            {
                printf("Matrices are not the same size.");
                return -1;
            }
            for(int i = 0; i < r2; i++)
            {
                seqSub[i] = 0;
                for(int j = 0; j < c2; j++)
                    seqSub[i * c2 + j] = A[i * c2 + j] - B[i * c2 + j];

            }
        }
    
        MPI_Bcast(A, c1 * r1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(B, c2 * r2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        double *localA = new double[r1 * c1 / commsz];
        double *localB = new double[r2 * c2 / commsz];
        double *localPar = new double[r1 * c1 / commsz];
    
        MPI_Scatter(A, c1 * r1 / commsz, MPI_DOUBLE, localA, c1 * r1 / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(B, c2 * r2 / commsz, MPI_DOUBLE, localA, c2 * r2 / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
        for (int i = 0; i < r1 / commsz; i++) 
        {
            parAdd[i] = 0;
            for(int j = 0; j < c1; j++)
                localPar[i] += localA[i * c1 + j] - localB[i * c1 + j];
        }
        
        MPI_Gather(parSub, r1 / commsz, MPI_DOUBLE, localPar, r1 / commsz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        MPI_Barrier(MPI_COMM_WORLD);
        
        if(rank == 0)
        {
            bool test = equalMatrices(seqSub, parSub, 0.00001, c1, r1);
            if(test == false)
                cout << "The matrices are not equal." << endl;
            if(test == true)
                cout << "The matrices are equal." << endl;
        }
    }
    
    delete[] A;
    delete[] B;
    delete[] seqAdd;
    delete[] parAdd;
    delete[] seqSub;
    delete[] parSub;
    delete[] seqMult;
    delete[] parMult;
    //delete[] seqTrans;
    //delete[] parTrans;

	MPI_Finalize();

	return 0;
}

