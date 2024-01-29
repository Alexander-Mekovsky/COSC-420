//Alexander Mekovsky and Gavin Beauchamp

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ctime>
#include <cmath>
#include <random>
#include <cstring>

using namespace std;

void check(double seq, double para){
    if((seq-para)<0.000001){
        printf("Serial and Parallel are the same\n");
        return;
    }
    
     printf("Serial and Parallel are not the same\n");
}

int main(void) {
	int commsz, rank, n; 
    double seqSum = 0;
    double par1Sum = 0;
    double par2Sum = 0;
    double starttime, endtime;
    
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> rng(-1.0, 1.0);

	/* Start up MPI */
	MPI_Init(NULL, NULL);

	/* Get the number of processes */
	MPI_Comm_size(MPI_COMM_WORLD, &commsz);

	/* Get my rank among all the processes */
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        
    srand(time(0) * time(0));
    double seqtime, para1time, para2time;
    

    
    if(rank == 0) 
    {
		cout << "Input size of vectors: ";
		cin >> n;
	}
    
    double *v = NULL;
    double *w = NULL;
    double *seqV = NULL;
    double *seqW = NULL;
    double *par1V = NULL;
    double *par1W = NULL;
    double *par2V = NULL;
    double *par2W = NULL;
    if(rank == 0)
    {
        v = new double[n];
        w = new double[n];
        seqV = new double[n];
        seqW = new double[n];
        par1V = new double[n];
        par1W = new double[n];
        par2V = new double[n];
        par2W = new double[n];
        
        for (int i = 0; i < n; i++)
        {
            v[i] = rng(gen);
            w[i] = rng(gen);
        }
    
        for(int i = 0; i < n; i++)
        {  
            seqV[i] = v[i];
            seqW[i] = w[i];
            par1V[i] = v[i];
            par1W[i] = w[i];
            par2V[i] = v[i];
            par2W[i] = w[i];
        }   
    }
    if(rank == 0)
    {
        double temp = 0;
        starttime = MPI_Wtime();
        for(int i = 0; i < n; i++)
        {
            temp = v[i] * w[i];
            seqSum += temp;
        }
        endtime = MPI_Wtime();
        seqtime = endtime-starttime;
        cout << "\n-----Sequential-----\n" << endl;
        printf("Dot Product sequentially = %.8f done in %.8f sec.\n", seqSum, seqtime);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //Send and Receive


    if(rank == 0)
    {
        double number;
        starttime = MPI_Wtime();
        for(int i = 1; i < commsz; i++)
        {
            MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
        }

        for(int i = 0; i < n; i++)
        {
            MPI_Send(&par1V[i], 1, MPI_DOUBLE, i % (commsz-1) + 1, 0, MPI_COMM_WORLD);
            MPI_Send(&par1W[i], 1, MPI_DOUBLE, i % (commsz-1) + 1, 0, MPI_COMM_WORLD);
        }
        for(int i = 0; i < n; i++)
        {
            MPI_Recv(&number, 1, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            par1Sum += number;
        }
        endtime = MPI_Wtime();
        para1time = endtime-starttime;
        cout << "\n-----Parallel Point to Point-----\n" << endl;
        printf("Dot Product parallel with point-to-point = %.8f done in %.8f sec.\n", par1Sum, para1time);
        printf("The speedup for parallel with point-to-point = %.8f\n", seqtime/para1time);
        printf("The efficiency for parallel with point-to-point = %.8f \n", (seqtime/para1time)/commsz);
        check(seqSum, par1Sum);
    
    }
    double temp = 0;
    double a, b;
    for(int i = 1; i < commsz; i++)
    {
        if(i == rank)
        {
            MPI_Recv(&n, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }        
    }
    for(int i = 0; i < n; i++)
    {
        if(i % (commsz-1)+1 == rank)
        {
            MPI_Recv(&a, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            MPI_Recv(&b, 1, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            temp = a * b;
            MPI_Send(&temp, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //Broadcasting
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    double *temp1 = new double[n];
    double *temp2 = new double[n];
    if(rank == 0)
    {
        for(int i = 0; i < n; i++)
        {
            temp1[i] = par2V[i];
            temp2[i] = par2W[i];
        }
        starttime = MPI_Wtime();
    }
    MPI_Bcast(temp1, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(temp2, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double tempAdd = 0;
    
    for(int i = rank; i < n; i+= commsz)
    {
        tempAdd += (temp1[i] * temp2[i]);
    }
    
    MPI_Reduce(&tempAdd, &par2Sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
    if(rank == 0) 
    {
        endtime = MPI_Wtime();
        para2time = endtime-starttime;
        cout << "\n-----Parallel Collective Communication-----\n" << endl;
        printf("Dot Product parallel with collective communication = %.8f done in %.8f sec.\n", par2Sum, para2time);
        printf("The speedup for parallel with collective communication = %.8f\n", seqtime/para2time);
        printf("The efficiency for parallel with collective communication = %.8f\n", (seqtime/para2time)/commsz);
        check(seqSum, par2Sum);
    }

    
	MPI_Finalize();
    
    delete[] v;
    delete[] w;
    delete[] seqV;
    delete[] seqW;
    delete[] par1V;
    delete[] par1W;
    delete[] par2V;
    delete[] par2W;
    delete[] temp1;
    delete[] temp2;
    
	return 0;
}
