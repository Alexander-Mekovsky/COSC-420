#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include <random>
#include <cstring>

using namespace std;

void seqRowBackSub(double **A, double *b, double *x, int s)
{
    for(int i = s - 1; i >= 0; i--)
    {
        x[i] = b[i];
        for(int j = i + 1; j <= s - 1; j++)
        {
            x[i] = x[i] - (A[i][j] * x[i]);
        }
        x[i] = x[i]/A[i][i];
    }
}

void parRowBackSub(double **A, double *b, double *x, int s, int numThreads)
{
#pragma omp parallel for num_threads(numThreads) shared(A, b, x, s)
    for(int i = s - 1; i >= 0; i--)
    {
        x[i] = b[i];
        for(int j = i + 1; j <= s - 1; j++)
        {
            x[i] = x[i] - (A[i][j] * x[i]);
        }
        x[i] = x[i]/A[i][i];
    }
}

void seqColBackSub(double **A, double *b, double *x, int s)
{
    for(int i = 0; i < s - 1; i++)
    {
        x[i] = b[i];
    }
    for(int i = s - 1; i >= 1; i--)
    {
        x[i] = x[i] / A[i][i];
        for(int j = 0; j < i - 1; j++)
        {
            x[j] = x[j] - (A[j][i] * x[j]);
        }
    }
}

void parColBackSub(double **A, double *b, double *x, int s, int numThreads)
{
    for(int i = 0; i < s - 1; i++)
    {
        x[i] = b[i];
    }
#pragma omp parallel for num_threads(numThreads) shared(A, b, x, s)
    for(int i = s - 1; i >= 1; i--)
    {
        x[i] = x[i] / A[i][i];
        for(int j = 0; j < i - 1; j++)
        {
            x[j] = x[j] - (A[j][i] * x[j]);
        }
    }
}

bool compareRuns(double **seq, double **par, int s, double tolerance)
{
    bool result = true;
    for(int i = 0; i < s; i++)
    {
        for(int j = 0; j < s; j++)
        {
            if(fabs(seq[i][j] - par[i][j]) > tolerance)
            {
                result = false;
                break;
            }
        }
    }
    return result;
}


int main()
{
    int numThreads, s;
	double start, end;
	srand(time(0));

	cout << "Input the coefficient square array size: ";
	cin >> s;

	cout << "Input number of threads: ";
	cin >> numThreads;

	if (numThreads <= 0) {
#pragma omp parallel
		numThreads = omp_get_num_threads();
	}

    double **A = new double*[s * s];

    double *b = new double[s];
    double *x = new double[s];

    double **seqRow = new double*[s * s];
    double **seqCol = new double*[s * s];
    double **parRow = new double*[s * s];
    double **parCol = new double*[s * s];

    for(int i = 0; i < s; i++)
    {
        A[i] = new double[s];
        seqRow[i] = new double[s];
        seqCol[i] = new double[s];
        parRow[i] = new double[s];
        parCol[i] = new double[s];
    }


    for(int i = 0; i < s; i++)
    {
        for(int j = 0; j < s; j++)
        {
            A[i][j] = (rand() % 1000) + 1;
        }
    }

    for(int i = 0; i < s; i++)
    {
        for(int j = 0; j < s; j++)
        {
            seqRow[i][j] = A[i][j];
            seqCol[i][j] = A[i][j];
            parRow[i][j] = A[i][j];
            parCol[i][j] = A[i][j];
        }
    }

	for(int i = 0; i < s; i++)
		b[i] = (rand() % 1000);

    //memcpy(seqRow, A, s * s * sizeof(double));
    //memcpy(seqCol, A, s * s * sizeof(double));
    //memcpy(parRow, A, s * s * sizeof(double));
    //memcpy(parCol, A, s * s * sizeof(double));

    start = omp_get_wtime();
    seqRowBackSub(seqRow, b, x, s);
    end = omp_get_wtime();
    double seqRowTime = end - start;
    cout << "-----Sequential Row-----" << endl;
    printf("Time (sec.) %.16f\n", seqRowTime);

    start = omp_get_wtime();
    seqColBackSub(seqCol, b, x, s);
    end = omp_get_wtime();
    double seqColTime = end - start;
    cout << "-----Sequential Column-----" << endl;
    printf("Time (sec.) %.16f\n", seqColTime);

    start = omp_get_wtime();
    parRowBackSub(parRow, b, x, s, numThreads);
    end = omp_get_wtime();
    double parRowTime = end - start;
    cout << "-----Parallel Row-----" << endl;
    printf("Time (sec.) %.16f\n", parRowTime);
    printf("---------------\n");
    printf("Speedup is %0.16f, Efficiency is %.16f \n", (seqRowTime/parRowTime), ((seqRowTime/parRowTime)/numThreads));

    start = omp_get_wtime();
    parColBackSub(parCol, b, x, s, numThreads);
    end = omp_get_wtime();
    double parColTime = end - start;
    cout << "-----Parallel Column-----" << endl;
    printf("Time (sec.) %.16f\n", parColTime);
    printf("---------------\n");
    printf("Speedup is %0.16f, Efficiency is %.16f \n", (seqColTime/parColTime), ((seqColTime/parColTime)/numThreads));

    bool rowResult = compareRuns(seqRow, parRow, s, 0.0000001);
    bool colResult = compareRuns(seqCol, parCol, s, 0.0000001);

    if(rowResult)
        cout << "Both row matrices are equal." << endl;
    else
        cout << "Both row matrices are not equal." << endl;

    if(colResult)
        cout << "Both col matrices are equal." << endl;
    else
        cout << "Both col matrices are not equal." << endl;

    delete[] A;
	delete[] seqRow;
	delete[] seqCol;
	delete[] parRow;
    delete[] parCol;

    return 0;
}