#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <omp.h>

using namespace std;

bool compareMatrices(double *y, double *yp, int n, double tol) {
	bool retval = true;
	for (int i = 0; i < n; i++)
		retval = retval && (fabs(y[i] - yp[i]) < tol);

	return retval;
}

int main() {
	int numThreads, r, c;
	double start, end;
	srand(time(0));

	cout << "Input the coefficient array size r c: ";
	cin >> r >> c;

	cout << "Input number of threads: ";
	cin >> numThreads;

	if (numThreads <= 0) {
#pragma omp parallel
		numThreads = omp_get_num_threads();
	}

	double *A = new double[r * c];
	double *x = new double[r * c];
	double *y = new double[r];
	double *yp = new double[r];

    double min = -10.0;
    double max = 10.0;

	// Populate A and x.
	for (int i = 0; i < r * c; i++)
		A[i] = min + static_cast<double>(rand()) / RAND_MAX * (max - min);

	for (int i = 0; i < r * c; i++)
		x[i] = min + static_cast<double>(rand()) / RAND_MAX * (max - min);

	start = omp_get_wtime();
	for (int i = 0; i < r; i++) {
		y[i] = 0.0;
		for (int j = 0; j < c; j++)
			y[i] += A[i * c + j] * x[i * c + j];
	}
	end = omp_get_wtime();
    double seqTime = end - start;
	cout << "-----Sequential-----" << endl;
    printf("Time (sec.) %.16f\n", seqTime);


	start = omp_get_wtime();
#pragma omp parallel for num_threads(numThreads) shared(A, x, yp, r, c)
	for (int i = 0; i < r; i++) {
		yp[i] = 0.0;
		for (int j = 0; j < c; j++)
			yp[i] += A[i * c + j] * x[i * c + j];
	}
	end = omp_get_wtime();
    double parTime = end - start;
	cout << "-----Parallel-----" << endl;
    printf("Time (sec.) %.16f\n", parTime);
    printf("---------------\n");
    printf("Speedup is %0.16f, Efficiency is %.16f \n", (seqTime/parTime), ((seqTime/parTime)/numThreads));

	bool result = compareMatrices(y, yp, r, 0.0000001);
	if (result)
		printf("Matrices are equal. \n");
	else
		printf("Matrices are not equal. \n");

	delete[] A;
	delete[] x;
	delete[] y;
	delete[] yp;

	return 0;
}
