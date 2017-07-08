#include "stdafx.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <random>
#include <math.h>
#include <immintrin.h>

using namespace std;

void populateMatrix(double **matA, int length_size);
void printMatrix(double **matX, int length_size);
void deAllocateMemory(double ** matX, int length_size);
void matrixMultiplicationParallelOptimized(double **matA, double **matB, double **matBT, double **matC, int length_size);
double fRand(double fMin, double fMax);
void printResults(double *time_taken_array, int iterations, int mat_size);
void transpose(double **matX, double **matXT, int length_size);
double getMatrixValue(double *vecX, double *vecY, int length_size);

int main()
{
	// This is calculated according to the initial mean and standard deviation
	// Denotes the number of times, the calculations should be sampled
	int iterations = 10;
	//parameters for matrix sizes
	int min_mat_size = 200;
	int max_mat_size = 2000;
	int increment = 200;

	//matrix A, B, C, BT
	double **matA, **matB, **matC, **matBT;
	//an array for storing excution time for mean and standard deviation calculation
	double *time_taken_array = new double[iterations];
	// varibale to facilitate the calcluation of time taken to matrix multiplication in each sample case
	double end_time = 0;
	double start_time = 0;

	printf("Matrix Size\tAvarage Time taken(s)\tStandard Deviation\tRequired Sample Size\t\n");

	for (int i = min_mat_size; i <= max_mat_size; i += increment)
	{
		// allocate memory for the matrixes
		matA = new double*[i];
		matB = new double*[i];
		matBT = new double*[i];
		matC = new double*[i];

		for (int j = 0; j < i; j++) {

			matA[j] = new double[i];
			matB[j] = new double[i];
			matBT[j] = new double[i];
			matC[j] = new double[i];
		}

		for (int k = 0; k < iterations; k++) {
			//populate the arrays
			populateMatrix(matA, i);
			populateMatrix(matB, i);

			//matrix multiplciation and time tracking
			start_time = omp_get_wtime();
			matrixMultiplicationParallelOptimized(matA, matB, matBT, matC, i);
			end_time = omp_get_wtime();
			time_taken_array[k] = end_time - start_time;

		}
		//print the results
		printResults(time_taken_array, iterations, i);

		//de-allocate the memoery used
		deAllocateMemory(matA, i);
		deAllocateMemory(matB, i);
		deAllocateMemory(matBT, i);
		deAllocateMemory(matC, i);

	}
	std::cin.get();
}

// Population of the given matrix using random values 
void populateMatrix(double **matA, int length_size)
{
	for (int c = 0; c < length_size; c++) {
		for (int r = 0; r < length_size; r++)
		{
			matA[c][r] = fRand(1, 1000);
		}
	}
}

// Generate a random double value
double fRand(double fMin, double fMax)
{
	//C++11 random function to generate the random number
	//seeding the random number generator 
	random_device rd;
	// random number generator 
	mt19937 mt(rd());
	// distribution of numbers 
	uniform_real_distribution<double> dist(fMin, fMax);

	return dist(mt);
}

// THis method is to print the generated matrix (testing only)
void printMatrix(double **matX, int length_size)
{
	for (int c = 0; c < length_size; c++) {
		for (int r = 0; r < length_size; r++)
		{
			cout << matX[c][r] << " ";
		}
		cout << endl;
	}
	cout << endl;
	cout << endl;
}

// Multiply the matrixs A and B and store the result in matrix C
void matrixMultiplicationParallelOptimized(double **matA, double **matB, double **matBT, double **matC, int length_size)
{
	//take the transpose of matrix B
	transpose(matB, matBT, length_size);
	//parellel code block
	#pragma omp parallel for
	for (int i = 0; i < length_size; i++) {
		for (int j = 0; j < length_size; j++) {
			matC[i][j] = getMatrixValue(matA[i], matBT[j], length_size);
		}
	}
 }

// Deallocation of memory
void deAllocateMemory(double ** matX, int length_size)
{
	for (int i = 0; i < length_size; ++i) {
		delete[] matX[i];
	}
}

//print the results 
void printResults(double *time_taken_array, int iterations, int mat_size)
{
	double mean_time = 0;
	double standard_deviation = 0;
	double sum_for_mean = 0;
	double sum_for_sd = 0;
	double required_samples = 0;

	for (int i = 0; i < iterations; ++i) {
		sum_for_mean += time_taken_array[i];
	}
	mean_time = sum_for_mean / iterations;

	for (int i = 0; i < iterations; ++i) {
		sum_for_sd += pow((time_taken_array[i] - mean_time), 2);
	}
	standard_deviation = sqrt(sum_for_sd / (iterations - 1));

	required_samples = pow((100 * 1.96 * standard_deviation) / (5 * mean_time), 2);

	printf("%d\t\t%4.6f\t\t%4.6f\t\t%4.2f\n", mat_size, mean_time, standard_deviation, required_samples);
}

// Take the transpose of the given matrix
void transpose(double **matX, double **matXT, int length_size) {
		//parallel code block
		#pragma omp parallel for
		for (int i = 0; i < length_size; i ++) {
			for (int j = 0; j < length_size; j ++) {
				matXT[i][j] = matX[j][i];
			}
		}
}

//get the sum of multiplicated values of given two rows
double getMatrixValue(double *vecX, double *vecY, int length_size) {
	double sum = 0.0;
	int block_size = 4;
	for (int i = 0; i < length_size; i+=block_size) {
		sum += vecX[i] * vecY[i] + vecX[i+1] * vecY[i+1] + vecX[i+2] * vecY[i+2] + vecX[i+3] * vecY[i+3];
	}
	return sum;
}
