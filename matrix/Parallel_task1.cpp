#include <iostream>
#include <omp.h>
#include <stdlib.h>  

using namespace std;

#define thrdsCount = 16

void generateMatrix(int N, int M, int mat[N][M]) {
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < M; j++) {
			mat[i][j] = rand() % 1000 + 1000;
		}
	}
}


int main() {
	//creating NxN matrix
	const int n = 100;
	//initializing by random numbers
	int matrix[n][n];
	double t_1, t_2, dt;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			matrix[i][j] = rand() % 1000 + 1000;
		}
	}
	int row_max, row_min, avg = 0;
	int matrix_max, matrix_min;
	int diag_max, diag_min, diag_avg = 0;
	//t_1 = omp_get_wtime();
	matrix_max = matrix[0][0];
	matrix_min = matrix[0][0];
	diag_max = matrix[0][0];
	diag_min = matrix[0][0];
	for (int i = 0; i < n; i++) {
		avg = matrix[i][0];
		row_max = matrix[i][0];
		row_min = matrix[i][0];
		for (int j = 1; j < n; j++) {
			if (matrix[i][j] > row_max) row_max = matrix[i][j];
			if (matrix[i][j] < row_min) row_min = matrix[i][j];
			if (i == j) {
				if (matrix[i][j] > diag_max) diag_max = matrix[i][j];
				if (matrix[i][j] < diag_min) diag_min = matrix[i][j];
				diag_avg += matrix[i][j];
			}
			avg += matrix[i][j];
		}
		if (row_max > matrix_max) matrix_max = row_max;
		if (row_min < matrix_min) matrix_min = row_min;
		cout << "For row " << i + 1 << " max = " << row_max << "; min = " << row_min << "; avg = " << avg / n << endl;
	}
	cout << "For main diog max = " << diag_max << "; min =" << diag_min << "; avg = " << diag_avg / n << endl;
	cout << "Matrix min = " << matrix_min << endl;
	cout << "Matrix max = " << matrix_max << endl;
}