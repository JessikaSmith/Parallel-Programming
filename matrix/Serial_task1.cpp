#include <iostream>
#include <omp.h>
#include <stdlib.h>  
#include <random>
#include <algorithm> 

using namespace std;

mt19937 mt(1729);
normal_distribution<float> dist(-1, 1);

float** generateMatrix(unsigned N) {
	float** matrix = 0;
	matrix = new float*[N];
	for (int i = 0; i < N; i++) {
		matrix[i] = new float[N];
		for (int j = 0; j < N; j++) {
			matrix[i][j] = dist(mt);
			//cout << matrix[i][j] << ' ';
		}
		//cout << endl;
	}
	return matrix;
}


int main() {
	const unsigned N = 10000;
	const int num_of_time_samples = 2;
	float time[num_of_time_samples];
	//initializing by random numbers
	float** matrix = generateMatrix(N);
	double t_1, t_2, dt = 0;
	float row_max, row_min, avg = 0;
	float* matrix_max, *matrix_min;
	float diag_max, diag_min, diag_avg = 0;
	float min_m[N];
	float max_m[N];
	float avg_m[N];
	for (int t = 0; t < num_of_time_samples; t++) {
		diag_max = matrix[0][0];
		diag_min = matrix[0][0];
		t_1 = omp_get_wtime();
		for (int i = 0; i < N; i++) {
			avg = matrix[i][0];
			row_max = matrix[i][0];
			row_min = matrix[i][0];
			for (int j = 1; j < N; j++) {
				if (matrix[i][j] > row_max) row_max = matrix[i][j];
				if (matrix[i][j] < row_min) row_min = matrix[i][j];
				if (i == j) {
					if (matrix[i][j] > diag_max) diag_max = matrix[i][j];
					if (matrix[i][j] < diag_min) diag_min = matrix[i][j];
					diag_avg += matrix[i][j];
				}
				avg += matrix[i][j];
			}
			max_m[i] = row_max;
			min_m[i] = row_min;
			avg_m[i] = avg / N;
		}
		matrix_max = max_element(max_m, max_m + N);
		matrix_min = min_element(min_m, min_m + N);
		t_2 = omp_get_wtime();
		time[t] = t_2 - t_1;
	}
	sort(time, time + num_of_time_samples);
	cout << "For main diog max = " << diag_max << "; min =" << diag_min << "; avg = " << diag_avg/N << endl;
	cout << "Matrix min = " << *matrix_min << endl;
	cout << "Matrix max = " << *matrix_max << endl;
	cout << "Computational time = " << time[1] << endl;
}