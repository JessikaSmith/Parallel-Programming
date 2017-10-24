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
		}
	}
	return matrix;
}

float* generateVector(unsigned N) {
	float* vec = new float[N];
	for (int i = 0; i < N; i++) {
		vec[i] = dist(mt);
	}
	return vec;
}

double computeParams(unsigned N, float** mat, unsigned num_of_time_samples,unsigned thrdsCount) {

	double t_1, t_2, dt = 0;
	double* time = new double[num_of_time_samples];

	float row_max, row_min, avg = 0;
	float* matrix_max, *matrix_min;
	float diag_max, diag_min, diag_avg = 0;
	float* min_m = new float[N];
	float* max_m = new float[N];
	float* avg_m = new float[N];
	for (int t = 0; t < num_of_time_samples; t++) {
		diag_max = mat[0][0];
		diag_min = mat[0][0];
		t_1 = omp_get_wtime();
		#pragma omp parallel for ordered schedule(static) shared(mat) \
			private(avg, row_max, row_min) num_threads(thrdsCount)
		for (int i = 0; i < N; i++) {
			avg = mat[i][0];
			row_max = mat[i][0];
			row_min = mat[i][0];
			for (int j = 1; j < N; j++) {
				if (mat[i][j] > row_max) {
					#pragma omp critical
					row_max = mat[i][j];
				}
				if (mat[i][j] < row_min) {
					#pragma omp critical
					row_min = mat[i][j];
				}
				if (i == j) {
					if (mat[i][j] > diag_max) {
						#pragma omp critical
						diag_max = mat[i][j];
					}
					if (mat[i][j] < diag_min) {
						#pragma omp critical
						diag_min = mat[i][j];
					}
					#pragma omp atomic
					diag_avg += mat[i][j];
				}
				#pragma omp atomic
				avg += mat[i][j];
			}
			max_m[i] = row_max;
			min_m[i] = row_min;
			avg_m[i] = avg / N;
			matrix_max = max_element(max_m, max_m + N);
			matrix_min = min_element(min_m, min_m + N);
		}
		t_2 = omp_get_wtime();
		time[t] = t_2 - t_1;	
		cout << "For main diog max = " << diag_max << "; min =" << diag_min << "; avg = " << diag_avg / N << endl;
		cout << "Matrix min = " << *matrix_min << endl;
		cout << "Matrix max = " << *matrix_max << endl;
	}

	double minim = time[0];
	for (int y = 1; y < num_of_time_samples; y++) {
		if (time[y] < minim) minim = time[y];
	}
	cout << "Min time " << minim << " for " << thrdsCount << " threads"<< endl;
	return minim;
}

double computeParamsVec(float* vec1, float* vec2, int n, unsigned num_of_time_samples, unsigned thrdsCount) {

	double t_1, t_2, dt = 0;
	double* time = new double[num_of_time_samples];

	for (int t = 0; t < num_of_time_samples; t++) {
		t_1 = omp_get_wtime();
		double res = 0;
		#pragma omp parallel for reduction(+:res) num_threads(thrdsCount)
		for (int i = 0; i < n; i++) {
			res += vec1[i] * vec2[i];
		}
		//cout << "Result: " << res << endl;
		t_2 = omp_get_wtime();
		time[t] = t_2 - t_1;
	}
	double minim = time[0];
		for (int y = 1; y < num_of_time_samples; y++) {
			if (time[y] < minim) minim = time[y];
		}
	cout << "Min time" << minim << " for " << thrdsCount << " threads" << endl;
	return minim;
}

int main() {
	const unsigned N = 10000;
	const unsigned num_of_time_samples = 9;

	//initializing matrix by random numbers
	cout << "==========  Speedup for matrix computations  ==========" << endl;
	float** mat = generateMatrix(N);
	const int max_thd = 16;
	double one = 1;
	for (int num_thd = 1; num_thd < max_thd+1; num_thd++) {
		if (num_thd == 1) {
			one = computeParams(N, mat, num_of_time_samples, num_thd);
			cout << 1 << endl;
		} else
			if (num_thd % 2 != 0) continue;
		//cout << num_thd << endl;
		cout << one/computeParams(N, mat, num_of_time_samples, num_thd) << endl;
	}
	//initialize 2 vectors
	unsigned vec_size = 100000000;
	float* vec1 = generateVector(vec_size);
	float* vec2 = generateVector(vec_size);
	cout << "==========  Speedup for vector computations  ==========" << endl;
	for (int num_thd = 1; num_thd < max_thd + 1; num_thd++) {
		if (num_thd == 1) {
			one = computeParamsVec(vec1, vec2, vec_size, num_of_time_samples, num_thd);
			cout << 1 << endl;
		}
		else
			if (num_thd % 2 != 0) continue;
		//cout << num_thd << endl;
		cout << one / computeParamsVec(vec1, vec2, vec_size, num_of_time_samples, num_thd) << endl;
	}
}
