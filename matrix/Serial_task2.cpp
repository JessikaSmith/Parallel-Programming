#include <iostream>
#include <omp.h>
#include <stdlib.h>  

using namespace std;

int main(){
	int i, chunk, iter_num;
	float a[1000], b[1000];
	float result = 0.0;
	int n = 100;
	double t_1, t_2, t_11;
    double dt = 0;
	chunk = 10;
	iter_num = 1000;
	for (i = 0; i < n; i++){
		a[i] = i * 1.0;
		b[i] = i * 2.0;
	}
	for (int r = 1; r < 17; r++) {
		#define thrdsCount = r
		dt = 0;
		for (int ti = 0; ti < iter_num; ti++) {
			result = 0.0;
			t_1 = omp_get_wtime();
			#pragma omp parallel for default(shared) private(i) schedule(static, chunk) reduction(+:result)
			for (i = 0; i < n; i++) result += (a[i] * b[i]);
			t_2 = omp_get_wtime();
			dt += t_2 - t_1;
		}
		if (r == 1) {
			t_11 = dt/iter_num;
		}
		//cout << "Measured time:" << dt/iter_num << endl;
		cout << "Relative SpeedUp, p = " << r <<": "<< t_11 / (dt/iter_num) << endl;
	}
}