#include <assert.h>
#include <stdio.h>

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions and utilities to work with CUDA
#include <helper_cuda.h>
#include <helper_functions.h>
#include "device_launch_parameters.h"
#define __blocksize 8
#define __chunk 8
__global__ void cholesky4(double* A, double* L, int n) {
	int _b = blockIdx.x;
	int _t = threadIdx.x;
	int _j = _b * __blocksize + _t;
	int _s = _j* __chunk;
	int _e = _s + __chunk;
	if (_s >= n)return;
	if (_e > n)_e = n;



	int _b2 = blockIdx.y;
	int _t2 = threadIdx.y;
	int _j2 = _b2 * __blocksize + _t2;
	int _s2 = _j2 * __chunk;
	int _e2 = _s2 + __chunk;
	if (_s2 >= n)return;
	if (_s2 >= _e2)return;
	for (int j = _s; j < _e; j++) {

		double s = 0;
		double* _ptr = &L[j * n + 0];
		for (int k = 0; k < j; k++) {
			s += *_ptr * *_ptr;
			_ptr++;
		}
		L[j * n + j] = sqrt(A[j * n + j] - s);
		int __s2 = 0;
		int __e2 = 0;
		if (_s2 >= j + 1)__s2 = _s2; else __s2 = j + 1;
		if (_e2 <= n)__e2 = _e2; else __e2 = n;
		if (__s2 >= __e2)return;
			for (int i = __s2; i < __e2; i++) {
				double s = 0;
				double* ptr = &L[i * n + 0];
				double* ptr2 = &L[j * n + 0];
				for (int k = 0; k < j; k++) {
					//s += L[i * n + k] * L[j * n + k];
					s += *ptr * *ptr2;
					ptr++;
					ptr2++;
				}
				L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
			}

	}
}


void kernel(double* A,double *work, int N,cudaStream_t stream) {
	
	dim3 threads(__blocksize, __blocksize);
	int ff = N / __blocksize / __chunk;
	if (ff == 0)ff = 1;
	dim3 grid(ff,ff);

	cholesky4<<<grid, threads,0,stream>>>(A, work, N);

}