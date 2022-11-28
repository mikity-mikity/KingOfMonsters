/*#include <assert.h>
#include <stdio.h>

// CUDA runtime
#include <cuda_runtime.h>

// Helper functions and utilities to work with CUDA
#include <helper_cuda.h>
#include <helper_functions.h>
#include "device_launch_parameters.h"
#define __blocksize 8
#define __chunk 128

#define __blocksize2 4
#define __chunk2 4

__global__ void __add(double* value, int* row, int* col, int N, int M, double* value2, int* index) {
	int _b = blockIdx.x;
	int _t = threadIdx.x;
	int _j = _b * __blocksize + _t;
	int _s = _j * __chunk;
	int _e = _s + __chunk;
	if (_s >= _e)return;
	if (_e > N)_e = N;
	int* _row = row + _s;
	for (int i = _s; i < _e; i++)
	{
		int S = *_row;// row[i];
		_row++;
		int E =* _row;// row[i + 1];
		int __row = i;
		int* _col = col + S;
		double* _val = value + S;
		for (int k = S; k < E; k++)
		{
			int __col = *_col;// col[k];
			double __val = *_val;;// value[k];
			int __index = index[__row * M + __col];
			value2[__index] += __val;
			_val++;
			_col++;
		}		
	}
}
__global__ void __vecmul(double* value, int* row, int* col, int N, int M, double* value2) {
	int _b = blockIdx.x;
	int _t = threadIdx.x;
	int _j = _b * __blocksize + _t;
	int _s = _j * __chunk;
	int _e = _s + __chunk;
	if (_s >= _e)return;
	if (_e > N)_e = N;
	int* _row = row + _s;
	double* _val2 = &value2[_s];
	for (int i = _s; i < _e; i++)
	{
		double val2 = sqrt(*_val2);
		int S = *_row;// row[i];
		_row++;
		int E = *_row;// row[i + 1];
		int __row = i;
		int* _col = col + S;
		double* _val = value + S;
		for (int k = S; k < E; k++)
		{
			int __col = *_col;// col[k];
			(*_val) *= val2;
			_val++;
			_col++;
		}
		_val2++;
	}
}
void kernel(double* value,int* row,int* col,int N,int M,double* value2,int* index,cudaStream_t stream) {
	
	dim3 threads(__blocksize);
	int ff = N / __blocksize / __chunk;
	ff++;
	dim3 grid(ff);
	cudaStreamSynchronize(stream);
	__add <<<grid, threads, 0, stream >> > (value, row, col, N, M, value2,index);

}

void kernel2(double* value, int* row, int* col, int N, int M, double* value2, cudaStream_t stream) {

	dim3 threads(__blocksize2);
	int ff = N / __blocksize2 / __chunk2;
	ff++;
	dim3 grid(ff);
	cudaStreamSynchronize(stream);
	__vecmul<< <grid, threads, 0, stream >> > (value, row, col, N, M, value2);

}
*/