
//#include "cuda_runtime.h"
//#include "device_launch_paraMeters.h"

#include<iostream>
#include <fstream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include <omp.h>
#include<cuda.h>
#include <cuda_runtime_api.h>
#include <chrono>
#include <vector>
using namespace std::chrono;

using namespace std;
int main() {
    int N = 2000;

    double* m;
    cuInit(0);
    cudaSetDevice(0);
    cudaMallocHost((void**)&m, sizeof(double) * N * N);
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            m[i * N + j] = 1;
        }
    }
#pragma omp parallel for
    for (int tt = 0; tt < 2; tt++)
    {
        cudaSetDevice(tt);
        double* m_gpu;
        cudaMalloc((void**)&m_gpu, sizeof(double) * N * N);
        auto start = high_resolution_clock::now();
        cudaMemcpy(m_gpu, m, sizeof(double) * N * N, cudaMemcpyHostToDevice);
        auto end = high_resolution_clock::now();
        auto duration = end - start;
        std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
        std::cout << d.count() << "ms" << std::endl;
        cudaFreeHost(m);
        cudaFree(m_gpu);
    }
    std::cin.get();

}