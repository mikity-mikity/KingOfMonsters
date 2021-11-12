#include "pch.h"
//#include "cuda_runtime.h"
//#include "device_launch_paraMeters.h"

#include<iostream>
#include <fstream>
#include<iomanip>
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include <omp.h>
/*#include <cusolverDn.h>
#include <cublas_v2.h>
#include <cuda_runtime_api.h>*/
/*
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/SparseQR"
#include "eigen-3.4.0/Eigen/SparseLU"
#include "eigen-3.4.0/Eigen/SparseCholesky"
#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT
*/
#include <chrono>
#include <vector>
using namespace std::chrono;
#include"eigenwrapper.h"
#include "eigen-3.4.0/Eigen/Dense"
using namespace std;
int main() {
    {

        int N = 2;
        double* c;
        double* m;
        double* e;
        cuInit(0);
        cudaSetDevice(0);
        c = (double*)malloc(sizeof(double) * N * N);
        cudaMallocHost((void**)&m, sizeof(double) * N * N);
        cudaMallocHost((void**)&e, sizeof(double) * N * N);
        Eigen::MatrixXd f(N, N);
        m[0] = 1;
        m[1] = 1;
        m[2] = 1;
        m[3] = 3;


        memcpy(f.data(), m, sizeof(double) * N * N);
        std::cout << f << std::endl;
//#pragma omp parallel for
        cudaStream_t stream;
        cudaStreamCreate(&stream);
        cudaSetDevice(0);
        double* m_gpu;
        cudaMallocAsync((void**)&m_gpu, sizeof(double) * N * N,stream);

        cudaMemcpyAsync(m_gpu, m, sizeof(double) * N * N, cudaMemcpyHostToDevice,stream);
        double* work;
        int work_size;
        int work_size1=0;
        int work_size2=0;
        cusolverDnHandle_t solver;
        cusolverDnCreate(&solver);
        cusolverDnSetStream(solver,stream);
        cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size1);
        cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size2);
        work_size = std::max(work_size1, work_size2);
        cudaMallocAsync(&work, sizeof(double)*work_size,stream);
        cudaMemsetAsync(work, 0, sizeof(double) * work_size1,stream);
        int *devInfo;
        cudaMallocAsync(&devInfo, sizeof(int),stream);
        cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work, work_size1, devInfo);
        cudaMemsetAsync(work, 0, sizeof(double) * work_size2,stream);
        cusolverDnDpotri(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work, work_size2, devInfo);


        cudaMemcpyAsync(e, m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost,stream);
        cudaDeviceSynchronize();

        cudaFreeAsync(work,stream);
        cudaFreeAsync(m_gpu, stream);
        cudaFreeAsync(devInfo, stream);

        memcpy(f.data(), e, sizeof(double) * N * N);
        std::cout << f << std::endl;

        cudaFreeHost(m);
        cudaFreeHost(e);
        free(c);
        std::cin.get();


    }

    /*

    int n = 1;
    int s = 1;
    int M = 2;
    int N = 2;

    for (int t = 0; t < 2; t++)
    {
        kingghidorah::cuda cuda(M);
        //std::cout << kingghidorah::_mySparse::_testopenmp() << std::endl;

        kingghidorah::_mySparse mat;
        mat.resize(M, N);
        //kingghidorah::_mySparse ret;
        mat._resize(M, N);
        mat.adddat(0, 0, 2);
        mat.adddat(0, 1, 1);
        mat.adddat(1, 0, 1);
        mat.adddat(1, 1, 4);
        mat.merge();
        mat.clearcoeff();
        double rhs[2000];
        for (int i = 0; i < N; i++)rhs[i] = 1;
        Eigen::MatrixXd f(N, N);
        f.setIdentity();
        Eigen::MatrixXd g(N, N);
        
        int NT = mat.ofAtA(&mat,true);

        kingghidorah::_mySparse mat2;
        kingghidorah::_mySparse mat3;
        mat3.resize(M, M);
        mat2._OfDuplicate(&mat);

        std::cout << cuda.valid() << std::endl;
        std::cout << "count:" << cuda.count() << std::endl;
        std::cout << "fastest" << cuda.fastest() << std::endl;
        std::cout << cuda.device_name() << std::endl;
        std::vector<double> x(N);
        for (int i = 0; i < N; i++)x[i] = 1;
        std::cout << "DN" << std::endl;
        kingghidorah::_mySparse ret;
        ret.resize(N, N);
        auto start = high_resolution_clock::now();
        mat._solveI_gpu(&cuda,&ret);
        auto end = high_resolution_clock::now();
        auto duration = end - start;
        auto d = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
        std::cout << d.count() << "ms" << std::endl;
        std::cout << ret._at(0, 0) << std::endl;
        std::cout << ret._at(0, 1) << std::endl;
        std::cout << ret._at(1, 0) << std::endl;
        std::cout << ret._at(1, 1) << std::endl;
        for (int i = 0; i < 2; i++)
        {
            start = high_resolution_clock::now();
            mat._solveI_gpu_omp(&cuda, &ret);
            std::cout << "MG" << std::endl;
            end = high_resolution_clock::now();
            duration = end - start;
            std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
            std::cout << d.count() << "ms" << std::endl;
            std::cout << ret._at(0, 0) << std::endl;
            std::cout << ret._at(0, 1) << std::endl;
            std::cout << ret._at(1, 0) << std::endl;
            std::cout << ret._at(1, 1) << std::endl;
        }

        //cuda.dispose();
        std::cin.get();
    }*/
    //std::cin.get();

}