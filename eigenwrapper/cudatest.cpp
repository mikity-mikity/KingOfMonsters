#include "pch.h"
#ifndef _CPU
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
using namespace KingOfMonsters;
#define _N 4
#define _M 4
int main() {
    {
        auto _cuda = new KingOfMonsters::cuda(_N);
        
        _mySparse m[_M];
        _mySparse M;
        for (int i = 0; i < _M; i++)
        {
            m[i].init(_N, _N);
        }
        M.init(_N, _N);
        M.begin_construct();
        int index[_N];
        for (int i = 0; i < _N; i++)index[i] = i;
        for (int tt = 0; tt < 1; tt++)
        {
            for (int kk = 0; kk <1; kk++)
            {
                //double f[_N];
                M.Clear();
                for (int ii = 0; ii < _M; ii++)
                {
                    m[ii].Clear();
                    for (int i = 0; i < _N; i++)
                    {
                        //memset(f, 0, sizeof(double) * _N);
                        for (int j = 0;j<_N;j++)
                        {
                            if (j == i)
                            {
                                m[ii].adddat(i, j, 1.0);
                            }
                            else
                                if(j>=0&&j<_N)
                                {
                                    m[ii].adddat(i, j, 0.1);
                                }
                        }
                        m[ii].addcoeff(1.0);
                    }
                    M.addmat(&m[ii]);
                }
                M.end_construct(_N);
                M.ofDat();
                M.clearcoeff();
                auto start = high_resolution_clock::now();
                std::cout << M._mat[0] << std::endl;
                std::cout << M._mat[1] << std::endl;
                std::cout << M._mat[2] << std::endl;
                std::cout << "start" << std::endl;
                //M.ofAtA(&M, true);
                std::cout << M._mat[0] << std::endl;
                M.ofAtA_gpu(_cuda , &M, true);
                std::cout << "end" << std::endl;
                std::cout << M._mat[0] << std::endl;
                //m._shrink(_N, false, true);
                //Eigen::VectorXd ret(_N);
                //m._solve0_gpu(&m, &rhs, &ret, ii);
                auto stop = high_resolution_clock::now();
                auto duration = duration_cast<microseconds>(stop - start);
                if(kk==0)std::cout << duration.count() << "microseconds"<<std::endl;
                
                //speed[ii] = duration.count();
            }
        }

        /*int N = 2000;
        double* c;
        double* m;
        double* e;
        cuInit(0);
        cudaSetDevice(0);
        c = (double*)malloc(sizeof(double) * N * N);
        cudaMallocHost((void**)&m, sizeof(double) * N * N);
        cudaMallocHost((void**)&e, sizeof(double) * N * N);
        Eigen::MatrixXd f(N, N);

        for (int i = 0; i < N; i++)
        {
            for (int j = i - 2; j < i + 3; j++)
            {
                if (j >= 0 && j < N)m[i * N + j] =1;
            }
            m[i * N + i] = 3;
        }
        memcpy(f.data(), m, sizeof(double) * N * N);
        std::cout << f.topLeftCorner(4,4) << std::endl << std::endl;
        int _count = 0;
        cudaGetDeviceCount(&_count);
        double** gpu_matrix = new double* [_count];
        double** gpu_work = new double* [_count];
        cudaStream_t* streams = new cudaStream_t[_count];
        cusolverDnHandle_t* solver=new cusolverDnHandle_t[_count];
        for (int i = 0; i < _count; i++)
        {
            cudaSetDevice(i);
            cusolverDnCreate(&solver[i]);
            cusolverDnSetStream(solver[i], streams[i]);

        }
#pragma omp parallel for
        for (int ii = 0; ii < _count; ii++)
        {
            cudaSetDevice(ii);
            cudaStreamCreate(&streams[ii]);
            cudaMallocAsync((void**)&gpu_matrix[ii], sizeof(double) * N * N, streams[ii]);
            cudaMemcpyAsync(gpu_matrix[ii], m, sizeof(double) * N * N, cudaMemcpyHostToDevice, streams[ii]);
            cudaMallocAsync((void**)&gpu_work[ii], sizeof(double) * N * N, streams[ii]);
        }

        for (int ii = 0; ii <_count; ii++)
        {
            auto _start = std::chrono::high_resolution_clock::now();
            cudaSetDevice(ii);
            cudaStream_t stream = streams[ii];
            double* m_gpu = gpu_matrix[ii];

            cudaMemcpyAsync(m_gpu, m, sizeof(double) * N * N, cudaMemcpyHostToDevice, stream);
            double* work = gpu_work[ii];
            int work_size1 = 0;
            int work_size2 = 0;
            cusolverDnDpotrf_bufferSize(solver[ii], CUBLAS_FILL_MODE_UPPER, N, m_gpu, N, &work_size1);
            cudaMemsetAsync(work, 0, sizeof(double) * N*N, stream);
            int* devInfo;
            cudaMallocAsync(&devInfo, sizeof(int), stream);
            cusolverDnSetStream(solver[ii], stream);
            if (ii == 0)
            {
                cusolverDnDpotrf(solver[ii], CUBLAS_FILL_MODE_UPPER, N, m_gpu, N, work, work_size1, devInfo);
                cudaMemcpyAsync(e, m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost, stream);
            }
            else {
                kernel(m_gpu, work, N, stream);
                cudaMemcpyAsync(e, work, sizeof(double) * N * N, cudaMemcpyDeviceToHost, stream);
            }
            cudaDeviceSynchronize();
            cudaFreeAsync(devInfo, stream);

            memcpy(f.data(), e, sizeof(double) * N * N);
            std::cout << f.topLeftCorner(4, 4) << std::endl << std::endl;
            std::cout << f.bottomRightCorner(4, 4) << std::endl << std::endl;

            auto _end = std::chrono::high_resolution_clock::now();
            std::cout << "potrf"<<std::chrono::duration_cast<std::chrono::milliseconds>(_end - _start).count() << "ms" << std::endl;
        }
        */
        /*
        auto start = std::chrono::high_resolution_clock::now();
#pragma omp parallel for
        for (int ii = 0; ii < _count; ii++)
        {
            auto _start = std::chrono::high_resolution_clock::now();
            cudaSetDevice(ii);
            cudaStream_t stream=streams[ii];
            double* m_gpu=gpu_matrix[ii];

            cudaMemcpyAsync(m_gpu, e, sizeof(double) * N * N, cudaMemcpyHostToDevice, stream);
            double* work=gpu_work[ii];
            int work_size2 = 0;
            cusolverDnDpotri_bufferSize(solver[ii], CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size2);
            int* devInfo;
            cudaMallocAsync(&devInfo, sizeof(int), stream);
            cudaMemsetAsync(work, 0, sizeof(double) * work_size2, stream);
            
            cusolverDnDpotri(solver[ii], CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work, work_size2, devInfo);
            cudaMemcpyAsync(e, m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost, stream);
            cudaDeviceSynchronize();

            cudaFreeAsync(work, stream);
            cudaFreeAsync(m_gpu, stream);
            cudaFreeAsync(devInfo, stream);

            memcpy(f.data(), e, sizeof(double) * N * N);
            auto _end = std::chrono::high_resolution_clock::now();
            std::cout << f.topLeftCorner(4,4)<< std::endl << std::endl;


            std::cout << "potri"<<std::chrono::duration_cast<std::chrono::milliseconds>(_end - _start).count() << "ms" << std::endl;
        }
        auto end = std::chrono::high_resolution_clock::now();

        std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
        */
        /*for (int i = 0; i < _count; i++)
        {
            cudaStreamDestroy(streams[i]);
            cusolverDnDestroy(solver[i]);
        }
        cudaFreeHost(m);
        cudaFreeHost(e);
        free(c);
        delete[] gpu_matrix;
        delete[] gpu_work;

        for (int i = 0; i < _count; i++)
        {
            cudaSetDevice(i);
            cudaDeviceReset();
        }*/
        std::cin.get();


    }

    /*

    int n = 1;
    int s = 1;
    int M = 2;
    int N = 2;

    for (int t = 0; t < 2; t++)
    {
        KingOfMonsters::cuda cuda(M);
        //std::cout << KingOfMonsters::_mySparse::_testopenmp() << std::endl;

        KingOfMonsters::_mySparse mat;
        mat.resize(M, N);
        //KingOfMonsters::_mySparse ret;
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

        KingOfMonsters::_mySparse mat2;
        KingOfMonsters::_mySparse mat3;
        mat3.resize(M, M);
        mat2._OfDuplicate(&mat);

        std::cout << cuda.valid() << std::endl;
        std::cout << "count:" << cuda.count() << std::endl;
        std::cout << "fastest" << cuda.fastest() << std::endl;
        std::cout << cuda.device_name() << std::endl;
        std::vector<double> x(N);
        for (int i = 0; i < N; i++)x[i] = 1;
        std::cout << "DN" << std::endl;
        KingOfMonsters::_mySparse ret;
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
#endif