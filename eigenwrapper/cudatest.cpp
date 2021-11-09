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
using namespace std;
int main() {
    int n = 1;
    int s = 1;
    int M = 2000;
    int N = 2000;
    
    for (int t = 0; t < 2; t++)
    {
        kingghidorah::cuda cuda(M);
        std::cout << kingghidorah::_mySparse::_testopenmp() << std::endl;

        kingghidorah::_mySparse mat;
        mat.resize(M, N);
        //kingghidorah::_mySparse ret;
        mat._resize(M, N);
        for (int i = 0; i < M; i++)
        {
            mat.adddat(i, i, i + 5+i%4);
        }
        mat.merge();
        mat.clearcoeff();
        double rhs[2000];
        for (int i = 0; i < N; i++)rhs[i] = 1;

        auto start = high_resolution_clock::now();

        int NT = mat.ofAtA(&mat);

        kingghidorah::_mySparse mat2;
        kingghidorah::_mySparse mat3;
        mat3.resize(M, M);
        mat2._OfDuplicate(&mat);
        
        std::cout << cuda.valid() << std::endl;
        std::cout << "count:" << cuda.count() << std::endl;
        std::cout << "fastest" << cuda.fastest() << std::endl;
        std::cout << cuda.device_name() << std::endl;
        auto now = std::chrono::high_resolution_clock::now();
        std::vector<double> x(N);
        for (int i = 0; i < N; i++)x[i] = 1;
        std::cout << "DN" << std::endl;
        kingghidorah::_mySparse ret;
        ret.resize(N, N);
        mat._solveI_gpu(&cuda, &ret);
        for (int i = 0; i < 2; i++)std::cout << ret._at(i, i) << std::endl;
        for (int i = 0; i < 2; i++)
        {
            mat._solveI_gpu_omp(&cuda, &ret);
            std::cout << "MG" << std::endl;
            for (int i = 0; i < 2; i++)std::cout << ret._at(i, i) << std::endl;
        }
        
        //cuda.dispose();
        std::cin.get();
    }

}