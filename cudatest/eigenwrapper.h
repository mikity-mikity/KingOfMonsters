#pragma once

//#define EIGEN_DONT_PARALLELIZE

#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/SparseQR"
#include "eigen-3.4.0/Eigen/SparseLU"
#include "eigen-3.4.0/Eigen/SparseCholesky"	
#include <cuda_runtime.h>
#include <device_launch_paraMeters.h>
#include <stdlib.h>
#include <stdio.h>
#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cublas_v2.h>
#include<cuda.h>
#include <cuda_runtime_api.h>
#include <chrono>
#include <vector>
//#define EIGEN_DONT_PARALLELIZE
//#define EIGEN_DONT_ALIGN

#define MAXDEVICE 4
using namespace std::chrono;
using std::vector;
using std::string;

//#define EIGEN_MALLOC_ALREADY_ALIGNED  0

namespace kingghidorah {
	class cuda {
	private:

		std::vector<std::string> rank;
		bool failed;
		bool initialized;
		int _count = 0;
		int _fastest = 0;
		cusolverDnHandle_t solver_handle[MAXDEVICE];
		cublasHandle_t cublas_handle[MAXDEVICE];
		double* __mgM[MAXDEVICE];
		double* __mgrhs[MAXDEVICE];
		double* __mgM2=0;
		double* __mgrhs2=0;
		double* __mgC[MAXDEVICE];
		double* __work[MAXDEVICE];
		int work_size[MAXDEVICE];
		double* _array_d_A[MAXDEVICE];
		double* _array_d_B[MAXDEVICE];
		double* _array_d_work[MAXDEVICE];
		int* __info[MAXDEVICE];
		int _deviceList[MAXDEVICE];
		//double* _L=0;
		std::vector<int> speed;

		cusolverMgHandle_t mg_solver = 0;


	public:
		int prevT_A = 0;
		int prevN = 0;
		int prevwn = 0;
		cuda(int N);
		~cuda();
		cusolverDnHandle_t& solver(int ii);
		cublasHandle_t& blas(int ii);
		cusolverMgHandle_t mgsolver();
		//double* L();
		bool valid();
		std::string device_name();
		int* info(int i);
		double* work_M(int i);
		double* work_rhs(int i);
		double* work_M2();
		double* work_rhs2();
		double* work_C(int i);
		double* work(int N,int i);
		int& count();
		int& fastest();
		void dispose();
		bool canpeeraccess(int i, int j);
		double** array_d_A();
		double** array_d_B();
		double** array_d_work();
		int* devicelist();
	};
	
	class _myPermutation {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
		_myPermutation(int* ptr, int N);
	};
	class _myLLT {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::LLT<Eigen::MatrixXd>* LLT;
	};
	class _mySparse {

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		vector<vector<double>> _coeff;
	private:
		std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> _mat;
		Eigen::MatrixXd mats;
		Eigen::MatrixXd _dmat;
		
		vector<Eigen::VectorXd> coeff;
		int _nt=0;
		int _mt = 0;
		int _dat_count = 0;
		//bool e_init = false;
		//Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> qr;
		//Eigen::SparseLU< Eigen::SparseMatrix<double>> lu;

		//Eigen::HouseholderQR<Eigen::MatrixXd> qr2;
		//Eigen::PartialPivLU<Eigen::MatrixXd> lu2;
		vector<vector<Eigen::Triplet<double>>> dat;
		//vector<Eigen::Triplet<double>> dat2;
	public:
		_mySparse();
		~_mySparse();
		void freeze();
		double L2Norm(double* ptr1, int N1, double* ptr2, int N2);
		Eigen::VectorXd Vector(double* ptr1, int N1);
		void plus(_mySparse* m, double sc);
		double at(int i, int ii);
		double _at(int i);
		double _at(int i, int j);
		int num_elem(int j);
		int cols();
		std::string info();
		void permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm);
		void shrink(int M);
		void _permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm);
		void _shrink(int M);
		void _permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm2);
		void _shrink(int M,int N);
		Eigen::VectorXd get_coeff(int ii);
		int rows();
		int _rows();
		int _cols();
		int __rows();
		void copycoefffrom(kingghidorah::_mySparse* mat);
		void init(int n, int m);
		int resize(int n, int m);
		void _resize(int n, int m);
		void reserve(int n);
		void addemptyrow(int ii);
		void addrow(int ii, int* ptr, double* data, double sc, int N);
		void addrow(int ii, int* ptr, double* data, int shift, double sc, int N, bool add);
		void adddat(int ii, int j, double value);
		void addcoeff(double sc);
		void addmat(_mySparse* mat);
		void OfDuplicate(_mySparse* mat);
		void _OfDuplicate(_mySparse* mat);
		void ofDat();
		void freezecoeff();
		int ofAtA(_mySparse* A);
		std::string _ofAtA(_mySparse* A);
		void ofAtB(_mySparse* B);
		void _ofAtB(_mySparse* B, _mySparse* C);
		Eigen::VectorXd _ofBtAB(_mySparse* B,double* ptr,int N, _mySparse* C);
		void _ofAtB_gpu(kingghidorah::cuda* cuda,_mySparse* B, _mySparse* C);
		Eigen::VectorXd Atb(double* ptr, int N);
		Eigen::VectorXd _Atb(double* ptr, int N);
		void merge();
		void computeQR();
		void computeLU();
		void computeLLT(Eigen::LLT<Eigen::MatrixXd>* LLT);
		int nonzeros();
		void Clear();
		void setmat(Eigen::SparseMatrix<double> mat,int ii);
		void setmat(const Eigen::MatrixXd &mat);
		void setmiddlecolum(Eigen::SparseMatrix<double> f, int start, int end);
		Eigen::VectorXd solve0(double* rhs, int N);
		Eigen::VectorXd _solve0(double* rhs, int N);
		Eigen::VectorXd _solve0_gpu(kingghidorah::cuda* cuda, double* rhs, int N, int device);
		Eigen::VectorXd _solve0_gpu_mg(kingghidorah::cuda* cuda, double* rhs, int N);
		Eigen::MatrixXd _solve0(_myLLT* LLT, _mySparse* rhs);
		void _solve0_gpu(kingghidorah::cuda* cuda, _mySparse* rhs, _mySparse* ret);
		void _solveI(_mySparse* ret);
		void _solveI_gpu(kingghidorah::cuda* cuda, _mySparse* ret);
		void _solveI_gpu_omp(kingghidorah::cuda* cuda, _mySparse* ret);
		void _solve0_gpu_mg(kingghidorah::cuda* cuda, _mySparse* rhs, _mySparse* ret);
		void _solveI_gpu_mg(kingghidorah::cuda* cuda,_mySparse* ret);
		Eigen::VectorXd __solve0(double* rhs, int N);
		Eigen::MatrixXd inv();
		Eigen::MatrixXd solve0(_mySparse* rhs);
		void minus(_mySparse* m);
		void clearcoeff();
		void addsmallidentity(double salt);
		void begin_construct();
		void end_construct(int c);
		int numBlocks();
		static std::string _testopenmp();

	};
}