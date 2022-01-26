#pragma once

//#define EIGEN_DONT_PARALLELIZE

#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/SparseQR"
#include "eigen-3.4.0/Eigen/SparseLU"
#include "eigen-3.4.0/Eigen/SparseCholesky"	

/*#include "eigen-3.3.8/Eigen/Sparse"
#include "eigen-3.3.8/Eigen/Dense"
#include "eigen-3.3.8/Eigen/SparseQR"
#include "eigen-3.3.8/Eigen/SparseLU"
#include "eigen-3.3.8/Eigen/SparseCholesky"
*/

#include <cuda_runtime.h>
#include <device_launch_paraMeters.h>
#include <stdlib.h>
#include <stdio.h>
#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cusolverSp.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>
#include<cuda.h>
#include <cuda_runtime_api.h>
#include <chrono>
#include <vector>
#include <map>
//#define EIGEN_DONT_PARALLELIZE
//#define EIGEN_MALLOC_ALREADY_ALIGNED  0
//#define EIGEN_DONT_ALIGN

#define MAXDEVICE 4
using namespace std::chrono;
using std::vector;
using std::string;

//void kernel(double* A, double* work, int N, cudaStream_t stream);
void kernel(double* value, int* row, int* col, int N, int M, double* value2, int* index, cudaStream_t stream);
namespace KingOfMonsters {
	class cuda {
	public:

	private:
		bool _canpeer = false;
		std::vector<std::string> rank;
		bool failed;
		bool initialized;
		int _count = 0;
		int _fastest = 0;
		std::vector<std::vector<cusolverDnHandle_t>> solver_handle;
		cublasHandle_t cublas_handle[MAXDEVICE];
		double* __mgM[MAXDEVICE];
		double* __mgrhs[MAXDEVICE];
		//double* __mgM2 = 0;
		//double* __mgrhs2 = 0;
		double* __mgC[MAXDEVICE];
		double* __work[MAXDEVICE];
		int work_size[MAXDEVICE];
		//double* _array_d_A[MAXDEVICE];
		//double* _array_d_B[MAXDEVICE];
		//double* _array_d_work[MAXDEVICE];
		int* __info[MAXDEVICE];
		int _deviceList[MAXDEVICE];
		//double* _L=0;
		std::vector<int> speed;

		//cusolverMgHandle_t mg_solver = 0;

		std::vector< std::vector<cudaStream_t>> _streams;

	public:
		//cudaLibMgMatrixDesc_t descrA;
		//cudaLibMgGrid_t gridA;
		//cusolverMgGridMapping_t mapping = CUDALIBMG_GRID_MAPPING_COL_MAJOR;

		//int prevT_A = 0;
		//int prevN = 0;
		//int prevwn = 0;
		static void disable();
		cuda(int N);
		~cuda();
		cusolverDnHandle_t& solver(int ii, int kk);
		cublasHandle_t& blas(int ii);
		//cusolverMgHandle_t mgsolver();
		//double* L();
		bool valid();
		bool canpeer();
		std::string device_name();
		int* info(int i);
		double* work_M(int i);
		double* work_rhs(int i);
		//double* work_M2();
		//double* work_rhs2();
		double* work_C(int i);
		double* work(int N, int i);
		double* work(int N, int i, cudaStream_t stream);
		int& count();
		int& fastest();
		void dispose();
		cudaStream_t& __streams(int i, int j);
		//bool canpeeraccess(int i, int j);
		//double** array_d_A();
		//double** array_d_B();
		//double** array_d_work();
		int* devicelist();
	};

	class _myPermutation {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm;
		_myPermutation(int* ptr, int N);
	};
	class _myDoubleArray
	{
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::VectorXd __v;

	};
	class _myLLT {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::LLT<Eigen::MatrixXd>* LLT;
	};
	struct spgemm {
	public:
		bool initialized = false;
		int* dA_csrOffsets = 0, * dA_columns = 0, * dB_csrOffsets = 0, * dB_columns = 0,
			* dC_csrOffsets = 0, * dC_columns = 0, * dD_csrOffsets = 0, * dD_columns = 0;
		double* dA_values = 0, * dB_values = 0, * dC_values = 0, * dD_values = 0;
		int* index = 0;

		int A_num_rows;
		int A_num_cols;
		int A_nnz;

		int B_num_rows;
		int B_num_cols;
		int B_nnz;
		//cusparseSpMatDescr_t matA, matB, matC;
		//cusparseSpGEMMDescr_t spgemmDesc;
		int C_num_rows, C_num_cols, C_nnz;
		int D_num_rows, D_num_cols, D_nnz;

		void* dBuffer1 = NULL, * dBuffer2 = NULL;
		size_t bufferSize1 = 0, bufferSize2 = 0;
	};
	class _mySparse {

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		vector<vector<double>> _coeff;
		std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor>> _mat;
		Eigen::MatrixXd _dmat;
	private:
		vector<Eigen::VectorXd> coeff;
		int space = 0;
		int _nt = 0;
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
		static double computeeigen(_mySparse* i1, _mySparse* i2, _mySparse* f1, int N);
		static double computeeigen2(_mySparse* i1, _mySparse* i2, _mySparse* f1, int N);
		static Eigen::VectorXd minilla(_mySparse* i1, _mySparse* i2, _mySparse* i3, _myDoubleArray* grad);
		_mySparse();
		~_mySparse();
		void freeze(bool _do);
		void _freeze();
		double L2Norm(Eigen::VectorXd* a, Eigen::VectorXd* b);
		Eigen::VectorXd Vector(double* ptr1, int N1);
		Eigen::VectorXd Vector(Eigen::VectorXd* a);
		void Vector(Eigen::VectorXd* a, Eigen::VectorXd* ret);
		void plus(_mySparse* m, double sc, bool dense, bool sparse);
		double at(int i, int ii);
		double _at(int i);
		double _at(int i, int j);
		int num_elem(int j);
		int cols();
		void join();
		std::string info();
		void permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm);
		void shrink(int M);
		void _permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm, bool sparse, bool dense);
		void _permuteCols(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm, bool sparse, bool dense);
		void _shrink(int M, bool sparse, bool dense);
		void _permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm2);
		void _shrink(int M, int N);
		Eigen::VectorXd get_coeff(int ii);
		double __at(int i, int j);
		int rows();
		int _rows();
		int _cols();
		int __rows();
		void copycoefffrom(KingOfMonsters::_mySparse* mat);
		void init(int n, int m);
		int resize(int n, int m);
		void _resize(int n, int m);
		void reserve(int n);
		void addemptyrow(int ii);
		void addrow(int ii, int* ptr, double* data, double sc, int N);
		void addrow(int ii, int* ptr, double* data, double sc, int N, double _coeff);
		void addrow(int ii, int* ptr, double* data, int shift, double sc, int N, bool add, double __coeff);
		void adddat(int ii, int j, double value);
		void addcoeff(double sc);
		void addmat(_mySparse* mat);
		void OfDuplicate(_mySparse* mat);
		void _OfDuplicate(_mySparse* mat);
		void ofDat();
		void freezecoeff();
		std::string ofAtA(_mySparse* A, bool sparse);
		std::string ofAtA_gpu(cuda* _cuda, _mySparse* A, bool sparse);
		std::string _ofAtA(_mySparse* A);
		//void ofAtB_gpu(_mySparse* B, bool sparse);
		void ofAtB(_mySparse* B, bool sparse);
		void _ofAtB(_mySparse* B, _mySparse* C);
		void _ofBtAB(_mySparse* B, Eigen::VectorXd* b, _mySparse* C, Eigen::VectorXd* ret);
		void _ofBtAB_qr(_mySparse* B, Eigen::VectorXd* b, _mySparse* C, Eigen::VectorXd* ret);
		Eigen::VectorXd Atb(double* ptr, int N);
		Eigen::VectorXd _Atb(double* ptr, int N);
		void Atb(double* ptr, double* ptr2, double sc, int N, Eigen::VectorXd* c);
		void merge();
		void computeQR();
		void computeLU();
		void computeLLT(Eigen::LLT<Eigen::MatrixXd>* LLT);
		int nonzeros();
		void Clear();
		void setmat(Eigen::SparseMatrix<double, Eigen::ColMajor>& mat, int ii);
		void setmat(const Eigen::MatrixXd& mat);
		void setmiddlecolum(Eigen::SparseMatrix<double, Eigen::ColMajor>& f, int start, int end);
		void solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		void _solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		void _solve0_lu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		std::string _solve0_gpu(KingOfMonsters::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int device);
		std::string _solveLU_gpu(KingOfMonsters::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int device);
		Eigen::MatrixXd _solve0(_myLLT* LLT, _mySparse* rhs);
		void _solve0_gpu(KingOfMonsters::cuda* cuda, _mySparse* rhs, _mySparse* ret);
		int _solveI(_mySparse* ret);
		std::string _solveI_gpu(KingOfMonsters::cuda* cuda, _mySparse* ret);
		std::string _solveI_gpu_omp(KingOfMonsters::cuda* cuda, _mySparse* ret);
		std::string _solveI_gpu_single(KingOfMonsters::cuda* cuda, _mySparse* ret);
		void plus(Eigen::SparseMatrix<double, Eigen::ColMajor>* m);
		//void _solveI_gpu_mg(KingOfMonsters::cuda* cuda, _mySparse* ret);
		void __solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		Eigen::MatrixXd inv();
		Eigen::MatrixXd solve0(_mySparse* rhs);
		void minus(_mySparse* m);
		void clearcoeff();
		void addsmallidentity(double salt, bool sparse, bool dense);
		void begin_construct();
		void end_construct(int c);
		int numBlocks();
		void _plus(int i, int j, double val);
		static std::string _testopenmp();
		//Eigen::SparseMatrix<double>* e = 0;
		//Eigen::SparseMatrix<double>* e2 = 0;

	};
}