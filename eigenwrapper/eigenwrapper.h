#pragma once

//#define EIGEN_DONT_PARALLELIZE
#ifdef _CPU
#define EIGEN_USE_MKL_ALL
#define EIGEN_USE_LAPACK
#include "eigen-3.4.0/Eigen/PardisoSupport"
#endif
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
#ifndef _CPU
#include <cuda_runtime.h>
#include <device_launch_paraMeters.h>
#include<cuda.h>
#include <cuda_runtime_api.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <cusolverDn.h>
#include <cusolverMg.h>
#include <cublas_v2.h>
#include <chrono>
#include <vector>
#include <map>
#include <string>
//#define EIGEN_DONT_PARALLELIZE
//#define EIGEN_MALLOC_ALREADY_ALIGNED  0
//#define EIGEN_DONT_ALIGN

#define MAXDEVICE 4
using namespace std::chrono;
using std::vector;
using std::string;
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int64_t;
//void kernel(double* A, double* work, int64_t N, cudaStream_t stream);
//void kernel(double* value, int64_t* row, int64_t* col, int64_t N, int64_t M, double* value2, int64_t* index, cudaStream_t stream);
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
		//std::vector<std::vector<cusolverSpHandle_t>> solver_handleSp;
		//std::vector<std::vector<cusparseHandle_t>> cusparse_handle;

		cublasHandle_t cublas_handle[MAXDEVICE];
		double* __mgM[MAXDEVICE];
		double* __mgM2[MAXDEVICE];
		double* __mgrhs[MAXDEVICE];
		//double* __mgM2 = 0;
		//double* __mgrhs2 = 0;
		double* __mgC[MAXDEVICE];
		double* __work[MAXDEVICE];
		int64_t work_size[MAXDEVICE];
		double* __work2[MAXDEVICE];
		int64_t work_size2[MAXDEVICE];
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

		//int64_t prevT_A = 0;
		//int64_t prevN = 0;
		//int64_t prevwn = 0;
		static void disable();
		cuda(int64_t N);
		~cuda();
		cusolverDnHandle_t& solver(int64_t ii, int64_t kk);
		//cusolverSpHandle_t& solverSp(int64_t ii, int64_t kk);
		cublasHandle_t& blas(int64_t ii);
		//cusolverMgHandle_t mgsolver();
		//double* L();
		bool valid();
		bool canpeer();
		std::string device_name();
#ifndef _CPU
		int* info(int i);

		double* work_M(int i);
		double* work_rhs(int i);
		double* work_C(int i);
		double* work(int64_t N, int i);
		double* work2(int64_t N, int i);
		double* work(int64_t N, int i, cudaStream_t stream);
#endif
		int& count();
		int& fastest();
		void dispose();
		cudaStream_t& __streams(int64_t i, int64_t j);
		//bool canpeeraccess(int64_t i, int64_t j);
		//double** array_d_A();
		//double** array_d_B();
		//double** array_d_work();
		int* devicelist();
	};

	class _myPermutation {
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t> perm;
		_myPermutation(int64_t* ptr, int64_t N);
	};
	class _myDoubleArray
	{
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		Eigen::VectorXd __v;
		void plus_useindex(double* ptr, double sc, int64_t N, int64_t* index);
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
		int64_t* dA_csrOffsets = 0, * dA_columns = 0, * dB_csrOffsets = 0, * dB_columns = 0,
			* dC_csrOffsets = 0, * dC_columns = 0, * dD_csrOffsets = 0, * dD_columns = 0;
		double* dA_values = 0, * dB_values = 0, * dC_values = 0, * dD_values = 0;
		int64_t* index = 0;

		int64_t A_num_rows;
		int64_t A_num_cols;
		int64_t A_nnz;

		int64_t B_num_rows;
		int64_t B_num_cols;
		int64_t B_nnz;
		//cusparseSpMatDescr_t matA, matB, matC;
		//cusparseSpGEMMDescr_t spgemmDesc;
		int64_t C_num_rows, C_num_cols, C_nnz;
		int64_t D_num_rows, D_num_cols, D_nnz;

		void* dBuffer1 = NULL, * dBuffer2 = NULL;
		size_t bufferSize1 = 0, bufferSize2 = 0;
	};
	class _myMicroMatrix {
	public:
		Eigen::MatrixXd _mat;
		Eigen::VectorXd _z, _w, _tmp;
		std::vector<int> indices;
		double _sc;
		double sqrsc;
	};
	class _mySparseVector
	{
	public:
		Eigen::SparseVector<double> _vec;
	};
	class _mySparse {

	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	public:
		vector<vector<double>> _coeff;
		std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> _mat;
		Eigen::MatrixXd _dmat;
		Eigen::MatrixXd _prevmat;

		vector<Eigen::VectorXd> coeff;
	private:
		int64_t space = 0;
	public:
		int64_t _nt = 0;
		int64_t _mt = 0;
		int64_t _dat_count = 0;
		//bool e_init = false;
		//Eigen::SparseQR< Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int64_t>> qr;
		//Eigen::SparseLU< Eigen::SparseMatrix<double>> lu;

		//Eigen::HouseholderQR<Eigen::MatrixXd> qr2;
		//Eigen::PartialPivLU<Eigen::MatrixXd> lu2;
	public:
		vector<vector<Eigen::Triplet<double>>> dat;
		//vector<Eigen::Triplet<double>> dat2;
	public:
		void makePattern();
		double mulboth(Eigen::VectorXd& u, Eigen::VectorXd& v);
		void mulleft(Eigen::VectorXd& u, Eigen::VectorXd& ret, double sc);
		void mulright(Eigen::VectorXd& u, Eigen::VectorXd& ret, double sc);
		void setzero(int row);
		static double computeeigen(_mySparse* i1, _mySparse* i2, _mySparse* f1, int64_t N);
		static void project(_mySparse* i1, _mySparse* i2, _mySparse* f1, _myDoubleArray* v1, _myDoubleArray* v2, _myDoubleArray* _ret1, _myDoubleArray* _ret2);
		_mySparse();
		~_mySparse();
		void freeze(bool _do);
		void freeze2();
		void _freeze();
		double L2Norm(Eigen::VectorXd* a, Eigen::VectorXd* b);
		Eigen::VectorXd Vector(double* ptr1, int64_t N1);
		Eigen::VectorXd Vector(Eigen::VectorXd* a);
		void Vector(Eigen::VectorXd* a, Eigen::VectorXd* ret);
		void plus(_mySparse* m, double sc, bool dense, bool sparse);
		double at(int64_t i, int64_t ii);
		double _at(int64_t i);
		double _at(int64_t i, int64_t j);
		int64_t num_elem(int64_t j);
		int64_t cols();
		//void join();
		std::string info();
		void permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm);
		void shrink(int64_t M);
		void _permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm, bool sparse, bool dense);
		void _permuteCols(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm, bool sparse, bool dense);
		void _permuteRows(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm, bool sparse, bool dense);
		void __permuteCols(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm);
		void __permuteRows(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm);
		void _shrink(int64_t M, bool sparse, bool dense);
		void _shrinkCols(int64_t M, bool sparse, bool dense);
		void _permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm2);
		void _shrink(int64_t M, int64_t N, bool sparse, bool dense);
		Eigen::VectorXd get_coeff(int64_t ii);
		double __at(int64_t i, int64_t j);
		void scale(int i, double sc);
		void scale(double sc);
		void _scale(double sc);
		int64_t rows();
		int64_t _rows();
		int64_t _cols();
		int64_t __rows();
		void copycoefffrom(KingOfMonsters::_mySparse* mat);
		void init(int64_t n, int64_t m);
		int64_t resize(int64_t n, int64_t m);
		void _resize(int64_t n, int64_t m);
		void reserve(int64_t n);
		void addemptyrow(int64_t ii);
		void addrow(int64_t ii, int64_t* ptr, double* data, double sc, int64_t N);
		void addrow(int64_t ii, int64_t* ptr, double* data, double sc, int64_t N, double _coeff);
		void addrow(int64_t ii, int64_t* ptr, double* data, double* data2, double sc, int64_t N, double c1, double c2);
		void addrow(int64_t ii, int64_t* ptr, double* data, int64_t shift, double sc, int64_t N, bool add, double __coeff);
		void adddat(int64_t ii, int64_t j, double value);
		void addcoeff(double sc);
		void addmat(_mySparse* mat);
		void OfDuplicate(_mySparse* mat);
		void _OfDuplicate(_mySparse* mat);
		void ofDat();
		void freezecoeff();
		std::string ofAtA(_mySparse* A, bool sparse);
		
		void add_usemap(int64_t i, int64_t j, double val);
		void add(int64_t i, int64_t j, double val);
		void set_usemap(int64_t i, int64_t j, double val);
		int find_location(int64_t i, int64_t j);
		void add_uselocation(int64_t location, double val);
		void set_uselocation(int64_t location, double val);
		std::string ofAtA_gpu(cuda* _cuda, _mySparse* A, bool sparse);
		std::string _ofAtA(_mySparse* A);
		std::string _ofAtA_sparse(_mySparse* A);
		//void ofAtB_gpu(_mySparse* B, bool sparse);
		void ofAtB(_mySparse* B, bool sparse,bool AorB);
		void _ofAtB(_mySparse* B, _mySparse* C);
		void _ofBtAB(_mySparse* B, /*Eigen::VectorXd* b, */_mySparse* C/*, Eigen::VectorXd* ret*/);
		void _ofCtAB(_mySparse* B, _mySparse* C, /*Eigen::VectorXd* b, */_mySparse* D/*, Eigen::VectorXd* ret*/);
		void _ofBtAB(_mySparse* B, _mySparse* B2, /*Eigen::VectorXd* b, */_mySparse* C/*, Eigen::VectorXd* ret*/);
		//void _ofBtAB2(_mySparse* B, _mySparse* C, _mySparse* Q, _mySparse* R, KingOfMonsters::cuda* cuda);
		void _ofCBtAB(_mySparse* B, _mySparse* C, _mySparse* D);
		void _ofCBtAB2(_mySparse* B, _mySparse* C, _mySparse* D, _mySparse* E);
		//void _ofBtAB_qr(_mySparse* B, Eigen::VectorXd* b, _mySparse* C, Eigen::VectorXd* ret);
		Eigen::VectorXd Atb(double* ptr, int64_t N);
		Eigen::VectorXd _Atb(double* ptr, int64_t N);
		void Atb(double* ptr, int64_t N, Eigen::VectorXd* c);
		void Atb(double* ptr, double* ptr2, double sc, int64_t N, Eigen::VectorXd* c);
		void merge();
		void computeQR();
		void computeLU();
		void computeLLT(Eigen::LLT<Eigen::MatrixXd>* LLT);
		int64_t nonzeros();
		void Clear();
		void setmat(Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>& mat, int64_t ii);
		void setmat(const Eigen::MatrixXd& mat);
		void setmiddlecolum(Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>& f, int64_t start, int64_t end);
		void solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		void LSsolve(Eigen::VectorXd* rhs, Eigen::VectorXd* ret,double,int mode);
		void Project(Eigen::VectorXd* rhs, Eigen::VectorXd* ret, double);
		//void _solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		//Eigen::MatrixXd _solve0(_myLLT* LLT, _mySparse* mat);
		std::string _solve0_lu_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int ordering);
		std::string _solve0_chol_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int ordering);
		void solve0_lu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		void _solve0_lu_cg(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		std::string _solve0_gpu(KingOfMonsters::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int64_t device);
		//std::string _QR_gpu(KingOfMonsters::cuda* cuda, Eigen::MatrixXd* Q, Eigen::MatrixXd* R, int64_t device);
		std::string _solveLU_gpu(KingOfMonsters::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int64_t device);
		std::string _solveLU_sparse_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		std::string _solveCG_sparse_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		std::string _solveLU_dense_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		void turnDense();
		void _solve0_gpu(KingOfMonsters::cuda* cuda, _mySparse* rhs, _mySparse* ret);
		int64_t _solveI(_mySparse* ret);
		std::string _solveI_dense(_mySparse* ret);
		std::string _solveI_gpu_sparse(KingOfMonsters::cuda* cuda, _mySparse* ret);
		std::string _solveI_gpu(KingOfMonsters::cuda* cuda, _mySparse* ret);
		std::string _solveI_cpu(_mySparse* ret);
		std::string _solveI_gpu_omp(KingOfMonsters::cuda* cuda, _mySparse* ret);
		//std::string AinvBA(KingOfMonsters::cuda* cuda, _mySparse* A, _mySparse* ret);
		std::string _solveI_gpu_single(KingOfMonsters::cuda* cuda, _mySparse* ret);
		void plus(Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>* m);
		//void _solveI_gpu_mg(KingOfMonsters::cuda* cuda, _mySparse* ret);
		void __solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret);
		
		Eigen::MatrixXd inv();
		Eigen::MatrixXd solve0(_mySparse* rhs);
		void minus(_mySparse* m);
		void clearcoeff();
		void addsmallidentity(double salt, bool sparse, bool dense, int m);
		void addsmallidentity(double salt, bool sparse, bool dense);
		void begin_construct();
		void end_construct(int64_t c);
		void end_construct2();
		int64_t numBlocks();
		void _plus(int64_t i, int64_t j, double val);
		static std::string _testopenmp();
		//Eigen::SparseMatrix<double>* e = 0;
		//Eigen::SparseMatrix<double>* e2 = 0;

	};
	class _helper {
	public:
		static double VarPro(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt);
		static double ALT(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt);
		static double Simple(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt);
		static double GN(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt, KingOfMonsters::cuda* cuda);
		static void write(Eigen::VectorXd* coeff, Eigen::VectorXd* phi0, Eigen::VectorXd* zz0, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt);
	};

}