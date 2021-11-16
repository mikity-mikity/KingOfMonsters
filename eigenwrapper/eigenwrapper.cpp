// eigenwrapper.cpp : Defines the functions for the static library.
//solveI_gpu
#include "pch.h"
#include "framework.h"
#include<iostream>
#include<stdlib.h>
#include<stdio.h>

#include "eigenwrapper.h"
#include "utill.h"
#include <omp.h> 
#include <iostream>
#include <fstream>
int previdentiyN = 0;
//std::vector<cudaStream_t> streams;
Eigen::MatrixXd I;
#define STRTREAMCOUNT 2
kingghidorah::cuda::cuda(int N) {
	I.resize(0, 0);
	omp_set_dynamic(false);
	omp_set_num_threads(16);
	prevT_A = 0;
	prevN = 0;
	prevwn = 0;
	CUresult res;
	res = cuInit(0);
	previdentiyN = 0;

	auto err = cuDeviceGetCount(&_count);

	if (_count > 0)
	{
		for (int i = 0; i < _count; i++)_deviceList[i] = i;
	}
	if (_count > 1)
	{
		_canpeer=enablePeerAccess(_count, _deviceList);
	}
	else {
		_canpeer = false;
	}
	solver_handle.resize(_count);
	_streams.resize(_count);
	for (int ii = 0; ii < _count; ii++)
	{
		solver_handle[ii].resize(STRTREAMCOUNT);
		_streams[ii].resize(STRTREAMCOUNT);
		cudaSetDevice(ii);
		for(int kk=0;kk< STRTREAMCOUNT;kk++)
		cudaStreamCreate(&_streams[ii][kk]);
	}

	//std::cout << err << "yay" << std::endl;
	for (int i = 0; i < _count; i++)
	{
		for (int j = 0; j < STRTREAMCOUNT; j++)
		{
			solver_handle[i][j] = 0;
		}
	}

	for (int i = 0; i < MAXDEVICE; i++)
	{
		cublas_handle[i] = 0;
	}


	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		cusolverStatus_t status;
		for(int j=0;j<STRTREAMCOUNT;j++)
			status = cusolverDnCreate(&solver_handle[ii][j]);
		auto status2 = cublasCreate(&cublas_handle[ii]);

		if (status == cusolverStatus_t::CUSOLVER_STATUS_SUCCESS)
		{
			initialized = true;
			failed = false;
		}
		else {
			initialized = false;
			failed = true;
			for(int j=0;j<STRTREAMCOUNT;j++)
			solver_handle[ii][j] = 0;
			return;
		}
		if (status2 == cublasStatus_t::CUBLAS_STATUS_SUCCESS)
		{
			initialized = true;
			failed = false;
		}
		else {
			initialized = false;
			failed = true;
			cublas_handle[ii] = 0;
			return;
		}
	}
	if (!initialized || failed)return;


	__mgM2 = 0;
	__mgrhs2 = 0;

	for (int i = 0; i < MAXDEVICE; i++)
	{
		__mgM[i] = 0;
		__mgrhs[i] = 0;
		__mgC[i] = 0;
		__work[i] = 0;
		work_size[i] = 0;
		__info[i] = 0;
	}
	_fastest = 0;

	int _N = 1000;
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		cudaMalloc(&__mgM[ii], sizeof(double) * _N * _N);
		cudaMalloc(&__mgrhs[ii], sizeof(double) * _N * _N);
		cudaMalloc(&__mgC[ii], sizeof(double) * _N * _N);
	}
	cudaMallocHost(&__mgM2, sizeof(double) * _N * _N);
	cudaMallocHost(&__mgrhs2, sizeof(double) * _N * _N);

	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		//cudaStreamCreate(&streams[ii]);
		//cusolverDnSetStream(solver_handle[ii], streams[ii]);
		//cublasSetStream(cublas_handle[ii], streams[ii]);
	}
	speed.resize(count());
	for (int i = 0; i < count(); i++)speed[i] = i;
	for (int kk = 0; kk < 3; kk++)
	{
		for (int ii = 0; ii < count(); ii++)
		{
			auto start = high_resolution_clock::now();
			_mySparse m;
			m.init(_N, _N);
			for (int i = 0; i < _N; i++)
			{
				for (int j = i - 5; j < i + 5; j++)
				{
					if (j >= 0 && j < _N)
						m.adddat(i, j, i);
				}
			}
			double* rhs;
			rhs = new double[_N];
			double* ret = new double[_N];
			for (int i = 0; i < _N; i++)rhs[i] = i;
			m.ofDat();
			m.clearcoeff();
			m._ofAtA(&m);
			m._solve0_gpu(this, rhs, _N, ret,ii);
			delete[] rhs;
			delete[] ret;
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			speed[ii] = duration.count();
		}
	}
	_fastest = std::distance(speed.begin(), std::min_element(speed.begin(), speed.end()));
	//_fastest = 0;
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		cudaFree(__mgM[ii]);
		cudaFree(__mgrhs[ii]);
		cudaFree(__mgC[ii]);
	}
	cudaFreeHost(__mgM2);
	cudaFreeHost(__mgrhs2);
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);		
		cudaMalloc(&__mgM[ii], sizeof(double) * N * N);
		cudaMalloc(&__mgrhs[ii], sizeof(double) * N * N);
		cudaMalloc(&__mgC[ii], sizeof(double) * N * N);
		cudaMalloc(&__info[ii], sizeof(int) * 10);
	}
	cudaMallocHost(&__mgM2, sizeof(double) * N * N);
	cudaMallocHost(&__mgrhs2, sizeof(double) * N * N);
	for (int i = 0; i < MAXDEVICE; i++)
	{
		_array_d_A[i] = 0;// = new double* [count()];
		_array_d_B[i] = 0;
		_array_d_work[i] = 0;
	}



	/*if (_count == 4)
	{
		_count = 2;
		_fastest = 0;
		cusolverStatus_t status = cusolverMgDeviceSelect(
			mg_solver,
			2,
			_deviceList); //seem like cannot be called twice. So call this only once here.
	}
	else */ {
		cusolverMgCreate(&mg_solver);
		cusolverStatus_t status = cusolverMgDeviceSelect(
			mg_solver,
			_count,
			_deviceList); //seem like cannot be called twice. So call this only once here.


		assert(CUSOLVER_STATUS_SUCCESS == status);
		{
			auto status = cusolverMgCreateDeviceGrid(&(gridA), 1, _count, _deviceList, mapping);
			assert(CUSOLVER_STATUS_SUCCESS == status);
		}
		const int IA = 1;
		const int JA = 1;
		const int T_A = N;// / nbGpus;
		const int lda = N;
		int NRHS = N;
		status = cusolverMgCreateMatrixDesc(
			&descrA,
			N,
			N,
			N,
			T_A,
			CUDA_R_64F,
			gridA);
		assert(CUSOLVER_STATUS_SUCCESS == status);
		createMat<double>(
			_count,
			_deviceList,
			N,
			T_A,
			lda,
			array_d_A()
			);
		int64_t lwork_potrf = 0;
		int64_t lwork_potri = 0;

		status = cusolverMgPotrf_bufferSize(
			mg_solver,
			CUBLAS_FILL_MODE_LOWER,
			N,
			(void**)array_d_A(),
			IA,
			JA,
			descrA,
			CUDA_R_64F,
			&lwork_potrf);
		assert(CUSOLVER_STATUS_SUCCESS == status);



		status = cusolverMgPotri_bufferSize(
			mg_solver,
			CUBLAS_FILL_MODE_LOWER,
			N,
			(void**)array_d_A(),
			IA,
			JA,
			descrA,
			CUDA_R_64F,
			&lwork_potri);
		assert(CUSOLVER_STATUS_SUCCESS == status);
		prevT_A = T_A;
		prevN = N;
		int lwork = (lwork_potrf > lwork_potri) ? lwork_potrf : lwork_potri;
		double** _array_d_work = array_d_work();
		//if (prevwn < lwork)
		{
			workspaceAlloc(
				_count,
				_deviceList,
				sizeof(double) * lwork,
				(void**)_array_d_work
			);
			prevwn = lwork;
		}
	}
	//_L = new double[N * N];
	/*if (CUSOLVER_STATUS_SUCCESS != status)
	{
		initialized = false;
		failed = true;
		return;
	}*/
}


int* kingghidorah::cuda::devicelist()
{
	return _deviceList;
}
double** kingghidorah::cuda::array_d_A()
{
	return this->_array_d_A;
}
double** kingghidorah::cuda::array_d_B()
{
	return this->_array_d_B;
}
double** kingghidorah::cuda::array_d_work()
{
	return this->_array_d_work;
}
bool kingghidorah::cuda::canpeeraccess(int i, int j)
{
	int canAccessPeer = 0;
	cudaDeviceCanAccessPeer(&canAccessPeer, i, j);
	return canAccessPeer;
}
/*double* kingghidorah::cuda::L()
{
	return this->_L;
}*/
bool kingghidorah::cuda::canpeer()
{
	return _canpeer;
}
cudaStream_t& kingghidorah::cuda::__streams(int i, int j)
{
	return _streams[i][j];
}
void kingghidorah::cuda::dispose() {
	if (valid())
	{
		if (prevwn != 0)
		{
			workspaceFree(_count, _deviceList, (void**)array_d_work());
			prevwn = 0;
		}
		if (prevN != 0)
		{
			destroyMat(
				_count,
				_deviceList,
				prevN,
				prevT_A,
				(void**)array_d_A());
		}
		for (int ii = 0; ii < _count; ii++)
		{
			cudaSetDevice(ii);
			for(int kk=0;kk< STRTREAMCOUNT;kk++)
			cudaStreamDestroy(_streams[ii][kk]);
		}
		previdentiyN = 0;
		/*if (_L != 0)
			delete[] _L;
		_L = 0;*/
		if (prevN != 0)
			destroyMat(
				_count,
				_deviceList,
				prevN,
				prevT_A,
				(void**)_array_d_A);
		if (prevwn != 0)
			workspaceFree(_count, _deviceList, (void**)_array_d_work);
		prevN = 0;
		prevT_A = 0;
		prevwn = 0;
		for (int i = 0; i < MAXDEVICE; i++)
		{
			_array_d_A[i] = 0;
			_array_d_B[i] = 0;
			_array_d_work[i] = 0;
		}
		if (__mgM2 != 0)
		{
			cudaFree(__mgM2);
		}
		if (__mgrhs2 != 0)
		{
			cudaFree(__mgrhs2);
		}
		__mgM2 = 0;
		__mgrhs2 = 0;
		for (int i = 0; i < count(); i++)
		{
			cudaSetDevice(i);
			if (__mgM[i] != 0)
			{
				cudaFree(__mgM[i]);
			}
			if (__mgrhs[i] != 0)
			{
				cudaFree(__mgrhs[i]);
			}
			if (__mgC[i] != 0)
			{
				cudaFree(__mgC[i]);
			}
			if (__work[i] != 0)
			{
				cudaFree(__work[i]);
			}
			if (__info[i] != 0)
			{
				cudaFree(__info[i]);
			}
			__info[i] = 0;
			__mgM[i] = 0;
			__mgrhs[i] = 0;
			__mgC[i] = 0;
			__work[i] = 0;
			work_size[i] = 0;
			for(int j=0;j<STRTREAMCOUNT;j++)
			if (solver_handle[i][j] != 0)
			{
				cusolverDnDestroy(solver_handle[i][j]);
				solver_handle[i][j] = 0;
			}
			if (cublas_handle != 0)
				cublasDestroy(cublas_handle[i]);
			cublas_handle[i] = 0;
		}
		if (mg_solver != 0)
		{
			cusolverMgDestroy(mg_solver);
		}
		mg_solver = 0;
		cudaDeviceReset();
	}
	initialized = false;
}
cusolverMgHandle_t kingghidorah::cuda::mgsolver() {
	return mg_solver;
}
int* kingghidorah::cuda::info(int i) {
	return __info[i];
}
kingghidorah::cuda::~cuda() {

	dispose();
}
double* kingghidorah::cuda::work_M(int i) {
	return __mgM[i];
}
double* kingghidorah::cuda::work_rhs(int i) {
	return __mgrhs[i];
}
double* kingghidorah::cuda::work_M2() {
	return __mgM2;
}
double* kingghidorah::cuda::work_rhs2() {
	return __mgrhs2;
}
double* kingghidorah::cuda::work_C(int i) {
	return __mgC[i];
}
double* kingghidorah::cuda::work(int N, int device) {
	if (N < work_size[device])
	{
		//do nothing
		return __work[device];
	}
	else {
		if (__work[device] != 0)
		{
			//cudaSetDevice(device);
			cudaFree(__work[device]);
			cudaMalloc(&__work[device], N * sizeof(double));
		}
		work_size[device] = N;
		return __work[device];
	}
}
double* kingghidorah::cuda::work(int N, int device,cudaStream_t stream) {
	if (N < work_size[device])
	{
		//do nothing
		return __work[device];
	}
	else {
		if (__work[device] != 0)
		{
			//cudaSetDevice(device);
			cudaFreeAsync(__work[device],stream);
			cudaMallocAsync(&__work[device], N * sizeof(double),stream);
		}
		work_size[device] = N;
		return __work[device];
	}
}
bool kingghidorah::cuda::valid() {
	//return true;
	return initialized && (!failed);
}
string  kingghidorah::cuda::device_name() {
	if (valid()) {
		char name[100];

		if (count() == 0) {
			failed = true;
			return "device not found";
		}
		std::stringstream ss;
		for (int i = 0; i < count(); i++)
		{
			CUdevice dev;
			cuDeviceGet(&dev, i);
			cuDeviceGetName(name, 100, dev);
			ss << name << "::" << "SPEED<<" << speed[i] << std::endl;
		}
		if (count() > 1)
		{
			ss << "multiple GPUs available!" << std::endl;
			for (int i = 0; i < count(); i++)
			{
				for (int j = i + 1; j < count(); j++)
				{
					ss << "peer access" << "(" << i << "," << j << ")" << canpeeraccess(i, j) << std::endl;
				}
			}

		}
		return ss.str();
	}
	return "invalid";
}
cusolverDnHandle_t& kingghidorah::cuda::solver(int ii,int kk) {
	return solver_handle[ii][kk];
}
cublasHandle_t& kingghidorah::cuda::blas(int ii) {
	return cublas_handle[ii];
}
inline int& kingghidorah::cuda::count() {
	return _count;
}
inline int& kingghidorah::cuda::fastest() {
	return _fastest;
}
kingghidorah::_mySparse::_mySparse()
{
	//_smat.resize(1);
	//_smat = new Eigen::SparseMatrix<double>();
}
kingghidorah::_mySparse::~_mySparse()
{
	if (__c * __r!= 0) {
		//cudaFreeHost(__dmat);
	}
	//delete _smat;
}

// TODO: This is an example of a library function
kingghidorah::_myPermutation::_myPermutation(int* ptr, int N)
{
	//Eigen::VectorXi indices(N);
	Eigen::VectorXi indices(N);
	for (int i = 0; i < N; i++)
	{
		indices[i] = *ptr;
		ptr++;
	}
	perm.indices() = indices;
}
std::string kingghidorah::_mySparse::_testopenmp()
{
	std::stringstream ss;
	int mt = omp_get_max_threads();
	ss << "num threads:" << mt << std::endl;
#pragma omp parallel for
	for (int i = 0; i < 100; i++)
	{
		int ct = omp_get_thread_num();
#pragma omp critical
		{
			ss << "current threads:" << ct << std::endl;
		}
		for (int t = 0; t < ct * 100; t++)
		{
			auto f = new double[10];
			delete[] f;
		}
	}
	return ss.str();
}
/*void kingghidorah::_mySparse::freeze(bool _do) {
	this->_dmat.setZero(this->rows(), this->cols());
	if(_do)
	for (int ii = 0; ii < _nt; ii++)
	{
		this->_dmat += coeff[ii].asDiagonal()*this->_mat[ii];
	}
	else
		for (int ii = 0; ii < _nt; ii++)
		{
			this->_dmat += this->_mat[ii];
		}
}*/
double kingghidorah::_mySparse::L2Norm(double* ptr1, int N1, double* ptr2, int N2) {
	auto a = Eigen::Map<Eigen::VectorXd, Eigen::Aligned128>(ptr1, N1);
	auto b = Eigen::Map<Eigen::VectorXd, Eigen::Aligned128>(ptr2, N2);
	return a.transpose() * this->_mat[0] * b;
}
void kingghidorah::_mySparse::Vector(Eigen::VectorXd *ptr1, int N1, Eigen::VectorXd* ptr2) {
	//auto a = Eigen::Map<Eigen::VectorXd, Eigen::Aligned128>(ptr1, N1);
	//auto b = Eigen::Map<Eigen::VectorXd, Eigen::Aligned128>(ptr1, N1);
	(*ptr2)=this->_mat[0] * (*ptr1);
}
void kingghidorah::_mySparse::plus(_mySparse* m, double sc,bool dense,bool sparse) {
	if (sparse)
		this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
	if (dense)
	{
		//Eigen::Map<Eigen::MatrixXd,Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
		_dmat += m->_mat[0] * sc;
	}
}
void kingghidorah::_mySparse::setmat(Eigen::SparseMatrix<double> mat, int ii) {
	this->_mat[ii] = mat;
}
void kingghidorah::_mySparse::setmat(const Eigen::MatrixXd& mat) {
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
	_dmat = mat;
}
double kingghidorah::_mySparse::at(int i, int ii) {
	return this->_mat[ii].data().value(i);
}
int kingghidorah::_mySparse::num_elem(int ii)
{
	return this->_mat[ii].data().size();
}
double kingghidorah::_mySparse::_at(int i) {
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
	return _dmat.data()[i];
}
double kingghidorah::_mySparse::_at(int i, int j) {
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
	return _dmat(i, j);
}
int kingghidorah::_mySparse::cols() {
	return _mat[0].cols();
}
void kingghidorah::_mySparse::_resize(int n, int m) {
	if (n < __c && m < __r) {
		__c = n;
		__r = m;
	}
	else {
		//if (__r * __c != 0)cudaFreeHost(__dmat);
		__c = n;
		__r = m;
		//cudaMallocHost(&__dmat, sizeof(double) * __c * __r*1.5);

	}
	this->_dmat.conservativeResize(n, m);
}
void kingghidorah::_mySparse::setmiddlecolum(Eigen::SparseMatrix<double> f, int start, int end) {
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
	_dmat.middleCols(start, end - start) = f;
}
void kingghidorah::_mySparse::permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm)
{
	for (int ii = 0; ii < _nt; ii++)
	{
		perm.transpose().applyThisOnTheRight(_mat[ii]);
	}
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm,bool sparse,bool dense)
{
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
	int nn = _dmat.rows();
	int numthreads = 0;
	numthreads = omp_get_max_threads();
	int S = nn / numthreads / 2;
	auto pt = perm.transpose();
	if(sparse)
	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			//_mat[0] = 
			perm.applyThisOnTheLeft(_mat[0]);
			perm.transpose().applyThisOnTheRight(_mat[0]);
		}
	}
	if (dense)
	{
#pragma omp parallel for
		for (int i = 0; i < nn; i += S)
		{
			int start = i;
			int end = i + S;
			if (end > nn)end = nn;
			//_dmat.middleRows(start, end - start) = ;//
			perm.transpose().applyThisOnTheRight(_dmat.middleRows(start, end - start));

		}

#pragma omp parallel for
		for (int i = 0; i < nn; i += S)
		{
			int start = i;
			int end = i + S;
			if (end > nn)end = nn;
			//_dmat.middleCols(start, end - start) = ;
			perm.applyThisOnTheLeft(_dmat.middleCols(start, end - start));
		}
	}
	/*_dmat.applyOnTheLeft(perm);
	_dmat.applyOnTheRight(perm.transpose());*/
}
void kingghidorah::_mySparse::shrink(int M)
{
	for (int ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii].leftCols(M);
}
void kingghidorah::_mySparse::_shrink(int M,bool sparse,bool dense)
{
	if (sparse)
	{
		if (_mat.size() >= 1)
		{
			//if ((_mat[0].rows() == _dmat.rows()) && (_mat[0].cols() == _dmat.cols()))
			{
				_mat[0].conservativeResize(M, M);
			}
		}
	}
	if (dense)
	{
		//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
		//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(this->__dmat, M, M);
		/*for (int i = 0; i < M; i++)
		{
			_dmat2.col(i) = _dmat.block(0, i, M, 1);
		}*/
		_dmat.conservativeResize(M, M);
		__c = M;
		__r = M;
	}
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm2)
{
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, __c, __r);
	int nn = _dmat.rows();
	int numthreads = omp_get_max_threads();

	int S = nn / numthreads / 2;
	auto pt = perm2.transpose();

	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0] = perm * (_mat[0]) * pt;
		}
	}
	/*
#pragma omp parallel for
	for (int i = 0; i < nn; i += S)
	{
		int start = i;
		int end = i + S;
		if (end > nn)end = nn;
		_dmat.middleRows(start, end - start) = _dmat.middleRows(start, end - start) * pt;

	}
#pragma omp parallel for
	for (int i = 0; i < nn; i += S)
	{
		int start = i;
		int end = i + S;
		if (end > nn)end = nn;
		_dmat.middleCols(start, end - start) = perm * _dmat.middleCols(start, end - start);

	}*/

}
void kingghidorah::_mySparse::_shrink(int M, int N)
{
	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0].conservativeResize(M, N);
		}
	}
	//_dmat.conservativeResize(M, N);//;// _dmat.topLeftCorner(M, N);// f;
}
int kingghidorah::_mySparse::rows() {
	int _ret = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		_ret += _mat[ii].rows();
	}
	return _ret;
}
int kingghidorah::_mySparse::_rows() {
	return __r;
}
int kingghidorah::_mySparse::_cols() {
	return __c;
}
int kingghidorah::_mySparse::__rows() {
	int _ret = 0;
	for (int ii = 0; ii < this->_coeff.size(); ii++)
	{
		_ret += _coeff[ii].size();
	}
	return _ret;
}
void kingghidorah::_mySparse::Clear() {
	if (_nt == 1)
	{
		this->dat[0].resize(0);
		this->_coeff[0].resize(0);
	}
	this->dat.resize(_nt);
	this->_coeff.resize(_nt);

	if (this->coeff.size() != _nt)this->coeff.resize(_nt);
	if (this->_mat.size() != _nt)this->_mat.resize(_nt);
#pragma omp parallel for
	for (int ii = 0; ii < _nt; ii++)
	{
		this->coeff[ii].setZero();
		this->_mat[ii].setZero();
	}
}

void kingghidorah::_mySparse::init(int n, int m)
{
	this->_nt = 1;

	resize(n, m);
	_resize(n, m);
}
int kingghidorah::_mySparse::resize(int n, int m) {
	this->_nt = this->dat.size();
	if (_nt == 0 || _nt == 1)
	{
		_nt = 1;
		_coeff.resize(_nt);
		dat.resize(_nt);
		_mat.resize(_nt);
		coeff.resize(_nt);
		_mat[0].resize(n, m);
	}

	return this->_nt;
}
void kingghidorah::_mySparse::reserve(int n) {
	_mat.reserve(n);
}
void kingghidorah::_mySparse::addemptyrow(int ii) {
	_coeff[0].push_back(1);
}
void kingghidorah::_mySparse::addrow(int ii, int* ptr, double* data, double sc, int N)
{
	addrow(ii, ptr, data, 0, sc, N, true);
}
void kingghidorah::_mySparse::addrow(int ii, int* ptr, double* data, int shift, double sc, int N, bool add)
{
	data += shift;
	ptr += shift;
	if (dat.size() != 1)dat.resize(1);
	if (_coeff.size() != 1)_coeff.resize(1);

	for (int i = 0; i < N; i++)
	{
		dat[0].push_back(Eigen::Triplet<double>(ii, *ptr, (*data)));
		ptr++;
		data++;
	}
	if (add)
	{
		_coeff[0].push_back(sc);
	}
}

void kingghidorah::_mySparse::adddat(int ii, int j, double value)
{
	dat[0].push_back(Eigen::Triplet<double>(ii, j, value));
}
void kingghidorah::_mySparse::addcoeff(double sc) {
	_coeff[0].push_back(sc);
}
Eigen::VectorXd kingghidorah::_mySparse::get_coeff(int ii) {
	return this->coeff[ii];
}
void  kingghidorah::_mySparse::copycoefffrom(kingghidorah::_mySparse* mat)
{
	this->_nt = mat->_nt;
	this->coeff.resize(_nt);
	for (int ii = 0; ii < _nt; ii++)
	{
		this->coeff[ii] = Eigen::VectorXd(mat->get_coeff(ii));
	}
}

void kingghidorah::_mySparse::begin_construct()
{
	_dat_count = 0;
}
void kingghidorah::_mySparse::end_construct(int cc)
{
	_nt = _dat_count;

	_mat.resize(_nt);
	coeff.resize(_nt);
	for (int ii = 0; ii < _nt; ii++)
	{
		_mat[ii].resize(_coeff[ii].size(), cc);
	}
}
void kingghidorah::_mySparse::addmat(_mySparse* mat)
{
	_dat_count++;
	if (dat.size() < _dat_count)dat.resize(_dat_count);
	if (_coeff.size() < _dat_count)_coeff.resize(_dat_count);
	dat[_dat_count - 1] = mat->dat[0];
	_coeff[_dat_count - 1] = (mat->_coeff[0]);
}
void kingghidorah::_mySparse::OfDuplicate(_mySparse* mat)
{
	this->_nt = mat->_nt;
	this->_mat.resize(_nt);
	this->_coeff.resize(_nt);
	this->coeff.resize(_nt);
	for (int ii = 0; ii < _nt; ii++)
	{
		//Eigen::SparseMatrix<double> M(mat->_mat[ii]);
		//std::vector<double> vec(mat->_coeff[ii]);
		
		this->_mat[ii] = mat->_mat[ii];// M;
		this->_coeff[ii] = mat->_coeff[ii];// vec;
	}
}
void kingghidorah::_mySparse::_OfDuplicate(_mySparse* mat)
{
	if (this->_mat.size() == 0)
	this->_mat.resize(1);
	this->_mat[0] = mat->_mat[0];
	//this->_dmat = this->_mat[0];// mat->_dmat;
}
void kingghidorah::_mySparse::ofDat()
{
	if (_mat.size() != _nt)_mat.resize(_nt);
#pragma omp parallel for
	for (int ii = 0; ii < _nt; ii++)
	{
		_mat[ii].reserve(dat[ii].size());
		_mat[ii].setFromTriplets(dat[ii].begin(), dat[ii].end());
	}
}
void kingghidorah::_mySparse::freezecoeff() {
#pragma omp parallel for
	for (int ii = 0; ii < _nt; ii++)
	{
		coeff[ii] = Eigen::Map<Eigen::VectorXd, Eigen::Aligned128>(_coeff[ii].data(), _coeff[ii].size());
	}
}
int kingghidorah::_mySparse::numBlocks()
{
	return this->dat.size();
}
std::vector<Eigen::SparseMatrix<double>> e;
//std::vector<Eigen::MatrixXd> e;
int kingghidorah::_mySparse::ofAtA(_mySparse* A,bool sparse)
{
	int nn = A->cols();
	int mt = omp_get_max_threads();
	_mt = mt;

	if (e.size() < mt)
	{
		e.resize(mt);
	}
	for (int i = 0; i < mt; i++) {
		e[i].resize(nn, nn);
		e[i].setZero();
	}
	this->_resize(nn, nn);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(this->__dmat, nn, nn);
	//Eigen::MatrixXd ___dmat(nn, nn);
//#pragma omp parallel sections
//	{
//#pragma omp section
		{
			//memset(this->__dmat, 0, sizeof(double) * nn * nn);
			if (!sparse)_dmat.setZero();
		}
//#pragma omp section
	//	{
#pragma omp parallel for
	for (int _ii = 0; _ii < mt; _ii++)
	{
		int S = 0;
		int E = 0;
		auto _e = e[_ii];
		S = _ii * _nt / mt;
		E = (_ii + 1) * _nt / mt;
		for (int ii = S; ii < E; ii++)
		{
			e[_ii] += (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
		}
	}
	if (sparse) {
		if (this->_mat.size() == 0)this->_mat.resize(1);
		this->_mat[0].resize(nn, nn);
		this->_mat[0].setZero();
		for (int i = 0; i < mt; i++) {
			this->_mat[0] += e[i];
		}
		//this->_dmat = this->_mat[0];
	}
	else {
		this->_dmat.setZero(nn, nn);
		for (int i = 0; i < mt; i++) {
			this->_dmat += e[i];
		}
	}
	//_dmat = ___dmat;
	//this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	return _nt;
}
void kingghidorah::_mySparse::_freeze() {
	//this->_dmat = this->_mat[0];
}
std::string kingghidorah::_mySparse::_ofAtA(_mySparse* A)
{
	std::stringstream ss;
#ifdef _DEBUG
	ss << _nt << std::endl;
	for (int ii = 0; ii < _nt; ii++)
	{
		ss << ii << "th" << std::endl;
		ss << A->_mat[ii].rows() << "," << A->_mat[ii].cols() << "," << coeff[ii].size() << std::endl;
	}
#endif
	this->_resize(A->cols(), A->cols());
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, A->cols(), A->cols());
	__c = A->cols();
	__r = A->cols();
	//this->_dmat.resize(A->cols(), A->cols());
	_dmat.setZero();
	for (int ii = 0; ii < _nt; ii++)
	{
		auto _ret = (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
		_dmat += _ret;
	}
	return ss.str();
}
std::string kingghidorah::_mySparse::info()
{
	return "viennacl has been abandoned";
}

/*void kingghidorah::_mySparse::_ofAtB_gpu(cuda* cuda, _mySparse* B, _mySparse* C)
{
	this->freeze(true);
	B->freeze(false);
	int nn = this->_cols();
	int mm = B->_cols();
	int kk = this->_rows();
	C->_dmat.resize(nn, mm);
	C->_dmat.setZero();
	int mt = omp_get_max_threads();

	int ss = mm / cuda->count() / 4;
	if (!cuda->valid())return;
	int job = 0;
	double alpha = 1.0;
	double beta = 1.0;
	int cc = cuda->count();
#pragma omp parallel for
	for (int ii = 0; ii < cc; ii++)
	{
		cudaSetDevice(ii);
		auto cublas = cuda->blas(ii);

		int S = 0;
		int E = 0;
		double* matrix_A_gpu = cuda->work_M(ii);
		double* matrix_B_gpu = cuda->work_rhs(ii);
		double* matrix_C_gpu = cuda->work_C(ii);
		cudaMemset(matrix_C_gpu, 0, mm * nn * sizeof(double));
		cudaMemcpy(matrix_A_gpu, this->_dmat.data(), nn * kk * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(matrix_B_gpu, B->_dmat.data(), mm * kk * sizeof(double), cudaMemcpyHostToDevice);
		cudaDeviceSynchronize();
		while (true)
		{
#pragma omp critical
			{
				S = job;
				E = job + ss;
				if (E >= mm)E = mm;
				job = E;
			}
			if (S >= mm)break;
			cublasDgemm(cublas, cublasOperation_t::CUBLAS_OP_T, cublasOperation_t::CUBLAS_OP_N, nn, E - S, kk, &alpha, matrix_A_gpu, kk, matrix_B_gpu + S * kk, kk, &beta, matrix_C_gpu + S * nn, nn);
			cudaMemcpy(C->_dmat.data() + S * nn, matrix_C_gpu + S * nn, (E - S) * nn * sizeof(double), cudaMemcpyDeviceToHost);
		}

	}
	if (C->_mat.size() == 0)C->_mat.resize(1);
	C->_mat[0] = C->_dmat.sparseView(1.0, 0.0000000000001);
}*/
/*void kingghidorah::_mySparse::_ofAtB(_mySparse* B, _mySparse* C)
{
	this->freeze(true);
	B->freeze(false);

	int nn = this->_cols();
	int mm = B->_cols();
	C->_dmat.resize(nn, mm);
	C->_dmat.setZero();

	int mt = omp_get_max_threads();

	int ss = mm / mt / 2;
	auto left = this->_dmat.transpose();
	auto right = B->_mat[0];

#pragma omp parallel for
	for (int ii = 0; ii < mm; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= mm)E = mm;
		C->_dmat.middleCols(S, E - S) = left * right.middleCols(S, E - S);
	}
	if (C->_mat.size() == 0)C->_mat.resize(1);
	C->_mat[0] = C->_dmat.sparseView(1.0, 0.0000000000001);
}*/
Eigen::VectorXd kingghidorah::_mySparse::_ofBtAB(_mySparse* B, double* ptr, int N, _mySparse* C)
{
	static Eigen::MatrixXd D;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	int nn = B->_mat[0].cols();
	int kk = _dmat.cols();
	C->_resize(nn, nn);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(C->__dmat, nn, nn);

	//C->_dmat.resize(nn, nn);

	int mt = omp_get_max_threads();

	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];

	D.resize(nn, kk);
	int ss = kk / mt / 2;
#pragma omp sections
	{
#pragma omp section
		{
			//_dmat2.setZero();
			C->_dmat.setZero();
			D.setZero();
		}
#pragma omp section
		{
#pragma omp parallel for
			for (int ii = 0; ii < kk; ii += ss)
			{
				int S = ii;
				int E = ii + ss;
				if (E >= kk)E = kk;
				D.middleCols(S, E - S) = left * mid.middleCols(S, E - S);
			}
		}
	}
	ss = nn / mt / 2;
#pragma omp parallel for
	for (int ii = 0; ii < nn; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= nn)E = nn;
		C->_dmat.middleCols(S, E - S) = D * right.middleCols(S, E - S);
	}
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> b(ptr, N);
	return D * b;
}

void kingghidorah::_mySparse::ofAtB(_mySparse* B, bool sparse)
{

	int nn = this->cols();
	int mm = B->cols();
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, nn,mm);
	__r = nn;
	__c = mm;
	//this->_dmat.resize(nn, mm);
	int mt = omp_get_max_threads();

	_mt = mt;
	if (mt > e.size())
		e.resize(mt);
	for (int i = 0; i < mt; i++) {
		e[i].resize(nn, mm);
		e[i].setZero();
	}
#pragma omp parallel for
	for (int _ii = 0; _ii < mt; _ii++)
	{
		int S = _ii * _nt / mt;
		int E = (_ii + 1) * _nt / mt;

		for (int ii = S; ii < E; ii++)
		{
			e[_ii] += this->_mat[ii].transpose() * coeff[ii].asDiagonal() * B->_mat[ii];
		}
	}
	if (_mat.size() == 0)_mat.resize(1);
	this->_mat[0].resize(nn, mm);
	this->_mat[0].setZero();

	if (sparse)
	{
		for (int i = 0; i < mt; i++) {
			this->_mat[0] += e[i];
		}
	}
	else {
		for (int i = 0; i < mt; i++) {
			this->_dmat += e[i];
		}
	}
	//this->_dmat = this->_mat[0];
}

Eigen::VectorXd kingghidorah::_mySparse::Atb(double* ptr, int N)
{
	Eigen::VectorXd ret(this->cols());
	ret.setZero();
	int offset = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		int ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> b(ptr + offset, ee);
		ret += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		offset += ee;
	}
	return ret;
}
void kingghidorah::_mySparse::Atb(double* ptr, int N,Eigen::VectorXd *ret)
{
	//static Eigen::VectorXd ret;
	//auto ret = Eigen::Map<Eigen::VectorXd, Eigen::Aligned128>(_ret, this->cols());
	ret->setZero();
	int offset = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		int ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		(*ret) += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		offset += ee;
	}
}

Eigen::VectorXd kingghidorah::_mySparse::_Atb(double* ptr, int N)
{
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> b(ptr, N);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r,__c);
	return _dmat.transpose() * b;
}
void kingghidorah::_mySparse::merge()
{
	this->ofDat();
}
void kingghidorah::_mySparse::computeQR()
{
	Eigen::HouseholderQR<Eigen::MatrixXd> qr;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	qr.compute(_dmat);
}
void kingghidorah::_mySparse::computeLU()
{
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	lu.compute(_dmat);
}
void kingghidorah::_mySparse::computeLLT(Eigen::LLT<Eigen::MatrixXd>* _LLT)
{
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	_LLT->compute(_dmat);
}
int kingghidorah::_mySparse::nonzeros() {
	int _ret = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		_ret += _mat[ii].nonZeros();
	}
	return _ret;
}
Eigen::VectorXd kingghidorah::_mySparse::solve0(double* rhs, int N) {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	lu.compute(_dmat);
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> b(rhs, N);
	Eigen::VectorXd x(_dmat.cols());
	x.setZero();
	x = lu.solve(b);
	return x;
}

void kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, double* rhs, int N, double *ret,int device) {
	//Eigen::VectorXd x(N);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> x(ret, N);
	this->_freeze();

	if (!cuda->valid())
	{
		x(0) = 10;
		return;
	}
	cudaSetDevice(device);
	auto solver = cuda->solver(device,0);
	auto blas = cuda->blas(device);
	//auto stream=streams[device];
	//cudaStreamCreate(&stream);
	//cusolverDnSetStream(solver, stream);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	x.setZero();


	double* gpu_rhs = cuda->work_rhs(device);
	double* gpu_matrix = cuda->work_M(device);
	cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_rhs, rhs, N * sizeof(double), cudaMemcpyHostToDevice);

	int work_size = 0;
	int* devInfo_on_gpu = cuda->info(device);
	//cudaMalloc(&devInfo_on_gpu, sizeof(int));

	// --- CUDA CHOLESKY initialization
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);

	// --- CUDA POTRF execution
	double* work = cuda->work(work_size, device);

	auto now = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
	int devInfo_on_cpu = 0;
	cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);


	if (0 != devInfo_on_cpu) {
		x(0) = devInfo_on_cpu;
		return;
	}


	cusolverDnDpotrs(solver, CUBLAS_FILL_MODE_LOWER, N, 1, gpu_matrix, N, gpu_rhs, N, devInfo_on_gpu);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Dn:" << duration.count() << "ms" << std::endl;

	cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);

	if (devInfo_on_cpu != 0) {
		x(0) = 24;
		return;
	}

	cudaMemcpy(x.data(), gpu_rhs, sizeof(double) * N, cudaMemcpyDeviceToHost);

	//cudaFree(work);

	//cudaFree(devInfo_on_gpu);

	cudaDeviceSynchronize();

}
Eigen::MatrixXd kingghidorah::_mySparse::_solve0(_myLLT* LLT, _mySparse* mat)
{
	//this function assumes that LLT decomposition has been done already
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(mat->__dmat, __r, __c);

	Eigen::MatrixXd ret(_dmat.cols(),_dmat.rows());
	int nn = mat->cols();
	int numthreads = 0;
#pragma omp parallel
	{
#pragma omp single
		numthreads = omp_get_num_threads();
	}
	int S = nn / numthreads / 2;

	//ret = LLT->LLT->solve(mat->_dmat);

#pragma omp parallel for
	for (int i = 0; i < nn; i += S)
	{
		int start = i;
		int end = i + S;
		if (end > nn)end = nn;
		//ret.middleRows(start, end - start) = LLT->LLT->solve(_dmat2.middleCols(start, end - start)).transpose();
		ret.middleRows(start, end - start) = LLT->LLT->solve(mat->_dmat.middleCols(start, end - start)).transpose();

	}
	return ret.transpose();
}

void kingghidorah::_mySparse::_solveI_gpu_mg(kingghidorah::cuda* cuda, _mySparse* ret)
{
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	int N = _dmat.cols();
	//Eigen::MatrixXd x(N,nn);
	ret->__c = N;
	ret->__r = N;
	//ret->_dmat.resize(N, N);
	if (!cuda->valid())return;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(ret->__dmat, N, N);

	//_dmat2.setZero();
	ret->_dmat.setZero();
	cudaSetDevice(cuda->fastest());

	cusolverMgHandle_t solver = cuda->mgsolver();
	//auto status = cusolverMgCreate(&solver);

	int nbGpus = cuda->count();
	//std::vector<int> _deviceList(nbGpus);
	//for (int i = 0; i < nbGpus; i++)_deviceList[i] = i;
	int* deviceList = cuda->devicelist();
	//enablePeerAccess(nbGpus, deviceList);

	const int IA = 1;
	const int JA = 1;
	const int T_A = N;// / nbGpus;
	const int lda = N;

	int  info = 0;

	cudaLibMgMatrixDesc_t descrA=cuda->descrA;

	cudaLibMgGrid_t gridA=cuda->gridA;

	cusolverMgGridMapping_t mapping = cuda->mapping;

	double** array_d_A = cuda->array_d_A();


	int64_t lwork_potrf = 0;
	int64_t lwork_potri = 0;
	int64_t lwork = 0;

	/*auto status = cusolverMgDeviceSelect(
		solver,
		nbGpus,
		deviceList);
	assert(CUSOLVER_STATUS_SUCCESS == status);
	{
	}*///always fails
	//auto status = cusolverMgCreateDeviceGrid(&gridA, 1, nbGpus, deviceList, mapping);
	//assert(CUSOLVER_STATUS_SUCCESS == status);


	//assert(CUSOLVER_STATUS_SUCCESS == status);



	/*array_d_A = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_A);
	array_d_B = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_B);
	*/
	/*if (cuda->prevN < N || cuda->prevT_A < T_A)
	{
		if (cuda->prevN != 0)
			destroyMat(
				nbGpus,
				deviceList,
				cuda->prevN,
				cuda->prevT_A,
				(void**)array_d_A);

	

		cuda->prevN = N;
		cuda->prevT_A = T_A;
	}*/
	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		N,

		//__dmat,
		_dmat.data(),
		lda,

		N,
		T_A,
		lda,
		array_d_A,
		IA,
		JA
		);


	
	auto cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	auto now = std::chrono::high_resolution_clock::now();
	auto status = cusolverMgPotrf(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		(void**)cuda->array_d_work(),
		lwork,
		&info
	);
	if (CUSOLVER_STATUS_SUCCESS != status)
	{
		assert(CUSOLVER_STATUS_ALLOC_FAILED != status);
		assert(CUSOLVER_STATUS_ARCH_MISMATCH != status);
		assert(CUSOLVER_STATUS_EXECUTION_FAILED != status);
		assert(CUSOLVER_STATUS_INTERNAL_ERROR != status);
		assert(CUSOLVER_STATUS_INVALID_LICENSE != status);
		assert(CUSOLVER_STATUS_INVALID_VALUE != status);
		assert(CUSOLVER_STATUS_INVALID_WORKSPACE != status);
		assert(CUSOLVER_STATUS_IRS_INFOS_NOT_DESTROYED != status);
		assert(CUSOLVER_STATUS_IRS_INFOS_NOT_INITIALIZED != status);
		assert(CUSOLVER_STATUS_IRS_INTERNAL_ERROR != status);
		assert(CUSOLVER_STATUS_IRS_MATRIX_SINGULAR != status);
		assert(CUSOLVER_STATUS_IRS_NOT_SUPPORTED != status);
		//_dmat2(0, 0) = info;
		ret->_dmat(0, 0) = info;
		return;
	}

	cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	assert(0 == info);



	status = cusolverMgPotri(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		(void**)cuda->array_d_work(),
		lwork,
		&info
	);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Mg" << duration.count() << "ms" << std::endl;

	assert(CUSOLVER_STATUS_SUCCESS == status);
	cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	assert(0 == info);

	printf("step 11: solution vector B \n");
	memcpyD2H<double>(
		nbGpus,
		deviceList,
		N,
		N,

		N,
		T_A,
		lda,
		array_d_A,
		IA,
		JA,

		//ret->__dmat,
		ret->_dmat.data(),
		lda
		);


#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			//_dmat2(j, i) = _dmat2(i, j);
			ret->_dmat(j, i) = ret->_dmat(i, j);
		}
	}


	//if (NULL != array_d_A) free(array_d_A);
	//if (NULL != array_d_B) free(array_d_B);
	//if (NULL != array_d_work) free(array_d_work);

}

void kingghidorah::_mySparse::_solveI(_mySparse* ret)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(__dmat, __r, __c);

	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);	
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
	llt.compute(this->_mat[0]);
	int nn = this->_mat[0].rows();
	//ret->_dmat.resize(nn, nn);
	ret->__c = nn;
	ret->__r = nn;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(ret->__dmat, nn, nn);
	if (I.rows() != nn || I.cols() != nn)
	{
		I.resize(nn, nn);
		I.setIdentity();
	}
	int mt = omp_get_max_threads();
	int ee = mt * 4;
#pragma omp parallel for
	for (int i = 0; i < ee; i++)
	{
		int S = i * nn / ee;
		int E = (i + 1) * nn / ee;
		//_dmat2.middleCols(S, E - S) = llt.solve(I.middleCols(S, E - S));
		ret->_dmat.middleCols(S, E - S) = llt.solve(I.middleCols(S, E - S));
	}
}

/*inline double inner_sum(double* li, double* lj, int n) {
	double s = 0;
	for (int i = 0; i < n; i++) {
		s += li[i] * lj[i];
	}
	return s;
}

inline double inner_sum3(double* li, double* lj, int n) {
	double s1 = 0, s2 = 0, s3 = 0;
	int i;
	for (i = 0; i < (n & (-3)); i += 3) {
		s1 += li[i] * lj[i + 0];
		s2 += li[i + 1] * lj[i + 1];
		s3 += li[i + 2] * lj[i + 2];
	}
	double sum = 0;
	for (; i < n; i++) sum += li[i] * lj[i + 0];
	sum += s1 + s2 + s3;
	return sum;
}

void cholesky4(double* A, double* L,int n) {

	if (L == NULL)
		exit(EXIT_FAILURE);

	for (int j = 0; j < n; j++) {

		double s = 0;
		for (int k = 0; k < j; k++) {
			s += L[j * n + k] * L[j * n + k];
		}
		L[j * n + j] = sqrt(A[j * n + j] - s);
#pragma omp parallel for
		for (int i = j + 1; i < n; i++) {
			double s = 0;
			for (int k = 0; k < j; k++) {
				s += L[i * n + k] * L[j * n + k];
			}
			L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
		}
	}
}

void cholesky5(double* A, double* L,int n) {

	for (int j = 0; j < n; j++) {
		double s = inner_sum(&L[j * n], &L[j * n], j);
		L[j * n + j] = sqrt(A[j * n + j] - s);
#pragma omp parallel for schedule(static, 8)
		for (int i = j + 1; i < n; i++) {
			double s = inner_sum(&L[j * n], &L[i * n], j);
			L[i * n + j] = (1.0 / L[j * n + j] * (A[i * n + j] - s));
		}
	}
}

void cholesky2(double* A, double* L,int n) {

	if (L == NULL)
		exit(EXIT_FAILURE);

	for (int i = 0; i < n; i++) {
		double s = 0;
		for (int k = 0; k < i; k++) {
			s += L[k * n + i] * L[k * n + i];
		}
		L[i * n + i] = sqrt(A[i * n + i] - s);
#pragma omp parallel for
		for (int j = i + 1; j < n; j++) {
			double s = 0;
			for (int k = 0; k < i; k++) {
				s += L[k * n + i] * L[k * n + j];

			}
			L[i * n + j] = (1.0 / L[i * n + i] * (A[i * n + j] - s));
		}
	}

}*/

void initidentiy(kingghidorah::cuda* cuda,int N) {
	double _a = 0;
	double _b = 1;
	//if (I.cols() != N || I.rows() != N)
	{
		//I.setIdentity(N, N);
#pragma omp parallel for
		for (int ii = 0; ii < cuda->count(); ii++)
		{
			cudaSetDevice(ii);
			double* gpu_rhs = cuda->work_C(ii);
			auto stream = cuda->__streams(ii, 0);
			//cudaMemcpyAsync(gpu_rhs, I.data(), sizeof(double) * N * N, cudaMemcpyHostToDevice, stream);
			cudaMemsetAsync(gpu_rhs, 0, N * N * sizeof(double), stream);
			for (int i = 0; i < N; i++)
			{
				cudaMemcpyAsync(gpu_rhs + (i+i*N), &_b, 1*sizeof(double), cudaMemcpyHostToDevice,stream);
			}
		}
	}
}
std::string kingghidorah::_mySparse::_solveI_gpu_single(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();
	std::stringstream sss;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);

	int N = _dmat.cols();

	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(ret->__dmat, N, N);
	ret->__r = N;
	ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int nn = cuda->count();
	initidentiy(cuda, N);
	//_dmat2.setZero();
	ret->_dmat.setZero();


	

		auto solver = cuda->solver(cuda->fastest(),0);

		cudaSetDevice(cuda->fastest());
		//auto stream = streams[cuda->fastest()];
		double* gpu_matrix = cuda->work_M(cuda->fastest());
		cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
		int* devInfo_on_gpu = cuda->info(cuda->fastest());

		int work_size = 0;

		// --- CUDA CHOLESKY initialization
		cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);
		// --- CUDA POTRF execution	
		double* work = cuda->work(work_size,cuda->fastest());

		cudaMemset(work, 0, work_size * sizeof(double));
		//cusolverDnSetStream(solver, stream);
		cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
		cudaDeviceSynchronize();
	



		double* gpu_rhs = cuda->work_C(cuda->fastest());


		int devInfo_on_cpu = 0;
		//int _S = 0;
		//int _E = 0;

		/*for (int kk = 0; kk < 16; kk++)
		{
			cudaStreamSynchronize(cuda->__streams(cuda->fastest(), kk));
		}*/
		bool exit = false;
#pragma omp parallel for
	for (int kk = 0; kk < STRTREAMCOUNT; kk++)
	{
		int S = kk * N / STRTREAMCOUNT;
		int E = (kk + 1) * N / STRTREAMCOUNT;
		cudaStream_t _stream = cuda->__streams(cuda->fastest(), kk);
		cusolverDnHandle_t _solver = cuda->solver(cuda->fastest(), kk);
		cusolverDnSetStream(_solver, _stream);

		cusolverDnDpotrs(_solver, CUBLAS_FILL_MODE_LOWER, N, E - S, gpu_matrix, N, gpu_rhs + S * N, N, devInfo_on_gpu);
		//cudaMemcpyAsync(ret->__dmat + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
		cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
	}

	cudaDeviceSynchronize();
	return sss.str();

}
std::string kingghidorah::_mySparse::_solveI_gpu_omp(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();
	std::stringstream sss;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);

	int N = _dmat.cols();
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(ret->__dmat,N, N);
	ret->__r = N;
	ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int nn = cuda->count();
	initidentiy(cuda, N);
	//_dmat2.setZero();
	ret->_dmat.setZero();


	if (cuda->canpeer())
	{

		auto solver = cuda->solver(cuda->fastest(), 0);

		cudaSetDevice(cuda->fastest());
		//auto stream = streams[cuda->fastest()];
		double* gpu_matrix = cuda->work_M(cuda->fastest());
		double* gpu_rhs = cuda->work_C(cuda->fastest());
		cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
		int* devInfo_on_gpu = cuda->info(cuda->fastest());

		int work_size = 0;

		// --- CUDA CHOLESKY initialization
		cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);
		// --- CUDA POTRF execution	
		double* work = cuda->work(work_size, cuda->fastest());

		cudaMemset(work, 0, work_size * sizeof(double));
		//cusolverDnSetStream(solver, stream);
		cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
		cudaDeviceSynchronize();
#pragma omp parallel for
		for (int ii = 0; ii < nn; ii++)
		{
			if (ii != cuda->fastest())
			{
				cudaStream_t _stream = cuda->__streams(ii, 0);
				//cudaStreamCreate(&_stream);
				cudaMemcpyAsync(cuda->work_M(ii), gpu_matrix, N * N * sizeof(double), cudaMemcpyDeviceToDevice, _stream);
				//cudaStreamDestroy(_stream);
			}
		}
	}
	else {
		auto solver = cuda->solver(cuda->fastest(), 0);

		cudaSetDevice(cuda->fastest());
		//auto stream = streams[cuda->fastest()];
		double* gpu_matrix = cuda->work_M(cuda->fastest());
		// --- CUDA CHOLESKY initialization
		int work_size = 0;
		cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);
		cudaMemcpyAsync(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice, cuda->__streams(cuda->fastest(), 0));
		int* devInfo_on_gpu = cuda->info(cuda->fastest());


		// --- CUDA POTRF execution	
		double* work = cuda->work(work_size, cuda->fastest(), cuda->__streams(cuda->fastest(), 1));

		cudaMemsetAsync(work, 0, work_size * sizeof(double), cuda->__streams(cuda->fastest(), 1));
		cuStreamSynchronize(cuda->__streams(cuda->fastest(), 0));
		cuStreamSynchronize(cuda->__streams(cuda->fastest(), 1));
		cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
		cudaDeviceSynchronize();
#pragma omp parallel for
		for (int ss = 0; ss < STRTREAMCOUNT; ss++)
		{
			cudaSetDevice(cuda->fastest());
			int S = ss * N / STRTREAMCOUNT;
			int E = (ss + 1) * N / STRTREAMCOUNT;
			cudaStream_t _streamX = cuda->__streams(cuda->fastest(), ss);
			//cudaMemcpyAsync(ret->__dmat + S * N, gpu_matrix + S * N, sizeof(double) * N * (E - S), cudaMemcpyDeviceToHost, _streamX);
			cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_matrix + S * N, sizeof(double) * N * (E - S), cudaMemcpyDeviceToHost, _streamX);
			cudaStreamSynchronize(_streamX);
			//cudaDeviceSynchronize();

#pragma omp parallel for
			for (int ii = 0; ii < nn; ii++)
			{
				if (ii != cuda->fastest())
				{
					cudaSetDevice(ii);
					cudaStream_t _stream = cuda->__streams(ii, ss);
					//cudaMemcpyAsync(cuda->work_M(ii) + S * N, ret->__dmat + S * N, sizeof(double) * N * (E - S), cudaMemcpyHostToDevice, _stream);
					cudaMemcpyAsync(cuda->work_M(ii) + S * N, ret->_dmat.data() + S * N, sizeof(double) * N * (E - S), cudaMemcpyHostToDevice, _stream);
				}
			}
		}
		//cudaStreamDestroy(_stream);
	}

	int job = 0;
	int ss = N / cuda->count() / STRTREAMCOUNT/4;
	if (ss == 0) ss = 1;
	for (int i = 0; i < nn; i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
	}
	std::vector<int>count(nn);
#pragma omp parallel for
	for (int ii = 0; ii < nn; ii++)
	{
		auto solver = cuda->solver(ii,0);
		
		cudaSetDevice(ii);
		//auto stream=streams[ii];
		//cusolverDnSetStream(solver, stream);

		double* gpu_matrix = cuda->work_M(ii);
		double* gpu_rhs = cuda->work_C(ii);
		int* devInfo_on_gpu = cuda->info(ii);

		int devInfo_on_cpu = 0;

		bool exit = false;
		while (true)
		{
#pragma omp parallel for
			for (int kk = 0; kk < STRTREAMCOUNT; kk++)
			{
				int S = 0;
				int E = 0;

#pragma critical
				{
					S = job;
					E = S + ss;
					if (E > N)E = N;
					job = E;
				}
				if (S >= N) {
					exit = true;
					break;
				}
				count[ii]++;
				//int S = _S + (_E - _S) * kk / 4;
				//int E = _S + (_E - _S) * (kk+1) / 4;
				cudaStream_t _stream = cuda->__streams(ii, kk);
				auto _solver = cuda->solver(ii, kk);
				cusolverDnSetStream(_solver, _stream);

				cusolverDnDpotrs(_solver, CUBLAS_FILL_MODE_LOWER, N, E - S, gpu_matrix, N, gpu_rhs + S * N, N, devInfo_on_gpu);
				//cudaMemcpyAsync(ret->__dmat + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
				cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_rhs + S * N, (E - S)* N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
			}
			if (exit)break;
		}

		//cusolverDnDestroy(_solver);
		//cusolverDnSetStream(solver,streams[ii]);
	}
	cudaDeviceSynchronize();
	sss << "jobs taken: ";
	for (int i = 0; i < nn; i++)
	{
		sss << i<<"-" << count[i] << "  ";
	}
	sss << std::endl;
	return sss.str();
	/*
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			ret->_dmat(j, i) = ret->_dmat(i, j);
		}
	}
	*/

}

void kingghidorah::_mySparse::_solveI_gpu(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);

	int N = _dmat.cols();

	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(ret->__dmat, N, N);
	ret->__r = N;
	ret->__c = N;
	//ret->_dmat.resize(N, N);
	//_dmat2.setZero();
	ret->_dmat.setZero();
	//initidentiy(cuda, N);
	if (!cuda->valid())return;
	
		int ii = cuda->fastest();
		auto solver = cuda->solver(ii,0);
		//cudaStream_t stream = streams[ii];
		//cusolverDnSetStream(solver, stream);
		cudaSetDevice(ii);
		double* m_gpu = cuda->work_M(ii);

		cudaMemcpyAsync(m_gpu, this->_dmat.data(), sizeof(double) * N * N, cudaMemcpyHostToDevice,cuda->__streams(cuda->fastest(),0));
		double* work;
		int work_size;
		int work_size1 = 0;
		int work_size2 = 0;
		cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size1);
		cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size2);
		work_size = std::max(work_size1, work_size2);
		work = cuda->work(work_size, ii, cuda->__streams(cuda->fastest(), 1));
		//cudaMallocAsync(&work, sizeof(double) * work_size, stream);
		cudaMemsetAsync(work, 0, sizeof(double) * work_size1, cuda->__streams(cuda->fastest(), 1));
		int* devInfo = cuda->info(ii);
		cudaDeviceSynchronize();
		cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work, work_size1, devInfo);
		double* work2 = cuda->work_rhs(ii);
		//cudaMalloc(&work2, sizeof(double) * work_size2);
		cudaMemset(work2, 0, sizeof(double) * work_size2);
		cusolverDnDpotri(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work2, work_size2, devInfo);
		//cudaFree(work2);

		//cudaMemcpy(ret->__dmat, m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
		cudaMemcpy(ret->_dmat.data(), m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();

	
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			//_dmat2(j, i) = _dmat2(i, j);
			ret->_dmat(j, i) = ret->_dmat(i, j);
		}
	}
	
}

/*void kingghidorah::_mySparse::_solveI_gpu(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();

	int N = this->_dmat.cols();

	ret->_dmat.resize(N, N);
	ret->_dmat.setZero();
	initidentiy(cuda, N);
	if (!cuda->valid())return;

	int ii = cuda->fastest();
	auto solver = cuda->solver(ii);
	auto blas = cuda->blas(ii);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	int job = 0;

	cudaSetDevice(ii);
	auto stream = streams[ii];
	cusolverDnSetStream(solver, stream);
	//cudaStreamCreate(&stream);
	//cusolverDnSetStream(solver, stream);

	double* gpu_matrix = cuda->work_M(ii);
	double* gpu_rhs = cuda->work_rhs(ii);
	cudaMemcpyAsync(gpu_matrix, this->_dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice, stream);
	int* devInfo_on_gpu = cuda->info(ii);

	int work_size = 0;
	int work_size1 = 0;
	int work_size2 = 0;

	// --- CUDA CHOLESKY initialization
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size1);
	//cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size2);
	work_size = work_size1;// std::max(work_size1, work_size2);
	// --- CUDA POTRF execution	
	double* work = cuda->work(work_size, ii);

	cudaMemsetAsync(work, 0, work_size1 * sizeof(double), stream);
	auto start = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size1, devInfo_on_gpu);
	int devInfo_on_cpu = 0;
	cudaMemcpyAsync(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int) * 1, cudaMemcpyDeviceToHost, stream);
	if (devInfo_on_cpu != 0)
	{
		ret->_dmat.data()[0] = devInfo_on_cpu;
		return;
	}
	auto end = high_resolution_clock::now();
	auto duration = end - start;
	std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
	std::cout << "regular:Dpotrf___AAA___:" << d.count() << "ms" << std::endl;

	//cudaStreamSynchronize(stream);
	//cudaDeviceSynchronize(); 
	//start = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrs(solver, CUBLAS_FILL_MODE_LOWER, N, N, gpu_matrix, N, gpu_rhs, N, devInfo_on_gpu);
	cudaMemcpyAsync(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int) * 1, cudaMemcpyDeviceToHost, stream);
	if (devInfo_on_cpu != 0)
	{
		ret->_dmat.data()[0] = devInfo_on_cpu;
		return;
	}
	cudaMemcpyAsync(ret->_dmat.data(), gpu_rhs, N * N * sizeof(double), cudaMemcpyDeviceToHost, stream);
	//cudaStreamDestroy(_stream);

	//cudaStreamDestroy(stream);

	cudaDeviceSynchronize();

}*/
void kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, _mySparse* mat, _mySparse* ret)
{
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(mat->__dmat, mat->__r, mat->__c);
	//int nn = _dmat2.cols();
	int nn = mat->_dmat.cols();
	int N = _dmat.cols();
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat3(ret->__dmat, N, nn);

	//Eigen::MatrixXd x(N,nn);
	/*if (ret->_dmat.cols() != nn || ret->_dmat.rows() != N)
	{
		ret->_dmat.resize(N, nn);
	}*/
	ret->__r = N;
	ret->__c = nn;
	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest(),0);
	auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	//_dmat3.setZero();
	ret->_dmat.setZero();
	int job = 0;

	cudaSetDevice(cuda->fastest());
	double* gpu_matrix = cuda->work_M(cuda->fastest());
	cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	int* devInfo_on_gpu;
	cudaMalloc(&devInfo_on_gpu, sizeof(int));
	double* work;
	int work_size = 0;

	// --- CUDA CHOLESKY initialization
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);

	// --- CUDA POTRF execution	
	cudaMalloc(&work, work_size * sizeof(double));
	auto now = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
	int devInfo_on_cpu = 0;
	cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);

	if (0 != devInfo_on_cpu) {
		//_dmat3(0, 0) = 2;
		ret->_dmat(0, 0) = 2;
		return;
	}
	cudaMemcpy(_dmat.data(), gpu_matrix, N * N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(work);
	cudaFree(devInfo_on_gpu);
	bool exit = false;

	int S = nn / cuda->count() / 5;

#pragma omp parallel for
	for (int i = 0; i < cuda->count(); i++) {
		cudaSetDevice(i);
		auto solver = cuda->solver(i,0);
		auto blas = cuda->blas(i);
		double* _gpu_matrix = cuda->work_M(i);
		double* gpu_rhs = cuda->work_rhs(i);
		cudaMemcpy(_gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
		//cudaMemcpy(gpu_rhs, _dmat2.data(), N * nn * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(gpu_rhs, mat->_dmat.data(), N * nn * sizeof(double), cudaMemcpyHostToDevice);
		int* _devInfo_on_gpu = 0;
		int _devInfo_on_cpu = 0;
		cudaMalloc(&_devInfo_on_gpu, sizeof(int));
		while (true)
		{
			int nextjob = -1;
#pragma omp critical
			{
				nextjob = job;
				job += S;
			}
			if (nextjob >= nn)break;
			int start = nextjob;
			int end = nextjob + S;
			if (end > nn)end = nn;
			cusolverDnDpotrs(solver, CUBLAS_FILL_MODE_LOWER, N, end - start, _gpu_matrix, N, gpu_rhs + nextjob * N, N, _devInfo_on_gpu);
			cudaMemcpy(&_devInfo_on_cpu, _devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
			if (_devInfo_on_cpu != 0) {
				exit = true;
				break;
			}
			//cudaMemcpy(_dmat3.data() + nextjob * N, gpu_rhs + nextjob * N, sizeof(double) * N * (end - start), cudaMemcpyDeviceToHost);
			cudaMemcpy(ret->_dmat.data() + nextjob * N, gpu_rhs + nextjob * N, sizeof(double) * N * (end - start), cudaMemcpyDeviceToHost);
		}
		cudaFree(_devInfo_on_gpu);
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Dn" << duration.count() << "ms" << std::endl;

	if (exit)
	{
		//_dmat3(0, 0) = 4;
		ret->_dmat(0, 0) = 4;
		return;
	}
}

void kingghidorah::_mySparse::_solve0(double* rhs, int N,double *ret) {
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> LLT;
	LLT.compute(_mat[0]);
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> b(rhs, N);
	//Eigen::VectorXd x(_mat[0].rows());
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> x(ret, _mat[0].rows());
	x.setZero();
	x = LLT.solve(b);
}
Eigen::VectorXd kingghidorah::_mySparse::__solve0(double* rhs, int N) {
	Eigen::LLT<Eigen::MatrixXd> LLT;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);

	LLT.compute(_dmat);
	Eigen::Map<Eigen::VectorXd, Eigen::Aligned128> b(rhs, N);
	Eigen::VectorXd x(_dmat.rows());
	x.setZero();
	x = LLT.solve(b);
	return x;
}
Eigen::MatrixXd kingghidorah::_mySparse::inv() {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);

	lu.compute(_dmat);
	Eigen::MatrixXd I(_dmat.rows(), _dmat.cols());
	I.setIdentity();
	return lu.solve(I);
}

Eigen::MatrixXd kingghidorah::_mySparse::solve0(_mySparse* rhs)
{
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;

	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(rhs->__dmat, rhs->__r, rhs->__c);

	//_mat.makeCompressed();
	lu.compute(_dmat);
	//Eigen::MatrixXd _x(_dmat.rows(), _dmat2.cols());
	//_x = lu.solve(_dmat2);
	Eigen::MatrixXd _x(_dmat.rows(), rhs->_dmat.cols());
	_x = lu.solve(rhs->_dmat);

	return _x;// .sparseView(0.000000000000001, 1.0);
}


void kingghidorah::_mySparse::minus(_mySparse* m) {
	this->_freeze();
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat2(m->__dmat, m->__r, m->__c);

	//_dmat = _dmat - _dmat2;
	_dmat = _dmat - m->_dmat;
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
}
void kingghidorah::_mySparse::clearcoeff() {
	for (int ii = 0; ii < _nt; ii++)
	{
		this->coeff[ii] = Eigen::VectorXd::Ones(this->_mat[ii].rows());
	}
}
Eigen::SparseMatrix<double> id;
void kingghidorah::_mySparse::addsmallidentity(double salt,bool sparse,bool dense) {

	if (dense)
	{
		//Eigen::Map<Eigen::MatrixXd, Eigen::Aligned128> _dmat(__dmat, __r, __c);

		id.conservativeResize(_dmat.rows(), _dmat.cols());
		id.setIdentity();
		_dmat += (id * salt);
	}
	if (sparse)

	{
		if (this->_mat.size() >= 1)
		{
			id.resize(this->_mat[0].rows(), this->_mat[0].cols());
			id.setIdentity();
			this->_mat[0] += (id * salt);
		}
	}
}
/*int main()
{
#ifdef VIENNACL_WITH_OPENCL
			std::cout << std::endl;
			std::cout << "----------------------------------------------" << std::endl;
			std::cout << "               Device Info" << std::endl;
			std::cout << "----------------------------------------------" << std::endl;
			std::cout << std::endl;
			std::cout << viennacl::ocl::current_device().info() << std::endl;
			std::cout << std::endl;
#endif
			int N = 1000;
			for (int t = 0; t < 10; t++) {
			viennacl::compressed_matrix<double, 1> dat(N,N);
			Eigen::SparseMatrix<double> M(N, N);
			for (int i = 0; i < N; i++)
			{
				for (int j = -8; j < 8; j++)
				{
					if (i + j >= 0 && i + j < N)
					{
						M.insert(i + j, i) = i + j;
					}
				}
			}

				viennacl::tools::timer timer;
				int time_previous = timer.get();

				dat = _copy(M);
				std::cout << "10:" << (timer.get() - time_previous) << "ms" << std::endl;
				time_previous = timer.get();

				viennacl::compressed_matrix<double, 1> ee = viennacl::linalg::prod(dat, dat);
				std::cout << "10:" << (timer.get() - time_previous) << "ms" << std::endl;
				time_previous = timer.get();
				M = _copy(ee);
				std::cout << "10:" << (timer.get() - time_previous) << "ms" << std::endl;
				time_previous = timer.get();
			}
			//viennacl::linalg::prod_impl()
			//std::cout << M<<std::endl;
		}*/