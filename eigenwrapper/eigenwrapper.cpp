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
#define STREAMCOUNT 2
bool __cuinit = false;
void kingghidorah::cuda::disable()
{
	__cuinit = false;
}
kingghidorah::cuda::cuda(int N) {
	I.resize(0, 0);
	omp_set_dynamic(false);
	//omp_set_num_threads(16);
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
		__cuinit = true;
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
		solver_handle[ii].resize(STREAMCOUNT);
		_streams[ii].resize(STREAMCOUNT);
		cudaSetDevice(ii);
		for(int kk=0;kk< STREAMCOUNT;kk++)
		cudaStreamCreate(&_streams[ii][kk]);
	}

	//std::cout << err << "yay" << std::endl;
	for (int i = 0; i < _count; i++)
	{
		for (int j = 0; j < STREAMCOUNT; j++)
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
		for(int j=0;j<STREAMCOUNT;j++)
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
			for(int j=0;j<STREAMCOUNT;j++)
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
			Eigen::VectorXd rhs(_N);
			for (int i = 0; i < _N; i++)rhs[i] = i;
			m.ofDat();
			m.clearcoeff();
			m._ofAtA(&m);
			Eigen::VectorXd ret(_N);
			m._solve0_gpu(this, &rhs,&ret,ii);
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
			for(int kk=0;kk< STREAMCOUNT;kk++)
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
			for(int j=0;j<STREAMCOUNT;j++)
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
	//delete _smat;
	if (___dmat != 0)
	{
		if (__cuinit)
		{
			cudaFreeHost(___dmat);
		}
		else {
			free(___dmat);
		}
		___dmat = 0;
		__r = 0;
		__c = 0;
	}
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
void kingghidorah::_mySparse::freeze(bool _do) {
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, this->rows(), this->cols());
	__r = this->rows(); __c = this->cols();
	_dmat.setZero();
	if(_do)
	for (int ii = 0; ii < _nt; ii++)
	{
		_dmat += coeff[ii].asDiagonal()*this->_mat[ii];
	}
	else
		for (int ii = 0; ii < _nt; ii++)
		{
			_dmat += this->_mat[ii];
		}
}
double kingghidorah::_mySparse::L2Norm(Eigen::VectorXd*a, Eigen::VectorXd* b) {
	return (*a).transpose() * this->_mat[0] * (*b);
}
Eigen::VectorXd kingghidorah::_mySparse::Vector(double* ptr1, int N1) {
	auto a = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr1, N1);
	return this->_mat[0] * a;
}
Eigen::VectorXd kingghidorah::_mySparse::Vector(Eigen::VectorXd* vec) {
	return this->_mat[0] * *vec;
}
void kingghidorah::_mySparse::Vector(Eigen::VectorXd* vec, Eigen::VectorXd *ret) {
	*ret=this->_mat[0] * *vec;
}
void kingghidorah::_mySparse::plus(_mySparse* m, double sc,bool dense,bool sparse) {
	if(sparse)
	this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
	if (dense)
	{
		Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r,__c);
		//this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
		_dmat += m->_mat[0] * sc;
	}
}
void kingghidorah::_mySparse::setmat(Eigen::SparseMatrix<double> mat, int ii) {
	this->_mat[ii] = mat;
}
void kingghidorah::_mySparse::setmat(const Eigen::MatrixXd& mat) {
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, mat.rows(), mat.cols());
	__r = mat.rows();
	__c = mat.cols();
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
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	return _dmat.data()[i];
}
double kingghidorah::_mySparse::_at(int i, int j) {
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	return _dmat(i, j);
}
int kingghidorah::_mySparse::cols() {
	return _mat[0].cols();
}
void kingghidorah::_mySparse::_resize(int n, int m) {
	Eigen::Map<Eigen::MatrixXd> map1(___dmat, __r, __c);
	Eigen::Map<Eigen::MatrixXd> map2(___dmat, n, m);

#pragma omp for schedule(static) ordered
	for (int i = 0; i < m; i++)
	{
#pragma omp ordered
		{
			map2.col(i) = map1.block(0, i, n, 1);
		}
	}

	__r = n;
	__c = m;
}
void kingghidorah::_mySparse::setmiddlecolum(Eigen::SparseMatrix<double> f, int start, int end) {
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	_dmat.middleCols(start, end - start) = f;
}
void kingghidorah::_mySparse::permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm)
{
	for (int ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii] * perm.transpose();
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm,bool sparse,bool dense)
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	int nn = __r;
	int numthreads = 0;
	numthreads = omp_get_max_threads();
	int S = nn / numthreads / 2;
	auto pt = perm.transpose();
	if(sparse)
	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0] = perm * (_mat[0]) * pt;
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
			_dmat.middleRows(start, end - start) = _dmat.middleRows(start, end - start) * pt;
			//perm.transpose().applyThisOnTheRight(_dmat.middleRows(start, end - start));
		}

#pragma omp parallel for
		for (int i = 0; i < nn; i += S)
		{
			int start = i;
			int end = i + S;
			if (end > nn)end = nn;
			_dmat.middleCols(start, end - start) = perm * _dmat.middleCols(start, end - start);
			//perm.applyThisOnTheLeft(_dmat.middleCols(start, end - start));
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
		_resize(M, M);
		//_dmat.conservativeResize(M, M);
	}
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm2)
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	int nn = __r;
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
	dat[0].reserve(dat[0].size() + N);
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
		coeff[ii] = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(_coeff[ii].data(), _coeff[ii].size());
	}
}
int kingghidorah::_mySparse::numBlocks()
{
	return this->dat.size();
}
//std::vector<Eigen::MatrixXd> e;
int kingghidorah::_mySparse::ofAtA(_mySparse* A,bool sparse)
{
	static std::vector<Eigen::SparseMatrix<double>> e;
	int nn = A->cols();
	int mt = omp_get_max_threads();
	_mt = mt*1;

	if (e.size() < _mt)
	{
		e.resize(_mt);
	}
	for (int i = 0; i < _mt; i++) {
		e[i].resize(nn, nn);
		e[i].setZero();
		e[i].reserve(nn * nn / 10);
	}
#pragma omp parallel for
	for (int _ii = 0; _ii < _mt; _ii++)
	{
		int S = 0;
		int E = 0;
		auto _e = e[_ii];
		S = _ii * _nt / _mt;
		E = (_ii + 1) * _nt / _mt;
		for (int ii = S; ii < E; ii++)
		{
			e[_ii] += (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
		}
	}
	if (sparse) {
		if (this->_mat.size() == 0)this->_mat.resize(1);
		this->_mat[0].resize(nn, nn);
		this->_mat[0].setZero();
		for (int i = 0; i < _mt; i++) {
			this->_mat[0] += e[i];
		}
		//this->_dmat = this->_mat[0];
	}
	else {
		//this->_dmat.setZero(nn, nn);
		if (__r == 0 || __c == 0)
		{
			if (__cuinit)
			{
				cudaMallocHost(&___dmat, sizeof(double) * nn * nn * 1.2);
			}
			else {
				___dmat=(double*)malloc(sizeof(double) * nn * nn * 1.2);
			}
		}
		__r = nn;
		__c = nn;
		Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, nn);
		_dmat.setZero();
		for (int i = 0; i < _mt; i++) {
			_dmat += e[i];
			//this->_dmat += e[i];
		}
		//this->_dmat = x;
	}
	//this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	return _nt;
}
void kingghidorah::_mySparse::_freeze() {
	//this->_dmat = this->_mat[0];
}
std::string kingghidorah::_mySparse::_ofAtA(_mySparse* A)
{
	if (__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&___dmat, sizeof(double) * A->cols() * A->cols());
		}
		else {
			___dmat=(double*)malloc(sizeof(double) * A->cols() * A->cols());
		}
	}
	__r = A->cols();
	__c = A->cols();
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, A->cols(), A->cols());
	std::stringstream ss;
#ifdef _DEBUG
	ss << _nt << std::endl;
	for (int ii = 0; ii < _nt; ii++)
	{
		ss << ii << "th" << std::endl;
		ss << A->_mat[ii].rows() << "," << A->_mat[ii].cols() << "," << coeff[ii].size() << std::endl;
	}
#endif

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

void kingghidorah::_mySparse::_ofAtB(_mySparse* B, _mySparse* C)
{
	this->freeze(true);
	B->freeze(false);

	int nn = this->_cols();
	int mm = B->_cols();
	if (C->__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&C->___dmat, sizeof(double) * nn * mm);
		}
		else {
			C->___dmat = (double*)malloc(sizeof(double) * nn * mm);
		}
	}
	C->__r = nn;
	C->__c = mm;
	//C->_dmat.resize(nn, mm);
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, mm);
	Eigen::Map<Eigen::MatrixXd> C_dmat(C->___dmat, nn, mm);
	C_dmat.setZero();

	int mt = omp_get_max_threads();

	int ss = mm / mt / 2;
	auto left = _dmat.transpose();
	auto right = B->_mat[0];

#pragma omp parallel for
	for (int ii = 0; ii < mm; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= mm)E = mm;
		C_dmat.middleCols(S, E - S) = left * right.middleCols(S, E - S);
	}
	if (C->_mat.size() == 0)C->_mat.resize(1);
	C->_mat[0] = C_dmat.sparseView(1.0, 0.0000000000001);
}
void kingghidorah::_mySparse::_ofBtAB(_mySparse* B,Eigen::VectorXd *b,_mySparse* C,Eigen::VectorXd* ret)
{

	static Eigen::MatrixXd D;

	int nn = B->_mat[0].cols();
	int kk = __c;


	if (C->__r == 0)
	{
		if (__cuinit) {
			cudaMallocHost(&C->___dmat, sizeof(double) * nn * nn);
		}
		else {
			C->___dmat = (double*)malloc(sizeof(double) * nn * nn);
		}
	}
	C->__r = nn;
	C->__c = nn;
	//C->_dmat.resize(nn, mm);
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::MatrixXd> C_dmat(C->___dmat, nn, nn);


	//C->_dmat.resize(nn, nn);
	C_dmat.setZero();

	int mt = omp_get_max_threads();

	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];

	D.resize(nn, kk);
	D.setZero();
	int ss = kk / mt / 2;

#pragma omp parallel for
	for (int ii = 0; ii < kk; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= kk)E = kk;
		D.middleCols(S, E - S) = left * mid.middleCols(S, E - S);
	}
	ss = nn / mt / 2;
#pragma omp parallel for
	for (int ii = 0; ii < nn; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= nn)E = nn;
		C_dmat.middleCols(S, E - S) = D * right.middleCols(S, E - S);
	}
	//Eigen::Map<Eigen::VectorXd> b(ptr, N);
	*ret = D * *b;
}

void kingghidorah::_mySparse::ofAtB(_mySparse* B, bool sparse)
{
	static std::vector<Eigen::SparseMatrix<double>> e;

	int nn = this->cols();
	int mm = B->cols();
	if (__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&___dmat, sizeof(double) * nn * mm);
		}
		else {
			___dmat = (double*)malloc(sizeof(double) * nn * mm);
		}
	}
	__r = nn;
	__c = mm;
	//C->_dmat.resize(nn, mm);
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, mm);

	//this->_dmat.resize(nn, mm);
	_dmat.setZero();
	int mt = omp_get_max_threads();

	_mt = mt*1;
	
	if (_mt > e.size())
		e.resize(_mt);
	for (int i = 0; i < _mt; i++) {
		e[i].resize(nn, mm);
		e[i].setZero();
		e[i].reserve(nn * mm / 10);
	}
#pragma omp parallel for
	for (int _ii = 0; _ii < _mt; _ii++)
	{
		int S = _ii * _nt / _mt;
		int E = (_ii + 1) * _nt / _mt;

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
		for (int i = 0; i < _mt; i++) {
			this->_mat[0] += e[i];
		}
	}
	else {
		for (int i = 0; i < _mt; i++) {
			_dmat += e[i];
		}
	}
	//this->_dmat = this->_mat[0];
}
void kingghidorah::_mySparse::Atb(double* ptr, int N,Eigen::VectorXd *c)
{
	c->resize(this->cols());
	c->setZero();
	int offset = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		int ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		*c += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		offset += ee;
	}
	//return ret;
}
Eigen::VectorXd kingghidorah::_mySparse::Atb(double* ptr, int N)
{
	Eigen::VectorXd ret(this->cols());
	ret.setZero();
	int offset = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		int ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		ret += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		offset += ee;
	}
	return ret;
}
Eigen::VectorXd kingghidorah::_mySparse::_Atb(double* ptr, int N)
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::VectorXd> b(ptr, N);
	return _dmat.transpose() * b;
}
void kingghidorah::_mySparse::merge()
{
	this->ofDat();
}
void kingghidorah::_mySparse::computeQR()
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::HouseholderQR<Eigen::MatrixXd> qr;
	qr.compute(_dmat);
}
void kingghidorah::_mySparse::computeLU()
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	lu.compute(_dmat);
}
void kingghidorah::_mySparse::computeLLT(Eigen::LLT<Eigen::MatrixXd>* _LLT)
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
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
void kingghidorah::_mySparse::solve0(Eigen::VectorXd* rhs,Eigen::VectorXd *ret) {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	lu.compute(_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_dmat.cols());
	//Eigen::VectorXd x(_dmat.cols());
	//x.setZero();
	*ret = lu.solve(*rhs);
}


void kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, Eigen::VectorXd *rhs, Eigen::VectorXd *ret,int device) {
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::VectorXd x = *ret;
	this->_freeze();
	int N = rhs->rows();
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
	x.resize(N);
	x.setZero();


	double* gpu_rhs = cuda->work_rhs(device);
	double* gpu_matrix = cuda->work_M(device);
	cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_rhs, rhs->data(), N * sizeof(double), cudaMemcpyHostToDevice);

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

	cudaMemcpy(ret->data(), gpu_rhs, sizeof(double) * N, cudaMemcpyDeviceToHost);

	//cudaFree(work);

	//cudaFree(devInfo_on_gpu);

	cudaDeviceSynchronize();
	//cudaStreamDestroy(stream);
	return;
}
Eigen::MatrixXd kingghidorah::_mySparse::_solve0(_myLLT* LLT, _mySparse* mat)
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//this function assumes that LLT decomposition has been done already
	Eigen::MatrixXd ret(this ->__c, mat->__c);
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
		ret.middleRows(start, end - start) = LLT->LLT->solve(_dmat.middleCols(start, end - start)).transpose();

	}
	return ret.transpose();
}

void kingghidorah::_mySparse::_solveI_gpu_mg(kingghidorah::cuda* cuda, _mySparse* ret)
{
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	int N = this->__c;
	//Eigen::MatrixXd x(N,nn);
	if (ret->__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N);
		}
	}
	ret->__r = N;
	ret->__c = N;
	Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	//ret->_dmat.resize(N, N);
	if (!cuda->valid())return;
	ret_dmat.setZero();
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


	//auto status = cusolverMgCreateDeviceGrid(&gridA, 1, nbGpus, deviceList, mapping);
	//assert(CUSOLVER_STATUS_SUCCESS == status);


	//assert(CUSOLVER_STATUS_SUCCESS == status);



	
	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		N,

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
		ret_dmat(0, 0) = info;
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

		ret_dmat.data(),
		lda
		);


#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			ret_dmat(j, i) = ret_dmat(i, j);
		}
	}


	//if (NULL != array_d_A) free(array_d_A);
	//if (NULL != array_d_B) free(array_d_B);
	//if (NULL != array_d_work) free(array_d_work);

}

void kingghidorah::_mySparse::_solveI(_mySparse* ret)
{
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);	
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
	llt.compute(this->_mat[0]);
	int nn = this->_mat[0].rows();
	if (ret->__r == 0)
	{
		if (__cuinit)
		{			
			cudaMallocHost(&ret->___dmat, sizeof(double) * nn * nn);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * nn * nn);
		}
	}
	ret->__r = nn;
	ret->__c = nn;
	Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, nn, nn);
	//ret->_dmat.resize(nn, nn);
	if (I.rows() != nn || I.cols() != nn)
	{
		I.resize(nn, nn);
		I.setIdentity();
	}
	int mt = omp_get_max_threads();
	int ee = mt * 2;
#pragma omp parallel for
	for (int i = 0; i < ee; i++)
	{
		int S = i * nn / ee;
		int E = (i + 1) * nn / ee;
		ret_dmat.middleCols(S, E - S) = llt.solve(I.middleCols(S, E - S));
	}
}

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
	int N = __c;
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

	if (ret->__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N);
		}
	}
	ret->__r = N;
	ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int nn = cuda->count();
	initidentiy(cuda, N);
	Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret_dmat.setZero();


	

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

	
		bool exit = false;
#pragma omp parallel for
	for (int kk = 0; kk < STREAMCOUNT; kk++)
	{
		int S = kk * N / STREAMCOUNT;
		int E = (kk + 1) * N / STREAMCOUNT;
		cudaStream_t _stream = cuda->__streams(cuda->fastest(), kk);
		cusolverDnHandle_t _solver = cuda->solver(cuda->fastest(), kk);
		cusolverDnSetStream(_solver, _stream);

		cusolverDnDpotrs(_solver, CUBLAS_FILL_MODE_LOWER, N, E - S, gpu_matrix, N, gpu_rhs + S * N, N, devInfo_on_gpu);
		cudaMemcpyAsync(ret_dmat.data() + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
	}

	cudaDeviceSynchronize();
	return sss.str();

}
std::string kingghidorah::_mySparse::_solveI_gpu_omp(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();
	std::stringstream sss;
	int N = this->__c;

	if (ret->__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N);
		}
	}
	ret->__r = N;
	ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int nn = cuda->count();
	initidentiy(cuda, N);
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret_dmat.setZero();


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
		for (int ss = 0; ss < STREAMCOUNT; ss++)
		{
			cudaSetDevice(cuda->fastest());
			int S = ss * N / STREAMCOUNT;
			int E = (ss + 1) * N / STREAMCOUNT;
			cudaStream_t _streamX = cuda->__streams(cuda->fastest(), ss);
			cudaMemcpyAsync(ret_dmat.data() + S * N, gpu_matrix + S * N, sizeof(double) * N * (E - S), cudaMemcpyDeviceToHost, _streamX);
			cudaStreamSynchronize(_streamX);
			//cudaDeviceSynchronize();

#pragma omp parallel for
			for (int ii = 0; ii < nn; ii++)
			{
				if (ii != cuda->fastest())
				{
					cudaSetDevice(ii);
					cudaStream_t _stream = cuda->__streams(ii, ss);
					cudaMemcpyAsync(cuda->work_M(ii) + S * N, ret_dmat.data() + S * N, sizeof(double) * N * (E - S), cudaMemcpyHostToDevice, _stream);
				}
			}
		}
		//cudaStreamDestroy(_stream);
	}

	int job = 0;
	int ss = N / cuda->count() / STREAMCOUNT/4;
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
			for (int kk = 0; kk < STREAMCOUNT; kk++)
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
				cudaMemcpyAsync(ret_dmat.data() + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
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
}

void kingghidorah::_mySparse::_solveI_gpu(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();

	int N = __c;
	if (ret->__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N);
		}
	}
	ret->__r = N;
	ret->__c = N;
	//ret->_dmat.resize(N, N);
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, N, N);
	Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret_dmat.setZero();
	//initidentiy(cuda, N);
	if (!cuda->valid())return;
	
		int ii = cuda->fastest();
		auto solver = cuda->solver(ii,0);
		//cudaStream_t stream = streams[ii];
		//cusolverDnSetStream(solver, stream);
		cudaSetDevice(ii);
		double* m_gpu = cuda->work_M(ii);

		cudaMemcpyAsync(m_gpu, _dmat.data(), sizeof(double) * N * N, cudaMemcpyHostToDevice,cuda->__streams(cuda->fastest(),0));
		double* work;
		int work_size;
		int work_size1 = 0;
		int work_size2 = 0;
		cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size1);
		cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size2);
		work_size = N * N;// std::max(work_size1, work_size2);
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

		cudaMemcpy(ret_dmat.data(), m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();

		ret_dmat.triangularView<Eigen::Upper>() = ret_dmat.triangularView<Eigen::Lower>().transpose();
/*#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			ret_dmat(j, i) = ret_dmat(i, j);
		}
	}*/
	
}

void kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, _mySparse* mat, _mySparse* ret)
{
	int nn = mat->__c;
	int N = this->__c;
	//Eigen::MatrixXd x(N,nn);

	if (ret->__r == 0)
	{
		if (__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * nn);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * nn);
		}
	}
	ret->__c = nn;
	ret->__r = N;
	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest(),0);
	auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, nn);
	Eigen::Map<Eigen::MatrixXd> mat_dmat(mat->___dmat, mat->__r, mat->__c);
	ret_dmat.setZero();
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
		ret_dmat(0, 0) = 2;
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
		cudaMemcpy(gpu_rhs, mat_dmat.data(), N * nn * sizeof(double), cudaMemcpyHostToDevice);
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
			cudaMemcpy(ret_dmat.data() + nextjob * N, gpu_rhs + nextjob * N, sizeof(double) * N * (end - start), cudaMemcpyDeviceToHost);
		}
		cudaFree(_devInfo_on_gpu);
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Dn" << duration.count() << "ms" << std::endl;

	if (exit)
	{
		ret_dmat(0, 0) = 4;
		return;
	}
}

void kingghidorah::_mySparse::_solve0(Eigen::VectorXd* rhs,Eigen::VectorXd *ret) {
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> LLT;
	LLT.compute(_mat[0]);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_mat[0].cols());
	ret->setZero();
	//Eigen::VectorXd x(_mat[0].rows());
	//x.setZero();
	*ret = LLT.solve(*rhs);
	//return x;
}
void kingghidorah::_mySparse::__solve0(Eigen::VectorXd* rhs, Eigen::VectorXd *ret) {
	Eigen::LLT<Eigen::MatrixXd> LLT;
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

	LLT.compute(_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->resize(_dmat.cols());
	ret->setZero();
	*ret = LLT.solve(*rhs);
}
Eigen::MatrixXd kingghidorah::_mySparse::inv() {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	lu.compute(_dmat);
	Eigen::MatrixXd I(_dmat.rows(), _dmat.cols());
	I.setIdentity();
	return lu.solve(I);
}

Eigen::MatrixXd kingghidorah::_mySparse::solve0(_mySparse* rhs)
{
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	Eigen::Map<Eigen::MatrixXd> rhs_dmat(rhs->___dmat, rhs->__r, rhs->__c);

	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//_mat.makeCompressed();
	lu.compute(_dmat);
	Eigen::MatrixXd _x(__c, rhs->__c);
	_x = lu.solve(rhs_dmat);

	return _x;// .sparseView(0.000000000000001, 1.0);
}


void kingghidorah::_mySparse::minus(_mySparse* m) {
	this->_freeze();
	Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::MatrixXd> m_dmat(m->___dmat, m->__r, m->__c);

	_dmat = _dmat - m_dmat;
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
		id.resize(__r,__c);
		id.setIdentity();
		Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

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