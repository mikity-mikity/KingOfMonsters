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
#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT
#define EIGEN_USE_LAPACKE
int previdentiyN = 0;
//std::vector<cudaStream_t> streams;
Eigen::MatrixXd I;
#define STREAMCOUNT 8
bool __cuinit = false;
//static std::vector<Eigen::SparseMatrix<double>> e;
std::vector<std::vector<cusparseHandle_t>> sp_handle;
std::map<std::tuple<kingghidorah::_mySparse*,int>, kingghidorah::spgemm> dict;
void kingghidorah::cuda::disable()
{
	__cuinit = false;
}
kingghidorah::cuda::cuda(int N) {
	dict.clear();
	Eigen::initParallel();
	Eigen::setNbThreads(omp_get_max_threads());
	//e.shrink_to_fit();
	//e.clear();
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
		_canpeer = enablePeerAccess(_count, _deviceList);
	}
	else {
		_canpeer = false;
	}
	solver_handle.resize(_count);
	sp_handle.resize(_count);
	_streams.resize(_count);
	for (int ii = 0; ii < _count; ii++)
	{
		solver_handle[ii].resize(STREAMCOUNT);
		sp_handle[ii].resize(STREAMCOUNT);
		_streams[ii].resize(STREAMCOUNT);
		cudaSetDevice(ii);
		for (int kk = 0; kk < STREAMCOUNT; kk++)
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
	for (int ii = 0; ii < count() * STREAMCOUNT; ii++)
	{
		cudaSetDevice(ii%count());
		cusparseCreate(&sp_handle[ii%count()][ii/count()]);
	}
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		cusolverStatus_t status;
		for (int j = 0; j < STREAMCOUNT; j++)
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
			for (int j = 0; j < STREAMCOUNT; j++)
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
	//cusparseCreate(&sp_handle);
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
	/*for (int kk = 0; kk < 3; kk++)
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
			//m._ofAtA(&m);
			Eigen::VectorXd ret(_N);
			//m._solve0_gpu(this, &rhs, &ret, ii);
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			speed[ii] = duration.count();
		}
	}*/
	_fastest = 0;// std::distance(speed.begin(), std::min_element(speed.begin(), speed.end()));
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
		/*cusolverMgCreate(&mg_solver);
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
		*/
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
			for (int kk = 0; kk < STREAMCOUNT; kk++)
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
			for (int j = 0; j < STREAMCOUNT; j++)
				if (solver_handle[i][j] != 0)
				{
					cusolverDnDestroy(solver_handle[i][j]);
					solver_handle[i][j] = 0;
				}
			for (int j = 0; j < STREAMCOUNT; j++)
			{
				if (sp_handle[i][j] != 0)
				{
					cusparseDestroy(sp_handle[i][j]);
					sp_handle[i][j] = 0;
				}
			}
			if (cublas_handle != 0)
			{
				cublasDestroy(cublas_handle[i]);
			}
			cublas_handle[i] = 0;
		}
		/*if (mg_solver != 0)
		{
			cusolverMgDestroy(mg_solver);
		}
		mg_solver = 0;*/
		//if (sp_handle != 0)cusparseDestroy(sp_handle);
		//sp_handle = 0;
		cudaDeviceReset();
	}
	initialized = false;
}
/*
cusolverMgHandle_t kingghidorah::cuda::mgsolver() {
	return mg_solver;
}*/
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
double* kingghidorah::cuda::work(int N, int device, cudaStream_t stream) {
	if (N < work_size[device])
	{
		//do nothing
		return __work[device];
	}
	else {
		if (__work[device] != 0)
		{
			//cudaSetDevice(device);
			cudaFreeAsync(__work[device], stream);
			cudaMallocAsync(&__work[device], N * sizeof(double), stream);
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
cusolverDnHandle_t& kingghidorah::cuda::solver(int ii, int kk) {
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
	dat.reserve(1000);
	coeff.reserve(1000);
	_coeff.reserve(1000);
	_mat.reserve(1000);
	_mt = omp_get_max_threads();
	//prevmat.resize(1, 1);
	//e = new Eigen::SparseMatrix<double>[200];
	//e2 = new Eigen::SparseMatrix<double>[200];
	//_smat.resize(1);
	//_smat = new Eigen::SparseMatrix<double>();
}
kingghidorah::_mySparse::~_mySparse()
{
	
	/*if (dict.contains(this))
	{
		spgemm __spgemm_dat=dict[this];

		cusparseSpGEMM_destroyDescr(__spgemm_dat.spgemmDesc);
		cusparseDestroySpMat(__spgemm_dat.matA);
		cusparseDestroySpMat(__spgemm_dat.matB);
		cusparseDestroySpMat(__spgemm_dat.matC);

		//--------------------------------------------------------------------------
		// device memory deallocation
		cudaFree(__spgemm_dat.dBuffer1);
		cudaFree(__spgemm_dat.dBuffer2);
		cudaFree(__spgemm_dat.dA_csrOffsets);
		cudaFree(__spgemm_dat.dA_columns);
		cudaFree(__spgemm_dat.dA_values);
		cudaFree(__spgemm_dat.dB_csrOffsets);
		cudaFree(__spgemm_dat.dB_columns);
		cudaFree(__spgemm_dat.dB_values);
		cudaFree(__spgemm_dat.dC_csrOffsets);
		cudaFree(__spgemm_dat.dC_columns);
		cudaFree(__spgemm_dat.dC_values);
		dict.erase(this);
	}*/
	//if (e != 0)delete[] e;
	//e = 0;
	//if (e2 != 0)delete[] e2;
	//e2 = 0;

	//delete _smat;
	/*if (___dmat != 0)
	{
		if (false)//__cuinit)
		{
			cudaFreeHost(___dmat);
		}
		else {
			free(___dmat);
		}
		___dmat = 0;
		__r = 0;
		__c = 0;
	}*/
	//__r = 0;
	//__c = 0;
	//_dmat.resize(0, 0);
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
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, this->rows(), this->cols());
	//__r = this->rows(); __c = this->cols();
	_dmat.setZero(this->rows(), this->cols());
	if (_do)
		for (int ii = 0; ii < _nt; ii++)
		{
			if(coeff[ii].size()!=0)
			_dmat += coeff[ii].asDiagonal() * this->_mat[ii];
		}
	else
		for (int ii = 0; ii < _nt; ii++)
		{
			if(this->_mat[ii].rows()!=0)
			_dmat += this->_mat[ii];
		}
}
double kingghidorah::_mySparse::L2Norm(Eigen::VectorXd* a, Eigen::VectorXd* b) {
	return (*a).transpose() * this->_mat[0] * (*b);
}
Eigen::VectorXd kingghidorah::_mySparse::Vector(double* ptr1, int N1) {
	auto a = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr1, N1);
	return this->_mat[0] * a;
}
Eigen::VectorXd kingghidorah::_mySparse::Vector(Eigen::VectorXd* vec) {
	return this->_mat[0] * *vec;
}
void kingghidorah::_mySparse::Vector(Eigen::VectorXd* vec, Eigen::VectorXd* ret) {
	*ret = this->_mat[0] * *vec;
}
void kingghidorah::_mySparse::plus(_mySparse* m, double sc, bool dense, bool sparse) {
	if (sparse)
		this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
	if (dense)
	{
		//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
		//this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
		_dmat += m->_mat[0] * sc;
	}
}
void kingghidorah::_mySparse::setmat(Eigen::SparseMatrix<double, Eigen::RowMajor> &mat, int ii) {
	this->_mat[ii] = mat;
}
void kingghidorah::_mySparse::setmat(const Eigen::MatrixXd& mat) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, mat.rows(), mat.cols());
	//__r = mat.rows();
	//__c = mat.cols();
	_dmat.resize(mat.rows(), mat.cols());
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
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	return _dmat.data()[i];
}
double kingghidorah::_mySparse::_at(int i, int j) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	return _dmat(i, j);
}
int kingghidorah::_mySparse::cols() {
	return _mat[0].cols();
}
void kingghidorah::_mySparse::_resize(int n, int m) {

	//Eigen::Map<Eigen::MatrixXd> map1(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> map2(___dmat, n, m);

/*#pragma omp for schedule(static) ordered
	for (int i = 0; i < m; i++)
	{
#pragma omp ordered
		{
			_dmat.col(i) = _dmat.block(0, i, n, 1);
		}
	}
	*/
	_dmat.conservativeResize(n, m);
	//__r = n;
	//__c = m;
}
void kingghidorah::_mySparse::setmiddlecolum(Eigen::SparseMatrix<double, Eigen::RowMajor> &f, int start, int end) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	_dmat.middleCols(start, end - start) = f;
}
void kingghidorah::_mySparse::permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm)
{
	for (int ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii] * perm.transpose();
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm, bool sparse, bool dense)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> _tmp(___dmat+__r*__c, __r, __c);
	int nn = _dmat.rows();// __r;
	int numthreads = 0;
	//_mt=omp_get_max_threads();
	int S = nn / _mt / 2;
	auto pt = perm.transpose();
	if (sparse||dense)
		if (_mat.size() >= 1)
		{
			//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
			{
				_mat[0] = perm * (_mat[0]) * pt;
			}
		}
	if (false/*dense*/)
	{
		//_tmp.noalias() = _dmat * pt;
		//_dmat.noalias() = perm * _tmp;
		//_dmat.applyOnTheLeft(perm);
		//_dmat.applyOnTheRight(perm.transpose());
		//prrm.transpose().applyThisOnTheRight(_dmat);
		//perm.applyThisOnTheLeft(_dmat);

#pragma omp parallel for
		for (int i = 0; i < nn; i += S)
		{
			int start = i;
			int end = i + S;
			if (end > nn)end = nn;
			_dmat.middleRows(start, end - start).applyOnTheRight(pt);// = _dmat.middleRows(start, end - start) * pt;
			//_tmp.middleRows(start, end - start).noalias()= _dmat.middleRows(start, end - start) * pt;
			//perm.transpose().applyThisOnTheRight(_dmat.middleRows(start, end - start));
		}

#pragma omp parallel for
		for (int i = 0; i < nn; i += S)
		{
			int start = i;
			int end = i + S;
			if (end > nn)end = nn;
			//_dmat.middleCols(start, end - start).noalias() = perm * _dmat.middleCols(start, end - start);
			_dmat.middleCols(start, end - start).applyOnTheLeft(perm);// = perm * _dmat.middleCols(start, end - start);
			//perm.applyThisOnTheLeft(_dmat.middleCols(start, end - start));
		}
	}

}
void kingghidorah::_mySparse::shrink(int M)
{
	for (int ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii].leftCols(M);
}
void kingghidorah::_mySparse::_shrink(int M, bool sparse, bool dense)
{
	if (true/*sparse*/)
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
		//__r = M;
		//__c = M;
		//Eigen::Map<Eigen::MatrixXd> map1(___dmat, __r, __c);
		_dmat.resize(M, M);
		_dmat = _mat[0];
		//_resize(M, M);
		//_dmat.conservativeResize(M, M);
	}
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic>& perm2)
{
	int nn = _dmat.rows();// __r;
	//int numthreads = omp_get_max_threads();

	int S = nn / _mt / 2;
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
	return _dmat.rows();// __r;
}
int kingghidorah::_mySparse::_cols() {
	return _dmat.cols();// __c;
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
	if (_nt == 0)_nt = 1;
	this->dat.resize(_nt);
	this->_coeff.resize(_nt);

	if (_nt == 1)
	{
		this->dat[0].clear();
		this->_coeff[0].clear();
	}


	if (this->coeff.size() != _nt)this->coeff.resize(_nt);
	if (this->_mat.size() != _nt)this->_mat.resize(_nt);

//#pragma omp parallel for
	//for (int ii = 0; ii < _nt; ii++)
	{
		//this->coeff[ii].setZero();
		//this->_mat[ii].setZero();
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
	eigen_assert(_coeff[0].size() == ii );
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
		//eigen_assert(*data != 0);
		//eigen_assert(*data < 10000);
		//eigen_assert(*data > -10000);
		//eigen_assert(*ptr < _mat[0].cols());
		//eigen_assert(*ptr >= 0);
		//if(*data!=0 && *data<100&&*data>-100)
			dat[0].push_back(Eigen::Triplet<double>(ii, *ptr, (*data)));
		ptr++;
		data++;
	}
	if (add)
	{
		eigen_assert(ii == _coeff[0].size());
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
	//_mat.clear();
	//_mat.shrink_to_fit();
	_mat.resize(_nt);
	coeff.resize(_nt);
#pragma omp parallel for
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
		this->_mat[ii].resize(mat->_mat[ii].rows(), mat->_mat[ii].cols());
		this->_mat[ii].reserve(mat->_mat[ii].nonZeros());
		this->_mat[ii] = mat->_mat[ii];// M;
		
		this->coeff[ii].resize(mat->coeff[ii].size());
		this->_coeff[ii].resize(mat->coeff[ii].size());

		this->coeff[ii] = mat->coeff[ii];
		this->_coeff[ii] = mat->_coeff[ii];

		/*for (int k = 0; k < mat->coeff[ii].size(); k++)
		{
			this->coeff[ii](k) = mat->coeff[ii](k);
			this->_coeff[ii][k] = mat->_coeff[ii][k];
		}*/
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
		_mat[ii].setZero();
		_mat[ii].reserve(dat[ii].size());
		if (dat[ii].size() > 0)
		{
			_mat[ii].setFromTriplets(dat[ii].begin(), dat[ii].end());
			_mat[ii].makeCompressed();
		}
		else {
			_mat[ii].setZero();
		}
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

std::string kingghidorah::_mySparse::ofAtA( _mySparse* A, bool sparse)
{
	static std::map<_mySparse*, Eigen::SparseMatrix<double, Eigen::RowMajor>> dict2;
	static std::map< Eigen::SparseMatrix<double, Eigen::RowMajor>*, std::vector<std::vector<int>>>  map;
	static std::vector<std::vector<int>> index;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e2;
	auto ss = std::stringstream();
	auto now = high_resolution_clock::now();
	int nn = this->cols();
	int mm = this->cols();
	int _mt = omp_get_max_threads();
	omp_set_num_threads(_mt);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	int __mt = _mt;// std::min(_nt / 10, _mt * 10);
	if (e.size() < __mt)
	{
		e.resize(__mt);
		e2.resize(__mt);
	}
	Eigen::initParallel();
	Eigen::setNbThreads(1);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	ss << _nt << ":nt:" << std::endl;
	ss << _mt << ":_mt:" << std::endl;
	Eigen::initParallel();
	int job = 0;
	int sss = _nt / __mt / 8;
	int __nt2 = _nt;
	if (sss == 0)sss = 1;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* prevmat;
	if (dict2.contains(this)) {
		prevmat = &dict2[this];
	}
	else {
		Eigen::SparseMatrix<double, Eigen::RowMajor> _prevmat(nn, nn);
		dict2[this] = _prevmat;
		prevmat = &dict2[this];
		prevmat->resize(nn, nn);
		prevmat->reserve(nn * nn / 10);
	}
	//prevmat->makeCompressed();
	std::vector<std::vector<int>>* _map = 0;
	if (map.contains(prevmat))
	{
		_map = &map[prevmat];
	}
#pragma omp parallel for
	for (int i = 0; i < __mt; i++) {
		e[i].resize(nn, mm);
		if (e[i].nonZeros() != prevmat->nonZeros())
		{
			e[i] = *prevmat;
			e[i].makeCompressed();
		}
		memset(e[i].valuePtr(), 0, sizeof(double) * prevmat->nonZeros());
		e2[i].resize(nn, mm);
		e2[i].reserve(nn * mm / 100);
	}
	index.resize(__mt);
#pragma omp parallel for schedule(static)
	for (int _ii = 0; _ii < __mt; _ii++)
	{
		int S = 0;
		int E = 0;
		//int K = 0;
		//auto _e = e[_ii];
		for (int tt = 0; tt < 40; tt++)
		{
#pragma omp critical
			{
				S = job;
				E = job + sss;
				if (E > _nt)E = _nt;
				job = E;
			}
			if (S >= _nt)break;
			for (int ii = S; ii < E; ii++)
			{
				/*Eigen::SparseMatrix<double> __coeff(coeff[ii].size(), coeff[ii].size());
				for (int k = 0; k < coeff[ii].size(); k++)
				{
					__coeff.insert(k, k) = coeff[ii](k);
				}*/
				//e[_ii] = 
//#pragma omp critical
				{

					e2[_ii] = this->_mat[ii].transpose() * coeff[ii].asDiagonal() * this->_mat[ii];
					index[_ii].resize(e2[_ii].nonZeros());
					if (_map == 0)
					{
					}
					else {
						int count = 0;
						for (int k = 0; k < nn; ++k) {
							for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(e2[_ii], k); it; ++it) {
								//e[0].coeffRef(it.row(), it.col()) += it.value();
								index[_ii][count] = (*_map)[it.row()][it.col()];
								count++;
							}
						}
						//e[_ii].makeCompressed();
					}

//#pragma omp critical
					{
						if (_map == 0)
						{
							e[_ii] += e2[_ii];
						}
						else {
							int* ptr = &index[_ii][0];
							for (int k = 0; k < nn; ++k) {
								for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(e2[_ii], k); it; ++it) {
									//e[0].coeffRef(it.row(), it.col()) += it.value();
									*(e[_ii].valuePtr() + *ptr) += it.value();
									ptr++;
								}
							}
							//e[_ii].makeCompressed();
						}
					}
				}
			}
		}
		//e[_ii].makeCompressed();
	}
	//this->_mat[0] = *prevmat;
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	Eigen::setNbThreads(_mt);

	for (int tt = 0; tt < 40; tt++)
	{
#pragma omp parallel for
		for (int i = 0; i < __mt; i += 2)
		{
			if (i + 1 < __mt) {
				if (_map == 0 || e[i].nonZeros() != e[i + 1].nonZeros())
				{
					e[i] += e[i + 1];
				}
				else {
					Eigen::Map<Eigen::VectorXd> map1(e[i].valuePtr(), e[i].nonZeros());
					Eigen::Map<Eigen::VectorXd> map2(e[i + 1].valuePtr(), e[i].nonZeros());
					map1 += map2;
				}
			}
		}
		int _ct = 0;
#pragma omp parallel for ordered schedule(dynamic)
		for (int i = 0; i < __mt; i += 2)
		{
#pragma omp ordered
			if (_map == 0 || e[i].nonZeros() != e[i / 2].nonZeros())
			{
				e[i / 2] = e[i];
			}
			else
			{
				memcpy(e[i / 2].valuePtr(), e[i].valuePtr(), sizeof(double) * e[i].nonZeros());
			}
#pragma omp atomic
			_ct++;
		}
		__mt = _ct;
		if (__mt == 1)break;
	}
	if (true) {
		if (this->_mat.size() == 0)this->_mat.resize(1);
		this->_mat[0].resize(nn, nn);
		this->_mat[0].reserve(nn * nn / 20);
		this->_mat[0] = e[0];
		//for (int i = 1; i < __mt; i++) {
		//	this->_mat[0] += e[i];
		//}
		this->_mat[0].makeCompressed();
		if (_map == 0 || prevmat->nonZeros() != this->_mat[0].nonZeros())
		{
			*prevmat = this->_mat[0];
			prevmat->makeCompressed();
			dict2[this] = *prevmat;
			//build map
			std::vector<std::vector<int>> __map;
			__map.resize(this->_mat[0].rows());
			for (int i = 0; i < this->_mat[0].rows(); i++)
			{
				__map[i].resize(this->_mat[0].cols());
			}
			for (int k = 0; k < prevmat->outerSize(); ++k) {
				for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(*prevmat, k); it; ++it) {
					int S = prevmat->outerIndexPtr()[k];
					int E = prevmat->outerIndexPtr()[k + 1];

					for (int tt = S; tt < E; tt++)
					{
						if (prevmat->innerIndexPtr()[tt] == it.col())
							__map[it.row()][it.col()] = tt;
					}
				}
			}
			map[prevmat] = __map;
		}
	}


	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	if (!sparse) {
		//this->_dmat.setZero(nn, nn);
		//this->_tmp.setZero(nn, nn);
		/*if (__r == 0 || __c == 0)
		{
			if (false)//__cuinit)
			{
				cudaMallocHost(&___dmat, sizeof(double) * nn * nn * 2);
			}
			else {
				___dmat = (double*)malloc(sizeof(double) * nn * nn * 2);
			}
		}*/
		//__r = nn;
		//__c = nn;
		/*Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, nn);
		_dmat = e[0];
		for (int i = 1; i < _mt; i ++) {
			_dmat += e[i];
		}*/
		//this->_dmat = x;
	}
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	//this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	ss << "sum:" << e[0].sum() << std::endl;
	return ss.str();
}




std::string kingghidorah::_mySparse::ofAtA_gpu(cuda* _cuda, _mySparse* A, bool sparse)
{
	static std::vector<Eigen::SparseMatrix<double,Eigen::RowMajor>> e;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> tmp;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> tmp2;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> res;
	auto ss = std::stringstream();
	auto now = high_resolution_clock::now();
	this->join();

	this->_nt = 1;
	if (e.size() < _mt)e.resize(_mt);
	tmp.resize(_mt);
	tmp2.resize(_mt);
	res.resize(_mt);
	int nn = this->cols();
	//cuInit(0);
	int _mt = omp_get_max_threads();
	for (int i = 0; i < _mt; i++)
	{
		e[i].resize(nn, nn);
		e[i].setZero();
		e[i].reserve(nn * nn / 20);
		e[i].makeCompressed();
		res[i].resize(nn, nn);
		res[i].setZero();
		res[i].reserve(nn * nn / 20);
		res[i].makeCompressed();
	}
	//std::vector <cusolverSpHandle_t> solversp_handles(_mt);
	//std::vector <CUstream> __streams(_mt);
	//cuInit(0);
	cudaSetDevice(0);

	//cusparseCreate(&sp_handle);
		
	//static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e;
	auto ALG = CUSPARSE_SPGEMM_CSR_ALG_DETERMINITIC;
	int __mt = 1;// STREAMCOUNT* _cuda->count();
	omp_set_num_threads(__mt);
		cudaSetDevice(0);
		int _ss = 20;// _nt / __mt / 4;
		//if (_ss == 0)_ss = 1;
		int job = 0;


		kingghidorah::spgemm ____spgemm_dat;
		auto key = std::tuple<kingghidorah::_mySparse*,int>(this,0);
		if (dict.contains(key))
		{
			____spgemm_dat = dict[key];
		}
		else {
			____spgemm_dat.initialized = false;
		}
		//if (!____spgemm_dat.initialized)
		//{

		
		int max_rows = 0;
		int max_cols = 0;
		int max_nnz = 0;
		for (int ii = 0; ii < _nt; ii++)
		{
			int num_rows = this->_mat[ii].rows();
			int num_cols = this->_mat[ii].cols();
			int nnz = this->_mat[ii].nonZeros();
			max_rows = std::max(max_rows, num_rows);
			max_cols = std::max(max_cols, num_cols);
			max_nnz = std::max(max_nnz, nnz);

		}
		if (____spgemm_dat.initialized == false)
		{
			for (int _ii = 0; _ii < __mt; _ii++)
			{
				int device = _ii % _cuda->count();
				cudaSetDevice(device);

				kingghidorah::spgemm __spgemm_dat;
				auto __key = std::tuple<kingghidorah::_mySparse*, int>(this, _ii);

				__spgemm_dat.A_num_rows = max_cols;
				__spgemm_dat.A_num_cols = max_rows;
				__spgemm_dat.A_nnz = max_nnz;

				__spgemm_dat.B_num_rows = max_rows;
				__spgemm_dat.B_num_cols = max_cols;
				__spgemm_dat.B_nnz = max_nnz;

				__spgemm_dat.C_num_rows = max_cols;
				__spgemm_dat.C_num_cols = max_cols;
				__spgemm_dat.C_nnz = max_cols * max_cols / 3+1000;
				auto err = cudaMalloc(&__spgemm_dat.dA_csrOffsets, sizeof(int) * (__spgemm_dat.A_num_rows + 1));
				err = cudaMalloc(&__spgemm_dat.dA_columns, sizeof(int) * __spgemm_dat.A_nnz);
				err = cudaMalloc(&__spgemm_dat.dA_values, sizeof(double) * __spgemm_dat.A_nnz);

				// allocate B
				err = cudaMalloc(&__spgemm_dat.dB_csrOffsets, (__spgemm_dat.B_num_rows + 1) * sizeof(int));
				err = cudaMalloc(&__spgemm_dat.dB_columns, __spgemm_dat.B_nnz * sizeof(int));
				err = cudaMalloc(&__spgemm_dat.dB_values, __spgemm_dat.B_nnz * sizeof(double));

				// allocate B
				err = cudaMalloc(&__spgemm_dat.dC_csrOffsets, (__spgemm_dat.C_num_rows + 1) * sizeof(int));
				err = cudaMalloc(&__spgemm_dat.dC_columns, __spgemm_dat.C_nnz * sizeof(int));
				err = cudaMalloc(&__spgemm_dat.dC_values, __spgemm_dat.C_nnz * sizeof(double));

				err = cudaMalloc((double**)&__spgemm_dat.dBuffer1, max_rows * max_cols *2 + 1000);
				err = cudaMalloc((double**)&__spgemm_dat.dBuffer2, max_rows * max_cols *2 + 1000);
				__spgemm_dat.bufferSize1 = max_rows * max_cols *2 + 1000;
				__spgemm_dat.bufferSize2 = max_rows * max_cols *2 + 1000;
				__spgemm_dat.initialized = true;
				dict[__key] = __spgemm_dat;
			}
		}
			
		//}

		cudaDeviceSynchronize();
		//for (int _tt = 0; _tt < _cuda->count(); _tt++)
		//{
#pragma omp parallel for
		for (int _ii = 0; _ii < __mt; _ii++)
			{
			int device =  _ii% _cuda->count();
			cudaSetDevice(device);
				kingghidorah::spgemm* __spgemm_dat;
				auto __key = std::tuple<kingghidorah::_mySparse*, int>(this, _ii);
				__spgemm_dat = &dict[__key];

				CUstream stream = _cuda->__streams(device, _ii/_cuda->count());
				cusparseHandle_t handle = sp_handle[device][_ii/ _cuda->count()];
				//CUstream stream;// = _cuda->__streams(device, _ii / _cuda->count());
				//cuStreamCreate(&stream,0);
				//cusparseHandle_t     handle;// = sp_handle;
				//cusparseCreate(&handle);
											//auto handle = _cuda->sp_handle[_ii % _cuda->count()][_ii / _cuda->count()];
				cusparseSetStream(handle, stream);
				
				
				double               alpha = 1.0;
				double               beta = 0.0;
				cusparseOperation_t opA = CUSPARSE_OPERATION_NON_TRANSPOSE;
				cusparseOperation_t opB = CUSPARSE_OPERATION_NON_TRANSPOSE;
				cudaDataType        computeType = CUDA_R_64F;
				cusparseSpMatDescr_t matA, matB, matC;
				cusparseSpGEMMDescr_t spgemmDesc;
				size_t bufferSize1 = 0, bufferSize2 = 0;
				int64_t C_num_rows1 = nn, C_num_cols1 = nn;

				//cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
	

				for (int kk = 0; kk < 100000; kk++)
				{
					
					//if (_ii == 0 )break;
					int S = 0;// _ii* _nt / __mt;
					int E = 0;// (_ii + 1)* _nt / __mt;

					//int S = 0;
					//int E = 0;
#pragma omp critical
					{
						S = job;
						E = S + _ss;
						if (E > _nt)E = _nt;
						job = E;
					}
					if (S >= _nt)break;
					




					for (int ii = S; ii < E; ii++)
					{
						int64_t C_nnz1 = 0;
						if (C_nnz1 == 0)
							cudaMemsetAsync(__spgemm_dat->dC_csrOffsets, 0, sizeof(int) * (nn + 1), stream);
						//int64_t* nnzTotalDevHostPtr = &C_nnz1;
						//csrgemm2Info_t info = NULL;
						//cusparseCreateCsrgemm2Info(&info);

						tmp[_ii] = (this->_mat[ii].transpose()) * coeff[ii].asDiagonal();
						tmp2[_ii] = (this->_mat[ii]);

						tmp2[_ii].makeCompressed();
						tmp[_ii].makeCompressed();
						int A_num_rows = tmp[_ii].rows();
						int A_num_cols = tmp[_ii].cols();
						int A_nnz = tmp[_ii].nonZeros();
						int B_num_rows = tmp2[_ii].rows();
						int B_num_cols = tmp2[_ii].cols();
						int B_nnz = tmp2[_ii].nonZeros();
						C_num_rows1 = A_num_rows;
						C_num_cols1 = B_num_cols;
						if (A_nnz * B_nnz == 0)continue;


						//std::cout << "A_nnz" << A_nnz << std::endl;
						//std::cout << "B_nnz" << B_nnz << std::endl;
						//--------------------------------------------------------------------------
						// Device memory management: Allocate and copy A, B


						// copy A
						if (A_nnz > 0)
						{
							auto err = cudaMemcpyAsync(__spgemm_dat->dA_csrOffsets, tmp[_ii].outerIndexPtr(), (A_num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);
							err = cudaMemcpyAsync(__spgemm_dat->dA_columns, tmp[_ii].innerIndexPtr(), A_nnz * sizeof(int), cudaMemcpyHostToDevice, stream);
							err = cudaMemcpyAsync(__spgemm_dat->dA_values, tmp[_ii].valuePtr(), A_nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
						}
						// copy B
						if (B_nnz > 0)
						{
							auto err = cudaMemcpyAsync(__spgemm_dat->dB_csrOffsets, tmp2[_ii].outerIndexPtr(), (B_num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);
							err = cudaMemcpyAsync(__spgemm_dat->dB_columns, tmp2[_ii].innerIndexPtr(), B_nnz * sizeof(int), cudaMemcpyHostToDevice, stream);
							err = cudaMemcpyAsync(__spgemm_dat->dB_values, tmp2[_ii].valuePtr(), B_nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
						}
						//--------------------------------------------------------------------------
						// CUSPARSE APIs



						// Create sparse matrix A in CSR format
						if (A_nnz == 0)
						{
							auto status = cusparseCreateCsr(&matA, A_num_rows, A_num_cols, 0,
								NULL, NULL, NULL,
								CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
								CUSPARSE_INDEX_BASE_ZERO, computeType);
						}
						else {
							auto status = cusparseCreateCsr(&matA, A_num_rows, A_num_cols, A_nnz,
								__spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns, __spgemm_dat->dA_values,
								CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
								CUSPARSE_INDEX_BASE_ZERO, computeType);
							eigen_assert(status == 0);
						}

						if (B_nnz == 0)
						{
							auto status = cusparseCreateCsr(&matB, B_num_rows, B_num_cols, 0,
								NULL, NULL, NULL,
								CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
								CUSPARSE_INDEX_BASE_ZERO, computeType);
							eigen_assert(status == 0);
						}
						else {
							auto status = cusparseCreateCsr(&matB, B_num_rows, B_num_cols, B_nnz,
								__spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns, __spgemm_dat->dB_values,
								CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
								CUSPARSE_INDEX_BASE_ZERO, computeType);
							eigen_assert(status == 0);
						}
						{
							auto status = cusparseCreateCsr(&matC, C_num_rows1, C_num_cols1, C_nnz1,
								__spgemm_dat->dC_csrOffsets, __spgemm_dat->dC_columns, __spgemm_dat->dC_values,
								CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
								CUSPARSE_INDEX_BASE_ZERO, computeType);
							eigen_assert(status == 0);
						}
						//--------------------------------------------------------------------------
						// SpGEMM Computation
						auto status = cusparseSpGEMM_createDescr(&spgemmDesc);


						// ask bufferSize1 bytes for external memory
						bufferSize1 = 0;
						bufferSize2 = 0;
						status = cusparseSpGEMM_workEstimation(handle, opA, opB,
							&alpha, matA, matB, &beta, matC,
							computeType, ALG,
							spgemmDesc, &bufferSize1, NULL);
						eigen_assert(status == 0);
						if (bufferSize1 > __spgemm_dat->bufferSize1)
						{
							cudaFree(__spgemm_dat->dBuffer1);
							auto err = cudaMalloc((double**)&__spgemm_dat->dBuffer1, bufferSize1 * 2);
							__spgemm_dat->bufferSize1 = bufferSize1 * 2;
						}
						status = cusparseSpGEMM_workEstimation(handle, opA, opB,
							&alpha, matA, matB, &beta, matC,
							computeType, ALG,
							spgemmDesc, &bufferSize1, __spgemm_dat->dBuffer1);
						eigen_assert(status == 0);





						status = cusparseSpGEMM_compute(handle, opA, opB,
							&alpha, matA, matB, &beta, matC,
							computeType, ALG,
							spgemmDesc, &bufferSize2, NULL);
						eigen_assert(status == 0);
						if (bufferSize2 > __spgemm_dat->bufferSize2)
						{
							cudaFree(__spgemm_dat->dBuffer2);
							auto err = cudaMalloc((double**)&__spgemm_dat->dBuffer2, bufferSize2 * 2);
							__spgemm_dat->bufferSize2 = bufferSize2 * 2;
						}

						// ask bufferSize2 bytes for external memory
						//cudaMemsetAsync(__spgemm_dat->dBuffer2, 0, bufferSize2);

						// compute the intermediate product of A * B
						status = cusparseSpGEMM_compute(handle, opA, opB,
							&alpha, matA, matB, &beta, matC,
							computeType, ALG,
							spgemmDesc, &bufferSize2, __spgemm_dat->dBuffer2);
						eigen_assert(status == 0);
						// get matrix C non-zero entries C_nnz1
						C_num_rows1 = 0;
						C_num_cols1 = 0;
						auto prevC_nnz1 = C_nnz1;
						C_nnz1 = 0;
						status = cusparseSpMatGetSize(matC, &C_num_rows1, &C_num_cols1,
							&C_nnz1);
						eigen_assert(status == 0);
						//eigen_assert(C_nnz1 == 0);
						//std::cout << "rows" << C_num_rows1 << "cols" << C_num_cols1 <<"nnz"<<C_nnz1<< std::endl;
						// allocate matrix C

						// update matC with the new pointers
						//if(prevC_nnz1==0)
							status = cusparseCsrSetPointers(matC, __spgemm_dat->dC_csrOffsets, __spgemm_dat->dC_columns, __spgemm_dat->dC_values);
						//eigen_assert(status == 0);
						// if beta != 0, cusparseSpGEMM_copy reuses/updates the values of dC_values
							eigen_assert(status == 0);
						// copy the final products to the matrix C
						status = cusparseSpGEMM_copy(handle, opA, opB,
							&alpha, matA, matB, &beta, matC,
							computeType, ALG, spgemmDesc);

						//cusparseSpCsrgeam()
						//cusparsecsrmv
						eigen_assert(status == 0);
						//cusparseDestroyCsrgemm2Info(info);
						if (C_nnz1 == 0) {
							res[_ii].setZero();
						}
						else {
							res[_ii].resize(C_num_rows1, C_num_cols1);
							res[_ii].setZero();
							res[_ii].reserve(C_nnz1);
							res[_ii].resizeNonZeros(C_nnz1);
							Eigen::SparseMatrix<double, Eigen::RowMajor>ff2 = res[_ii];
							auto err = cudaMemcpyAsync(res[_ii].outerIndexPtr(), __spgemm_dat->dC_csrOffsets, (C_num_rows1 + 1) * sizeof(int), cudaMemcpyDeviceToHost, stream);
							err = cudaMemcpyAsync(res[_ii].innerIndexPtr(), __spgemm_dat->dC_columns, C_nnz1 * sizeof(int), cudaMemcpyDeviceToHost, stream);
							err = cudaMemcpyAsync(res[_ii].valuePtr(), __spgemm_dat->dC_values, C_nnz1 * sizeof(double), cudaMemcpyDeviceToHost, stream);
							e[_ii] += res[_ii];
							Eigen::SparseMatrix<double, Eigen::RowMajor>ff = res[_ii];
							Eigen::SparseMatrix<double, Eigen::RowMajor>gg = e[_ii];
						}
					}
				}

				
				//cuStreamDestroy(stream);
				//cusparseDestroy(handle);
#pragma omp critical
				{
				//dict.insert_or_assign(key, __spgemm_dat);
			//	dict[key] = __spgemm_dat;
				}

			}
		
		cudaDeviceSynchronize();
		//this->_mat[0] = e[0];
		//return ss.str();
		//int __mt = _mt;
		for (int tt = 0; tt < 40; tt++)
		{
#pragma omp parallel for
			for (int i = 0; i < __mt; i += 2)
			{
				if (i + 1 < __mt) {
					e[i] += e[i + 1];
				}
			}
			int _ct = 0;
#pragma omp parallel for ordered schedule(dynamic)
			for (int i = 0; i < __mt; i += 2)
			{
#pragma omp ordered
				if(i!=0)
					e[i / 2] = e[i];
#pragma omp atomic
				_ct++;
			}
			__mt = _ct;
			if (__mt == 1)break;
		}
		//eigen_assert(e[0].sum() < 1000 && e[0].sum() > -1000);
		ss<<"sum:" << e[0].sum() << std::endl;
		if (true) {
			if (this->_mat.size() == 0)this->_mat.resize(1);
			this->_mat[0].resize(nn, nn);
			this->_mat[0].reserve(nn * nn / 20);
			this->_mat[0] = e[0];
			//for (int i = 1; i < __mt; i++) {
			//	this->_mat[0] += e[i];
			//}
			this->_mat[0].makeCompressed();
			//this->_dmat = this->_mat[0];
		}
		auto end = high_resolution_clock::now();
		auto duration = duration_cast<milliseconds>(now - end);
		ss << "gpu" << duration.count() << "ms" << std::endl;
		omp_set_num_threads(_mt);
	return ss.str();
}


/*std::string kingghidorah::_mySparse::ofAtA_gpu(cuda* _cuda, _mySparse* A, bool sparse)
{
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> tmp;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> tmp2;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> res;
	auto ss = std::stringstream();
	auto now = high_resolution_clock::now();
	//this->join();
	if (e.size() < _mt)e.resize(_mt);
	tmp.resize(_mt);
	tmp2.resize(_mt);
	res.resize(_mt);
	int nn = this->cols();
	//cuInit(0);
	int _mt = omp_get_max_threads();
	for (int i = 0; i < _mt; i++)
	{
		e[i].resize(nn, nn);
		e[i].setZero();
		e[i].reserve(nn * nn / 20);
		e[i].makeCompressed();
		res[i].resize(nn, nn);
		res[i].setZero();
		res[i].reserve(nn * nn / 20);
		res[i].makeCompressed();
	}
	//std::vector <cusolverSpHandle_t> solversp_handles(_mt);
	//std::vector <CUstream> __streams(_mt);
	//cuInit(0);
	cudaSetDevice(0);

	//cusparseCreate(&sp_handle);

	//static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e;

	int __mt = STREAMCOUNT * _cuda->count();
	omp_set_num_threads(__mt);
	cudaSetDevice(0);
	int _ss = 20;// _nt / __mt / 4;
	//if (_ss == 0)_ss = 1;
	int job = 0;


	kingghidorah::spgemm ____spgemm_dat;
	auto key = std::tuple<kingghidorah::_mySparse*, int>(this, 0);
	if (dict.contains(key))
	{
		____spgemm_dat = dict[key];
	}
	else {
		____spgemm_dat.initialized = false;
	}
	//if (!____spgemm_dat.initialized)
	//{


	int max_rows = 0;
	int max_cols = 0;
	int max_nnz = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		int num_rows = this->_mat[ii].rows();
		int num_cols = this->_mat[ii].cols();
		int nnz = this->_mat[ii].nonZeros();
		max_rows = std::max(max_rows, num_rows);
		max_cols = std::max(max_cols, num_cols);
		max_nnz = std::max(max_nnz, nnz);

	}
	if (____spgemm_dat.initialized == false)
	{
		for (int _ii = 0; _ii < __mt; _ii++)
		{
			int device = _ii % _cuda->count();
			cudaSetDevice(device);

			kingghidorah::spgemm __spgemm_dat;
			auto __key = std::tuple<kingghidorah::_mySparse*, int>(this, _ii);

			__spgemm_dat.A_num_rows = max_cols;
			__spgemm_dat.A_num_cols = max_rows;
			__spgemm_dat.A_nnz = max_nnz;

			__spgemm_dat.B_num_rows = max_rows;
			__spgemm_dat.B_num_cols = max_cols;
			__spgemm_dat.B_nnz = max_nnz;

			__spgemm_dat.C_num_rows = max_cols;
			__spgemm_dat.C_num_cols = max_cols;
			__spgemm_dat.C_nnz = max_cols * max_cols / 20;

			__spgemm_dat.D_num_rows = max_cols;
			__spgemm_dat.D_num_cols = max_cols;
			__spgemm_dat.D_nnz = max_cols * max_cols / 20;

			auto err = cudaMalloc(&__spgemm_dat.dA_csrOffsets, sizeof(int) * (__spgemm_dat.A_num_rows + 1));
			err = cudaMalloc(&__spgemm_dat.dA_columns, sizeof(int) * __spgemm_dat.A_nnz);
			err = cudaMalloc(&__spgemm_dat.dA_values, sizeof(double) * __spgemm_dat.A_nnz);

			// allocate B
			err = cudaMalloc(&__spgemm_dat.dB_csrOffsets, (__spgemm_dat.B_num_rows + 1) * sizeof(int));
			err = cudaMalloc(&__spgemm_dat.dB_columns, __spgemm_dat.B_nnz * sizeof(int));
			err = cudaMalloc(&__spgemm_dat.dB_values, __spgemm_dat.B_nnz * sizeof(double));

			// allocate C
			err = cudaMalloc(&__spgemm_dat.dC_csrOffsets, (__spgemm_dat.C_num_rows + 1) * sizeof(int));
			err = cudaMalloc(&__spgemm_dat.dC_columns, __spgemm_dat.C_nnz * sizeof(int));
			err = cudaMalloc(&__spgemm_dat.dC_values, __spgemm_dat.C_nnz * sizeof(double));

			// allocate D
			err = cudaMalloc(&__spgemm_dat.dD_csrOffsets, (__spgemm_dat.D_num_rows + 1) * sizeof(int));
			err = cudaMalloc(&__spgemm_dat.dD_columns, __spgemm_dat.D_nnz * sizeof(int));
			err = cudaMalloc(&__spgemm_dat.dD_values, __spgemm_dat.D_nnz * sizeof(double));

			err = cudaMalloc((double**)&__spgemm_dat.dBuffer1, max_rows * max_cols * 2 + 1000);
			err = cudaMalloc((double**)&__spgemm_dat.dBuffer2, max_rows * max_cols * 2 + 1000);
			__spgemm_dat.bufferSize1 = max_rows * max_cols * 2 + 1000;
			__spgemm_dat.bufferSize2 = max_rows * max_cols * 2 + 1000;
			__spgemm_dat.initialized = true;
			dict[__key] = __spgemm_dat;
		}
	}

	//}

	cudaDeviceSynchronize();
	//for (int _tt = 0; _tt < _cuda->count(); _tt++)
	//{
#pragma omp parallel for
	for (int _ii = 0; _ii < __mt; _ii++)
	{
		int device = _ii % _cuda->count();
		cudaSetDevice(device);
		kingghidorah::spgemm* __spgemm_dat;
		auto __key = std::tuple<kingghidorah::_mySparse*, int>(this, _ii);
		__spgemm_dat = &dict[__key];

		CUstream stream = _cuda->__streams(device, _ii / _cuda->count());
		cusparseHandle_t handle = sp_handle[device][_ii / _cuda->count()];
		//CUstream stream;// = _cuda->__streams(device, _ii / _cuda->count());
		//cuStreamCreate(&stream,0);
		//cusparseHandle_t     handle;// = sp_handle;
		//cusparseCreate(&handle);
									//auto handle = _cuda->sp_handle[_ii % _cuda->count()][_ii / _cuda->count()];
		cusparseSetStream(handle, stream);


		double               alpha = 1.0;
		double               beta = 1.0;
		
		cusparseOperation_t opA = CUSPARSE_OPERATION_NON_TRANSPOSE;
		cusparseOperation_t opB = CUSPARSE_OPERATION_NON_TRANSPOSE;
		size_t bufferSize1 = 0, bufferSize2 = 0;
		int64_t C_num_rows = nn, C_num_cols = nn;
		int64_t D_num_rows = nn, D_num_cols = nn;

		cusparseSetPointerMode(handle, CUSPARSE_POINTER_MODE_HOST);
		int which = -1;
		int D_nnz = 0;
		int C_nnz = 0;
		cusparseMatDescr_t matA, matB, matC, matD;
		cusparseCreateMatDescr(&matA);
		cusparseSetMatType(matA, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(matA, CUSPARSE_INDEX_BASE_ZERO);
		cusparseCreateMatDescr(&matB);
		cusparseSetMatType(matB, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(matB, CUSPARSE_INDEX_BASE_ZERO);
		cusparseCreateMatDescr(&matC);
		cusparseSetMatType(matC, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(matC, CUSPARSE_INDEX_BASE_ZERO);
		cusparseCreateMatDescr(&matD);
		cusparseSetMatType(matD, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(matD, CUSPARSE_INDEX_BASE_ZERO);
		csrgemm2Info_t info = NULL;
		cusparseCreateCsrgemm2Info(&info);

		for (int kk = 0; kk < 100000; kk++)
		{
			//if (_ii == 0 )break;
			int S = 0;// _ii* _nt / __mt;
			int E = 0;// (_ii + 1)* _nt / __mt;

			//int S = 0;
			//int E = 0;
#pragma omp critical
			{
				S = job;
				E = S + _ss;
				if (E > _nt)E = _nt;
				job = E;
			}
			if (S >= _nt)break;


			//int counter = 0;


			for (int ii = S; ii < E; ii++)
			{
				int* nnzTotalDevHostPtr = 0;
				if (which == 0)
				{
					nnzTotalDevHostPtr = &D_nnz;
				}
				else {
					nnzTotalDevHostPtr = &C_nnz;
				}

				tmp[_ii] = (this->_mat[ii].transpose());// *coeff[ii].asDiagonal();
				tmp2[_ii] = (this->_mat[ii]);

				tmp2[_ii].makeCompressed();
				tmp[_ii].makeCompressed();

				for (int i = 0; i < tmp2[_ii].rows(); i++)
				{
					int S = tmp2[_ii].outerIndexPtr()[i];
					int E = tmp2[_ii].outerIndexPtr()[i+1];
					for (int k = S; k < E; k++) {
						*(tmp2[_ii].valuePtr()+k) *= coeff[ii][i];
					}
				}
				int A_num_rows = tmp[_ii].rows();
				int A_num_cols = tmp[_ii].cols();
				int A_nnz = tmp[_ii].nonZeros();
				int B_num_rows = tmp2[_ii].rows();
				int B_num_cols = tmp2[_ii].cols();
				int B_nnz = tmp2[_ii].nonZeros();
				C_num_rows = A_num_rows;
				C_num_cols = B_num_cols;
				D_num_rows = A_num_rows;
				D_num_cols = B_num_cols;
				if (A_nnz * B_nnz == 0)continue;
				which++;
				which = which % 2;

				//counter++;


				//std::cout << "A_nnz" << A_nnz << std::endl;
				//std::cout << "B_nnz" << B_nnz << std::endl;
				//--------------------------------------------------------------------------
				// Device memory management: Allocate and copy A, B


				// copy A
				if (A_nnz > 0)
				{
					auto err = cudaMemcpyAsync(__spgemm_dat->dA_csrOffsets, tmp[_ii].outerIndexPtr(), (A_num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);
					err = cudaMemcpyAsync(__spgemm_dat->dA_columns, tmp[_ii].innerIndexPtr(), A_nnz * sizeof(int), cudaMemcpyHostToDevice, stream);
					err = cudaMemcpyAsync(__spgemm_dat->dA_values, tmp[_ii].valuePtr(), A_nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
				}
				// copy B
				if (B_nnz > 0)
				{
					auto err = cudaMemcpyAsync(__spgemm_dat->dB_csrOffsets, tmp2[_ii].outerIndexPtr(), (B_num_rows + 1) * sizeof(int), cudaMemcpyHostToDevice, stream);
					err = cudaMemcpyAsync(__spgemm_dat->dB_columns, tmp2[_ii].innerIndexPtr(), B_nnz * sizeof(int), cudaMemcpyHostToDevice, stream);
					err = cudaMemcpyAsync(__spgemm_dat->dB_values, tmp2[_ii].valuePtr(), B_nnz * sizeof(double), cudaMemcpyHostToDevice, stream);
				}
				//--------------------------------------------------------------------------
				// CUSPARSE APIs
				if(C_nnz==0)
					cudaMemsetAsync(__spgemm_dat->dC_csrOffsets, 0, sizeof(int) * (nn + 1),stream);
				if(D_nnz==0)
					cudaMemsetAsync(__spgemm_dat->dD_csrOffsets, 0, sizeof(int) * (nn + 1),stream);

	
				// Create sparse matrix A in CSR format
				size_t bufferSize = 0;
					
		


				
				//__spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns, __spgemm_dat->dA_values,
					if (which == 0)
					{
						cusparseDcsrgemm2_bufferSizeExt(handle, nn, nn, A_num_cols, &alpha,
							matA, A_nnz, __spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns,
							matB, B_nnz, __spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns,
							&beta,
							matC, C_nnz, __spgemm_dat->dC_csrOffsets, __spgemm_dat->dC_columns,
							info,
							&bufferSize);
					}
					else {
						cusparseDcsrgemm2_bufferSizeExt(handle, nn, nn, A_num_cols, &alpha,
							matA, A_nnz, __spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns,
							matB, B_nnz, __spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns,
							&beta,
							matD, D_nnz, __spgemm_dat->dD_csrOffsets, __spgemm_dat->dD_columns,
							info,
							&bufferSize);

					}
				if (bufferSize > __spgemm_dat->bufferSize1)
				{
					cudaFree(__spgemm_dat->dBuffer1);
					auto err = cudaMalloc((double**)&__spgemm_dat->dBuffer1, bufferSize * 2);
					__spgemm_dat->bufferSize1 = bufferSize * 2;
				}
				//cudaMalloc(&buffer, bufferSize);

				// step 3: compute csrRowPtrC
				//cudaMalloc((void**)&csrRowPtrC, sizeof(int) * (m + 1));
				if (which == 0)
				{
					cusparseXcsrgemm2Nnz(handle, nn, nn, A_num_cols,
						matA, A_nnz, __spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns,
						matB, B_nnz, __spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns,
						matC, C_nnz, __spgemm_dat->dC_csrOffsets, __spgemm_dat->dC_columns,
						matD, __spgemm_dat->dD_csrOffsets, nnzTotalDevHostPtr,
						info, __spgemm_dat->dBuffer1);
				}
				else {
					cusparseXcsrgemm2Nnz(handle, nn, nn, A_num_cols,
						matA, A_nnz, __spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns,
						matB, B_nnz, __spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns,
						matD, D_nnz, __spgemm_dat->dD_csrOffsets, __spgemm_dat->dD_columns,
						matC, __spgemm_dat->dC_csrOffsets, nnzTotalDevHostPtr,
						info, __spgemm_dat->dBuffer1);

				}
				if (which == 0)
				{
					int D_base = 0;
					if (NULL != nnzTotalDevHostPtr) {
						D_nnz = *nnzTotalDevHostPtr;
					}
					else {
						cudaMemcpyAsync(&D_nnz, __spgemm_dat->dD_csrOffsets + nn, sizeof(int), cudaMemcpyDeviceToHost,stream);
						cudaMemcpyAsync(&D_base, __spgemm_dat->dD_csrOffsets, sizeof(int), cudaMemcpyDeviceToHost, stream);
						D_nnz -= D_base; //compute nnz for D
					}
				}
				else {
					int C_base = 0;
					if (NULL != nnzTotalDevHostPtr) {
						C_nnz = *nnzTotalDevHostPtr;
					}
					else {
						cudaMemcpyAsync(&C_nnz, __spgemm_dat->dC_csrOffsets + nn, sizeof(int), cudaMemcpyDeviceToHost, stream);
						cudaMemcpyAsync(&C_base, __spgemm_dat->dC_csrOffsets, sizeof(int), cudaMemcpyDeviceToHost, stream);
						C_nnz -= C_base; //compute nnz for D
					}
				}

				// step 4: finish sparsity pattern and value of C	
				//cudaMalloc((void**)&csrColIndC, sizeof(int) * nnzC);
				//cudaMalloc((void**)&csrValC, sizeof(double) * nnzC);
				// Remark: set csrValC to null if only sparsity pattern is required.
				if (which == 0)
				{
					cusparseDcsrgemm2(handle, nn, nn, A_num_cols, &alpha,
						matA, A_nnz, __spgemm_dat->dA_values, __spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns,
						matB, B_nnz, __spgemm_dat->dB_values, __spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns,
						&beta,
						matC, C_nnz, __spgemm_dat->dC_values, __spgemm_dat->dC_csrOffsets, __spgemm_dat->dC_columns,
						matD, __spgemm_dat->dD_values, __spgemm_dat->dD_csrOffsets, __spgemm_dat->dD_columns,
						info, __spgemm_dat->dBuffer1);
				}
				else {
					cusparseDcsrgemm2(handle, nn, nn, A_num_cols, &alpha,
						matA, A_nnz, __spgemm_dat->dA_values, __spgemm_dat->dA_csrOffsets, __spgemm_dat->dA_columns,
						matB, B_nnz, __spgemm_dat->dB_values, __spgemm_dat->dB_csrOffsets, __spgemm_dat->dB_columns,
						&beta,
						matD, D_nnz, __spgemm_dat->dD_values, __spgemm_dat->dD_csrOffsets, __spgemm_dat->dD_columns,
						matC, __spgemm_dat->dC_values, __spgemm_dat->dC_csrOffsets, __spgemm_dat->dC_columns,
						info, __spgemm_dat->dBuffer1);

				}

				//cudaMemcpyAsync(__spgemm_dat->dC_csrOffsets, __spgemm_dat->dD_csrOffsets, sizeof(int) * (nn+1), cudaMemcpyDeviceToDevice, stream);
				//cudaMemcpyAsync(__spgemm_dat->dC_columns, __spgemm_dat->dD_columns, sizeof(int) * D_nnz, cudaMemcpyDeviceToDevice, stream);
				//cudaMemcpyAsync(__spgemm_dat->dC_values, __spgemm_dat->dD_values, sizeof(double) * D_nnz, cudaMemcpyDeviceToDevice, stream);
				//C_nnz = D_nnz;
				//eigen_assert(status == 0);
			
			}
		}
		cusparseDestroyCsrgemm2Info(info);

		if (which == 0)
		{
			if (D_nnz == 0) {
				res[_ii].setZero();
			}
			else {
				res[_ii].resize(D_num_rows, D_num_cols);
				res[_ii].setZero();
				res[_ii].reserve(D_nnz);
				res[_ii].resizeNonZeros(D_nnz);
				auto err = cudaMemcpyAsync(res[_ii].outerIndexPtr(), __spgemm_dat->dD_csrOffsets, (D_num_rows + 1) * sizeof(int), cudaMemcpyDeviceToHost, stream);
				err = cudaMemcpyAsync(res[_ii].innerIndexPtr(), __spgemm_dat->dD_columns, D_nnz * sizeof(int), cudaMemcpyDeviceToHost, stream);
				err = cudaMemcpyAsync(res[_ii].valuePtr(), __spgemm_dat->dD_values, D_nnz * sizeof(double), cudaMemcpyDeviceToHost, stream);
				cudaStreamSynchronize(stream);
				e[_ii] += res[_ii];
			}
		}
		else if(which==1){
			if (C_nnz == 0) {
				res[_ii].setZero();
			}
			else {
				res[_ii].resize(C_num_rows, C_num_cols);
				res[_ii].setZero();
				res[_ii].reserve(C_nnz);
				res[_ii].resizeNonZeros(C_nnz);
				auto err = cudaMemcpyAsync(res[_ii].outerIndexPtr(), __spgemm_dat->dC_csrOffsets, (C_num_rows + 1) * sizeof(int), cudaMemcpyDeviceToHost, stream);
				err = cudaMemcpyAsync(res[_ii].innerIndexPtr(), __spgemm_dat->dC_columns, C_nnz * sizeof(int), cudaMemcpyDeviceToHost, stream);
				err = cudaMemcpyAsync(res[_ii].valuePtr(), __spgemm_dat->dC_values, C_nnz * sizeof(double), cudaMemcpyDeviceToHost, stream);
				cudaStreamSynchronize(stream);
				e[_ii] += res[_ii];
			}
		}

		//cuStreamDestroy(stream);
		//cusparseDestroy(handle);
//#pragma omp critical
	//	{
			//dict.insert_or_assign(key, __spgemm_dat);
		//	dict[key] = __spgemm_dat;
		//}

	}

	cudaDeviceSynchronize();

	//int __mt = _mt;
	for (int tt = 0; tt < 40; tt++)
	{
#pragma omp parallel for
		for (int i = 0; i < __mt; i += 2)
		{
			if (i + 1 < __mt) {
				e[i] += e[i + 1];
			}
		}
		int _ct = 0;
#pragma omp parallel for ordered schedule(dynamic)
		for (int i = 0; i < __mt; i += 2)
		{
#pragma omp ordered
			if (i != 0)
				e[i / 2] = e[i];
#pragma omp atomic
			_ct++;
		}
		__mt = _ct;
		if (__mt == 1)break;
	}
	//eigen_assert(e[0].sum() < 1000 && e[0].sum() > -1000);
	ss << "sum:" << e[0].sum() << std::endl;
	if (true) {
		if (this->_mat.size() == 0)this->_mat.resize(1);
		this->_mat[0].resize(nn, nn);
		this->_mat[0].reserve(nn * nn / 20);
		this->_mat[0] = e[0];
		//for (int i = 1; i < __mt; i++) {
		//	this->_mat[0] += e[i];
		//}
		this->_mat[0].makeCompressed();
		//this->_dmat = this->_mat[0];
	}
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(now - end);
	ss << "gpu" << duration.count() << "ms" << std::endl;
	omp_set_num_threads(_mt);
	return ss.str();
}*/

void kingghidorah::_mySparse::_freeze() {
	//this->_dmat = this->_mat[0];
}
std::string kingghidorah::_mySparse::_ofAtA(_mySparse* A)
{
	/*if (__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&___dmat, sizeof(double) * A->cols() * A->cols()*2);
		}
		else {
			___dmat = (double*)malloc(sizeof(double) * A->cols() * A->cols()*2);
		}
	}*/
	//__r = A->cols();
	//__c = A->cols();
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, A->cols(), A->cols());
	std::stringstream ss;
#ifdef _DEBUG
	ss << _nt << std::endl;
	for (int ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{
			ss << ii << "th" << std::endl;
			ss << A->_mat[ii].rows() << "," << A->_mat[ii].cols() << "," << coeff[ii].size() << std::endl;

		}
	}
#endif

	//this->_dmat.resize(A->cols(), A->cols());
	_dmat.setZero(A->cols(), A->cols());
	//_tmp.setZero(__r, __c);
	for (int ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{

			auto _ret = (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
			_dmat += _ret;
		}
	}
	return ss.str();
}
std::string kingghidorah::_mySparse::info()
{
	return "viennacl has been abandoned";
}
void kingghidorah::_mySparse::join()
{
	//eigen_assert(this->_mat[0].nonZeros() > 0);
	int __rows = this->rows();
	int offset = this->_mat[0].rows();
	int offsetnonzeros = this->_mat[0].nonZeros();;
	this->_mat[0].conservativeResize(__rows,this->_mat[0].cols());
	this->coeff[0].conservativeResize(__rows, 1);
	int nnz = 0;
	//for (int i = 0; i < _nt; i++)
	//{
	//	nnz += this->_mat[i].nonZeros();
	//}
	this->_mat[0].resizeNonZeros(nnz);
	this->_mat[0].reserve(nnz);
	this->_mat[0].makeCompressed();

	for (int i = 1; i < _nt; i++)
	{
		this->_mat[0].middleRows(offset, this->_mat[i].rows()) = this->_mat[i];
		//this->_mat[i].makeCompressed();
		//int nnz0 = this->_mat[i].nonZeros();
		//memcpy(this->_mat[0].valuePtr() + offsetnonzeros, this->_mat[i].valuePtr(), sizeof(double) * nnz0);
		//memcpy(this->_mat[0].innerIndexPtr() + offsetnonzeros, this->_mat[i].innerIndexPtr(), sizeof(double) * nnz0);
		
		this->coeff[0].middleRows(offset, this->_mat[i].rows()) = this->coeff[i];
		offset += this->_mat[i].rows();
		//for (int k = 1; k < this->_mat[i].rows()+1; k++)
		//{
		//	this->_mat[0].outerIndexPtr()[offset + k] = offsetnonzeros + this->_mat[i].outerIndexPtr()[k];
		//}
		//offsetnonzeros += nnz0;
		
	}
}
void kingghidorah::_mySparse::_ofAtB(_mySparse* B, _mySparse* C)
{
	this->freeze(true);
	B->freeze(false);

	int nn = this->_cols();
	int mm = B->_cols();
	/*f(C->__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&C->___dmat, sizeof(double) * nn * mm*2);
		}
		else {
			C->___dmat = (double*)malloc(sizeof(double) * nn * mm*2);
		}
	}*/
	//C->__r = nn;
	//C->__c = mm;
	C->_dmat.resize(nn, mm);
	//C->_tmp.resize(nn, mm);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, mm);
	//Eigen::Map<Eigen::MatrixXd> C_dmat(C->___dmat, nn, mm);
	//C_dmat.setZero();

	//int mt = omp_get_max_threads();

	int ss = mm / _mt;
	auto left = _dmat.transpose();
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
	C->_mat[0].makeCompressed();
}
void kingghidorah::_mySparse::_ofBtAB_qr(_mySparse* B, Eigen::VectorXd* b, _mySparse* C, Eigen::VectorXd* ret)
{
	this->join();
	B->join();

	static Eigen::SparseMatrix<double, Eigen::ColMajor> _a;
	static Eigen::SparseMatrix<double,Eigen::ColMajor> q;
	static Eigen::SparseMatrix<double, Eigen::ColMajor > tmp;
	this->_mat[0].makeCompressed();
	_a = this->_mat[0];
	Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int>> qr;
	qr.setPivotThreshold(0.00000000001);
	qr.compute(_a);
	q = qr.m_Q;
	tmp = q * B->_mat[0];
	C->_dmat = tmp.transpose() * tmp;
}
void kingghidorah::_mySparse::_ofBtAB(_mySparse* B, Eigen::VectorXd* b, _mySparse* C, Eigen::VectorXd* ret)
{
	Eigen::MatrixXd D;

	int nn = B->_mat[0].cols();
	int kk = _dmat.cols();// __c;


	/*if (C->__r == 0)
	{
		if (false){//__cuinit) {
			cudaMallocHost(&C->___dmat, sizeof(double) * nn * nn*2);
		}
		else {
			C->___dmat = (double*)malloc(sizeof(double) * nn * nn*2);
		}
	}*/
	//C->__r = nn;
	//C->__c = nn;
	//C->_dmat.resize(nn, mm);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> C_dmat(C->___dmat, nn, nn);


	//C->_dmat.resize(nn, nn);
	C->_dmat.setZero(nn,nn);
	//C->_tmp.setZero(nn,nn);

	//int mt = omp_get_max_threads();

	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];

	D.resize(nn, kk);
	int ss = kk / _mt / 2;
	
#pragma omp parallel for
	for (int ii = 0; ii < kk; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= kk)E = kk;
		D.middleCols(S, E - S).noalias() = left * mid.middleCols(S, E - S);
	}
	//D.noalias() = left * mid;
	ss = nn / _mt / 2;
#pragma omp parallel for
	for (int ii = 0; ii < nn; ii += ss)
	{
		int S = ii;
		int E = ii + ss;
		if (E >= nn)E = nn;
		C->_dmat.middleCols(S, E - S).noalias() = D * right.middleCols(S, E - S);
	}
	//C_dmat.noalias() = D * right;
	//Eigen::Map<Eigen::VectorXd> b(ptr, N);
	*ret = D * *b;
}

void kingghidorah::_mySparse::ofAtB(_mySparse* B, bool sparse)
{
	static std::map<_mySparse*, Eigen::SparseMatrix<double, Eigen::RowMajor>> dict2;
	static std::map< Eigen::SparseMatrix<double, Eigen::RowMajor>*, std::vector<std::vector<int>>>  map;
	static std::vector<std::vector<int>> index;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e;
	static std::vector<Eigen::SparseMatrix<double, Eigen::RowMajor>> e2;
	auto ss = std::stringstream();
	auto now = high_resolution_clock::now();
	int nn = this->cols();
	int mm = B->cols();
	int _mt = omp_get_max_threads();
	omp_set_num_threads(_mt);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	int __mt = _mt;// std::min(_nt / 10, _mt * 10);
	if (e.size() < __mt)
	{
		e.resize(__mt);
		e2.resize(__mt);
	}
	Eigen::initParallel();
	Eigen::setNbThreads(1);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	ss << _nt << ":nt:" << std::endl;
	ss << _mt << ":_mt:" << std::endl;
	Eigen::initParallel();
	int job = 0;
	int sss = _nt / __mt / 8;
	int __nt2 = _nt;
	if (sss == 0)sss = 1;
	Eigen::SparseMatrix<double, Eigen::RowMajor>* prevmat;
	if (dict2.contains(this)) {
		prevmat = &dict2[this];
	}
	else {
		Eigen::SparseMatrix<double, Eigen::RowMajor> _prevmat(nn, nn);
		dict2[this] = _prevmat;
		prevmat = &dict2[this];
		prevmat->resize(nn, nn);
		prevmat->reserve(nn * nn / 10);
	}
	//prevmat->makeCompressed();
	std::vector<std::vector<int>>* _map = 0;
	if (map.contains(prevmat))
	{
		_map = &map[prevmat];
	}
#pragma omp parallel for
	for (int i = 0; i < __mt; i++) {
		e[i].resize(nn, mm);
		if (e[i].nonZeros() != prevmat->nonZeros())
		{
			e[i] = *prevmat;
			e[i].makeCompressed();
		}
		memset(e[i].valuePtr(), 0, sizeof(double) * prevmat->nonZeros());
		e2[i].resize(nn, mm);
		e2[i].reserve(nn * mm / 100);
	}
	index.resize(__mt);
#pragma omp parallel for schedule(static)
	for (int _ii = 0; _ii < __mt; _ii++)
	{
		int S = 0;
		int E = 0;
		//int K = 0;
		//auto _e = e[_ii];
		for (int tt = 0; tt < 40; tt++)
		{
#pragma omp critical
			{
				S = job;
				E = job + sss;
				if (E > _nt)E = _nt;
				job = E;
			}
			if (S >= _nt)break;
			for (int ii = S; ii < E; ii++)
			{
				/*Eigen::SparseMatrix<double> __coeff(coeff[ii].size(), coeff[ii].size());
				for (int k = 0; k < coeff[ii].size(); k++)
				{
					__coeff.insert(k, k) = coeff[ii](k);
				}*/
				//e[_ii] = 
//#pragma omp critical
				{

					e2[_ii] = this->_mat[ii].transpose() * coeff[ii].asDiagonal() *B->_mat[ii];
					index[_ii].resize(e2[_ii].nonZeros());
					if (_map == 0)
					{
					}
					else {
						int count = 0;
						for (int k = 0; k < nn; ++k) {
							for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(e2[_ii], k); it; ++it) {
								//e[0].coeffRef(it.row(), it.col()) += it.value();
								index[_ii][count] = (*_map)[it.row()][it.col()];
								count++;
							}
						}
						//e[_ii].makeCompressed();
					}

//#pragma omp critical
					{
						if (_map == 0)
						{
							e[_ii] += e2[_ii];
						}
						else {
							int* ptr = &index[_ii][0];
							for (int k = 0; k < nn; ++k) {
								for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(e2[_ii], k); it; ++it) {
									//e[0].coeffRef(it.row(), it.col()) += it.value();
									*(e[_ii].valuePtr()+*ptr) += it.value();
									ptr++;
								}
							}
							//e[_ii].makeCompressed();
						}
					}
				}
			}
		}
		//e[_ii].makeCompressed();
	}
	//this->_mat[0] = *prevmat;
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	Eigen::setNbThreads(_mt);

	for (int tt = 0; tt < 40; tt++)
	{
#pragma omp parallel for
		for (int i = 0; i < __mt; i += 2)
		{
			if (i + 1 < __mt) {
				if (_map == 0 || e[i].nonZeros() != e[i + 1].nonZeros())
				{
					e[i] += e[i + 1];
				}
				else {
					Eigen::Map<Eigen::VectorXd> map1(e[i].valuePtr(), e[i].nonZeros());
					Eigen::Map<Eigen::VectorXd> map2(e[i + 1].valuePtr(), e[i].nonZeros());
					map1 += map2;
				}
			}
		}
		int _ct = 0;
#pragma omp parallel for ordered schedule(dynamic)
		for (int i = 0; i < __mt; i += 2)
		{
#pragma omp ordered
			if (_map == 0 || e[i].nonZeros() != e[i / 2].nonZeros())
			{
				e[i / 2] = e[i];
			}
			else
			{
				memcpy(e[i / 2].valuePtr(), e[i].valuePtr(), sizeof(double) * e[i].nonZeros());
			}
#pragma omp atomic
			_ct++;
		}
		__mt = _ct;
		if (__mt == 1)break;
	}
	if (true) {
		if (this->_mat.size() == 0)this->_mat.resize(1);
		this->_mat[0].resize(nn, nn);
		this->_mat[0].reserve(nn * nn / 20);
		this->_mat[0] = e[0];
		//for (int i = 1; i < __mt; i++) {
		//	this->_mat[0] += e[i];
		//}
		this->_mat[0].makeCompressed();
		if (_map == 0 || prevmat->nonZeros() != this->_mat[0].nonZeros())
		{
			*prevmat = this->_mat[0];
			prevmat->makeCompressed();
			dict2[this] = *prevmat;
			//build map
			std::vector<std::vector<int>> __map;
			__map.resize(this->_mat[0].rows());
			for (int i = 0; i < this->_mat[0].rows(); i++)
			{
				__map[i].resize(this->_mat[0].cols());
			}
			for (int k = 0; k < prevmat->outerSize(); ++k) {
				for (Eigen::SparseMatrix<double, Eigen::RowMajor>::InnerIterator it(*prevmat, k); it; ++it) {
					int S = prevmat->outerIndexPtr()[k];
					int E = prevmat->outerIndexPtr()[k + 1];

					for (int tt = S; tt < E; tt++)
					{
						if (prevmat->innerIndexPtr()[tt] == it.col())
							__map[it.row()][it.col()] = tt;
					}
				}
			}
			map[prevmat] = __map;
		}
	}


	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	if (!sparse) {
		//this->_dmat.setZero(nn, nn);
		//this->_tmp.setZero(nn, nn);
		/*if (__r == 0 || __c == 0)
		{
			if (false)//__cuinit)
			{
				cudaMallocHost(&___dmat, sizeof(double) * nn * nn * 2);
			}
			else {
				___dmat = (double*)malloc(sizeof(double) * nn * nn * 2);
			}
		}*/
		//__r = nn;
		//__c = nn;
		/*Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, nn);
		_dmat = e[0];
		for (int i = 1; i < _mt; i ++) {
			_dmat += e[i];
		}*/
		//this->_dmat = x;
	}
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	//this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	ss << "sum:" << e[0].sum() << std::endl;
	return;// ss.str();
}
void kingghidorah::_mySparse::Atb(double* ptr, int N, Eigen::VectorXd* c)
{
	c->resize(this->cols());
	c->setZero();
	int offset = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		int ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{

			/*Eigen::SparseMatrix<double> __coeff(coeff[ii].size(), coeff[ii].size());
			for (int k = 0; k < coeff[ii].size(); k++)
			{
				__coeff.insert(k, k) = coeff[ii](k);
			}*/
			*c += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		}
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
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{
			/*Eigen::SparseMatrix<double> __coeff(coeff[ii].size(), coeff[ii].size());
			for (int k = 0; k < coeff[ii].size(); k++)
			{
				__coeff.insert(k, k) = coeff[ii](k);
			}*/

			ret += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		}
		offset += ee;
	}
	return ret;
}
Eigen::VectorXd kingghidorah::_mySparse::_Atb(double* ptr, int N)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::VectorXd> b(ptr, N);
	return _dmat.transpose() * b;
}
void kingghidorah::_mySparse::merge()
{
	this->ofDat();
}
void kingghidorah::_mySparse::computeQR()
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::HouseholderQR<Eigen::MatrixXd> qr;
	qr.compute(_dmat);
}
void kingghidorah::_mySparse::computeLU()
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	lu.compute(_dmat);
}
void kingghidorah::_mySparse::computeLLT(Eigen::LLT<Eigen::MatrixXd>* _LLT)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	_LLT->compute(_dmat);
}
int kingghidorah::_mySparse::nonzeros() {
	int _ret = 0;
	for (int ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)

		_ret += _mat[ii].nonZeros();
	}
	return _ret;
}
void kingghidorah::_mySparse::solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	lu.compute(_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_dmat.cols());
	//Eigen::VectorXd x(_dmat.cols());
	//x.setZero();
	*ret = lu.solve(*rhs);
}


void kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int device) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::VectorXd x = *ret;
	//this->_freeze();
	int N = rhs->rows();
	if (!cuda->valid())
	{
		x(0) = 10;
		return;
	}
	cudaSetDevice(device);
	auto solver = cuda->solver(device, 0);
	auto blas = cuda->blas(device);
	//auto stream=streams[device];
	//cudaStreamCreate(&stream);
	//cusolverDnSetStream(solver, stream);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	x.resize(N);
	//x.setZero();


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

	//auto now = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
	//int devInfo_on_cpu = 0;
	//cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);


	//if (0 != devInfo_on_cpu) {
	//	x(0) = devInfo_on_cpu;
	//	return;
	//}


	cusolverDnDpotrs(solver, CUBLAS_FILL_MODE_LOWER, N, 1, gpu_matrix, N, gpu_rhs, N, devInfo_on_gpu);
	// auto end = std::chrono::high_resolution_clock::now();
	//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	//std::cout << "Dn:" << duration.count() << "ms" << std::endl;

	//cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);

	//if (devInfo_on_cpu != 0) {
	//	x(0) = 24;
	//	return;
	//}

	cudaMemcpy(ret->data(), gpu_rhs, sizeof(double) * N, cudaMemcpyDeviceToHost);

	//cudaFree(work);

	//cudaFree(devInfo_on_gpu);

	cudaDeviceSynchronize();
	//cudaStreamDestroy(stream);
	return;
}
Eigen::MatrixXd kingghidorah::_mySparse::_solve0(_myLLT* LLT, _mySparse* mat)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//this function assumes that LLT decomposition has been done already
	Eigen::MatrixXd ret(this->_dmat.cols(), mat->_dmat.cols());
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
/*
void kingghidorah::_mySparse::_solveI_gpu_mg(kingghidorah::cuda* cuda, _mySparse* ret)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	int N = this->_dmat.cols();
	//Eigen::MatrixXd x(N,nn);
	
	//ret->__r = N;
	//ret->__c = N;
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	//ret->_dmat.resize(N, N);
	if (!cuda->valid())return;
	ret->_dmat.setZero(N, N);
	//ret->_tmp.setZero(N, N);
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

	cudaLibMgMatrixDesc_t descrA = cuda->descrA;

	cudaLibMgGrid_t gridA = cuda->gridA;

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

		ret->_dmat.data(),
		lda
		);


#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			ret->_dmat(j, i) = ret->_dmat(i, j);
		}
	}


	//if (NULL != array_d_A) free(array_d_A);
	//if (NULL != array_d_B) free(array_d_B);
	//if (NULL != array_d_work) free(array_d_work);

}
*/
void kingghidorah::_mySparse::_solveI(_mySparse* ret)
{
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);	
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::RowMajor> llt;
	llt.compute(this->_mat[0]);
	int nn = this->_mat[0].rows();
	/*if(ret->__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * nn * nn*2);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * nn * nn * 2);
		}
	}*/
	//ret->__r = nn;
	//ret->__c = nn;
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, nn, nn);
	ret->_dmat.resize(nn, nn);
	if (I.rows() != nn || I.cols() != nn)
	{
		I.resize(nn, nn);
		I.setIdentity();
	}
	//int mt = omp_get_max_threads();
	int ee = _mt * 2;
#pragma omp parallel for
	for (int i = 0; i < ee; i++)
	{
		int S = i * nn / ee;
		int E = (i + 1) * nn / ee;
		ret->_dmat.middleCols(S, E - S) = llt.solve(I.middleCols(S, E - S));
	}
}

void initidentiy(kingghidorah::cuda* cuda, int N) {
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
				cudaMemcpyAsync(gpu_rhs + (i + i * N), &_b, 1 * sizeof(double), cudaMemcpyHostToDevice, stream);
			}
		}
	}
}
std::string kingghidorah::_mySparse::_solveI_gpu_single(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();
	std::stringstream sss;
	int N = _dmat.cols();// __c;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

	/*if (ret->__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N * 2);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N * 2);
		}
	}*/
	//ret->__r = N;
	//ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int nn = cuda->count();
	initidentiy(cuda, N);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret->_dmat.setZero(N, N);
	//ret->_tmp.setZero(N, N);




	auto solver = cuda->solver(cuda->fastest(), 0);

	cudaSetDevice(cuda->fastest());
	//auto stream = streams[cuda->fastest()];
	double* gpu_matrix = cuda->work_M(cuda->fastest());
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
		cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
	}

	cudaDeviceSynchronize();
	return sss.str();

}
std::string kingghidorah::_mySparse::_solveI_gpu_omp(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();
	std::stringstream sss;
	int N = this->_dmat.cols();

	/*if (ret->__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N * 2);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N * 2);
		}
	}*/
	//ret->__r = N;
	//ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int nn = cuda->count();
	initidentiy(cuda, N);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret->_dmat.setZero(N, N);
	//ret->_tmp.setZero(N, N);


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
					cudaMemcpyAsync(cuda->work_M(ii) + S * N, ret->_dmat.data() + S * N, sizeof(double) * N * (E - S), cudaMemcpyHostToDevice, _stream);
				}
			}
		}
		//cudaStreamDestroy(_stream);
	}

	int job = 0;
	int ss = N / cuda->count() / STREAMCOUNT / 4;
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
		auto solver = cuda->solver(ii, 0);

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

#pragma omp critical
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
				cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_rhs + S * N, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
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
		sss << i << "-" << count[i] << "  ";
	}
	sss << std::endl;
	return sss.str();
}

void kingghidorah::_mySparse::_solveI_gpu(kingghidorah::cuda* cuda, _mySparse* ret)
{
	//this->_freeze();

	int N = _dmat.cols();// __c;
	/*if (ret->__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * N * 2);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * N * 2);
		}
	}*/
	//ret->__r = N;
	//ret->__c = N;
	//ret->_dmat.resize(N, N);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, N, N);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret->_dmat.setZero(N, N);
	//ret->_tmp.setZero(N, N);
	//initidentiy(cuda, N);
	if (!cuda->valid())return;

	int ii = cuda->fastest();
	auto solver = cuda->solver(ii, 0);
	//cudaStream_t stream = streams[ii];
	//cusolverDnSetStream(solver, stream);
	cudaSetDevice(ii);
	double* m_gpu = cuda->work_M(ii);

	cudaMemcpyAsync(m_gpu, _dmat.data(), sizeof(double) * N * N, cudaMemcpyHostToDevice, cuda->__streams(cuda->fastest(), 0));
	double* work;
	int work_size;
	int work_size1 = 0;
	int work_size2 = 0;
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size1);
	cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size2);
	work_size = N * N;// std::max(work_size1, work_size2);
	work = cuda->work(work_size, ii, cuda->__streams(cuda->fastest(), 1));
	//cudaMallocAsync(&work, sizeof(double) * work_size, stream);
	cudaMemset(work, 0, sizeof(double) * work_size1);
	int* devInfo = cuda->info(ii);

	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work, work_size1, devInfo);
	double* work2 = cuda->work_rhs(ii);
	//cudaMalloc(&work2, sizeof(double) * work_size2);
	cudaMemset(work2, 0, sizeof(double) * work_size2);
	cusolverDnDpotri(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work2, work_size2, devInfo);
	//cudaFree(work2);

	cudaMemcpy(ret->_dmat.data(), m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	ret->_dmat.triangularView<Eigen::Upper>() = ret->_dmat.triangularView<Eigen::Lower>().transpose();
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
	int nn = mat->_dmat.cols();// __c;
	int N = this->_dmat.cols();// __c;
	//Eigen::MatrixXd x(N,nn);

	/*if (ret->__r == 0)
	{
		if (false)//__cuinit)
		{
			cudaMallocHost(&ret->___dmat, sizeof(double) * N * nn * 2);
		}
		else {
			ret->___dmat = (double*)malloc(sizeof(double) * N * nn * 2);
		}
	}*/
	//ret->__c = nn;
	//ret->__r = N;
	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest(), 0);
	auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, nn);
	//Eigen::Map<Eigen::MatrixXd> mat_dmat(mat->___dmat, mat->__r, mat->__c);
	ret->_dmat.setZero(N, nn);
	//ret->_tmp.setZero(N, N);
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
		auto solver = cuda->solver(i, 0);
		auto blas = cuda->blas(i);
		double* _gpu_matrix = cuda->work_M(i);
		double* gpu_rhs = cuda->work_rhs(i);
		cudaMemcpy(_gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
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
			cudaMemcpy(ret->_dmat.data() + nextjob * N, gpu_rhs + nextjob * N, sizeof(double) * N * (end - start), cudaMemcpyDeviceToHost);
		}
		cudaFree(_devInfo_on_gpu);
	}
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Dn" << duration.count() << "ms" << std::endl;

	if (exit)
	{
		ret->_dmat(0, 0) = 4;
		return;
	}
}

void kingghidorah::_mySparse::_solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::RowMajor> LLT;
	LLT.compute(_mat[0]);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_mat[0].cols());
	ret->setZero();
	//Eigen::VectorXd x(_mat[0].rows());
	//x.setZero();
	*ret = LLT.solve(*rhs);
	//return x;
}
void kingghidorah::_mySparse::__solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	Eigen::LLT<Eigen::MatrixXd> LLT;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

	LLT.compute(_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->resize(_dmat.cols());
	ret->setZero();
	*ret = LLT.solve(*rhs);
}
Eigen::MatrixXd kingghidorah::_mySparse::inv() {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	lu.compute(_dmat);
	Eigen::MatrixXd I(_dmat.rows(), _dmat.cols());
	I.setIdentity();
	return lu.solve(I);
}

Eigen::MatrixXd kingghidorah::_mySparse::solve0(_mySparse* rhs)
{
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd> rhs_dmat(rhs->___dmat, rhs->__r, rhs->__c);

	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//_mat.makeCompressed();
	lu.compute(_dmat);
	Eigen::MatrixXd _x(_dmat.cols(), rhs->_dmat.cols());
	_x = lu.solve(rhs->_dmat);

	return _x;// .sparseView(0.000000000000001, 1.0);
}


void kingghidorah::_mySparse::minus(_mySparse* m) {
	this->_freeze();
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> m_dmat(m->___dmat, m->__r, m->__c);

	_dmat = _dmat - m->_dmat;
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
}
void kingghidorah::_mySparse::clearcoeff() {
	for (int ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{
			this->coeff[ii] = Eigen::VectorXd::Ones(this->_mat[ii].rows());
		}
		else {
			this->coeff[ii].resize(0);
		}
	}
}
Eigen::SparseMatrix<double, Eigen::RowMajor> id;
void kingghidorah::_mySparse::addsmallidentity(double salt, bool sparse, bool dense) {

	if (dense)
	{
		id.resize(_dmat.rows(), _dmat.cols());
		id.setIdentity();
		//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

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