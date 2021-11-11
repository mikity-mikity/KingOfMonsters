// eigenwrapper.cpp : Defines the functions for the static library.
//
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
double* identity = 0;
int previdentiyN = 0;
kingghidorah::cuda::cuda(int N) {
	omp_set_dynamic(false);
	omp_set_num_threads(16);
	prevT_A = 0;
	prevN = 0;
	prevwn = 0;
	CUresult res;
	res = cuInit(0);
	if (identity != 0)cudaFreeHost(identity);
	identity = 0;
	previdentiyN = 0;

	auto err = cuDeviceGetCount(&_count);
	//mg_solver = 0;
	//cusolverMgCreate(&mg_solver);
	if (_count > 0)
	{
		for (int i = 0; i < _count; i++)_deviceList[i] = i;
	}
	if (_count > 1)
	{
		assert(0 == enablePeerAccess(_count, _deviceList));
	}




	//std::cout << err << "yay" << std::endl;
	for (int i = 0; i < MAXDEVICE; i++)
	{
		solver_handle[i] = 0;
		cublas_handle[i] = 0;
	}


	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		auto status = cusolverDnCreate(&solver_handle[ii]);
		auto status2 = cublasCreate(&cublas_handle[ii]);

		if (status == cusolverStatus_t::CUSOLVER_STATUS_SUCCESS)
		{
			initialized = true;
			failed = false;
		}
		else {
			initialized = false;
			failed = true;
			solver_handle[ii] = 0;
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
			for (int i = 0; i < _N; i++)rhs[i] = i;
			m.ofDat();
			m.clearcoeff();
			m._ofAtA(&m);
			m._solve0_gpu(this, rhs, _N, ii);
			delete[] rhs;
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			speed[ii] = duration.count();
		}
	}
	_fastest = std::distance(speed.begin(), std::min_element(speed.begin(), speed.end()));
	//_fastest = 1;
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
	/*cusolverStatus_t status = cusolverMgDeviceSelect(
		mg_solver,
		_count,
		_deviceList);*/ //seem like cannot be called twice. So call this only once here.
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
void kingghidorah::cuda::dispose() {
	if (valid())
	{
		if (identity != 0)cudaFreeHost(identity);
		identity = 0;
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
			if (solver_handle[i] != 0)
				cusolverDnDestroy(solver_handle[i]);
			solver_handle[i] = 0;
			if (cublas_handle != 0)
				cublasDestroy(cublas_handle[i]);
			cublas_handle[i] = 0;
		}
		if (mg_solver != 0)
		{
			//cusolverMgDestroy(mg_solver);
		}
		mg_solver = 0;
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
			cudaFree(__work[device]);
		cudaMalloc(&__work[device], N * sizeof(double));
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
cusolverDnHandle_t& kingghidorah::cuda::solver(int ii) {
	return solver_handle[ii];
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
}
double kingghidorah::_mySparse::L2Norm(double* ptr1, int N1, double* ptr2, int N2) {
	auto a = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr1, N1);
	auto b = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr2, N2);
	return a.transpose() * this->_mat[0] * b;
}
Eigen::VectorXd kingghidorah::_mySparse::Vector(double* ptr1, int N1) {
	auto a = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr1, N1);
	return this->_mat[0] * a;
}
void kingghidorah::_mySparse::plus(_mySparse* m, double sc) {
	this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
	this->_dmat = this->_mat[0];// this->_dmat + m->_dmat * sc;
}
void kingghidorah::_mySparse::setmat(Eigen::SparseMatrix<double> mat, int ii) {
	this->_mat[ii] = mat;
}
void kingghidorah::_mySparse::setmat(const Eigen::MatrixXd& mat) {
	this->_dmat = mat;
}
double kingghidorah::_mySparse::at(int i, int ii) {
	return this->_mat[ii].data().value(i);
}
int kingghidorah::_mySparse::num_elem(int ii)
{
	return this->_mat[ii].data().size();
}
double kingghidorah::_mySparse::_at(int i) {
	return this->_dmat.data()[i];
}
double kingghidorah::_mySparse::_at(int i, int j) {
	return this->_dmat(i, j);
}
int kingghidorah::_mySparse::cols() {
	return _mat[0].cols();
}
void kingghidorah::_mySparse::_resize(int n, int m) {
	this->_dmat.resize(n, m);
}
void kingghidorah::_mySparse::setmiddlecolum(Eigen::SparseMatrix<double> f, int start, int end) {
	this->_dmat.middleCols(start, end - start) = f;
}
void kingghidorah::_mySparse::permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm)
{
	for (int ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii] * perm.transpose();
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm,bool sparse,bool dense)
{

	int nn = _dmat.rows();
	int numthreads = 0;
	numthreads = omp_get_max_threads();
	int S = nn / numthreads / 2;
	auto pt = perm.transpose();
	if(sparse)
	if (_mat.size() >= 1)
	{
		if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
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

		}

#pragma omp parallel for
		for (int i = 0; i < nn; i += S)
		{
			int start = i;
			int end = i + S;
			if (end > nn)end = nn;
			_dmat.middleCols(start, end - start) = perm * _dmat.middleCols(start, end - start);
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
			if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
			{
				_mat[0].conservativeResize(M, M);
			}
		}
	}
	if (dense)
	{
		_dmat.conservativeResize(M, M);
	}
}
void kingghidorah::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic> perm2)
{
	int nn = _dmat.rows();
	int numthreads = omp_get_max_threads();

	int S = nn / numthreads / 2;
	auto pt = perm2.transpose();

	if (_mat.size() >= 1)
	{
		if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
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
		if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
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
	return _dmat.rows();
}
int kingghidorah::_mySparse::_cols() {
	return _dmat.cols();
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
	if (mat->_mat.size() > 0)
	{
		if(this->_mat.size()==0)
			this->_mat.resize(1);
		this->_mat[0] = mat->_mat[0];
	}
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
std::vector<Eigen::SparseMatrix<double>> e;
//std::vector<Eigen::MatrixXd> e;
int kingghidorah::_mySparse::ofAtA(_mySparse* A)
{
	int nn = A->cols();
	this->_dmat.setZero(nn, nn);
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
	this->_mat[0].setZero();
	for (int i = 0; i < mt; i++) {
		this->_dmat += e[i];
	}
	this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	//this->_dmat = this->_mat[0];
	//this->_dmat = this->_mat[0];
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

	this->_dmat.resize(A->cols(), A->cols());
	this->_dmat.setZero();
	for (int ii = 0; ii < _nt; ii++)
	{
		auto _ret = (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
		this->_dmat += _ret;
	}
	return ss.str();
}
std::string kingghidorah::_mySparse::info()
{
	return "viennacl has been abandoned";
}

void kingghidorah::_mySparse::_ofAtB_gpu(cuda* cuda, _mySparse* B, _mySparse* C)
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
}
void kingghidorah::_mySparse::_ofAtB(_mySparse* B, _mySparse* C)
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
}
Eigen::VectorXd kingghidorah::_mySparse::_ofBtAB(_mySparse* B, double* ptr, int N, _mySparse* C)
{
	static Eigen::MatrixXd D;
	//std::stringstream sss;
	//sss << B->_mat[0].rows() << "," << B->_mat[0].cols() << "," << this->_dmat.rows() << "," << this->_dmat.cols() << std::endl;

	//auto ff = B->_mat[0].transpose() * this->_dmat * B->_mat[0];
	//C->_dmat = ff;
	//return;
	int nn = B->_dmat.cols();
	int kk = this->_dmat.cols();

	C->_dmat.resize(nn, nn);
	C->_dmat.setZero();

	int mt = omp_get_max_threads();

	auto left = B->_mat[0].transpose();
	auto mid = this->_dmat;
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
		C->_dmat.middleCols(S, E - S) = D * right.middleCols(S, E - S);
	}
	Eigen::Map<Eigen::VectorXd> b(ptr, N);
	return D * b;
}

void kingghidorah::_mySparse::ofAtB(_mySparse* B)
{

	int nn = this->cols();
	int mm = B->cols();
	this->_dmat.resize(nn, mm);
	this->_dmat.setZero();
	int mt = omp_get_max_threads();

	_mt = mt;
	if(mt>e.size())
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
	this->_mat[0].resize(nn,mm);
	this->_mat[0].setZero();

	for (int i = 0; i < mt; i++) {
		this->_dmat += e[i];
	}
	/*for (int i = 0; i < mt; i++) {
		this->_dmat += e[i];
		this->_mat[0] += e[i];
	}*/
	this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
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
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		ret += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		offset += ee;
	}
	return ret;
}
Eigen::VectorXd kingghidorah::_mySparse::_Atb(double* ptr, int N)
{
	Eigen::Map<Eigen::VectorXd> b(ptr, N);
	return this->_dmat.transpose() * b;
}
void kingghidorah::_mySparse::merge()
{
	this->ofDat();
}
void kingghidorah::_mySparse::computeQR()
{
	Eigen::HouseholderQR<Eigen::MatrixXd> qr;
	qr.compute(_dmat);
}
void kingghidorah::_mySparse::computeLU()
{
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	lu.compute(_dmat);
}
void kingghidorah::_mySparse::computeLLT(Eigen::LLT<Eigen::MatrixXd>* _LLT)
{
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
	lu.compute(_dmat);
	Eigen::Map<Eigen::VectorXd> b(rhs, N);
	Eigen::VectorXd x(_dmat.cols());
	x.setZero();
	x = lu.solve(b);
	return x;
}

Eigen::VectorXd kingghidorah::_mySparse::_solve0_gpu_mg(kingghidorah::cuda* cuda, double* rhs, int N) {
	this->_freeze();
	Eigen::VectorXd x(N);
	x.setZero();


	if (!cuda->valid())
	{
		x(0) = 10;
		return x;
	}
	cusolverStatus_t status;
	auto solver = cuda->mgsolver();
	int nbGpus = cuda->count();
	int* deviceList = cuda->devicelist();

	const int NRHS = 1;

	const int IA = 1;
	const int JA = 1;
	const int T_A = std::max(1, N / nbGpus / 3);
	const int lda = N;


	const int IB = 1;
	const int JB = 1;
	const int T_B = std::max(1, NRHS / nbGpus / 3);
	const int ldb = N;


	int  info = 0;

	cudaLibMgMatrixDesc_t descrA;
	cudaLibMgMatrixDesc_t descrB;
	cudaLibMgGrid_t gridA;
	cudaLibMgGrid_t gridB;
	cusolverMgGridMapping_t mapping = CUDALIBMG_GRID_MAPPING_COL_MAJOR;

	//double** array_d_A=0;
	//double** array_d_B=0;

	int64_t lwork_potrf = 0;
	int64_t lwork_potrs = 0;
	int64_t lwork = 0;
	//double** array_d_work = NULL;
	printf("step 1\n");


	status = cusolverMgCreateDeviceGrid(&gridA, 1, nbGpus, deviceList, mapping);
	if (CUSOLVER_STATUS_SUCCESS != status)
	{
		x[0] = 13;
		return x;
	}
	status = cusolverMgCreateDeviceGrid(&gridB, 1, nbGpus, deviceList, mapping);
	if (CUSOLVER_STATUS_SUCCESS != status)
	{
		x[0] = 14;
		return x;
	}
	printf("step 2\n");

	status = cusolverMgCreateMatrixDesc(
		&descrA,
		N,
		N,
		N,
		T_A,
		CUDA_R_64F,
		gridA);
	if (CUSOLVER_STATUS_SUCCESS != status)
	{
		x[0] = 15;
		return x;
	}
	status = cusolverMgCreateMatrixDesc(
		&descrB,
		N,
		1,
		N,
		T_B,
		CUDA_R_64F,
		gridB);
	printf("step 3\n");
	if (CUSOLVER_STATUS_SUCCESS != status)
	{
		x[0] = 16;
		return x;
	}
	double** array_d_A = cuda->array_d_A();
	double** array_d_B = cuda->array_d_B();

	createMat<double>(
		nbGpus,
		deviceList,
		N,
		T_A,
		lda,
		array_d_A
		);
	createMat<double>(
		nbGpus,
		deviceList,
		1,
		T_B,
		ldb,
		array_d_B
		);
	printf("step 3\n");
	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		N,

		this->_dmat.data(),
		lda,

		N,
		T_A,
		lda,
		array_d_A,
		IA,
		JA
		);
	printf("step 4\n");

	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		NRHS,

		rhs,
		ldb,

		1,
		T_B,
		ldb,
		array_d_B,
		IB,
		JB
		);

	printf("step 5\n");

	status = cusolverMgPotrf_bufferSize(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA, //base-1
		JA,//base-1
		descrA,
		CUDA_R_64F,
		&lwork_potrf);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	printf("step 6\n");


	status = cusolverMgPotrs_bufferSize(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		NRHS,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		(void**)array_d_B,
		IB,
		JB,
		descrB,
		CUDA_R_64F,
		&lwork_potrs);
	assert(CUSOLVER_STATUS_SUCCESS == status);
	printf("step 7\n");

	lwork = (lwork_potrf > lwork_potrs) ? lwork_potrf : lwork_potrs;
	printf("\tallocate device workspace, lwork = %lld \n", (long long)lwork);
	printf("step 8\n");

	double** array_d_work = cuda->array_d_work();//(double**)malloc(sizeof(double*)* nbGpus);
	//assert(NULL != array_d_work);
	printf("step 9\n");


	workspaceAlloc(
		nbGpus,
		deviceList,
		sizeof(double) * lwork,
		(void**)array_d_work
	);
	printf("step 10\n");

	auto cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	printf("step 11\n");
	auto now = std::chrono::high_resolution_clock::now();
	status = cusolverMgPotrf(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		(void**)array_d_work,
		lwork,
		&info
	);
	printf("step 12\n");
	//if(CUSOLVER_STATUS_SUCCESS != status)assert(CUSOLVER_STATUS_SUCCESS);
	assert(CUSOLVER_STATUS_SUCCESS == status);
	cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	assert(0 == info);
	printf("step 13\n");



	status = cusolverMgPotrs(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		NRHS,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		(void**)array_d_B,
		IB,
		JB,
		descrB,
		CUDA_R_64F,
		(void**)array_d_work,
		lwork,
		&info
	);
	printf("step 14\n");
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Mg" << duration.count() << "ms" << std::endl;

	assert(CUSOLVER_STATUS_SUCCESS == status);
	cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	assert(0 == info);

	printf("step 15\n");
	for (int i = 0; i < 4; i++)printf("%.3f\n", x(i));
	memcpyD2H<double>(
		nbGpus,
		deviceList,
		N,
		1,

		1,
		T_B,
		ldb,
		array_d_B,
		IB,
		JB,

		x.data(),
		ldb
		);
	printf("step 16\n");

	for (int i = 0; i < 4; i++)printf("%.3f\n", x(i));
	destroyMat(
		nbGpus,
		deviceList,
		N,
		T_A,
		(void**)array_d_A);
	printf("step 17\n");

	destroyMat(
		nbGpus,
		deviceList,
		1,
		T_B,
		(void**)array_d_B);
	printf("step 18\n");

	workspaceFree(nbGpus, deviceList, (void**)array_d_work);
	printf("step 19\n");


	//if (NULL != array_d_A) free(array_d_A);
	//if (NULL != array_d_B) free(array_d_B);

	//if (NULL != array_d_work) free(array_d_work);



	return x;
}
Eigen::VectorXd kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, double* rhs, int N, int device) {
	Eigen::VectorXd x(N);
	this->_freeze();

	if (!cuda->valid())
	{
		x(0) = 10;
		return x;
	}
	cudaSetDevice(device);
	auto solver = cuda->solver(device);
	auto blas = cuda->blas(device);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	x.setZero();


	double* gpu_rhs = cuda->work_rhs(device);
	double* gpu_matrix = cuda->work_M(device);
	cudaMemcpy(gpu_matrix, this->_dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_rhs, rhs, N * sizeof(double), cudaMemcpyHostToDevice);

	int work_size = 0;
	int* devInfo_on_gpu = cuda->info(0);
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
		x(0) = 22;
		return x;
	}


	cusolverDnDpotrs(solver, CUBLAS_FILL_MODE_LOWER, N, 1, gpu_matrix, N, gpu_rhs, N, devInfo_on_gpu);
	auto end = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
	std::cout << "Dn:" << duration.count() << "ms" << std::endl;

	cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);

	if (devInfo_on_cpu != 0) {
		x(0) = 24;
		return x;
	}

	cudaMemcpy(x.data(), gpu_rhs, sizeof(double) * N, cudaMemcpyDeviceToHost);

	//cudaFree(work);

	//cudaFree(devInfo_on_gpu);

	cudaDeviceSynchronize();
	return x;
}
Eigen::MatrixXd kingghidorah::_mySparse::_solve0(_myLLT* LLT, _mySparse* mat)
{
	//this function assumes that LLT decomposition has been done already
	Eigen::MatrixXd ret(mat->_dmat.cols(), this->_dmat.rows());
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
		ret.middleRows(start, end - start) = LLT->LLT->solve(mat->_dmat.middleCols(start, end - start)).transpose();

	}
	return ret.transpose();
}
void kingghidorah::_mySparse::_solve0_gpu_mg(kingghidorah::cuda* cuda, _mySparse* mat, _mySparse* ret)
{
	int nn = mat->_dmat.cols();
	int N = this->_dmat.cols();
	int NRHS = nn;
	//Eigen::MatrixXd x(N,nn);
	if (ret->_dmat.cols() != nn || ret->_dmat.rows() != N)
	{
		ret->_dmat.resize(N, nn);
	}
	if (!cuda->valid())return;
	ret->_dmat.setZero();


	cusolverMgHandle_t solver = cuda->mgsolver();
	//auto status = cusolverMgCreate(&solver);

	int nbGpus = cuda->count();
	//std::vector<int> _deviceList(nbGpus);
	//for (int i = 0; i < nbGpus; i++)_deviceList[i] = i;
	int* deviceList = cuda->devicelist();
	enablePeerAccess(nbGpus, deviceList);

	const int IA = 1;
	const int JA = 1;
	const int T_A = N;// / nbGpus;
	const int lda = N;


	const int IB = 1;
	const int JB = 1;
	const int T_B = NRHS / nbGpus;
	const int ldb = N;


	int  info = 0;

	cudaLibMgMatrixDesc_t descrA;
	cudaLibMgMatrixDesc_t descrB;
	cudaLibMgGrid_t gridA;
	cudaLibMgGrid_t gridB;
	cusolverMgGridMapping_t mapping = CUDALIBMG_GRID_MAPPING_COL_MAJOR;

	double** array_d_A = cuda->array_d_A();
	double** array_d_B = cuda->array_d_B();

	int64_t lwork_potrf = 0;
	int64_t lwork_potrs = 0;
	int64_t lwork = 0;
	double** array_d_work = cuda->array_d_work();
	auto status = cusolverMgDeviceSelect(
		solver,
		nbGpus,
		deviceList);
	assert(CUSOLVER_STATUS_SUCCESS == status);
	/*int canAccessPeer = 0;
	if (nbGpus > 1)
	{
		assert(0 == enablePeerAccess(nbGpus, deviceList));
	}*/
	/*if (nbGpus == 4)
	{
		status = cusolverMgCreateDeviceGrid(&gridA, 2, 2, deviceList, mapping);
		assert(CUSOLVER_STATUS_SUCCESS == status);
		status = cusolverMgCreateDeviceGrid(&gridB, 2,2, deviceList, mapping);
		assert(CUSOLVER_STATUS_SUCCESS == status);
	}
	else */ {
		status = cusolverMgCreateDeviceGrid(&gridA, 1, nbGpus, deviceList, mapping);
		assert(CUSOLVER_STATUS_SUCCESS == status);
		status = cusolverMgCreateDeviceGrid(&gridB, 1, nbGpus, deviceList, mapping);
		assert(CUSOLVER_STATUS_SUCCESS == status);
	}
	status = cusolverMgCreateMatrixDesc(
		&descrA,
		N,
		N,
		N,
		T_A,
		CUDA_R_64F,
		gridA);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	status = cusolverMgCreateMatrixDesc(
		&descrB,
		N,
		NRHS,
		N,
		T_B,
		CUDA_R_64F,
		gridB);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	/*array_d_A = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_A);
	array_d_B = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_B);
	*/
	createMat<double>(
		nbGpus,
		deviceList,
		N,
		T_A,
		lda,
		array_d_A
		);
	createMat<double>(
		nbGpus,
		deviceList,
		NRHS,
		T_B,
		ldb,
		array_d_B
		);
	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		N,

		this->_dmat.data(),
		lda,

		N,
		T_A,
		lda,
		array_d_A,
		IA,
		JA
		);
	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		NRHS,

		mat->_dmat.data(),
		ldb,

		NRHS,
		T_B,
		ldb,
		array_d_B,
		IB,
		JB
		);


	status = cusolverMgPotrf_bufferSize(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		&lwork_potrf);
	assert(CUSOLVER_STATUS_SUCCESS == status);



	status = cusolverMgPotrs_bufferSize(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		NRHS,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		(void**)array_d_B,
		IB,
		JB,
		descrB,
		CUDA_R_64F,
		&lwork_potrs);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	lwork = (lwork_potrf > lwork_potrs) ? lwork_potrf : lwork_potrs;
	printf("\tallocate device workspace, lwork = %lld \n", (long long)lwork);

	//array_d_work = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_work);

	workspaceAlloc(
		nbGpus,
		deviceList,
		sizeof(double) * lwork,
		(void**)array_d_work
	);
	auto cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	auto now = std::chrono::high_resolution_clock::now();
	status = cusolverMgPotrf(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		(void**)array_d_work,
		lwork,
		&info
	);
	assert(CUSOLVER_STATUS_SUCCESS == status);
	cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	assert(0 == info);



	status = cusolverMgPotrs(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		NRHS,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		(void**)array_d_B,
		IB,
		JB,
		descrB,
		CUDA_R_64F,
		(void**)array_d_work,
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
		NRHS,

		NRHS,
		T_B,
		ldb,
		array_d_B,
		IB,
		JB,

		ret->_dmat.data(),
		ldb
		);

	destroyMat(
		nbGpus,
		deviceList,
		N,
		T_A,
		(void**)array_d_A);
	destroyMat(
		nbGpus,
		deviceList,
		1,
		T_B,
		(void**)array_d_B);

	workspaceFree(nbGpus, deviceList, (void**)array_d_work);


	//if (NULL != array_d_A) free(array_d_A);
	//if (NULL != array_d_B) free(array_d_B);
	//if (NULL != array_d_work) free(array_d_work);

}

void kingghidorah::_mySparse::_solveI_gpu_mg(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();

	int N = this->_dmat.cols();
	//Eigen::MatrixXd x(N,nn);
	if (ret->_dmat.cols() != N || ret->_dmat.rows() != N)
	{
		ret->_dmat.resize(N, N);
	}
	if (!cuda->valid())return;
	ret->_dmat.setZero();


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

	cudaLibMgMatrixDesc_t descrA;

	cudaLibMgGrid_t gridA;

	cusolverMgGridMapping_t mapping = CUDALIBMG_GRID_MAPPING_COL_MAJOR;

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
	auto status = cusolverMgCreateDeviceGrid(&gridA, 1, nbGpus, deviceList, mapping);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	status = cusolverMgCreateMatrixDesc(
		&descrA,
		N,
		N,
		N,
		T_A,
		CUDA_R_64F,
		gridA);
	assert(CUSOLVER_STATUS_SUCCESS == status);



	/*array_d_A = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_A);
	array_d_B = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_B);
	*/
	if (cuda->prevN < N || cuda->prevT_A < T_A)
	{
		if (cuda->prevN != 0)
			destroyMat(
				nbGpus,
				deviceList,
				cuda->prevN,
				cuda->prevT_A,
				(void**)array_d_A);

		createMat<double>(
			nbGpus,
			deviceList,
			N,
			T_A,
			lda,
			array_d_A
			);

		cuda->prevN = N;
		cuda->prevT_A = T_A;
	}
	memcpyH2D<double>(
		nbGpus,
		deviceList,
		N,
		N,

		this->_dmat.data(),
		lda,

		N,
		T_A,
		lda,
		array_d_A,
		IA,
		JA
		);


	status = cusolverMgPotrf_bufferSize(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		&lwork_potrf);
	assert(CUSOLVER_STATUS_SUCCESS == status);



	status = cusolverMgPotri_bufferSize(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		&lwork_potri);
	assert(CUSOLVER_STATUS_SUCCESS == status);

	lwork = (lwork_potrf > lwork_potri) ? lwork_potrf : lwork_potri;
	printf("\tallocate device workspace, lwork = %lld \n", (long long)lwork);
	double** array_d_work = cuda->array_d_work();
	//array_d_work = (double**)malloc(sizeof(double*) * nbGpus);
	assert(NULL != array_d_work);
	if (cuda->prevwn < lwork)
	{
		if (cuda->prevwn != 0)
			workspaceFree(nbGpus, deviceList, (void**)array_d_work);
		workspaceAlloc(
			nbGpus,
			deviceList,
			sizeof(double) * lwork,
			(void**)array_d_work
		);
		cuda->prevwn = lwork;
	}
	auto cudaStat = cudaDeviceSynchronize();
	assert(cudaSuccess == cudaStat);
	auto now = std::chrono::high_resolution_clock::now();
	status = cusolverMgPotrf(
		solver,
		CUBLAS_FILL_MODE_LOWER,
		N,
		(void**)array_d_A,
		IA,
		JA,
		descrA,
		CUDA_R_64F,
		(void**)array_d_work,
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
		(void**)array_d_work,
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
Eigen::MatrixXd I;
void kingghidorah::_mySparse::_solveI(_mySparse* ret)
{
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt;
	llt.compute(this->_mat[0]);
	int nn = this->_mat[0].rows();
	ret->_dmat.resize(nn, nn);
	I.resize(nn, nn);
	I.setIdentity();
	int mt = omp_get_max_threads();
	int ee = mt * 4;
#pragma omp parallel for
	for (int i = 0; i < ee; i++)
	{
		int S = i * nn / ee;
		int E = (i + 1) * nn / ee;
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


void kingghidorah::_mySparse::_solveI_gpu_omp(kingghidorah::cuda* cuda, _mySparse* ret)
{
	//static Eigen::MatrixXd x;
	this->_freeze();
	int N = this->_dmat.cols();
	if (identity == 0)
	{
		cudaMallocHost(&identity,sizeof(double) * N * N);
		previdentiyN = N;
		memset(identity, 0, sizeof(double) * N * N);
		for (int i = 0; i < N; i++)identity[i + i * N] = 1;
	}
	if (N > previdentiyN)
	{
		if (identity != 0)cudaFreeHost(identity);
		cudaMallocHost(&identity, sizeof(double) * N * N);
		previdentiyN = N;
		memset(identity, 0, sizeof(double) * N * N);
		for (int i = 0; i < N; i++)identity[i + i * N] = 1;
	}
	//I.setIdentity(N, N);
	ret->_dmat.resize(N, N);

	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest());
	auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->_dmat.setZero();
	int job = 0;

	cudaSetDevice(cuda->fastest());

	//double* tmp = cuda->work_M2();
	//memcpy(tmp, this->_dmat.data(), N * N * sizeof(double));


	double* gpu_matrix = cuda->work_M(cuda->fastest());
	cudaMemcpy(gpu_matrix, this->_dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	int* devInfo_on_gpu = cuda->info(cuda->fastest());
	//cudaMalloc(&devInfo_on_gpu, sizeof(int));
	int work_size = 0;
	//int work_size1 = 0;
	//int work_size2 = 0;

	// --- CUDA CHOLESKY initialization
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);
	//cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size2);
	//work_size = work_size1;// std::max(work_size1, work_size2);
	// --- CUDA POTRF execution	
	//cudaMalloc(&work, work_size * sizeof(double));
	double* work = cuda->work(work_size, cuda->fastest());

	cudaMemset(work, 0, work_size * sizeof(double));
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
	if (cuda->count() > 1)
	{
		cudaSetDevice(cuda->fastest());
		for (int i = 0; i < cuda->count(); i++)
		{
			if (i != cuda->fastest())
			{
				cudaMemcpy(cuda->work_M(i), cuda->work_M(cuda->fastest()), sizeof(double) * N * N,cudaMemcpyDeviceToDevice);
			}
		}
	}

	int devInfo_on_cpu = 0;

#pragma omp parallel for
	for (int i = 0; i < cuda->count(); i++)
	{
		cudaSetDevice(i);
		int _S = i * N / cuda->count();
		int _E = (i + 1) * N / cuda->count();
		double* gpu_rhs = cuda->work_rhs(i);
		double* gpu_matrix = cuda->work_M(i);
		cudaStream_t stream[3];
		cudaStreamCreate(&stream[0]);
		cudaStreamCreate(&stream[1]);
		cudaStreamCreate(&stream[2]);
		for (int j = 0; j < 2; j++)
		{
			int S = _S + (_E - _S) * j / 2;
			int E = _S + (_E - _S) * (j+1) / 2;

			cudaMemcpy(gpu_rhs + S * N, identity + S * N, (E - S) * N * sizeof(double), cudaMemcpyHostToDevice);
			cusolverDnDpotrs(solver, CUBLAS_FILL_MODE_LOWER, N, E - S, gpu_matrix, N, gpu_rhs + S * N, N, devInfo_on_gpu);
			cudaMemcpy(ret->_dmat.data() + N * S, gpu_rhs + N * S, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost);
		}
		cudaStreamDestroy(stream[0]);
		cudaStreamDestroy(stream[1]);
		cudaStreamDestroy(stream[2]);

	}
	


}
void kingghidorah::_mySparse::_solveI_gpu(kingghidorah::cuda* cuda, _mySparse* ret)
{
	this->_freeze();

	int N = this->_dmat.cols();

	ret->_dmat.resize(N, N);

	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest());
	auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->_dmat.setZero();
	int job = 0;

	cudaSetDevice(cuda->fastest());

	double* tmp = cuda->work_M2();
	memcpy(tmp, this->_dmat.data(), N * N * sizeof(double));


	double* gpu_matrix = cuda->work_M(cuda->fastest());
	cudaMemcpy(gpu_matrix, tmp, N * N * sizeof(double), cudaMemcpyHostToDevice);
	int* devInfo_on_gpu = cuda->info(cuda->fastest());
	//cudaMalloc(&devInfo_on_gpu, sizeof(int));
	int work_size = 0;
	int work_size1 = 0;
	int work_size2 = 0;

	// --- CUDA CHOLESKY initialization
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size1);
	cusolverDnDpotri_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size2);
	work_size = std::max(work_size1, work_size2);
	// --- CUDA POTRF execution	
	//cudaMalloc(&work, work_size * sizeof(double));
	double* work = cuda->work(work_size, cuda->fastest());

	cudaMemset(work, 0, work_size1 * sizeof(double));
	auto start = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size1, devInfo_on_gpu);

	auto end = high_resolution_clock::now();
	auto duration = end - start;
	std::chrono::milliseconds d = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
	std::cout << "regular:Dpotrf:" << d.count() << "ms" << std::endl;

	int devInfo_on_cpu = 0;
	/*cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);

	if (0 != devInfo_on_cpu) {
		ret->_dmat(0, 0) = 34;
		return;
	}*/

	cudaMemset(work, 0, work_size2 * sizeof(double));
	start = std::chrono::high_resolution_clock::now();
	cusolverDnDpotri(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size2, devInfo_on_gpu);
	end = high_resolution_clock::now();
	duration = end - start;
	d = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
	std::cout << "regular:Dpotri:" << d.count() << "ms" << std::endl;

	/*cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
	if (devInfo_on_cpu != 0) {
		ret->_dmat(0, 0) = 12;
		return;
	}*/
	cudaMemcpy(ret->_dmat.data(), gpu_matrix, N * N * sizeof(double), cudaMemcpyDeviceToHost);
	//cudaFree(work);
	//cudaFree(devInfo_on_gpu);
	//ret->_dmat.triangularView<Eigen::Upper>() = ret->_dmat.triangularView<Eigen::Lower>();
	//memcpy(tmp, ret->_dmat.data(), N * N * sizeof(double));
#pragma omp parallel for
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < i; j++)
		{
			ret->_dmat(j, i) = ret->_dmat(i, j);
		}
	}
}
void kingghidorah::_mySparse::_solve0_gpu(kingghidorah::cuda* cuda, _mySparse* mat, _mySparse* ret)
{
	int nn = mat->_dmat.cols();
	int N = this->_dmat.cols();
	//Eigen::MatrixXd x(N,nn);
	if (ret->_dmat.cols() != nn || ret->_dmat.rows() != N)
	{
		ret->_dmat.resize(N, nn);
	}
	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest());
	auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->_dmat.setZero();
	int job = 0;

	cudaSetDevice(cuda->fastest());
	double* gpu_matrix = cuda->work_M(cuda->fastest());
	cudaMemcpy(gpu_matrix, this->_dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
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
	cudaMemcpy(this->_dmat.data(), gpu_matrix, N * N * sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(work);
	cudaFree(devInfo_on_gpu);
	bool exit = false;

	int S = nn / cuda->count() / 5;

#pragma omp parallel for
	for (int i = 0; i < cuda->count(); i++) {
		cudaSetDevice(i);
		auto solver = cuda->solver(i);
		auto blas = cuda->blas(i);
		double* _gpu_matrix = cuda->work_M(i);
		double* gpu_rhs = cuda->work_rhs(i);
		cudaMemcpy(_gpu_matrix, this->_dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
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

Eigen::VectorXd kingghidorah::_mySparse::_solve0(double* rhs, int N) {
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> LLT;
	LLT.compute(_mat[0]);
	Eigen::Map<Eigen::VectorXd> b(rhs, N);
	Eigen::VectorXd x(_mat[0].rows());
	x.setZero();
	x = LLT.solve(b);
	return x;
}
Eigen::VectorXd kingghidorah::_mySparse::__solve0(double* rhs, int N) {
	Eigen::LLT<Eigen::MatrixXd> LLT;
	LLT.compute(_dmat);
	Eigen::Map<Eigen::VectorXd> b(rhs, N);
	Eigen::VectorXd x(_dmat.rows());
	x.setZero();
	x = LLT.solve(b);
	return x;
}
Eigen::MatrixXd kingghidorah::_mySparse::inv() {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	lu.compute(_dmat);
	Eigen::MatrixXd I(_dmat.rows(), _dmat.cols());
	I.setIdentity();
	return lu.solve(I);
}

Eigen::MatrixXd kingghidorah::_mySparse::solve0(_mySparse* rhs)
{
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//_mat.makeCompressed();
	lu.compute(_dmat);
	Eigen::MatrixXd _x(_dmat.rows(), rhs->_dmat.cols());
	_x = lu.solve(rhs->_dmat);

	return _x;// .sparseView(0.000000000000001, 1.0);
}


void kingghidorah::_mySparse::minus(_mySparse* m) {
	this->_freeze();

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
void kingghidorah::_mySparse::addsmallidentity(double salt) {
	id.resize(this->_dmat.rows(), this->_dmat.cols());
	id.setIdentity();
	this->_dmat += (id * salt);
	this->_mat[0] += (id * salt);
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