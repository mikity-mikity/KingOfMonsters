#pragma once
#pragma once
#include <cmath>
#include<vector>
#include <cstring>
#include <string>
#include<iostream>
#include"eigenwrapper.h"
#include <omp.h> 
using std::vector;
using std::string;
using namespace System;
using namespace System::Threading::Tasks;
#define EIGEN_DONT_ALIGN
namespace kingghidorah {
	public ref class myDoubleArray {
	public:
		//double* _arr = 0;
		_myDoubleArray *_arr=0;
		int _N = 0;
	public:
		void setzero(int S, int N)
		{
			_arr->__v.middleRows(S, N).setZero();
		}
		myDoubleArray(int N)
		{
			_arr = new _myDoubleArray();
			//_arr = new double[N];
			_N = N;
			_arr->__v.resize(N);
			_arr->__v.setZero();
		}
		/*inline double* data()
		{
			return _arr;
		}*/
		inline int size() {
			return _N;
		}
		void copyfrom(myDoubleArray^ arr, int N)
		{
			_arr->__v = arr->_arr->__v;
			//System::Runtime::InteropServices::Marshal::Copy( arr,0, (System::IntPtr)_arr->data(), N);
		}
		void minus() {
			_arr->__v.noalias() = -(_arr->__v);
		}
		void plus(myDoubleArray^ a, double sc) {
			_arr->__v.noalias() += _arr->__v * sc;
		}

		!myDoubleArray()
		{
			if (_arr != 0)
			{
				delete _arr;
			}
			_arr = 0;
			_N = 0;
		}
		~myDoubleArray()
		{
			if (_arr != 0) {
				delete _arr;
			}
			_arr = 0;
			_N = 0;
		}
		void set(int i, double val)
		{
			(_arr->__v)(i) = val;
		}
		double at(int i)
		{
			return _arr->__v(i);
		}
		void scale(double sc)
		{
			_arr->__v *= sc;
		}
		void resize(int N)
		{
			_arr->__v.conservativeResize(N);
			if (N > _N)
			{
				_arr->__v.middleRows(_N, N - _N).setZero();
			}
			_N = N;
		}
		void reset(int N)
		{
			_N = N;
			_arr->__v.resize(_N);
			_arr->__v.setZero();
		}
		void minus(myDoubleArray^ b)
		{
			_arr->__v.noalias() -= b->_arr->__v;
		}
		double L2Norm() {
			return _arr->__v.norm();
		}
	};
	public ref class myIntArray {
	public:
		int* _arr = 0;
		int _N = 0;
	public:
		void setzero(int S, int N)
		{
			memset(_arr + S, 0, sizeof(int) * N);
		}
		myIntArray(int N)
		{
			_arr = new int[N];
			_N = N;
		}
		inline int* data()
		{
			return _arr;
		}
		inline int size()
		{
			return _N;
		}
		!myIntArray()
		{
			if (_arr != 0)
			{
				delete[] _arr;
			}
			_arr = 0;
			_N = 0;
		}
		~myIntArray()
		{
			if (_arr != 0) {
				delete[] _arr;
			}
			_arr = 0;
			_N = 0;
		}
		void set(int i, int val)
		{
			_arr[i] = val;
		}
	};
	public ref class myCuda {
	private:
		cuda* _cuda = 0;
	public:
		static void disable() {
			cuda::disable();
		}
		myCuda(int N) {
			_cuda = new kingghidorah::cuda(N);
		}
		!myCuda() {
			if (_cuda != 0)
				delete _cuda;
			_cuda = 0;
		}
		~myCuda() {
			if (_cuda != 0)
				delete _cuda;
			_cuda = 0;
		}
		void dispose() {
			if (_cuda != 0)
			{
				_cuda->dispose();
				delete _cuda;
			}
			_cuda = 0;
		}
		kingghidorah::cuda* cuda() {
			return _cuda;
		}
		bool isValid() {
			return _cuda->valid();
		}
		System::String^ device_name() {
			std::string s = _cuda->device_name();
			System::String^ ret = gcnew System::String(s.c_str(), 0, s.size());
			return ret;
		}
		int fastest() {
			return _cuda->fastest();
		}
		int count() {
			return _cuda->count();
		}
		bool canpeeraccess(int i, int j)
		{
			return _cuda->canpeeraccess(i, j);
		}
	};
	public ref class myLLT {
	public:
		_myLLT* LLT;
		myLLT()
		{
			LLT = new _myLLT();
		}
		!myLLT() {
			if (LLT != 0)
				delete LLT;
			LLT = 0;
		}
		~myLLT() {
			if (LLT != 0)
				delete LLT;
			LLT = 0;
		}
	};
	public ref class myPermutation {
	public:
		_myPermutation* p;
	public:
		myPermutation(array<int>^ index)
		{

			pin_ptr<int> ptr = &index[0];

			p = new _myPermutation(ptr, index->Length);

			ptr = nullptr;
		}
		!myPermutation() {
			if (p != 0)
				delete(p);
			p = 0;
		}
		~myPermutation() {
			if (p != 0)
				delete(p);
			p = 0;
		}
		void perm(array<double>^ vec) {
			pin_ptr<double> ptr = &vec[0];
			int N = vec->Length;
			Eigen::Map<Eigen::VectorXd> b(ptr, N);

			b.applyOnTheLeft(this->p->perm);
			//array<double>^ ret = gcnew array<double>(vec->Length);
			//pin_ptr<double> ptr2 = &ret[0];
			//Eigen::Map<Eigen::VectorXd>(ptr2, N) = b;
			//Eigen::VectorXd _b = b;
			//memcpy(ptr2, _b.data(), sizeof(double) * N);
			ptr = nullptr;
			//ptr2 = nullptr;
			//return ret;
		}
		/*void  permback(array<double>^ vec) {
			pin_ptr<double> ptr = &vec[0];
			int N = vec->Length;
			Eigen::Map<Eigen::VectorXd> b(ptr, N);
			//auto _ret = this->p->perm.transpose() * b;
			b.applyOnTheLeft(this->p->perm.transpose());

			//array<double>^ ret = gcnew array<double>(N);
			//pin_ptr<double> ptr2 = &ret[0];

			//Eigen::VectorXd _b = _ret;
			//memcpy(ptr2, _b.data(), sizeof(double) * N);
			//Eigen::Map<Eigen::VectorXd>(ptr2, N) = _ret;

			ptr = nullptr;
			//ptr2 = nullptr;
			//return ret;
		}*/
		void perm(myDoubleArray^ vec) {
			//double* ptr = vec->data();
			int N = vec->size();
			//Eigen::Map<Eigen::VectorXd> b(ptr, N);
			vec->_arr->__v.applyOnTheLeft(this->p->perm);
		}
		void permback(myDoubleArray^ vec) {
			//double* ptr = vec->data();
			//int N = vec->size();
			//Eigen::Map<Eigen::VectorXd> b(ptr, N);
			//b.applyOnTheLeft(this->p->perm.transpose());
			vec->_arr->__v.applyOnTheLeft(this->p->perm.transpose());
		}
	};
	public ref class mySparse {
	public:
		_mySparse* dat;
		mySparse() {
			dat = 0;
			dat = new _mySparse();
			dat->init(0, 0);
		}
		mySparse(int n, int m)
		{
			dat = 0;
			dat = new _mySparse();
			dat->init(n, m);
		}
		mySparse(mySparse^ m)
		{
			dat = 0;
			dat = new _mySparse();
			dat->init(m->rows(), m->cols());
			this->dat->OfDuplicate(m->dat);
			this->dat->copycoefffrom(m->dat);
		}
		void ofDuplicate(mySparse^ m)
		{
			dat->init(m->rows(), m->cols());
			this->dat->OfDuplicate(m->dat);
			this->dat->copycoefffrom(m->dat);
		}
		void _ofDuplicate(mySparse^ m)
		{
			dat->init(m->rows(), m->cols());
			this->dat->_OfDuplicate(m->dat);
		}
		void freeze(bool _do) {
			this->dat->freeze(_do);
		}
		mySparse^ Duplicate()
		{
			return gcnew mySparse(this);
		}
		~mySparse()
		{
			delete(dat);
			dat = 0;
		}
		!mySparse()
		{
			delete(dat);
			dat = 0;
		}
		void ofDat()
		{
			dat->ofDat();
		}
		void merge()
		{
			dat->merge();
		}

		void computeQR() {
			dat->computeQR();
		}
		void computeLU() {
			dat->computeLU();
		}
		int resize(int n, int m)
		{
			return dat->resize(n, m);
		}
		void reserve(int n) {
			dat->reserve(n);
		}
		void addmat(mySparse^ m)
		{
			dat->addmat(m->dat);
		}
		int cols() {
			return dat->cols();
		}
		int rows() {
			return dat->rows();
		}
		int _cols() {
			return dat->_cols();
		}
		int _rows() {
			return dat->_rows();
		}
		int __rows() {
			return dat->__rows();
		}
		void shrink(int M)
		{
			dat->shrink(M);
		}
		void permute(myPermutation^ p)
		{
			dat->permute(p->p->perm);
		}
		void _shrink(int M, bool sparse, bool dense)
		{
			dat->_shrink(M, sparse, dense);
		}
		void _permute(myPermutation^ p, bool sparse, bool dense)
		{
			dat->_permute(p->p->perm, sparse, dense);
		}
		void _shrink(int M, int N)
		{
			dat->_shrink(M, N);
		}
		void _permute(myPermutation^ p, myPermutation^ q)
		{
			dat->_permute(p->p->perm, q->p->perm);
		}
		int ofAtA(mySparse^ m, bool sparse) {
			return dat->ofAtA(m->dat, sparse);
		}
		System::String^ _ofAtA(mySparse^ m) {
			auto str = dat->_ofAtA(m->dat);
			auto ret = gcnew System::String(str.c_str());
			return ret;
		}
		void ofAtB(mySparse^ m, bool sparse) {
			dat->ofAtB(m->dat, sparse);
		}
		void _ofAtB(mySparse^ A, mySparse^ B)
		{
			A->dat->_ofAtB(B->dat, this->dat);
		}
		void _ofBtAB(mySparse^ A, mySparse^ B, myDoubleArray^ b, myDoubleArray^ ret)
		{
			//pin_ptr<double> ptr = &b[0];

			A->dat->_ofBtAB(B->dat, &b->_arr->__v, this->dat, &ret->_arr->__v);
			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			//ptr = nullptr;
			//return ret;
		}
		void addemptyrow(int ii) {
			dat->addemptyrow(ii);
		}
		System::String^ info() {
			auto str = this->dat->info();
			System::String^ s = gcnew System::String(str.data(), 0, str.size());
			return s;
		}
		void addrow(int ii, array<int>^ index, array<double>^ data, double sc, int N) {
			pin_ptr<int> ptr = &index[0];
			pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, ptr, dat, sc, N);
			ptr = nullptr;
			dat = nullptr;
		}
		void addrow(int ii, array<int>^ index, array<double>^ data, int shift, double sc, int N, bool add) {
			pin_ptr<int> ptr = &index[0];
			pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, ptr, dat, shift, sc, N, add);
			ptr = nullptr;
			dat = nullptr;
		}
		void addrow(int ii, myIntArray^ index, myDoubleArray^ data, double sc, int N) {
			//pin_ptr<int> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), sc, N);
			//ptr = nullptr;
			//dat = nullptr;
		}
		void addrow(int ii, myIntArray^ index, myDoubleArray^ data, int shift, double sc, int N, bool add) {
			//pin_ptr<int> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), shift, sc, N, add);
			//ptr = nullptr;
			//dat = nullptr;
		}
		array<double>^ Atb(array<double>^ b)
		{
			pin_ptr<double> ptr = &b[0];
			Eigen::VectorXd rhs = dat->Atb(ptr, b->Length);
			array<double>^ ret = gcnew array<double>(rhs.rows());
			System::Runtime::InteropServices::Marshal::Copy((IntPtr)rhs.data(), ret, 0, rhs.rows());
			ptr = nullptr;
			return ret;
		}
		void Atb(array<double>^ b, kingghidorah::myDoubleArray^ ret)
		{
			pin_ptr<double> ptr = &b[0];
			dat->Atb(ptr, b->Length, &ret->_arr->__v);
			//array<double>^ ret = gcnew array<double>(rhs.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)rhs.data(), ret, 0, rhs.rows());
			ptr = nullptr;
			//return ret;
		}
		array<double>^ _Atb(array<double>^ b)
		{
			pin_ptr<double> ptr = &b[0];
			Eigen::VectorXd rhs = dat->_Atb(ptr, b->Length);
			array<double>^ ret = gcnew array<double>(rhs.rows());
			System::Runtime::InteropServices::Marshal::Copy((IntPtr)rhs.data(), ret, 0, rhs.rows());
			ptr = nullptr;
			return ret;
		}
		void solve0(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->solve0(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}
		void _solve0(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->_solve0(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}
		void __solve0(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->__solve0(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}

		void _solve0_gpu(myCuda^ gpu, myDoubleArray^ rhs, myDoubleArray^ ret, int device) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->_solve0_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}

		mySparse^ solve0(mySparse^ rhs) {
			Eigen::LLT<Eigen::MatrixXd>* _LLT = new Eigen::LLT<Eigen::MatrixXd>();
			dat->computeLLT(_LLT);
			myLLT^ LLT = gcnew myLLT();
			LLT->LLT->LLT = _LLT;

			Eigen::MatrixXd _ret = this->dat->_solve0(LLT->LLT, rhs->dat);
			mySparse^ ret = gcnew mySparse(rhs->rows(), rhs->cols());

			return ret;
		}

		mySparse^ solve0_gpu(myCuda^ gpu, mySparse^ rhs, mySparse^ ret) {
			this->dat->_solve0_gpu(gpu->cuda(), rhs->dat, ret->dat);
			return ret;
		}
		mySparse^ solveI(mySparse^ ret) {
			this->dat->_solveI(ret->dat);
			return ret;
		}
		mySparse^ solveI_gpu(myCuda^ gpu, mySparse^ ret) {
			this->dat->_solveI_gpu(gpu->cuda(), ret->dat);
			return ret;
		}
		System::String^ solveI_gpu_omp(myCuda^ gpu, mySparse^ ret) {
			std::string ss = this->dat->_solveI_gpu_omp(gpu->cuda(), ret->dat);
			auto ee = gcnew System::String(ss.c_str());
			return ee;
		}
		System::String^ solveI_gpu_single(myCuda^ gpu, mySparse^ ret) {
			std::string ss = this->dat->_solveI_gpu_single(gpu->cuda(), ret->dat);
			auto ee = gcnew System::String(ss.c_str());
			return ee;
		}
		mySparse^ solveI_gpu_mg(myCuda^ gpu, mySparse^ ret) {
			this->dat->_solveI_gpu_mg(gpu->cuda(), ret->dat);
			return ret;
		}
		void freezecoeff()
		{
			dat->freezecoeff();
		}
		void addcoeff(double sc) {
			dat->addcoeff(sc);
		}
		double at_each(int i, int ii) {
			return dat->at(i, ii);
		}
		void begin_construct()
		{
			dat->begin_construct();
		}
		void end_construct(int c)
		{
			dat->end_construct(c);
		}
		int num_elem(int ii)
		{
			return dat->num_elem(ii);
		}
		double _at(int i) {
			return dat->_at(i);
		}
		double _at(int i, int j) {
			return dat->_at(i, j);
		}
		int nonZeros() {
			return dat->nonzeros();
		}
		void adddat(int i, int j, double val) {
			dat->adddat(i, j, val);
		}

		void minus(mySparse^ m)
		{
			dat->minus(m->dat);
		}
		void clearcoeff()
		{
			dat->clearcoeff();
		}
		double L2Norm(myDoubleArray^ a, myDoubleArray^ b) {
			double ret = dat->L2Norm(&a->_arr->__v, &b->_arr->__v);

			return ret;
		}
		array<double>^ vector(array<double>^ a) {
			pin_ptr<double> ptr1 = &a[0];
			auto _ret = dat->Vector(ptr1, a->Length);
			ptr1 = nullptr;
			array<double>^ ret = gcnew array<double>(_ret.rows());
			System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			return ret;
		}
		array<double>^ vector(myDoubleArray^ a) {
			auto _ret = dat->Vector(&a->_arr->__v);
			array<double>^ ret = gcnew array<double>(_ret.rows());
			System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			return ret;
		}
		void vector(myDoubleArray^ a, myDoubleArray^ b) {
			b->resize(dat->cols());
			dat->Vector(&a->_arr->__v, &b->_arr->__v);
		}
		void plus(mySparse^ m, double a, bool dense, bool sparse) {
			this->dat->plus(m->dat, a, dense, sparse);
		}
		void addsmallidentity(double salt, bool sparse, bool dense) {
			this->dat->addsmallidentity(salt, sparse, dense);
		}
		void Clear() {
			this->dat->Clear();
		}
		int numBlocks()
		{
			return this->dat->numBlocks();
		}
		static System::String^ _testopenmp()
		{
			auto _str = kingghidorah::_mySparse::_testopenmp();
			return gcnew System::String(_str.c_str());
		}
		static System::String^ testopenmp()
		{
			auto str = gcnew System::String("");
			int mt = omp_get_max_threads();
			str += "num threads:" + mt.ToString() + "\n";
#pragma omp parallel for
			for (int i = 0; i < 100; i++)
			{
				int ct = omp_get_thread_num();
#pragma omp critical
				{
					str += "current threads:" + ct.ToString() + "\n";
				}
				for (int t = 0; t < ct * 100; t++)
				{
					auto f = new double[10];
				}
			}
			return str;
		}
	};
}