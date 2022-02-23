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
using std::basic_string;
using namespace System;
using namespace System::Threading::Tasks;
//#define EIGEN_DONT_ALIGN
namespace KingOfMonsters {
	public ref class mySparseVector {
	public:
		_mySparseVector* _vec=0;
		mySparseVector()
		{
			this->_vec = new _mySparseVector();
		}
		mySparseVector(Int64 N)
		{
			this->_vec = new _mySparseVector();
			this->_vec->_vec.resize(N);
			this->_vec->_vec.setZero();
		}
		void plus(Int64 i, double val)
		{
			this->_vec->_vec.coeffRef(i) += val;
		}
		void Zeros() {
			this->_vec->_vec.setZero();
		}
		~mySparseVector()
		{
			if (_vec != 0)
			{
				delete _vec;
			}
			_vec = 0;
		}
		!mySparseVector()
		{
			if (_vec != 0)
			{
				delete _vec;
			}
			_vec = 0;
		}
	};
	public ref class myDoubleArray {
	public:
		//double* _arr = 0;
		_myDoubleArray *_arr=0;
		Int64 _N = 0;
	public:
		double dot(myDoubleArray^ v)
		{
			return this->_arr->__v.dot(v->_arr->__v);
		}
		myDoubleArray^ subVector(Int64 i, Int64 N)
		{
			myDoubleArray^ ret = gcnew myDoubleArray(N);
			ret->_arr->__v = this->_arr->__v.middleRows(i, N);
			return ret;
		}
		void setzero(Int64 S, Int64 N)
		{
			_arr->__v.middleRows(S, N).setZero();
		}
		myDoubleArray(Int64 N)
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
		inline Int64 size() {
			return _N;
		}
		void copyfrom(myDoubleArray^ arr, Int64 N)
		{
			_arr->__v = arr->_arr->__v;
			//System::Runtime::InteropServices::Marshal::Copy( arr,0, (System::IntPtr)_arr->data(), N);
		}
		void minus() {
			_arr->__v = -(_arr->__v);
		}
		void plus(myDoubleArray^ a, double sc) {
			_arr->__v += a->_arr->__v * sc;
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
		void ofJoin(myDoubleArray^ a, myDoubleArray^ b)
		{
			this->_arr->__v.topRows(a->_arr->__v.size()) = a->_arr->__v;
			this->_arr->__v.bottomRows(b->_arr->__v.size()) = b->_arr->__v;
		}
		void set(Int64 i, double val)
		{
			(_arr->__v)(i) = val;
		}
		void ofAplusB(myDoubleArray^ A, myDoubleArray^ B)
		{
			this->_arr ->__v = A->_arr->__v + B->_arr->__v;
		}
		void plus(Int64 i, double val)
		{
			(_arr->__v)(i) += val;
		}
		void set(myDoubleArray^ f)
		{
			_arr->__v = f->_arr->__v;
			//System::Runtime::InteropServices::Marshal::Copy( arr,0,(IntPtr) _arr->__v.data(), arr->Length);

		}
		double at(Int64 i)
		{
			return _arr->__v(i);
		}
		void scale(double sc)
		{
			_arr->__v *= sc;
		}
		void resize(Int64 N)
		{	
 			Eigen::VectorXd v(N);
			v.setZero();
			if (_N < N)
			{
				v.middleRows(0, _N) = this->_arr->__v;
				this->_arr->__v.resize(N);
				this->_arr->__v = v;
			}
			else if(N<_N){
				v = this->_arr->__v.middleRows(0, N);
				this->_arr->__v.resize(N);
				this->_arr->__v = v;
			}
			//_arr->__v.conservativeResize(N);
			/*if (N > _N)
			{
				_arr->__v.middleRows(_N, N - _N).setZero();
			}*/
			_N = N;
		}
		void reset(Int64 N)
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
		Int64* _arr = 0;
		Int64 _N = 0;
	public:
		void setzero(Int64 S, Int64 N)
		{
			memset(_arr + S, 0, sizeof(Int64) * N);
		}
		myIntArray(Int64 N)
		{
			_arr = new Int64[N];
			_N = N;
		}
		inline Int64* data()
		{
			return _arr;
		}
		inline Int64 size()
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
		void set(Int64 i, Int64 val)
		{
			_arr[i] = val;
		}
		Int64 at(Int64 i)
		{
			return _arr[i];
		}
	};
	public ref class myCuda {
	private:
		cuda* _cuda = 0;
	public:
		static void disable() {
			cuda::disable();
		}
		myCuda(Int64 N) {
			_cuda = new KingOfMonsters::cuda(N);
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
		KingOfMonsters::cuda* cuda() {
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
		Int64 fastest() {
			return _cuda->fastest();
		}
		Int64 count() {
			return _cuda->count();
		}
		/*bool canpeeraccess(Int64 i, Int64 j)
		{
			return _cuda->canpeeraccess(i, j);
		}*/
	};
	public ref class myLLT {
	public:
		_myLLT* LLT=0;
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
		myPermutation(array<Int64>^ index)
		{

			pin_ptr<Int64> ptr = &index[0];

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
			Int64 N = vec->Length;
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
			Int64 N = vec->Length;
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
			Int64 N = vec->size();
			//Eigen::Map<Eigen::VectorXd> b(ptr, N);
			vec->_arr->__v.applyOnTheLeft(this->p->perm);
		}
		void permback(myDoubleArray^ vec) {
			//double* ptr = vec->data();
			//Int64 N = vec->size();
			//Eigen::Map<Eigen::VectorXd> b(ptr, N);
			//b.applyOnTheLeft(this->p->perm.transpose());
			vec->_arr->__v.applyOnTheLeft(this->p->perm.transpose());
		}
	};
	
	public ref class mySparse {
	public:
		_mySparse* dat = 0;
		mySparse^ computeKernel()
		{
			Eigen::SparseMatrix<double> ff = (this->dat->_mat[0] * this->dat->_mat[0].transpose());
			Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
			lu.compute(ff);
			auto I = Eigen::SparseMatrix<double>(ff.rows(), ff.cols());
			Eigen::MatrixXd ret=lu.solve(I);
			auto II = Eigen::SparseMatrix<double>(this->dat->_mat[0].cols(), this->dat->_mat[0].cols());

			mySparse^ kernel = gcnew mySparse();
			kernel->dat->_dmat= II - this->dat->_mat[0].transpose() * ret * this->dat->_mat[0];
			return kernel;
		}
		void setFromList(System::Collections::Generic::List < System::Collections::Generic::List<System::Tuple<Int64, Int64>^>^>^ tt)
		{
			std::vector<Eigen::Triplet<double>> dat;
			long count = 0;
			for (int i = 0; i < tt->Count; i++)
			{
				count += tt[i]->Count;
			}
			dat.resize(count);
			int current = 0;
			for (int i = 0; i < tt->Count; i++)
			{
				auto ss = tt[i];
				for (int j = 0; j < ss->Count; j++)
				{
					auto ee = ss[j];
					dat[current] = Eigen::Triplet<double>(ee->Item1, ee->Item2, 0);
					current++;
				}
			}
			this->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
		}
		myDoubleArray^ contract()
		{
			myDoubleArray^ vec = gcnew myDoubleArray(this->dat->_mat[0].cols());
			for (Int64 i = 0; i < this->dat->_mat[0].cols(); i++)
			{
				double val = this->dat->_mat[0].col(i).sum();
				vec->_arr->__v(i) = val;
			}
			return vec;
		}
		void ofStack(mySparse^ A, mySparse^ B)
		{
			std::vector<Eigen::Triplet<double>> dat;
			for (Int64 k = 0; k < A->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor,Int64>::InnerIterator it(A->dat->_mat[0], k); it; ++it)
				{
					dat.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				}
			}
			for (Int64 k = 0; k < B->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(B->dat->_mat[0], k); it; ++it)
				{
					dat.push_back(Eigen::Triplet<double>(it.row() + A->dat->_mat[0].rows(), it.col(), it.value()));
				}
			}
			for (Int64 k = 0; k < B->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(B->dat->_mat[0], k); it; ++it)
				{
					dat.push_back(Eigen::Triplet<double>(it.col(), it.row()+ A->dat->_mat[0].cols(), it.value()));
				}
			}

			this->dat->_mat[0].setZero();
			this->dat->_mat[0].reserve(dat.size());
			this->dat->_mat[0].resize(A->dat->_mat[0].rows() + B->dat->_mat[0].rows(), A->dat->_mat[0].cols() + B->dat->_mat[0].rows());
			this->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
		}
		mySparse^ trimMatrix(Int64 L1, Int64 L2)
		{
			/*public mySparse trimMatrix(mySparse D, Int64 L1, Int64 L2) {
				//return D.subMatrix(0, L1, 0, L1);

				var F = new mySparse(L1 * 2, L2 * 2);
				for (Int64 i = 0; i < L1; i++) {
					F._set(i, i, 1);
					F._set(i + L1, i + L2, 1);
				}
				return (F.Multiply(D) as mySparse).Multiply(F.Transpose()) as mySparse; 
			}*/
			mySparse^ ret = gcnew mySparse(L1*2, L1*2);
			std::vector<Eigen::Triplet<double>> dat;
			for (Int64 k = 0; k < this->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(this->dat->_mat[0], k); it; ++it)
				{
					if (it.col() < L1 && it.row() < L1)
					{
						dat.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
					}
					if (it.col()>L2&&it.col() < L1+L2 && it.row()>L2&&it.row() < L1+L2)
					{
						dat.push_back(Eigen::Triplet<double>(it.row()-(L2-L1), it.col() - (L2 - L1), it.value()));
					}
				}
			}
			ret->dat->_mat[0].setZero();
			ret->dat->_mat[0].reserve(dat.size());
			ret->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
			return ret;
		}
		mySparse^ zz()
		{
			mySparse^ ret = gcnew mySparse(this->dat->_mat[0].rows() / 3, this->dat->_mat[0].cols() / 3);
			std::vector<Eigen::Triplet<double>> dat;
			for (Int64 k = 0; k < this->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(this->dat->_mat[0], k); it; ++it)
				{
					if (it.row() % 3 == 2 && it.col() % 3 == 2)
					{
						dat.push_back(Eigen::Triplet<double>(it.row()/3, it.col()/3, it.value()));
					}
				}
			}
			ret->dat->_mat[0].setZero();
			ret->dat->_mat[0].reserve(dat.size());
			ret->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
			return ret;
		}
		mySparse^ xy()
		{
			mySparse^ ret = gcnew mySparse((this->dat->_mat[0].rows() / 3) * 2, (this->dat->_mat[0].cols() / 3) * 2);
			std::vector<Eigen::Triplet<double>> dat;
			for (Int64 k = 0; k < this->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(this->dat->_mat[0], k); it; ++it)
				{
					if (it.row() % 3 == 0 && it.col() % 3 == 0)
					{
						dat.push_back(Eigen::Triplet<double>(it.row() / 3, it.col() / 3, it.value()));
					}
					if (it.row() % 3 == 1 && it.col() % 3 == 1)
					{
						dat.push_back(Eigen::Triplet<double>(it.row() / 3+1, it.col() / 3+1, it.value()));
					}
				}
			}
			ret->dat->_mat[0].setZero();
			ret->dat->_mat[0].reserve(dat.size());
			ret->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
			return ret;
		}
		void ofStack2(mySparse^ A, mySparse^ B)
		{
			std::vector<Eigen::Triplet<double>> dat;
			for (Int64 k = 0; k < A->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(A->dat->_mat[0], k); it; ++it)
				{
					dat.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
				}
			}
			for (Int64 k = 0; k < B->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(B->dat->_mat[0], k); it; ++it)
				{
					dat.push_back(Eigen::Triplet<double>(it.row() + A->dat->_mat[0].rows(), it.col(), it.value()));
				}
			}

			this->dat->_mat[0].setZero();
			this->dat->_mat[0].reserve(dat.size());
			this->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
		}
		mySparse() {
			dat = 0;
			dat = new _mySparse();
			dat->init(0, 0);
		}
		mySparse(Int64 n, Int64 m)
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
		double sum() {
			if(this->dat->_mat[0].nonZeros()>0)
			return this->dat->_mat[0].cwiseAbs() .sum();
		}
		double _sum() {
			return this->dat->_dmat.cwiseAbs().sum();
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
		void ofAplusB(double alpha, mySparse^ A, double beta, mySparse^ B)
		{
			this->dat->init(B->rows(), B->cols());
			this->dat->_mat[0] = alpha * A->dat->_mat[0] + beta * B->dat->_mat[0];
		}
		mySparse^ subMatrix(Int64 i,Int64 N,Int64 j,Int64 M)
		{
			auto ret = gcnew mySparse(N, M);
			ret->dat->_mat[0] = this->dat->_mat[0].block(i, j, N, M);
			return ret;
		}
		void plus(mySparse^ B,double sc)
		{			
			this->dat->_mat[0] = this->dat->_mat[0]+B->dat->_mat[0];// *sc;
		}
		void plus(mySparseVector ^vec, double sc)
		{
			this->dat->_mat[0] += vec->_vec->_vec * vec->_vec->_vec.transpose() * sc;
		}
		void plus(mySparseVector^ vec, mySparseVector^ vec2,double sc)
		{
			this->dat->_mat[0] += vec->_vec->_vec * vec2->_vec->_vec.transpose() * sc;
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
			if(dat!=0)
				delete(dat);
			dat = 0;
		}
		!mySparse()
		{
			if(dat!=0)
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
		Int64 resize(Int64 n, Int64 m)
		{
			return dat->resize(n, m);
		}
		void reserve(Int64 n) {
			dat->reserve(n);
		}
		void addmat(mySparse^ m)
		{
			dat->addmat(m->dat);
		}
		Int64 cols() {
			return dat->cols();
		}
		Int64 rows() {
			return dat->rows();
		}
		Int64 _cols() {
			return dat->_cols();
		}
		Int64 _rows() {
			return dat->_rows();
		}
		Int64 __rows() {
			return dat->__rows();
		}
		void shrink(Int64 M)
		{
			dat->shrink(M);
		}
		void permute(myPermutation^ p)
		{
			dat->permute(p->p->perm);
		}
		void _shrink(Int64 M, bool sparse, bool dense)
		{
			dat->_shrink(M, sparse, dense);
		}
		void _permute(myPermutation^ p, bool sparse, bool dense)
		{
			dat->_permute(p->p->perm, sparse, dense);
		}
		void _permuteCols(myPermutation^ p, bool sparse, bool dense)
		{
			dat->_permuteCols(p->p->perm, sparse, dense);
		}
		void _shrink(Int64 M, Int64 N)
		{
			dat->_shrink(M, N);
		}
		void _permute(myPermutation^ p, myPermutation^ q)
		{
			dat->_permute(p->p->perm, q->p->perm);
		}
		System::String ^ ofAtA(mySparse^ m, bool sparse) {
			return gcnew System::String(dat->ofAtA(m->dat, sparse).c_str());
		}
		System::String^ ofAtA_gpu(myCuda ^cuda,mySparse^ m, bool sparse) {
			return gcnew System::String( dat->ofAtA_gpu(cuda->cuda(), m->dat, sparse).c_str());
		}
		System::String^ _ofAtA(mySparse^ m) {
			auto str = dat->_ofAtA(m->dat);
			auto ret = gcnew System::String(str.c_str());
			return ret;
		}
		void ofAtB(mySparse^ m, bool sparse) {
			dat->ofAtB(m->dat, sparse);
		}
		/*void ofAtB_gpu(mySparse^ m, bool sparse) {
			dat->ofAtB_gpu(m->dat, sparse);
		}*/
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
		
		
		void addemptyrow(Int64 ii) {
			dat->addemptyrow(ii);
		}
		System::String^ info() {
			auto str = this->dat->info();
			System::String^ s = gcnew System::String(str.data(), 0, str.size());
			return s;
		}
		void addrow(Int64 ii, array<Int64>^ index, array<double>^ data, double sc, Int64 N) {
			pin_ptr<Int64> ptr = &index[0];
			pin_ptr<double> dat = &data[0];

			Int64* _ptr = ptr;
			double* _dat = dat;
			this->dat->addrow(ii, _ptr, _dat, sc, N);
			ptr = nullptr;
			dat = nullptr;
		}
		void addrow(Int64 ii, array<Int64>^ index, array<double>^ data, Int64 shift, double sc, Int64 N, bool add) {
			pin_ptr<Int64> ptr = &index[0];
			pin_ptr<double> dat = &data[0];
			Int64* _ptr = ptr;
			double* _dat = dat;

			this->dat->addrow(ii, _ptr, _dat, shift, sc, N, add,1.0);
			ptr = nullptr;
			dat = nullptr;
		}
		void addrow(Int64 ii, myIntArray^ index, myDoubleArray^ data, double sc, Int64 N) {
			//pin_ptr<Int64> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), sc, N);
			//ptr = nullptr;
			//dat = nullptr;
		}
		void addrow(Int64 ii, myIntArray^ index, myDoubleArray^ data, double sc, Int64 N,double coeff) {
			//pin_ptr<Int64> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), sc, N, coeff);
			//ptr = nullptr;
			//dat = nullptr;
		}
		void addrow(Int64 ii, myIntArray^ index, myDoubleArray^ data, Int64 shift, double sc, Int64 N, bool add) {
			//pin_ptr<Int64> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), shift, sc, N, add,1.0);
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
		void Atb(array<double>^ b, array<double>^ c,double sc, KingOfMonsters::myDoubleArray^ ret)
		{
			pin_ptr<double> ptr = &b[0];
			pin_ptr<double> ptr2 = &c[0];
			dat->Atb(ptr, ptr2,sc,b->Length, &ret->_arr->__v);
			//array<double>^ ret = gcnew array<double>(rhs.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)rhs.data(), ret, 0, rhs.rows());
			ptr = nullptr;
			ptr2 = nullptr;
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
		void _solve0_lu(myDoubleArray^ rhs, myDoubleArray^ ret,int ordering,bool meh) {
			//pin_ptr<double> ptr = &rhs[0];
			//dat->_mat[0].setIdentity();
			if (meh) {
				rhs->_arr->__v = this->dat->_mat[0].transpose() * rhs->_arr->__v;
				this->dat->_mat[0] = this->dat->_mat[0].transpose() * this->dat->_mat[0];
			}
			double nn = 0.00000000001;
			bool allocerr = false;
			std::string _str = "";
			for (int i = 0; i < 50; i++)
			{
				std::string str = dat->_solve0_lu(&rhs->_arr->__v, &ret->_arr->__v, ordering);
				_str += str;
				if (str.find("PIVOT") != string::npos) {

					nn *= 100;
					this->dat->addsmallidentity(nn, true, false);
				}
				else {
					if (str.find("ALLOCERR") != string::npos) {
						allocerr = true;
					}
					break;

				}
			}
			if (allocerr)
			{
				_str += std::string("\n") + std::string("swithing to CPU") + std::string("\n");
				for (int i = 0; i < 50; i++)
				{
					std::string str = dat->_solve0_lu_cpu(&rhs->_arr->__v, &ret->_arr->__v, ordering);
					_str += str;
					if (str.find("PIVOT") != string::npos) {

						nn *= 100;
						this->dat->addsmallidentity(nn, true, false);
					}
					else {
						break;
					}
				}
			}
			if (_str == "")_str = "success";
			System::Console::WriteLine(gcnew System::String(_str.c_str()));
		}
		void _solve0_lu_cpu(myDoubleArray^ rhs, myDoubleArray^ ret, int ordering, bool meh) {
			if (meh) {
				rhs->_arr->__v = this->dat->_mat[0].transpose() * rhs->_arr->__v;
				this->dat->_mat[0] = this->dat->_mat[0].transpose() * this->dat->_mat[0];
			}
			double nn = 0.00000000001;
			bool allocerr = false;
			std::string _str = "";
			
			for (int i = 0; i < 10; i++)
			{
				std::string str = dat->_solve0_lu_cpu(&rhs->_arr->__v, &ret->_arr->__v, ordering);
				
				//lu.compute(this->dat->_mat[0]);
				_str += str;
				if(str.find("SUCCESS") == string::npos) 
				{
					nn *= 100;
					this->dat->addsmallidentity(nn, true, false);
				}
				else {					
					break;
				}
			}
			if (_str == "")_str = "success";
			System::Console::WriteLine(gcnew System::String(_str.c_str()));
		}
		void solve0_lu(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->solve0_lu(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}
		void _solve0_lu_cg(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->_solve0_lu_cg(&rhs->_arr->__v, &ret->_arr->__v);

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

		System::String^ _solve0_gpu(myCuda^ gpu, myDoubleArray^ rhs, myDoubleArray^ ret, Int64 device) {
			//pin_ptr<double> ptr = &rhs[0];

			auto ss = dat->_solve0_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);
			System::String^ ee = gcnew System::String(ss.c_str());
			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
			return ee;
		}
		System::String^ _solveLU_gpu(myCuda^ gpu, myDoubleArray^ rhs, myDoubleArray^ ret, Int64 device) {
			//pin_ptr<double> ptr = &rhs[0];

			auto ss = dat->_solveLU_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);
			System::String^ ee = gcnew System::String(ss.c_str());
			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
			return ee;
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

		void solve0_gpu(myCuda^ gpu, mySparse^ rhs, mySparse^ ret) {
			this->dat->_solve0_gpu(gpu->cuda(), rhs->dat, ret->dat);
		}
		Int64 solveI(mySparse^ ret) {
			return this->dat->_solveI(ret->dat);
		}
		System::String^ solveI_gpu(myCuda^ gpu, mySparse^ ret) {
			auto _ss=this->dat->_solveI_gpu(gpu->cuda(), ret->dat);
			auto ss = gcnew System::String(_ss.c_str());

			return ss;
		}
		System::String^ solveI_gpu_sparse(myCuda^ gpu, mySparse^ ret) {
			auto _ss = this->dat->_solveI_gpu_sparse(gpu->cuda(), ret->dat);
			auto ss = gcnew System::String(_ss.c_str());

			return ss;
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
		/*mySparse^ solveI_gpu_mg(myCuda^ gpu, mySparse^ ret) {
			this->dat->_solveI_gpu_mg(gpu->cuda(), ret->dat);
			return ret;
		}*/
		void freezecoeff()
		{
			dat->freezecoeff();
		}
		void addcoeff(double sc) {
			dat->addcoeff(sc);
		}
		double at_each(Int64 i, Int64 ii) {
			return dat->at(i, ii);
		}
		void begin_construct()
		{
			dat->begin_construct();
		}
		void end_construct(Int64 c)
		{
			dat->end_construct(c);
		}
		Int64 num_elem(Int64 ii)
		{
			return dat->num_elem(ii);
		}
		double _at(Int64 i) {
			return dat->_at(i);
		}
		double _at(Int64 i, Int64 j) {
			return dat->_at(i, j);
		}
		double at(Int64 i, Int64 j) {
			return dat->__at(i, j);
		}
		/*void _plus(Int64 i, Int64 j, double val) {
			this->dat->_mat[0].coeffRef(i,j) += val;
		}*/
		void _plus(Int64 i, Int64 j, double val) {
			this->dat->_plus(i, j, val);
			//this->dat->_mat[0].coeffRef(i, j) += val;
		}
		void _set(Int64 i, Int64 j, double val) {
			this->dat->_mat[0].coeffRef(i, j) = val;
		}
		void increaseCapacityBy(Int64 nn)
		{
			if (this->dat->_mat[0].data().allocatedSize() - this->dat->_mat[0].nonZeros() < nn)
			{
				this->dat->_mat[0].reserve(this->dat->_mat[0].nonZeros() + nn * 2);
			}
		}
		Int64 nonZeros() {
			return dat->nonzeros();
		}
		void adddat(Int64 i, Int64 j, double val) {
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
		myDoubleArray^ vector(myDoubleArray^ a) {
			auto _ret = dat->Vector(&a->_arr->__v);
			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			myDoubleArray^ ret = gcnew myDoubleArray(_ret.rows());
			ret->_arr->__v = _ret;
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
		Int64 numBlocks()
		{
			return this->dat->numBlocks();
		}
		static System::String^ _testopenmp()
		{
			auto _str = KingOfMonsters::_mySparse::_testopenmp();
			return gcnew System::String(_str.c_str());
		}
		static System::String^ testopenmp()
		{
			auto str = gcnew System::String("");
			Int64 mt = omp_get_max_threads();
			str += "num threads:" + mt.ToString() + "\n";
#pragma omp parallel for
			for (Int64 i = 0; i < 100; i++)
			{
				Int64 ct = omp_get_thread_num();
#pragma omp critical
				{
					str += "current threads:" + ct.ToString() + "\n";
				}
				for (Int64 t = 0; t < ct * 100; t++)
				{
					auto f = new double[10];
				}
			}
			return str;
		}
	};

	public ref class helper {
	public:
		static double computeeigen(mySparse^ i1, mySparse^ i2, mySparse^ t1, Int64 N)
		{
			double maxval = _mySparse::computeeigen(i1->dat, i2->dat, t1->dat, N);


			return maxval;
			//auto vv = mat.eigenvalues();
		}
		static double computeeigen2(mySparse^ i1, mySparse^ i2, mySparse^ t1, Int64 N)
		{
			double maxval = _mySparse::computeeigen2(i1->dat, i2->dat, t1->dat, N);


			return maxval;
			//auto vv = mat.eigenvalues();
		}
		static void minilla(mySparse^ i1, mySparse^ i2, mySparse^ i3, myDoubleArray^ grad,myDoubleArray^ ret1, myDoubleArray^ ret2)
		{
			Eigen::VectorXd v = _mySparse::minilla(i1->dat, i2->dat, i3->dat, grad->_arr);

			ret1->_arr->__v = v.topRows(i1->dat->cols());
			ret2->_arr->__v = v.bottomRows(i2->dat->cols());

			//auto vv = mat.eigenvalues();
		}

	};
}