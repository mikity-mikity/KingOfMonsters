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
	
	public ref class denseMatrix {
	private:
		Eigen::MatrixXd* mat = 0;
	public:
		void genEigen(denseMatrix ^a,denseMatrix ^ b, [Runtime::InteropServices::Out]double%  l1, [Runtime::InteropServices::Out]double% l1i,[Runtime::InteropServices::Out]double% l2, [Runtime::InteropServices::Out]double% l2i){
			Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solve(* (a->mat), * (b->mat), true);
			
			/*if (std::abs(a->mat->determinant()) > std::abs(b->mat->determinant()))
			{
				Eigen::MatrixXd M=a->mat->inverse() * *b->mat;
				Eigen::EigenSolver<Eigen::MatrixXd> solve(M);
				*this->mat = solve.eigenvectors().real();
				l1 = solve.eigenvalues()(0).real();
				l2 = solve.eigenvalues()(1).real();
				l1i = solve.eigenvalues()(0).imag();
				l2i = solve.eigenvalues()(1).imag();
			}
			else {
				Eigen::MatrixXd M = b->mat->inverse() * *a->mat;
				Eigen::EigenSolver<Eigen::MatrixXd> solve(M);
				*this->mat = solve.eigenvectors().real();
				l1 = solve.eigenvalues()(0).real();
				l2 = solve.eigenvalues()(1).real();
				l1i = solve.eigenvalues()(0).imag();
				l2i = solve.eigenvalues()(1).imag();
			}*/

			* this->mat = solve.eigenvectors().real();
			l1 = solve.eigenvalues()(0).real();
			l2 = solve.eigenvalues()(1).real();
			l1i = solve.eigenvalues()(0).imag();
			l2i = solve.eigenvalues()(1).imag();

		}
		Eigen::MatrixXd& get()
		{
			return *mat;
		}
		void set(Eigen::MatrixXd _m)
		{
			*mat = _m;
		}
		void set(int i, int j, double val)
		{
			(*mat)(i, j) = val;
		}
		double get(int i, int j)
		{
			return (*mat)(i, j);
		}
		void setzero() {
			mat->setZero();
		}
		void resize(int n, int m)
		{
			mat->resize(n, m);
		}
		int rows() {
			return mat->rows();
		}
		int cols() {
			return mat->cols();
		}
		denseMatrix(int n, int m)
		{
			if (mat != 0)del();
			mat = new Eigen::MatrixXd();
			mat->resize(n, m);
			mat->setZero();
		}
		denseMatrix()
		{
			if (mat != 0)del();
			mat = new Eigen::MatrixXd();
		}
		!denseMatrix()
		{
			del();
		}
		~denseMatrix()
		{
			del();
		}
		void del() {
			if (mat != 0)
			{
				delete mat;
			}
			mat = 0;
		}
	};
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
		void scale(double sc)
		{
			this->_vec->_vec *= sc;
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
		_myDoubleArray* _arr = 0;
		Int64 _N = 0;
	public:
		void assemble(myDoubleArray^ v, myDoubleArray^ w)
		{
			this->_arr->__v.resize(v->_arr->__v.rows() + w->_arr->__v.rows());
			this->_arr->__v.topRows(v->_arr->__v.rows()) = v->_arr->__v;
			this->_arr->__v.bottomRows(w->_arr->__v.rows()) = w->_arr->__v;

		}
		void split(myDoubleArray^ v, myDoubleArray^ w, int N)
		{
			v->_arr->__v = this->_arr->__v.topRows(N);
			w->_arr->__v = this->_arr->__v.bottomRows(this->_arr->__v.rows()-N);

		}
		void transpose()
		{
			this->_arr->__v.transposeInPlace();
		}
		
		myDoubleArray^ duplicate()
		{
			myDoubleArray^ _new = gcnew myDoubleArray(this->_arr->__v.rows());
			_new->_arr->__v = this->_arr->__v;
			return _new;
		}
		void _plus(mySparseVector^ v, double sc)
		{
			this->_arr->__v += sc * v->_vec->_vec;
		}
		void plus_useindex(myDoubleArray^ vec, double sc, int N, array<int>^index)
		{
			for (int i = 0; i < N; i++)
			{
				this->_arr->__v(index[i]) += sc * vec->_arr->__v(i);
			}
		}
		void addResidual(myDoubleArray^ r, myDoubleArray^ ret)
		{
			ret->_arr->__v += this->_arr->__v * r->_arr->__v;
		}
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
		void setzero()
		{
			_arr->__v.setZero();
		}
		void setconstant(double val)
		{
			_arr->__v.setConstant(val);
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
 			//Eigen::VectorXd v(N);
			//v.setZero();
			if (_N < N)
			{
				//v.middleRows(0, _N) = this->_arr->__v;
				this->_arr->__v.conservativeResize(N);
				this->_arr->__v.bottomRows(N - _N).setZero();
				//this->_arr->__v = v;
			}
			else if(N<_N){
				this->_arr->__v.conservativeResize(N);
				//v = this->_arr->__v.middleRows(0, N);
				//this->_arr->__v.resize(N);
				//this->_arr->__v = v;
			}
			//_arr->__v.conservativeResize(N);
			/*if (N > _N)
			{
				_arr->__v.middleRows(_N, N - _N).setZero();
			}*/
			_N = N;
		}
		void split(int N,myDoubleArray ^ ret)
		{
			ret->_arr->__v = this->_arr->__v.bottomRows(this->_arr->__v.size() - N);
			this->_arr->__v = this->_arr->__v.topRows(N);
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
		void minus(myDoubleArray^ b,double sc)
		{
			_arr->__v.noalias() -= b->_arr->__v*sc;
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
	public ref class myIntArray2 {
	public:
		Int64* _arr = 0;
		Int64 _N = 0;  //rows
		Int64 _M = 0;  //cols
	public:		
		myIntArray2(Int64 N,Int64 M)
		{
			_arr = new Int64[N*M];
			_N = N;
		}
		inline Int64* data()
		{
			return _arr;
		}
		inline Int64 size()
		{
			return _N*_M;
		}
		!myIntArray2()
		{
			if (_arr != 0)
			{
				delete[] _arr;
			}
			_arr = 0;
			_N = 0;
		}
		~myIntArray2()
		{
			if (_arr != 0) {
				delete[] _arr;
			}
			_arr = 0;
			_N = 0;
		}
		void set(Int64 i, int j,Int64 val)
		{
			_arr[i*_M+j] = val;
		}
		Int64 at(Int64 i,Int64 j)
		{
			return _arr[i*_M+j];
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
		myPermutation(int N,System::Collections::Generic::List<long long>^ index)
		{
			int L1 = N-1;
			array<long long>^ ii = gcnew array<long long>(N);
			for (int i = 0; i < index->Count; i++)
			{
				ii[index[i]] = L1;
				L1--;
			}
			int L2 = 0;
			for (int i = 0; i < N; i++)
			{
				if (index->Contains(i)) {}
				else {
					ii[i] = L2;
					L2++;
				}
				
			}
			pin_ptr<Int64> ptr = &ii[0];
			p = new _myPermutation(ptr, N);
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
			//ptr2 = nullptr;es
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
		int size()
		{
			return this->p->perm.size();
		}
	};
	
	public ref class mySparse {
	public:
		_mySparse* dat = 0;
		System::String^ tostring()
		{
			System::String^ str = gcnew System::String("");
			auto ee = Eigen::SparseMatrix<double,0,int64_t>(dat->_mat[0].transpose());
			double ss = (dat->_mat[0] - ee).squaredNorm();
			//double ss = ff.sum();
			int count = 0;
			for (int i = 0; i < dat->_mat[0].cols(); i++)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(dat->_mat[0], i); it; ++it) 
				{
					str = str + it.row().ToString() + "," + it.col().ToString() + "," + it.value().ToString() +"::" +dat->_mat[0].coeffRef(it.col(), it.row()).ToString() + "\n";
					count++;
					if (count > 40)break;
				}
				if (count > 40)break;
			}
			return str;
		}
		mySparse^ AtA() {
			mySparse^ newMat = gcnew mySparse();
			newMat->dat->_mat.resize(1);
			newMat->dat->_mat[0] = this->dat->_mat[0].transpose() * this->dat->_mat[0];
			return newMat;
		}
		mySparse^ AtBC(mySparse ^B,mySparse ^C)
		{
			mySparse^ newMat = gcnew mySparse();
			newMat->dat->_mat.resize(1);
			newMat->dat->_mat[0] = this->dat->_mat[0].transpose() * B->dat->_mat[0]*C->dat->_mat[0];
			return newMat;
		}
		void addResidual(myDoubleArray^ r, myDoubleArray^ ret)
		{
			ret->_arr->__v += this->dat->_mat[0].transpose() * r->_arr->__v;
		}
		void scale(double sc)
		{
			this->dat->scale(sc);
		}
		void _scale(double sc)
		{
			this->dat->_scale(sc);
		}
		void plus(mySparse^ m, bool sparse)
		{
			if (sparse)
			{
				this->dat->_mat[0] += m->dat->_mat[0];
			}
			else {
				this->dat->_dmat += m->dat->_mat[0];
			}
		}
		void plus(mySparse^ m, double sc)
		{
				this->dat->_mat[0] += m->dat->_mat[0]*sc;
		}
		void plusvvt(myDoubleArray^ v,double sc)
		{
			this->dat->_dmat += sc*v->_arr->__v * v->_arr->__v.transpose();
		}
		void _plusvvt(mySparseVector ^ v, double sc)
		{
			this->dat->_mat[0] += sc * v->_vec->_vec * v->_vec->_vec.transpose();
		}
		
		void multiply(myDoubleArray^ v, myDoubleArray^ ret)
		{
			ret->_arr->__v = this->dat->_mat[0] * v->_arr->__v;
			if (this->dat->coeff[0].size() == this->dat->_mat[0].rows())
			{
				ret->_arr->__v = this->dat->coeff[0].asDiagonal() * ret->_arr->__v;
			}
		}
		void _multiply(myDoubleArray^ v, myDoubleArray^ ret)
		{
			ret->_arr->__v = this->dat->_mat[0] * v->_arr->__v;
		}
		void leftmultiply(myDoubleArray^ v, myDoubleArray^ ret)
		{
			ret->_arr->__v = v->_arr->__v.transpose() * this->dat->_mat[0];
		}
		void ofAtAsimple(mySparse^ m,bool sparse)
		{
			if (sparse)
			{
				this->dat->_mat.resize(1);
				this->dat->_mat[0] = m->dat->_mat[0].transpose() * m->dat->_mat[0];
			}
			else {
				this->dat->_dmat = m->dat->_mat[0].transpose() * m->dat->_mat[0];
			}
		}

		mySparse^ computeKernel()
		{
			Eigen::SparseMatrix<double> ff = (this->dat->_mat[0] * this->dat->_mat[0].transpose());
			Eigen::SparseLU<Eigen::SparseMatrix<double>> lu;
			lu.compute(ff);
			auto I = Eigen::SparseMatrix<double>(ff.rows(), ff.cols());
			Eigen::MatrixXd ret = lu.solve(I);
			auto II = Eigen::SparseMatrix<double>(this->dat->_mat[0].cols(), this->dat->_mat[0].cols());

			mySparse^ kernel = gcnew mySparse();
			kernel->dat->_dmat = II - this->dat->_mat[0].transpose() * ret * this->dat->_mat[0];
			return kernel;
		}
		void resetData()
		{
			this->dat->dat.resize(1);
			this->dat->dat[0].clear();
			this->dat->_nt = 1;
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
		void setZero() {
			this->dat->_mat[0].setZero();
		}
		void setZero(int row)
		{
			this->dat->setzero(row);
		}
		void setOne(int N)
		{
			std::vector<Eigen::Triplet<double>> dat;
			dat.resize(N * N);
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					dat.push_back(Eigen::Triplet<double>(i, j, 1));
				}
			}
			this->dat->_mat[0].resize(N, N);
			this->dat->_mat[0].setFromTriplets(dat.begin(), dat.end());
		}
		
		
		void ofStack(System::Collections::Generic::List<mySparse^> ^jacobians)
		{
			int M = this->dat->_mat[0].rows();
			int N = this->dat->_mat[0].cols();
			for (int i = 0; i < jacobians->Count; i++)
			{
				M += jacobians[i]->dat->_mat[0].rows();
			}
			this->dat->_dmat.resize(M, N);

			this->dat->_dmat.topRows(this->dat->_mat[0].rows()) = this->dat->_mat[0];
			M = this->dat->_mat[0].rows();
			for (int i = 0; i < jacobians->Count; i++)
			{
				this->dat->_dmat.middleRows(M, jacobians[i]->dat->_mat[0].rows()) = jacobians[i]->dat->_mat[0];
				M += jacobians[i]->dat->_mat[0].rows();
			}
			return;
		}
		void assemble(mySparse^ A, mySparse^ B, mySparse^ C)
		{
			this->dat->_dmat.resize(A->dat->_dmat.cols() + C->dat->_dmat.cols(), A->dat->_dmat.cols() + C->dat->_dmat.cols());
			this->dat->_dmat.setZero();
			this->dat->_dmat.topLeftCorner(A->dat->_dmat.cols(), A->dat->_dmat.cols()) = A->dat->_dmat;
			this->dat->_dmat.bottomRightCorner(C->dat->_dmat.cols(), C->dat->_dmat.cols()) = C->dat->_dmat;
			if (B != nullptr)
			{
				this->dat->_dmat.topRightCorner(A->dat->_dmat.rows(), C->dat->_dmat.cols()) = B->dat->_dmat;
				this->dat->_dmat.bottomLeftCorner(C->dat->_dmat.cols(), A->dat->_dmat.cols()) = B->dat->_dmat.transpose();
			}

		}
		void ofStack(mySparse^ A, mySparse^ B)
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
			for (Int64 k = 0; k < B->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(B->dat->_mat[0], k); it; ++it)
				{
					dat.push_back(Eigen::Triplet<double>(it.col(), it.row() + A->dat->_mat[0].cols(), it.value()));
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
			mySparse^ ret = gcnew mySparse(L1 * 2, L1 * 2);
			std::vector<Eigen::Triplet<double>> dat;
			for (Int64 k = 0; k < this->dat->_mat[0].outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, Int64>::InnerIterator it(this->dat->_mat[0], k); it; ++it)
				{
					if (it.col() < L1 && it.row() < L1)
					{
						dat.push_back(Eigen::Triplet<double>(it.row(), it.col(), it.value()));
					}
					if (it.col() > L2 && it.col() < L1 + L2 && it.row() > L2 && it.row() < L1 + L2)
					{
						dat.push_back(Eigen::Triplet<double>(it.row() - (L2 - L1), it.col() - (L2 - L1), it.value()));
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
						dat.push_back(Eigen::Triplet<double>(it.row() / 3, it.col() / 3, it.value()));
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
						dat.push_back(Eigen::Triplet<double>(it.row() / 3 + 1, it.col() / 3 + 1, it.value()));
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
		void ofvv(myDoubleArray^ v,double sc)
		{
			int n = v->_arr->__v.size();
			this->dat->_dmat.resize(n, n);
			
			this->dat->_dmat = sc*(v->_arr->__v * v->_arr->__v.transpose());
		}
		double sum() {
			if (this->dat->_mat[0].nonZeros() > 0)
				return this->dat->_mat[0].cwiseAbs().sum();
			return -1;
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
		mySparse^ subMatrix(Int64 i, Int64 N, Int64 j, Int64 M)
		{
			auto ret = gcnew mySparse(N, M);
			ret->dat->_mat[0] = this->dat->_mat[0].block(i, j, N, M);
			return ret;
		}
		void plus(mySparseVector^ vec, double sc)
		{
			this->dat->_mat[0] += vec->_vec->_vec * vec->_vec->_vec.transpose() * sc;
		}
		void plus(mySparseVector^ vec, mySparseVector^ vec2, double sc)
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
		mySparse^ _duplicate()
		{
			mySparse^ newmat = gcnew mySparse();
			newmat->dat->_mat.push_back(this->dat->_mat[0]);
			return newmat;
		}
		~mySparse()
		{
			if (dat != 0)
				delete(dat);
			dat = 0;
		}
		!mySparse()
		{
			if (dat != 0)
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
		void setZero(bool sparse) {
			if (sparse)
			{
				this->dat->_mat[0].setZero();
			}
			else {
				this->dat->_dmat.setZero();
			}
		}
		void reserve(Int64 n) {
			dat->reserve(n);
		}
		void addmat(mySparse^ m)
		{
			dat->addmat(m->dat);
		}
		void scale(int i, double sc)
		{
			this->dat->scale(i, sc);
		}
		Int64 cols() {
			return dat->cols();
		}
		void _resize(int n, int m)
		{
			Eigen::SparseMatrix<double> mm(n, m);
			for (int64_t k = 0; k < this->dat->_mat[0].cols(); ++k) {
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(this->dat->_mat[0], k); it; ++it) {
					mm.coeffRef(it.row(), it.col()) = it.value();
				}
			}
			this->dat->_mat[0].resize(n, m);
			this->dat->_mat[0].setZero();
			this->dat->_mat[0] = mm;
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
		void permute_dense(myPermutation^ p)
		{
			dat->_dmat = p->p->perm * dat->_dmat * p->p->perm.transpose();
		}
		void _shrink_dense(Int64 N)
		{
			this->dat->_dmat.conservativeResize(N, N);
		}
		void _shrink(Int64 M, bool sparse, bool dense)
		{
			dat->_shrink(M, sparse, dense);
		}
		void _shrinkCols(Int64 M, bool sparse, bool dense)
		{
			dat->_shrinkCols(M, sparse, dense);
		}
		void memory()
		{
			this->dat->_prevmat = this->dat->_dmat;
		}
		double approximate(mySparse^ A,double sc)
		{
			Eigen::MatrixXd deltaA = A->dat->_dmat - A->dat->_prevmat;
			deltaA *= sc;
			auto I = Eigen::MatrixXd::Identity(A->dat->_dmat.rows(), A->dat->_dmat.cols());
			Eigen::MatrixXd C = I - A->dat->_prevmat * this->dat->_dmat;
			Eigen::MatrixXd deltaB=this->dat->_dmat*(C - deltaA * this->dat->_dmat);
			this->dat->_dmat += deltaB;
			return (Eigen::MatrixXd::Identity(A->dat->_dmat.rows(), A->dat->_dmat.cols()) - A->dat->_dmat * this->dat->_dmat).sum();

		}
		void _permute(myPermutation^ p, bool sparse, bool dense)
		{
			dat->_permute(p->p->perm, sparse, dense);
		}
		void _permuteCols(myPermutation^ p, bool sparse, bool dense)
		{
			dat->_permuteCols(p->p->perm, sparse, dense);
		}
		void _shrink(Int64 M, Int64 N,bool sparse,bool dense)
		{
			dat->_shrink(M, N,sparse,dense);
		}
		void _permute(myPermutation^ p, myPermutation^ q)
		{
			

				dat->_permute(p->p->perm, q->p->perm);
			
			
		}
		void _permuteL(myPermutation^ p, myPermutation^ q)
		{
			
				
				{
					dat->__permuteRows(p->p->perm);
				}
				

		}
		void _permuteR(myPermutation^ p, myPermutation^ q)
		{

			{
				dat->__permuteCols(q->p->perm);
			}

		}
		System::String^ ofAtA(mySparse^ m, bool sparse) {
			return gcnew System::String(dat->ofAtA(m->dat, sparse).c_str());
		}
		System::String^ ofAtA_gpu(myCuda^ cuda, mySparse^ m, bool sparse) {
			return gcnew System::String(dat->ofAtA_gpu(cuda->cuda(), m->dat, sparse).c_str());
		}
		System::String^ _ofAtA(mySparse^ m) {
			auto str = dat->_ofAtA(m->dat);
			auto ret = gcnew System::String(str.c_str());
			return ret;
		}
		System::String^ _ofAtA_sparse(mySparse^ m) {
			auto str = dat->_ofAtA_sparse(m->dat);
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

			A->dat->_ofBtAB(B->dat, this->dat);
			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			//ptr = nullptr;
			//return ret;
		}
		void _ofBtAB(mySparse^ A, mySparse^ B, mySparse^ B2, myDoubleArray^ b, myDoubleArray^ ret)
		{
			//pin_ptr<double> ptr = &b[0];

			A->dat->_ofBtAB(B->dat, B2->dat,this->dat);
			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			//ptr = nullptr;
			//return ret;
		}

		double _trace()
		{
			return this->dat->_dmat.trace();
		}
		void _ofCBtAB(mySparse^ C, mySparse^ A, mySparse^ B)
		{
			//pin_ptr<double> ptr = &b[0];
			//singularvalues = gcnew myDoubleArray(0);
			A->dat->_ofCBtAB(B->dat, C->dat, this->dat);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			//ptr = nullptr;
			//return ret;
		}
		void _ofCBtAB2(mySparse^ C, mySparse^ A, mySparse^ B, mySparse^ D)
		{
			//pin_ptr<double> ptr = &b[0];
			//singularvalues = gcnew myDoubleArray(0);
			A->dat->_ofCBtAB2(B->dat, C->dat, this->dat, D->dat);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());
			//ptr = nullptr;
			//return ret;
		}
		void mult(myDoubleArray^ a, myDoubleArray^ b)
		{
			b->_arr->__v = this->dat->_dmat * a->_arr->__v;
		}
		void enabledense()
		{
			this->dat->_dmat.resize(this->dat->_mat[0].rows(), this->dat->_mat[0].cols());
			this->dat->_dmat.setZero();
		}
		void makePattern()
		{
			this->dat->makePattern();
		}
		void plusmat_usemap(denseMatrix^ a, double sc, int N, array<int>^ index)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					this->dat->add_usemap(index[i], index[j], sc * (a->get())(i,j) );
				}
			}

		}
		void plusdyad_usemap(myDoubleArray^ a, myDoubleArray^ b, double sc,int N,array<int> ^index)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					this->dat->add_usemap(index[i], index[j], sc * a->_arr->__v(i) * b->_arr->__v(j));
				}
			}
		}

		void plusmat_uselocation(denseMatrix^ a, double sc, int N, array<int,2>^ location)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					this->dat->add_uselocation(location[i,j], sc * (a->get())(i, j));
				}
			}

		}
		void plusdyad_uselocation(myDoubleArray^ a, myDoubleArray^ b, double sc, int N, array<int,2>^ location)
		{
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					this->dat->add_uselocation(location[i,j], sc * a->_arr->__v(i) * b->_arr->__v(j));
				}
			}
		}
		void plus_usemap(Int64 i, Int64 j, double value)
		{
			this->dat->add_usemap(i, j, value);
		}
		void set_usemap(Int64 i, Int64 j, double value)
		{
			this->dat->set_usemap(i, j, value);
		}
		int find_location(Int64 i, Int64 j)
		{
			return this->dat->find_location(i, j);
		}
		void plus_uselocation(Int64 location, double value)
		{
			this->dat->add_uselocation(location, value);
		}
		void set_useloaction(Int64 location, double value)
		{
			this->dat->set_uselocation(location, value);
		}

		void mult(mySparse^ a, mySparse^ b)
		{
			b->dat->_dmat = this->dat->_dmat * a->dat->_dmat;
		}
		void _mult(myDoubleArray^ a, myDoubleArray^ b)
		{
			b->_arr->__v = this->dat->_mat[0] * a->_arr->__v;
		}
		void _transposemultplus(myDoubleArray^ a, myDoubleArray^ b)
		{
			b->_arr->__v += this->dat->_mat[0].transpose() * a->_arr->__v;
		}
		void _mult(mySparse^ a, mySparse^ b)
		{
			b->dat->_dmat = this->dat->_mat[0] * a->dat->_dmat;
		}
		void _mult_dense(myDoubleArray^ a, myDoubleArray^ ret)
		{
			ret->_arr->__v = a->_arr->__v.transpose() * this->dat->_dmat;
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

			this->dat->addrow(ii, _ptr, _dat, shift, sc, N, add, 1.0);
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
		void addrow(Int64 ii, myIntArray^ index, myDoubleArray^ data, double sc, Int64 N, double coeff) {
			//pin_ptr<Int64> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), sc, N, coeff);
			//ptr = nullptr;
			//dat = nullptr;
		}
		void addrow(Int64 ii, myIntArray^ index, myDoubleArray^ data, Int64 shift, double sc, Int64 N, bool add) {
			//pin_ptr<Int64> ptr = &index[0];
			//pin_ptr<double> dat = &data[0];

			this->dat->addrow(ii, index->_arr, data->_arr->__v.data(), shift, sc, N, add, 1.0);
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
		void Atb(array<double>^ b,  KingOfMonsters::myDoubleArray^  ret)
		{
			pin_ptr<double> ptr = &b[0];
			dat->Atb(ptr, b->Length,&ret->_arr->__v);
			//array<double>^ ret = gcnew array<double>(rhs.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)rhs.data(), ret, 0, rhs.rows());
			ptr = nullptr;
			//return ret;
		}
		void Atb(array<double>^ b, array<double>^ c, double sc, KingOfMonsters::myDoubleArray^ ret)
		{
			pin_ptr<double> ptr = &b[0];
			pin_ptr<double> ptr2 = &c[0];
			dat->Atb(ptr, ptr2, sc, b->Length, &ret->_arr->__v);
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
		void fillZeros()
		{
			memset((this->dat->_mat[0]).valuePtr(), 0, sizeof(double) * this->dat->_mat[0].nonZeros());
		}
		void solve0(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->solve0(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}
		void LSsolve(myDoubleArray^ rhs, myDoubleArray^ ret,double salt,int mode) {

			dat->LSsolve(&rhs->_arr->__v, &ret->_arr->__v,salt,mode);

		}
		void Project(myDoubleArray^ rhs, myDoubleArray^ ret, double salt) {

			dat->Project(&rhs->_arr->__v, &ret->_arr->__v,salt);

		}
		void trimCols(System::Collections::Generic::List<int> ^ index)
		{
			int N = this->dat->_mat[0].cols();
			int L1 = N;
			array<long long>^ ii = gcnew array<long long>(N);
			for (int i=0;i<index->Count;i++)
			{
				ii[index[i]] = L1;
				L1--;
			}
			int L2 = 0;
			for (int i = 0; i < N; i++)
			{
				if (index->Contains(i)) {}
				else {
					ii[index[i]] = L2;
				}
				L2++;
			}
			myPermutation ^perm = gcnew myPermutation(ii);
			this->dat->_mat[0] = this->dat->_mat[0] * perm->p->perm.transpose();
			this->dat->_mat[0] = this->dat->_mat[0].leftCols(L2);			
		}
		void trimRows(System::Collections::Generic::List<int>^ index)
		{
			int N = this->dat->_mat[0].rows();
			int L1 = N;
			array<long long>^ ii = gcnew array<long long>(N);
			for (int i = 0; i < index->Count; i++)
			{
				ii[index[i]] = L1;
				L1--;
			}
			int L2 = 0;
			for (int i = 0; i < N; i++)
			{
				if (index->Contains(i)) {}
				else {
					ii[index[i]] = L2;
				}
				L2++;
			}
			myPermutation^ perm = gcnew myPermutation(ii);
			this->dat->_mat[0] = perm->p->perm*this->dat->_mat[0] ;
			this->dat->_mat[0] = this->dat->_mat[0].topRows(L2);
		}

		/*System::String^ _solve0_lu(myDoubleArray^ rhs, myDoubleArray^ ret, int ordering, bool meh) {
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
			return gcnew System::String(_str.c_str());
		}*/
		void _solve0_lu_cpu(myDoubleArray^ rhs, myDoubleArray^ ret, int ordering, bool meh,double nnn) {
			mySparse^ m = nullptr;
			myDoubleArray^ v = nullptr;
			if (meh) {

				//rhs->_arr->__v = this->dat->_mat[0].transpose() * rhs->_arr->__v;
				
				Eigen::VectorXd _v = this->dat->_mat[0].transpose() * rhs->_arr->__v;
				Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> _m = this->dat->_mat[0].transpose() * this->dat->_mat[0];

				//this->dat->_mat[0] = this->dat->_mat[0].transpose() * this->dat->_mat[0];
				m = gcnew mySparse(_m.rows(), _m.cols());
				m->dat->_mat[0] = _m;
				v = gcnew myDoubleArray(_v.size());
				v->_arr->__v = _v;
			}
			else {
				m = this;
				v = rhs;
			}
			double nn = 0.00000000001;
			bool allocerr = false;
			std::string _str = "";
			if (nnn != 0)
			{
				m->dat->addsmallidentity(nnn, true, false);
				nn = nnn;
			}

			for (int i = 0; i < 10; i++)
			{
				std::string str = m->dat->_solve0_lu_cpu(&v->_arr->__v, &ret->_arr->__v, ordering);

				//lu.compute(this->dat->_mat[0]);
				_str += str;
				if (str.find("SUCCESS") == string::npos)
				{
					nn *= 100;
					m->dat->addsmallidentity(nn, true, false);
				}
				else {
					break;
				}
			}
			if (_str == "")_str = "success";
			System::Console::WriteLine(gcnew System::String(_str.c_str()));
		}
		void _solve0_lu_cpu(myDoubleArray^ rhs, myDoubleArray^ ret, int ordering, bool meh) {
			mySparse^ m = nullptr;
			double nnn = 0;
			myDoubleArray^ v = nullptr;
			if (meh) {

				//rhs->_arr->__v = this->dat->_mat[0].transpose() * rhs->_arr->__v;

				Eigen::VectorXd _v = this->dat->_mat[0].transpose() * rhs->_arr->__v;
				Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> _m = this->dat->_mat[0].transpose() * this->dat->_mat[0];

				//this->dat->_mat[0] = this->dat->_mat[0].transpose() * this->dat->_mat[0];
				m = gcnew mySparse(_m.rows(), _m.cols());
				m->dat->_mat[0] = _m;
				v = gcnew myDoubleArray(_v.size());
				v->_arr->__v = _v;
			}
			else {
				m = this;
				v = rhs;
			}
			double nn = 0.00000000001;
			bool allocerr = false;
			std::string _str = "";
			if (nnn != 0)
			{
				m->dat->addsmallidentity(nnn, true, false);
				nn = nnn;
			}

			for (int i = 0; i < 10; i++)
			{
				std::string str = m->dat->_solve0_lu_cpu(&v->_arr->__v, &ret->_arr->__v, ordering);

				//lu.compute(this->dat->_mat[0]);
				_str += str;
				if (str.find("SUCCESS") == string::npos)
				{
					nn *= 100;
					m->dat->addsmallidentity(nn, true, false);
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
		/*void _solve0(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->_solve0(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}*/
		void __solve0(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//pin_ptr<double> ptr = &rhs[0];

			dat->__solve0(&rhs->_arr->__v, &ret->_arr->__v);

			//array<double>^ ret = gcnew array<double>(_ret.rows());
			//System::Runtime::InteropServices::Marshal::Copy((IntPtr)_ret.data(), ret, 0, _ret.rows());

			//ptr = nullptr;
			//return ret;
		}

		System::String^ _solve_gpu(myCuda^ gpu, myDoubleArray^ rhs, myDoubleArray^ ret, Int64 device) {
			//pin_ptr<double> ptr = &rhs[0];

			auto ss = dat->_solveLU_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);
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
		System::String^ _solveLU_sparse_cpu(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//auto ss = dat->_solveLU_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);
			//System::String^ ee = gcnew System::String(ss.c_str());
			auto ss = dat->_solveLU_sparse_cpu(&rhs->_arr->__v, &ret->_arr->__v);

			return gcnew System::String(ss.c_str());
		}
		System::String^ _solveCG_sparse_cpu(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//auto ss = dat->_solveLU_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);
			//System::String^ ee = gcnew System::String(ss.c_str());
			auto ss = dat->_solveCG_sparse_cpu(&rhs->_arr->__v, &ret->_arr->__v);

			return gcnew System::String(ss.c_str());
		}
		System::String^ _solveLU_dense_cpu(myDoubleArray^ rhs, myDoubleArray^ ret) {
			//auto ss = dat->_solveLU_gpu(gpu->cuda(), &rhs->_arr->__v, &ret->_arr->__v, device);
			//System::String^ ee = gcnew System::String(ss.c_str());
			auto ss = dat->_solveLU_dense_cpu(&rhs->_arr->__v, &ret->_arr->__v);

			return gcnew System::String(ss.c_str());
		}
		void turnDense()
		{
			this->dat->turnDense();
		}
		/*mySparse^ solve0(mySparse^ rhs) {
			Eigen::LLT<Eigen::MatrixXd>* _LLT = new Eigen::LLT<Eigen::MatrixXd>();
			dat->computeLLT(_LLT);
			myLLT^ LLT = gcnew myLLT();
			LLT->LLT->LLT = _LLT;

			Eigen::MatrixXd _ret = this->dat->_solve0(LLT->LLT, rhs->dat);
			mySparse^ ret = gcnew mySparse(rhs->rows(), rhs->cols());

			return ret;
		}*/

		void solve0_gpu(myCuda^ gpu, mySparse^ rhs, mySparse^ ret) {
			this->dat->_solve0_gpu(gpu->cuda(), rhs->dat, ret->dat);
		}
		Int64 solveI(mySparse^ ret) {
			return this->dat->_solveI(ret->dat);
		}
		Int64 solveI_dense(mySparse^ ret) {
			return this->dat->_solveI_dense(ret->dat);
		}
		System::String^ solveI_gpu(myCuda^ gpu, mySparse^ ret) {
			auto _ss = this->dat->_solveI_gpu_single(gpu->cuda(), ret->dat);
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
		/*System::String^ AinvBA(myCuda^ gpu, mySparse^ A, mySparse^ ret) {
			std::string ss = this->dat->AinvBA(gpu->cuda(), A->dat, ret->dat);
			auto ee = gcnew System::String(ss.c_str());
			return ee;
		}*/
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
		void setcoeff(myDoubleArray^ a)
		{
			this->dat->coeff[0].resize(a->size());
			int N = a->size();
			for (int i = 0; i < N; i++)
			{
				this->dat->coeff[0](i) = a->at(i);
			}
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
		void plus_dense(mySparse^ m, double a) {
			int _mt;
			int n = Eigen::nbThreads();
			Eigen::setNbThreads(1);
			Eigen::initParallel();
#pragma omp parallel
			{
#pragma omp single
				_mt = omp_get_num_threads();
			}
#pragma omp parallel for
			for (int nn = 0; nn < _mt; nn++)
			{
				int S = (m->dat->_dmat.cols() * nn) / _mt;
				int E = (m->dat->_dmat.cols() * (nn+1)) / _mt;
				this->dat->_dmat.middleCols(S,E-S) += a * m->dat->_dmat.middleCols(S,E-S);
			}
			Eigen::setNbThreads(0);

		}

		void plus_AtA_dense(myDoubleArray^ m, double a) {
			int _mt;
			int n = Eigen::nbThreads();
			Eigen::setNbThreads(1);
			Eigen::initParallel();
#pragma omp parallel
			{
#pragma omp single
				_mt = omp_get_num_threads();
			}
#pragma omp parallel for
			for (int nn = 0; nn < _mt; nn++)
			{
				int S = (this->dat->_dmat.cols() * nn) / _mt;
				int E = (this->dat->_dmat.cols() * (nn + 1)) / _mt;
				this->dat->_dmat.middleCols(S, E - S) += a * m->_arr->__v * m->_arr->__v.middleRows(S, E - S).transpose();
			}
			Eigen::setNbThreads(0);

		}
		void addsmallidentity(double salt, bool sparse, bool dense) {
			this->dat->addsmallidentity(salt, sparse, dense);
		}
		void addsmallidentity(double salt, bool sparse, bool dense,int m) {
			this->dat->addsmallidentity(salt, sparse, dense,m);
		}
		void addsmallones(double salt)
		{
			this->dat->_dmat += salt * Eigen::MatrixXd::Ones(this->dat->_dmat.rows(), this->dat->_dmat.cols());
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
	
	public ref class sparseMatrix {
	private:
		Eigen::SparseMatrix<double>* mat = 0;
	public:
		inline Eigen::SparseMatrix<double>& get()
		{
			return *mat;
		}
		inline void set(Eigen::SparseMatrix<double> _m)
		{
			*mat = _m;
		}
		void resize(int n, int m)
		{
			mat->resize(n, m);
		}
		int rows() {
			return mat->rows();
		}
		int cols() {
			return mat->cols();
		}
		sparseMatrix(int n, int m)
		{
			if (mat != 0)del();
			mat = new Eigen::SparseMatrix<double>();
			mat->resize(n, m);
			mat->setZero();
		}
		sparseMatrix()
		{
			if (mat != 0)del();
			mat = new Eigen::SparseMatrix<double>();
		}
		!sparseMatrix()
		{
			del();
		}
		~sparseMatrix()
		{
			del();
		}
		void del() {
			if (mat != 0)
			{
				delete mat;
			}
			mat = 0;
		}
	};
	template<typename Scalar, typename StorageIndex = typename Eigen::SparseMatrix<Scalar>::StorageIndex >
	class _Triplet
	{
	public:
		_Triplet() : m_row(0), m_col(0), m_value(0) {}

		_Triplet(const StorageIndex& i, const StorageIndex& j, const Scalar& v = Scalar(0))
			: m_row(i), m_col(j), m_value(v)
		{}

		/** \returns the row index of the element */
		const StorageIndex& row()const { return m_row; }

		/** \returns the column index of the element */
		const StorageIndex& col() const { return m_col; }

		/** \returns the value of the element */
		const Scalar& value() const { return m_value; }
		/** \returns the row index of the element */
		StorageIndex& _row() { return m_row; }

		/** \returns the column index of the element */
		StorageIndex& _col() { return m_col; }

		/** \returns the value of the element */
		Scalar& _value() { return m_value; }
	public:
		StorageIndex m_row, m_col;
		Scalar m_value;
	};
	public ref class workspace {
	public:
		std::vector<_Triplet<double>>* _dat=0;
		workspace()
		{
			if (_dat != 0)del();
			_dat = new std::vector<_Triplet<double>>();
			_dat->clear();
		}
		System::String^ tostring()
		{
			System::String^ str = gcnew System::String("");

			for (auto dd : *_dat)
			{
				str = str + dd.row().ToString() + "," + dd.col().ToString() + "," + dd.value().ToString()+"\n";
			}
			return str;
		}
		Eigen::SparseMatrix<double> _tosparse(int n)
		{
			Eigen::SparseMatrix<double> newmat(n, n);
			for (auto tr : *_dat)
			{
				newmat.coeffRef(tr.row(), tr.col()) += tr.value();
			}
			return newmat;
		}
		mySparse^ __tosparse(int n)
		{
			auto mat = _tosparse(n);
			mySparse^ ret = gcnew mySparse(n, n);
			ret->dat->_mat.resize(1);
			ret->dat->_mat[0] = mat;
			return ret;
		}
		
		Eigen::SparseMatrix<double> tosparse(int n, bool flip, bool both)
		{
			Eigen::SparseMatrix<double> newmat(n, n);
			if (both)
			{
				if ((*_dat)[0].col() == -1)
				{
					newmat.resize(n, 1);
					newmat.setZero();
					newmat.reserve(_dat->size());
					
					for (auto tr : *_dat)
					{
						newmat.coeffRef(tr.row(),0) = tr.value();
					}
					
				}
				else if ((*_dat)[0].row() == -1)
				{
					newmat.resize(1, n);
					newmat.setZero();
					newmat.reserve(_dat->size());
					for (auto tr : *_dat)
					{
						newmat.coeffRef(0, tr.col()) = tr.value();
					}
				}
				else {
					newmat.setZero();
					newmat.reserve(_dat->size());
					for (auto tr : *_dat)
					{
						newmat.coeffRef(tr.row(), tr.col()) = tr.value();
					}
				}
			}
			else {
				if (flip)//U,oV
				{
					if ((*_dat)[0].col() == -1)
					{
						newmat.resize(1, n);
						newmat.setZero();
						newmat.reserve(_dat->size());
						for (auto tr : *_dat)
						{
							newmat.coeffRef(0, tr.row()) = tr.value();
						}
					}
					else if ((*_dat)[0].row() == -1)
					{
						newmat.resize(0, 0);
						newmat.setZero();
					}
					else {
						newmat.setZero();
						newmat.reserve(_dat->size());
						for (auto tr : *_dat)
						{
							newmat.coeffRef(tr.col(), tr.row()) = tr.value();
						}
					}
				}
				else {//U,oV
					if ((*_dat)[0].col() == -1)
					{
						newmat.resize(0, 0);
						newmat.setZero();
					}
					else if ((*_dat)[0].row() == -1)
					{
						newmat.resize(1, n);
						newmat.setZero();
						newmat.reserve(_dat->size());
						for (auto tr : *_dat)
						{
							newmat.coeffRef(0, tr.col()) = tr.value();
						}
					}
					else {
						newmat.setZero();
						newmat.reserve(_dat->size());
						for (auto tr : *_dat)
						{
							newmat.coeffRef(tr.row(), tr.col()) = tr.value();
						}
					}
				}
			}
			return newmat;
		}
		Eigen::MatrixXd contract(denseMatrix^ U, denseMatrix^ V, int n,int m,bool flip,bool both)
		{
			auto _V = V->get();
			auto _U = U->get();
			Eigen::MatrixXd newmat(n,n);
			if (both)
			{
				if ((*_dat)[0].col() == -1)
				{
					newmat.resize(n, 1);
					newmat.setZero();
					for (int j = 0; j < n; j++)
					{
						for (auto tr : *_dat)
						{
							newmat.coeffRef(j, 0) += _U.coeffRef(tr.row(), j) * tr.value();
						}
					}
				}
				else if ((*_dat)[0].row() == -1)
				{
					newmat.resize(1, n);
					newmat.setZero();
					for (int j = 0; j < n; j++)
					{
						for (auto tr : *_dat)
						{
							newmat.coeffRef(0, j) += _V.coeffRef(tr.col(), j) * tr.value();
						}
					}
				}
				else {
					newmat.setZero();
					for (int i = 0; i < n; i++)
					{
						for (int j = 0; j < n; j++)
						{
							for (auto tr : *_dat)
							{
								newmat.coeffRef(i, j) += _V.coeffRef(tr.col(), j) * _U.coeffRef(tr.row(), i) * tr.value();
							}
						}
					}
				}
			}
			else {
				if (flip)//U,oV
				{
					if ((*_dat)[0].col() == -1)
					{
						newmat.resize(1, n);
						newmat.setZero();
						for (int j = 0; j < n; j++)
						{
							for (auto tr : *_dat)
							{
								newmat.coeffRef(0, j) += _U.coeffRef(tr.row(), j) * tr.value();
							}
						}
					}
					else if ((*_dat)[0].row() == -1)
					{
						newmat.resize(0, 0);
						newmat.setZero();
					}
					else {
						newmat.setZero();
						for (int i = 0; i < n; i++)
						{
							for (int j = 0; j < n; j++)
							{
								for (auto tr : *_dat)
								{
									newmat.coeffRef(j, i) += _V.coeffRef(tr.col(), j) * _U.coeffRef(tr.row(), i) * tr.value();
								}
							}
						}
					}
				}
				else {//U,oV
					if ((*_dat)[0].col() == -1)
					{
						newmat.resize(0, 0);
						newmat.setZero();
					}
					else if ((*_dat)[0].row() == -1)
					{
						newmat.resize(1, n);
						newmat.setZero();
						for (int j = 0; j < n; j++)
						{
							for (auto tr : *_dat)
							{
								newmat.coeffRef(0, j) += _V.coeffRef(tr.col(), j) * tr.value();
							}
						}
					}
					else {
						newmat.setZero();
						for (int i = 0; i < n; i++)
						{
							for (int j = 0; j < n; j++)
							{
								for (auto tr : *_dat)
								{
									newmat.coeffRef(i, j) += _V.coeffRef(tr.col(), j) * _U.coeffRef(tr.row(), i) * tr.value();
								}
							}
						}
					}
				}
			}
			/*if (U == nullptr)
			{
				auto _V = V->get();
				int nv = _V.rows();
				if ((*_dat)[0].row() == -1)
				{
					newmat.resize(1, n);
					newmat.setZero();
					for (int j = 0; j < n; j++)
					{
						for (auto tr : *_dat)
						{
							//if (tr.col() < nv)
							{
								newmat.coeffRef(0, j) += _V.coeffRef(tr.col(), j) * tr.value();
							}
						}

					}
				}else if ((*_dat)[0].col() == -1)
				{
					newmat.resize(0, 0);
				}
				else {
					newmat.resize(m, n);
					newmat.setZero();
					int nv = _V.rows();
					for (int i = 0; i < n; i++)
					{
						for (auto tr : *_dat)
						{
							//if (tr.col() < nv)
							{
								newmat.coeffRef(tr.row(), i) += _V.coeffRef(tr.col(), i) * tr.value();
							}
						}
					}
				}
			}
			else if (V == nullptr) {
				auto _U = U->get();
				if ((*_dat)[0].col() == -1)
				{
					newmat.resize(1,n);
					newmat.setZero();

					int nu = _U.rows();
					for (int i = 0; i < n; i++)
					{
						for (auto tr : *_dat)
						{
							//if (tr.row() < nu)
							{
								newmat.coeffRef(0, i) += _U.coeffRef(tr.row(), i) * tr.value();
							}
						}
					}
				}else if ((*_dat)[0].row() == -1)
				{
					newmat.resize(0, 0);
				}
				else {
					newmat.resize(m,n);
					newmat.setZero();
					int nu = _U.rows();
					for (int i = 0; i < n; i++)
					{
						for (auto tr : *_dat)
						{
							//if (tr.row() < nu)
							{
								newmat.coeffRef(tr.col(), i) += _U.coeffRef(tr.row(), i) * tr.value();
							}
						}
					}
				}
			}
			*/
			return newmat;

		}
		void pushforward(myPermutation^ mU, myPermutation^ mV,int nU,int nV,workspace^ newone)
		{
			if ((*_dat).empty())return;
			auto pU = mU->p->perm;
			auto pV = mV->p->perm;
			(*newone->_dat).resize((*_dat).size());

			if ((*_dat)[0].row() == -1)
			{
				for (auto& triplet : (*newone->_dat))
				{
					triplet.m_col = pV.indices()(triplet.col());
				}
			}
			else if((*_dat)[0].col() == -1)
			{
				for (auto& triplet : (*newone->_dat))
				{
					triplet.m_row = pU.indices()(triplet.row());
				}
			}
			else {
				for (auto& triplet : (*newone->_dat))
				{
					triplet.m_col = pV.indices()(triplet.col());
					triplet.m_row = pU.indices()(triplet.row());
				}
			}
			
			/*for (auto it = (*newone->_dat).begin(); it != (*newone->_dat).end();) {
				if ((*it).row()>=nU|| (*it).col() >= nV) {
					it = (*newone->_dat).erase(it);
				}
				else {
					++it;
				}
			}*/
		}
		
		void mulright(Eigen::VectorXd* v, Eigen::VectorXd* ret,double sc,int jj)
		{
			if ((*_dat).empty())return;
			if ((*_dat)[0].row() == -1)return;
			if ((*_dat)[0].col() == -1)
			{
				//if (jj % 2 == 1)
				{
					for (const auto& triplet : *_dat)
					{
						ret->coeffRef(triplet.row()) += sc * (triplet.value());
					}
				}
			}
			else {
				//if (jj % 2 == 0)
				{
					for (const auto& triplet : *_dat)
					{
						ret->coeffRef(triplet.row()) += sc * ((*v)(triplet.col()) * triplet.value());
					}
				}
			}
		}
		void mulleft(Eigen::VectorXd* u, Eigen::VectorXd* ret,double sc,int jj)
		{
			if ((*_dat).empty())return;
			if ((*_dat)[0].col() == -1)return;
			if ((*_dat)[0].row() == -1)
			{
				//if (jj % 2 == 1)
				{
					for (const auto& triplet : *_dat)
					{
						ret->coeffRef(triplet.col()) += sc * (triplet.value());
					}
				}
			}
			else {
				//if (jj % 2 == 0)
				{
					for (const auto& triplet : *_dat)
					{
						ret->coeffRef(triplet.col()) += sc * ((*u)(triplet.row()) * triplet.value());
					}
				}
			}
		}
		double mulboth(Eigen::VectorXd* u, Eigen::VectorXd* v,int jj)
		{
			double val = 0;
			if ((*_dat).empty())return 0;
			if ((*_dat)[0].row() == -1)
			{
				//if (jj % 2 == 1)
				{
					for (const auto& triplet : *_dat)
					{
						val += ((*v)(triplet.col())) * triplet.value();
					}
				}
			}
			else if ((*_dat)[0].col() == -1)
			{
				//if (jj % 2 == 1)
				{
					for (const auto& triplet : *_dat)
					{
						val += ((*u)(triplet.row())) * triplet.value();
					}
				}
			}
			else {
				//if (jj % 2 == 0)
				{
					for (const auto& triplet : *_dat)
					{
						val += ((*u)(triplet.row())) * ((*v)(triplet.col())) * triplet.value();
					}
				}
			}
			return val;
		}
		~workspace() {
			if (_dat != 0)
				del();
			_dat = 0;
		}
		!workspace() {
			if (_dat != 0)
				del();
			_dat = 0;
		}
		void del() {
			if (_dat != 0)delete _dat;
			_dat = 0;
		}
	};
	public ref class helper {
	public:
		static void orthogonalize(myDoubleArray^ _a, myDoubleArray^ _b, myDoubleArray^ ret)
		{
			double norm = _a->_arr->__v.dot(_b->_arr->__v);
			double normB = _b->_arr->__v.squaredNorm();
			if (normB < 0.00000000001)
			{
				ret->_arr->__v = _a->_arr->__v;
			}
			else {
				ret->_arr->__v = _a->_arr->__v - norm / normB * _b->_arr->__v;
			}

			
		}
		static void find(myDoubleArray^ dx0,myDoubleArray^ _b, mySparse^ _F,myDoubleArray^ ret,double alpha)
		{
			int N = _b->_arr->__v.size();
			Eigen::MatrixXd F = _F->dat->_dmat;
			Eigen::MatrixXd FIFI = F;// 0.5 * (F.transpose() + F);
			Eigen::MatrixXd M = F + alpha * Eigen::MatrixXd::Identity(N, N);


			Eigen::VectorXd dx = _b->_arr->__v;
			Eigen::VectorXd prev_dx = _b->_arr->__v;
			Eigen::VectorXd b = _b->_arr->__v;
			double normb = b.norm();
			Eigen::VectorXd c;
			if (_b->L2Norm() < 0.000000000001)return;
			for (int i = 0; i < 50; i++)
			{
				prev_dx = dx;
				dx = M * dx;
				c = b - 2.0 * dx;
				//solve (dx+lambda c)(dx+lambda c)=(dx+lambda c)b
				double A = c.dot(c);
				double B = 2 * dx.dot(c)-c.dot(b);
				double C = dx.dot(dx) - dx.dot(b);
				double lambda1 = (-B + sqrt(B * B - 4 * A * C)) / (2 * A);
				double lambda2 = (-B - sqrt(B * B - 4 * A * C)) / (2 * A);
				if (abs(lambda1) < abs(lambda2))
				{
					dx = dx + lambda1 * c;
				}
				else {
					dx = dx + lambda2 * c;
				}
				//double norm = dx.norm();
				//dx = dx+b/normb*norm*0.1;
				if (dx.norm() < 0.000000000001)break;
				if ((prev_dx - dx).norm() < 0.00000000001)break;
			}
			ret->_arr->__v = dx;
		}
		static double computeeigen(mySparse^ i1, mySparse^ i2, mySparse^ t1, Int64 N)
		{
			double maxval = _mySparse::computeeigen(i1->dat, i2->dat, t1->dat, N);


			return maxval;
			//auto vv = mat.eigenvalues();
		}
		static void project(mySparse^ i1, mySparse^ i2, mySparse^ t1, myDoubleArray^ v1, myDoubleArray^ v2, myDoubleArray^ ret1, myDoubleArray^ ret2)
		{
			_mySparse::project(i1->dat, i2->dat, t1->dat, v1->_arr, v2->_arr, ret1->_arr, ret2->_arr);
		}
		static void GN(mySparse^ mat1, mySparse^ mat2, mySparse^ mat3, myDoubleArray^ rhs1, myDoubleArray^ rhs2, myDoubleArray^ ret1, myDoubleArray^ ret2,int L1phi,int L1Z,myCuda ^cuda)
		{
			Eigen::MatrixXd M(L1phi + L1Z, L1phi + L1Z);
			M.topLeftCorner(L1phi, L1phi) = mat1->dat->_mat[0];
			M.bottomRightCorner(L1Z, L1Z) = mat2->dat->_mat[0];
			M.topRightCorner(L1phi, L1Z) = mat3->dat->_mat[0];
			M.bottomLeftCorner(L1Z, L1phi) = mat3->dat->_mat[0].transpose();

			M += Eigen::MatrixXd::Identity(L1phi + L1Z, L1phi + L1Z) * 0.000000000000001;
			Eigen::VectorXd rhs(L1phi + L1Z);
			rhs.topRows(L1phi) = rhs1->_arr->__v;
			rhs.bottomRows(L1Z) = rhs2->_arr->__v;

			mySparse^ m = gcnew mySparse();
			m->dat->_dmat = M;
			myDoubleArray^ _rhs = gcnew myDoubleArray(L1phi+L1Z);
			_rhs->_arr->__v = rhs;
			myDoubleArray^ _ret = gcnew myDoubleArray(L1phi + L1Z);
			m->_solveLU_gpu(cuda,_rhs, _ret,cuda->fastest());
			ret1->_arr->__v = _ret->_arr->__v.topRows(L1phi);
			ret2->_arr->__v = _ret->_arr->__v.bottomRows(L1Z);
			/*Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
			auto ret=lu.solve(rhs);

			ret1->_arr->__v = ret.topRows(L1phi);
			ret2->_arr->__v = ret.bottomRows(L1Z);
			*/

		}
		static void GN2(mySparse^ mat1, mySparse^ mat2, mySparse^ mat3, myDoubleArray^ rhs1, myDoubleArray^ rhs2, myDoubleArray^ ret1, myDoubleArray^ ret2, int L1phi, int L1Z, myCuda^ cuda)
		{
			Eigen::MatrixXd M(L1phi + L1Z, L1phi + L1Z);
			M.topLeftCorner(L1phi, L1phi) = mat1->dat->_mat[0];
			M.topRightCorner(L1phi, L1Z) = mat3->dat->_mat[0];

			//M.bottomLeftCorner(L1Z, L1phi) = mat3->dat->_mat[0].transpose();
			//M.bottomRightCorner(L1Z, L1Z) = mat2->dat->_mat[0];
			M += Eigen::MatrixXd::Identity(L1phi + L1Z, L1phi + L1Z) * 0.000000000000001;
			Eigen::VectorXd rhs(L1phi+L1Z);
			rhs.topRows(L1phi) = rhs1->_arr->__v;
			rhs.bottomRows(L1Z) = rhs2->_arr->__v;

			/*mySparse^ m = gcnew mySparse();
			m->dat->_dmat = M;
			myDoubleArray^ _rhs = gcnew myDoubleArray(L1phi + L1Z);
			_rhs->_arr->__v = rhs;
			myDoubleArray^ _ret = gcnew myDoubleArray(L1phi + L1Z);
			m->_solveLU_gpu(cuda, _rhs, _ret, cuda->fastest());
			ret1->_arr->__v = _ret->_arr->__v.topRows(L1phi);
			ret2->_arr->__v = _ret->_arr->__v.bottomRows(L1Z);
			*/
			Eigen::FullPivLU<Eigen::MatrixXd> lu(M);
			auto ret=lu.solve(rhs);
			ret1->_arr->__v = ret.topRows(L1phi);
			ret2->_arr->__v = ret.bottomRows(L1Z);
			

		}
		//helper.pushforward(_mats,mZ,mphi,L1Z,L1phi);
		//helper.computeKrylovSubspace(_mats, __U, __V, __W, _C, 200);

		/*static System::Collections::Generic::List<workspace^>^ pushforward(System::Collections::Generic::List<workspace^>^ _mats, myPermutation^ mU, myPermutation^ mV, int nU, int nV)
		{
			System::Collections::Generic::List<workspace^>^ newmats = gcnew System::Collections::Generic::List<workspace^>();
			int _mt = omp_get_max_threads();
			int m = _mats->Count;
			for (int i = 0; i < m; i++)
			{
				newmats->Add(nullptr);
			}

#pragma omp parallel
			{
#pragma omp single
				_mt = omp_get_num_threads();
			}
#pragma omp parallel for
			for (int ii = 0; ii < _mt; ii++)
			{
				int S = ii * m / _mt;
				int E = (ii + 1) * m / _mt;
				for (int i = S; i < E; i++)
				{
					auto M = _mats[i];
					auto newworkspace = gcnew workspace();
					M->pushforward(mU, mV,nU,nV,newworkspace);
					newmats[i] = newworkspace;
				}
			}
			return newmats;
		}*/
		//note that _mats is the one before pushforward
		//The one used for Krylov decompoision was used to compute U,V,W and now they are pulled back. It is safe to multiply the original _mats and U,V,W.
		static void contract(System::Collections::Generic::List<workspace^>^ _mats, denseMatrix^ U, denseMatrix^ V, denseMatrix^ W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, int _C);
		
		/*static void pullback(denseMatrix^ _U, denseMatrix^ _V, myPermutation^ mU, myPermutation^ mV, int nU, int nV)
		{
			int nu0 = _U->get().rows();
			int nv0 = _V->get().rows();

			_U->get().conservativeResize(nU, _U->get().cols());
			_V->get().conservativeResize(nV, _V->get().cols());
			_U->get().bottomRows(nU - nu0).setZero();
			_V->get().bottomRows(nV - nv0).setZero();

			_U->get().applyOnTheLeft(mU->p->perm.transpose());
			_V->get().applyOnTheLeft(mV->p->perm.transpose());
			
		}*/
		static double VarPro(System::Collections::Generic::List<double>^ __coeff, myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt);
		static double GN(System::Collections::Generic::List<double>^ __coeff, myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt,myCuda ^ cuda);
		static double ALT(System::Collections::Generic::List<double>^ __coeff, myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt);
		static double simple(System::Collections::Generic::List<double>^ __coeff, myDoubleArray^ phi, myDoubleArray^ zz, denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt);
		static void write(System::Collections::Generic::List<double>^ __coeff, myDoubleArray^ phi0, myDoubleArray^ zz0, myDoubleArray^ phi, myDoubleArray^ zz,denseMatrix^ __U, denseMatrix^ __V, denseMatrix^ __W, array<sparseMatrix^>^ _mats1, array<sparseMatrix^>^ _mats2, array<sparseMatrix^>^ _mats3, myDoubleArray^ _r1, myDoubleArray^ _r2, double dt, int tt);

		static void computeKrylovSubspace(System::Collections::Generic::List<double>^ __coeff,System::Collections::Generic::List<workspace^>^ _mats, denseMatrix^ _U, denseMatrix^ _V, denseMatrix^ _W, int nU, int nV, int r, myPermutation^ mphi, myPermutation^ mZ, myDoubleArray^ phi, myDoubleArray^ zz, System::Collections::Generic::List<Tuple<int, int>^>^bb1, System::Collections::Generic::List<Tuple<int, int>^>^ bb2);
		

	};

	public ref class myMicroMatrix {
		_myMicroMatrix* _dat = 0;
		int N;
	public:
		myMicroMatrix()
		{
			_dat = new _myMicroMatrix();
			_dat->_sc = 1.0;
		}
		~myMicroMatrix()
		{
			if (_dat != 0)
			{
				delete _dat;
			}
			_dat = 0;
		}
		!myMicroMatrix()
		{
			if (_dat != 0)
			{
				delete _dat;
			}
			_dat = 0;
		}
		double check() {
			return (_dat->_mat/* - _dat->_mat.transpose()*/).squaredNorm();

		}
		void scale(double sc)
		{
			_dat->_sc = sc;
			_dat->sqrsc = sqrt(sc);

		}
		void init(int I, System::Collections::Generic::List<int>^ star)
		{
			this->_dat->indices.push_back(I);
			for (int i = 0; i < star->Count; i++)
			{
				this->_dat->indices.push_back(star[i]);
			}
			N = star->Count + 1;
			this->_dat->_mat.resize(N, N);
			this->_dat->_mat.setZero();
			this->_dat->_w.resize(N);
			this->_dat->_w.setZero();
			this->_dat->_z.resize(N);
			this->_dat->_z.setZero();
		}
		void set(int i, int j, double val)
		{
			this->_dat->_mat(i, j) = val;
		}
		void updateNodes(myDoubleArray^ z, myDoubleArray^ w)
		{
			for (int i = 0; i < N; i++)
			{
				this->_dat->_w(i) = w->_arr->__v(this->_dat->indices[i]);
				this->_dat->_z(i) = z->_arr->__v(this->_dat->indices[i]);
			}
		}

		void computeJacobian(mySparse^ M, int I, bool z,bool firsttime)
		{
			
			if (z)
			{
				this->_dat->_tmp = ((this->_dat->_w.transpose() * this->_dat->_mat).transpose()) * _dat->sqrsc;
			}
			else {
				this->_dat->_tmp = ((this->_dat->_mat * this->_dat->_z)) * _dat->sqrsc;
			}
			if (firsttime)
			{
				for (int i = 0; i < N; i++)
				{
					M->dat->dat[0].push_back(Eigen::Triplet<double>(I, this->_dat->indices[i], this->_dat->_tmp(i, 0)));
				}
			}
			else {
				for (int i = 0; i < N; i++)
				{
					M->dat->_mat[0].coeffRef(I, this->_dat->indices[i]) = this->_dat->_tmp(i, 0);
				}
			}
		}
		double computeResidual(double rho)
		{
			return _dat->_sc * ((this->_dat->_w.transpose() * this->_dat->_mat *this->_dat->_z)(0,0)+rho);
		}
		double getscale() {
			return _dat->_sc;
		}


	};
}