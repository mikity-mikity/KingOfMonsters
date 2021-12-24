#pragma once
#include <cmath>
#include<vector>

using namespace System;
#include <cstring>
#include <string>
//using namespace Alea;
//using namespace Alea::Parallel;
using std::vector;
using std::string;
#include "mySparseLibrary.h"
namespace KingOfMonsters {
	
	public class _buffer {
	public
		double mem[12000];
	};
	public ref class buffer {
	public:
		_buffer* _buf=0;

		buffer() {
			_buf = new _buffer();
		}
		~buffer() {
			if(_buf!=0)
				delete _buf;
			_buf = 0;
		}
		!buffer() {
			if(_buf!=0)
				delete _buf;
			_buf = 0;
		}
	};
	public class _memS_ref {
	public:
		int _nNode;
		int dim[2]{ 0,0 };
		int _uDim, _vDim;
		double* node=0;
		double* def = 0;
		double* buf_z = 0;
		double* buf_phi = 0;
		double* buf_b = 0;
		double* buf_D = 0;

		double _gi[6];
		double _Gi[6];
		double _gij[4];
		double _Gij[4];
		double _bij[12];
		double _Sij[4];
		double _Gammaijk[8];
		double* __mat = 0;
		double** M[2]{ 0,0 };
		int** dd=0;
		double* d0=0;
		double* d1[2]{ 0,0 };
		double* d2[4]{ 0,0,0,0 };
		double* d2_star[4]{ 0,0,0,0 };
		const int star2[4]{ 1,-1,-1,1 };
		const int _star[4]{ 3,1,2,0 };
		/*star2[0] = 1;
		star2[1] = -1;
		star2[2] = -1;
		star2[3] = 1;

		_star[0] = 3;
		_star[1] = 1;
		_star[2] = 2;
		_star[3] = 0;*/
		double* B[4]{ 0,0,0,0 };
		double* tt0[2]{ 0,0 }, * hh0[2]{ 0,0 }, * tt1[4]{ 0,0,0,0 }, * hh1[4]{ 0,0,0,0 }, * tt2[8]{ 0,0,0,0,0,0,0,0 }, * hh2[8]{ 0,0,0,0,0,0,0,0 };
		const int ___ll[4]{ 0,3,6,9 };
	public:
		bool initialized=false;
		double refDv=0,_refDv=0;
		double _x=0, _y=0, _z=0, __z=0, Z=0, _Z=0;
		inline void set_z(double &z) {
			this->_z = z;
		}
		inline void set__z(double &z) {
			this->__z = z;
		}
		inline void set_buffer(double* buf)
		{
			node = buf;
			def = &buf[2000];
			buf_z = &buf[4000];
			buf_phi = &buf[6000];
			buf_b = &buf[8000];
			buf_D = &buf[10000];
		}
		inline void set_node(const int &i, const int &s, const double val) {
			node[___ll[i] + s] = val;
		}
		inline void set_buf_z(const int &i, const double &val) {
			buf_z[i] = val;
		}
		inline void set_node(double* ptr, const int &N) {
			//buf_z[i] = val;
			memcpy(node, ptr, sizeof(double) * N);
		}
		inline void set_buf_z(double* ptr, const int &N) {
			//buf_z[i] = val;
			memcpy(buf_z, ptr, sizeof(double) * N);
		}
		inline void set_buf_phi(double* ptr, const int &N) {
			//buf_z[i] = val;
			memcpy(buf_phi, ptr, sizeof(double) * N);
		}
		inline void set_def(double* ptr, const int &N) {
			//buf_z[i] = val;
			memcpy(def, ptr, sizeof(double) * N);
		}
		inline void set_buf_phi(const int &i, const double val) {
			buf_phi[i] = val;
		}
		inline void set_def(const int &i, const int &s, const double &val) {
			def[___ll[i]  + s] = val;
		}
		inline double get_node(const int &i, const int &s) {
			return node[___ll[i] + s];
		}
		inline double get__gi(const int &i, const int &s) {
			return _gi[___ll[i] + s];
		}
		inline double get__Gi(const int &i, const int &s) {
			return _Gi[___ll[i] + s];
		}
		inline double get__gij(const int &i, const int &j) {
			return _gij[(i<<1) + j];
		}
		inline double get__Gij(const int &i, const int &j) {
			return _Gij[(i<<1) + j];
		}
		inline double get__bij(const int &i, const int &j, const int &s) {
			return _bij[___ll[((i<<1) + j)] + s];
		}
		inline double get__Gammaijk(const int &i, const int &j, const int &k) {
			return _Gammaijk[(((i<<1) + j) <<1) + k];
		}
	public:
		_memS_ref() {
			initialized = false;
			node = 0;
			def = 0;
			_nNode = 0;
			__z = -10000000;

			M[0] = 0;
			M[1] = 0;
			dd = 0;

			d0 = 0;
			d1[0] = 0;
			d1[1] = 0;
			d2[0] = 0;
			d2[1] = 0;
			d2[2] = 0;
			d2[3] = 0;
			d2_star[0] = 0;
			d2_star[1] = 0;
			d2_star[2] = 0;
			d2_star[3] = 0;



			B[0] = 0;
			B[1] = 0;
			B[2] = 0;
			B[3] = 0;
			tt0[0] = 0;
			tt0[1] = 0;
			hh0[0] = 0;
			hh0[1] = 0;
			tt1[0] = 0;
			tt1[1] = 0;
			tt1[2] = 0;
			tt1[3] = 0;
			hh1[0] = 0;
			hh1[1] = 0;
			hh1[2] = 0;
			hh1[3] = 0;
			tt2[0] = 0;
			tt2[1] = 0;
			tt2[2] = 0;
			tt2[3] = 0;
			tt2[4] = 0;
			tt2[5] = 0;
			tt2[6] = 0;
			tt2[7] = 0;
			hh2[0] = 0;
			hh2[1] = 0;
			hh2[2] = 0;
			hh2[3] = 0;
			hh2[4] = 0;
			hh2[5] = 0;
			hh2[6] = 0;
			hh2[7] = 0;
		}
		~_memS_ref() {
			del();
		}
		
		void del() {		
			if (__mat != 0)delete[] __mat;
			__mat = 0;
			if (dd != 0)
			{
				for (int i = 0; i < _nNode; i++) {
					delete[] dd[i];
				}
				delete[] dd;
				dd = 0;
			}
			if (d0 != 0) {
				delete[] d0;
				d0 = 0;
			}
			if (d1[0] != 0)
			{
				delete[] d1[0];
				delete[] d1[1];
				d1[0] = 0;
				d1[1] = 0;
			}
			if (d2[0] != 0)
			{
				delete[] d2[0];
				delete[] d2[1];
				delete[] d2[2];
				delete[] d2[3];
				d2[0] = 0;
				d2[1] = 0;
				d2[2] = 0;
				d2[3] = 0;
				delete[] d2_star[0];
				delete[] d2_star[1];
				delete[] d2_star[2];
				delete[] d2_star[3];
				d2_star[0] = 0;
				d2_star[1] = 0;
				d2_star[2] = 0;
				d2_star[3] = 0;
			}

			if (B[0] != 0) {
				delete[] B[0];
				delete[] B[1];
				delete[] B[2];
				delete[] B[3];
				B[0] = 0;
				B[1] = 0;
				B[2] = 0;
				B[3] = 0;

			}
			if (tt0[0] != 0) {
				delete[] tt0[0];
				tt0[0] = 0;
				delete[] tt0[1];
				tt0[1] = 0;
				delete[] hh0[0];
				delete[] hh0[1];
				hh0[0] = 0;
				hh0[1] = 0;
				delete[] tt1[0];
				delete[] tt1[1];
				delete[] tt1[2];
				delete[] tt1[3];
				tt1[0] = 0;
				tt1[1] = 0;
				tt1[2] = 0;
				tt1[3] = 0;
				delete[] hh1[0];
				delete[] hh1[1];
				delete[] hh1[2];
				delete[] hh1[3];
				hh1[0] = 0;
				hh1[1] = 0;
				hh1[2] = 0;
				hh1[3] = 0;
				delete[] tt2[0];
				delete[] tt2[1];
				delete[] tt2[2];
				delete[] tt2[3];
				delete[] tt2[4];
				delete[] tt2[5];
				delete[] tt2[6];
				delete[] tt2[7];
				tt2[0] = 0;
				tt2[1] = 0;
				tt2[2] = 0;
				tt2[3] = 0;
				tt2[4] = 0;
				tt2[5] = 0;
				tt2[6] = 0;
				tt2[7] = 0;

				delete[] hh2[0];
				delete[] hh2[1];
				delete[] hh2[2];
				delete[] hh2[3];
				delete[] hh2[4];
				delete[] hh2[5];
				delete[] hh2[6];
				delete[] hh2[7];
				hh2[0] = 0;
				hh2[1] = 0;
				hh2[2] = 0;
				hh2[3] = 0;
				hh2[4] = 0;
				hh2[5] = 0;
				hh2[6] = 0;
				hh2[7] = 0;

			}
			if (M[0] != 0) {
				for (int i = 0; i < _uDim; i++) {
					delete[] M[0][i];
				}
				for (int i = 0; i < _vDim; i++) {
					delete[] M[1][i];
				}
				delete[] M[0];
				delete[] M[1];
				M[0] = 0;
				M[1] = 0;
			}
			/*_nNode = 0;
			__z = -10000000;
			_nNode = 0;
			_uDim = -1;
			_vDim = -1;*/
		}
		void update(int nNode, int uDim, int vDim) {
			if (_nNode != nNode || _uDim != uDim || _vDim != vDim) {
				initialized = false;

				del();
				_nNode = nNode;
				_uDim = uDim;
				_vDim = vDim;
				dim[0] = _uDim;
				dim[1] = _vDim;
				//__mat = new double[nNode * nNode];
				d0 = new double[nNode];

				d1[0] = new double[nNode];
				d1[1] = new double[nNode];

				d2[0] = new double[nNode];
				d2[1] = new double[nNode];
				d2[2] = new double[nNode];
				d2[3] = new double[nNode];

				d2_star[0] = new double[nNode];
				d2_star[1] = new double[nNode];
				d2_star[2] = new double[nNode];
				d2_star[3] = new double[nNode];

				B[0] = new double[nNode * nNode];
				B[1] = new double[nNode * nNode];
				B[2] = new double[nNode * nNode];
				B[3] = new double[nNode * nNode];

				tt0[0] = new double[uDim];
				tt0[1] = new double[vDim];

				hh0[0] = new double[uDim];
				hh0[1] = new double[vDim];

				tt1[0] = new double[uDim];
				tt1[1] = new double[vDim];
				tt1[2] = new double[uDim];
				tt1[3] = new double[vDim];
				hh1[0] = new double[uDim];
				hh1[1] = new double[vDim];
				hh1[2] = new double[uDim];
				hh1[3] = new double[vDim];

				tt2[0] = new double[uDim];
				tt2[1] = new double[vDim];
				tt2[2] = new double[uDim];
				tt2[3] = new double[vDim];
				tt2[4] = new double[uDim];
				tt2[5] = new double[vDim];
				tt2[6] = new double[uDim];
				tt2[7] = new double[vDim];
				hh2[0] = new double[uDim];
				hh2[1] = new double[vDim];
				hh2[2] = new double[uDim];
				hh2[3] = new double[vDim];
				hh2[4] = new double[uDim];
				hh2[5] = new double[vDim];
				hh2[6] = new double[uDim];
				hh2[7] = new double[vDim];


				M[0] = new double* [uDim];
				M[1] = new double* [vDim];
				for (int i = 0; i < uDim; i++) {
					M[0][i] = new double[uDim];
				}
				for (int i = 0; i < vDim; i++) {
					M[1][i] = new double[vDim];
				}
				dd = new int* [nNode];
				for (int i = 0; i < nNode; i++) {
					dd[i] = new int[2];
				}
				

			}
			else {
				initialized = false;

			}
		}
	};
	const int ___ee[2]{ 0,1 };
	public class _memS {
	public:
		std::string mode;
		_memS_ref* _ref=0;
		double N[3];
		double bodyF;
		double BCF[2];
		double BCF2[2];
		double BCF6[2];
	public:
		double* __grad_z=0;
		double* __grad_phi = 0;
		int _nNode;
	private:
		//double* __mat = 0;
		double* _F3 = 0;
		double _F3_phi[2];
		double _F3_z[2];
		double _SLOPE_phi[2];
		double _SLOPE_z[2];
		double* __grad_C_z = 0;
		double* __grad_C_phi = 0;
		double* __grad_D_z = 0;
		double* __grad_D_phi = 0;
		double _K_phi[2];
		double* _K = 0;
		double lo[2]{ 0,0, };
		//double** M[2];
		//int** dd;
		int dim[2];
		int _uDim, _vDim;
		double gi[6];
		double Gi[6];
		double gij[4];
		double Gij[4];
		double bij[12];
		double Gammaijk[8];
		double _ss[4];
		double _Sij[4];
		double sc;
		//double* d0;
		//double* d1[2];
		//double* d2[4];
		//double* d2_star[4];
		//double star2[4];
		//int _star[4];
		//double* B[4];
		//double* tt0[2], * hh0[2], * tt1[4], * hh1[4], * tt2[8], * hh2[8];
		double* gradN[3]{ 0,0,0 };
		double* gradG = 0;
		const int ___ll[4]{ 0,3,6,9 };
	public:
		double x, y, z, Z, _Z, phi;
		double dv, _dv;
	public:
		inline void set_lo(double L1, double L2) {
			lo[0] = L1;
			lo[1] = L2;
		}
		inline void set_M(int i, int j, int k, double val) {
			_ref->M[i][j][k] = val;
		}
		inline void set_dd(int i, int j, int val) {
			_ref->dd[i][j] = val;
		}
		inline double component(double* G1, double* G2) {
			return G1[0] * _ss[0] * G2[0] + G1[1] * _ss[1] * G2[0] + G1[0] * _ss[1] * G2[1] + G1[1] * _ss[3] * G2[1];
		}
		inline void set_Sij(const double &A, const double &B, const double &C) {
			_ss[0] = C;
			_ss[1] = -B;
			_ss[2] = -B;
			_ss[3] = A;

			_Sij[0] = component(Gi, Gi);
			_Sij[3] = component(&(Gi[3]), &(Gi[3]));
			_Sij[1] = component(&(Gi[3]), Gi);
			_Sij[2] = component(Gi, &(Gi[3]));
		}

		inline double get_gi(const int &i, const int &s) {
			return gi[___ll[i] + s];
		}
		inline double get_Gi(const int &i, const int &s) {
			return Gi[___ll[i] + s];
		}
		inline double get_gij(const int &i, const int &j) {
			return gij[(i<<1) + j];
		}
		inline double get_Gij(const int &i, const int &j) {
			return Gij[(i<<1) + j];
		}
		inline double get_bij(const int &i, const int &j, const int &s) {
			return bij[___ll[((i << 1) + j)] + s];
		}
		inline double get_Gammaijk(const int &i, const int &j, const int &k) {
			return Gammaijk[(((i<<1) + j) <<1)+ k];
		}

		inline double get_tt0(const int &i, const int &s) {
			return _ref->tt0[i][s];
		}
		inline double get_hh0(const int &i, const int &s) {
			return _ref->hh0[i][s];
		}
		inline double get_tt1(const int &i, const int &j, const int &s) {
			return _ref->tt1[(i<<1) + j][s];
		}
		inline double get_hh1(const int &i, const int &j, const int &s) {
			return _ref->hh1[(i<<1) + j][s];
		}
		inline double get_tt2(const int &i, const int &j, const int &k, const int &s) {
			return _ref->tt2[(((i<<1) + j) <<1) + k][s];
		}
		inline double get_hh2(const int &i, const int &j, const int &k, const int &s) {
			return _ref->hh2[(((i<<1) + j) <<1) + k][s];
		}

	public:
		inline double _pow(const double &f, const int &k) {
			double val = 1;
			for (int i = 0;i < k;i++)
			{
				val *= f;
			}
			return val;
		}
		double __hh0(const int &j, const int &k) {
			double t = lo[j];
			return pow(t, (dim[j] - k - 1));
		}
		double __tt0(const int &j, const int &k) {
			double val = 0;
			double t = lo[j];
			for (int l = 0; l < dim[j]; l++)
			{
				val += get_hh0(j, l) * _ref->M[j][l][k];
			}
			return val;
		}
		
		double __hh1(const int &m, const int &j, const int &k) {
			double t = lo[j];
			if (j != m)
			{
				return pow(t, (dim[j] - k - 1));
			}
			else if (j == m && k < dim[j] - 1)
			{
				return (dim[j] - k - 1) * pow(t, (dim[j] - k - 2));
			}
			else {
				return 0;
			}
		}
		double __tt1(const int &m, const int &j, const int &k) {
			double val = 0;
			double t = lo[j];

			for (int l = 0; l < dim[j]; l++)
			{

				val += get_hh1(m, j, l) * _ref->M[j][l][k];
			}
			return val;
		}
		double __hh2(int n, int m, int j, int k) {
			double t = lo[j];
			if (j != m && j != n)
			{
				return pow(t, (dim[j] - k - 1));
			}
			if ((j != m && j == n) || (j == m && j != n))
			{
				if (k < dim[j] - 1)
				{
					return (dim[j] - k - 1) * pow(t, (dim[j] - k - 2));
				}
				else {
					return 0;
				}
			}
			if (j == m && j == n)
			{
				if (k < dim[j] - 2)
					return (dim[j] - k - 1) * (dim[j] - k - 2) * pow(t, (dim[j] - k - 3));
				else return  0;
			}
			return 0;
		}
		double __tt2(int n, int m, int j, int k) {
			double val = 0;
			for (int l = 0; l < dim[j]; l++)
			{
				val += get_hh2(n, m, j, l) * _ref->M[j][l][k];
			}
			return val;
		}
		double _shape(int k)
		{
			double shape = 1.0;
			for (int j = 0; j < 2; j++)
			{
				shape *= get_tt0(j, _ref->dd[k][j]);
			}
			//return shape[k];
			return shape;
		}
		double _C(int m, int k)
		{
			double C = 1.0;
			for (int j = 0; j < 2; j++)
			{
				C *= get_tt1(m, j, _ref->dd[k][j]);
			}
			//return C[m, k];
			return C;
		}
		double _D(int m, int n, int k)
		{
			double D = 1.0;
			for (int j = 0; j < 2; j++)
			{
				D *= get_tt2(m, n, j, _ref->dd[k][j]);
			}
			return D;
		}
		double _B(int i, int j, int u, int v) {
			//対称化
			double val = _ref->d1[i][u] * _ref->d1[j][v];
			val += _ref->d1[j][u] * _ref->d1[i][v];
			return val;
		}
	public:
		_memS(string ultimate) {
			this->mode = ultimate;
			if (mode=="U") {
				__grad_z = 0;
				__grad_phi = 0;
				__grad_C_z = 0;
				__grad_C_phi = 0;
				__grad_D_z = 0;
				__grad_D_phi = 0;
				//__mat = 0;
				_K = 0;
				//__mat2 = 0;
			}
			_F3 = 0;
			_nNode = 0;
			_uDim = -1;
			_vDim = -1;
			//M[0] = 0;
			//M[1] = 0;
			//dd = 0;
			gradN[0] = 0;
			gradN[1] = 0;
			gradN[2] = 0;
			gradG = 0;
			/*d0 = 0;
			d1[0] = 0;
			d1[1] = 0;
			d2[0] = 0;
			d2[1] = 0;
			d2[2] = 0;
			d2[3] = 0;
			d2_star[0] = 0;
			d2_star[1] = 0;
			d2_star[2] = 0;
			d2_star[3] = 0;

			star2[0] = 1;
			star2[1] = -1;
			star2[2] = -1;
			star2[3] = 1;

			_star[0] = 3;
			_star[1] = 1;
			_star[2] = 2;
			_star[3] = 0;

			B[0] = 0;
			B[1] = 0;
			B[2] = 0;
			B[3] = 0;
			tt0[0] = 0;
			tt0[1] = 0;
			hh0[0] = 0;
			hh0[1] = 0;
			tt1[0] = 0;
			tt1[1] = 0;
			tt1[2] = 0;
			tt1[3] = 0;
			hh1[0] = 0;
			hh1[1] = 0;
			hh1[2] = 0;
			hh1[3] = 0;
			tt2[0] = 0;
			tt2[1] = 0;
			tt2[2] = 0;
			tt2[3] = 0;
			tt2[4] = 0;
			tt2[5] = 0;
			tt2[6] = 0;
			tt2[7] = 0;
			hh2[0] = 0;
			hh2[1] = 0;
			hh2[2] = 0;
			hh2[3] = 0;
			hh2[4] = 0;
			hh2[5] = 0;
			hh2[6] = 0;
			hh2[7] = 0;*/

		}
		~_memS() {
			del();
		}
		void update2() {
			if (!_ref->initialized)
			{
				_ref->initialized = true;

				for (auto const& j:___ee)
				{
					for (int k = 0; k < dim[j]; k++)
					{
						_ref->hh0[j][k] = __hh0(j, k);
					}
				}
				for (auto const& j : ___ee)
				{
					for (int k = 0; k < dim[j]; k++) {
						_ref->tt0[j][k] = __tt0(j, k);
					}
				}
				for (auto const& m : ___ee)
				{
					for (int j = 0; j < 2; j++)
					{
						for (int k = 0; k < dim[j]; k++) {
							_ref->hh1[m * 2 + j][k] = __hh1(m, j, k);
						}
					}
				}
				double** ptr;
				ptr = _ref->tt1;
				for (auto const& m : ___ee)
				{
					for (auto const& j:___ee)
					{
						for (int k = 0; k < dim[j]; k++) {
							*(*ptr + k) = __tt1(m, j, k);
							//tt1[m * 2 + j][k] = __tt1(m, j, k);
						}
						ptr++;
					}
				}
				ptr = _ref->hh2;
				for (auto const& n:___ee)
				{
					for (auto const& m : ___ee)
					{
						for (auto const& j : ___ee)
						{
							for (int k = 0; k < dim[j]; k++) {
								*(*ptr + k) = __hh2(n, m, j, k);
								//hh2[(n * 2 + m) * 2 + j][k] = __hh2(n, m, j, k);
							}
							ptr++;
						}
					}
				}
				ptr = _ref->tt2;
				for (auto const& n : ___ee)
				{
					for (auto const& m : ___ee)
					{
						for (auto const& j : ___ee)
						{
							for (int k = 0; k < dim[j]; k++) {
								*(*ptr + k) = __tt2(n, m, j, k);
								//tt2[(n * 2 + m) * 2 + j][k] = __tt2(n, m, j, k);
							}
							ptr++;
						}
					}
				}
				double* ptr2;
				ptr2 = _ref->d0;
				for (int j = 0; j < _nNode; j++) {
					*ptr2 = _shape(j);
					ptr2++;
				}
				ptr = _ref->d1;
				for (auto const& i : ___ee) {
					ptr2 = *ptr;
					for (int j = 0; j < _nNode; j++) {
						//d1[i][j] = _C(i, j);
						*ptr2 = _C(i, j);
						ptr2++;
					}
					ptr++;
				}
				ptr = _ref->d2;
				for (auto const& i : ___ee) {
					for (auto const& ii : ___ee) {
						ptr2 = *ptr;
						for (int j = 0; j < _nNode; j++) {
							//d2[i * 2 + ii][j] = _D(i, ii, j);
							*ptr2 = _D(i, ii, j);
							ptr2++;
						}
						ptr++;
					}
				}

				ptr = _ref->B;
				for (auto const& i : ___ee) {
					for (auto const& ii : ___ee) {
						ptr2 = *ptr;
						for (int j = 0; j < _nNode; j++) {
							for (int jj = 0; jj < _nNode; jj++) {
								//B[i * 2 + ii][j * _nNode + jj] = _B(i, ii, j, jj);
								*ptr2 = _B(i, ii, j, jj);
								ptr2++;
							}
						}
						ptr++;
					}
				}
			}
			double X = 0, Y = 0, ZZ = 0;
			double *ptr2 = _ref->d0;
			double* ptr3;
			ptr3 = _ref->node;
			for (int i = 0; i < _nNode; i++)
			{
				X += *ptr2 * *ptr3;
				ptr3++;
				Y += *ptr2 * *ptr3;
				ptr3++;
				ZZ += *ptr2 * *ptr3;
				ptr3++;
				ptr2++;
			}
			x = X;
			y = Y;
			z = ZZ;

			ZZ = 0;
			ptr2 = _ref->d0;
			ptr3 = _ref->buf_z;
			for (int i = 0; i < _nNode; i++)
			{
				ZZ += *ptr2 * *ptr3;
				ptr3++;
				ptr2++;
			}
			this->Z = ZZ;
			double* pptr1 = &_ref->buf_phi[0];
			double* pptr2 = &_ref->d0[0];
			double val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			this->phi = val;
			//covariant base vectors

			for (auto const& j : ___ee) {
				double fx = 0, fy = 0, fz = 0;
				ptr2 = _ref->d1[j];
				ptr3 = _ref->node;
				for (int i = 0; i < _nNode; i++)
				{
					fx += *ptr2 * *ptr3;
					ptr3++;
					fy += *ptr2 * *ptr3;
					ptr3++;
					fz += *ptr2 * *ptr3;
					ptr3++;
					ptr2++;
				}
				gi[j * 3 + 0] = fx;
				gi[j * 3 + 1] = fy;
				gi[j * 3 + 2] = fz;
			}
			gij[0] = gi[0] * gi[0] + gi[1] * gi[1];
			gij[1] = gi[0] * gi[3] + gi[1] * gi[4];
			gij[2] = gij[1];
			gij[3] = gi[3] * gi[3] + gi[4] * gi[4];


			_dv = sqrt(_det2(gij));
			//sc = 1 / _dv / _dv;
			sc = 1.0 / _det2(gij);
			gij[0] = gi[0] * gi[0] + gi[1] * gi[1] + gi[2] * gi[2];
			gij[1] = gi[0] * gi[3] + gi[1] * gi[4] + gi[2] * gi[5];
			gij[2] = gij[1];
			gij[3] = gi[3] * gi[3] + gi[4] * gi[4] + gi[5] * gi[5];

			_inv2(gij, Gij);
			dv = sqrt(_det2(gij));

			//contravatiant base vectors
			double Fx = 0, Fy = 0, Fz = 0;
			Fx += get_gi(0, 0) * Gij[0] + get_gi(1, 0) * Gij[1];
			Fy += get_gi(0, 1) * Gij[0] + get_gi(1, 1) * Gij[1];
			Fz += get_gi(0, 2) * Gij[0] + get_gi(1, 2) * Gij[1];
			Gi[0] = Fx;
			Gi[1] = Fy;
			Gi[2] = Fz;
			Fx = 0; Fy = 0; Fz = 0;
			Fx += get_gi(0, 0) * Gij[2] + get_gi(1, 0) * Gij[3];
			Fy += get_gi(0, 1) * Gij[2] + get_gi(1, 1) * Gij[3];
			Fz += get_gi(0, 2) * Gij[2] + get_gi(1, 2) * Gij[3];
			Gi[3] = Fx;
			Gi[4] = Fy;
			Gi[5] = Fz;

			double gx = get_gi(0, 0);
			double gy = get_gi(0, 1);
			double gz = get_gi(0, 2);
			double hx = get_gi(1, 0);
			double hy = get_gi(1, 1);
			double hz = get_gi(1, 2);
			double Nx = gy * hz - gz * hy;
			double Ny = gz * hx - gx * hz;
			double Nz = gx * hy - gy * hx;
			double  norm = sqrt(Nx * Nx + Ny * Ny + Nz * Nz);
			Nx = Nx / norm;
			Ny = Ny / norm;
			Nz = Nz / norm;
			N[0] = Nx;
			N[1] = Ny;
			N[2] = Nz;

			double **ptr = _ref->d2;
			for (auto const& j : ___ee) {
				for (auto const& k : ___ee) {
					double fx = 0, fy = 0, fz = 0;
					int e = j * 2 + k;
					ptr2 = *ptr;
					ptr3 = _ref->node;
					for (int i = 0; i < _nNode; i++)
					{

						fx += *ptr2 * *ptr3;
						ptr3++;
						fy += *ptr2 * *ptr3;
						ptr3++;
						fz += *ptr2 * *ptr3;
						ptr3++;
						ptr2++;
					}
					ptr++;
					e *= 3;
					bij[e + 0] = fx;
					bij[e + 1] = fy;
					bij[e + 2] = fz;
				}
			}
			int ccc = 0;
			static const int fff[8]{ 0,2,1,3,4,6,5,7 };
			//static const int ggg[8]{ 0,1,4,5,2,3,6,7 };
			//static const int hhh[8]{ 0,4,2,6,1,5,3,7 };
			int eee = 0;
			for (auto const& i : ___ee) {
				for (auto const& j : ___ee) {
					for (auto const& k : ___ee) {
						double val = 0;
						
						val += bij[eee + 0] * get_Gi(k, 0);
						val += bij[eee + 1] * get_Gi(k, 1);
						val += bij[eee + 2] * get_Gi(k, 2);
						Gammaijk[ccc] = val;
						ccc++;
					}
					eee += 3;
				}
			}


			/*ccc = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						if (i <= j) {
							if (j <= k)
							{
								if (i <= k)
								{
								}
								else {
									//Gammaijk[ccc] = -Gammaijk[hhh[ccc]];
								}
							}
							else {
								Gammaijk[ccc] = -Gammaijk[ggg[ccc]];   //swap 1th and 2th bits
							}
						}
						//else {
							//Gammaijk[ccc] = Gammaijk[fff[ccc]];   //swap 1th and 2th bits
						//}
						ccc++;
					}
				}
			}*/

			if (mode == "U") {

				double* ptr4;
				/*for (auto const& j : ___ee) {
					double fx = 0, fy = 0, fz = 0;
					ptr2 = _ref->d1[j];
					ptr3 = _ref->node;
					ptr4 = _ref->buf_z;
					for (int i = 0; i < _nNode; i++)
					{
						fx += *ptr2 * *ptr3;
						ptr3++;
						fy += *ptr2 * *ptr3;
						ptr3++;
						fz += *ptr2 * *ptr4;
						ptr3++;
						ptr2++;
						ptr4++;
					}
					gi[j * 3 + 0] = fx;
					gi[j * 3 + 1] = fy;
					gi[j * 3 + 2] = fz;
				}
				gij[0] = gi[0] * gi[0] + gi[1] * gi[1]+ gi[2] * gi[2];
				gij[1] = gi[0] * gi[3] + gi[1] * gi[4]+ gi[2] * gi[5];
				gij[2] = gij[1];
				gij[3] = gi[3] * gi[3] + gi[4] * gi[4]+ gi[5] * gi[5];


				dv = sqrt(_det2(gij));*/


				double val = 0;
				double val2 = 0;
				double val3 = 0;
				double* pptr = 0;
				double* pptr1 = &_ref->buf_z[0];
				double* pptr2 = &_ref->d0[0];
				/*for (int i = 0; i < _nNode; i++)
				{
					//val += _ref->d0[i] * _ref->buf_z[i];
					val += (*pptr1) * (*pptr2);
					pptr1++;
					pptr2++;
				}
				this->Z = val;
				val = 0;
				pptr1 = &_ref->buf_phi[0];
				pptr2 = &_ref->d0[0];
				for (int i = 0; i < _nNode; i++)
				{
					//val += _ref->d0[i] * _ref->buf_phi[i];
					val += (*pptr1) * (*pptr2);
					pptr1++;
					pptr2++;
				}
				this->phi = val;*/
				if (_ref->__mat == 0)
				{
					_ref->__mat = new double[_nNode * _nNode];
					pptr = &_ref->__mat[0];
					static const int sss[4]{ 1,2,0,3 };
					for (int i = 0; i < _nNode; i++)
					{
						for (int j = 0; j < _nNode; j++) {
							double val = 0;
							for (auto const& k : sss) {
								if (k == 2)
								{
									val *= 2.0;
								}
								else {
									//double tmp = ;
										val += _ref->star2[k] * (_ref->d2[_ref->_star[k]][i] - Gammaijk[((_ref->_star[k]) << 1) + 0] * _ref->d1[0][i] - Gammaijk[((_ref->_star[k]) << 1) + 1] * _ref->d1[1][i]) *
										(_ref->d2[k][j] - Gammaijk[((k) << 1) + 0] * _ref->d1[0][j] - Gammaijk[((k) << 1) + 1] * _ref->d1[1][j]);
								}
							}
							*pptr = val;
							pptr++;
						}
					}
				}
	
				pptr = &_ref->__mat[0];
				pptr1 = &_ref->buf_phi[0];
				static const int sss[4]{ 1,2,0,3 };
				for (int i = 0; i < _nNode; i++)
				{
					pptr2 = &_ref->buf_z[0];
					for (int j = 0; j < _nNode; j++) {
						val3 += *pptr * (*pptr2) */*_ref->buf_z[j] **/ (*pptr1);//_ref->buf_phi[i];
						pptr++;
						pptr2++;
					}
					pptr1++;
				}
				bodyF = val3;
				pptr = &_ref->__mat[0];
				pptr2 = &__grad_z[0];
				for (int i = 0; i < _nNode; i++) {
					val2 = 0;
					pptr1 = &_ref->buf_z[0];
					for (int j = 0; j < _nNode; j++) {
						val2 += (*pptr)/*__mat[i * _nNode + j]*/ * (*pptr1)/*_ref->buf_z[j]*/;
						pptr1++;
						pptr++;
					}
					//__grad_z[i] = val2*sc;
					//eigen_assert(val2<1000 && val2>-1000);
					*pptr2 = val2 * sc;
					pptr2++;
				}
				pptr = _ref->__mat;
				pptr1 = &_ref->buf_phi[0];
				pptr2 = &__grad_phi[0];
				for (int j = 0; j < _nNode; j++) {
					*pptr2 = 0;
					pptr2++;
				}
				for (int i = 0; i < _nNode; i++) {
					//val2 = 0;
					pptr2 = &__grad_phi[0];
					for (int j = 0; j < _nNode; j++) {
						val2 = (*pptr)/*__mat[i * _nNode + j]*/ * (*pptr1)/* _ref->buf_phi[i]*/;
						pptr++;
						*pptr2 += val2 * sc;
						pptr2++;
					}
					//__grad_phi[j] = val2 * sc;
					pptr1++;
				}
			}
			if ( mode == "SLOPE") {
				for (auto const& l : ___ee) {
					_SLOPE_phi[l] = __SLOPE_phi(l);
				}
				for (auto const& l : ___ee) {
					_SLOPE_z[l] = __SLOPE_z(l);
				}
				for (int i = 0; i < _nNode; i++) {
					for (auto const& l : ___ee) {
						__grad_C_z[l * _nNode + i] = this->__F4_z(l, i);
						__grad_C_phi[l * _nNode + i] = this->__F4_phi(l, i);
					}
				}
				for (int i = 0; i < _nNode; i++) {
					for (auto const& l : ___ee) {
						__grad_D_z[l * _nNode + i] = this->__F5_z(l, i);
						__grad_D_phi[l * _nNode + i] = this->__F5_phi(l, i);
					}				
				}
				for (int i = 0; i < _nNode; i++) {
					for (auto const& l : ___ee) {
						_K[l * _nNode + i] = __K(l, i);
					}
				}
				for (auto const& l : ___ee)
				{
					_K_phi[l] = __K_phi(l);
				}


				BCF[0] = 0;
				BCF[1] = 0;
				BCF2[0] = 0;
				BCF2[1] = 0;
				BCF6[0] = 0;
				BCF6[1] = 0;
				for (int j = 0; j < _nNode; j++) {
					BCF[0] += __grad_C_z[j] * _ref->buf_phi[j];
					BCF[1] += __grad_C_z[_nNode + j] * _ref->buf_phi[j];
					BCF2[0] += __grad_D_phi[j] * _ref->buf_z[j];
					BCF2[1] += __grad_D_phi[_nNode + j] * _ref->buf_z[j];
				}
				for (int j = 0; j < _nNode; j++) {
					BCF6[0] += this->F6_z(0, j) * _ref->buf_phi[j];
					BCF6[1] += this->F6_z(1, j) * _ref->buf_phi[j];
				}
			}

			set_Sij(1, 0, 1);
			for (int i = 0; i < _nNode; i++) {
				for (auto const& l : ___ee) {
					_F3[l * _nNode + i] = __F3(l, i);
				}
			}
			for (auto const& l : ___ee) {
				_F3_phi[l] = __F3_phi(l);
			}
			for (auto const& l : ___ee) {
				_F3_z[l] = __F3_z(l);
			}


		}
		void del() {
			if (mode=="U" ) {
				if (__grad_z != 0)delete[] __grad_z;
				if (__grad_phi != 0)delete[] __grad_phi;
				//if (__mat != 0)delete[] __mat;
				if (__grad_C_z != 0)delete[] __grad_C_z;
				if (__grad_C_phi != 0)delete[] __grad_C_phi;
				if (__grad_D_z != 0)delete[] __grad_D_z;
				if (__grad_D_phi != 0)delete[] __grad_D_phi;
				if (_K != 0)delete[] _K;
				__grad_z = 0;
				__grad_phi = 0;
				//__mat = 0;
				__grad_C_z = 0;
				__grad_C_phi = 0;
				__grad_D_z = 0;
				__grad_D_phi = 0;
				_K = 0;
				//if (__mat2 != 0)delete[] __mat2;
			}
			if (_F3 != 0)delete[] _F3;
			_F3 = 0;
			/*
			if (_ref->M[0] != 0) {
				for (int i = 0; i < _uDim; i++) {
					delete[] M[0][i];
				}
				for (int i = 0; i < _vDim; i++) {
					delete[] M[1][i];
				}
				delete[] M[0];
				delete[] M[1];
				M[0] = 0;
				M[1] = 0;
			}

			if (dd != 0)
			{
				for (int i = 0; i < _nNode; i++) {
					delete[] dd[i];
				}
				delete[] dd;
				dd = 0;
			}
			if (d0 != 0) {
				delete[] d0;
				d0 = 0;
			}
			if (d1[0] != 0)
			{
				delete[] d1[0];
				delete[] d1[1];
				d1[0] = 0;
				d1[1] = 0;
			}

			if (d2[0] != 0)
			{
				delete[] d2[0];
				delete[] d2[1];
				delete[] d2[2];
				delete[] d2[3];
				d2[0] = 0;
				d2[1] = 0;
				d2[2] = 0;
				d2[3] = 0;
				delete[] d2_star[0];
				delete[] d2_star[1];
				delete[] d2_star[2];
				delete[] d2_star[3];
				d2_star[0] = 0;
				d2_star[1] = 0;
				d2_star[2] = 0;
				d2_star[3] = 0;
			}

			if (B[0] != 0) {
				delete[] B[0];
				delete[] B[1];
				delete[] B[2];
				delete[] B[3];
				B[0] = 0;
				B[1] = 0;
				B[2] = 0;
				B[3] = 0;

			}*/
			if (gradN[0] != 0) {
				delete[] gradN[0];
				delete[] gradN[1];
				delete[] gradN[2];
				gradN[0] = 0;
				gradN[1] = 0;
				gradN[2] = 0;

			}
			if (gradG != 0)delete[] gradG;
			/*if (tt0[0] != 0) {
				delete[] tt0[0];
				tt0[0] = 0;
				delete[] tt0[1];
				tt0[1] = 0;
				delete[] hh0[0];
				delete[] hh0[1];
				hh0[0] = 0;
				hh0[1] = 0;
				delete[] tt1[0];
				delete[] tt1[1];
				delete[] tt1[2];
				delete[] tt1[3];
				tt1[0] = 0;
				tt1[1] = 0;
				tt1[2] = 0;
				tt1[3] = 0;
				delete[] hh1[0];
				delete[] hh1[1];
				delete[] hh1[2];
				delete[] hh1[3];
				hh1[0] = 0;
				hh1[1] = 0;
				hh1[2] = 0;
				hh1[3] = 0;
				delete[] tt2[0];
				delete[] tt2[1];
				delete[] tt2[2];
				delete[] tt2[3];
				delete[] tt2[4];
				delete[] tt2[5];
				delete[] tt2[6];
				delete[] tt2[7];
				tt2[0] = 0;
				tt2[1] = 0;
				tt2[2] = 0;
				tt2[3] = 0;
				tt2[4] = 0;
				tt2[5] = 0;
				tt2[6] = 0;
				tt2[7] = 0;

				delete[] hh2[0];
				delete[] hh2[1];
				delete[] hh2[2];
				delete[] hh2[3];
				delete[] hh2[4];
				delete[] hh2[5];
				delete[] hh2[6];
				delete[] hh2[7];
				hh2[0] = 0;
				hh2[1] = 0;
				hh2[2] = 0;
				hh2[3] = 0;
				hh2[4] = 0;
				hh2[5] = 0;
				hh2[6] = 0;
				hh2[7] = 0;
			}*/
		}

		void update(int nNode, int uDim, int vDim) {
			if (!_ref->initialized)
			{
				_ref->update(nNode, uDim, vDim);
			}
			if (_nNode != nNode || _uDim != uDim || _vDim != vDim) {

				del();


				_nNode = nNode;
				_uDim = uDim;
				_vDim = vDim;
				dim[0] = _uDim;
				dim[1] = _vDim;
				if (mode=="U"||mode=="SLOPE")
				{
					__grad_z = new double[nNode];
					__grad_phi = new double[nNode];
					//__mat = new double[nNode * nNode];
					//__mat2 = new double[nNode * nNode];
					__grad_C_phi = new double[nNode * 2];
					__grad_C_z = new double[nNode * 2];
					__grad_D_phi = new double[nNode * 2];
					__grad_D_z = new double[nNode * 2];
					_K = new double[nNode*2];
				}
				_F3 = new double[nNode * 2];
				if (mode == "SENSITIVITY" || mode == "SHELL")
				{
					gradN[0] = new double[3 * _nNode];
					gradN[1] = new double[3 * _nNode];
					gradN[2] = new double[3 * _nNode];
					gradG = new double[8 * _nNode];
				}
				/*d0 = new double[nNode];

				d1[0] = new double[nNode];
				d1[1] = new double[nNode];

				d2[0] = new double[nNode];
				d2[1] = new double[nNode];
				d2[2] = new double[nNode];
				d2[3] = new double[nNode];

				d2_star[0] = new double[nNode];
				d2_star[1] = new double[nNode];
				d2_star[2] = new double[nNode];
				d2_star[3] = new double[nNode];

				B[0] = new double[nNode * nNode];
				B[1] = new double[nNode * nNode];
				B[2] = new double[nNode * nNode];
				B[3] = new double[nNode * nNode];

				tt0[0] = new double[uDim];
				tt0[1] = new double[vDim];

				hh0[0] = new double[uDim];
				hh0[1] = new double[vDim];

				tt1[0] = new double[uDim];
				tt1[1] = new double[vDim];
				tt1[2] = new double[uDim];
				tt1[3] = new double[vDim];
				hh1[0] = new double[uDim];
				hh1[1] = new double[vDim];
				hh1[2] = new double[uDim];
				hh1[3] = new double[vDim];

				tt2[0] = new double[uDim];
				tt2[1] = new double[vDim];
				tt2[2] = new double[uDim];
				tt2[3] = new double[vDim];
				tt2[4] = new double[uDim];
				tt2[5] = new double[vDim];
				tt2[6] = new double[uDim];
				tt2[7] = new double[vDim];
				hh2[0] = new double[uDim];
				hh2[1] = new double[vDim];
				hh2[2] = new double[uDim];
				hh2[3] = new double[vDim];
				hh2[4] = new double[uDim];
				hh2[5] = new double[vDim];
				hh2[6] = new double[uDim];
				hh2[7] = new double[vDim];
				gradN[0] = new double[3 * _nNode];
				gradN[1] = new double[3 * _nNode];
				gradN[2] = new double[3 * _nNode];
				gradG = new double[8 * _nNode];



				M[0] = new double* [uDim];
				M[1] = new double* [vDim];
				for (int i = 0; i < uDim; i++) {
					M[0][i] = new double[uDim];
				}
				for (int i = 0; i < vDim; i++) {
					M[1][i] = new double[vDim];
				}
				dd = new int* [nNode];
				for (int i = 0; i < nNode; i++) {
					dd[i] = new int[2];
				}*/

			}
		}
	private:
		inline void _inv2(double* From, double* to)
		{
			double det = _det2(From);
			if (det == 0)det = 0.0000000000001;
			to[0] = From[3] / det;
			to[3] = From[0] / det;
			to[1] = -From[1] / det;
			to[2] = -From[2] / det;

		}
		inline double _det2(double* m)
		{
			return m[0] * m[3] - m[1] * m[2];
		}
	public:
		inline double __SLOPE_phi(int l) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_phi[i] * get_Gij(k, l);
				}
			}
			return val;
		}
		inline double __SLOPE_z(int l) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_z[i] * get_Gij(k, l);
				}
			}
			return val;

		}
		inline double SLOPE_phi(int l) {
			/*double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_phi[i] * get_Gij(k, l);
				}
			}
			return val;*/
			return _SLOPE_phi[l];// val;
		}
		inline double SLOPE_z(int l) {
			/*double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_z[i] * get_Gij(k, l);
				}
			}
			return val;*/
			return _SLOPE_z[l];// val;
		}
		inline double SLOPE_phi(double dcdtstar0, double dcdtstar1,double sc){

			return sc*(_SLOPE_phi[0] * dcdtstar0 + _SLOPE_phi[1] * dcdtstar1);// val;
		}
		inline double SLOPE_z(double dcdtstar0,double dcdtstar1,double sc) {

			return sc*(_SLOPE_z[0]*dcdtstar0+_SLOPE_z[1]*dcdtstar1);// val;
		}
		inline double SLOPE(int l,int i) {
			return _ref->d1[0][i]*get_Gij(0, l)+ _ref->d1[1][i] * get_Gij(1, l);

		}
		inline void SLOPE(double* ptr, double dcdt1, double dcdt2, double sc) {
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1= ((_ref->d1[0][i] * get_Gij(0, 0) + _ref->d1[1][i] * get_Gij(1, 0))*dcdt1+
					(_ref->d1[0][i] * get_Gij(0, 1) + _ref->d1[1][i] * get_Gij(1, 1)) * dcdt2)*sc;
				ptr1++;

			}

		}
		inline double Dc(double a, double b,int l) {
			double f1 = a * this->get_gi(0, 0) + b * this->get_gi(0, 1);
			double f2 = a * this->get_gi(1, 0) + b * this->get_gi(1, 1);
			return f1* get_Gij(0, l) + f2 * get_Gij(1, l);
		}
		inline double __K(int l, int I) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				val += _ref->d1[k][I] * this->get_Gij(k, l);
			}
			return val;
		}
		inline double K(int l,int I)
		{
			return _K[l * _nNode + I];// val;
		}
		//this is essentially The distance between two points in the force diagram
		double __K_phi(int l)
		{
			double val = 0;
			for (int I = 0; I < _nNode; I++) {
				double val2 = 0;
				for (auto k:___ee)
				{
					val2 += _ref->d1[k][I] * this->get_Gij(k, l);
				}
				val += val2 * _ref->buf_phi[I];
			}
			return val;
		}
		double K_phi(int l)
		{
			return _K_phi[l];// val;
		}
		//stress function
		double F(int i, int j) {
			int e = i * _nNode + j;
			return 0.5*(_Sij[0] * _ref->B[0][e] + 2 * _Sij[1] * _ref->B[1][e] + _Sij[3] * _ref->B[3][e]) * _ref->refDv;
		}
		//stress function L2
		double F2(int i, int j) {
			return
				//return  (_Sij[0] * (d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * d1[1][i]) + 2 * _Sij[1] * (d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * d1[1][i]) + _Sij[3] * (d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * d1[1][i]))*
				(_Sij[0] * (_ref->d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][i]) + 2 * _Sij[1] * (_ref->d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) + _Sij[3] * (_ref->d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][i]))
				* (_Sij[0] * (_ref->d2[0][j] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][j] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][j]) + 2 * _Sij[1] * (_ref->d2[1][j] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][j] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][j]) + _Sij[3] * (_ref->d2[3][j] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][j] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][j])) * _ref->refDv;// / _ref->refDv / _ref->refDv;

		}
		double __F3(int i, int I) {
			double val = 0;
			for (int ii = 0; ii < 2; ii++)
			{
				val += _Sij[(ii<<1) + i] * _ref->d1[ii][I];
			}
			return val;// *_ref->refDv;
		}
		double F3(int i, int I) {
			double val = 0;
			for (int ii = 0; ii < 2; ii++)
			{
				val += _Sij[(ii<<1) + i] * _ref->d1[ii][I];
			}
			return val;// _F3[i * _nNode + I];// val;// *_ref->refDv;
		}

		
		double __F3_z(int i) {
			double val = 0;
			for (int I = 0; I < _nNode; I++)
			{
				/*double val2 = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					val2 += _Sij[ii * 2 + i] * d1[ii][I];
				}*/
				val += /*val2 */ _F3[i * _nNode + I] * _ref->buf_z[I];
			}
			return val;// *_ref->refDv;
		}
		double F3_z(int i) {
			double val = 0;
			for (int I = 0; I < _nNode; I++)
			{
				double val2 = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					val2 += _Sij[(ii<<1) + i] * _ref->d1[ii][I];
				}
				val += val2 * _ref->buf_z[I];
			}
			return val;// _F3_z[i];// val;// *_ref->refDv;
		}
		double __F3_phi(int i) {
			double val = 0;
			for (int I = 0; I < _nNode; I++)
			{
				/*double val2 = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					val2 += _Sij[ii * 2 + i] * d1[ii][I];
				}*/
				val += _F3[i*_nNode+I]/* val2*/ * _ref->buf_phi[I];
			}
			return val;// *_ref->refDv;
		}
		double F3_phi(int i) {
			double val = 0;
			for (int I = 0; I < _nNode; I++)
			{
				double val2 = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					val2 += _Sij[(ii<<1) + i] * _ref->d1[ii][I];
				}
				val += val2 * _ref->buf_phi[I];
			}
			return val;// _F3_phi[i];// val;// *_ref->refDv;
		}
		double F4(int i, int I,int J) {
			double val = 0;
			double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			for (int ii = 0; ii < 2; ii++)
			{
				int ij = ii * 2 + i;
				val += _ref->star2[ij]*(_ref->d2[_ref->_star[ij]][I]- Gammaijk[(_ref->_star[ij])*2+0]* _ref->d1[0][I] - Gammaijk[(_ref->_star[ij]) * 2 + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
			}
			return sc*val;
		}
		double F6(int k, int I, int J) {
			double val = 0;
			double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			for (int i = 0; i < 2; i++) {
				for (int ii = 0; ii < 2; ii++)
				{
					int ij = ii * 2 + i;
					for (int l = 0; l < 2; l++)
					{
						val += (_ref->d2[ij][I] - Gammaijk[((ij)<<1) + 0] * _ref->d1[0][I] - Gammaijk[((ij)<<1)+ 1] * _ref->d1[1][I]) * this->get_Gij(ii, k) * this->get_Gij(i, l) * _ref->d1[l][J];
					}
				}
			}
			return sc * val;
		}
		double F6_z(int k, int I) {
			double val = 0;
			double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			for (int J = 0; J < _nNode; J++)
			{
				double val2 = 0;
				for (int i = 0; i < 2; i++) {
					for (int ii = 0; ii < 2; ii++)
					{
						int ij = ii * 2 + i;
						for (int l = 0; l < 2; l++)
						{
							val2 += (_ref->d2[ij][I] - Gammaijk[((ij) <<1) + 0] * _ref->d1[0][I] - Gammaijk[((ij)<<1) + 1] * _ref->d1[1][I]) * this->get_Gij(ii, k) * this->get_Gij(i, l) * _ref->d1[l][J];
						}
					}
				}
				val += val2 * _ref->buf_z[J];
			}
			return sc * val;
		}
		double F6_phi(int k, int J) {
			double val = 0;
			double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			for (int I = 0; I < _nNode; I++)
			{
				double val2 = 0;
				for (int i = 0; i < 2; i++) {
					for (int ii = 0; ii < 2; ii++)
					{
						int ij = ii * 2 + i;
						for (int l = 0; l < 2; l++)
						{
							val += (_ref->d2[ij][I] - Gammaijk[((ij) <<1)+ 0] * _ref->d1[0][I] - Gammaijk[((ij) <<1) + 1] * _ref->d1[1][I]) * this->get_Gij(ii, k) * this->get_Gij(i, l) * _ref->d1[l][J];
						}
					}
				}
				val += val2 * _ref->buf_phi[I];
			}
			return sc * val;
		}
	private:

		double __F4_z(int i, int I) {
			double val = 0;
			double val2 = 0;
			
			for (int J = 0; J < _nNode; J++) {
				val = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					int ij = ii * 2 + i;
					val += _ref->star2[ij] * (_ref->d2[_ref->_star[ij]][I] - Gammaijk[((_ref->_star[ij]) <<1) + 0] * _ref->d1[0][I] - Gammaijk[((_ref->_star[ij]) <<1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
				}
				val2 += val * _ref->buf_z[J];
			}
			return val2;
		}
		double __F4_phi(int i, int J) {
			double val = 0;
			double val2 = 0;

			for (int I = 0; I < _nNode; I++) {
				val = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					int ij = (ii<<1) + i;
					val += _ref->star2[ij] * (_ref->d2[_ref->_star[ij]][I] - Gammaijk[((_ref->_star[ij]) <<1) + 0] * _ref->d1[0][I] - Gammaijk[((_ref->_star[ij]) <<1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
				}
				val2 += val * _ref->buf_phi[I];
			}
			return val2;
		}
		double __F5_z(int i, int J) {
			double val = 0;
			double val2 = 0;

			for (int I = 0; I < _nNode; I++) {
				val = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					int ij = (ii<<1) + i;
					val += _ref->star2[ij] * (_ref->d2[_ref->_star[ij]][I] - Gammaijk[((_ref->_star[ij]) <<1)+ 0] * _ref->d1[0][I] - Gammaijk[((_ref->_star[ij]) <<1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
				}
				val2 += val * _ref->buf_z[I];
			}
			return val2;
		}
		double __F5_phi(int i, int I) {
			double val = 0;
			double val2 = 0;

			for (int J = 0; J < _nNode; J++) {
				val = 0;
				for (int ii = 0; ii < 2; ii++)
				{
					int ij = (ii<<1) + i;
					val += _ref->star2[ij] * (_ref->d2[_ref->_star[ij]][I] - Gammaijk[((_ref->_star[ij]) <<1) + 0] * _ref->d1[0][I] - Gammaijk[((_ref->_star[ij]) <<1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
				}
				val2 += val * _ref->buf_phi[J];
			}
			return val2;
		}
	public:
		double F4_z(int i, int I) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;

			return sc * __grad_C_z[i * _nNode + I];
		}
		double F4_phi(int i, int J) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;

			return sc *  __grad_C_phi[i * _nNode + J];
		}
		void F4_z(double *ptr,double dcdtstar0,double dcdtstar1,double sc2) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			double* ptr1 = ptr;
			double _sc1 = sc * sc2* dcdtstar0;
			double _sc2 = sc * sc2* dcdtstar1;

			double* ptr2 = __grad_C_z;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = _sc1* (*ptr2) +_sc2*(*(ptr2+_nNode));
				ptr1++;
				ptr2++;
			}

			//return sc * __grad_C_z[i * _nNode + I];
		}
		void F4_phi(double* ptr, double dcdtstar0, double dcdtstar1, double sc2) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			double* ptr1 = ptr;
			double _sc1 = sc * sc2*dcdtstar0;
			double _sc2 = sc * sc2*dcdtstar1;

			double* ptr2 = __grad_C_phi;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = _sc1*(*ptr2) +_sc2*(*(ptr2+_nNode)) ;
				ptr1++;
				ptr2++;
			}

			//return sc * __grad_C_z[i * _nNode + I];
		}
		double F5_z(int i, int I) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;

			return sc *  __grad_D_z[i * _nNode + I];
		}
		double F5_phi(int i, int J) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;

			return sc * __grad_D_phi[i * _nNode + J];
		}
		double error() {
			double bij00 = this->get_bij(0, 0, 2) - this->get_Gammaijk(0, 0, 0) * this->get_gi(0, 2) - this->get_Gammaijk(0, 0, 1) * this->get_gi(1, 2);
			double bij01 = this->get_bij(0, 1, 2) - this->get_Gammaijk(0, 1, 0) * this->get_gi(0, 2) - this->get_Gammaijk(0, 1, 1) * this->get_gi(1, 2);
			double bij11 = this->get_bij(1, 1, 2) - this->get_Gammaijk(1, 1, 0) * this->get_gi(0, 2) - this->get_Gammaijk(1, 1, 1) * this->get_gi(1, 2);

			return _ref->_Sij[0] * bij00 +
				2 * _ref->_Sij[1] * bij01
				+ _ref->_Sij[3] * bij11;

		}
		//body force term reference
		double G(int i) {
			return _ref->d0[i] * _ref->refDv;
		}
		//body force term projected
		double G4(int i) {
			return _ref->d0[i] * _dv;
		}
		//body force term accurate element area
		double G2(int i) {
			return _ref->d0[i] * dv;
		}
		double _d0(int i) {
			return _ref->d0[i];
		}
		double _d2(int l,int i) {
			return _ref->d2[l][i] - _ref->_Gammaijk[l * 2 + 0] * _ref->d1[0][i] - _ref->_Gammaijk[l * 2 + 1] * _ref->d1[1][i];
		}
		double G3(int i) {
			//return d0[i] * _ref->refDv;
			return (_Sij[0] * (_ref->d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][i]) + 2 * _Sij[1] * (_ref->d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) + _Sij[3] * (_ref->d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][i])) * _ref->refDv;// / _ref->refDv / _ref->refDv;

		}
		//Hessian
		double _B(int i, int j) {
			double val = 0;
			for (int l = 0; l < 2; l++)
			{
				for (int m = 0; m < 2; m++)
				{

					int n = 2;

					val += 0.5 * _ref->d0[i] * _ref->d1[l][j] * this->get_gi(m, n) * this->dv * this->get_Gij(l, m);
					val += 0.5 * _ref->d0[i] * _ref->d1[m][j] * this->get_gi(l, n) * this->dv * this->get_Gij(l, m);

				}
			}
			return val;
		}
		//linear spring term
		double R(int i) {
			return -_ref->d0[i] * _ref->_z;
		}
		double R2(int i) {
			return -_ref->d0[i] * pow((_ref->_z - _ref->__z), 3) * _ref->_refDv;
		}
		double R3(int i) {
			return _ref->d0[i] * _ref->_refDv;
		}
		double MASS(int i, int j) {
			return _ref->d0[i] * _ref->d0[j] * _ref->_refDv;
		}
		double SMOOTH(int I, int J) {
			double val = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++)
						{
							val+=(_ref->d2[i * 2 + j][I]-this->get_Gammaijk(i,j,0)* _ref->d1[0][I]- this->get_Gammaijk(i, j, 1) * _ref->d1[1][I]) * this->get_Gij(i, k)* this->get_Gij(j, l)*(_ref->d2[k * 2 + l][J] - this->get_Gammaijk(k, l, 0) * _ref->d1[0][J] - this->get_Gammaijk(k, l, 1) * _ref->d1[1][J]);
						}
					}
				}
			}
			return val * _ref->_refDv;
		}
		double SMOOTH_z(int I)
		{
			double val = 0;
			for (int J = 0; J < _nNode; J++)
			{
				double val2 = 0;
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++)
					{
						for (int k = 0; k < 2; k++) {
							for (int l = 0; l < 2; l++)
							{
								val2 += (_ref->d2[i * 2 + j][I] - this->get_Gammaijk(i, j, 0) * _ref->d1[0][I] - this->get_Gammaijk(i, j, 1) * _ref->d1[1][I]) * this->get_Gij(i, k) * this->get_Gij(j, l) * (_ref->d2[k * 2 + l][J] - this->get_Gammaijk(k, l, 0) * _ref->d1[0][J] - this->get_Gammaijk(k, l, 1) * _ref->d1[1][J]);
							}
						}
					}
				}
				val += val2 * _ref->buf_z[J];
			}
			return val * _ref->_refDv;
		}
		double SMOOTH_phi(int I)
		{
			double val = 0;
			for (int J = 0; J < _nNode; J++)
			{
				double val2 = 0;
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++)
					{
						for (int k = 0; k < 2; k++) {
							for (int l = 0; l < 2; l++)
							{
								val2 += (_ref->d2[i * 2 + j][I] - this->get_Gammaijk(i, j, 0) * _ref->d1[0][I] - this->get_Gammaijk(i, j, 1) * _ref->d1[1][I]) * this->get_Gij(i, k) * this->get_Gij(j, l) * (_ref->d2[k * 2 + l][J] - this->get_Gammaijk(k, l, 0) * _ref->d1[0][J] - this->get_Gammaijk(k, l, 1) * _ref->d1[1][J]);
							}
						}
					}
				}
				val += val2 * _ref->buf_phi[J];
			}
			return val * _ref->_refDv;
		}
		double area() {
			return _ref->_refDv;
		}
		double A() {
			return this->dv;
		}
		double dA(int i) {
			double da = 0;
			for (int l = 0; l < 2; l++)
			{
				for (int m = 0; m < 2; m++)
				{
					da += 0.5 * get_Gij(l, m) * (get_gi(l, 2) * _ref->d1[m][i] + get_gi(m, 2) * _ref->d1[l][i]) * dv;
				}
			}
			return da;
		}

		void fM(double _la,double _mu,double*a,double*b,double*c,double*d)
		{
			//vector<double> ret;
			double S[4]; //covariant
			S[0] = 0;
			S[1] = 0;
			S[2] = 0;
			S[3] = 0;
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					double val = 0;
					for (int m = 0; m < 2; m++) {
						for (int n = 0; n < 2; n++) {


							double A = _la * _ref->get__Gij(l, k) * _ref->get__Gij(n, m) +2*_mu* _ref->get__Gij(l, n) * _ref->get__Gij(k, m);
							//A = 1.0;

							double D = 0;
							//double E = 0;
							for (int s = 0; s < 3; s++)
							{
								auto ff = get_gi(n, s);
								auto gg = _ref->get__gi(n, s);

								D += _ref->get__gi(m, s) * (get_gi(n, s)-_ref->get__gi(n, s));
								D += _ref->get__gi(n, s) * (get_gi(m, s)-_ref->get__gi(m, s));
							}
							//double D2 = /*get_gij(n, m) - */ _ref->get__gij(n, m);
							val += A * D;
							
							//S[(m << 1) + n] = D;

						}
					}
					S[(k << 1) + l] = val;// get_gij(k, l);// -_ref->get__gij(k, l);
				}
			}
			*a = S[0] * _ref->get__gi(0, 0) * _ref->get__gi(0, 0) + S[1] * _ref->get__gi(0, 0) * _ref->get__gi(1, 0) + S[2] * _ref->get__gi(1, 0) * _ref->get__gi(0, 0) + S[3] * _ref->get__gi(1, 0) * _ref->get__gi(1, 0);
			*b = S[0] * _ref->get__gi(0, 0) * _ref->get__gi(0, 1) + S[1] * _ref->get__gi(0, 0) * _ref->get__gi(1, 1) + S[2] * _ref->get__gi(1, 0) * _ref->get__gi(0, 1) + S[3] * _ref->get__gi(1, 0) * _ref->get__gi(1, 1);
			*c = S[0] * _ref->get__gi(0, 1) * _ref->get__gi(0, 0) + S[1] * _ref->get__gi(0, 1) * _ref->get__gi(1, 0) + S[2] * _ref->get__gi(1, 1) * _ref->get__gi(0, 0) + S[3] * _ref->get__gi(1, 1) * _ref->get__gi(1, 0);
			*d = S[0] * _ref->get__gi(0, 1) * _ref->get__gi(0, 1) + S[1] * _ref->get__gi(0, 1) * _ref->get__gi(1, 1) + S[2] * _ref->get__gi(1, 1) * _ref->get__gi(0, 1) + S[3] * _ref->get__gi(1, 1) * _ref->get__gi(1, 1);

			

			/*S[0] = get_gij(0, 0) - _ref->get__gij(0, 0);
			S[1] = get_gij(1, 0) - _ref->get__gij(1, 0);
			S[2] = get_gij(0, 1) - _ref->get__gij(0, 1);
			S[3] = get_gij(1, 1) - _ref->get__gij(1, 1);*/

			//return (S[0] * S[3] - S[1] * S[2]);// *(_ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1));
		}
	    double eM(double _la, double _mu) {
			double membrane = 0;
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					
					for (int m = 0; m < 2; m++) {
						for (int n = 0; n < 2; n++) {
							

							double A = _la * _ref->get__Gij(l, k) * _ref->get__Gij(n, m) + 2*_mu*_ref->get__Gij(l, n) * _ref->get__Gij(k, m);


							double D = 0;
							double E = 0;
							for (int s = 0; s < 3; s++)
							{
								D += _ref->get__gi(m, s) * (get_gi(n, s) - _ref->get__gi(n, s));
								D += _ref->get__gi(n, s) * (get_gi(m, s) - _ref->get__gi(m, s));
								E += _ref->get__gi(k, s) * (get_gi(l, s) - _ref->get__gi(l, s));
								E += _ref->get__gi(l, s) * (get_gi(k, s) - _ref->get__gi(k, s));
							}
							membrane += 0.25 * A * D * E * _ref->refDv;
							//Print((bij - _bij).ToString() + "," +  (tup.Gij[m, n] - tup._Gij[m, n]).ToString());
						}
					}
				}
			}
			return membrane;
		}

		double eK(double _la,double _mu) {
			double bending = 0;
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					double __bij2 = get_bij(k, l, 0) * N[0] + get_bij(k, l, 1) * N[1] + get_bij(k, l, 2) * N[2];
					double ___bij2 = _ref->get__bij(k, l, 0) * N[0] + _ref->get__bij(k, l, 1) * N[1] + _ref->get__bij(k, l, 2) * N[2];
					for (int m = 0; m < 2; m++) {
						for (int n = 0; n < 2; n++) {
							double __bij = get_bij(m, n, 0) * N[0] + get_bij(m, n, 1) * N[1] + get_bij(m, n, 2) * N[2];
							double ___bij = _ref->get__bij(m, n, 0) * N[0] + _ref->get__bij(m, n, 1) * N[1] + _ref->get__bij(m, n, 2) * N[2];

							double A = _la * _ref->get__Gij(l, k) * _ref->get__Gij(n, m) + 2*_mu*_ref->get__Gij(l, n) * _ref->get__Gij(k, m);

							double D = 0;
							double E = 0;
							for (int s = 0; s < 3; s++)
							{
								D += N[s] * ((get_bij(m, n, s) - _ref->get__bij(m, n, s)) - _ref->get__Gammaijk(m, n, 0) * (get_gi(0, s) - _ref->get__gi(0, s)) - _ref->get__Gammaijk(m, n, 1) * (get_gi(1, s) - _ref->get__gi(1, s)));
								E += N[s] * ((get_bij(k, l, s) - _ref->get__bij(k, l, s)) - _ref->get__Gammaijk(k, l, 0) * (get_gi(0, s) - _ref->get__gi(0, s)) - _ref->get__Gammaijk(k, l, 1) * (get_gi(1, s) - _ref->get__gi(1, s)));
								//E += 0.5 * N[s] * ((tup.bij[l, k, s] - tup._bij[l, l, s]) - tup._Gammaijk[l, k, 0] * (tup.gi[0, s] - tup._gi[0, s]) - tup._Gammaijk[l, m, 1] * (tup.gi[1, s] - tup._gi[1, s]));
							}
							bending += A * D * E * _ref->refDv;

						}
					}
				}
			}
			return bending;
		}
		void computeGrads() {
			for (int i = 0; i < _nNode; i++) {
				for (int s = 0; s < 3; s++) {
					double valx = 0;
					double valy = 0;
					double valz = 0;
					for (int k = 0; k < 2; k++) {
						valx += -_ref->d1[k][i] * (N[0]) * get_Gi(k, s);
						valy += -_ref->d1[k][i] * (N[1]) * get_Gi(k, s);
						valz += -_ref->d1[k][i] * (N[2]) * get_Gi(k, s);

					}
					gradN[s][i * 3 + 0] = valx;
					gradN[s][i * 3 + 1] = valy;
					gradN[s][i * 3 + 2] = valz;
				}
			}
			//gradient of Gammaijk
			for (int A = 0; A < _nNode; A++) {
				for (int i = 0; i < 2; i++) {
					for (int j = 0; j < 2; j++) {
						for (int k = 0; k < 2; k++) {
							double val = 0;
							val = _ref->d2[i*2+ j][A] * get_Gi(k, 2);
							for (int l = 0; l < 2; l++) {
								val += get_bij(i, j, 2) * get_Gij(k, l) * _ref->d1[l][A];
								for (int m = 0; m < 2; m++) {
									
									for (int nn = 0; nn < 2; nn++) {
										double AA = -get_Gij(m, k) * (_ref->d1[nn][A] * get_gi(m, 2) + _ref->d1[m][A] * get_gi(nn, 2)) * get_Gij(nn, l);
										for (int oo = 0; oo < 3; oo++) {
											val += get_bij(i, j, oo) * get_gi(l, oo) * AA;
										}
									}
								}
							}
							gradG[((i*2+j)*2+k)*_nNode+A] = val;
						}
					}
				}
			}
		}
		//gradient of the bending matrix
		double L(int alpha, double _la, double _mu)
		{

			double ddv = 0;
			for (int n = 0; n < 2; n++) {
				for (int m = 0; m < 2; m++) {
					ddv += 0.5 * get_Gij(n, m) * (get_gi(n, 2) * _ref->d1[m][alpha] + get_gi(m, 2) * _ref->d1[n][alpha]);
				}
			}


			double coeff = dv;// * tup.area * tup.K;
			double membrane = 0;
			double bending = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++)
						{


							double ik = 0;
							double jl = 0;
							double ij = 0;
							double kl = 0;
							for (int n = 0; n < 2; n++) {
								for (int m = 0; m < 2; m++) {
									ik += -get_Gij(i, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(k, m);
									jl += -get_Gij(j, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(l, m);
									ij += -get_Gij(i, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(j, m);
									kl += -get_Gij(k, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(l, m);
								}
							}
							double ijkl = _la * (ij * get_Gij(k, l) + kl * get_Gij(i, j)) + 2 * _mu * (ik * get_Gij(j, l) + get_Gij(i, k) * jl);

							double ff = get_Gij(k, l);

							double ijkl2 = _la * get_Gij(i, j) * get_Gij(k, l) + 2 * _mu * get_Gij(i, k) * get_Gij(j, l);

							for (int I = 0; I < _nNode; I++) {

								for (int J = 0; J < _nNode; J++) {

									//if(fx[I] ||fx[J])continue;
									double _E1 = _ref->d2[i * 2 + j][I] - get_Gammaijk(i, j, 0) * _ref->d1[0][I] - get_Gammaijk(i, j, 1) * _ref->d1[1][I];
									double _E2 = _ref->d2[k * 2 + l][J] - get_Gammaijk(k, l, 0) * _ref->d1[0][J] - get_Gammaijk(k, l, 1) * _ref->d1[1][J];


									for (int _i = 0; _i < 3; _i++) {
										for (int _j = 0; _j < 3; _j++) {
											double E1 = N[_i] * _E1;
											double E2 = N[_j] * _E2;
											double AA = _ref->def[I * 3 + _i] * _ref->def[J * 3 + _j] * E1 * E2 * coeff; //hh3 term
											double CC = ijkl + ijkl2 * ddv;
											bending += AA * CC;  //hh3 term
										}
									}
								}
							}
						}
					}
				}
			}

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++)
						{
							double ijkl = _la * get_Gij(i, j) * get_Gij(k, l) + 2 * _mu * get_Gij(i, k) * get_Gij(j, l);

							for (int I = 0; I < _nNode; I++) {

								for (int J = 0; J < _nNode; J++) {

									//if(fx[I] ||fx[J])continue;
									double _E1 = (_ref->d2[k * 2 + l][J] - get_Gammaijk(k, l, 0) * _ref->d1[0][J] - get_Gammaijk(k, l, 1) * _ref->d1[1][J]);
									double _E2 = (_ref->d2[i * 2 + j][I] - get_Gammaijk(i, j, 0) * _ref->d1[0][I] - get_Gammaijk(i, j, 1) * _ref->d1[1][I]);
									double A1 = _ref->d2[i * 2 + j][I];
									double A2 = _ref->d2[k * 2 + l][J];



									double _B11 = 0, _B21 = 0;
									double _B12 = 0, _B22 = 0;
									for (int kk = 0; kk < 2; kk++) {
										_B11 += -gradG[((i * 2 + j) * 2 + kk) * _nNode + alpha] * (_ref->d1[kk][I]);

										_B21 += -get_Gammaijk(i, j, kk) * _ref->d1[kk][I];
										_B12 += -gradG[((k * 2 + l) * 2 + kk) * _nNode + alpha] * (_ref->d1[kk][J]);
										_B22 += -get_Gammaijk(k, l, kk) * _ref->d1[kk][J];
									}

									for (int mm = 0; mm < 3; mm++) {
										for (int nn = 0; nn < 3; nn++) {

											double E1 = N[mm] * _E1;
											double E2 = N[nn] * _E2;
											double C1 = 0;
											double C2 = 0;

											C1 = ((A1 + _B21) * gradN[nn][alpha * 3 + 2] + _B11 * N[nn]);
											C2 = ((A2 + _B22) * gradN[mm][alpha * 3 + 2] + _B12 * N[mm]);

											bending += _ref->def[I * 3 + nn] * _ref->def[J * 3 + mm] * ijkl * (C1 * E1 + C2 * E2) * coeff;
										}
									}

								}
							}
						}
					}
				}
			}
			return bending;
		}

		//gradient of the membrane stiffness matrix
		double _M(int alpha,double _la, double _mu)
		{

			double ddv = 0;
			for (int n = 0; n < 2; n++) {
				for (int m = 0; m < 2; m++) {
					ddv += 0.5 * get_Gij(n, m) * (get_gi(n, 2) * _ref->d1[m][alpha] + get_gi(m, 2) * _ref->d1[n][alpha]);
				}
			}


			double coeff = dv;// * tup.area * tup.K;
			double membrane = 0;
			double bending = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++)
						{


							double ik = 0;
							double jl = 0;
							double ij = 0;
							double kl = 0;
							for (int n = 0; n < 2; n++) {
								for (int m = 0; m < 2; m++) {
									ik += -get_Gij(i, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(k, m);
									jl += -get_Gij(j, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(l, m);
									ij += -get_Gij(i, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(j, m);
									kl += -get_Gij(k, n) * (_ref->d1[n][alpha] * get_gi(m, 2) + _ref->d1[m][alpha] * get_gi(n, 2)) * get_Gij(l, m);
								}
							}
							double ijkl = _la * (ij * get_Gij(k, l) + kl * get_Gij(i, j)) + 2*_mu*(ik * get_Gij(j, l) + get_Gij(i, k) * jl);



							double ijkl2 = _la * get_Gij(i, j) * get_Gij(k, l) + 2*_mu * (get_Gij(i, k) * get_Gij(j, l));

							for (int I = 0; I < _nNode; I++) {

								for (int J = 0; J < _nNode; J++) {

									//if(fx[I] ||fx[J])continue;

									for (int _i = 0; _i < 3; _i++) {
										for (int _j = 0; _j < 3; _j++) {

											double FF = 0.5 * (_ref->d1[i][I] * get_gi(j, _i) + _ref->d1[j][I] * get_gi(i, _i));
											double GG = 0.5 * (_ref->d1[k][J] * get_gi(l, _j) + _ref->d1[l][J] * get_gi(k, _j));

											double BB = _ref->def[I * 3 + _i] * _ref->def[J * 3 + _j] * FF * GG * coeff;  //hh term
											double CC = ijkl + ijkl2 * ddv;

											membrane += BB * CC;  //hh term
										}
									}
								}
							}
						}
					}
				}
			}

			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++) {
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++)
						{
							double ijkl = _la *get_Gij(i, j)* get_Gij(k, l) + 2*_mu*get_Gij(i, k) * get_Gij(j, l);

							for (int I = 0; I < _nNode; I++) {

								for (int J = 0; J < _nNode; J++) {

									//if(fx[I] ||fx[J])continue;

									double _F1 = _ref->d1[i][alpha] * _ref->d1[j][I] + _ref->d1[j][alpha] * _ref->d1[i][I];//=zero when x or y
									double _F2 = _ref->d1[k][alpha] * _ref->d1[l][J] + _ref->d1[l][alpha] * _ref->d1[k][J];//=zero when x or y
									for (int mm = 0; mm < 3; mm++) {
										for (int nn = 0; nn < 3; nn++) {
											double F1 = 0.25 * _F1 * (get_gi(k, mm) * _ref->d1[l][J] + get_gi(l, mm) * _ref->d1[k][J]);
											double F2 = 0.25 * _F2 * (get_gi(i, nn) * _ref->d1[j][I] + get_gi(j, nn) * _ref->d1[i][I]);
											if (nn == 2)membrane += _ref->def[I * 3 + 2] * _ref->def[J * 3 + mm] * ijkl * (F1)*coeff;
											if (mm == 2)membrane += _ref->def[I * 3 + nn] * _ref->def[J * 3 + 2] * ijkl * (F2)*coeff;
										}
									}

								}
							}
						}
					}
				}
			}
		return membrane;
		}
		
		//bending term
		double K(int i, int k2, int j, int k,double _la, double _mu)
		{
			double _val3 = 0;
			double _val4 = 0;
			double _val5 = 0;
			double _val6 = 0;
			for (int l = 0; l < 2; l++)
			{
				for (int m = 0; m < 2; m++)
				{
					for (int g = 0; g < 2; g++)
					{
						for (int h = 0; h < 2; h++)
						{
							double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));
							double D = (_ref->d2[g * 2 + h][j] *N[k]+ this->get_bij(g, h, k) * gradN[k][j]); // -Gammaijk[(g * 2 + h) * 2 + 0] * _ref->d1[0][j] - Gammaijk[(g * 2 + h) * 2 + 1] * _ref->d1[1][j]);
							double E = (_ref->d2[l * 2 + m][i] * N[2] + this->get_bij(l, m, k2) * gradN[k2][i]);// - Gammaijk[(l * 2 + m) * 2 + 0] * _ref->d1[0][i] - Gammaijk[(l * 2 + m) * 2 + 1] * _ref->d1[1][i]);
							_val3 += A *  (D) * (E);
							//_val4 += A * this->gradN[k][j] * get_bij(g, h, k) * this->gradN[k2][i] * get_bij(l, m, k2);
							//_val5 += A * this->gradN[k][j] * get_bij(g, h, k) * N[k2] * E;
							//_val6 += A * this->gradN[k2][i] * get_bij(l, m, k2) * N[k] * D;
						}
					}
				}
			}
			return (_val4 + _val3 + _val5 + _val6) * this->dv;
		}
		double dK(int i, int k, double _la, double _mu)
		{
			double _val3 = 0;
			for (int l = 0; l < 2; l++)
			{
				for (int m = 0; m < 2; m++)
				{
					for (int g = 0; g < 2; g++)
					{
						for (int h = 0; h < 2; h++)
						{
							double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));
							double D = _ref->get__bij(g,h,k);
							double E = (_ref->d2[l * 2 + m][i] - Gammaijk[(l * 2 + m) * 2 + 0] * _ref->d1[0][i] - Gammaijk[(l * 2 + m) * 2 + 1] * _ref->d1[1][i]);
							_val3 += A * N[k]* (D) * (E);
						}
					}
				}
			}
			return (_val3) * this->dv;
		}
		void K(_mySparse* M,_mySparse* mat, int* _index, double _la, double _mu, double __sc)
		{
			const static int kk[3]{ 0,1,2 };
			const static int ll[2]{ 0,1 };
			//double _sc = sc * _ref->refDv;
			//Eigen::SparseMatrix<double, Eigen::ColMajor> _mat(mat->_mat[0].rows(), mat->_mat[0].cols());
			std::vector<Eigen::Triplet<double>> dat;//
			//dat.resize(_nNode * 3 * _nNode * 3);
			//dat.clear();
			dat.reserve(_nNode * 3 * _nNode * 3);
			Eigen::SparseMatrix<double, Eigen::ColMajor>* _mat = &mat->_mat[0];
			_mat->setZero();
			_mat->makeCompressed();
			//_mat.reserve(_nNode * 3 * _nNode * 3);
			for (int i = 0; i < _nNode; i++)
			{
				int I = _index[i] * 3;

				for (int j = 0; j < _nNode; j++)
				{
					int J = _index[j] * 3;
					for (int k = 0; k < 3; k++)
					{
						for (int k2 = 0; k2 < 3; k2++)
						{
							double _val3 = 0;
							/*for (int l = 0; l < 2; l++)
							{
								for (int m = 0; m < 2; m++)
								{
									for (int g = 0; g < 2; g++)
									{
										for (int h = 0; h < 2; h++)
										{
											double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));
											double D = (_ref->d2[(g*2) + h][j] - Gammaijk[(((g * 2) + h) * 2) + 0] * _ref->d1[0][j] - Gammaijk[(((g * 2)+ h) * 2)+ 1]* _ref->d1[1][j]);
											double E = (_ref->d2[(l * 2) + m][i] - Gammaijk[(((l * 2) + m) * 2) + 0] * _ref->d1[0][i] - Gammaijk[(((l * 2)+ m) * 2)+ 1] * _ref->d1[1][i]);
											_val3 += A * N[k] * N[k2] * (D) * (E);
										}
									}
								}
							}*/
							_val3 = K(i, k, j, k2, _la, _mu);
							//dat[(i*3+k) * (_nNode * 3) + (j*3 + k2)] = Eigen::Triplet<double>(I + k, J + k2, _val3 * sc);

							dat.push_back(Eigen::Triplet<double>(I+k, J+k2, _val3 * sc));
						}
					}
				}
			}
			_mat->setFromTriplets(dat.begin(),dat.end());
			M->_mat[0] += *_mat;
			//_mat->makeCompressed();
			//_mat->finalize();
			//mat->plus(&_mat);
		}
		//ultimate term
		int star(int i) {
			if (i == 0)return 3;
			if (i == 1)return 2;
			if (i == 2)return 1;
			if (i == 3)return 0;
			return 0;
		}

		double U(int I, int J) {
			double val = 0;
			//double sc = 1.0 / _ref->_refDv / _ref->_refDv;
			/*for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					int ij = i * 2 + j;
					val += sc*star2[ij]*(d2[star(ij)][I]-Gammaijk[(star(ij))*2+0]*d1[0][I] - Gammaijk[(star(ij)) * 2 + 1] * d1[1][I]) * 
						(d2[ij][J] - Gammaijk[(ij) * 2 + 0] * d1[0][J] - Gammaijk[(ij) * 2 + 1] * d1[1][J]);
				}
			}*/
			return sc*_ref->__mat[I * _nNode + J];
		}
		/*double U2(int I, int J) {
			double sc = 1 / _ref->_refDv / _ref->_refDv;
			return sc * sc * __mat2[I * _nNode + J];
		}*/
		void U_z(double* ptr)
		{
			/*for (int i = 0;i < _nNode;i++) {
				*(ptr + i) = sc*__grad_z[i];
			}*/
			memcpy(ptr, __grad_z, sizeof(double) * _nNode);
		}
		void U_phi(double* ptr)
		{
			/*for (int i = 0;i < _nNode;i++) {
				*(ptr + i) = sc*__grad_phi[i];
			}*/
			memcpy(ptr, __grad_phi, sizeof(double) * _nNode);
		}

		double U_z(int I) {
			//double val = 0;
			//double sc = 1.0 / _ref->_refDv / _ref->_refDv;
			return __grad_z[I];
		}
		double U_phi(int J) {
			//double val = 0;
			//double sc = 1.0 / _ref->_refDv / _ref->_refDv;
			return  __grad_phi[J];
		}
		//membrane boundary term
		double MB(int i,int k)
		{
			double _val3 = 0;
			for (int j = 0; j < 2; j++) {
				_val3 += _Sij[(i<<1)+j] * this->get_gi(j, k);
			}	
			return _val3 * _ref->_refDv;
		}
		//membrane boundary term
		double MB3(int i,int I)
		{
			double _val3 = 0;
			for (int j = 0; j < 2; j++) {
				_val3 += _Sij[(i <<1) + j] * _ref->d1[j][I];
			}
			return _val3;
		}
		//membrane boundary term
		double MB2(int i, int k)
		{
			double _val3 = 0;
			for (int j = 0; j < 2; j++) {
				_val3 += _Sij[(i <<1) + j] * this->get_gi(j, k);
			}
			return _val3;
		}
		double S(int i, int j) {
			return _Sij[(i<<1) + j];
		}
		//bending boundary term
		double KB(int j, int k,int l,int m, double _la, double _mu)
		{
			double _val3 = 0;

			
			for (int g = 0; g < 2; g++)
			{
				for (int h = 0; h < 2; h++)
				{
					double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2*_mu*_ref->get__Gij(h, m) * _ref->get__Gij(g, l));
					double D = (_ref->d2[g * 2 + h][j] - Gammaijk[(g * 2 + h) * 2 + 0] * _ref->d1[0][j] - Gammaijk[(g * 2 + h) * 2 + 1] * _ref->d1[1][j]);
					_val3 += A * N[k] * (D);
				}
			}

			return _val3 * _ref->refDv * _ref->refDv;
		}
		//angle term
		double T(int i, int s,int m) {
			double val = 0;
			for (int l = 0; l < 2; l++)
			{

				val+= _ref->d1[l][i]*N[s] * _ref->get__Gij(l, m);
			}
			return val* _ref->refDv;
		}
		//membrane boundary handle
		double T2(int l, int i) {
			double val = 0;
			for (int m = 0; m < 2; m++) {
				val += _ref->d1[l][i] * get_Gij(l, m);
			}
			return val;
		}
		//membrane term
		double H(int i, int k2, int j, int k, double _la, double _mu)
		{
			const static int kk[3]{ 0,1,2 };
			const static int ll[2]{ 0,1};

			double _val4 = 0;
			
			for (auto l:ll)
			{
				for (auto m : ll)
				{
					for (auto g:ll)
					{
						for (auto h:ll)
						{

							double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) +2*_mu* _ref->get__Gij(h, m) * _ref->get__Gij(g, l));

							double FF = (_ref->d1[g][j] * get_gi(h, k) + _ref->d1[h][j] * get_gi(g, k));
							double GG = (_ref->d1[l][i] * get_gi(m, k2) + _ref->d1[m][i] * get_gi(l, k2));
							_val4 += A*FF * GG;
						}
					}
				}
			}
			return _val4 * _ref->refDv * 0.25;
		}
		void H(_mySparse* M,_mySparse * mat,int* _index,double _la, double _mu,double sc)
		{
			const static int kk[3]{ 0,1,2 };
			const static int ll[2]{ 0,1 };
			//double _sc = sc * _ref->refDv * 0.25;
			std::vector<Eigen::Triplet<double>> dat;//
			dat.resize(_nNode * 3 * _nNode * 3);
			//dat.clear();
			Eigen::SparseMatrix<double,Eigen::ColMajor>*_mat = &mat->_mat[0];
			_mat->setZero();
			_mat->makeCompressed();
			//_mat.reserve(_nNode * 3 * _nNode * 3);
			for (int i = 0; i < _nNode; i++)
			{
				int I = _index[i]*3;
				for (int j = 0; j < _nNode; j++)
				{
					int J = _index[j] * 3;
					for (int k=0;k<3;k++)
					{
						for (int k2 = 0; k2 < 3; k2++)
						{
							double _val4 = 0;
							/*for (int l = 0; l<2; l++)
							{
								for (int m = 0; m < 2; m++)
								{
									for (int g = 0; g < 2; g++)
									{
										for (int h = 0; h < 2; h++)
										{
											double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));

											double FF = (_ref->d1[g][j] * get_gi(h, k2) + _ref->d1[h][j] * get_gi(g, k2));
											double GG = (_ref->d1[l][i] * get_gi(m, k) + _ref->d1[m][i] * get_gi(l, k));
											_val4 += A * FF * GG;
										}
									}
								}
							}*/
							_val4 = H(i, k, j, k2, _la, _mu);
							//_mat.insert(I + k, J + k2) = _val4 * _sc;
							dat[(i*3+k)*(_nNode*3)+(j*3+k2)]=Eigen::Triplet<double>(I + k, J + k2, _val4 * sc);
							//mat->_plus(I+k, J+k2, _val4 * _sc);
						}
					}
				}
			}
			_mat->setFromTriplets(dat.begin(),dat.end());
			M->_mat[0] += *_mat;
			//_mat->makeCompressed();
			//_mat->finalize();
			//mat->plus(&_mat);
			//return _val4 * _ref->refDv * 0.25;
		}
		void memory(_memS_ref* __mem) {
			if (__mem->__z < -1000) {
				__mem->__z = z;
			}
			__mem->_x = x;
			__mem->_y = y;
			__mem->_z = z;
			__mem->Z = Z;
			__mem->_Z = Z;
			__mem->refDv = dv;
			__mem->_refDv = _dv;
			std::memcpy(__mem->_Sij, _Sij, sizeof(double) * 4);
			std::memcpy(__mem->_gi, gi, sizeof(double) * 6);
			std::memcpy(__mem->_Gi, Gi, sizeof(double) * 6);
			std::memcpy(__mem->_gij, gij, sizeof(double) * 4);
			std::memcpy(__mem->_Gij, Gij, sizeof(double) * 4);
			//std::memset(__mem->_gij, 0, sizeof(double) * 4);
			//std::memset(__mem->_Gij, 0, sizeof(double) * 4);
			std::memcpy(__mem->_bij, bij, sizeof(double) * 12);
			std::memcpy(__mem->_Gammaijk, Gammaijk, sizeof(double) * 8);

		}
	};

public ref class memS_ref {
public:
	_memS_ref* __mem=0;
	void setbuffer(buffer^ buf) {
		__mem->set_buffer(buf->_buf->mem);
	}
	double _x() {
		return __mem->_x;
	}
	double _y() {
		return __mem->_y;
	}
	double _z() {
		return __mem->_z;
	}
	double _Z() {
		return __mem->_Z;
	}
	memS_ref() {
		__mem = new _memS_ref();
	}
	void dispose() {
		if(__mem!=0)
		delete __mem;
		__mem = 0;
	}
	~memS_ref() {
		dispose();
	}
	!memS_ref() {
		dispose();
	}
	void set_z(double z)
	{
		__mem->set_z(z);
	}
	void update_z_phi(int nNode, KingOfMonsters::myDoubleArray^ Z, KingOfMonsters::myDoubleArray^ phi) {
		if (Z != nullptr)
		{
			//eigen_assert(Z->_arr->__v.norm() < 10000 && Z->_arr->__v.norm() > -100000);
			__mem->set_buf_z(Z->_arr->__v.data(), nNode);
			/*for (int i = 0; i < nNode; i++) {
				__mem->set_buf_z(i, (*Z->_arr)(i));
			}*/
		}
		if (phi != nullptr)
		{
			//eigen_assert(phi->_arr->__v.norm() < 10000 && phi->_arr->__v.norm() > -100000);
			__mem->set_buf_phi(phi->_arr->__v.data(), nNode);
			/*for (int i = 0; i < nNode; i++) {
				__mem->set_buf_phi(i, (*phi->_arr)(i));
			}*/
		}
	}
	void update3(int nNode, KingOfMonsters::myDoubleArray^ node, KingOfMonsters::myDoubleArray^ def,bool ignoreZ) {
		if (node != nullptr) {
			__mem->set_node(node->_arr->__v.data(), nNode * 3);
			/*for (int i = 0; i < nNode; i++) {
				int e = i * 3;
				__mem->set_node(i, 0, (*node->_arr)(e + 0));
				__mem->set_node(i, 1, (*node->_arr)(e + 1));
				if (!ignoreZ) {
					__mem->set_node(i, 2, (*node->_arr)(e + 2));
				}
				else{
					__mem->set_node(i, 2, 0);
				}
			}*/
		}
		if (def != nullptr)
		{
			__mem->set_def(def->_arr->__v.data(), nNode * 3);

			/*for (int i = 0; i < nNode; i++) {
				int e = i * 3;
				__mem->set_def(i, 0, def[e + 0]);
				__mem->set_def(i, 1, def[e + 1]);
				if (!ignoreZ) {
					__mem->set_def(i, 2, def[e + 2]);
				}
				else {
					__mem->set_def(i, 2, 0);
				}
			}*/
		}
	}
	void update3(int nNode, KingOfMonsters::myDoubleArray^ node, KingOfMonsters::myDoubleArray^  def) {
		update3(nNode, node, def, false);
	}
	void update(int nNode, int uDim, int vDim) {
		__mem->update(nNode, uDim, vDim);
	}
};
public ref class memS
{
public:
	_memS* __mem=0;
public:
	double x, y, z,Z,phi;
	double _x, _y, _z,_Z,_phi;
	double refDv;
	double _refDv;
	double dv;
public:
	double orient(memS^ another) {
		double dot = another->__mem->N[0] * this->__mem->N[0] + another->__mem->N[1] * this->__mem->N[1] + another->__mem->N[2] * this->__mem->N[2];
		if (dot < 0)return -1;
		return 1;
	}
	void setRef(memS_ref^ _mem_) {
		__mem->_ref = _mem_->__mem;
	}
	void computeGrads() {
		__mem->computeGrads();
	}
	double B(int i, int j) {
		return __mem->_B(i, j);
	}
	//body force accurate area accounting for the slope
	double G2(int i) {
		return __mem->G2(i);
	}
	//body force least square
	double G3(int i) {
		return __mem->G3(i);
	}
	//body force projected area
	double G4(int i) {
		return __mem->G4(i);
	}
	double R(int i) {
		return __mem->R(i);
	}
	double A() {
		return __mem->A();
	}
	double dA(int i) {
		return __mem->dA(i);
	}
	double R2(int i) {
		return __mem->R2(i);
	}
	double R3(int i) {
		return __mem->R3(i);
	}
	double MASS(int i, int j) {
		return __mem->MASS(i, j);
	}
	double SMOOTH(int i, int j) {
		return __mem->SMOOTH(i, j);
	}
	double SMOOTH_z(int i) {
		return __mem->SMOOTH_z(i);
	}
	double SMOOTH_phi(int i) {
		return __mem->SMOOTH_phi(i);
	}
	double area() {
		return __mem->area();
	}
	double U(int I, int J) {
		return __mem->U(I, J);
	}
	double U_z(int I) {
		return __mem->U_z(I);
	}
	double U_phi(int J) {
		return __mem->U_phi(J);
	}
	void d0(mySparse^ mat, int ii, myIntArray^ index, double sc)
	{
		mat->dat->addrow(ii, index->_arr, __mem->_ref->d0, sc, __mem->_nNode);
	}
	void U_z(mySparse^ mat, int ii, myIntArray^ index,double sc)
	{

		mat->dat->addrow(ii, index->_arr, __mem->__grad_z, sc,__mem->_nNode);

	}
	void U_phi(mySparse^ mat, int ii, myIntArray^ index,double sc)
	{

		mat->dat->addrow(ii, index->_arr, __mem->__grad_phi, sc, __mem->_nNode);

	}
	void U_z(array<double>^ X) {
		System::Runtime::InteropServices::Marshal::Copy((IntPtr)__mem->__grad_z, X, 0, __mem->_nNode);
	}
	void U_phi(array<double>^ x) {
		System::Runtime::InteropServices::Marshal::Copy((IntPtr)__mem->__grad_phi, x, 0, __mem->_nNode);
	}
	void U_z(myDoubleArray^ X) {
		memcpy(X->_arr->__v.data(),__mem->__grad_z,sizeof(double)* __mem->_nNode);
	}
	void U_phi(myDoubleArray^ x) {
		memcpy(x->_arr->__v.data(), __mem->__grad_phi, sizeof(double) * __mem->_nNode);
	}
	array<double>^ fM(double _la, double _mu)
	{
		double a, b, c, d;
		__mem->fM(_la,_mu,&a,&b,&c,&d);
		array<double>^ arr = gcnew array<double>(4);
		arr[0] = a;
		arr[1] = b;
		arr[2] = c;
		arr[3] = d;
		return arr;
	}
		/*void fM(array<double>^ ret) {
			vector<double> _ret = __mem->fM();
			ret[0] = _ret[0];
			ret[1] = _ret[1];
			ret[2] = _ret[2];
			ret[3] = _ret[3];
			ret[4] = _ret[4];
			ret[5] = _ret[5];

		}*/
	double SLOPE_phi(int l) {
		return __mem->SLOPE_phi(l);
	}
	double SLOPE_z(int l) {
		return __mem->SLOPE_z(l);
	}
	double SLOPE_phi(double dcdtstar0,double dcdtstar1,double sc) {
		return __mem->SLOPE_phi(dcdtstar0,dcdtstar1,sc);
	}
	double SLOPE_z(double dcdtstar0, double dcdtstar1,double sc) {
		return __mem->SLOPE_z(dcdtstar0, dcdtstar1,sc);
	}
	double SLOPE(int l, int i) {
		return __mem->SLOPE(l, i);
	}

	void SLOPE(myDoubleArray^ arr,double dcdt1,double dcdt2,double sc) {
		__mem->SLOPE(arr->_arr->__v.data(), dcdt1, dcdt2, sc);
		//return __mem->SLOPE(l, i);
	}
	void SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc,int shift) {
		__mem->SLOPE(arr->_arr->__v.data()+shift, dcdt1, dcdt2, sc);
		//return __mem->SLOPE(l, i);
	}
	double eM(double _la, double _mu) {
			return __mem->eM(_la,_mu);
		}
		double eK(double _la, double _mu) {
			return __mem->eK(_la,_mu);
		}
		double L(int a, double _la, double _mu) {
			return __mem->L(a, _la, _mu);
		}
		double M(int a, double _la, double _mu) {
			return __mem->_M(a, _la, _mu);
		}
		double d0(int i) {
			return __mem->_d0(i);
		}
		double d2(int l,int i) {
			return __mem->_d2(l,i);
		}
		double gi(int i,int s){
			return __mem->get_gi(i, s);
		}
		///error norm
		double error() {
			return __mem->error();
		}
		///stress function
		double F(int i,int j) {
			return __mem->F(i,j);
		}
		//least squares
		double F2(int i, int j) {
			return __mem->F2(i, j);
		}
		double F3(int i, int I) {
			return __mem->F3(i, I);
		}
		double F3_phi(int i) {
			return __mem->F3_phi(i);
		}
		double F3_z(int i) {
			return __mem->F3_z(i);
		}
		double Dc(double a, double b,int l) {
			return __mem->Dc(a, b, l);
		}
		double F4(int i, int I, int J) {
			return __mem->F4(i, I, J);
		}
		double F4_z(int i, int I)
		{
			return __mem->F4_z(i, I);
		}
		double F4_phi(int i, int J)
		{
			return __mem->F4_phi(i, J);
		}
		void F4_z(myDoubleArray^ arr,double dcdtstar0,double dcdtstar1,double sc)
		{
			__mem->F4_z(arr->_arr->__v.data(),dcdtstar0,dcdtstar1,sc);
		}
		void F4_phi(myDoubleArray^ arr, double dcdtstar0, double dcdtstar1, double sc)
		{
			__mem->F4_phi(arr->_arr->__v.data(), dcdtstar0, dcdtstar1, sc);
		}
		double F5_z(int i, int I)
		{
			return __mem->F5_z(i, I);
		}
		double F5_phi(int i, int J)
		{
			return __mem->F5_phi(i, J);
		}

		double F6(int i, int I, int J) {
			return __mem->F6(i, I, J);
		}
		double F6_z(int i, int I)
		{
			return __mem->F6_z(i, I);
		}
		double F6_phi(int i, int J)
		{
			return __mem->F6_phi(i, J);
		}

		double S(int i, int j) {
			return __mem->S(i, j);
		}
		
		//body force reference area
		double G(int i) {
			return __mem->G(i);
		}
		//bending term
		double K(int i, int s, int j, int k,double _la,double _mu) {
			return __mem->K(i, s, j, k,_la,_mu);
		}
		///membrane boundary term
		double MB(int i, int k)
		{
			return __mem->MB(i, k);
		}
		///membrane boundary term2
		double MB2(int i, int k)
		{
			return __mem->MB2(i, k);
		}
		double MB3(int l,int i)
		{
			return __mem->MB3(l,i);
		}
		///<summary>
		///bending boundary term
		///</summary>
		double KB(int j, int k, int l, int m, double _la, double _mu) {
			return __mem->KB(j, k, l, m,_la,_mu);
		}
		//rotation angle at boundary
		double T(int i, int s, int m) {
			return __mem->T(i, s, m);
		}
		//membrane boundary term handle
		double T2(int l, int i) {
			return __mem->T2(l, i);
		}
		//membrane term
		double H(int i, int s, int j, int k, double _la, double _mu) {
			return __mem->H(i, s, j, k, _la, _mu);
		}
		double K(int l, int I)
		{
			return __mem->K(l, I);
		}
		double K_phi(int j)
		{
			return __mem->K_phi(j);
		}
		void H(mySparse^ M, mySparse^ mat, myIntArray^ index, double _la, double _mu, double sc)
		{
			__mem->H(M->dat,mat->dat, index->data(), _la, _mu, sc);
		}
		void K(mySparse^ M,mySparse^ mat, myIntArray^ index, double _la, double _mu, double sc)
		{
			__mem->K(M->dat,mat->dat, index->data(), _la, _mu, sc);
		}
		void S(double val1, double val2, double val3) {
			__mem->set_Sij(val1, val2, val3);
		}

		double bodyF() {
			double sc = 1.0/__mem->_ref->_refDv/ __mem->_ref->_refDv;
			return sc*__mem->bodyF;
		}
		double BCF(int l) {
			double sc = 1.0 / __mem->_ref->_refDv / __mem->_ref->_refDv;
			return  sc*__mem->BCF[l];
		}
		double BCF2(int l) {
			double sc = 1.0 / __mem->_ref->_refDv / __mem->_ref->_refDv;
			return sc*__mem->BCF2[l];
		}
		double BCF6(int l) {
			return __mem->BCF6[l];
		}
		double dcdtstar(double x, double y, int i)
		{
			return x* __mem->get_gi(i, 0) + y * __mem->get_gi(i, 1);
		}
		void compute() {
			__mem->update2();
			x = __mem->x;
			y = __mem->y;
			z = __mem->z;
			Z = __mem->Z;
			phi = __mem->phi;
			dv = __mem->dv;
			_refDv = __mem->_dv;
			refDv = __mem->dv;
		}
		void update_lo(array<double>^ lo) {
			__mem->set_lo(lo[0], lo[1]);
		}
		void update_elem(int nNode, int uDim, int vDim,array<double, 3>^ M, array<int, 2>^ dd) {
			__mem->update(nNode, uDim, vDim);
			int i = 0;
			for (int j = 0; j < uDim; j++) {
				for (int k = 0; k < uDim; k++) {
					__mem->set_M(i, j, k, M[i, j, k]);
				}
			}
			i = 1;
			for (int j = 0; j < vDim; j++) {
				for (int k = 0; k < vDim; k++) {
					__mem->set_M(i, j, k, M[i, j, k]);
				}
			}
			for (int i = 0; i < nNode; i++) {
				for (int j = 0; j < 2; j++) {
					__mem->set_dd(i, j, dd[i, j]);
				}
			}
		}
		void memory(memS_ref^ _mem_) {
			__mem->memory(_mem_->__mem );
			_x = __mem->_ref->_x;
			_y = __mem->_ref->_y;
			_z = __mem->_ref->_z;
			Z = __mem->_ref->Z;
			//phi = __mem->phi;
		}
		memS() {
			__mem = new _memS(std::string(""));
		}
		void static func(int i) {

		}
		memS(System::String ^ultimate) {
			
			if (ultimate == "U")
				__mem = new _memS("U");
			else if (ultimate == "SLOPE")
				__mem = new _memS("SLOPE");
			else if (ultimate == "SHELL")
				__mem = new _memS("SHELL");
			else if (ultimate == "SENSITIVITY")
				__mem = new _memS("SENSITIVITY");
			else
				__mem = new _memS("");
		}
		void dispose() {
			if(__mem!=0)
			delete __mem ;
			__mem = 0;
		}
		~memS() {
			dispose();
		}
		!memS() {
			dispose();
		}
	};


}
