#pragma once
#include <cmath>
#include<vector>
#include<math.h>
using namespace System;
#include <cstring>
#include <string>_SLOPE_z
using std::vector;
using std::string;

#include "mySparseLibrary.h"
//#define EIGEN_DONT_PARALLELIZE
//#define EIGEN_MALLOC_ALREADY_ALIGNED  0

namespace KingOfMonsters {

	public class _buffer {
	public:
		double mem[12000];
	};
	public ref class buffer {
	public:
		_buffer* _buf = 0;

		buffer() {
			_buf = new _buffer();
		}
		~buffer() {
			if (_buf != 0)
				delete _buf;
			_buf = 0;
		}
		!buffer() {
			if (_buf != 0)
				delete _buf;
			_buf = 0;
		}
	};
	const int ___ll[4]{ 0,3,6,9 };
	enum _RAM {
		SAVE = 1,
		MAX = 2
	};

	public class _memS_ref {
	public:
		_RAM RAM;
		int _nNode;
		int dim[2]{ 0,0 };
		int _uDim, _vDim;
		double* node = 0;
		double* def = 0;
		double* buf_z = 0;
		double* buf_phi = 0;
		double* buf_xi = 0;
		double* buf_eta = 0;
		double* buf_mu = 0;
		double* buf_nu = 0;

		double _gi[6];
		double _Gi[6];
		double _gij[4];
		double _Gij[4];
		double _bij[12];
		double _Sij[4];
		double _Gammaijk[8];

		double* __mat = 0;
		double** M[2]{ 0,0 };
		int** dd = 0;
		double* d0 = 0;


		//double* d1[2]{ 0,0 };
		//double* d2[4]{ 0,0,0,0 };
		//double* d2_star[4]{ 0,0,0,0 };

		//double* B[4]{ 0,0,0,0 };
		//double* tt0[2]{ 0,0 }, * hh0[2]{ 0,0 }, * tt1[4]{ 0,0,0,0 }, * hh1[4]{ 0,0,0,0 }, * tt2[8]{ 0,0,0,0,0,0,0,0 }, * hh2[8]{ 0,0,0,0,0,0,0,0 };

		double** d1 = 0; //2
		double** d2 = 0; //4
		double** d2_star = 0;  //4

		double** B = 0; //4
		double** tt0 = 0, ** hh0 = 0;//2
		double** tt1 = 0, ** hh1 = 0; //4
		double** tt2 = 0, ** hh2 = 0;//8

	public:
		bool initialized = false;
		double refDv = 0, _refDv = 0;
		double _x = 0, _y = 0, _z = 0, __z = 0, Z = 0, _Z = 0;
		double xi = 0, eta = 0,mu=0;
		//double w = 0;
		inline void set_z(double& z) {
			this->_z = z;
		}
		inline void set__z(double& z) {
			this->__z = z;
		}
		inline void set_buffer(double* buf)
		{
			node = buf;
			def = &buf[1500];
			buf_z = &buf[3000];
			buf_phi = &buf[4500];
			buf_xi = &buf[6000];
			buf_eta = &buf[7500];
			buf_mu = &buf[9000];
			buf_nu = &buf[10500];
			//buf_chi = &buf[9000];
			//buf_b = &buf[12000];
			//buf_D = &buf[14000];
			//buf_W = &buf[8000];
		}
		inline void set_node(const int& i, const int& s, const double val) {
			node[___ll[i] + s] = val;
		}
		inline void set_node(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(node, ptr, sizeof(double) * N);
		}
		inline void set_buf_z(const int& i, const double& val) {
			buf_z[i] = val;
		}
		inline void set_buf_z(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_z, ptr, sizeof(double) * N);
		}
		inline void set_buf_phi(const int& i, const double val) {
			buf_phi[i] = val;
		}
		inline void set_buf_phi(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_phi, ptr, sizeof(double) * N);
		}
		inline void set_def(const int& i, const int& s, const double& val) {
			def[___ll[i] + s] = val;
		}
		inline void set_def(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(def, ptr, sizeof(double) * N);
		}
		inline void set_buf_xi(const int& i, const double& val) {
			buf_xi[i] = val;
		}
		inline void set_buf_xi(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_xi, ptr, sizeof(double) * N);
		}
		inline void set_buf_eta(const int& i, const double& val) {
			buf_eta[i] = val;
		}
		inline void set_buf_eta(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_eta, ptr, sizeof(double) * N);
		}
		inline void set_buf_mu(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_mu, ptr, sizeof(double) * N);
		}
		inline void set_buf_nu(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_nu, ptr, sizeof(double) * N);
		}
		inline double& get_node(const int& i, const int& s) {
			return node[___ll[i] + s];
		}
		inline double& get__gi(const int& i, const int& s) {
			return _gi[___ll[i] + s];
		}
		inline double& get__Gi(const int& i, const int& s) {
			return _Gi[___ll[i] + s];
		}
		inline double& get__gij(const int& i, const int& j) {
			return _gij[(i << 1) + j];
		}
		inline double& get__Gij(const int& i, const int& j) {
			return _Gij[(i << 1) + j];
		}
		inline double& get__bij(const int& i, const int& j, const int& s) {
			return _bij[___ll[((i << 1) + j)] + s];
		}
		inline double& get__Gammaijk(const int& i, const int& j, const int& k) {
			return _Gammaijk[(((i << 1) + j) << 1) + k];
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
			d1 = 0;
			d2 = 0;
			d2_star = 0;
			/*d1[0] = 0;
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
			B[3] = 0;*/

			B = 0;
			tt0 = 0;
			hh0 = 0;
			tt1 = 0;
			hh1 = 0;
			tt2 = 0;
			hh2 = 0;

			/*
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
			*/
		}
		~_memS_ref() {
			del();
		}

		void del() {
			if (RAM == MAX)
			{
				if (dd != 0)
				{
					delete[] M[0];
					delete[] M[1];
					M[0] = 0;
					M[1] = 0;
					//}
					//if (dd != 0)
					//{
					for (int i = 0; i < _nNode; i++) {
						delete[] dd[i];
					}
					delete[] dd;
					dd = 0;
				}
				if (__mat != 0)delete[] __mat;
				__mat = 0;
				if (d0 != 0) {
					delete[] d0;
					d0 = 0;
				}
				if (d1 != 0)
				{
					delete[] d1[0];
					delete[] d1[1];
					//d1[0] = 0;
					//d1[1] = 0;
					delete[] d1;
					d1 = 0;

				}
				if (d2 != 0)
				{
					delete[] d2[0];
					delete[] d2[1];
					delete[] d2[2];
					delete[] d2[3];
					//d2[0] = 0;
					//d2[1] = 0;
					//d2[2] = 0;
					//d2[3] = 0;
					delete[] d2;
					d2 = 0;
					delete[] d2_star[0];
					delete[] d2_star[1];
					delete[] d2_star[2];
					delete[] d2_star[3];
					//d2_star[0] = 0;
					//d2_star[1] = 0;
					//d2_star[2] = 0;
					//d2_star[3] = 0;
					delete[] d2_star;
					d2_star = 0;
				}

				if (B != 0) {
					delete[] B[0];
					delete[] B[1];
					delete[] B[2];
					delete[] B[3];
					//B[0] = 0;
					//B[1] = 0;
					//B[2] = 0;
					//B[3] = 0;
					delete[] B;
					B = 0;
				}
				if (tt0 != 0) {
					delete[] tt0[0];
					//tt0[0] = 0;
					delete[] tt0[1];
					//tt0[1] = 0;
					delete[] tt0;
					tt0 = 0;

					delete[] hh0[0];
					delete[] hh0[1];
					//hh0[0] = 0;
					//hh0[1] = 0;
					delete[] hh0;
					hh0 = 0;

					delete[] tt1[0];
					delete[] tt1[1];
					delete[] tt1[2];
					delete[] tt1[3];
					//tt1[0] = 0;
					//tt1[1] = 0;
					//tt1[2] = 0;
					//tt1[3] = 0;
					delete[] tt1;
					tt1 = 0;

					delete[] hh1[0];
					delete[] hh1[1];
					delete[] hh1[2];
					delete[] hh1[3];
					//hh1[0] = 0;
					//hh1[1] = 0;
					//hh1[2] = 0;
					//hh1[3] = 0;
					delete[] hh1;
					hh1 = 0;

					delete[] tt2[0];
					delete[] tt2[1];
					delete[] tt2[2];
					delete[] tt2[3];
					delete[] tt2[4];
					delete[] tt2[5];
					delete[] tt2[6];
					delete[] tt2[7];
					//tt2[0] = 0;
					//tt2[1] = 0;
					//tt2[2] = 0;
					//tt2[3] = 0;
					//tt2[4] = 0;
					//tt2[5] = 0;
					//tt2[6] = 0;
					//tt2[7] = 0;
					delete[] tt2;
					tt2 = 0;

					delete[] hh2[0];
					delete[] hh2[1];
					delete[] hh2[2];
					delete[] hh2[3];
					delete[] hh2[4];
					delete[] hh2[5];
					delete[] hh2[6];
					delete[] hh2[7];
					//hh2[0] = 0;
					//hh2[1] = 0;
					//hh2[2] = 0;
					//hh2[3] = 0;
					//hh2[4] = 0;
					//hh2[5] = 0;
					//hh2[6] = 0;
					//hh2[7] = 0;
					delete[] hh2;
					hh2 = 0;
					//}
				}
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
				if (RAM == MAX)
				{
					__mat = new double[_nNode * _nNode];
					d0 = new double[nNode];

					d1 = new double* [2];
					d1[0] = new double[nNode];
					d1[1] = new double[nNode];

					d2 = new double* [4];
					d2[0] = new double[nNode];
					d2[1] = new double[nNode];
					d2[2] = new double[nNode];
					d2[3] = new double[nNode];

					d2_star = new double* [4];
					d2_star[0] = new double[nNode];
					d2_star[1] = new double[nNode];
					d2_star[2] = new double[nNode];
					d2_star[3] = new double[nNode];

					B = new double* [4];
					B[0] = new double[nNode * nNode];
					B[1] = new double[nNode * nNode];
					B[2] = new double[nNode * nNode];
					B[3] = new double[nNode * nNode];

					tt0 = new double* [2];
					tt0[0] = new double[uDim];
					tt0[1] = new double[vDim];

					hh0 = new double* [2];
					hh0[0] = new double[uDim];
					hh0[1] = new double[vDim];

					tt1 = new double* [4];
					tt1[0] = new double[uDim];
					tt1[1] = new double[vDim];
					tt1[2] = new double[uDim];
					tt1[3] = new double[vDim];
					hh1 = new double* [4];
					hh1[0] = new double[uDim];
					hh1[1] = new double[vDim];
					hh1[2] = new double[uDim];
					hh1[3] = new double[vDim];

					tt2 = new double* [8];
					tt2[0] = new double[uDim];
					tt2[1] = new double[vDim];
					tt2[2] = new double[uDim];
					tt2[3] = new double[vDim];
					tt2[4] = new double[uDim];
					tt2[5] = new double[vDim];
					tt2[6] = new double[uDim];
					tt2[7] = new double[vDim];
					hh2 = new double* [8];
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
		_RAM RAM;

		_memS_ref* _ref = 0;
		double N[3];
		double bodyF;
		double BCF[2];
		double BCF2[2];
		double BCF6[2];
	public:
		double* __grad_z = 0;
		double* __grad_phi = 0;
		double* __grad = 0;
		double* __grad2 = 0;
		double* __grad_z_tmp = 0;
		double* __grad_phi_tmp = 0;
		double* __mix_phi = 0;
		double* __mix_z = 0;

		double* __grad_xi = 0;
		double* __grad_eta = 0;

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
		double* __grad_C_chi = 0;
		double* __grad_C_psi = 0;
		double* __grad_D_z = 0;
		double* __grad_D_phi = 0;
		double* __grad_D_chi = 0;
		double* __grad_D_psi = 0;
		double _K_phi[2];
		double _K_chi[2];
		double _K_psi[2];
		double* _K = 0;
		double lo[2]{ 0,0, };


		int dim[2];
		int _uDim, _vDim;
	public:
		double gi[6];
		double gi2[6];
		double Gi[6];
		double Gi2[6];
		double gij[4];
		double gij2[4];
		double Gij2[4];
		double Gij[4];
		double bij[12];
		double Gammaijk[8];
		double Gammaijk2[8];
		double _ss[4];
		double _Sij[4];
	public:
		double sc;
	public:
		///////shared memory//////
		double** M[2];
		int** dd;
		double* __mat = 0;
		double* d0;
		double* d1[2];
		double* d2[4];
		double* d2_star[4];
		const int star2[4]{ 1,-1,-1,1 };
		const int _star[4]{ 3,1,2,0 };
		double* B[4];
		double* tt0[2], * hh0[2], * tt1[4], * hh1[4], * tt2[8], * hh2[8];
		///////shared memory//////

		double* gradN[3]{ 0,0,0 };
		double* gradG = 0;
		const int ___ll[4]{ 0,3,6,9 };
	public:
		double x, y, z, Z, _Z, phi, eta, xi,mu,nu;/// , chi;
		double dv, _dv;
	public:
		inline void set_lo(double L1, double L2) {
			lo[0] = L1;
			lo[1] = L2;
		}
		inline void set_M(int i, int j, int k, double val) {
			_ref->M[i][j][k] = val;
		}
		inline void set_dd_save(int i, int j, int val) {
			dd[i][j] = val;
		}
		inline void set_M_save(int i, int j, int k, double val) {
			M[i][j][k] = val;
		}
		inline void set_dd(int i, int j, int val) {
			_ref->dd[i][j] = val;
		}
		inline double component(double* G1, double* G2) {
			return G1[0] * _ss[0] * G2[0] + G1[1] * _ss[1] * G2[0] + G1[0] * _ss[1] * G2[1] + G1[1] * _ss[3] * G2[1];
		}
		inline void set_Sij(const double& A, const double& B, const double& C) {
			_ss[0] = C;
			_ss[1] = -B;
			_ss[2] = -B;
			_ss[3] = A;

			_Sij[0] = component(Gi, Gi);
			_Sij[3] = component(&(Gi[3]), &(Gi[3]));
			_Sij[1] = component(&(Gi[3]), Gi);
			_Sij[2] = component(Gi, &(Gi[3]));
		}

		inline double& get_gi(const int& i, const int& s) {
			return gi[___ll[i] + s];
		}
		inline double& get_gi2(const int& i, const int& s) {
			return gi2[___ll[i] + s];
		}
		inline double& get_Gi(const int& i, const int& s) {
			return Gi[___ll[i] + s];
		}
		inline double& get_Gi2(const int& i, const int& s) {
			return Gi2[___ll[i] + s];
		}
		inline double& get_gij(const int& i, const int& j) {
			return gij[(i << 1) + j];
		}
		inline double& get_gij2(const int& i, const int& j) {
			return gij2[(i << 1) + j];
		}
		inline double& get_Gij(const int& i, const int& j) {
			return Gij[(i << 1) + j];
		}
		inline double& get_Gij2(const int& i, const int& j) {
			return Gij2[(i << 1) + j];
		}
		inline double& get_bij(const int& i, const int& j, const int& s) {
			return bij[___ll[((i << 1) + j)] + s];
		}
		inline double& get_Gammaijk(const int& i, const int& j, const int& k) {
			return Gammaijk[(((i << 1) + j) << 1) + k];
		}
		inline double& get_Gammaijk2(const int& i, const int& j, const int& k) {
			return Gammaijk2[(((i << 1) + j) << 1) + k];
		}

		inline double& get_tt0(const int& i, const int& s) {
			return _ref->tt0[i][s];
		}
		inline double& get_hh0(const int& i, const int& s) {
			return _ref->hh0[i][s];
		}
		inline double& get_tt1(const int& i, const int& j, const int& s) {
			return _ref->tt1[(i << 1) + j][s];
		}
		inline double& get_hh1(const int& i, const int& j, const int& s) {
			return _ref->hh1[(i << 1) + j][s];
		}
		inline double& get_tt2(const int& i, const int& j, const int& k, const int& s) {
			return _ref->tt2[(((i << 1) + j) << 1) + k][s];
		}
		inline double& get_hh2(const int& i, const int& j, const int& k, const int& s) {
			return _ref->hh2[(((i << 1) + j) << 1) + k][s];
		}

	public:
		inline double _pow(const double& f, const int& k) {
			double val = 1;
			for (int i = 0; i < k; i++)
			{
				val *= f;
			}
			return val;
		}
		double __hh0(const int& j, const int& k) {
			double t = lo[j];
			return pow(t, (dim[j] - k - 1));
		}
		double __tt0(const int& j, const int& k) {
			double val = 0;
			double t = lo[j];
			for (int l = 0; l < dim[j]; l++)
			{
				val += get_hh0(j, l) * _ref->M[j][l][k];
			}
			return val;
		}

		double __hh1(const int& m, const int& j, const int& k) {
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
		double __tt1(const int& m, const int& j, const int& k) {
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
		_memS(string ultimate, _RAM RAM) {
			this->mode = ultimate;
			this->RAM = RAM;
			if (mode == "U") {
				__grad = 0;
				__grad_z = 0;
				__grad_phi = 0;
				__grad_z_tmp = 0;
				__grad_phi_tmp = 0;
				__grad_xi = 0;
				__grad_eta = 0;
				__mix_phi = 0;
				__mix_z = 0;
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

			if (RAM == SAVE)
			{
				__mat = 0;
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

		}
		~_memS() {
			del();
		}
		void update2() {
			if (!_ref->initialized || RAM == SAVE)
			{


				for (auto const& j : ___ee)
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
							_ref->hh1[(m << 1) + j][k] = __hh1(m, j, k);
						}
					}
				}
				double** ptr;
				ptr = _ref->tt1;
				for (auto const& m : ___ee)
				{
					for (auto const& j : ___ee)
					{
						for (int k = 0; k < dim[j]; k++) {
							*(*ptr + k) = __tt1(m, j, k);
							//tt1[m * 2 + j][k] = __tt1(m, j, k);
						}
						ptr++;
					}
				}
				ptr = _ref->hh2;
				for (auto const& n : ___ee)
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
				/*_ref->w = 0;
				double W = 0;
				for (int j = 0; j < _nNode; j++) {
					W += _ref->d0[j] * _ref->buf_W[j];
				}
				_ref->w = W;*/

				//for (int i = 0; i < _nNode; i++)
				//{
					//_ref->d0[i] *= _ref->buf_W[i] / _ref->w;
				//}


				ptr = _ref->d1;
				for (auto const& i : ___ee) {
					ptr2 = *ptr;
					for (int j = 0; j < _nNode; j++) {
						//d1[i][j] = _C(i, j);
						//*ptr2 = _ref->buf_W[j] * _C(i, j) / _ref->w;
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
							//*ptr2 = _ref->buf_W[j] * _D(i, ii, j) / _ref->w;
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
								//*ptr2 = _ref->buf_W[j] * _ref->buf_W[jj] * _B(i, ii, j, jj) / _ref->w / _ref->w;
								*ptr2 = _B(i, ii, j, jj);
								ptr2++;
							}
						}
						ptr++;
					}
				}
			}
			double X = 0, Y = 0, ZZ = 0, W = 0;
			double* ptr2 = _ref->d0;
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



			pptr1 = &_ref->buf_xi[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			this->xi = val;


			pptr1 = &_ref->buf_eta[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			this->eta = val;


			pptr1 = &_ref->buf_mu[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			this->mu = val;


			/*pptr1 = &_ref->buf_nu[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			this->nu = val;*/
			/*pptr1 = &_ref->buf_chi[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			this->chi = val;*/
			//covariant base vectors
			if (mode == "U" || mode == "SLOPE")
			{
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
						fz += *ptr2 * _ref->buf_z[i];
						ptr3++;
						ptr2++;
					}
					gi[j * 3 + 0] = fx;
					gi[j * 3 + 1] = fy;
					gi[j * 3 + 2] = fz;
				}
			}
			else {
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
			}
			//sc = 1 / _dv / _dv;

			if (mode == "U") {
				gij[0] = gi[0] * gi[0] + gi[1] * gi[1] + gi[2] * gi[2];
				gij[1] = gi[0] * gi[3] + gi[1] * gi[4] + gi[2] * gi[5];
				gij[2] = gij[1];
				gij[3] = gi[3] * gi[3] + gi[4] * gi[4] + gi[5] * gi[5];
				dv = sqrt(_det2(gij));
				gi2[0] = gi[0];
				gi2[1] = gi[1];
				gi2[2] = gi[2];
				gi2[3] = gi[3];
				gi2[4] = gi[4];
				gi2[5] = gi[5];
				
				_inv2(gij, Gij2);
				double Fx = 0, Fy = 0, Fz = 0;
				Fx += get_gi2(0, 0) * Gij2[0] + get_gi2(1, 0) * Gij2[1];
				Fy += get_gi2(0, 1) * Gij2[0] + get_gi2(1, 1) * Gij2[1];
				Fz += get_gi2(0, 2) * Gij2[0] + get_gi2(1, 2) * Gij2[1];
				Gi2[0] = Fx;
				Gi2[1] = Fy;
				Gi2[2] = Fz;
				Fx = 0; Fy = 0; Fz = 0;
				Fx += get_gi2(0, 0) * Gij2[2] + get_gi2(1, 0) * Gij2[3];
				Fy += get_gi2(0, 1) * Gij2[2] + get_gi2(1, 1) * Gij2[3];
				Fz += get_gi2(0, 2) * Gij2[2] + get_gi2(1, 2) * Gij2[3];
				Gi2[3] = Fx;
				Gi2[4] = Fy;
				Gi2[5] = Fz;


				gij2[0] = gij[0];
				gij2[1] = gij[1];
				gij2[2] = gij[2];
				gij2[3] = gij[3];

				gi[2] = 0;
				gi[5] = 0;

				gij[0] = gi[0] * gi[0] + gi[1] * gi[1];
				gij[1] = gi[0] * gi[3] + gi[1] * gi[4];
				gij[2] = gij[1];
				gij[3] = gi[3] * gi[3] + gi[4] * gi[4];
				_dv = sqrt(_det2(gij));

				sc = 1.0 / _det2(gij);
				_inv2(gij, Gij);
			}
			else {

				gij2[0] = gi[0] * gi[0] + gi[1] * gi[1];
				gij2[1] = gi[0] * gi[3] + gi[1] * gi[4];
				gij2[2] = gij2[1];
				gij2[3] = gi[3] * gi[3] + gi[4] * gi[4];
				_dv = sqrt(_det2(gij2));
				_inv2(gij2, Gij2);
				double Fx = 0, Fy = 0;
				Fx += get_gi(0, 0) * Gij2[0] + get_gi(1, 0) * Gij2[1];
				Fy += get_gi(0, 1) * Gij2[0] + get_gi(1, 1) * Gij2[1];
				Gi2[0] = Fx;
				Gi2[1] = Fy;
				Gi2[2] = 0;
				Fx = 0; Fy = 0;
				Fx += get_gi(0, 0) * Gij2[2] + get_gi(1, 0) * Gij2[3];
				Fy += get_gi(0, 1) * Gij2[2] + get_gi(1, 1) * Gij2[3];
				Gi2[3] = Fx;
				Gi2[4] = Fy;
				Gi2[5] = 0;

				gij[0] = gi[0] * gi[0] + gi[1] * gi[1] + gi[2] * gi[2];
				gij[1] = gi[0] * gi[3] + gi[1] * gi[4] + gi[2] * gi[5];
				gij[2] = gij[1];
				gij[3] = gi[3] * gi[3] + gi[4] * gi[4] + gi[5] * gi[5];
				dv = sqrt(_det2(gij));

				sc = 1.0 / _det2(gij);
				_inv2(gij, Gij);

			}
			//contravatiant base vectors
			double Fx = 0, Fy = 0, Fz = 0;
			
			Fx = 0, Fy = 0, Fz = 0;
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
			if (mode == "SLOPE")
			{
				double** ptr = _ref->d2;
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
						bij[e + 2] = 0;
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

							val += bij[eee + 0] * get_Gi2(k, 0);
							val += bij[eee + 1] * get_Gi2(k, 1);
							//val += bij[eee + 2] * get_Gi2(k, 2);
							Gammaijk[ccc] = val;
							ccc++;
						}
						eee += 3;
					}
				}
			}
			else {


				double** ptr = _ref->d2;
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
							val = 0;
							val += bij[eee + 0] * get_Gi2(k, 0);
							val += bij[eee + 1] * get_Gi2(k, 1);
							val += bij[eee + 2] * get_Gi2(k, 2);
							Gammaijk2[ccc] = val;
							ccc++;
						}
						eee += 3;
					}
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
				static const int sss[4]{ 0,1,2,3 };
				/*const int star2[4]{1,-1,-1,1};
				const int _star[4]{ 3,2,1,0 };
				*/
				if (!_ref->initialized || RAM == SAVE)
				{
					pptr = &_ref->__mat[0];
					double* _pt00 = &_ref->d1[0][0];
					double* _pt01 = &_ref->d1[1][0];

					/*for (int i = 0; i < _ref->_nNode; i++)
					{
						a += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
						b += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
						c += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
					}*/
					for (int i = 0; i < _nNode; i++)
					{
						double* _pt0 = &_ref->d1[0][0];
						double* _pt1 = &_ref->d1[1][0];
						for (int j = 0; j < _nNode; j++) {
							double val = 0;
							for (auto const& k : sss) {
								{
									int mm = (_star[k]) << 1;

									int nn = (k) << 1;

									val += star2[k] * (_ref->d2[_star[k]][i] - Gammaijk[mm] * (*_pt00) - Gammaijk[mm + 1] * (*_pt01)) *
										(_ref->d2[k][j] - Gammaijk[nn] * (*_pt0) - Gammaijk[nn + 1] * (*_pt1));
								}
							}
							*pptr = val;
							pptr++;
							_pt0++;
							_pt1++;
						}
						_pt00++;
						_pt01++;
					}
				}
				/*for (int i = 0; i < _nNode; i++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						if (j >= i)continue;
						_ref->__mat[i * _nNode + j] = (_ref->__mat[i * _nNode + j] + _ref->__mat[i * _nNode + j]) / 2.0;
						_ref->__mat[j * _nNode + i] = _ref->__mat[i * _nNode + j];
					}
				}*/
				pptr = &_ref->__mat[0];
				pptr1 = &_ref->buf_z[0];
				val3 = 0;
				for (int i = 0; i < _nNode; i++)
				{
					pptr2 = &_ref->buf_phi[0];
					for (int j = 0; j < _nNode; j++) {
						val3 += *pptr * (*pptr2) */*_ref->buf_z[j] **/ (*pptr1);//_ref->buf_phi[i];
						pptr++;
						pptr2++;
					}
					pptr1++;
				}
				bodyF = val3;
				pptr = &_ref->__mat[0];
				pptr2 = &__grad_phi[0];
				for (int i = 0; i < _nNode; i++) {
					val2 = 0;
					pptr1 = &_ref->buf_phi[0];
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
				pptr1 = &_ref->buf_z[0];
				pptr2 = &__grad_z[0];
				for (int j = 0; j < _nNode; j++) {
					*pptr2 = 0;
					pptr2++;
				}
				for (int i = 0; i < _nNode; i++) {
					//val2 = 0;
					pptr2 = &__grad_z[0];
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
			_ref->initialized = true;
			if (mode == "SLOPE") {
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
			if (mode == "U") {
				if (__grad_z != 0)delete[] __grad_z;
				if (__grad_phi != 0)delete[] __grad_phi;
				if (__grad_z_tmp != 0)delete[] __grad_z_tmp;
				if (__grad_phi_tmp != 0)delete[] __grad_phi_tmp;
				if (__grad_xi != 0)delete[] __grad_xi;
				if (__grad_eta != 0)delete[] __grad_eta;
				if (__mix_phi != 0)delete[] __mix_phi;
				if (__mix_z != 0)delete[] __mix_z;


				//if (__mat != 0)delete[] __mat;
				if (__grad_C_z != 0)delete[] __grad_C_z;
				if (__grad_C_phi != 0)delete[] __grad_C_phi;
				if (__grad_D_z != 0)delete[] __grad_D_z;
				if (__grad_D_phi != 0)delete[] __grad_D_phi;
				if (__grad != 0)delete[] __grad;
				if (__grad2 != 0)delete[] __grad2;
				if (_K != 0)delete[] _K;
				__grad = 0;
				__grad2 = 0;
				__grad_z = 0;
				__grad_phi = 0;
				__grad_z_tmp = 0;
				__grad_phi_tmp = 0;
				__grad_xi = 0;
				__grad_eta = 0;
				__mix_phi = 0;
				__mix_z = 0;
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
			if (RAM == SAVE)
			{
				if (tt0[0] != 0) {
					delete[] __mat;
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

					for (int i = 0; i < _nNode; i++) {
						delete[] dd[i];
					}
					delete[] dd;
					dd = 0;


				}
			}
		}
		void update_ref(int nNode, int uDim, int vDim) {
			if (!_ref->initialized)
			{
				_ref->update(nNode, uDim, vDim);
			}
		}
		void update(int nNode, int uDim, int vDim) {
			if (_nNode != nNode || _uDim != uDim || _vDim != vDim) {

				del();


				_nNode = nNode;
				_uDim = uDim;
				_vDim = vDim;
				dim[0] = _uDim;
				dim[1] = _vDim;
				if (mode == "U" || mode == "SLOPE")
				{
					__grad = new double[nNode];
					__grad2 = new double[nNode];
					__grad_z = new double[nNode];
					__grad_z_tmp = new double[nNode];
					__grad_phi_tmp = new double[nNode];
					__grad_xi = new double[nNode];
					__grad_eta = new double[nNode];
					__mix_z = new double [nNode];
					__mix_phi = new double[ nNode];

					__grad_phi = new double[nNode];
					//__mat = new double[nNode * nNode];
					//__mat2 = new double[nNode * nNode];
					__grad_C_phi = new double[nNode * 2];
					__grad_C_z = new double[nNode * 2];
					__grad_D_phi = new double[nNode * 2];
					__grad_D_z = new double[nNode * 2];
					_K = new double[nNode * 2];
				}
				_F3 = new double[nNode * 2];
				if (mode == "SENSITIVITY" || mode == "SHELL")
				{
					gradN[0] = new double[3 * _nNode];
					gradN[1] = new double[3 * _nNode];
					gradN[2] = new double[3 * _nNode];
					gradG = new double[8 * _nNode];
				}
				if (RAM == SAVE)
				{
					__mat = new double[nNode * nNode];
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
				/*
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
		void getMat_slope(int64_t* index, std::vector<_Triplet<double>>* _dat, double dcdt1, double dcdt2, bool airy)
		{
			_dat->clear();
			_dat->reserve(_nNode);
			if (airy)
			{
				int64_t* ptr2 = &index[0];
				for (int i = 0; i < _nNode; i++)
				{
					double val = 0;
					for (int k = 0; k < 2; k++)
					{
						val += _ref->d1[k][i] * get_Gij(k, 0) * dcdt1;
						val += _ref->d1[k][i] * get_Gij(k, 1) * dcdt2;
					}
					_dat->push_back(_Triplet<double>(*ptr2, -1, val));
					ptr2++;
				}
			}
			else {
				int64_t* ptr2 = &index[0];
				for (int i = 0; i < _nNode; i++)
				{
					double val = 0;
					for (int k = 0; k < 2; k++)
					{
						val += _ref->d1[k][i] * get_Gij(k, 0) * dcdt1;
						val += _ref->d1[k][i] * get_Gij(k, 1) * dcdt2;
					}
					_dat->push_back(_Triplet<double>(-1, *ptr2, val));
					ptr2++;
				}
			}
		}
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
		inline double SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {

			return sc * (_SLOPE_phi[0] * dcdtstar0 + _SLOPE_phi[1] * dcdtstar1);// val;
		}
		inline double SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {

			return sc * (_SLOPE_z[0] * dcdtstar0 + _SLOPE_z[1] * dcdtstar1);// val;
		}
	

		inline void constant_SLOPE(double* ptr, double dcdtstar0, double dcdtstar1, double sc, bool accurate)
		{

			double U1 = dcdtstar1;
			double U2 = -dcdtstar0;
			double v1 = dcdtstar0;
			double v2 = dcdtstar1;
			double* ptr1 = ptr;
			if (accurate)
			{
				for (int l = 0; l < _nNode; l++)
				{
					double H11 = 0, H12 = 0, H22 = 0;
					H11 += (_ref->d2[0][l] - _ref->_Gammaijk[0] * _ref->d1[0][l] - _ref->_Gammaijk[1] * _ref->d1[1][l]);
					H22 += (_ref->d2[3][l] - _ref->_Gammaijk[6] * _ref->d1[0][l] - _ref->_Gammaijk[7] * _ref->d1[1][l]);
					H12 += (_ref->d2[1][l] - _ref->_Gammaijk[2] * _ref->d1[0][l] - _ref->_Gammaijk[3] * _ref->d1[1][l]);
					double val = 0;

					val += H11 * U1 * (v1 * this->get_Gij(0, 0) + v2 * this->get_Gij(1, 0));
					val += H12 * U2 * (v1 * this->get_Gij(0, 0) + v2 * this->get_Gij(1, 0));
					val += H12 * U1 * (v1 * this->get_Gij(0, 1) + v2 * this->get_Gij(1, 1));
					val += H22 * U2 * (v1 * this->get_Gij(0, 1) + v2 * this->get_Gij(1, 1));
					*ptr1 = val;
					ptr1++;
				}
			}
			else {
				for (int l = 0; l < _nNode; l++)
				{
					double H11 = 0, H12 = 0, H22 = 0;
					H11 += (_ref->d2[0][l] - _ref->_Gammaijk[0] * _ref->d1[0][l] - _ref->_Gammaijk[1] * _ref->d1[1][l]);
					H22 += (_ref->d2[3][l] - _ref->_Gammaijk[6] * _ref->d1[0][l] - _ref->_Gammaijk[7] * _ref->d1[1][l]);
					H12 += (_ref->d2[1][l] - _ref->_Gammaijk[2] * _ref->d1[0][l] - _ref->_Gammaijk[3] * _ref->d1[1][l]);
					double val = 0;

					val += H11 * U1 * (v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0));
					val += H12 * U2 * (v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0));
					val += H12 * U1 * (v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1));
					val += H22 * U2 * (v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1));
					*ptr1 = val;
					ptr1++;
				}
			}
		}


		inline double constant_SLOPE_z(double dcdtstar0, double dcdtstar1, bool accurate) {

			double H11 = 0, H12 = 0, H22 = 0;
			double* pptr1 = &_ref->buf_z[0];

			double U1 = dcdtstar1;
			double U2 = -dcdtstar0;
			double v1 = dcdtstar0;
			double v2 = dcdtstar1;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				H11 += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				H22 += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				H12 += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			double val = 0;
			if (accurate)
			{
				val += H11 * U1 * (v1 * this->get_Gij(0, 0) + v2 * this->get_Gij(1, 0));
				val += H12 * U2 * (v1 * this->get_Gij(0, 0) + v2 * this->get_Gij(1, 0));
				val += H12 * U1 * (v1 * this->get_Gij(0, 1) + v2 * this->get_Gij(1, 1));
				val += H22 * U2 * (v1 * this->get_Gij(0, 1) + v2 * this->get_Gij(1, 1));
			}
			else {
				val += H11 * U1 * (v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0));
				val += H12 * U2 * (v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0));
				val += H12 * U1 * (v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1));
				val += H22 * U2 * (v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1));
			}
			return val;
		}
		inline double constant_SLOPE_phi(double dcdtstar0, double dcdtstar1, bool accurate) {

			double H11 = 0, H12 = 0, H22 = 0;
			double* pptr1 = &_ref->buf_phi[0];

			double U1 = dcdtstar1;
			double U2 = -dcdtstar0;
			double v1 = dcdtstar0;
			double v2 = dcdtstar1;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				H11 += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				H22 += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				H12 += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			double val = 0;
			if (accurate)
			{
				val += H11 * U1 * (v1 * this->get_Gij(0, 0) + v2 * this->get_Gij(1, 0));
				val += H12 * U2 * (v1 * this->get_Gij(0, 0) + v2 * this->get_Gij(1, 0));
				val += H12 * U1 * (v1 * this->get_Gij(0, 1) + v2 * this->get_Gij(1, 1));
				val += H22 * U2 * (v1 * this->get_Gij(0, 1) + v2 * this->get_Gij(1, 1));
			}
			else {
				val += H11 * U1 * (v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0));
				val += H12 * U2 * (v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0));
				val += H12 * U1 * (v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1));
				val += H22 * U2 * (v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1));
			}
			return val;
		}
		inline void shear(double* ptr, int uv)
		{
			double* ptr1 = ptr;
			if (uv == 0)
			{
				for (int i = 0; i < _ref->_nNode; i++)
				{
					*ptr1 = 0;
					*ptr1 += _ref->_Gij[2] * _ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i];
					//*ptr1 += _ref->_Gij[2] * _ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i];
					*ptr1 += _ref->_Gij[3] * _ref->d2[2][i] - _ref->_Gammaijk[4] * _ref->d1[0][i] - _ref->_Gammaijk[5] * _ref->d1[1][i];
					//*ptr1 += _ref->_Gij[3] * _ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i];
				}
			}
			else {
				for (int i = 0; i < _ref->_nNode; i++)
				{
					*ptr1 = 0;
					//*ptr1 += _ref->_Gij[0] * _ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i];
					*ptr1 += _ref->_Gij[0] * _ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i];
					//*ptr1 += _ref->_Gij[1] * _ref->d2[2][i] - _ref->_Gammaijk[4] * _ref->d1[0][i] - _ref->_Gammaijk[5] * _ref->d1[1][i];
					*ptr1 += _ref->_Gij[1] * _ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i];
				}
			}
		}
		inline double shear_z(int uv)
		{
			double val = 0;
			if (uv == 0)
			{
				for (int i = 0; i < _ref->_nNode; i++)
				{
					val += (_ref->_Gij[2] * _ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
					//val += (_ref->_Gij[2] * _ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
					val += (_ref->_Gij[3] * _ref->d2[2][i] - _ref->_Gammaijk[4] * _ref->d1[0][i] - _ref->_Gammaijk[5] * _ref->d1[1][i]) * _ref->buf_z[i];
					//val += (_ref->_Gij[3] * _ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
				}
			}
			else {
				for (int i = 0; i < _ref->_nNode; i++)
				{
					//val += (_ref->_Gij[0] * _ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
					val += (_ref->_Gij[0] * _ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
					//val += (_ref->_Gij[1] * _ref->d2[2][i] - _ref->_Gammaijk[4] * _ref->d1[0][i] - _ref->_Gammaijk[5] * _ref->d1[1][i]) * _ref->buf_z[i];
					val += (_ref->_Gij[1] * _ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
				}
			}
			return val;
		}
		inline double shear_phi(int uv)
		{
			double val = 0;
			if (uv == 0)
			{
				for (int i = 0; i < _ref->_nNode; i++)
				{
					val += (_ref->_Gij[2] * _ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
					//val += (_ref->_Gij[2] * _ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
					val += (_ref->_Gij[3] * _ref->d2[2][i] - _ref->_Gammaijk[4] * _ref->d1[0][i] - _ref->_Gammaijk[5] * _ref->d1[1][i]) * _ref->buf_z[i];
					//val += (_ref->_Gij[3] * _ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
				}
			}
			else
			{
				for (int i = 0; i < _ref->_nNode; i++)
				{
					//val += (_ref->_Gij[0] * _ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
					val += (_ref->_Gij[0] * _ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
					//val += (_ref->_Gij[1] * _ref->d2[2][i] - _ref->_Gammaijk[4] * _ref->d1[0][i] - _ref->_Gammaijk[5] * _ref->d1[1][i]) * _ref->buf_z[i];
					val += (_ref->_Gij[1] * _ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
				}
			}
			return val;
		}
		inline double SLOPE(int l, int i) {
			return _ref->d1[0][i] * get_Gij(0, l) + _ref->d1[1][i] * get_Gij(1, l);

		}
		inline void SLOPE(double* ptr, double dcdt1, double dcdt2, double sc) {
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = ((_ref->d1[0][i] * get_Gij(0, 0) + _ref->d1[1][i] * get_Gij(1, 0)) * dcdt1 +
					(_ref->d1[0][i] * get_Gij(0, 1) + _ref->d1[1][i] * get_Gij(1, 1)) * dcdt2) * sc;
				ptr1++;

			}

		}
		inline void _SLOPE(double* ptr, double dcdt1, double dcdt2,double sc) {
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = ((_ref->d1[0][i] * _ref->get__Gij(0, 0) + _ref->d1[1][i] * _ref->get__Gij(1, 0)) * dcdt1 +
					(_ref->d1[0][i] * _ref->get__Gij(0, 1) + _ref->d1[1][i] * _ref->get__Gij(1, 1)) * dcdt2)*sc;
				ptr1++;

			}

		}
		inline double SLOPE_xi(double dcdtstar0, double dcdtstar1) {

			double xi_u = 0, xi_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_xi[i];
				xi_v += _ref->d1[1][i] * _ref->buf_xi[i];
			}
			double f_U = xi_u * _ref->get__Gij(0, 0) + xi_v * _ref->get__Gij(1, 0);
			double f_V = xi_u * _ref->get__Gij(0, 1) + xi_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
		inline double SLOPE_eta(double dcdtstar0, double dcdtstar1) {

			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				eta_u += _ref->d1[0][i] * _ref->buf_eta[i];
				eta_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
			double f_U = eta_u * _ref->get__Gij(0, 0) + eta_v * _ref->get__Gij(1, 0);
			double f_V = eta_u * _ref->get__Gij(0, 1) + eta_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}

		inline double SLOPE_mu(double dcdtstar0, double dcdtstar1) {

			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_mu[i];
				f_v += _ref->d1[1][i] * _ref->buf_mu[i];
			}
			double f_U = f_u * this->get_Gij(0, 0) + f_v * this->get_Gij(1, 0);
			double f_V = f_u * this->get_Gij(0, 1) + f_v * this->get_Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}

		inline double Dc(double a, double b, int l) {
			double f1 = a * this->get_gi(0, 0) + b * this->get_gi(0, 1);
			double f2 = a * this->get_gi(1, 0) + b * this->get_gi(1, 1);
			return f1 * get_Gij(0, l) + f2 * get_Gij(1, l);
		}
		inline double __a(int l)
		{
			return this->get_gi(0, 0) * get_Gij(0, l) + this->get_gi(1, 0) * get_Gij(1, l);
		}
		inline double __b(int l)
		{
			return this->get_gi(0, 1) * get_Gij(0, l) + this->get_gi(1, 1) * get_Gij(1, l);
		}
		inline double __K(int l, int I) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				val += _ref->d1[k][I] * this->get_Gij(k, l);
			}
			return val;
		}
		inline double K(int l, int I)
		{
			return _K[l * _nNode + I];// val;
		}
		//this is essentially The distance between two points in the force diagram
		double __K_phi(int l)
		{
			double val = 0;
			for (int I = 0; I < _nNode; I++) {
				double val2 = 0;
				for (auto k : ___ee)
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
			return 0.5 * (_Sij[0] * _ref->B[0][e] + 2 * _Sij[1] * _ref->B[1][e] + _Sij[3] * _ref->B[3][e]) * _ref->refDv;
		}
		void F(_mySparse* mat, int64_t* index, double sc) {
			for (int i = 0; i < _nNode; i++)
			{
				int I = index[i];
				for (int j = 0; j < _nNode; j++)
				{
					int J = index[j];
					double val = F(i, j) * sc;
					mat->_mat[0].coeffRef(i * 3 + 0, j * 3 + 0) += val;
					mat->_mat[0].coeffRef(i * 3 + 0, j * 3 + 0) += val;
					mat->_mat[0].coeffRef(i * 3 + 0, j * 3 + 0) += val;
				}
			}
		}
		void getMat(int64_t* index, std::vector<_Triplet<double>>* _dat)
		{
			//mat->_mat[0].setZero();
			//mat->_mat[0].reserve(_nNode * _nNode);
			double* ptr = &_ref->__mat[0];
			int64_t* ptr2 = &index[0];
			_dat->clear();
			_dat->reserve(_nNode * _nNode);
			double norm = 0;
			for (int i = 0; i < _nNode; i++)
			{
				int64_t* ptr3 = &index[0];
				for (int j = 0; j < _nNode; j++)
				{
					//int I = index[i];
					//int J = index[j];
					//mat->_mat[0].coeffRef(*ptr2, *ptr3) = *ptr;// _ref->__mat[i * _nNode + j];
					_dat->push_back(_Triplet<double>(*ptr2, *ptr3, sc * _ref->__mat[i * _nNode + j]));
					norm += (sc * _ref->__mat[i * _nNode + j]) * (sc * _ref->__mat[i * _nNode + j]);
					//if (i == j)
					//	_dat->push_back(_Triplet<double>(*ptr2, *ptr3, _ref->d2[0][i]* _ref->d2[0][j]+ _ref->d2[3][i] * _ref->d2[3][j]));
					ptr++;
					ptr3++;
				}
				ptr2++;
			}
			norm = std::sqrt(norm);
			//mat->_mat[0].setFromTriplets(_dat->begin(), _dat->end());
		}
		void getMat_Galerkin(int64_t* index, std::vector<_Triplet<double>>* _dat, int I, double w)
		{

			int kk[4]{ 3,2,1,0 };
			int kk1[4]{ 1,1,0,0 };
			int kk2[4]{ 1,0,1,0 };

			int tt[4]{ 1,-1,-1,1 };
			int ss[4]{ 0,0,1,1, };
			int uu[4]{ 0,1,0,1, };
			int II = -1;
			for (int i = 0; i < _nNode; i++)
			{
				if (index[i] == I)II = i;
			}
			if (II == -1)return;

			_dat->reserve(_dat->size() + _nNode * _nNode);
			for (int i = 0; i < _nNode; i++)
			{
				for (int j = 0; j < _nNode; j++)
				{
					double val = 0;
					for (int k = 0; k < 4; k++)
					{
						val += w * _ref->refDv * sc * tt[k] * (_ref->d2[kk[k]][i] - _ref->get__Gammaijk(kk1[k], kk2[k], 0) * _ref->d1[0][i] - _ref->get__Gammaijk(kk1[k], kk2[k], 1) * _ref->d1[1][i]) * (_ref->d1[ss[k]][j] * _ref->d1[uu[k]][II] + _ref->d1[ss[k]][II] * _ref->d1[uu[k]][j]) * 0.5;
					}
					_dat->push_back(_Triplet<double>(index[i], index[j], val));
				}

			}

			//mat->_mat[0].setZero();
			//mat->_mat[0].reserve(_nNode * _nNode);
			/*double* ptr = &_ref->__mat[0];
			int64_t* ptr2 = &index[0];
			_dat->clear();
			_dat->reserve(_nNode * _nNode);
			double norm = 0;
			for (int i = 0; i < _nNode; i++)
			{
				int64_t* ptr3 = &index[0];
				for (int j = 0; j < _nNode; j++)
				{
					//int I = index[i];
					//int J = index[j];
					//mat->_mat[0].coeffRef(*ptr2, *ptr3) = *ptr;// _ref->__mat[i * _nNode + j];
					_dat->push_back(_Triplet<double>(*ptr2, *ptr3, sc * _ref->__mat[i * _nNode + j]));
					norm += (sc * _ref->__mat[i * _nNode + j]) * (sc * _ref->__mat[i * _nNode + j]);
					//if (i == j)
					//	_dat->push_back(_Triplet<double>(*ptr2, *ptr3, _ref->d2[0][i]* _ref->d2[0][j]+ _ref->d2[3][i] * _ref->d2[3][j]));
					ptr++;
					ptr3++;
				}
				ptr2++;
			}
			norm = std::sqrt(norm);*/
			//mat->_mat[0].setFromTriplets(_dat->begin(), _dat->end());
		}
		void getMatG(int64_t* index, int64_t* index2, std::vector<_Triplet<double>>* _dat, double sc)
		{
			//double* ptr = &_ref->__mat[0];
			int64_t* ptr2 = &index[0];
			//_dat->clear();
			//_dat->reserve(_nNode * _nNode);
			for (int i = 0; i < _nNode; i++)
			{
				int64_t* ptr3 = &index2[0];
				for (int j = 0; j < _nNode; j++)
				{
					//int I = index[i];
					//int J = index[j];
					//mat->_mat[0].coeffRef(*ptr2, *ptr3) = *ptr;// _ref->__mat[i * _nNode + j];
					//_dat->push_back(_Triplet<double>(*ptr2, *ptr3, sc * _ref->__mat[i * _nNode + j]));
					//if (i == j)
					_dat->push_back(_Triplet<double>(*ptr2, *ptr3, (
						(_ref->d2[3][i] - _ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][i] - _ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][i]) * (_ref->d2[0][j] - _ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][j] - _ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][j]) +
						(_ref->d2[0][i] - _ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][i] - _ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][i]) * (_ref->d2[3][j] - _ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][j] - _ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][j]) -
						(_ref->d2[1][i] - _ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - _ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) * (_ref->d2[2][j] - _ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][j] - _ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][j]) -
						(_ref->d2[2][i] - _ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - _ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) * (_ref->d2[1][j] - _ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][j] - _ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][j])) * sc));

					//ptr++;
					ptr3++;
				}
				ptr2++;
			}
			//mat->_mat[0].setFromTriplets(_dat->begin(), _dat->end());
		}
		void getMat_d0(int64_t* index, std::vector<_Triplet<double>>* _dat)
		{
			//mat->_mat[0].setZero();
			//mat->_mat[0].reserve(_nNode * _nNode);
			double* ptr = &_ref->d0[0];
			int64_t* ptr2 = &index[0];
			_dat->clear();
			_dat->reserve(_nNode * _nNode);
			for (int i = 0; i < _nNode; i++)
			{
				_dat->push_back(_Triplet<double>(-1, *ptr2, *ptr));
				ptr++;
				ptr2++;
			}
			//mat->_mat[0].setFromTriplets(_dat->begin(), _dat->end());
		}
		double detZ()
		{
			double* pptr = &_ref->__mat[0];
			double* pptr1 = &_ref->buf_z[0];
			double val3 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				double* pptr2 = &_ref->buf_z[0];
				for (int j = 0; j < _ref->_nNode; j++) {
					val3 += *pptr * (*pptr1) */*_ref->buf_z[j] **/ (*pptr2);//_ref->buf_phi[i];
					pptr++;
					pptr2++;
				}
				pptr1++;
			}
			return val3;
		}

		double _detZ()
		{

			double* pptr1 = &_ref->buf_z[0];
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				a += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				c += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				b += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			return a * c - b * b;
		}
		double _detZ_a()
		{

			double* pptr1 = &_ref->buf_z[0];
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				a += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				c += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				b += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			return a;
		}
		double _detZ_b()
		{

			double* pptr1 = &_ref->buf_z[0];
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				a += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				c += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				b += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			return b;
		}
		double _detZ_c()
		{

			double* pptr1 = &_ref->buf_z[0];
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				a += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				c += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				b += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			return c;
		}
		double _detphi()
		{

			double* pptr1 = &_ref->buf_phi[0];
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				a += pptr1[i] * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]);
				c += pptr1[i] * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]);
				b += pptr1[i] * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]);
			}
			return a * c - b * b;
		}
		double detphi()
		{
			double* pptr = &_ref->__mat[0];
			double* pptr1 = &_ref->buf_phi[0];
			double val3 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				double* pptr2 = &_ref->buf_phi[0];
				for (int j = 0; j < _ref->_nNode; j++) {
					val3 += *pptr * (*pptr1) */*_ref->buf_z[j] **/ (*pptr2);//_ref->buf_phi[i];
					pptr++;
					pptr2++;
				}
				pptr1++;
			}
			return val3;
		}
		//stress function L2
		double F2(int i, int j) {
			return
				//return  (_Sij[0] * (d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * d1[1][i]) + 2 * _Sij[1] * (d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * d1[1][i]) + _Sij[3] * (d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * d1[1][i]))*
				(_Sij[0] * (_ref->d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][i]) + 2 * _Sij[1] * (_ref->d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) + _Sij[3] * (_ref->d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][i]))
				* (_Sij[0] * (_ref->d2[0][j] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][j] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][j]) + 2 * _Sij[1] * (_ref->d2[1][j] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][j] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][j]) + _Sij[3] * (_ref->d2[3][j] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][j] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][j])) * _ref->refDv;// / _ref->refDv / _ref->refDv;

		}
		std::string F2(_mySparse* mat, int64_t* index, double sc) {
			Eigen::initParallel();
			std::stringstream ss;
			auto start = high_resolution_clock::now();

			for (int i = 0; i < _nNode; i++)
			{
				int I = index[i];
				for (int j = 0; j < _nNode; j++)
				{
					int J = index[j];
					mat->_mat[0].coeffRef(I, J) += F2(i, j) * sc;
				}
			}
			return ss.str();
		}
		double F2A(int i, int j) {
			return
				//return  (_Sij[0] * (d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * d1[1][i]) + 2 * _Sij[1] * (d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * d1[1][i]) + _Sij[3] * (d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * d1[1][i]))*
				(_Sij[0] * (_ref->d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][i]) + 2 * _Sij[1] * (_ref->d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) + _Sij[3] * (_ref->d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][i]))
				* (_ref->d2[0][j]) * _ref->refDv;// / _ref->refDv / _ref->refDv;
		}
		std::string F2A(_mySparse* mat, int64_t* index, double sc) {
			Eigen::initParallel();
			std::stringstream ss;
			auto start = high_resolution_clock::now();

			for (int i = 0; i < _nNode; i++)
			{
				int I = index[i];
				for (int j = 0; j < _nNode; j++)
				{
					int J = index[j];
					mat->_mat[0].coeffRef(I, J) += F2A(i, j) * sc;
				}
			}
			return ss.str();
		}
		double __F3(int i, int I) {
			double val = 0;
			for (int ii = 0; ii < 2; ii++)
			{
				val += _Sij[(ii << 1) + i] * _ref->d1[ii][I];
			}
			return val;// *_ref->refDv;
		}
		double F3(int i, int I) {
			double val = 0;
			for (int ii = 0; ii < 2; ii++)
			{
				val += _Sij[(ii << 1) + i] * _ref->d1[ii][I];
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
					val2 += _Sij[(ii << 1) + i] * _ref->d1[ii][I];
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
				val += _F3[i * _nNode + I]/* val2*/ * _ref->buf_phi[I];
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
					val2 += _Sij[(ii << 1) + i] * _ref->d1[ii][I];
				}
				val += val2 * _ref->buf_phi[I];
			}
			return val;// _F3_phi[i];// val;// *_ref->refDv;
		}
		double F4(int i, int I, int J) {
			double val = 0;
			double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			for (int ii = 0; ii < 2; ii++)
			{
				int ij = ii * 2 + i;
				val += star2[ij] * (_ref->d2[_star[ij]][I] - Gammaijk[(_star[ij]) * 2 + 0] * _ref->d1[0][I] - Gammaijk[(_star[ij]) * 2 + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
			}
			return sc * val;
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
						val += (_ref->d2[ij][I] - Gammaijk[((ij) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((ij) << 1) + 1] * _ref->d1[1][I]) * this->get_Gij(ii, k) * this->get_Gij(i, l) * _ref->d1[l][J];
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
							val2 += (_ref->d2[ij][I] - Gammaijk[((ij) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((ij) << 1) + 1] * _ref->d1[1][I]) * this->get_Gij(ii, k) * this->get_Gij(i, l) * _ref->d1[l][J];
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
							val += (_ref->d2[ij][I] - Gammaijk[((ij) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((ij) << 1) + 1] * _ref->d1[1][I]) * this->get_Gij(ii, k) * this->get_Gij(i, l) * _ref->d1[l][J];
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
					val += star2[ij] * (_ref->d2[_star[ij]][I] - Gammaijk[((_star[ij]) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((_star[ij]) << 1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
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
					int ij = (ii << 1) + i;
					val += star2[ij] * (_ref->d2[_star[ij]][I] - Gammaijk[((_star[ij]) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((_star[ij]) << 1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
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
					int ij = (ii << 1) + i;
					val += star2[ij] * (_ref->d2[_star[ij]][I] - Gammaijk[((_star[ij]) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((_star[ij]) << 1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
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
					int ij = (ii << 1) + i;
					val += star2[ij] * (_ref->d2[_star[ij]][I] - Gammaijk[((_star[ij]) << 1) + 0] * _ref->d1[0][I] - Gammaijk[((_star[ij]) << 1) + 1] * _ref->d1[1][I]) * _ref->d1[ii][J];
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

			return sc * __grad_C_phi[i * _nNode + J];
		}
		void F4_z(double* ptr, double dcdtstar0, double dcdtstar1, double sc2) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			double* ptr1 = ptr;
			double _sc1 = sc * sc2 * dcdtstar0;
			double _sc2 = sc * sc2 * dcdtstar1;

			double* ptr2 = __grad_C_z;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = _sc1 * (*ptr2) + _sc2 * (*(ptr2 + _nNode));
				ptr1++;
				ptr2++;
			}

			//return sc * __grad_C_z[i * _nNode + I];
		}
		void F4_phi(double* ptr, double dcdtstar0, double dcdtstar1, double sc2) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;
			double* ptr1 = ptr;
			double _sc1 = sc * sc2 * dcdtstar0;
			double _sc2 = sc * sc2 * dcdtstar1;

			double* ptr2 = __grad_C_phi;
			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = _sc1 * (*ptr2) + _sc2 * (*(ptr2 + _nNode));
				ptr1++;
				ptr2++;
			}

			//return sc * __grad_C_z[i * _nNode + I];
		}
		double F5_z(int i, int I) {
			//double sc = 1 / this->_ref->_refDv / this->_ref->_refDv;

			return sc * __grad_D_z[i * _nNode + I];
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
		double _d2(int l, int i) {
			return _ref->d2[l][i] - _ref->_Gammaijk[l * 2 + 0] * _ref->d1[0][i] - _ref->_Gammaijk[l * 2 + 1] * _ref->d1[1][i];
		}
		double G3(int i) {
			//return d0[i] * _ref->refDv;
			return (_Sij[0] * (_ref->d2[0][i] - this->_ref->get__Gammaijk(0, 0, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 0, 1) * _ref->d1[1][i]) + 2 * _Sij[1] * (_ref->d2[1][i] - this->_ref->get__Gammaijk(0, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(0, 1, 1) * _ref->d1[1][i]) + _Sij[3] * (_ref->d2[3][i] - this->_ref->get__Gammaijk(1, 1, 0) * _ref->d1[0][i] - this->_ref->get__Gammaijk(1, 1, 1) * _ref->d1[1][i])) * _ref->refDv;// / _ref->refDv / _ref->refDv;

		}
		void G3(_myDoubleArray* vec, int64_t* index, double sc) {
			for (int i = 0; i < _nNode; i++)
			{
				int I = index[i];
				vec->__v.coeffRef(I) += G3(i) * sc;
			}

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
		void MASS(_mySparse* mat, int64_t* index, double sc) {
			for (int i = 0; i < _nNode; i++)
			{
				int I = index[i];
				for (int j = 0; j < _nNode; j++)
				{
					int J = index[j];
					//(*dat)[i + j * _nNode] = Eigen::Triplet<double>(I, J, MASS(i, j) * sc);
					mat->_mat[0].coeffRef(I, J) += MASS(i, j) * sc;
				}
			}
		}
		double SMOOTH(int I, int J) {
			double val = 0;
			for (int i = 0; i < 2; i++) {
				for (int j = 0; j < 2; j++)
				{
					for (int k = 0; k < 2; k++) {
						for (int l = 0; l < 2; l++)
						{
							val += (_ref->d2[i * 2 + j][I] - this->get_Gammaijk(i, j, 0) * _ref->d1[0][I] - this->get_Gammaijk(i, j, 1) * _ref->d1[1][I]) * this->get_Gij(i, k) * this->get_Gij(j, l) * (_ref->d2[k * 2 + l][J] - this->get_Gammaijk(k, l, 0) * _ref->d1[0][J] - this->get_Gammaijk(k, l, 1) * _ref->d1[1][J]);
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

		void fM(double _la, double _mu, double* a, double* b, double* c, double* d)
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
					for (int m = 0; m < 2; m++)
					{
						for (int n = 0; n < 2; n++)
						{


							double A = _la * _ref->get__Gij(l, k) * _ref->get__Gij(n, m) + 2 * _mu * _ref->get__Gij(l, n) * _ref->get__Gij(k, m);
							//A = 1.0;

							double D = 0;
							//double E = 0;
							for (int s = 0; s < 3; s++)
							{
								auto ff = get_gi(n, s);
								auto gg = _ref->get__gi(n, s);

								D += _ref->get__gi(m, s) * (get_gi(n, s) - _ref->get__gi(n, s));
								D += _ref->get__gi(n, s) * (get_gi(m, s) - _ref->get__gi(m, s));
							}
							//double D2 = /*get_gij(n, m) - */ _ref->get__gij(n, m);
							val += A * D;

							//S[(m << 1) + n] = D;

						}
					}
					S[(k << 1) + l] = val;// get_gij(k, l);// -_ref->get__gij(k, l);
				}
			}
			double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
			double gg11 = _ref->get__gi(0, 0) * _ref->get__gi(0, 0) + _ref->get__gi(0, 1) * _ref->get__gi(0, 1);
			double gg21 = _ref->get__gi(0, 0) * _ref->get__gi(1, 0) + _ref->get__gi(0, 1) * _ref->get__gi(1, 1);
			double gg12 = _ref->get__gi(0, 0) * _ref->get__gi(1, 0) + _ref->get__gi(0, 1) * _ref->get__gi(1, 1);
			double gg22 = _ref->get__gi(1, 0) * _ref->get__gi(1, 0) + _ref->get__gi(1, 1) * _ref->get__gi(1, 1);
			double det2 = gg11 * gg22 - gg12 * gg21;
			double J = sqrt(det / det2);
		
			*a = S[0] * _ref->get__gi(0, 0) * _ref->get__gi(0, 0) + S[1] * _ref->get__gi(0, 0) * _ref->get__gi(1, 0) + S[2] * _ref->get__gi(1, 0) * _ref->get__gi(0, 0) + S[3] * _ref->get__gi(1, 0) * _ref->get__gi(1, 0);
			*b = S[0] * _ref->get__gi(0, 0) * _ref->get__gi(0, 1) + S[1] * _ref->get__gi(0, 0) * _ref->get__gi(1, 1) + S[2] * _ref->get__gi(1, 0) * _ref->get__gi(0, 1) + S[3] * _ref->get__gi(1, 0) * _ref->get__gi(1, 1);
			*c = S[0] * _ref->get__gi(0, 1) * _ref->get__gi(0, 0) + S[1] * _ref->get__gi(0, 1) * _ref->get__gi(1, 0) + S[2] * _ref->get__gi(1, 1) * _ref->get__gi(0, 0) + S[3] * _ref->get__gi(1, 1) * _ref->get__gi(1, 0);
			*d = S[0] * _ref->get__gi(0, 1) * _ref->get__gi(0, 1) + S[1] * _ref->get__gi(0, 1) * _ref->get__gi(1, 1) + S[2] * _ref->get__gi(1, 1) * _ref->get__gi(0, 1) + S[3] * _ref->get__gi(1, 1) * _ref->get__gi(1, 1);
			*a *= J;
			*b *= J;
			*c *= J;
			*d *= J;


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


							double A = _la * _ref->get__Gij(l, k) * _ref->get__Gij(n, m) + 2 * _mu * _ref->get__Gij(l, n) * _ref->get__Gij(k, m);


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

		double eK(double _la, double _mu) {
			double bending = 0;
			for (int k = 0; k < 2; k++) {
				for (int l = 0; l < 2; l++) {
					double __bij2 = get_bij(k, l, 0) * N[0] + get_bij(k, l, 1) * N[1] + get_bij(k, l, 2) * N[2];
					double ___bij2 = _ref->get__bij(k, l, 0) * N[0] + _ref->get__bij(k, l, 1) * N[1] + _ref->get__bij(k, l, 2) * N[2];
					for (int m = 0; m < 2; m++) {
						for (int n = 0; n < 2; n++) {
							double __bij = get_bij(m, n, 0) * N[0] + get_bij(m, n, 1) * N[1] + get_bij(m, n, 2) * N[2];
							double ___bij = _ref->get__bij(m, n, 0) * N[0] + _ref->get__bij(m, n, 1) * N[1] + _ref->get__bij(m, n, 2) * N[2];

							double A = _la * _ref->get__Gij(l, k) * _ref->get__Gij(n, m) + 2 * _mu * _ref->get__Gij(l, n) * _ref->get__Gij(k, m);

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
							val = _ref->d2[i * 2 + j][A] * get_Gi(k, 2);
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
							gradG[((i * 2 + j) * 2 + k) * _nNode + A] = val;
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
		double _M(int alpha, double _la, double _mu)
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



							double ijkl2 = _la * get_Gij(i, j) * get_Gij(k, l) + 2 * _mu * (get_Gij(i, k) * get_Gij(j, l));

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
							double ijkl = _la * get_Gij(i, j) * get_Gij(k, l) + 2 * _mu * get_Gij(i, k) * get_Gij(j, l);

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
		double K(int i, int k2, int j, int k, double _la, double _mu)
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
							double D = _ref->d2[g * 2 + h][j]/**N[k]+_ref->get__bij(g,h,k)*gradN[k][j];//*/ - Gammaijk[(g * 2 + h) * 2 + 0] * _ref->d1[0][j] - Gammaijk[(g * 2 + h) * 2 + 1] * _ref->d1[1][j];
							double E = _ref->d2[l * 2 + m][i] /** N[k2] + _ref->get__bij(l, m, k2) * gradN[k2][i];//*/ - Gammaijk[(l * 2 + m) * 2 + 0] * _ref->d1[0][i] - Gammaijk[(l * 2 + m) * 2 + 1] * _ref->d1[1][i];
							_val3 += A * (D)*N[k] * N[k2] * (E);
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
							double D = _ref->get__bij(g, h, k);
							double E = (_ref->d2[l * 2 + m][i] - Gammaijk[(l * 2 + m) * 2 + 0] * _ref->d1[0][i] - Gammaijk[(l * 2 + m) * 2 + 1] * _ref->d1[1][i]);
							_val3 += A * N[k] * (D) * (E);
						}
					}
				}
			}
			return (_val3)*this->dv;
		}
		void K(_mySparse* M, int64_t* _index, double _la, double _mu, double __sc)
		{
			const static int kk[3]{ 0,1,2 };
			const static int ll[2]{ 0,1 };
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
							_val3 = K(i, k, j, k2, _la, _mu);

							M->_mat[0].coeffRef(I + k, J + k2) += _val3 * __sc;
						}
					}
				}
			}
		}
		//ultimate term
		int star(int i) {
			if (i == 0)return 3;
			if (i == 1)return 2;
			if (i == 2)return 1;
			if (i == 3)return 0;
			return 0;
		}
		double Gamma000()
		{
			return Gammaijk[0];
		}
		double Gamma001()
		{
			return Gammaijk[1];
		}
		double Gamma010()
		{
			return Gammaijk[2];
		}
		double Gamma011()
		{
			return Gammaijk[3];
		}
		double Gamma110()
		{
			return Gammaijk[6];
		}
		double Gamma111()
		{
			return Gammaijk[7];
		}
		double f2uu()
		{
			double val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				val += _ref->d2[0][i] * _ref->buf_z[i];
			}
			return val;

		}
		double f2uv()
		{
			double val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				val += _ref->d2[1][i] * _ref->buf_z[i];
			}
			return val;

		}
		double f2vv()
		{
			double val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				val += _ref->d2[3][i] * _ref->buf_z[i];
			}
			return val;
		}
		double f1u()
		{
			double val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				val += _ref->d1[0][i] * _ref->buf_z[i];
			}
			return val;
		}
		double f1v()
		{
			double val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				val += _ref->d1[1][i] * _ref->buf_z[i];
			}
			return val;
		}

		void conjugate(double v1, double v2, double* s1, double* s2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;


			*s1 = (S12 * v1 + S22 * v2);
			*s2 = -(S11 * v1 + S12 * v2);
		}
		void transpose(double v1, double v2, double* s1, double* s2, bool accurate)
		{

			if (accurate)
			{
				*s1 = (this->get_Gij(0, 0) * v2 - this->get_Gij(0, 1) * v1);
				*s2 = (this->get_Gij(1, 0) * v2 - this->get_Gij(1, 1) * v1);
			}
			else {
				*s1 = (_ref->get__Gij(0, 0) * v2 - _ref->get__Gij(0, 1) * v1);
				*s2 = (_ref->get__Gij(1, 0) * v2 - _ref->get__Gij(1, 1) * v1);
			}
		}double st(double v1, double v2, bool accurate)
		{
			double length = 0;

			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
			}

			double S21 = S12;
			double val = 0;


			val = (S11 * v1 * v1 * this->get_gij(0, 1) + S11 * v1 * v2 * this->get_gij(1, 1) + S12 * v2 * v1 * this->get_gij(0, 1) + S12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
			val += -(S21 * v1 * v1 * this->get_gij(0, 0) + S21 * v1 * v2 * this->get_gij(1, 0) + S22 * v2 * v1 * this->get_gij(0, 0) + S22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;

			return val;
		}


		void st_z(double* ptr, double v1, double v2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double length = 0;

			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}

			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
			}

			double _S21 = _S12;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				double S21 = S12;

	

				val = (S11 * v1 * v1 * this->get_gij(0, 1) + S11 * v1 * v2 * this->get_gij(1, 1) + S12 * v2 * v1 * this->get_gij(0, 1) + S12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
				val += -(S21 * v1 * v1 * this->get_gij(0, 0) + S21 * v1 * v2 * this->get_gij(1, 0) + S22 * v2 * v1 * this->get_gij(0, 0) + S22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;
				

				*ptr1 = val;
				ptr1++;
			}
		}
		void st_phi(double* ptr, double v1, double v2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double length = 0;

			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}

			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
			}

			double _S21 = _S12;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				double S21 = S12;

				double g11 = 0, g12 = 0, g21 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				g21 = g12;

				val = (_S11 * v1 * v1 * g12 + _S11 * v1 * v2 * g22 + _S12 * v2 * v1 * g12 + _S12 * v2 * v2 * g22) / _ref->_refDv;
				val += -(_S21 * v1 * v1 * g11 + _S21 * v1 * v2 * g21 + _S22 * v2 * v1 * g11 + _S22 * v2 * v2 * g21) / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}
		}

		double st2(double v1, double v2, bool accurate)
		{
			double length = 0;

			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;
				
				v1 = v1 / length;
				v2 = v2 / length;
				
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			
			double S21 = S12;
			double val = 0;


			val =  (S11 * v1 * v1 * this->get_gij(0, 1) + S11 * v1 * v2 * this->get_gij(1, 1) + S12 * v2 * v1 * this->get_gij(0, 1) + S12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
			val += -(S21 * v1 * v1 * this->get_gij(0, 0) + S21 * v1 * v2 * this->get_gij(1, 0) + S22 * v2 * v1 * this->get_gij(0, 0) + S22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;

			return val;
		}


		void st2_phi(double* ptr, double v1, double v2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double length = 0;

			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}

			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}

			double _S21 = _S12;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				double S21 = S12;

				double g11 = 0, g12 = 0, g21 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				g21 = g12;

				val = (S11 * v1 * v1 * this->get_gij(0, 1) + S11 * v1 * v2 * this->get_gij(1, 1) + S12 * v2 * v1 * this->get_gij(0, 1) + S12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
				val += -(S21 * v1 * v1 * this->get_gij(0, 0) + S21 * v1 * v2 * this->get_gij(1, 0) + S22 * v2 * v1 * this->get_gij(0, 0) + S22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;
				//val += (_S11 * v1 * v1 * g12 + _S11 * v1 * v2 * g22 + _S12 * v2 * v1 * g12 + _S12 * v2 * v2 * g22) / _ref->_refDv;
				//val += -(_S21 * v1 * v1 * g11 + _S21 * v1 * v2 * g21 + _S22 * v2 * v1 * g11 + _S22 * v2 * v2 * g21) / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}
		}





		double tt(double v1, double v2, bool accurate)
		{
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
			}
			double s21 = s12;
			double val = 0;
			val += v1 * s11 * v1;
			val += v1 * s12 * v2;
			val += v2 * s21 * v1;
			val += v2 * s22 * v2;
			double norm = 0;
			if (accurate)
			{
				norm += v1 * this->get_gij(0, 0) * v1;
				norm += v1 * this->get_gij(0, 1) * v2;
				norm += v2 * this->get_gij(1, 0) * v1;
				norm += v2 * this->get_gij(1, 1) * v2;
			}
			else
			{
				norm += v1 * _ref->get__gij(0, 0) * v1;
				norm += v1 * _ref->get__gij(0, 1) * v2;
				norm += v2 * _ref->get__gij(1, 0) * v1;
				norm += v2 * _ref->get__gij(1, 1) * v2;
			}
			return val/ norm;
		}
		void tt_z(double* ptr, double v1, double v2, bool accurate)
		{
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double norm = 0;
			if (accurate)
			{
				norm += v1 * this->get_gij(0, 0) * v1;
				norm += v1 * this->get_gij(0, 1) * v2;
				norm += v2 * this->get_gij(1, 0) * v1;
				norm += v2 * this->get_gij(1, 1) * v2;
			}
			else
			{
				norm += v1 * _ref->get__gij(0, 0) * v1;
				norm += v1 * _ref->get__gij(0, 1) * v2;
				norm += v2 * _ref->get__gij(1, 0) * v1;
				norm += v2 * _ref->get__gij(1, 1) * v2;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				s11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				
				
				double s21 = s12;
				val += v1 * s11 * v1;
				val += v1 * s12 * v2;
				val += v2 * s21 * v1;
				val += v2 * s22 * v2;
				*ptr1 = val/ norm;
				ptr1++;
			}
		}

		

		double tt_2(double v1, double v2, bool accurate)
		{
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double s21 = s12;
			double val = 0;
			val += v1 * s11 * v1;
			val += v1 * s12 * v2;
			val += v2 * s21 * v1;
			val += v2 * s22 * v2;
			double norm = 0;
			if (accurate)
			{
				norm += v1 * this->get_gij(0, 0) * v1;
				norm += v1 * this->get_gij(0, 1) * v2;
				norm += v2 * this->get_gij(1, 0) * v1;
				norm += v2 * this->get_gij(1, 1) * v2;
			}
			else
			{
				norm += v1 * _ref->get__gij(0, 0) * v1;
				norm += v1 * _ref->get__gij(0, 1) * v2;
				norm += v2 * _ref->get__gij(1, 0) * v1;
				norm += v2 * _ref->get__gij(1, 1) * v2;
			}
			return val / norm;
		}
		void tt2_phi(double* ptr, double v1, double v2, bool accurate)
		{
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double norm = 0;
			if (accurate)
			{
				norm += v1 * this->get_gij(0, 0) * v1;
				norm += v1 * this->get_gij(0, 1) * v2;
				norm += v2 * this->get_gij(1, 0) * v1;
				norm += v2 * this->get_gij(1, 1) * v2;
			}
			else
			{
				norm += v1 * _ref->get__gij(0, 0) * v1;
				norm += v1 * _ref->get__gij(0, 1) * v2;
				norm += v2 * _ref->get__gij(1, 0) * v1;
				norm += v2 * _ref->get__gij(1, 1) * v2;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				s11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);


				double s21 = s12;
				val += v1 * s11 * v1;
				val += v1 * s12 * v2;
				val += v2 * s21 * v1;
				val += v2 * s22 * v2;
				*ptr1 = val / norm;
				ptr1++;
			}
		}
		double curvature(double v1, double v2, bool accurate)
		{
			double val = 0;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}


			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
			double sdet = sqrt(det);
			//val += (_ref->get__Gammaijk(0, 0, 0) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 0) * v1 * v2 + _ref->get__Gammaijk(1, 1, 0) * v2 * v2) * S1;
			//val += (_ref->get__Gammaijk(0, 0, 1) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 1) * v1 * v2 + _ref->get__Gammaijk(1, 1, 1) * v2 * v2) * S2;
			double gamma111=0, gamma112=0, gamma121=0, gamma122=0, gamma221=0, gamma222=0;

			double _gamma111 = _ref->get__Gammaijk(0, 0, 0) * this->get_gij2(0, 0) + _ref->get__Gammaijk(0, 0, 1) * this->get_gij2(1, 0);
			double _gamma112 = _ref->get__Gammaijk(0, 0, 0) * this->get_gij2(0, 1) + _ref->get__Gammaijk(0, 0, 1) * this->get_gij2(1, 1);
			double _gamma121 = _ref->get__Gammaijk(0, 1, 0) * this->get_gij2(0, 0) + _ref->get__Gammaijk(0, 1, 1) * this->get_gij2(1, 0);
			double _gamma122 = _ref->get__Gammaijk(0, 1, 0) * this->get_gij2(0, 1) + _ref->get__Gammaijk(0, 1, 1) * this->get_gij2(1, 1);
			double _gamma221 = _ref->get__Gammaijk(1, 1, 0) * this->get_gij2(0, 0) + _ref->get__Gammaijk(1, 1, 1) * this->get_gij2(1, 0);
			double _gamma222 = _ref->get__Gammaijk(1, 1, 0) * this->get_gij2(0, 1) + _ref->get__Gammaijk(1, 1, 1) * this->get_gij2(1, 1);

			val = (_gamma111 * v1 * v1 + 2 * _gamma121 * v1 * v2 + _gamma221 * v2 * v2) * s1/sdet;
			val += (_gamma112 * v1 * v1 + 2 * _gamma122 * v1 * v2 + _gamma222 * v2 * v2) * s2 / sdet;

			//{ab,c}=ac,b + bc,a -ab,c
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					//gamma111=g11,1
					gamma111 += 2*_ref->d2[0][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t* 2+0];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0]
						+ 2*_ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0]
						- 2*_ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0];
					//gamma121=g11,2
					gamma121 += 2*_ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0];
					//gamma122=g22,1
					gamma122 += 2*_ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2*_ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0]
						+ 2*_ref->d2[3][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0]
						- 2*_ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0];
					//gamma222=g22,2
					gamma222 += 2*_ref->d2[3][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0];

					gamma111 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->buf_xi[t]
						+ 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1]
						- 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->node[s *2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1]
						+ 2 * _ref->d2[3][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1]
						- 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];

					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t]
						+ 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t]
						- 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t]
						+ 2 * _ref->d2[3][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t]
						- 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];

				}
			}
			val -= (gamma111 * v1 * v1 + 2 * gamma121 * v1 * v2 + gamma221 * v2 * v2) * s1 / sdet;
			val -= (gamma112 * v1 * v1 + 2 * gamma122 * v1 * v2 + gamma222 * v2 * v2) * s2 / sdet;


			return val;
		}
		double curvature_Z(double *ptr,double v1, double v2, bool accurate)
		{
			double val = 0;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}


			//val += (_ref->get__Gammaijk(0, 0, 0) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 0) * v1 * v2 + _ref->get__Gammaijk(1, 1, 0) * v2 * v2) * S1;
			//val += (_ref->get__Gammaijk(0, 0, 1) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 1) * v1 * v2 + _ref->get__Gammaijk(1, 1, 1) * v2 * v2) * S2;
			double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;

			double _gamma111 = _ref->get__Gammaijk(0, 0, 0) * this->get_gij2(0, 0) + _ref->get__Gammaijk(0, 0, 1) * this->get_gij2(1, 0);
			double _gamma112 = _ref->get__Gammaijk(0, 0, 0) * this->get_gij2(0, 1) + _ref->get__Gammaijk(0, 0, 1) * this->get_gij2(1, 1);
			double _gamma121 = _ref->get__Gammaijk(0, 1, 0) * this->get_gij2(0, 0) + _ref->get__Gammaijk(0, 1, 1) * this->get_gij2(1, 0);
			double _gamma122 = _ref->get__Gammaijk(0, 1, 0) * this->get_gij2(0, 1) + _ref->get__Gammaijk(0, 1, 1) * this->get_gij2(1, 1);
			double _gamma221 = _ref->get__Gammaijk(1, 1, 0) * this->get_gij2(0, 0) + _ref->get__Gammaijk(1, 1, 1) * this->get_gij2(1, 0);
			double _gamma222 = _ref->get__Gammaijk(1, 1, 0) * this->get_gij2(0, 1) + _ref->get__Gammaijk(1, 1, 1) * this->get_gij2(1, 1);
			
			//{ab,c}=ac,b + bc,a -ab,c
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0]
						+ 2 * _ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0]
						- 2 * _ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0]
						+ 2 * _ref->d2[3][s] * _ref->node[s * 2 + 0] * _ref->d1[0][t] * _ref->node[t * 2 + 0]
						- 2 * _ref->d2[1][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->node[s * 2 + 0] * _ref->d1[1][t] * _ref->node[t * 2 + 0];

					gamma111 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->buf_xi[t]
						+ 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1]
						- 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1]
						+ 2 * _ref->d2[3][s] * _ref->node[s * 2 + 1] * _ref->d1[0][t] * _ref->node[t * 2 + 1]
						- 2 * _ref->d2[1][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->node[s * 2 + 1] * _ref->d1[1][t] * _ref->node[t * 2 + 1];

					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t]
						+ 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t]
						- 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t]
						+ 2 * _ref->d2[3][s] * _ref->buf_z[s] * _ref->d1[0][t] * _ref->buf_z[t]
						- 2 * _ref->d2[1][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];
				}
			}
			double* ptr1 = ptr;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
			double sdet = sqrt(det);
			double _s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			double _s2 = this->get_gij2(0,0) * v1 + this->get_gij2(0, 1) * v2;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double g11 = 0, g12 = 0, g21 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				g21 = g12;

				s1 = -g21 * v1 - g22 * v2;
				s2 = g11 * v1 + g12 * v2;
				double dsdet = 0.5 * (g11 * this->get_Gij2(0, 0) + 2 * g12 * this->get_Gij2(0, 1) + g22 * this->get_Gij2(1, 1)) * sdet;

				double __gamma111 = 0, __gamma112 = 0, __gamma121 = 0, __gamma122 = 0, __gamma221 = 0, __gamma222 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					//for (int t = 0; t < _ref->_nNode; t++)
					{

						__gamma111 += 2 * _ref->d2[0][t] * _ref->buf_z[t] * _ref->d1[0][s] ;
						//gamma112=2*g12,1-g11,2
						__gamma112 += 2 * _ref->d2[0][t] * _ref->buf_z[t] * _ref->d1[1][s]
							+ 2 * _ref->d2[1][t] * _ref->buf_z[t] * _ref->d1[1][s] 
							- 2 * _ref->d2[1][t] * _ref->buf_z[t] * _ref->d1[0][s];
						//gamma121=g11,2
						__gamma121 += 2 * _ref->d2[1][t] * _ref->buf_z[t] * _ref->d1[0][s];
						//gamma122=g22,1
						__gamma122 += 2 * _ref->d2[1][t] * _ref->buf_z[t] * _ref->d1[1][s];
						//gamma221=2*g12,2-g22,1
						__gamma221 += 2 * _ref->d2[1][t] * _ref->buf_z[t] * _ref->d1[1][s]
							+ 2 * _ref->d2[3][t] * _ref->buf_z[t] * _ref->d1[0][s] 
							- 2 * _ref->d2[1][t] * _ref->buf_z[t] * _ref->d1[1][s];
						//gamma222=g22,2
						__gamma222 += 2 * _ref->d2[3][t] * _ref->buf_z[t] * _ref->d1[1][s];



						__gamma111 += 2 * _ref->d2[0][s] * _ref->d1[0][t] * _ref->buf_z[t];
						//gamma112=2*g12,1-g11,2
						__gamma112 += 2 * _ref->d2[0][s]  * _ref->d1[1][t] * _ref->buf_z[t]
							+ 2 * _ref->d2[1][s]  * _ref->d1[1][t] * _ref->buf_z[t]
							- 2 * _ref->d2[1][s]  * _ref->d1[0][t] * _ref->buf_z[t];
						//gamma121=g11,2
						__gamma121 += 2 * _ref->d2[1][s]  * _ref->d1[0][t] * _ref->buf_z[t];
						//gamma122=g22,1
						__gamma122 += 2 * _ref->d2[1][s]  * _ref->d1[1][t] * _ref->buf_z[t];
						//gamma221=2*g12,2-g22,1
						__gamma221 += 2 * _ref->d2[1][s] * _ref->d1[1][t] * _ref->buf_z[t]
							+ 2 * _ref->d2[3][s] * _ref->d1[0][t] * _ref->buf_z[t]
							- 2 * _ref->d2[1][s] * _ref->d1[1][t] * _ref->buf_z[t];
						//gamma222=g22,2
						__gamma222 += 2 * _ref->d2[3][s] * _ref->buf_z[s] * _ref->d1[1][t] * _ref->buf_z[t];
					}
				}

				val = (_gamma111 * v1 * v1 + 2 * _gamma121 * v1 * v2 + _gamma221 * v2 * v2) * (s1 / sdet-_s1/sdet/sdet*dsdet);
				val += (_gamma112 * v1 * v1 + 2 * _gamma122 * v1 * v2 + _gamma222 * v2 * v2) * (s2 / sdet - _s2 / sdet / sdet * dsdet);
				val -= (gamma111 * v1 * v1 + 2 * gamma121 * v1 * v2 + gamma221 * v2 * v2) * (s1 / sdet - _s1 / sdet / sdet * dsdet);
				val -= (gamma112 * v1 * v1 + 2 * gamma122 * v1 * v2 + gamma222 * v2 * v2) * (s2 / sdet - _s2 / sdet / sdet * dsdet);

				val -= (__gamma111 * v1 * v1 + 2 * __gamma121 * v1 * v2 + __gamma221 * v2 * v2) * (_s1 / sdet);
				val -= (__gamma112 * v1 * v1 + 2 * __gamma122 * v1 * v2 + __gamma222 * v2 * v2) * (_s2 / sdet);

				*ptr1 = val;
				ptr1++;
			}
			return val;
		}
		void curvature_xi(double *ptr,double v1, double v2, bool accurate)
		{

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			
			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);

			double* ptr1 = ptr;
			double sdet = sqrt(det);

			//{ab,c}=ac,b + bc,a -ab,c
			for (int t = 0; t < _ref->_nNode; t++)
				{
				double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
							//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[1][t] 
						+ 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] ;
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] 
						+ 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[0][t] 
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[1][t] ;


					gamma111 += 2 * _ref->d2[0][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][t]  * _ref->d1[1][s] * _ref->buf_xi[s]
						+ 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						- 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][t]  * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][t]  * _ref->d1[1][s] * _ref->buf_xi[s];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][t]  * _ref->d1[1][s] * _ref->buf_xi[s]
						+ 2 * _ref->d2[3][t] * _ref->d1[0][s] * _ref->buf_xi[s]
						- 2 * _ref->d2[1][t]* _ref->d1[1][s] * _ref->buf_xi[s];
					//gamma222=g22,2
				}
				double val = -(gamma111 * v1 * v1 + 2 * gamma121 * v1 * v2 + gamma221 * v2 * v2) * s1/sdet;
				val -= (gamma112 * v1 * v1 + 2 * gamma122 * v1 * v2 + gamma222 * v2 * v2) * s2/sdet;

				*ptr1 = val;
				ptr1++;
			}
	
		}
		void curvature_eta(double* ptr, double v1, double v2, bool accurate)
		{
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double* ptr1 = ptr;
			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);
			double sdet = sqrt(det);
			//{ab,c}=ac,b + bc,a -ab,c
			for (int t = 0; t < _ref->_nNode; t++)
			{
				double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[0][t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[1][t];


					gamma111 += 2 * _ref->d2[0][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						+ 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						- 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						+ 2 * _ref->d2[3][t] * _ref->d1[0][s] * _ref->buf_eta[s]
						- 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s];
					//gamma222=g22,2
				}
				double val = -(gamma111 * v1 * v1 + 2 * gamma121 * v1 * v2 + gamma221 * v2 * v2) * s1/sdet;
				val -= (gamma112 * v1 * v1 + 2 * gamma122 * v1 * v2 + gamma222 * v2 * v2) * s2/sdet;
				*ptr1 = val;
				ptr1++;
			}
		}
		double curvature1(double v1, double v2, bool accurate)
		{
			double val = 0;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			val += (_ref->get__Gammaijk(0, 0, 0) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 0) * v1 * v2 + _ref->get__Gammaijk(1, 1, 0) * v2 * v2) * S1;
			val += (_ref->get__Gammaijk(0, 0, 1) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 1) * v1 * v2 + _ref->get__Gammaijk(1, 1, 1) * v2 * v2) * S2;
			double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;

			//{ab,c}=ac,b + bc,a -ab,c
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t]
						+ 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t]
						+ 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t];

					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t]
						+ 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t]
						+ 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t];

				}
			}
			//val -= (gamma111 * v1 * v1 + 2 * gamma121 * v1 * v2 + gamma221 * v2 * v2) * s1;
			//val -= (gamma112 * v1 * v1 + 2 * gamma122 * v1 * v2 + gamma222 * v2 * v2) * s2;
			val = gamma111 * v1 * v1 - _ref->get__Gammaijk(0, 0, 0) * v1 * v1;
			val +=2*( gamma121 * v1 * v2 - _ref->get__Gammaijk(0, 1, 0) * v1 * v2);
			val += gamma221 * v2 * v2 - _ref->get__Gammaijk(1, 1, 0) * v2 * v2;


			return val;
		}
		void curvature1_xi(double* ptr, double v1, double v2, bool accurate)
		{

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			double* ptr1 = ptr;


			//{ab,c}=ac,b + bc,a -ab,c
			for (int t = 0; t < _ref->_nNode; t++)
			{
				double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[0][t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[1][t];


					gamma111 += 2 * _ref->d2[0][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						+ 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						- 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						+ 2 * _ref->d2[3][t] * _ref->d1[0][s] * _ref->buf_xi[s]
						- 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s];
					//gamma222=g22,2
				}
				double val = gamma111 * v1 * v1 ;
				val += 2 * (gamma121 * v1 * v2);
				val += gamma221 * v2 * v2 ;


				*ptr1 = val;
				ptr1++;
			}

		}
		void curvature1_eta(double* ptr, double v1, double v2, bool accurate)
		{
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			double* ptr1 = ptr;


			//{ab,c}=ac,b + bc,a -ab,c
			for (int t = 0; t < _ref->_nNode; t++)
			{
				double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[0][t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[1][t];


					gamma111 += 2 * _ref->d2[0][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						+ 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						- 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						+ 2 * _ref->d2[3][t] * _ref->d1[0][s] * _ref->buf_eta[s]
						- 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s];
					//gamma222=g22,2
				}
				double val = gamma111 * v1 * v1;
				val += 2 * (gamma121 * v1 * v2);
				val += gamma221 * v2 * v2;

				*ptr1 = val;
				ptr1++;
			}
		}

		double curvature2(double v1, double v2, bool accurate)
		{
			double val = 0;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			val += (_ref->get__Gammaijk(0, 0, 0) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 0) * v1 * v2 + _ref->get__Gammaijk(1, 1, 0) * v2 * v2) * S1;
			val += (_ref->get__Gammaijk(0, 0, 1) * v1 * v1 + 2 * _ref->get__Gammaijk(0, 1, 1) * v1 * v2 + _ref->get__Gammaijk(1, 1, 1) * v2 * v2) * S2;
			double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;

			//{ab,c}=ac,b + bc,a -ab,c
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t]
						+ 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t]
						+ 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[0][t] * _ref->buf_xi[t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[1][t] * _ref->buf_xi[t];

					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t]
						+ 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t]
						+ 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[0][t] * _ref->buf_eta[t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[1][t] * _ref->buf_eta[t];

				}
			}
			//val -= (gamma111 * v1 * v1 + 2 * gamma121 * v1 * v2 + gamma221 * v2 * v2) * s1;
			//val -= (gamma112 * v1 * v1 + 2 * gamma122 * v1 * v2 + gamma222 * v2 * v2) * s2;
			val = gamma112 * v1 * v1 - _ref->get__Gammaijk(0, 0, 1) * v1 * v1;
			val += 2 * (gamma122 * v1 * v2 - _ref->get__Gammaijk(0, 1, 1) * v1 * v2);
			val += gamma222 * v2 * v2 - _ref->get__Gammaijk(1, 1, 1) * v2 * v2;


			return val;
		}
		void curvature2_xi(double* ptr, double v1, double v2, bool accurate)
		{

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			double* ptr1 = ptr;


			//{ab,c}=ac,b + bc,a -ab,c
			for (int t = 0; t < _ref->_nNode; t++)
			{
				double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[0][t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[0][t]
						- 2 * _ref->d2[1][s] * _ref->buf_xi[s] * _ref->d1[1][t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_xi[s] * _ref->d1[1][t];


					gamma111 += 2 * _ref->d2[0][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						+ 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						- 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_xi[s];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s]
						+ 2 * _ref->d2[3][t] * _ref->d1[0][s] * _ref->buf_xi[s]
						- 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_xi[s];
					//gamma222=g22,2
				}
				double val = gamma112 * v1 * v1;
				val += 2 * (gamma122 * v1 * v2);
				val += gamma222 * v2 * v2;


				*ptr1 = val;
				ptr1++;
			}

		}
		void curvature2_eta(double* ptr, double v1, double v2, bool accurate)
		{
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			double* ptr1 = ptr;


			//{ab,c}=ac,b + bc,a -ab,c
			for (int t = 0; t < _ref->_nNode; t++)
			{
				double gamma111 = 0, gamma112 = 0, gamma121 = 0, gamma122 = 0, gamma221 = 0, gamma222 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					//gamma111=g11,1
					gamma111 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[0][t];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t]
						+ 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[0][t]
						- 2 * _ref->d2[1][s] * _ref->buf_eta[s] * _ref->d1[1][t];
					//gamma222=g22,2
					gamma222 += 2 * _ref->d2[3][s] * _ref->buf_eta[s] * _ref->d1[1][t];


					gamma111 += 2 * _ref->d2[0][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma112=2*g12,1-g11,2
					gamma112 += 2 * _ref->d2[0][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						+ 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						- 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma121=g11,2
					gamma121 += 2 * _ref->d2[1][t] * _ref->d1[0][s] * _ref->buf_eta[s];
					//gamma122=g22,1
					gamma122 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s];
					//gamma221=2*g12,2-g22,1
					gamma221 += 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s]
						+ 2 * _ref->d2[3][t] * _ref->d1[0][s] * _ref->buf_eta[s]
						- 2 * _ref->d2[1][t] * _ref->d1[1][s] * _ref->buf_eta[s];
					//gamma222=g22,2
				}
				double val = gamma112 * v1 * v1;
				val += 2 * (gamma122 * v1 * v2);
				val += gamma222 * v2 * v2;

				*ptr1 = val;
				ptr1++;
			}
		}
		double guideBC(double v1, double v2, bool accurate, int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double e21 = e12;
			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

			double tr2 = 1;// _ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			//s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			//s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;

			//double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0,1);
			//double sdet = sqrt(det);

			//scale = 1;
			val = scale * (e11 * v1 * v1 * this->get_gij(0, 1) + e11 * v1 * v2 * this->get_gij(1, 1) + e12 * v2 * v1 * this->get_gij(0, 1) + e12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
			val += -scale * (e21 * v1 * v1 * this->get_gij(0, 0) + e21 * v1 * v2 * this->get_gij(1, 0) + e22 * v2 * v1 * this->get_gij(0, 0) + e22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;
			return val;
		}

		void guideBC_xi(double* ptr, double v1, double v2, bool accurate, int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr2 = 1;// _ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double _e11 = e11, _e12 = e12, _e22 = e22, _e21 = e12;
			double* ptr1 = ptr;
			s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;
			//double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);
			//double sdet = sqrt(det);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				double e21 = e12;
				//scale = 1;
				val = scale * (e11 * v1 * v1 * this->get_gij(0, 1) + e11 * v1 * v2 * this->get_gij(1, 1) + e12 * v2 * v1 * this->get_gij(0, 1) + e12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
				val += -scale * (e21 * v1 * v1 * this->get_gij(0, 0) + e21 * v1 * v2 * this->get_gij(1, 0) + e22 * v2 * v1 * this->get_gij(0, 0) + e22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;

				//double dtr=( e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
				double dtr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

				val += -1. / tr / tr * dtr * (_e11 * v1 * v1 * this->get_gij(0, 1) + _e11 * v1 * v2 * this->get_gij(1, 1) + _e12 * v2 * v1 * this->get_gij(0, 1) + _e12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
				val -= -1. / tr / tr * dtr * (_e21 * v1 * v1 * this->get_gij(0, 0) + _e21 * v1 * v2 * this->get_gij(1, 0) + _e22 * v2 * v1 * this->get_gij(0, 0) + _e22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}

		}
		void guideBC_eta(double* ptr, double v1, double v2, bool accurate, int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1)
			{
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

			double tr2 = 1;
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double _e11 = e11, _e12 = e12, _e22 = e22, _e21 = e12;
			double* ptr1 = ptr;
			s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;
			//double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);
			//double sdet = sqrt(det);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}

				double e21 = e12;
				//scale = 1;
				val = scale * (e11 * v1 * v1 * this->get_gij(0, 1) + e11 * v1 * v2 * this->get_gij(1, 1) + e12 * v2 * v1 * this->get_gij(0, 1) + e12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
				val += -scale * (e21 * v1 * v1 * this->get_gij(0, 0) + e21 * v1 * v2 * this->get_gij(1, 0) + e22 * v2 * v1 * this->get_gij(0, 0) + e22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;

				//double dtr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
				double dtr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

				val += -1. / tr / tr * dtr * (_e11 * v1 * v1 * this->get_gij(0, 1) + _e11 * v1 * v2 * this->get_gij(1, 1) + _e12 * v2 * v1 * this->get_gij(0, 1) + _e12 * v2 * v2 * this->get_gij(1, 1)) / _ref->_refDv;
				val -= -1. / tr / tr * dtr * (_e21 * v1 * v1 * this->get_gij(0, 0) + _e21 * v1 * v2 * this->get_gij(1, 0) + _e22 * v2 * v1 * this->get_gij(0, 0) + _e22 * v2 * v2 * this->get_gij(1, 0)) / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}

		}
		double guide(double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;
			
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double e21 = e12;
			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

			double tr2 = 1;// _ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;

			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0,1);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;
			//scale = 1;
			//val = scale * (e11 * v1 * v1 * this->get_gij2(0, 1) + e11 * v1 * v2 * this->get_gij2(1, 1) + e12 * v2 * v1 * this->get_gij2(0, 1) + e12 * v2 * v2 * this->get_gij2(1, 1)) / _ref->_refDv;
			//val += -scale * (e21 * v1 * v1 * this->get_gij2(0, 0) + e21 * v1 * v2 * this->get_gij2(1, 0) + e22 * v2 * v1 * this->get_gij2(0, 0) + e22 * v2 * v2 * this->get_gij2(1, 0)) / _ref->_refDv;
			
			val = e11 * s1 * v1 + e12 * s1 * v2 + e21 * s2 * v1 + e22 * s2 * v2;
			return val;
		}
		
		void guide_xi(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr2 = 1;// _ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double _e11 = e11, _e12 = e12, _e22 = e22, _e21 = e12;
			double* ptr1 = ptr;
			s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				double e21 = e12;
				//scale = 1;
				/*val = scale * (e11 * v1 * v1 * this->get_gij2(0, 1) + e11 * v1 * v2 * this->get_gij2(1, 1) + e12 * v2 * v1 * this->get_gij2(0, 1) + e12 * v2 * v2 * this->get_gij2(1, 1)) / _ref->_refDv;
				val += -scale * (e21 * v1 * v1 * this->get_gij2(0, 0) + e21 * v1 * v2 * this->get_gij2(1, 0) + e22 * v2 * v1 * this->get_gij2(0, 0) + e22 * v2 * v2 * this->get_gij2(1, 0)) / _ref->_refDv;

				//double dtr=( e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
				double dtr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

				val +=-1./tr/tr*dtr* (_e11 * v1 * v1 * this->get_gij2(0, 1) + _e11 * v1 * v2 * this->get_gij2(1, 1) + _e12 * v2 * v1 * this->get_gij2(0, 1) + _e12 * v2 * v2 * this->get_gij2(1, 1)) / _ref->_refDv;
				val -= -1. / tr / tr * dtr * (_e21 * v1 * v1 * this->get_gij2(0, 0) + _e21 * v1 * v2 * this->get_gij2(1, 0) + _e22 * v2 * v1 * this->get_gij2(0, 0) + _e22 * v2 * v2 * this->get_gij2(1, 0)) / _ref->_refDv;
				*/
				val = e11 * s1 * v1 +e12 * s1 * v2 + e21 * s2 * v1 + e22 * s2 * v2;

				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_eta(double* ptr, double v1, double v2, bool accurate, int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1)
			{
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));
			
			double tr2 = 1;
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double _e11 = e11, _e12 = e12, _e22 = e22, _e21 = e12;
			double* ptr1 =ptr;
			s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				
				double e21 = e12;
				//scale = 1;
				//val = scale * (e11 * v1 * v1 * this->get_gij2(0, 1) + e11 * v1 * v2 * this->get_gij2(1, 1) + e12 * v2 * v1 * this->get_gij2(0, 1) + e12 * v2 * v2 * this->get_gij2(1, 1)) / _ref->_refDv;
				//val += -scale * (e21 * v1 * v1 * this->get_gij2(0, 0) + e21 * v1 * v2 * this->get_gij2(1, 0) + e22 * v2 * v1 * this->get_gij2(0, 0) + e22 * v2 * v2 * this->get_gij2(1, 0)) / _ref->_refDv;

				//double dtr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
				//double dtr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

				//val += -1. / tr / tr * dtr * (_e11 * v1 * v1 * this->get_gij2(0, 1) + _e11 * v1 * v2 * this->get_gij2(1, 1) + _e12 * v2 * v1 * this->get_gij2(0, 1) + _e12 * v2 * v2 * this->get_gij2(1, 1)) / _ref->_refDv;
				//val -= -1. / tr / tr * dtr * (_e21 * v1 * v1 * this->get_gij2(0, 0) + _e21 * v1 * v2 * this->get_gij2(1, 0) + _e22 * v2 * v1 * this->get_gij2(0, 0) + _e22 * v2 * v2 * this->get_gij2(1, 0)) / _ref->_refDv;
				
				val = e11 * s1 * v1 + e12 * s1 * v2 + e21 * s2 * v1 + e22 * s2 * v2;

				*ptr1 = val;
				ptr1++;
			}

		}
		
		void guide_Z(double* ptr, double v1, double v2, bool accurate)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double e21 = e12;
			//double tr = (e11 + e22) / (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));
			double tr2 = 1;
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double* ptr1 = ptr;
			//double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(0, 1);
			//double sdet = sqrt(det);
			double _s1 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
			double _s2 = -this->get_gij2(0, 0) * v1 - this->get_gij2(0, 1) * v2;
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
				double g11 = 0, g12 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * _ref->d1[0][t] * _ref->buf_z[t] * _ref->d1[0][s];
					g12 += _ref->d1[0][t] * _ref->buf_z[t] * _ref->d1[1][s] + _ref->d1[1][t] * _ref->buf_z[t] * _ref->d1[0][s];
					g22 += 2 * _ref->d1[1][t] * _ref->buf_z[t] * _ref->d1[1][s];
				}
				double g21 = g12;

				//scale = 1;
				
				val = scale * (e11 * v1 * v1 * g12 + e11 * v1 * v2 * g22 + e12 * v2 * v1 * g12 + e12 * v2 * v2 * g22) / _ref->_refDv;
				val += -scale * (e21 * v1 * v1 * g11 + e21 * v1 * v2 * g21 + e22 * v2 * v1 *g11 + e22 * v2 * v2 * g21) / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}

		}
		double guide4(double v1, double v2, bool accurate)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			val = (this->get_gij2(0, 0) * v1 * v1 * _ref->get__gij(0, 1) + this->get_gij2(0, 0) * v1 * v2 * _ref->get__gij(1, 1) + this->get_gij2(0, 1) * v2 * v1 * _ref->get__gij(0, 1) + this->get_gij2(0, 1) * v2 * v2 * _ref->get__gij(1, 1)) / _ref->_refDv;
			val +=  (this->get_gij2(1, 0) * v1 * v1 * _ref->get__gij(0, 0) + this->get_gij2(1, 0) * v1 * v2 * _ref->get__gij(1, 0) + this->get_gij2(1, 1) * v2 * v1 * _ref->get__gij(0, 0) + this->get_gij2(1, 1)  * v2 * v2 * _ref->get__gij(1, 0)) / _ref->_refDv;
			return val;
		}
		void guide4_Z(double *ptr,double v1, double v2, bool accurate)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double g11 = 0, g12 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * _ref->d1[0][t] * _ref->buf_z[t] * _ref->d1[0][s];
					g12 += _ref->d1[0][t] * _ref->buf_z[t] * _ref->d1[1][s] + _ref->d1[1][t] * _ref->buf_z[t] * _ref->d1[0][s];
					g22 += 2 * _ref->d1[1][t] * _ref->buf_z[t] * _ref->d1[1][s];
				}
				double g21 = g12;
				val = (g11 * v1 * v1 * _ref->get__gij(0, 1) + g11 * v1 * v2 * _ref->get__gij(1, 1) + g12 * v2 * v1 * _ref->get__gij(0, 1) + g12 * v2 * v2 * _ref->get__gij(1, 1)) / _ref->_refDv;
				val += (g21 * v1 * v1 * _ref->get__gij(0, 0) + g21 * v1 * v2 * _ref->get__gij(1, 0) + g22 * v2 * v1 * _ref->get__gij(0, 0) + g22 * v2 * v2 * _ref->get__gij(1, 0)) / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}
		}

		double guide5(double v1, double v2, bool accurate)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12; 
			double det = (S11 * S22 - S12 * S12);
			double tr = S11 + S22;
			return (tr*tr-2*det)*sc;
		}
		void guide5_Z(double* ptr, bool accurate)
		{
			double val = 0;

	
			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double det = (S11 * S22 - S12 * S12);
			double tr = S11 + S22;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = 0, _S12 = 0, _S22 = 0;
				
				_S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) ;
				_S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) ;
				_S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) ;
				
				double _S21 = _S12;
				
				double S21 = S12;
				double ddet = (_S11 * S22 - _S12 * S12)+ (S11 * _S22 - S12 * _S12);
				double dtr = _S11 + _S22;
				*ptr1=(2*dtr * tr - 2 * ddet) * sc;
				ptr1++;
			}
		}
		double guide6( bool accurate)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double det = (S11 * S22 - S12 * S12);
			double tr = (S11 + S22);
			//return tr/_ref->_refDv;
			return (tr * tr - 2 * det) * sc;
		}
		void guide6_Z(double* ptr, bool accurate)
		{
			double val = 0;


			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double det = (S11 * S22 - S12 * S12);
			double tr = S11 + S22;
			val = (tr * tr - 2 * det);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = 0, _S12 = 0, _S22 = 0;

				_S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				_S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				_S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				double _S21 = _S12;

				double S21 = S12;
				double ddet = (_S11 * S22 - _S12 * S12) + (S11 * _S22 - S12 * _S12);
				double dtr = _S11 + _S22;
					
				*ptr1 = 2*tr*dtr*sc-2*ddet*sc;
				//*ptr1 = dtr / _ref->_refDv;
				ptr1++;
			}
		}
		double guide7(bool accurate)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double det = (S11 * S22 - S12 * S12);
			double tr = (S11 + S22);
			return tr/_ref->_refDv;
			//return (tr * tr - 2 * det) * sc;
		}
		void guide7_Z(double* ptr, bool accurate)
		{
			double val = 0;


			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double det = (S11 * S22 - S12 * S12);
			double tr = S11 + S22;
			val = (tr * tr - 2 * det);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = 0, _S12 = 0, _S22 = 0;

				_S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				_S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				_S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				double _S21 = _S12;

				double S21 = S12;
				double ddet = (_S11 * S22 - _S12 * S12) + (S11 * _S22 - S12 * _S12);
				double dtr = _S11 + _S22;

				//*ptr1 = 2 * tr * dtr * sc - 2 * ddet * sc;
				*ptr1 = dtr / _ref->_refDv;
				ptr1++;
			}
		}
		double vec(double V1,double V2)
		{
			

			double __mu = 0, __nu = 0;

			double v1 = V1 * _ref->get__Gi(0, 0)+ V2 * _ref->get__Gi(0, 1);
			double v2 = V1 * _ref->get__Gi(1, 0)+ V2 * _ref->get__Gi(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{
				__mu += _ref->d0[s] * _ref->buf_mu[s];
				__nu += _ref->d0[s] * _ref->buf_nu[s];
			}


			return __mu*v1 + __nu*v2;

		}
		void vec_mu(double* ptr, double V1, double V2)
		{

			double __mu = 0, __nu = 0;
			double* ptr1 = ptr;
			double v1 = V1 * _ref->get__Gi(0, 0) + V2 * _ref->get__Gi(0, 1);
			double v2 = V1 * _ref->get__Gi(1, 0) + V2 * _ref->get__Gi(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{
				__mu = _ref->d0[s];


				double val = __mu * v1;
				*ptr1 = val;
				ptr1++;
			}
		}
		void vec_nu(double* ptr, double V1, double V2)
		{
			double __mu = 0, __nu = 0;
			double* ptr1 = ptr;
			double v1 = V1 * _ref->get__Gi(0, 0) + V2 * _ref->get__Gi(0, 1);
			double v2 = V1 * _ref->get__Gi(1, 0) + V2 * _ref->get__Gi(1, 1);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				__nu = _ref->d0[s];


				double val = __nu * v2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double unit()
		{
			double __mu = 0, __nu = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				__mu += _ref->d0[i] * _ref->buf_mu[i];
				__nu += _ref->d0[i] * _ref->buf_nu[i];
				
			}

			return __mu * __mu * _ref->get__Gij(0, 0) + 2 * __mu * __nu * _ref->get__Gij(0, 1) + __nu * __nu * _ref->get__Gij(1, 1)-1;

		}
		void unit_mu(double *ptr)
		{
			double __mu = 0, __nu = 0;
			double *ptr1 = ptr;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				__mu += _ref->d0[i] * _ref->buf_mu[i];
				__nu += _ref->d0[i] * _ref->buf_nu[i];

			}

			for (int i = 0; i < _ref->_nNode; i++)
			{
				double __dmu = _ref->d0[i];
				double __dnu = 0;// _ref->d0[i] * _ref->buf_nu[i];
				double val=2*__dmu * __mu * _ref->get__Gij(0, 0) + 2 * (__mu * __dnu+__dmu * __nu) * _ref->get__Gij(0, 1) + 2*(__dnu * __nu) * _ref->get__Gij(1, 1) ;
				*ptr1 = val;
				ptr1++;
			}


		}
		void unit_nu(double* ptr)
		{
			double __mu = 0, __nu = 0;
			double* ptr1 = ptr;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				__mu += _ref->d0[i] * _ref->buf_mu[i];
				__nu += _ref->d0[i] * _ref->buf_nu[i];

			}

			for (int i = 0; i < _ref->_nNode; i++)
			{
				double __dmu = 0;// _ref->d0[i];
				double __dnu =  _ref->d0[i];
				double val = 2 * __dmu * __mu * _ref->get__Gij(0, 0) + 2 * (__mu * __dnu + __dmu * __nu) * _ref->get__Gij(0, 1) + 2 * (__dnu * __nu) * _ref->get__Gij(1, 1);
				*ptr1 = val;
				ptr1++;
			}


		}
		double get_v1()
		{
			double f_u = 0,f_v=0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f_u += (_ref->d1[0][s]) * _ref->buf_mu[s];
				f_v += (_ref->d1[1][s]) * _ref->buf_mu[s];
			}
			return f_v;
		}
		double get_v2()
		{
			double f_u = 0, f_v = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f_u += (_ref->d1[0][s]) * _ref->buf_mu[s];
				f_v += (_ref->d1[1][s]) * _ref->buf_mu[s];
			}
			return -f_u;

		}
		double curl_free()
		{
			double xi_y = 0, eta_x = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xi_y += (_ref->d1[0][s]) * _ref->buf_xi[s] * _ref->get__Gi(0, 1) + (_ref->d1[1][s]) * _ref->buf_xi[s] * _ref->get__Gi(1, 1);
				eta_x += (_ref->d1[0][s]) * _ref->buf_eta[s] * _ref->get__Gi(0, 0) + (_ref->d1[1][s]) * _ref->buf_eta[s] * _ref->get__Gi(1, 0);

			}
			return xi_y - eta_x;
		}
		void  curl_free_xi(double* ptr)
		{
			double xi_y = 0, eta_x = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xi_y += (_ref->d1[0][s]) * _ref->buf_xi[s] * _ref->get__Gi(0, 1) + (_ref->d1[1][s]) * _ref->buf_xi[s] * _ref->get__Gi(1, 1);
				eta_x += (_ref->d1[0][s]) * _ref->buf_eta[s] * _ref->get__Gi(0, 0) + (_ref->d1[1][s]) * _ref->buf_eta[s] * _ref->get__Gi(1, 0);

			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi_y = (_ref->d1[0][s]) * _ref->get__Gi(0, 1) + (_ref->d1[1][s]) * _ref->get__Gi(1, 1);
				double _eta_x = (_ref->d1[0][s]) * _ref->get__Gi(0, 0) + (_ref->d1[1][s]) * _ref->get__Gi(1, 0);
				*ptr1 = _xi_y;
				ptr1++;
			}
		}
		void  curl_free_eta(double* ptr)
		{
			double xi_y = 0, eta_x = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xi_y += (_ref->d1[0][s]) * _ref->buf_xi[s] * _ref->get__Gi(0, 1) + (_ref->d1[1][s]) * _ref->buf_xi[s] * _ref->get__Gi(1, 1);
				eta_x += (_ref->d1[0][s]) * _ref->buf_eta[s] * _ref->get__Gi(0, 0) + (_ref->d1[1][s]) * _ref->buf_eta[s] * _ref->get__Gi(1, 0);

			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi_y = (_ref->d1[0][s]) * _ref->get__Gi(0, 1) + (_ref->d1[1][s]) * _ref->get__Gi(1, 1);
				double _eta_x = (_ref->d1[0][s]) * _ref->get__Gi(0, 0) + (_ref->d1[1][s]) * _ref->get__Gi(1, 0);
				*ptr1 =  - _eta_x;
				ptr1++;
			}
		}
		double fair()
		{
			double f_uu = 0, f_uv = 0, f_vv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f_uu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_mu[s];
				f_uv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_mu[s];
				f_vv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_mu[s];
			}

			//return mu_u * _ref->get__Gi(0, 0) + mu_v * _ref->get__Gi(1, 0)
			//	+ nu_u * _ref->get__Gi(0, 1) + nu_v * _ref->get__Gi(1, 1);
			return f_uu * _ref->get__Gij(0, 0) + 2 * f_uv * _ref->get__Gij(0, 1) + f_vv * _ref->get__Gij(1, 1);

		}

		void fair_mu(double* ptr)
		{
			double* ptr1 = ptr;

			double f_uu = 0, f_uv = 0, f_vv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f_uu = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				f_uv = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				f_vv = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			

				double val= f_uu * _ref->get__Gij(0, 0) + 2 * f_uv * _ref->get__Gij(0, 1) + f_vv * _ref->get__Gij(1, 1);
				*ptr1 = val;
				ptr1 ++;
			}


		}
		
		double guide8(bool accurate,double v1,double v2,int mode)
		{
			double val = 0;
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double e11 = 0, e12 = 0, e22= 0;
			double g111 = 0, g112 = 0, g121 = 0, g122 = 0, g221 = 0, g222 = 0,g211=0,g212=0;
			
			for (int i = 0; i < _ref->_nNode; i++)
			{
				for (int j = 0; j < _ref->_nNode; j++)
				{
					g111 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					g211 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					g112 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g212 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g122 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g222 += (_ref->d2[3][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);

					g111 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[0][j] * _ref->buf_xi[j]);
					g211 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g112 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g212 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[3][j] * _ref->buf_xi[j]);
					g122 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g222 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[3][j] * _ref->buf_xi[j]);
					
					g111 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					g211 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					g112 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g212 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g122 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g222 += (_ref->d2[3][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);

					g111 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[0][j] * _ref->buf_eta[j]);
					g211 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g112 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g212 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[3][j] * _ref->buf_eta[j]);
					g122 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g222 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[3][j] * _ref->buf_eta[j]);

					e11 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					e12 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					e22 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					e11 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					e12 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					e22 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
				}

			}
			double e21 = e12;
			g111 += -this->get_Gammaijk(0, 0, 0) * e11 - this->get_Gammaijk(0, 0, 1) * e12 - this->get_Gammaijk(0, 0, 0) * e11 - this->get_Gammaijk(0, 0, 1) * e12;
			g211 += -this->get_Gammaijk(1, 0, 0) * e11 - this->get_Gammaijk(1, 0, 1) * e12 - this->get_Gammaijk(1, 0, 0) * e11 - this->get_Gammaijk(1, 0, 1) * e12;
			g112 += -this->get_Gammaijk(0, 0, 0) * e21 - this->get_Gammaijk(0, 0, 1) * e22 - this->get_Gammaijk(0, 1, 0) * e11 - this->get_Gammaijk(0, 1, 1) * e12;
			g212 += -this->get_Gammaijk(1, 0, 0) * e21 - this->get_Gammaijk(1, 0, 1) * e22 - this->get_Gammaijk(1, 1, 0) * e11 - this->get_Gammaijk(1, 1, 1) * e12;
			g122 += -this->get_Gammaijk(0, 1, 0) * e21 - this->get_Gammaijk(0, 1, 1) * e22 - this->get_Gammaijk(0, 1, 0) * e21 - this->get_Gammaijk(0, 1, 1) * e22;
			g222 += -this->get_Gammaijk(1, 1, 0) * e21 - this->get_Gammaijk(1, 1, 1) * e22 - this->get_Gammaijk(1, 1, 0) * e21 - this->get_Gammaijk(1, 1, 1) * e22;
			g121 = g112;
			g221 = g212;
			
			double det = get_gij2(0, 0) * get_gij2(1, 1) - get_gij2(0, 1) * get_gij2(0, 1);
			double s1 = (v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1))/sqrt(det);
			double s2 = ( - v1 * get_gij2(0, 0) - v2 * get_gij2(1, 0)) / sqrt(det);
			double t1 = s1, t2 = s2;
			if (mode == 1) {
				t1 = v1; t2 = v2;
			}
			//double tr2 = ((g111 * t1+ g211 * t2) * _ref->get__Gij(0, 0) + 2 * (g112 * t1+ g212 * t2) * _ref->get__Gij(0, 1) + (g122 * t1+ g222 * t2) * _ref->get__Gij(1, 1));
			//tr = 1;
			double tr = g111 * t1 * _ref->get__Gij(0, 0) + g112 * t1 * _ref->get__Gij(0, 1) + g121 * t1 * _ref->get__Gij(1, 0) + g122 * t1 * _ref->get__Gij(1, 1) +
				g211 * t2 * _ref->get__Gij(0, 0) + g212 * t2 * _ref->get__Gij(0, 1) + g221 * t2 * _ref->get__Gij(1, 0) + g222 * t2 * _ref->get__Gij(1, 1);

			//val = (g111 * t1 * v1 * s1 + g112 * t1 * v1 * s2 + g121 * t1 * v2 * s1 + g122 * t1 * v2 * s2 +
			//	g211 * t2 * v1 * s1 + g212 * t2 * v1 * s2 + g221 * t2 * v2 * s1 + g222 * t2 * v2 * s2);
			val = ((g111 * get_Gij(0, 0) * e12 + g111 * get_Gij(0, 1) * e22 + g112 * get_Gij(1, 0) * e12 + g112 * get_Gij(1, 1) * e22) * t1 +
				(g211 * get_Gij(0, 0) * e12 + g211 * get_Gij(0, 1) * e22 + g212 * get_Gij(1, 0) * e21 + g212 * get_Gij(1, 1) * e22) * t2 -
				(g121 * get_Gij(0, 0) * e11 + g121 * get_Gij(0, 1) * e21 + g122 * get_Gij(1, 0) * e11 + g122 * get_Gij(1, 1) * e21) * t1 -
				(g221 * get_Gij(0, 0) * e11 + g221 * get_Gij(0, 1) * e21 + g222 * get_Gij(1, 0) * e11 + g222 * get_Gij(1, 1) * e21) * t2) / _ref->_refDv;

			return val;
		}
		
		void guide8_xi(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;
			double length = 0;
			if (mode == 1)
			{
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
			}
			double e11 = 0, e12 = 0, e22 = 0;
			double g111 = 0, g112 = 0, g121 = 0, g122 = 0, g221 = 0, g222 = 0, g211 = 0, g212 = 0;

			for (int i = 0; i < _ref->_nNode; i++)
			{
				for (int j = 0; j < _ref->_nNode; j++)
				{
					g111 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					g211 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					g112 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g212 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g122 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g222 += (_ref->d2[3][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);

					g111 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[0][j] * _ref->buf_xi[j]);
					g211 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g112 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g212 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[3][j] * _ref->buf_xi[j]);
					g122 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g222 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[3][j] * _ref->buf_xi[j]);

					g111 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					g211 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					g112 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g212 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g122 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g222 += (_ref->d2[3][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);

					g111 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[0][j] * _ref->buf_eta[j]);
					g211 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g112 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g212 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[3][j] * _ref->buf_eta[j]);
					g122 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g222 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[3][j] * _ref->buf_eta[j]);

					e11 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					e12 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					e22 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					e11 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					e12 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					e22 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
				}

			}
			double e21 = e12;
			g111 += -this->get_Gammaijk(0, 0, 0) * e11 - this->get_Gammaijk(0, 0, 1) * e12 - this->get_Gammaijk(0, 0, 0) * e11 - this->get_Gammaijk(0, 0, 1) * e12;
			g211 += -this->get_Gammaijk(1, 0, 0) * e11 - this->get_Gammaijk(1, 0, 1) * e12 - this->get_Gammaijk(1, 0, 0) * e11 - this->get_Gammaijk(1, 0, 1) * e12;
			g112 += -this->get_Gammaijk(0, 0, 0) * e21 - this->get_Gammaijk(0, 0, 1) * e22 - this->get_Gammaijk(0, 1, 0) * e11 - this->get_Gammaijk(0, 1, 1) * e12;
			g212 += -this->get_Gammaijk(1, 0, 0) * e21 - this->get_Gammaijk(1, 0, 1) * e22 - this->get_Gammaijk(1, 1, 0) * e11 - this->get_Gammaijk(1, 1, 1) * e12;
			g122 += -this->get_Gammaijk(0, 1, 0) * e21 - this->get_Gammaijk(0, 1, 1) * e22 - this->get_Gammaijk(0, 1, 0) * e21 - this->get_Gammaijk(0, 1, 1) * e22;
			g222 += -this->get_Gammaijk(1, 1, 0) * e21 - this->get_Gammaijk(1, 1, 1) * e22 - this->get_Gammaijk(1, 1, 0) * e21 - this->get_Gammaijk(1, 1, 1) * e22;
			g121 = g112;
			g221 = g212;

			double* ptr1 = ptr;
			
			double det = get_gij2(0, 0) * get_gij2(1, 1) - get_gij2(0, 1) * get_gij2(0, 1);
			double s1 = (v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / sqrt(det);
			double s2 = (-v1 * get_gij2(0, 0) - v2 * get_gij2(1, 0)) / sqrt(det);
			double t1 = s1, t2 = s2;
			if (mode == 1) {
				t1 = v1; t2 = v2;
			}

			double tr = g111 * t1 * _ref->get__Gij(0, 0) + g112 * t1 * _ref->get__Gij(0, 1) + g121 * t1 * _ref->get__Gij(1, 0) + g122 * t1 * _ref->get__Gij(1, 1) +
				g211 * t2 * _ref->get__Gij(0, 0) + g212 * t2 * _ref->get__Gij(0, 1) + g221 * t2 * _ref->get__Gij(1, 0) + g222 * t2 * _ref->get__Gij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g111 = 0, _g112 = 0, _g121 = 0, _g122 = 0, _g221 = 0, _g222 = 0, _g211 = 0, _g212 = 0;
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int i = 0; i < _ref->_nNode; i++)
				{

					_g111 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][s]);
					_g211 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[0][s]);
					_g112 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][s]);
					_g212 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][s]);
					_g122 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][s]);
					_g222 += (_ref->d2[3][i] * _ref->buf_xi[i]) * (_ref->d1[1][s]);

					_g111 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[0][s] );
					_g211 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][s] );
					_g112 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][s] );
					_g212 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[3][s] );
					_g122 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[1][s]);
					_g222 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[3][s]);

					_g111 += (_ref->d2[0][s] ) * (_ref->d1[0][i] * _ref->buf_xi[i]);
					_g211 += (_ref->d2[1][s] ) * (_ref->d1[0][i] * _ref->buf_xi[i]);
					_g112 += (_ref->d2[0][s] ) * (_ref->d1[1][i] * _ref->buf_xi[i]);
					_g212 += (_ref->d2[1][s] ) * (_ref->d1[1][i] * _ref->buf_xi[i]);
					_g122 += (_ref->d2[1][s] ) * (_ref->d1[1][i] * _ref->buf_xi[i]);
					_g222 += (_ref->d2[3][s] ) * (_ref->d1[1][i] * _ref->buf_xi[i]);

					_g111 += (_ref->d1[0][s] ) * (_ref->d2[0][i] * _ref->buf_xi[i]);
					_g211 += (_ref->d1[0][s] ) * (_ref->d2[1][i] * _ref->buf_xi[i]);
					_g112 += (_ref->d1[0][s] ) * (_ref->d2[1][i] * _ref->buf_xi[i]);
					_g212 += (_ref->d1[0][s] ) * (_ref->d2[3][i] * _ref->buf_xi[i]);
					_g122 += (_ref->d1[1][s] ) * (_ref->d2[1][i] * _ref->buf_xi[i]);
					_g222 += (_ref->d1[1][s] ) * (_ref->d2[3][i] * _ref->buf_xi[i]);



					_e11 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][s] );
					_e12 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][s]);
					_e22 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][s]);

					_e11 += (_ref->d1[0][s] ) * (_ref->d1[0][i] * _ref->buf_xi[i]);
					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][i] * _ref->buf_xi[i]);
					_e22 += (_ref->d1[1][s]) * (_ref->d1[1][i] * _ref->buf_xi[i]);
				}
				double _e21 = _e12;

				_g111 += -this->get_Gammaijk(0, 0, 0) * _e11 - this->get_Gammaijk(0, 0, 1) * _e12 - this->get_Gammaijk(0, 0, 0) * _e11 - this->get_Gammaijk(0, 0, 1) * _e12;
				_g211 += -this->get_Gammaijk(1, 0, 0) * _e11 - this->get_Gammaijk(1, 0, 1) * _e12 - this->get_Gammaijk(1, 0, 0) * _e11 - this->get_Gammaijk(1, 0, 1) * _e12;
				_g112 += -this->get_Gammaijk(0, 0, 0) * _e21 - this->get_Gammaijk(0, 0, 1) * _e22 - this->get_Gammaijk(0, 1, 0) * _e11 - this->get_Gammaijk(0, 1, 1) * _e12;
				_g212 += -this->get_Gammaijk(1, 0, 0) * _e21 - this->get_Gammaijk(1, 0, 1) * _e22 - this->get_Gammaijk(1, 1, 0) * _e11 - this->get_Gammaijk(1, 1, 1) * _e12;
				_g122 += -this->get_Gammaijk(0, 1, 0) * _e21 - this->get_Gammaijk(0, 1, 1) * _e22 - this->get_Gammaijk(0, 1, 0) * _e21 - this->get_Gammaijk(0, 1, 1) * _e22;
				_g222 += -this->get_Gammaijk(1, 1, 0) * _e21 - this->get_Gammaijk(1, 1, 1) * _e22 - this->get_Gammaijk(1, 1, 0) * _e21 - this->get_Gammaijk(1, 1, 1) * _e22;

				_g121 = _g112;
				_g221 = _g212;
				double dtr = _g111 * t1 * _ref->get__Gij(0, 0) + _g112 * t1 * _ref->get__Gij(0, 1) + _g121 * t1 * _ref->get__Gij(1, 0) + _g122 * t1 * _ref->get__Gij(1, 1) +
					_g211 * t2 * _ref->get__Gij(0, 0) + _g212 * t2 * _ref->get__Gij(0, 1) + _g221 * t2 * _ref->get__Gij(1, 0) + _g222 * t2 * _ref->get__Gij(1, 1);

				//val = (_g111 * t1 * v1 * s1 + _g112 * t1 * v1 * s2 + _g121 * t1 * v2 * s1 + _g122 * t1 * v2 * s2 +
				//	_g211 * t2 * v1 * s1 + _g212 * t2 * v1 * s2 + _g221 * t2 * v2 * s1 + _g222 * t2 * v2 * s2);
				val = ((_g111 * get_Gij(0, 0) * e12 + _g111 * get_Gij(0, 1) * e22 + _g112 * get_Gij(1, 0) * e12 + _g112 * get_Gij(1, 1) * e22) * t1 +
					(_g211 * get_Gij(0, 0) * e12 + _g211 * get_Gij(0, 1) * e22 + _g212 * get_Gij(1, 0) * e21 + _g212 * get_Gij(1, 1) * e22) * t2 -
					(_g121 * get_Gij(0, 0) * e11 + _g121 * get_Gij(0, 1) * e21 + _g122 * get_Gij(1, 0) * e11 + _g122 * get_Gij(1, 1) * e21) * t1 -
					(_g221 * get_Gij(0, 0) * e11 + _g221 * get_Gij(0, 1) * e21 + _g222 * get_Gij(1, 0) * e11 + _g222 * get_Gij(1, 1) * e21) * t2) / _ref->_refDv;

				val+= ((g111 * get_Gij(0, 0) * _e12 + g111 * get_Gij(0, 1) * _e22 + g112 * get_Gij(1, 0) * _e12 + g112 * get_Gij(1, 1) * _e22) * t1 +
					(g211 * get_Gij(0, 0) * _e12 + g211 * get_Gij(0, 1) * _e22 + g212 * get_Gij(1, 0) * _e21 + g212 * get_Gij(1, 1) * _e22) * t2 -
					(g121 * get_Gij(0, 0) * _e11 + g121 * get_Gij(0, 1) * _e21 + g122 * get_Gij(1, 0) * _e11 + g122 * get_Gij(1, 1) * _e21) * t1 -
					(g221 * get_Gij(0, 0) * _e11 + g221 * get_Gij(0, 1) * _e21 + g222 * get_Gij(1, 0) * _e11 + g222 * get_Gij(1, 1) * _e21) * t2) / _ref->_refDv;
					
				*ptr1 = val;
				ptr1++;
			}
			
		}
		void guide8_eta(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;
			double length = 0;
			if (mode == 1)
			{
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
			}
			double e11 = 0, e12 = 0, e22 = 0;
			double g111 = 0, g112 = 0, g121 = 0, g122 = 0, g221 = 0, g222 = 0, g211 = 0, g212 = 0;

			for (int i = 0; i < _ref->_nNode; i++)
			{
				for (int j = 0; j < _ref->_nNode; j++)
				{
					g111 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					g211 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					g112 += (_ref->d2[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g212 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g122 += (_ref->d2[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					g222 += (_ref->d2[3][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);

					g111 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[0][j] * _ref->buf_xi[j]);
					g211 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g112 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g212 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d2[3][j] * _ref->buf_xi[j]);
					g122 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[1][j] * _ref->buf_xi[j]);
					g222 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d2[3][j] * _ref->buf_xi[j]);

					g111 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					g211 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					g112 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g212 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g122 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					g222 += (_ref->d2[3][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);

					g111 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[0][j] * _ref->buf_eta[j]);
					g211 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g112 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g212 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[3][j] * _ref->buf_eta[j]);
					g122 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[1][j] * _ref->buf_eta[j]);
					g222 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[3][j] * _ref->buf_eta[j]);

					e11 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[0][j] * _ref->buf_xi[j]);
					e12 += (_ref->d1[0][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					e22 += (_ref->d1[1][i] * _ref->buf_xi[i]) * (_ref->d1[1][j] * _ref->buf_xi[j]);
					e11 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][j] * _ref->buf_eta[j]);
					e12 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
					e22 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][j] * _ref->buf_eta[j]);
				}

			}
			double e21 = e12;
			g111 += -this->get_Gammaijk(0, 0, 0) * e11 - this->get_Gammaijk(0, 0, 1) * e12 - this->get_Gammaijk(0, 0, 0) * e11 - this->get_Gammaijk(0, 0, 1) * e12;
			g211 += -this->get_Gammaijk(1, 0, 0) * e11 - this->get_Gammaijk(1, 0, 1) * e12 - this->get_Gammaijk(1, 0, 0) * e11 - this->get_Gammaijk(1, 0, 1) * e12;
			g112 += -this->get_Gammaijk(0, 0, 0) * e21 - this->get_Gammaijk(0, 0, 1) * e22 - this->get_Gammaijk(0, 1, 0) * e11 - this->get_Gammaijk(0, 1, 1) * e12;
			g212 += -this->get_Gammaijk(1, 0, 0) * e21 - this->get_Gammaijk(1, 0, 1) * e22 - this->get_Gammaijk(1, 1, 0) * e11 - this->get_Gammaijk(1, 1, 1) * e12;
			g122 += -this->get_Gammaijk(0, 1, 0) * e21 - this->get_Gammaijk(0, 1, 1) * e22 - this->get_Gammaijk(0, 1, 0) * e21 - this->get_Gammaijk(0, 1, 1) * e22;
			g222 += -this->get_Gammaijk(1, 1, 0) * e21 - this->get_Gammaijk(1, 1, 1) * e22 - this->get_Gammaijk(1, 1, 0) * e21 - this->get_Gammaijk(1, 1, 1) * e22;
			g121 = g112;
			g221 = g212;
			double det = get_gij2(0, 0) * get_gij2(1, 1) - get_gij2(0, 1) * get_gij2(0, 1);
			double s1 = (v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / sqrt(det);
			double s2 = (-v1 * get_gij2(0, 0) - v2 * get_gij2(1, 0)) / sqrt(det);
			double t1 = s1, t2 = s2;
			if (mode == 1) {
				t1 = v1; t2 = v2;
			}
			
			double* ptr1 = ptr;
			double tr = g111 * t1 * _ref->get__Gij(0, 0) + g112 * t1 * _ref->get__Gij(0, 1) + g121 * t1 * _ref->get__Gij(1, 0) + g122 * t1 * _ref->get__Gij(1, 1) +
				g211 * t2 * _ref->get__Gij(0, 0) + g212 * t2 * _ref->get__Gij(0, 1) + g221 * t2 * _ref->get__Gij(1, 0) + g222 * t2 * _ref->get__Gij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g111 = 0, _g112 = 0, _g121 = 0, _g122 = 0, _g221 = 0, _g222 = 0, _g211 = 0, _g212 = 0;
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int i = 0; i < _ref->_nNode; i++)
				{

					_g111 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][s]);
					_g211 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[0][s]);
					_g112 += (_ref->d2[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][s]);
					_g212 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][s]);
					_g122 += (_ref->d2[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][s]);
					_g222 += (_ref->d2[3][i] * _ref->buf_eta[i]) * (_ref->d1[1][s]);

					_g111 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[0][s]);
					_g211 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][s]);
					_g112 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[1][s]);
					_g212 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d2[3][s]);
					_g122 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[1][s]);
					_g222 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d2[3][s]);

					_g111 += (_ref->d2[0][s]) * (_ref->d1[0][i] * _ref->buf_eta[i]);
					_g211 += (_ref->d2[1][s]) * (_ref->d1[0][i] * _ref->buf_eta[i]);
					_g112 += (_ref->d2[0][s]) * (_ref->d1[1][i] * _ref->buf_eta[i]);
					_g212 += (_ref->d2[1][s]) * (_ref->d1[1][i] * _ref->buf_eta[i]);
					_g122 += (_ref->d2[1][s]) * (_ref->d1[1][i] * _ref->buf_eta[i]);
					_g222 += (_ref->d2[3][s]) * (_ref->d1[1][i] * _ref->buf_eta[i]);

					_g111 += (_ref->d1[0][s]) * (_ref->d2[0][i] * _ref->buf_eta[i]);
					_g211 += (_ref->d1[0][s]) * (_ref->d2[1][i] * _ref->buf_eta[i]);
					_g112 += (_ref->d1[0][s]) * (_ref->d2[1][i] * _ref->buf_eta[i]);
					_g212 += (_ref->d1[0][s]) * (_ref->d2[3][i] * _ref->buf_eta[i]);
					_g122 += (_ref->d1[1][s]) * (_ref->d2[1][i] * _ref->buf_eta[i]);
					_g222 += (_ref->d1[1][s]) * (_ref->d2[3][i] * _ref->buf_eta[i]);


					_e11 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[0][s]);
					_e12 += (_ref->d1[0][i] * _ref->buf_eta[i]) * (_ref->d1[1][s]);
					_e22 += (_ref->d1[1][i] * _ref->buf_eta[i]) * (_ref->d1[1][s]);

					_e11 += (_ref->d1[0][s]) * (_ref->d1[0][i] * _ref->buf_eta[i]);
					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][i] * _ref->buf_eta[i]);
					_e22 += (_ref->d1[1][s]) * (_ref->d1[1][i] * _ref->buf_eta[i]);
				}

				double _e21 = _e12;
				_g111 += -this->get_Gammaijk(0, 0, 0) * _e11 - this->get_Gammaijk(0, 0, 1) * _e12 - this->get_Gammaijk(0, 0, 0) * _e11 - this->get_Gammaijk(0, 0, 1) * _e12;
				_g211 += -this->get_Gammaijk(1, 0, 0) * _e11 - this->get_Gammaijk(1, 0, 1) * _e12 - this->get_Gammaijk(1, 0, 0) * _e11 - this->get_Gammaijk(1, 0, 1) * _e12;
				_g112 += -this->get_Gammaijk(0, 0, 0) * _e21 - this->get_Gammaijk(0, 0, 1) * _e22 - this->get_Gammaijk(0, 1, 0) * _e11 - this->get_Gammaijk(0, 1, 1) * _e12;
				_g212 += -this->get_Gammaijk(1, 0, 0) * _e21 - this->get_Gammaijk(1, 0, 1) * _e22 - this->get_Gammaijk(1, 1, 0) * _e11 - this->get_Gammaijk(1, 1, 1) * _e12;
				_g122 += -this->get_Gammaijk(0, 1, 0) * _e21 - this->get_Gammaijk(0, 1, 1) * _e22 - this->get_Gammaijk(0, 1, 0) * _e21 - this->get_Gammaijk(0, 1, 1) * _e22;
				_g222 += -this->get_Gammaijk(1, 1, 0) * _e21 - this->get_Gammaijk(1, 1, 1) * _e22 - this->get_Gammaijk(1, 1, 0) * _e21 - this->get_Gammaijk(1, 1, 1) * _e22;

				_g121 = _g112;
				_g221 = _g212;
				//val = (_g111 * t1 * v1 * s1 + _g112 * t1 * v1 * s2 + _g121 * t1 * v2 * s1 + _g122 * t1 * v2 * s2 +
				//	_g211 * t2 * v1 * s1 + _g212 * t2 * v1 * s2 + _g221 * t2 * v2 * s1 + _g222 * t2 * v2 * s2);

				val = ((_g111 * get_Gij(0, 0) * e12 + _g111 * get_Gij(0, 1) * e22 + _g112 * get_Gij(1, 0) * e12 + _g112 * get_Gij(1, 1) * e22) * t1 +
					(_g211 * get_Gij(0, 0) * e12 + _g211 * get_Gij(0, 1) * e22 + _g212 * get_Gij(1, 0) * e21 + _g212 * get_Gij(1, 1) * e22) * t2 -
					(_g121 * get_Gij(0, 0) * e11 + _g121 * get_Gij(0, 1) * e21 + _g122 * get_Gij(1, 0) * e11 + _g122 * get_Gij(1, 1) * e21) * t1 -
					(_g221 * get_Gij(0, 0) * e11 + _g221 * get_Gij(0, 1) * e21 + _g222 * get_Gij(1, 0) * e11 + _g222 * get_Gij(1, 1) * e21) * t2) / _ref->_refDv;

				val += ((g111 * get_Gij(0, 0) * _e12 + g111 * get_Gij(0, 1) * _e22 + g112 * get_Gij(1, 0) * _e12 + g112 * get_Gij(1, 1) * _e22) * t1 +
					(g211 * get_Gij(0, 0) * _e12 + g211 * get_Gij(0, 1) * _e22 + g212 * get_Gij(1, 0) * _e21 + g212 * get_Gij(1, 1) * _e22) * t2 -
					(g121 * get_Gij(0, 0) * _e11 + g121 * get_Gij(0, 1) * _e21 + g122 * get_Gij(1, 0) * _e11 + g122 * get_Gij(1, 1) * _e21) * t1 -
					(g221 * get_Gij(0, 0) * _e11 + g221 * get_Gij(0, 1) * _e21 + g222 * get_Gij(1, 0) * _e11 + g222 * get_Gij(1, 1) * _e21) * t2) / _ref->_refDv;
				*ptr1 = val;
				ptr1++;
			}

		}
		
		double guide2(double v1, double v2, bool accurate)
		{
			double val = 0;


			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;

			double k11 = S11 * _ref->get__Gij(0, 0) * S11 + S11 * _ref->get__Gij(0, 1) * S21 + S12 * _ref->get__Gij(1, 0) * S11 + S12 * _ref->get__Gij(1, 1) * S21;
			double k12 = S11 * _ref->get__Gij(0, 0) * S12 + S11 * _ref->get__Gij(0, 1) * S22 + S12 * _ref->get__Gij(1, 0) * S12 + S12 * _ref->get__Gij(1, 1) * S22;
			double k22 = S21 * _ref->get__Gij(0, 0) * S12 + S21 * _ref->get__Gij(0, 1) * S22 + S22 * _ref->get__Gij(1, 0) * S12 + S22 * _ref->get__Gij(1, 1) * S22;
			double k21 = k12;
			val = (k11 * v1 * v1 * _ref->get__gij(0, 1) + k11 * v1 * v2 * _ref->get__gij(1, 1) + k12 * v2 * v1 * _ref->get__gij(0, 1) + k12 * v2 * v2 * _ref->get__gij(1, 1)) / _ref->_refDv;
			val -= (k21 * v1 * v1 * _ref->get__gij(0, 0) + k21 * v1 * v2 * _ref->get__gij(1, 0) + k22 * v2 * v1 * _ref->get__gij(0, 0) + k22 * v2 * v2 * _ref->get__gij(1, 0)) / _ref->_refDv;
			return val / _ref->_refDv;
		}
		void guide2_z(double* ptr, double v1, double v2, bool accurate)
		{
			double val = 0;


			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;

			double k11 = S11 * _ref->get__Gij(0, 0) * S11 + S11 * _ref->get__Gij(0, 1) * S21 + S12 * _ref->get__Gij(1, 0) * S11 + S12 * _ref->get__Gij(1, 1) * S21;
			double k12 = S11 * _ref->get__Gij(0, 0) * S12 + S11 * _ref->get__Gij(0, 1) * S22 + S12 * _ref->get__Gij(1, 0) * S12 + S12 * _ref->get__Gij(1, 1) * S22;
			double k22 = S21 * _ref->get__Gij(0, 0) * S12 + S21 * _ref->get__Gij(0, 1) * S22 + S22 * _ref->get__Gij(1, 0) * S12 + S22 * _ref->get__Gij(1, 1) * S22;
			double k21 = k12;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = 0, _S12 = 0, _S22 = 0;
				_S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				_S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				_S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;

				double _k11 = _S11 * _ref->get__Gij(0, 0) * S11 + _S11 * _ref->get__Gij(0, 1) * S21 + _S12 * _ref->get__Gij(1, 0) * S11 + _S12 * _ref->get__Gij(1, 1) * S21;
				double _k12 = _S11 * _ref->get__Gij(0, 0) * S12 + _S11 * _ref->get__Gij(0, 1) * S22 + _S12 * _ref->get__Gij(1, 0) * S12 + _S12 * _ref->get__Gij(1, 1) * S22;
				double _k22 = _S21 * _ref->get__Gij(0, 0) * S12 + _S21 * _ref->get__Gij(0, 1) * S22 + _S22 * _ref->get__Gij(1, 0) * S12 + _S22 * _ref->get__Gij(1, 1) * S22;
				_k11 += S11 * _ref->get__Gij(0, 0) * _S11 + S11 * _ref->get__Gij(0, 1) * _S21 + S12 * _ref->get__Gij(1, 0) * _S11 + S12 * _ref->get__Gij(1, 1) * _S21;
				_k12 += S11 * _ref->get__Gij(0, 0) * _S12 + S11 * _ref->get__Gij(0, 1) * _S22 + S12 * _ref->get__Gij(1, 0) * _S12 + S12 * _ref->get__Gij(1, 1) * _S22;
				_k22 += S21 * _ref->get__Gij(0, 0) * _S12 + S21 * _ref->get__Gij(0, 1) * _S22 + S22 * _ref->get__Gij(1, 0) * _S12 + S22 * _ref->get__Gij(1, 1) * _S22;
				double _k21 = _k12;
				val = (_k11 * v1 * v1 * _ref->get__gij(0, 1) + _k11 * v1 * v2 * _ref->get__gij(1, 1) + _k12 * v2 * v1 * _ref->get__gij(0, 1) + _k12 * v2 * v2 * _ref->get__gij(1, 1)) / _ref->_refDv;
				val -= (_k21 * v1 * v1 * _ref->get__gij(0, 0) + _k21 * v1 * v2 * _ref->get__gij(1, 0) + _k22 * v2 * v1 * _ref->get__gij(0, 0) + _k22 * v2 * v2 * _ref->get__gij(1, 0)) / _ref->_refDv;

				*ptr1 = val / _ref->_refDv;
				ptr1++;
			}

		}
		double guide3(double v1, double v2, bool accurate)
		{
			double val = 0;


			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			val = S11 * v1 * v1 * _ref->get__gij(0, 1);
			val += S11 * v1 * v2 * _ref->get__gij(1, 1);
			val += S12 * v2 * v1 * _ref->get__gij(0, 1);
			val += S12 * v2 * v2 * _ref->get__gij(1, 1);
			val -= S21 * v1 * v1 * _ref->get__gij(0, 0);
			val -= S21 * v1 * v2 * _ref->get__gij(1, 0);
			val -= S22 * v2 * v1 * _ref->get__gij(0, 0);
			val -= S22 * v2 * v2 * _ref->get__gij(1, 0);


			return val;
		}
		void guide3_z(double* ptr, double v1, double v2, bool accurate)
		{
			double val = 0;


			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			double S11 = 0, S12 = 0, S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = 0, _S12 = 0, _S22 = 0;
				_S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				_S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				_S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;

				double S21 = S12;
				val = _S11 * v1 * v1 * _ref->get__gij(0, 1);
				val += _S11 * v1 * v2 * _ref->get__gij(1, 1);
				val += _S12 * v2 * v1 * _ref->get__gij(0, 1);
				val += _S12 * v2 * v2 * _ref->get__gij(1, 1);
				val -= _S21 * v1 * v1 * _ref->get__gij(0, 0);
				val -= _S21 * v1 * v2 * _ref->get__gij(1, 0);
				val -= _S22 * v2 * v1 * _ref->get__gij(0, 0);
				val -= _S22 * v2 * v2 * _ref->get__gij(1, 0);


				*ptr1 = val;
				ptr1++;
			}			

		}
		

		

		double guide_supported(double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;

			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double tr = (e11*_ref->get__Gij(0, 0) + e22* _ref->get__Gij(1, 1)+2*e12* _ref->get__Gij(0, 1));
			double tr2 = 1;// _ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			//scale = 1;
			//tr = 1;
			val = (e11 * v1 * v1+ e12 * v1 * v2 + e12 * v2 * v1 + e22 * v2 * v2)/tr;
			return val;
		}
		void guide_supported_xi(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;


			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;

			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double tr = (e11 * _ref->get__Gij(0, 0) + e22 * _ref->get__Gij(1, 1) + 2 * e12 * _ref->get__Gij(0, 1));

		
			double _e11 = e11, _e12 = e12, _e22 = e22;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				//scale = 1;
				double _tr = (e11 + e22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
				//tr = 1;

				val = (e11 * v1 * v1 + e12 * v1 * v2 + e12 * v2 * v1 + e22 * v2 * v2)/tr;
				//val += -(_e11 * v1 * v1 + _e12 * v1 * v2 + _e12 * v2 * v1 + _e22 * v2 * v2) / tr / tr * _tr;

				*ptr1 = val;
				ptr1++;
			}

		}

		void guide_supported_eta(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;


			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;

			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;

			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
		
				}
			}
			double tr = (e11 * _ref->get__Gij(0, 0) + e22 * _ref->get__Gij(1, 1) + 2 * e12 * _ref->get__Gij(0, 1));

			double _e11 = e11, _e12 = e12, _e22 = e22;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				//scale = 1;
				//tr = 1;

				double _tr = (e11 + e22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));

				val = (e11 * v1 * v1 + e12 * v1 * v2 + e12 * v2 * v1 + e22 * v2 * v2)/tr;
				//val += -(_e11 * v1 * v1 + _e12 * v1 * v2 + _e12 * v2 * v1 + _e22 * v2 * v2) / tr / tr * _tr;

				*ptr1 = val;
				ptr1++;
			}

		}
	

		double guide_free(double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;

			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;

			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double e21 = e12;
			double s1 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
			double s2 = -this->get_gij(0, 0) * v1 - this->get_gij(0, 1) * v2;
			double det = this->get_gij(0, 0) * this->get_gij(1, 1) - this->get_gij(0, 1) * this->get_gij(0, 1);
			double  tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));
			//tr = 1;
			val = (e11 * v1 * v1 + 2 * e12 * v1 * v2 + e22 * v2 * v2) / tr;// / _ref->_refDv;
			val -= (e11 * s1 * s1 + 2 * e12 * s1 * s2 + e22 * s2 * s2)/tr/det;

			return val;
		}
		void guide_free_xi(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;

			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;

			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				if (length == 0)length = 1;

				v1 = v1 / length;
				v2 = v2 / length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			//double tr = (e11+e22)/(_ref->get__Gij(0,0)+_ref->get__Gij(1,1));// e11* _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1);
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

			double tr2 = 1;//_ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double _e11 = e11, _e12 = e12, _e22 = e22,_e21=_e12;
			double* ptr1 = ptr;
			double det = this->get_gij(0, 0) * this->get_gij(1, 1) - this->get_gij(0, 1) * this->get_gij(1, 0);

			double s1 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
			double s2 = -this->get_gij(0, 0) * v1 - this->get_gij(0, 1) * v2;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				double e21 = e12;
				/*double g11 = 0, g12 = 0, g21 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				g21 = g12;*/


				//double _tr = (_e11 + _e22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
				double _tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

				//tr = 1;
				val = (e11 * v1 * v1 + 2 * e12 * v1 * v2 + e22 * v2 * v2) / tr;// / _ref->_refDv;
				val -= (e11 * s1 * s1 + 2 * e12 * s1 * s2 + e22 * s2 * s2) / tr / det;// / _ref->_refDv;
				//val += -(_e11 * v1 * v1 + 2 * _e12 * v1 * v2 + _e22 * v2 * v2) / tr/tr*_tr;// / _ref->_refDv;
				//val -= -(_e11 * s1 * s1 + 2 * _e12 * s1 * s2 + _e22 * s2 * s2) / tr/tr*_tr / det;//  / _ref->_refDv;


				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_free_z(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;


			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;
			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			//double tr = (e11 + e22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));// e11* _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1);
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));
			double tr2 = 1;//_ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}

			double s1 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
			double s2 = -this->get_gij(0, 0) * v1 - this->get_gij(0, 1) * v2;
			double det = this->get_gij(0, 0) * this->get_gij(1, 1) - this->get_gij(0, 1) * this->get_gij(0, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double g11 = 0, g12 = 0, g22 = 0, g21 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2*(_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				g21 = g12;
				double _s1 = g12 * v1 +g22 * v2;
				double _s2 = g11 * v1 -g12 * v2;
				double _det = this->get_gij(0, 0) *g22 - this->get_gij(0, 1) * g12+g11 * this->get_gij(1, 1) - g12 * this->get_gij(0, 1);
				
				tr = 1;
				val = -((e11 * _s1 * s1 + 2 * e12 * _s1 * s2 + e22 * _s2 * s2)/tr / det - (e11 * s1 * _s1 + 2 * e12 * s1 * _s2 + e22 * s2 * _s2))/tr / det;// / _ref->_refDv;
				
				//val -= -((e11 * s1 * s1 + 2 * e12 * s1 * s2 + e22 * s2 * s2))/tr / det/det*_det;// / _ref->_refDv;

				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_free_eta(double* ptr, double v1, double v2, bool accurate,int mode)
		{
			double val = 0;

			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			
			double length = 0;
			if (mode == 1) {
				double __mu = 0;
				double __nu = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					__mu += (_ref->d0[s]) * _ref->buf_mu[s];
					__nu += (_ref->d0[s]) * _ref->buf_nu[s];
				}
				v1 = __nu;
				v2 = -__mu;

			}
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e11 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_xi[s] * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += (_ref->d1[1][s]) * _ref->buf_eta[s] * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			//double tr = (e11 + e22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1)); // e11* _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1);
			double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));
			double tr2 = 1;// _ref->get__gij(0, 0) + _ref->get__gij(1, 1);
			double scale = tr2 / tr;
			if (tr == 0) {
				scale = 1;
				tr = 1;
			}
			double _e11 = e11, _e12 = e12, _e22 = e22;
			double* ptr1 = ptr;

			double s1 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
			double s2 = -this->get_gij(0, 0) * v1 - this->get_gij(0, 1) * v2;
			double det = this->get_gij(0, 0) * this->get_gij(1, 1) - this->get_gij(0, 1) * this->get_gij(0, 1);
			tr = 1;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				//double _tr = (e11 + e22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
				double _tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 1) + e22 * _ref->get__Gij(1, 1));

				//tr = 1;
				val = (e11 * v1 * v1 + 2 * e12 * v1 * v2 + e22 * v2 * v2) / tr;///_ref->_refDv;
				val -= (e11 * s1 * s1 + 2 * e12 * s1 * s2 + e22 * s2 * s2)/tr / det;/// _ref->_refDv;
				//val += -(_e11 * v1 * v1 + 2 * _e12 * v1 * v2 + _e22 * v2 * v2) / tr / tr * _tr;
				//val -= -(_e11 * s1 * s1 + 2 * _e12 * s1 * s2 + _e22 * s2 * s2) / tr / tr * _tr / det;

				*ptr1 = val;
				ptr1++;
			}

		}
		
		double mix2(double v1, double v2, double w1, double w2, bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			E11 = w1 * v1 * v1;
			E12 = w1 * v1 * v2;
			E22 = w1 * v2 * v2;

			E22 += this->get_gij2(0, 0) * v1 * v1 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * v1 * v2 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * v2 * v1 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * v2 * v2 * this->get_gij2(1, 0);
			E11 += this->get_gij2(1, 0) * v1 * v1 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * v1 * v2 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * v2 * v1 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * v2 * v2 * this->get_gij2(1, 1);
			E12 -= this->get_gij2(0, 0) * v1 * v1 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * v1 * v2 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * v2 * v1 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * v2 * v2 * this->get_gij2(1, 1);
			E21 = E12;

			double S21 = S12;
			double _s21 = _s12;

			double s11 = 0;	double s12 = 0;	double s22 = 0;	double s21 = 0;
			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc; s12 *= sc; s22 *= sc; s21 *= sc;
			double TRACE = 1;

			val = (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22);// / sc;
			val -= (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21);// / sc;

			return val;

		}
		void mix2_phi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}


			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
				
			}
			double _S21 = _S12;
			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;


			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc; s12 *= sc; s22 *= sc; s21 *= sc;
			double E11 = w1 * v1 * v1;
			double E12 = w1 * v1 * v2;
			double E22 = w1 * v2 * v2;

			E22 += this->get_gij2(0, 0) * v1 * v1 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * v1 * v2 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * v2 * v1 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * v2 * v2 * this->get_gij2(1, 0);
			E11 += this->get_gij2(1, 0) * v1 * v1 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * v1 * v2 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * v2 * v1 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * v2 * v2 * this->get_gij2(1, 1);
			E12 -= this->get_gij2(0, 0) * v1 * v1 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * v1 * v2 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * v2 * v1 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * v2 * v2 * this->get_gij2(1, 1);
			double E21 = E12;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double S11_z = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double S12_z = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double S21_z = S12_z;
				double S22_z = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _g11 = 0, _g12 = 0, _g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				double _g21 = _g12;
				double _E22 = _g11 * v1 * v1 * this->get_gij2(0, 0) + _g11 * v1 * v2 * this->get_gij2(1, 0) + _g12 * v2 * v1 * this->get_gij2(0, 0) + _g12 * v2 * v2 * this->get_gij2(1, 0);
				double _E11 = _g21 * v1 * v1 * this->get_gij2(0, 1) + _g21 * v1 * v2 * this->get_gij2(1, 1) + _g22 * v2 * v1 * this->get_gij2(0, 1) + _g22 * v2 * v2 * this->get_gij2(1, 1);
				double _E12 = -(_g22 * v1 * v1 * this->get_gij2(0, 1) + _g11 * v1 * v2 * this->get_gij2(1, 1) + _g12 * v2 * v1 * this->get_gij2(0, 1) + _g12 * v2 * v2 * this->get_gij2(1, 1));
				_E22 += this->get_gij2(0, 0) * v1 * v1 * _g11 + this->get_gij2(0, 0) * v1 * v2 * _g21 + this->get_gij2(0, 1) * v2 * v1 * _g11 + this->get_gij2(0, 1) * v2 * v2 * _g21;
				_E11 += this->get_gij2(1, 0) * v1 * v1 * _g12 + this->get_gij2(1, 0) * v1 * v2 * _g22 + this->get_gij2(1, 1) * v2 * v1 * _g12 + this->get_gij2(1, 1) * v2 * v2 * _g22;
				_E12 += -(this->get_gij2(0, 0) * v1 * v1 * _g12 + this->get_gij2(0, 0) * v1 * v2 * _g22 + this->get_gij2(0, 1) * v2 * v1 *_g12+ this->get_gij2(0, 1) * v2 * v2 * +_g22);

				double _E21 = _E12;

				val = (s11 * E11 * S12_z + s11 * E12 * S22_z + s12 * E21 * S12_z + s12 * E22 * S22_z);
				val -= (s21 * E11 * S11_z + s21 * E12 * S21_z + s22 * E21 * S11_z + s22 * E22 * S21_z);
				val += (s11 * _E11 * _S12 + s11 * _E12 * _S22 + s12 * _E21 * _S12 + s12 * _E22 * _S22);
				val -= (s21 * _E11 * _S11 + s21 * _E12 * _S21 + s22 * _E21 * _S11 + s22 * _E22 * _S21);
				*ptr1 = val;
				ptr1++;
			}
		}
		void mix2_Z(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
		
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
	
			}
			/**V1 = 1; V2 = 0;
			v1 = 1; v2 = 0;
			s1 = 0; s2 = 1;
			S1 = 0; S2 = 1;*/
		
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];

			}

			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;


			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc; s12 *= sc; s22 *= sc; s21 *= sc;
			double E11 = w1 * v1 * v1;
			double E12 = w1 * v1 * v2;
			double E22 = w1 * v2 * v2;

			E22 += this->get_gij2(0, 0) * v1 * v1 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * v1 * v2 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * v2 * v1 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * v2 * v2 * this->get_gij2(1, 0);
			E11 += this->get_gij2(1, 0) * v1 * v1 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * v1 * v2 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * v2 * v1 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * v2 * v2 * this->get_gij2(1, 1);
			E12 -= this->get_gij2(0, 0) * v1 * v1 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * v1 * v2 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * v2 * v1 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * v2 * v2 * this->get_gij2(1, 1);
			double E21 = E12;


			double tr = (E11 + E22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
			if (tr == 0)tr = 1;

			double _S21 = _S12;
			//double TRACE = (_ref->get__gij(0, 0) + _ref->get__gij(1, 1));
			double TRACE = 1;
			//if (TRACE == 0)TRACE = 1;
			/*
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double __s22_phi = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double __s12_phi = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double __s21_phi = __s12_phi;
				double __s11_phi = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				__s11_phi *= sc; __s12_phi *= sc; __s21_phi *= sc; __s22_phi *= sc;
				double s11_phi = 0, s12_phi = 0, s21_phi = 0, s22_phi = 0;

				if (accurate)
				{
					s11_phi = this->get_gij2(0, 0) * __s11_phi * this->get_gij2(0, 0) + this->get_gij2(0, 0) * __s12_phi * this->get_gij2(1, 0) + this->get_gij2(0, 1) * __s21_phi * this->get_gij2(0, 0) + this->get_gij2(0, 1) * __s22_phi * this->get_gij2(1, 0);
					s12_phi = this->get_gij2(0, 0) * __s11_phi * this->get_gij2(0, 1) + this->get_gij2(0, 0) * __s12_phi * this->get_gij2(1, 1) + this->get_gij2(0, 1) * __s21_phi * this->get_gij2(0, 1) + this->get_gij2(0, 1) * __s22_phi * this->get_gij2(1, 1);
					s22_phi = this->get_gij2(1, 0) * __s11_phi * this->get_gij2(0, 1) + this->get_gij2(1, 0) * __s12_phi * this->get_gij2(1, 1) + this->get_gij2(1, 1) * __s21_phi * this->get_gij2(0, 1) + this->get_gij2(1, 1) * __s22_phi * this->get_gij2(1, 1);
					s21_phi = s12_phi;
				}
				else {
					s11_phi = _ref->get__gij(0, 0) * __s11_phi * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * __s12_phi * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * __s21_phi * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * __s22_phi * _ref->get__gij(1, 0);
					s12_phi = _ref->get__gij(0, 0) * __s11_phi * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * __s12_phi * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * __s21_phi * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * __s22_phi * _ref->get__gij(1, 1);
					s22_phi = _ref->get__gij(1, 0) * __s11_phi * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * __s12_phi * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * __s21_phi * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * __s22_phi * _ref->get__gij(1, 1);
					s21_phi = s12_phi;
				}
				//val += w2 * ((s12 * v1 + s22 * v2) * (-(_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * v1 - (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * v2));
				//val -= w2 * ((-s11 * v1 - s12 * v2) * ((_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * v1 + (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * v2));
				//val += w1 * ((s21 * V1 + s22 * V2) * (-S11 * v1 - S12 * v2) - (-s11 * V1 - s12 * V2) * (S21 * v1 + S22 * v2));
				//val += w2 * ((s21 * S1 + s22 * S2) * (-S11 * s1 - S12 * s2) - (-s11 * S1 - s12 * S2) * (S21 * s1 + S22 * s2));


				//TRACE = 1;
				val =  (s11_phi * E11 * _S12 + s11_phi * E12 * _S22 + s12_phi * E21 * _S12 + s12_phi * E22 * _S22);// / sc;
				val -= (s21_phi * E11 * _S11 + s21_phi * E12 * _S21 + s22_phi * E21 * _S11 + s22_phi * E22 * _S21);// / sc;

				*ptr1 = val;
				ptr1++;
			}

		}
		double e11(double cx,double cy,double L)
		{
			double xi = 0, eta = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
					xi += _ref->d0[s] * _ref->buf_xi[s];
					eta += _ref->d0[s] * _ref->buf_eta[s];

				
			}
			return (xi - cx) * (xi - cx) + (eta - cy) * (eta - cy) - L * L;
		}
		void e11_xi(double cx, double cy, double L,double *ptr)
		{
			double xi = 0, eta = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				xi += _ref->d0[s] * _ref->buf_xi[s];
				eta += _ref->d0[s] * _ref->buf_eta[s];


			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi = _ref->d0[s] ;
				double val = 2*(xi - cx) * _xi;
				*ptr1 = val;
				ptr1++;
			}

		}
		void e11_eta(double cx, double cy, double L, double* ptr)
		{
			double xi = 0, eta = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				xi += _ref->d0[s] * _ref->buf_xi[s];
				eta += _ref->d0[s] * _ref->buf_eta[s];


			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _eta = _ref->d0[s] ;
				double val = 2 * (eta - cy) * _eta;
				*ptr1 = val;
				ptr1++;
			}

		}
	
		double mix(double v1, double v2, double w1, double w2, bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0,E21;


			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;

			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;

			}

			double e1x=0, e1y=0, e2x=0, e2y=0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
	
				}
			}


			E11 = e22;
			E22 = e11;
			E12 = -e12;
			E21 = E12;


			//double tr = (E11 + E22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
			//double tr = (e11 * _ref->get__Gij(0, 0) + 2 * e12 * _ref->get__Gij(0, 2) + e22 * _ref->get__Gij(1, 1));
			double tr = (E11 * _ref->get__gij(0, 0) + 2 * E12 * _ref->get__gij(0, 1) + E22 * _ref->get__gij(1, 1));



			double S21 = S12;
			double _s21 = _s12;

			double s11 = 0;	double s12 = 0;	double s22 = 0;	double s21 = 0;
			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}

			double TRACE = 1;

			//tr = 1;
			val = 1. / tr / TRACE * (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22)*sc / _ref->_refDv;
			val -= 1. / tr / TRACE * (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21) * sc / _ref->_refDv;
			
			return val;

		}
		void remove3(int N, double* __ptr, double* __ptr2, double* __ptr3 )
		{
			Eigen::MatrixXd J(2, N);
			Eigen::VectorXd v(N);
			for (int i = 0; i < N; i++)
			{
				J(0, i) = __ptr2[i];
				J(1, i) = __ptr3[i];
				v(i) = __ptr[i];
			}
			v=v-J.transpose()*(J * J.transpose()).inverse()* J* v;
			for (int i = 0; i < N; i++)
			{
				__ptr[i] = v(i);
			}

		}

		void remove2(int N, double* __ptr, double* __ptr2, long long* index, double sc)
		{
		
			/*double __dot = 0, __norm = 0;
			double* ptr1 = __ptr;
			double* ptr2 = __ptr2;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				__dot += *ptr1 * *ptr2;
				__norm += *ptr2 * *ptr2;
				ptr1++;
				ptr2++;
			}
			ptr1 = __ptr3;
			ptr2 = __ptr4;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				__dot += *ptr1 * *ptr2;
				__norm += *ptr2 * *ptr2;
				ptr1++;
				ptr2++;
			}

			ptr2 = __ptr2;
			ptr1 = __ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				*ptr1 -= *ptr2 * __dot / __norm;
				ptr1++;
				ptr2++;
			}
			ptr2 = __ptr4;
			ptr1 = __ptr3;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				*ptr1 -= *ptr2 * __dot / __norm;
				ptr1++;
				ptr2++;
			}*/

		}
		void mix_phi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}

			}
			double _S21 = _S12;
	
			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;
		

			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}

			double E11 = e22;
			double E22 = e11;
			double E12 = -e12;
			double E21 = E12;


			
			
			//double tr = (E11 + E22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
			double tr = (E11 * _ref->get__gij(0, 0) + 2 * E12 * _ref->get__gij(0, 1) + E22 * _ref->get__gij(1, 1));

			double TRACE = 1;

			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double S11_z = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double S12_z = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double S21_z = S12_z;
				double S22_z = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _g11 = 0, _g12 = 0, _g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				double _g21 = _g12;



				double s11_z = _g11 * _s11 * this->get_gij2(0, 0) + _g11 * _s12 * this->get_gij2(1, 0) + _g12 * _s21 * this->get_gij2(0, 0) + _g12 * _s22 * this->get_gij2(1, 0);
				double s12_z = _g11 * _s11 * this->get_gij2(0, 1) + _g11 * _s12 * this->get_gij2(1, 1) + _g12 * _s21 * this->get_gij2(0, 1) + _g12 * _s22 * this->get_gij2(1, 1);
				double s22_z = _g21 * _s11 * this->get_gij2(0, 1) + _g21 * _s12 * this->get_gij2(1, 1) + _g22 * _s21 * this->get_gij2(0, 1) + _g22 * _s22 * this->get_gij2(1, 1);
				s11_z += this->get_gij2(0, 0) * _s11 * _g11 + this->get_gij2(0, 0) * _s12 * _g12 + this->get_gij2(0, 1) * _s21 * _g11 + this->get_gij2(0, 1) * _s22 * _g12;
				s12_z += this->get_gij2(0, 0) * _s11 * _g12 + this->get_gij2(0, 0) * _s12 * _g22 + this->get_gij2(0, 1) * _s21 * _g12 + this->get_gij2(0, 1) * _s22 * _g22;
				s22_z += this->get_gij2(1, 0) * _s11 * _g12 + this->get_gij2(1, 0) * _s12 * _g22 + this->get_gij2(1, 1) * _s21 *_g12 + this->get_gij2(1, 1) * _s22 * _g22;

				double s21_z = s12_z;
				//tr = 1;
				val = 1. / tr / TRACE * (s11 * E11 * S12_z + s11 * E12 * S22_z + s12 * E21 * S12_z + s12 * E22 * S22_z) * sc / _ref->_refDv;;// / sc;
				val -= 1. / tr / TRACE * (s21 * E11 * S11_z + s21 * E12 * S21_z + s22 * E21 * S11_z + s22 * E22 * S21_z) * sc / _ref->_refDv;// / sc;
				//val += 1. / tr / TRACE * (s11_z * E11 * _S12 + s11_z * E12 * _S22 + s12_z * E21 * _S12 + s12_z * E22 * _S22) * sc / _ref->_refDv;// / sc;
				//val -= 1. / tr / TRACE * (s21_z * E11 * _S11 + s21_z * E12 * _S21 + s22_z * E21 * _S11 + s22_z * E22 * _S21) * sc / _ref->_refDv;// / sc;


				*ptr1 = val;
				ptr1++;
			}
			
		}
		void mix_Z(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double length = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
	
			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;


			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}

			double E11 = e22;
			double E22 = e11;
			double E12 = -e12;
			double E21 = E12;


			//double tr = (E11 + E22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
			double tr = (E11 * _ref->get__gij(0, 0) + 2 * E12 * _ref->get__gij(0,1) + E22 * _ref->get__gij(1, 1));

			if (tr == 0)tr = 1;

			double _S21 = _S12;

			double TRACE = 1;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double __s22_phi = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double __s12_phi = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double __s21_phi = __s12_phi;
				double __s11_phi = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				//__s11_phi *= sc; __s12_phi *= sc; __s21_phi *= sc; __s22_phi *= sc;
				double s11_phi = 0, s12_phi = 0, s21_phi = 0, s22_phi = 0;

				if (accurate)
				{
					s11_phi = this->get_gij2(0, 0) * __s11_phi * this->get_gij2(0, 0) + this->get_gij2(0, 0) * __s12_phi * this->get_gij2(1, 0) + this->get_gij2(0, 1) * __s21_phi * this->get_gij2(0, 0) + this->get_gij2(0, 1) * __s22_phi * this->get_gij2(1, 0);
					s12_phi = this->get_gij2(0, 0) * __s11_phi * this->get_gij2(0, 1) + this->get_gij2(0, 0) * __s12_phi * this->get_gij2(1, 1) + this->get_gij2(0, 1) * __s21_phi * this->get_gij2(0, 1) + this->get_gij2(0, 1) * __s22_phi * this->get_gij2(1, 1);
					s22_phi = this->get_gij2(1, 0) * __s11_phi * this->get_gij2(0, 1) + this->get_gij2(1, 0) * __s12_phi * this->get_gij2(1, 1) + this->get_gij2(1, 1) * __s21_phi * this->get_gij2(0, 1) + this->get_gij2(1, 1) * __s22_phi * this->get_gij2(1, 1);
					s21_phi = s12_phi;
				}
				else {
					s11_phi = _ref->get__gij(0, 0) * __s11_phi * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * __s12_phi * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * __s21_phi * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * __s22_phi * _ref->get__gij(1, 0);
					s12_phi = _ref->get__gij(0, 0) * __s11_phi * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * __s12_phi * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * __s21_phi * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * __s22_phi * _ref->get__gij(1, 1);
					s22_phi = _ref->get__gij(1, 0) * __s11_phi * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * __s12_phi * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * __s21_phi * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * __s22_phi * _ref->get__gij(1, 1);
					s21_phi = s12_phi;
				}
				//tr = 1;
				val = 1. / tr / TRACE * (s11_phi * E11 * _S12 + s11_phi * E12 * _S22 + s12_phi * E21 * _S12 + s12_phi * E22 * _S22) * sc / _ref->_refDv;
				val -= 1. / tr / TRACE * (s21_phi * E11 * _S11 + s21_phi * E12 * _S21 + s22_phi * E21 * _S11 + s22_phi * E22 * _S21) * sc / _ref->_refDv;

				*ptr1 =val;
				ptr1++;
			}
			
		}
		void __U_z(double* ptr, bool accurate)
		{
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double _S11 = 0, _S12 = 0, _S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];

			}
			double _s21 = _s12;


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s22 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s11 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);



				double val = (_S11 * s11 + _S22 * s22 + 2 * _S12 * s12) * sc;

				if (accurate)
				{
		

					val = val;
					
					*ptr1 = val;
				}
				else {
					*ptr1 = val;
				}
				ptr1++;
			}
		}
		void __U_phi(double *ptr,double load,bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double _S11 = 0, _S12 = 0, _S22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];

			}
			double _s21 = _s12;

		
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);



				double val = (S11 * _s11 + S22 * _s22 + 2 * S12 * _s12)*sc;
				
				if (accurate)
				{
					double _g11 = 0, _g12 = 0, _g22 = 0;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					}
					double _g21 = _g12;
					//__mem->sc * __mem->bodyF-load* __mem->dv/ __mem->_ref->_refDv;
					
					//double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					//val += -load*ddv / _ref->_refDv;
					*ptr1 = val; 
				}
				else {
					*ptr1 = val;
				}
				ptr1++;
			}
		}
		void mix_xi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
		
			double e1x = 0, e1y = 0, e2x = 0, e2y = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
	
				}
			}
			e11 = e1x * e1x + e1y * e1y;
			e12 = e1x * e2x + e1y * e2y;
			e22 = e2x * e2x + e2y * e2y;
			double E11 = e22;
			double E22 = e11;
			double E12 = -e12;
			double E21 = E12;

			//double tr = (E11 + E22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
			double tr = (E11 * _ref->get__gij(0, 0) + 2 * E12 * _ref->get__gij(0, 1) + E22 * _ref->get__gij(1, 1));

			if (tr == 0) {
				
				tr = 1;
			}
			double S21 = S12;
			double _s21 = _s12;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;

			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}

			double TRACE = 1;
			double* ptr1 = ptr;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val;
				double e11_xi = 0, e12_xi = 0, e22_xi = 0;

				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11_xi += 2*(_ref->d1[0][s])  * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12_xi += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12_xi += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22_xi += 2*(_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				
				double _E11_xi = e22_xi, _E12_xi = -e12_xi, _E22_xi = e11_xi, _E21_xi = _E12_xi;
				//TRACE = 1;
				//tr = 1;
				val = 1./tr /TRACE* (s11 * _E11_xi * S12 + s11 * _E12_xi * S22 + s12 * _E21_xi * S12 + s12 * _E22_xi * S22) * sc / _ref->_refDv;;
				val -=1./tr/TRACE * (s21 * _E11_xi * S11 + s21 * _E12_xi * S21 + s22 * _E21_xi * S11 + s22 * _E22_xi * S21) * sc / _ref->_refDv;;
				//double dtr = (_E11_xi + _E22_xi) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
				double dtr = (_E11_xi * _ref->get__gij(0, 0) + 2 * _E12_xi * _ref->get__gij(0, 1) + _E22_xi * _ref->get__gij(1, 1));

				val += -1. / tr / tr /TRACE* (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22) * sc*dtr / _ref->_refDv;
				val -= -1. / tr / tr /TRACE* (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21) *sc* dtr / _ref->_refDv;
				
				
				

				
				*ptr1 = val;
				ptr1++;
			}
		}
		
		void mix_eta( double *ptr,double v1, double v2, double w1, double w2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double e1x = 0, e1y = 0, e2x = 0, e2y = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			e11 = e1x * e1x + e1y * e1y;
			e12 = e1x * e2x + e1y * e2y;
			e22 = e2x * e2x + e2y * e2y;
			double E11 = e22;
			double E22 = e11;
			double E12 = -e12;
			double E21 = E12;



			//double tr = (E11 + E22) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
			double tr = (E11 * _ref->get__gij(0, 0) + 2 * E12 * _ref->get__gij(0, 1) + E22 * _ref->get__gij(1, 1));
			if (tr == 0) {
				
				tr = 1;
			}
			double S21 = S12;
			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;
			

			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}


			double TRACE = 1;
			for(int s = 0; s < _ref->_nNode; s++)
			{
				double val;
				double e11_eta = 0, e12_eta = 0, e22_eta = 0;


				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11_eta += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12_eta += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12_eta += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22_eta += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				double _E11_eta = e22_eta, _E12_eta = -e12_eta, _E22_eta = e11_eta, _E21_eta = _E12_eta;
				//tr = 1;

				val = 1. / tr / TRACE * (s11 * _E11_eta * S12 + s11 * _E12_eta * S22 + s12 * _E21_eta * S12 + s12 * _E22_eta * S22) * sc / _ref->_refDv;// / sc;
				val -= 1. / tr / TRACE * (s21 * _E11_eta * S11 + s21 * _E12_eta * S21 + s22 * _E21_eta * S11 + s22 * _E22_eta * S21) * sc / _ref->_refDv;// / sc;

				//double dtr = (_E11_eta + _E22_eta) / (_ref->get__Gij(0, 0) + _ref->get__Gij(1, 1));
				double dtr = (_E11_eta * _ref->get__gij(0, 0) + 2 * _E12_eta * _ref->get__gij(0, 1) + _E22_eta * _ref->get__gij(1, 1));

				val += -1. / tr / tr / TRACE * (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22) *sc* dtr / _ref->_refDv;// / sc;
				val -= -1. / tr / tr / TRACE * (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21) *sc* dtr / _ref->_refDv;// / sc;

				*ptr1 = val;
				ptr1++;
			}
		}
	

		double mix_guide(double v1, double v2, double w1, double w2, bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			double E11 = 0, E12 = 0, E22 = 0, E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			double e1x = 0, e1y = 0, e2x = 0, e2y = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}

			/*for (int s = 0; s < _ref->_nNode; s++)
			{
				s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}*/

			

			double S21 = S12;
			double _s21 = _s12;

			double s11 = 0;	double s12 = 0;	double s22 = 0;	double s21 = 0;
			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc; s12 *= sc; s22 *= sc; s21 *= sc;

			E11 = w1 * v1 * v1 + w2 * s1 * s1;
			E12 = w1 * v1 * v2 + w2 * s1 * s2;
			E21 = E12;
			E22 = w1 * v2 * v2 + w2 * s2 * s2;
			val =  (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22);// / sc;
			val -=  (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21);// / sc;

			return val;

		}
		void mix_guide_phi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double e11 = 0, e12 = 0, e22 = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];

				//}
			}
			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;


			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc; s12 *= sc; s22 *= sc; s21 *= sc;

			double E11 = w1 * v1 * v1 + w2 * s1 * s1;
			double E12 = w1 * v1 * v2 + w2 * s1 * s2;
			double E21 = E12;
			double E22 = w1 * v2 * v2 + w2 * s2 * s2;



			double _S21 = _S12;
			double TRACE = 1;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double S11_z = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double S12_z = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double S21_z = S12_z;
				double S22_z = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _g11 = 0, _g12 = 0, _g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				double _g21 = _g12;
				double s11_z = _g11 * _s11 * this->get_gij2(0, 0) + _g11 * _s12 * this->get_gij2(1, 0) + _g12 * _s21 * this->get_gij2(0, 0) + _g12 * _s22 * this->get_gij2(1, 0);
				double s12_z = _g11 * _s11 * this->get_gij2(0, 1) + _g11 * _s12 * this->get_gij2(1, 1) + _g12 * _s21 * this->get_gij2(0, 1) + _g12 * _s22 * this->get_gij2(1, 1);
				double s22_z = _g21 * _s11 * this->get_gij2(0, 1) + _g21 * _s12 * this->get_gij2(1, 1) + _g22 * _s21 * this->get_gij2(0, 1) + _g22 * _s22 * this->get_gij2(1, 1);
				s11_z += this->get_gij2(0, 0) * _s11 * _g11 + this->get_gij2(0, 0) * _s12 * _g12 + this->get_gij2(0, 1) * _s21 * _g11 + this->get_gij2(0, 1) * _s22 * _g12;
				s12_z += this->get_gij2(0, 0) * _s11 * _g12 + this->get_gij2(0, 0) * _s12 * _g22 + this->get_gij2(0, 1) * _s21 * _g12 + this->get_gij2(0, 1) * _s22 * _g22;
				s22_z += this->get_gij2(1, 0) * _s11 * _g12 + this->get_gij2(1, 0) * _s12 * _g22 + this->get_gij2(1, 1) * _s21 * _g12 + this->get_gij2(1, 1) * _s22 * _g22;

				double s21_z = s12_z;
				s11_z *= sc; s12_z *= sc; s21_z *= sc; s22_z *= sc;

				val = (s11 * E11 * S12_z + s11 * E12 * S22_z + s12 * E21 * S12_z + s12 * E22 * S22_z);
				val -= (s21 * E11 * S11_z + s21 * E12 * S21_z + s22 * E21 * S11_z + s22 * E22 * S21_z);
				val += (s11_z * E11 * _S12 + s11_z * E12 * _S22 + s12_z * E21 * _S12 + s12_z * E22 * _S22);
				val -= (s21_z * E11 * _S11 + s21_z * E12 * _S21 + s22_z * E21 * _S11 + s22_z * E22 * _S21);


				*ptr1 = val;
				ptr1++;
			}

		}
		void mix_guide_Z(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			double _S11 = 0;
			double _S12 = 0;
			double _S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			/**V1 = 1; V2 = 0;
			v1 = 1; v2 = 0;
			s1 = 0; s2 = 1;
			S1 = 0; S2 = 1;*/
			double e1x = 0, e1y = 0, e2x = 0, e2y = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				_S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}

			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;


			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc; s12 *= sc; s22 *= sc; s21 *= sc;
			double E11 = w1 * v1 * v1 + w2 * s1 * s1;
			double E12 = w1 * v1 * v2 + w2 * s1 * s2;
			double E21 = E12;
			double E22 = w1 * v2 * v2 + w2 * s2 * s2;




			double _S21 = _S12;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double __s22_phi = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double __s12_phi = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double __s21_phi = __s12_phi;
				double __s11_phi = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				__s11_phi *= sc; __s12_phi *= sc; __s21_phi *= sc; __s22_phi *= sc;
				double s11_phi = 0, s12_phi = 0, s21_phi = 0, s22_phi = 0;

				if (accurate)
				{
					s11_phi = this->get_gij2(0, 0) * __s11_phi * this->get_gij2(0, 0) + this->get_gij2(0, 0) * __s12_phi * this->get_gij2(1, 0) + this->get_gij2(0, 1) * __s21_phi * this->get_gij2(0, 0) + this->get_gij2(0, 1) * __s22_phi * this->get_gij2(1, 0);
					s12_phi = this->get_gij2(0, 0) * __s11_phi * this->get_gij2(0, 1) + this->get_gij2(0, 0) * __s12_phi * this->get_gij2(1, 1) + this->get_gij2(0, 1) * __s21_phi * this->get_gij2(0, 1) + this->get_gij2(0, 1) * __s22_phi * this->get_gij2(1, 1);
					s22_phi = this->get_gij2(1, 0) * __s11_phi * this->get_gij2(0, 1) + this->get_gij2(1, 0) * __s12_phi * this->get_gij2(1, 1) + this->get_gij2(1, 1) * __s21_phi * this->get_gij2(0, 1) + this->get_gij2(1, 1) * __s22_phi * this->get_gij2(1, 1);
					s21_phi = s12_phi;
				}
				else {
					s11_phi = _ref->get__gij(0, 0) * __s11_phi * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * __s12_phi * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * __s21_phi * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * __s22_phi * _ref->get__gij(1, 0);
					s12_phi = _ref->get__gij(0, 0) * __s11_phi * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * __s12_phi * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * __s21_phi * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * __s22_phi * _ref->get__gij(1, 1);
					s22_phi = _ref->get__gij(1, 0) * __s11_phi * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * __s12_phi * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * __s21_phi * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * __s22_phi * _ref->get__gij(1, 1);
					s21_phi = s12_phi;
				}
				val =  (s11_phi * E11 * _S12 + s11_phi * E12 * _S22 + s12_phi * E21 * _S12 + s12_phi * E22 * _S22);// / sc;
				val -= (s21_phi * E11 * _S11 + s21_phi * E12 * _S21 + s22_phi * E21 * _S11 + s22_phi * E22 * _S21);// / sc;

				*ptr1 = val;
				ptr1++;
			}

		}
		double E22(double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double tr2 = 1;
			double tr = e11 + e22;
			double scale = tr2 / tr;
			e21 = e12;

			val =  scale*(e11 * s1 * s1 + e12 * s1 * s2 + e21 * s2 * s1 + e22 * s2 * s2);
			return val;
		}
		void E22_xi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double tr2 = 1;
			double tr = _e11 + _e22;
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				e21 = e12;
				val = scale * ((e11 * s1 * s1 + e12 * s1 * s2 + e21 * s2 * s1 + e22 * s2 * s2));
				val += -tr2/tr/tr * ((_e11 * s1 * s1 + _e12 * s1 * s2 + _e21 * s2 * s1 + _e22 * s2 * s2))*(e11+e22);

				*ptr1 = val;
				ptr1++;
			}
		}
		void E22_eta(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double tr2 = 1;
			double tr = _e11 + _e22;
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				e21 = e12;
				val = scale * ((e11 * s1 * s1 + e12 * s1 * s2 + e21 * s2 * s1 + e22 * s2 * s2));
				val += -tr2 / tr / tr * ((_e11 * s1 * s1 + _e12 * s1 * s2 + _e21 * s2 * s1 + _e22 * s2 * s2)) * (e11 + e22);

				*ptr1 = val;
				ptr1++;
			}
		}
		double testequal(double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			double tr2 = 1;
			double tr = e11 + e22;
			double scale = tr2 / tr;
			e21 = e12;
			scale = 1;
			val = scale * ((e11 * v1 * v1 + e12 * v1 * v2 + e21 * v2 * v1 + e22 * v2 * v2) - (e11 * s1 * s1 + e12 * s1 * s2 + e21 * s2 * s1 + e22 * s2 * s2));
			return val;
		}
		void testequal_xi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double tr2 = 1;
			double tr = _e11 + _e22;
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				e21 = e12;
				scale = 1;
				val = scale * ((e11 * v1 * v1 + e12 * v1 * v2 + e21 * v2 * v1 + e22 * v2 * v2) - (e11 * s1 * s1 + e12 * s1 * s2 + e21 * s2 * s1 + e22 * s2 * s2));
				//val += -tr2/tr/tr* ((_e11 * v1 * v1 + _e12 * v1 * v2 + _e21 * v2 * v1 + _e22 * v2 * v2) - (_e11 * s1 * s1 + _e12 * s1 * s2 + _e21 * s2 * s1 + _e22 * s2 * s2))*(e11+e22);

				*ptr1 = val;
				ptr1++;
			}
		}
		void testequal_eta(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double tr2 = 1;
			double tr = _e11 + _e22;
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				e21 = e12;
				scale = 1;
				val = scale * ((e11 * v1 * v1 + e12 * v1 * v2 + e21 * v2 * v1 + e22 * v2 * v2) - (e11 * s1 * s1 + e12 * s1 * s2 + e21 * s2 * s1 + e22 * s2 * s2));
				//val += -tr2 / tr / tr * ((_e11 * v1 * v1 + _e12 * v1 * v2 + _e21 * v2 * v1 + _e22 * v2 * v2) - (_e11 * s1 * s1 + _e12 * s1 * s2 + _e21 * s2 * s1 + _e22 * s2 * s2)) * (e11 + e22);

				*ptr1 = val;
				ptr1++;
			}
		}
		double setEij_A(double v1, double v2, double w1, double w2, bool accurate)
		{
			
			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;
			double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
			}
			/*double E11 = e22, E22 = e11, E12 = -e12, E21 = E12;
			double k11 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 0);
			double k12 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 1);
			double k22 = this->get_gij2(1, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * E22 * this->get_gij2(1, 1);
			double k21 = k12;
			*/
			double tr2 = 1;
			double tr = e11 + e22;
			double scale = tr2 / tr;
			e21 = e12;
			val = scale * ((e11 * v1 * v1) + (e12 * v1 * v2) + (e21 * v2 * v1) + (e22 * v2 * v2)) * w1;
 		    val-= scale *  ((e11 * s1 * s1) + (e12 * s1 * s2) + (e21 * s2 * s1) + (e22 * s2 * s2))*w2;
			return val;
		}
		void setEij_A_xi(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}

			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			/*double E11 = _e22, E22 = _e11, E12 = -_e12, E21 = E12;
			double _k11 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 0);
			double _k12 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 1);
			double _k22 = this->get_gij2(1, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * E22 * this->get_gij2(1, 1);
			double _k21 = _k12;*/

			double tr2 = 1;
			double tr = _e11 + _e22;
			
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				e21 = e12;
				/*double E11 = _e22, E22 = _e11, E12 = -_e12, E21 = E12;
				double k11 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 0);
				double k12 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 1);
				double k22 = this->get_gij2(1, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * E22 * this->get_gij2(1, 1);
				double k21 = k12;*/
				val = scale * ((e11 * v1 * v1) + (e12 * v1 * v2) + (e21 * v2 * v1) + (e22 * v2 * v2)) * w1;
				val += -scale * ((e11 * s1 * s1) + (e12 * s1 * s2) + (e21 * s2 * s1) + (e22 * s2 * s2))*w2;
				val += -w1 * tr2 / tr / tr * ((_e11 * v1 * v1) + (_e12 * v1 * v2) + (_e21 * v2 * v1) + (_e22 * v2 * v2)) * (e11 + e22);
				val += w2 *tr2 / tr / tr * ((_e11 * s1 * s1) + (_e12 * s1 * s2) + (_e21 * s2 * s1) + (_e22 * s2 * s2)) * (e11 + e22);


				*ptr1 = val;
				ptr1++;
			}


		}
		void  setEij_A_eta(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;

			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			/*double E11 = _e22, E22 = _e11, E12 = -_e12, E21 = E12;
			double _k11 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 0);
			double _k12 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 1);
			double _k22 = this->get_gij2(1, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * E21 * this->get_gij2(0, 1) + this->get_gij(21, 1) * E22 * this->get_gij2(1, 1);
			double _k21 = _k12;
			*/
			double tr2 = 1;
			double tr = _e11 + _e22;
			
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				e21 = e12;
				/*double E11 = _e22, E22 = _e11, E12 = -_e12, E21 = E12;
				double k11 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 0);
				double k12 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 1);
				double k22 = this->get_gij2(1, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * E22 * this->get_gij2(1, 1);
				double k21 = k12;
				*/
				val = scale * ((e11 * v1 * v1) + (e12 * v1 * v2) + (e21 * v2 * v1) + (e22 * v2 * v2)) * w1;
				val += -w2 * scale * ((e11 * s1 * s1) + (e12 * s1 * s2) + (e21 * s2 * s1) + (e22 * s2 * s2));
				val += -w1 *  tr2 / tr / tr * ((_e11 * v1 * v1) + (_e12 * v1 * v2) + (_e21 * v2 * v1) + (_e22 * v2 * v2)) * (e11 + e22);
				val += w2 * tr2 / tr / tr * ((_e11 * s1 * s1) + (_e12 * s1 * s2) + (_e21 * s2 * s1) + (_e22 * s2 * s2)) * (e11 + e22);
				*ptr1 = val;
				ptr1++;
			}

		}
		void  setEij_A_Z(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{

			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
			}
			s1 = -this->get_gij2(1, 0) * v1 - this->get_gij2(1, 1) * v2;
			s2 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
			double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
			double sdet = sqrt(det);
			s1 /= sdet;
			s2 /= sdet;
			double _s1 = s1;
			double _s2 = s2;
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double E11 = _e22, E22 = _e11, E12 = -_e12, E21 = E12;
			double _k11 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 0) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 0) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 0);
			double _k12 = this->get_gij2(0, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * E22 * this->get_gij2(1, 1);
			double _k22 = this->get_gij2(1, 0) * E11 * this->get_gij2(0, 1) + this->get_gij2(1, 0) * E12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * E21 * this->get_gij2(0, 1) + this->get_gij2(1, 1) * E22 * this->get_gij2(1, 1);
			double _k21 = _k12;

			double tr2 = 1;
			double tr = _e11 + _e22;

			double scale = tr2 / tr;
			double* ptr1 = ptr;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				val = 0;
				double g11 = 0, g12 = 0, g21 = 0, g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];
					g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}
				g21 = g12;
				/*double E11 = _e22, E22 = _e11, E12 = -_e12, E21 = E12;
				double k11 = g11 * E11 * this->get_gij2(0, 0) + g11 * E12 * this->get_gij2(1, 0) + g12 * E21 * this->get_gij2(0, 0) + g12 * E22 * this->get_gij2(1, 0);
				double k12 = g11 * E11 * this->get_gij2(0, 1) + g11 * E12 * this->get_gij2(1, 1) + g12 * E21 * this->get_gij2(0, 1) + g12 * E22 * this->get_gij2(1, 1);
				double k22 = g21 * E11 * this->get_gij2(0, 1) + g21 * E12 * this->get_gij2(1, 1) + g22 * E21 * this->get_gij2(0, 1) + g22 * E22 * this->get_gij2(1, 1);
				k11 += this->get_gij2(0, 0) * E11 * g11 + this->get_gij2(0, 0) * E12 * g21 + this->get_gij2(0, 1) * E21 * g11 + this->get_gij2(0, 1) * E22 * g21;
				k12 += this->get_gij2(0, 0) * E11 * g12 + this->get_gij2(0, 0) * E12 * g22 + this->get_gij2(0, 1) * E21 * g12 + this->get_gij2(0, 1) * E22 * g22;
				k22 += this->get_gij2(1, 0) * E11 * g12 + this->get_gij2(1, 0) * E12 * g22 + this->get_gij2(1, 1) * E21 * g12 + this->get_gij2(1, 1) * E22 * g22;

				double k21 = k12;
				*/
				s1 = -g21 * v1 - g22 * v2;
				s2 = g11 * v1 + g12 * v2;
				//double det = this->get_gij2(0, 0) * this->get_gij2(1, 1) - this->get_gij2(0, 1) * this->get_gij2(1, 0);
				//double sdet = sqrt(det);
				double dsdet = 0.5*(g11 * this->get_Gij2(0, 0) + 2*g12 * this->get_Gij2(0, 1)+ g22 * this->get_Gij2(1, 1))*sdet;
				//s1 /= sdet;
				//s2 /= sdet;
				s1 = s1 / sdet - _s1 / sdet / sdet * dsdet;
				s1 = s1 / sdet - _s1 / sdet / sdet * dsdet;

				val = -w2 * scale * ((_e11 * _s1 * s1) + (_e12 * _s1 * s2) + (_e21 * _s2 * s1) + (_e22 * _s2 * s2));
				val += -w2 * scale * ((_e11 * s1 * _s1) + (_e12 * s1 * _s2) + (_e21 * s2 * _s1) + (_e22 * s2 * _s2));
				*ptr1 = val;
				ptr1++;
			}
		}
		double setEij_B(double v1, double v2, double w1, double w2, bool accurate)
		{
			
			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			/*double _E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double _E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double _E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double _E22 = w1 * (v2 * v2) + w2 * (s2 * s2);*/
			double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					
				}
			}
			double tr2 = 1;
			double tr = e11 + e22;
			double scale = tr2 / tr;
			e21 = e12;
			val = scale*(((e11 * v1 * s1) + (e12 * v1 * s2) + (e21 * v2 * s1) + (e22 * v2 * s2)));
			return val;
		}

		
		void setEij_B_xi(double *ptr,double v1, double v2, double w1, double w2, bool accurate)
		{
			
			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double tr2 = 1;
			double tr = _e11 + _e22;
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				val = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					e12 +=(_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				e21 = e12;
				val = scale * ((e11 * v1 * s1) + (e12 * v1 * s2) + (e21 * v2 * s1) + (e22 * v2 * s2));
				val += -tr2/tr/tr * ((_e11 * v1 * s1) + (_e12 * v1 * s2) + (_e21 * v2 * s1) + (_e22 * v2 * s2))*(e11+e22);
				*ptr1=val;
				ptr1++;
			}
						
		}

	
		void setEij_B_eta(double* ptr, double v1, double v2, double w1, double w2, bool accurate)
		{
			
			double _E11 = 0, _E12 = 0, _E22 = 0, _E21;

			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			double val;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij2(0, 0) + v2 * v1 * this->get_gij2(1, 0) + v1 * v2 * this->get_gij2(0, 1) + v2 * v2 * this->get_gij2(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij2(0, 0) * v1 + this->get_gij2(0, 1) * v2;
				V2 = this->get_gij2(1, 0) * v1 + this->get_gij2(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij2(0, 0) + s2 * s1 * this->get_gij2(1, 0) + s1 * s2 * this->get_gij2(0, 1) + s2 * s2 * this->get_gij2(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij2(0, 0) * s1 + this->get_gij2(0, 1) * s2;
				S2 = this->get_gij2(1, 0) * s1 + this->get_gij2(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}


			double _e11 = 0, _e22 = 0, _e12 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];
					_e22 += _ref->buf_xi[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];

					_e11 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					_e22 += _ref->buf_eta[s] * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];

				}
			}
			double _e21 = _e12;
			double tr2 = 1;
			double tr = _e11 + _e22;
			double scale = tr2 / tr;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double e11 = 0, e12 = 0, e21 = 0, e22 = 0;
				val = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					e12 +=  (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];
					e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				e21 = e12;
				val = scale * ((e11 * v1 * s1) + (e12 * v1 * s2) + (e21 * v2 * s1) + (e22 * v2 * s2));
				val += -tr2 / tr / tr * ((_e11 * v1 * s1) + (_e12 * v1 * s2) + (_e21 * v2 * s1) + (_e22 * v2 * s2)) * (e11 + e22);
				*ptr1 = val;
				ptr1++;
			}

		}

		
		double theta(bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double s21 = s12;
			s11 *= sc;s12 *= sc;s21 *= sc;s22 *= sc;

			if (accurate)
			{
				return this->get_gij2(0, 0) * (s11 * S12 + s12 * S22)
					+ this->get_gij2(0, 1) * (s21 * S12 + s22 * S22)
					- this->get_gij2(1, 0) * (s11 * S11 + s12 * S21)
					- this->get_gij2(1, 1) * (s21 * S11 + s22 * S21);
			}
			else {
				return _ref->get__gij(0, 0) * (s11 * S12 + s12 * S22)
					+ _ref->get__gij(0, 1) * (s21 * S12 + s22 * S22)
					- _ref->get__gij(1, 0) * (s11 * S11 + s12 * S21)
					- _ref->get__gij(1, 1) * (s21 * S11 + s22 * S21);
			}


		}
		void theta_phi(double* ptr, bool accurate)
		{
			//double S11 = 0;
			//double S12 = 0;
			//double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				//S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				//S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double s21 = s12;
			s11 *= sc;
			s12 *= sc;
			s21 *= sc;
			s22 *= sc;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double S21 = S12;
				double S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

				if (accurate)
				{
					//double s21 = s12;

					/*double a11 = s11 * this->get_Gij2(0, 0) * S11 + s11 * this->get_Gij2(0, 1) * S21 + s12 * this->get_Gij2(1, 0) * S11 + s12 * this->get_Gij2(1, 1) * S21;
					double a21 = s21 * this->get_Gij2(0, 0) * S11 + s21 * this->get_Gij2(0, 1) * S21 + s22 * this->get_Gij2(1, 0) * S11 + s22 * this->get_Gij2(1, 1) * S21;
					double a12 = s11 * this->get_Gij2(0, 0) * S12 + s11 * this->get_Gij2(0, 1) * S22 + s12 * this->get_Gij2(1, 0) * S12 + s12 * this->get_Gij2(1, 1) * S22;
					double a22 = s21 * this->get_Gij2(0, 0) * S12 + s21 * this->get_Gij2(0, 1) * S22 + s22 * this->get_Gij2(1, 0) * S12 + s22 * this->get_Gij2(1, 1) * S22;

					double c11 = _ref->get__gij(0, 0) * this->get_Gij2(0, 0) + _ref->get__gij(0, 1) * this->get_Gij2(1, 0);
					double c21 = _ref->get__gij(1, 0) * this->get_Gij2(0, 0) + _ref->get__gij(1, 1) * this->get_Gij2(1, 0);
					double c12 = _ref->get__gij(0, 0) * this->get_Gij2(0, 1) + _ref->get__gij(0, 1) * this->get_Gij2(1, 1);
					double c22 = _ref->get__gij(1, 0) * this->get_Gij2(0, 1) + _ref->get__gij(1, 1) * this->get_Gij2(1, 1);

					double d11 = c11 * a11 + c12 * a21;
					double d21 = c21 * a11 + c22 * a21;
					double d12 = c11 * a12 + c12 * a22;
					double d22 = c21 * a12 + c22 * a22;
					val= a12 - a21;*/
					val += this->get_gij2(0, 0) * (s11 * S12 + s12 * S22);
					val += this->get_gij2(0, 1) * (s21 * S12 + s22 * S22);
					val -= this->get_gij2(1, 0) * (s11 * S11 + s12 * S21);
					val -= this->get_gij2(1, 1) * (s21 * S11 + s22 * S21);
					//val += s11 * this->get_Gij2(0, 0) * (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
					//val += s11 * this->get_Gij2(0, 1) * (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
					//val += s12 * this->get_Gij2(1, 0) * (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
					//val += s12 * this->get_Gij2(1, 1) * (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

					//val -= (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * this->get_Gij2(0, 0) * s12;
					//val -= (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * this->get_Gij2(0, 1) * s22;
					//val -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * this->get_Gij2(1, 0) * s12;
					//val -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * this->get_Gij2(1, 1) * s22;
				}
				else {
					val += _ref->get__gij(0, 0) * (s11 * S12 + s12 * S22);
					val += _ref->get__gij(0, 1) * (s21 * S12 + s22 * S22);
					val -= _ref->get__gij(1, 0) * (s11 * S11 + s12 * S21);
					val -= _ref->get__gij(1, 1) * (s21 * S11 + s22 * S21);
					//val += s11 * _ref->get__Gij(0, 0) * (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
					//val += s11 * _ref->get__Gij(0, 1) * (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
					//val += s12 * _ref->get__Gij(1, 0) * (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
					//val += s12 * _ref->get__Gij(1, 1) * (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

					//val -= (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->get__Gij(0, 0) * s12;
					//val -= (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->get__Gij(0, 1) * s22;
					//val -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->get__Gij(1, 0) * s12;
					//val -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->get__Gij(1, 1) * s22;
				}
				*ptr1 = val;
				ptr1++;
			}
		}
		void theta_Z(double* ptr, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			//double s11 = 0;
			//double s12 = 0;
			//double s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double s22 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double s12 = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double s21 = s12;
				double s11 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				s11 *= sc;
				s12 *= sc;
				s21 *= sc;
				s22 *= sc;

				if (accurate)
				{
					//double S21 = S12;


					//double a11 = s11 * this->get_Gij2(0, 0) * S11 + s11 * this->get_Gij2(0, 1) * S21 + s12 * this->get_Gij2(1, 0) * S11 + s12 * this->get_Gij2(1, 1) * S21;
					//double a21 = s21 * this->get_Gij2(0, 0) * S11 + s21 * this->get_Gij2(0, 1) * S21 + s22 * this->get_Gij2(1, 0) * S11 + s22 * this->get_Gij2(1, 1) * S21;
					//double a12 = s11 * this->get_Gij2(0, 0) * S12 + s11 * this->get_Gij2(0, 1) * S22 + s12 * this->get_Gij2(1, 0) * S12 + s12 * this->get_Gij2(1, 1) * S22;
					//double a22 = s21 * this->get_Gij2(0, 0) * S12 + s21 * this->get_Gij2(0, 1) * S22 + s22 * this->get_Gij2(1, 0) * S12 + s22 * this->get_Gij2(1, 1) * S22;

					/*double c11 = _ref->get__gij(0, 0) * this->get_Gij2(0, 0) + _ref->get__gij(0, 1) * this->get_Gij2(1, 0);
					double c21 = _ref->get__gij(1, 0) * this->get_Gij2(0, 0) + _ref->get__gij(1, 1) * this->get_Gij2(1, 0);
					double c12 = _ref->get__gij(0, 0) * this->get_Gij2(0, 1) + _ref->get__gij(0, 1) * this->get_Gij2(1, 1);
					double c22 = _ref->get__gij(1, 0) * this->get_Gij2(0, 1) + _ref->get__gij(1, 1) * this->get_Gij2(1, 1);

					double d11 = c11 * a11 + c12 * a21;
					double d21 = c21 * a11 + c22 * a21;
					double d12 = c11 * a12 + c12 * a22;
					double d22 = c21 * a12 + c22 * a22;*/
					val += this->get_gij2(0, 0) * (s11 * S12 + s12 * S22);
					val += this->get_gij2(0, 1) * (s21 * S12 + s22 * S22);
					val -= this->get_gij2(1, 0) * (s11 * S11 + s12 * S21);
					val -= this->get_gij2(1, 1) * (s21 * S11 + s22 * S21);



				}
				else {
					val += _ref->get__gij(0, 0) * (s11 * S12 + s12 * S22);
					val += _ref->get__gij(0, 1) * (s21 * S12 + s22 * S22);
					val -= _ref->get__gij(1, 0) * (s11 * S11 + s12 * S21);
					val -= _ref->get__gij(1, 1) * (s21 * S11 + s22 * S21);
					//val -= S11 * _ref->get__Gij(0, 0) * (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
					//val -= S11 * _ref->get__Gij(0, 1) * (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
					//val -= S12 * _ref->get__Gij(1, 0) * (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
					//val -= S12 * _ref->get__Gij(1, 1) * (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);

					//val += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->get__Gij(0, 0) * S12;
					//val += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->get__Gij(0, 1) * S22;
					//val += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->get__Gij(1, 0) * S12;
					//val += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->get__Gij(1, 1) * S22;
				}
				*ptr1 = val;
				ptr1++;
			}
		}
		double CC(double v1, double v2, bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 += -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double s11 = 0, s12 = 0, s22 = 0;
			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + 2 * this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s12 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + 2 * this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + 2 * _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s12 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + 2 * _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
			}
			s11 *= sc;
			s12 *= sc;
			s22 *= sc;

			return (s12 * v1 + s22 * v2) * (-S11 * v1 - S12 * v2) - (-s11 * v1 - s12 * v2) * (S12 * v1 + S22 * v2);
		}
		void CC_phi(double* ptr, double v1, double v2, bool accurate)
		{
			//double S11 = 0;
			//double S12 = 0;
			//double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 += -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				//S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				//S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double s11 = 0, s12 = 0, s22 = 0;
			if (accurate)
			{
				s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + 2 * this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
				s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s12 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
				s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + 2 * this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
			}
			else
			{
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + 2 * _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s12 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + 2 * _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
			}
			s11 *= sc;
			s12 *= sc;
			s22 *= sc;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double S21 = S12;
				val += (s12 * v1 + s22 * v2) * (-S11 * v1 - S12 * v2);
				val -= (-s11 * v1 - s12 * v2) * (S12 * v1 + S22 * v2);
				*ptr1 = val;
				ptr1++;
			}
		}
		void CC_Z(double* ptr, double v1, double v2, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			//double s11 = 0;
			//double s12 = 0;
			//double s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double _s22 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _s12 = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _s11 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double s11 = 0, s12 = 0, s22 = 0;
				if (accurate)
				{
					s11 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 0) + 2 * this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 0) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 0);
					s12 = this->get_gij2(0, 0) * _s11 * this->get_gij2(0, 1) + this->get_gij2(0, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(0, 1) * _s12 * this->get_gij2(0, 1) + this->get_gij2(0, 1) * _s22 * this->get_gij2(1, 1);
					s22 = this->get_gij2(1, 0) * _s11 * this->get_gij2(0, 1) + 2 * this->get_gij2(1, 0) * _s12 * this->get_gij2(1, 1) + this->get_gij2(1, 1) * _s22 * this->get_gij2(1, 1);
				}
				else
				{
					s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + 2 * _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
					s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s12 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
					s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + 2 * _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				}
				s11 *= sc;
				s12 *= sc;
				s22 *= sc;

				val += (s12 * v1 + s22 * v2) * (-S11 * v1 - S12 * v2);
				val -= (-s11 * v1 - s12 * v2) * (S12 * v1 + S22 * v2);

				*ptr1 = val;
				ptr1++;
			}
		}
		double BCEQ(double _t1/*up*/, double _t2/*up*/,double stt, double a,double b,bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double length = 0;
			double t1 = 0;//up
			double t2 = 0;//up
			double T1 = 0;//down
			double T2 = 0;//down
			double n1 = 0;//up
			double n2 = 0;//up
			double N1 = 0;//down
			double N2 = 0;//down
			double gamma = 0;
			double eta = 0;
			if (accurate)
			{
				length = sqrt(_t1 * _t1 * this->get_gij(0, 0) + _t2 * _t1 * this->get_gij(1, 0) + _t1 * _t2 * this->get_gij(0, 1) + _t2 * _t2 * this->get_gij(1, 1));
				gamma = length* length;
				t1 = _t1/length;
				t2 = _t2/length;
				T1 = this->get_gij(0, 0) * t1 + this->get_gij(0, 1) * t2;
				T2 = this->get_gij(1, 0) * t1 + this->get_gij(1, 1) * t2;
				n1 = T2;
				n2 = -T1;
				length = sqrt(n1 * n1 * this->get_gij(0, 0) + n2 * n1 * this->get_gij(1, 0) + n1 * n2 * this->get_gij(0, 1) + n2 * n2 * this->get_gij(1, 1));
				n1 /= length;
				n2 /= length;
				eta = length * length;
				N1 = this->get_gij(0, 0) * n1 + this->get_gij(0, 1) * n2;
				N2 = this->get_gij(1, 0) * n1 + this->get_gij(1, 1) * n2;
			}
			else {
				length = sqrt(_t1 * _t1 * _ref->get__gij(0, 0) + _t2 * _t1 * _ref->get__gij(1, 0) + _t1 * _t2 * _ref->get__gij(0, 1) + _t2 * _t2 * _ref->get__gij(1, 1));
				gamma = length * length;
				t1 = _t1 / length;
				t2 = _t2 / length;
				T1 = _ref->get__gij(0, 0) * t1 + _ref->get__gij(0, 1) * t2;
				T2 = _ref->get__gij(1, 0) * t1 + _ref->get__gij(1, 1) * t2;
				n1 = T2;
				n2 = -T1;
				length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + n2 * n1 * _ref->get__gij(1, 0) + n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
				eta = length * length;
				n1 /= length;
				n2 /= length;
				N1 = _ref->get__gij(0, 0) * n1 + _ref->get__gij(0, 1) * n2;
				N2 = _ref->get__gij(1, 0) * n1 + _ref->get__gij(1, 1) * n2;
			}
			double d1 = 0, d2 = 0, D1 = 0, D2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				d1 += _ref->d1[0][s] * _ref->buf_phi[s];
				d2 += _ref->d1[1][s] * _ref->buf_phi[s];
				D1 += _ref->d1[0][s] * _ref->buf_z[s];
				D2 += _ref->d1[1][s] * _ref->buf_z[s];
				s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double s21 = s12;
			//s11 = 0;
			// 
	
			//if (stt < -100)
			{
				val += (s11 * t1 * t1 + 2 * s12 * t1 * t2 + s22 * t2 * t2) * (D1 * n1 + D2 * n2) ; //no unit
			}
			/*else {
				//val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2);
				val += stt * (D1 * t1 + D2 * t2) * sc;
			}*/
			double f1 = a * this->get_gi(0, 0) + b * this->get_gi(0, 1);
			double f2 = a * this->get_gi(1, 0) + b * this->get_gi(1, 1);
			double dc = f1 * n1 + f2 * n2;

			if (stt < -100)
			{

				val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2 - dc);
			}
			else {
				val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * stt;

			}
			return val;
		}
		void BCEQ_z(double* ptr, double _t1, double _t2,double stt,double a,double b ,bool accurate)
		{


			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double length = 0;
			double t1 = 0;//up
			double t2 = 0;//up
			double T1 = 0;//down
			double T2 = 0;//down
			double n1 = 0;//up
			double n2 = 0;//up
			double N1 = 0;//down
			double N2 = 0;//down
			double gamma = 0;
			double eta = 0;
			if (accurate)
			{
				length = sqrt(_t1 * _t1 * this->get_gij(0, 0) + _t2 * _t1 * this->get_gij(1, 0) + _t1 * _t2 * this->get_gij(0, 1) + _t2 * _t2 * this->get_gij(1, 1));
				gamma = length * length;
				t1 = _t1 / length;
				t2 = _t2 / length;
				T1 = this->get_gij(0, 0) * t1 + this->get_gij(0, 1) * t2;
				T2 = this->get_gij(1, 0) * t1 + this->get_gij(1, 1) * t2;
				n1 = T2;
				n2 = -T1;
				length = sqrt(n1 * n1 * this->get_gij(0, 0) + n2 * n1 * this->get_gij(1, 0) + n1 * n2 * this->get_gij(0, 1) + n2 * n2 * this->get_gij(1, 1));

				n1 /= length;
				n2 /= length;
				N1 = this->get_gij(0, 0) * n1 + this->get_gij(0, 1) * n2;
				N2 = this->get_gij(1, 0) * n1 + this->get_gij(1, 1) * n2;
			}
			else {
				length = sqrt(_t1 * _t1 * _ref->get__gij(0, 0) + _t2 * _t1 * _ref->get__gij(1, 0) + _t1 * _t2 * _ref->get__gij(0, 1) + _t2 * _t2 * _ref->get__gij(1, 1));
				gamma = length * length;
				t1 = _t1 / length;
				t2 = _t2 / length;
				T1 = _ref->get__gij(0, 0) * t1 + _ref->get__gij(0, 1) * t2;
				T2 = _ref->get__gij(1, 0) * t1 + _ref->get__gij(1, 1) * t2;
				n1 = T2;
				n2 = -T1;
				length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + n2 * n1 * _ref->get__gij(1, 0) + n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
				n1 /= length;
				n2 /= length;
				N1 = _ref->get__gij(0, 0) * n1 + _ref->get__gij(0, 1) * n2;
				N2 = _ref->get__gij(1, 0) * n1 + _ref->get__gij(1, 1) * n2;
			}
			double d1 = 0, d2 = 0, D1 = 0, D2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				D1 += _ref->d1[0][s] * _ref->buf_z[s];
				D2 += _ref->d1[1][s] * _ref->buf_z[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double* ptr1 = ptr;
			double f1 = a * this->get_gi(0, 0) + b * this->get_gi(0, 1);
			double f2 = a * this->get_gi(1, 0) + b * this->get_gi(1, 1);
			double dc = f1 * n1 + f2 * n2;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				d1 = _ref->d1[0][s];
				d2 = _ref->d1[1][s];
				s11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double S21 = S12;
				double s21 = s12;
				//s11 = 0;
				{
					val += (s11 * t1 * t1 + 2 * s12 * t1 * t2 + s22 * t2 * t2) * (D1 * n1 + D2 * n2);
				}
				if (stt < -100)
				{
					val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2-dc);

				}
				else {
					//val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (stt);
				}

				

				*ptr1 = val;
				ptr1 ++;
			}


		}
		void BCEQ_phi(double* ptr, double _t1, double _t2, double stt, double a, double b, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double length = 0;
			double t1 = 0;//up
			double t2 = 0;//up
			double T1 = 0;//down
			double T2 = 0;//down
			double n1 = 0;//up
			double n2 = 0;//up
			double N1 = 0;//down
			double N2 = 0;//down
			double gamma = 0;
			//double eta = 0;
			//double det = 0;
			if (accurate)
			{
				length = sqrt(_t1 * _t1 * this->get_gij(0, 0) + _t2 * _t1 * this->get_gij(1, 0) + _t1 * _t2 * this->get_gij(0, 1) + _t2 * _t2 * this->get_gij(1, 1));
				gamma = length * length;
				t1 = _t1 / length;
				t2 = _t2 / length;
				T1 = this->get_gij(0, 0) * t1 + this->get_gij(0, 1) * t2;
				T2 = this->get_gij(1, 0) * t1 + this->get_gij(1, 1) * t2;
				n1 = T2;
				n2 = -T1;
				length = sqrt(n1 * n1 * this->get_gij(0, 0) + n2 * n1 * this->get_gij(1, 0) + n1 * n2 * this->get_gij(0, 1) + n2 * n2 * this->get_gij(1, 1));
				//eta = length * length;
				n1 /= length;
				n2 /= length;
				N1 = this->get_gij(0, 0) * n1 + this->get_gij(0, 1) * n2;
				N2 = this->get_gij(1, 0) * n1 + this->get_gij(1, 1) * n2;
			}
			else {
				length = sqrt(_t1 * _t1 * _ref->get__gij(0, 0) + _t2 * _t1 * _ref->get__gij(1, 0) + _t1 * _t2 * _ref->get__gij(0, 1) + _t2 * _t2 * _ref->get__gij(1, 1));
				gamma = length * length;
				t1 = _t1 / length;
				t2 = _t2 / length;
				T1 = _ref->get__gij(0, 0) * t1 + _ref->get__gij(0, 1) * t2;
				T2 = _ref->get__gij(1, 0) * t1 + _ref->get__gij(1, 1) * t2;
				n1 = T2;
				n2 = -T1;
				length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + n2 * n1 * _ref->get__gij(1, 0) + n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
				//eta = length * length;
				n1 /= length;
				n2 /= length;
				N1 = _ref->get__gij(0, 0) * n1 + _ref->get__gij(0, 1) * n2;
				N2 = _ref->get__gij(1, 0) * n1 + _ref->get__gij(1, 1) * n2;
			}
			//det = eta * gamma;
			double d1 = 0, d2 = 0, D1 = 0, D2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				d1 += _ref->d1[0][s] * _ref->buf_phi[s];
				d2 += _ref->d1[1][s] * _ref->buf_phi[s];
				s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
			}
			double* ptr1 = ptr;
			double f1 = a * this->get_gi(0, 0) + b * this->get_gi(0, 1);//down
			double f2 = a * this->get_gi(1, 0) + b * this->get_gi(1, 1);//down
			double dc = f1 * n1 + f2 * n2;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				D1 = _ref->d1[0][s];
				D2 = _ref->d1[1][s];
				S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double S21 = S12;
				double s21 = s12;
				//s11 = 0;
				//if (stt < -100)
				{
					val += (s11 * t1 * t1 + 2 * s12 * t1 * t2 + s22 * t2 * t2) * (D1 * n1 + D2 * n2); //no unit
				}
				/*
				else {
					val += stt * (D1 * t1 + D2 * t2) * sc;
				}*/
				if (stt < -100)
				{
					val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2-dc);
				}
				else {
					val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (stt);
				}

				
				*ptr1 = val;
				ptr1++;
			}
		}
		




		double mix_BC(double v1, double v2, bool accurate)
		{
			double val = 0;

			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}

			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double _s21 = _s12;



			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double s21 = 0;
			if (accurate)
			{
				s11 = this->get_gij(0, 0) * _s11 * this->get_gij(0, 0) + this->get_gij(0, 0) * _s12 * this->get_gij(1, 0) + this->get_gij(0, 1) * _s21 * this->get_gij(0, 0) + this->get_gij(0, 1) * _s22 * this->get_gij(1, 0);
				s12 = this->get_gij(0, 0) * _s11 * this->get_gij(0, 1) + this->get_gij(0, 0) * _s12 * this->get_gij(1, 1) + this->get_gij(0, 1) * _s21 * this->get_gij(0, 1) + this->get_gij(0, 1) * _s22 * this->get_gij(1, 1);
				s22 = this->get_gij(1, 0) * _s11 * this->get_gij(0, 1) + this->get_gij(1, 0) * _s12 * this->get_gij(1, 1) + this->get_gij(1, 1) * _s21 * this->get_gij(0, 1) + this->get_gij(1, 1) * _s22 * this->get_gij(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc;
			s12 *= sc;
			s22 *= sc;
			s21 *= sc;
			//val +=  ((s21 * V1 + s22 * V2) * (-S11 * v1 - S12 * v2) - (-s11 * V1 - s12 * V2) * (S21 * v1 + S22 * v2));
			double E11 = 1 * (v1 * v1);// +w2 * (S1 * s1);
			double E12 = 1 * (v1 * v2);// + w2 * (S1 * s2);
			double E21 = 1 * (v2 * v1);// + w2 * (S2 * s1);
			double E22 = 1 * (v2 * v2);// + w2 * (S2 * s2);
			val += (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22);
			val -= (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21);
			return val;
		}
		void mix_BC_phi(double* ptr, double v1, double v2, bool accurate)
		{
			//double S11 = 0;
			//double S12 = 0;
			//double S22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_s22 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s12 -= (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				_s11 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				//S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				//S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double _s21 = _s12;
			double* ptr1 = ptr;
			double s11 = 0;
			double s12 = 0;
			double s21 = 0;
			double s22 = 0;

			if (accurate)
			{
				s11 = this->get_gij(0, 0) * _s11 * this->get_gij(0, 0) + this->get_gij(0, 0) * _s12 * this->get_gij(1, 0) + this->get_gij(0, 1) * _s21 * this->get_gij(0, 0) + this->get_gij(0, 1) * _s22 * this->get_gij(1, 0);
				s12 = this->get_gij(0, 0) * _s11 * this->get_gij(0, 1) + this->get_gij(0, 0) * _s12 * this->get_gij(1, 1) + this->get_gij(0, 1) * _s21 * this->get_gij(0, 1) + this->get_gij(0, 1) * _s22 * this->get_gij(1, 1);
				s22 = this->get_gij(1, 0) * _s11 * this->get_gij(0, 1) + this->get_gij(1, 0) * _s12 * this->get_gij(1, 1) + this->get_gij(1, 1) * _s21 * this->get_gij(0, 1) + this->get_gij(1, 1) * _s22 * this->get_gij(1, 1);
				s21 = s12;
			}
			else {
				s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
				s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
				s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
				s21 = s12;
			}
			s11 *= sc;
			s12 *= sc;
			s22 *= sc;
			s21 *= sc;

			double E11 = 1 * (v1 * v1);// +w2 * (S1 * s1);
			double E12 = 1 * (v1 * v2);// + w2 * (S1 * s2);
			double E21 = 1 * (v2 * v1);// + w2 * (S2 * s1);
			double E22 = 1 * (v2 * v2);// + w2 * (S2 * s2);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double S21 = S12;
				double S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);


				//val += w2 * ((s12 * v1 + s22 * v2) * (-(_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * v1 - (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * v2));
				//val -= w2 * ((-s11 * v1 - s12 * v2) * ((_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * v1 + (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * v2));
				//val += ((s21 * V1 + s22 * V2) * (-S11 * v1 - S12 * v2) - (-s11 * V1 - s12 * V2) * (S21 * v1 + S22 * v2));
				val += (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22);
				val -= (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21);


				*ptr1 = val;
				ptr1++;
			}
		}
		void mix_BC_Z(double* ptr, double v1, double v2,bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			//double s11 = 0;
			//double s12 = 0;
			//double s22 = 0;
			double s1 = 0, s2 = 0;//up
			double S1 = 0, S2 = 0;//down
			double V1 = 0, V2 = 0;//down
			double length = 0;
			if (accurate)
			{
				length = sqrt(v1 * v1 * this->get_gij(0, 0) + v2 * v1 * this->get_gij(1, 0) + v1 * v2 * this->get_gij(0, 1) + v2 * v2 * this->get_gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = this->get_gij(0, 0) * v1 + this->get_gij(0, 1) * v2;
				V2 = this->get_gij(1, 0) * v1 + this->get_gij(1, 1) * v2;
				s1 = V2;
				s2 = -V1;
				length = sqrt(s1 * s1 * this->get_gij(0, 0) + s2 * s1 * this->get_gij(1, 0) + s1 * s2 * this->get_gij(0, 1) + s2 * s2 * this->get_gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = this->get_gij(0, 0) * s1 + this->get_gij(0, 1) * s2;
				S2 = this->get_gij(1, 0) * s1 + this->get_gij(1, 1) * s2;
			}
			else {
				length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + v2 * v1 * _ref->get__gij(1, 0) + v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
				v1 /= length;
				v2 /= length;
				V1 = _ref->get__gij(0, 0) * v1 + _ref->get__gij(0, 1) * v2;
				V2 = _ref->get__gij(1, 0) * v1 + _ref->get__gij(1, 1) * v2;
				S1 = v2;
				S2 = -v1;
				length = sqrt(s1 * s1 * _ref->get__gij(0, 0) + s2 * s1 * _ref->get__gij(1, 0) + s1 * s2 * _ref->get__gij(0, 1) + s2 * s2 * _ref->get__gij(1, 1));
				s1 /= length;
				s2 /= length;
				S1 = _ref->get__gij(0, 0) * s1 + _ref->get__gij(0, 1) * s2;
				S2 = _ref->get__gij(1, 0) * s1 + _ref->get__gij(1, 1) * s2;
			}

			for (int s = 0; s < _ref->_nNode; s++)
			{
				//s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				//s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];
			}
			double S21 = S12;
			double* ptr1 = ptr;
			double E11 = 1 * (v1 * v1);// +w2 * (S1 * s1);
			double E12 = 1 * (v1 * v2);// + w2 * (S1 * s2);
			double E21 = 1 * (v2 * v1);// + w2 * (S2 * s1);
			double E22 = 1 * (v2 * v2);// + w2 * (S2 * s2);

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;
				double _s22 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _s12 = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _s21 = _s12;
				double _s11 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double s11 = 0, s12 = 0, s22 = 0, s21 = 0;
				if (accurate)
				{
					s11 = this->get_gij(0, 0) * _s11 * this->get_gij(0, 0) + this->get_gij(0, 0) * _s12 * this->get_gij(1, 0) + this->get_gij(0, 1) * _s21 * this->get_gij(0, 0) + this->get_gij(0, 1) * _s22 * this->get_gij(1, 0);
					s12 = this->get_gij(0, 0) * _s11 * this->get_gij(0, 1) + this->get_gij(0, 0) * _s12 * this->get_gij(1, 1) + this->get_gij(0, 1) * _s21 * this->get_gij(0, 1) + this->get_gij(0, 1) * _s22 * this->get_gij(1, 1);
					s22 = this->get_gij(1, 0) * _s11 * this->get_gij(0, 1) + this->get_gij(1, 0) * _s12 * this->get_gij(1, 1) + this->get_gij(1, 1) * _s21 * this->get_gij(0, 1) + this->get_gij(1, 1) * _s22 * this->get_gij(1, 1);
					s21 = s12;
				}
				else {
					s11 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 0) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 0) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 0) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 0);
					s12 = _ref->get__gij(0, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(0, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(0, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(0, 1) * _s22 * _ref->get__gij(1, 1);
					s22 = _ref->get__gij(1, 0) * _s11 * _ref->get__gij(0, 1) + _ref->get__gij(1, 0) * _s12 * _ref->get__gij(1, 1) + _ref->get__gij(1, 1) * _s21 * _ref->get__gij(0, 1) + _ref->get__gij(1, 1) * _s22 * _ref->get__gij(1, 1);
					s21 = s12;
				}
				s11 *= sc;
				s12 *= sc;
				s22 *= sc;
				s21 *= sc;
				//val += ((s21 * V1 + s22 * V2) * (-S11 * v1 - S12 * v2) - (-s11 * V1 - s12 * V2) * (S21 * v1 + S22 * v2));
				

				val += (s11 * E11 * S12 + s11 * E12 * S22 + s12 * E21 * S12 + s12 * E22 * S22);
				val -= (s21 * E11 * S11 + s21 * E12 * S21 + s22 * E21 * S11 + s22 * E22 * S21);


				*ptr1 = val;
				ptr1++;
			}
		}
		double D1()
		{
			double S11 = 0, S12 = 0, S22 = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				S11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_phi[i];
				s11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
				s12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
				s22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
			}
			double b11 = S11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[3] + S22 * _ref->_Gi[3] * _ref->_Gi[3];
			double b12 = S11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[4] + S22 * _ref->_Gi[3] * _ref->_Gi[4];
			double b22 = S11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[1] * _ref->_Gi[4] + S22 * _ref->_Gi[4] * _ref->_Gi[4];
			double c11 = s11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[3] + s22 * _ref->_Gi[3] * _ref->_Gi[3];
			double c12 = s11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[4] + s22 * _ref->_Gi[3] * _ref->_Gi[4];
			double c22 = s11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[1] * _ref->_Gi[4] + s22 * _ref->_Gi[4] * _ref->_Gi[4];

			return b12 * c22 - b22 * c12;
		}
		double D2()
		{
			double S11 = 0, S12 = 0, S22 = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				S11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_phi[i];
				s11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
				s12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
				s22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
			}
			double b11 = S11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[3] + S22 * _ref->_Gi[3] * _ref->_Gi[3];
			double b12 = S11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[4] + S22 * _ref->_Gi[3] * _ref->_Gi[4];
			double b22 = S11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[1] * _ref->_Gi[4] + S22 * _ref->_Gi[4] * _ref->_Gi[4];
			double c11 = s11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[3] + s22 * _ref->_Gi[3] * _ref->_Gi[3];
			double c12 = s11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[4] + s22 * _ref->_Gi[3] * _ref->_Gi[4];
			double c22 = s11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[1] * _ref->_Gi[4] + s22 * _ref->_Gi[4] * _ref->_Gi[4];

			return b12 * c11 - b11 * c12;
		}
		void D1_phi(double* ptr)
		{
			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				S11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_phi[i];
				s11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
				s12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
				s22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
			}
			double b11 = S11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[3] + S22 * _ref->_Gi[3] * _ref->_Gi[3];
			double b12 = S11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[4] + S22 * _ref->_Gi[3] * _ref->_Gi[4];
			double b22 = S11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[1] * _ref->_Gi[4] + S22 * _ref->_Gi[4] * _ref->_Gi[4];
			double c11 = s11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[3] + s22 * _ref->_Gi[3] * _ref->_Gi[3];
			double c12 = s11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[4] + s22 * _ref->_Gi[3] * _ref->_Gi[4];
			double c22 = s11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[1] * _ref->_Gi[4] + s22 * _ref->_Gi[4] * _ref->_Gi[4];
			//b12* c22 - b22 * c12;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				double val = 0;
				val = 0;
				val += b12 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[1] * _ref->_Gi[1];
				val += b12 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[1] * _ref->_Gi[4];
				val += b12 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[4] * _ref->_Gi[4];

				val += -b22 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[1];
				val += -b22 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[4];
				val += -b22 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[3] * _ref->_Gi[4];
				*ptr1 = val;
				ptr1++;
			}


		}
		void D1_z(double* ptr)
		{
			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				S11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_phi[i];
				s11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
				s12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
				s22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
			}
			double b11 = S11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[3] + S22 * _ref->_Gi[3] * _ref->_Gi[3];
			double b12 = S11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[4] + S22 * _ref->_Gi[3] * _ref->_Gi[4];
			double b22 = S11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[1] * _ref->_Gi[4] + S22 * _ref->_Gi[4] * _ref->_Gi[4];
			double c11 = s11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[3] + s22 * _ref->_Gi[3] * _ref->_Gi[3];
			double c12 = s11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[4] + s22 * _ref->_Gi[3] * _ref->_Gi[4];
			double c22 = s11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[1] * _ref->_Gi[4] + s22 * _ref->_Gi[4] * _ref->_Gi[4];
			//b12* c22 - b22 * c12;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				double val = 0;
				val = 0;
				val += c22 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[1];
				val += c22 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[4];
				val += c22 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[3] * _ref->_Gi[4];

				val += -c12 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[1] * _ref->_Gi[1];
				val += -c12 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[1] * _ref->_Gi[4];
				val += -c12 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[4] * _ref->_Gi[4];
				*ptr1 = val;
				ptr1++;
			}
		}


		void D2_phi(double* ptr)
		{
			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				S11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_phi[i];
				s11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
				s12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
				s22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
			}
			//b12 * c11 - b11 * c12;
			double b11 = S11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[3] + S22 * _ref->_Gi[3] * _ref->_Gi[3];
			double b12 = S11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[4] + S22 * _ref->_Gi[3] * _ref->_Gi[4];
			double b22 = S11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[1] * _ref->_Gi[4] + S22 * _ref->_Gi[4] * _ref->_Gi[4];
			double c11 = s11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[3] + s22 * _ref->_Gi[3] * _ref->_Gi[3];
			double c12 = s11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[4] + s22 * _ref->_Gi[3] * _ref->_Gi[4];
			double c22 = s11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[1] * _ref->_Gi[4] + s22 * _ref->_Gi[4] * _ref->_Gi[4];

			for (int i = 0; i < _ref->_nNode; i++)
			{
				double val = 0;
				val = 0;
				val += b12 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[0];
				val += b12 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[3];
				val += b12 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[3] * _ref->_Gi[3];

				val += -b11 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[1];
				val += -b11 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[4];
				val += -b11 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[3] * _ref->_Gi[4];
				*ptr1 = val;
				ptr1++;
			}


		}
		void D2_z(double* ptr)
		{
			double* ptr1 = ptr;
			double S11 = 0, S12 = 0, S22 = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				S11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_phi[i];
				S22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_phi[i];
				s11 += (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->buf_z[i];
				s12 += (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->buf_z[i];
				s22 += (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->buf_z[i];
			}
			double b11 = S11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[3] + S22 * _ref->_Gi[3] * _ref->_Gi[3];
			double b12 = S11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[0] * _ref->_Gi[4] + S22 * _ref->_Gi[3] * _ref->_Gi[4];
			double b22 = S11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * S12 * _ref->_Gi[1] * _ref->_Gi[4] + S22 * _ref->_Gi[4] * _ref->_Gi[4];
			double c11 = s11 * _ref->_Gi[0] * _ref->_Gi[0] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[3] + s22 * _ref->_Gi[3] * _ref->_Gi[3];
			double c12 = s11 * _ref->_Gi[0] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[0] * _ref->_Gi[4] + s22 * _ref->_Gi[3] * _ref->_Gi[4];
			double c22 = s11 * _ref->_Gi[1] * _ref->_Gi[1] + 2 * s12 * _ref->_Gi[1] * _ref->_Gi[4] + s22 * _ref->_Gi[4] * _ref->_Gi[4];
			//return b12 * c11 - b11 * c12;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				double val = 0;
				val = 0;
				val += c11 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[1];
				val += c11 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[4];
				val += c11 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[3] * _ref->_Gi[4];

				val += -c12 * (_ref->d2[0][i] - _ref->_Gammaijk[0] * _ref->d1[0][i] - _ref->_Gammaijk[1] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[0];
				val += -c12 * 2 * (_ref->d2[1][i] - _ref->_Gammaijk[2] * _ref->d1[0][i] - _ref->_Gammaijk[3] * _ref->d1[1][i]) * _ref->_Gi[0] * _ref->_Gi[3];
				val += -c12 * (_ref->d2[3][i] - _ref->_Gammaijk[6] * _ref->d1[0][i] - _ref->_Gammaijk[7] * _ref->d1[1][i]) * _ref->_Gi[3] * _ref->_Gi[3];
				*ptr1 = val;
				ptr1++;
			}
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
			return sc * _ref->__mat[I * _nNode + J];
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
		void U_phi_z(_mySparse* mat, long long* index, int N, long long* map, double R)
		{
			if (map == 0)
			{
				for (int i = 0; i < _nNode; i++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						int I = index[i];
						int J = index[j];
						double val = __grad_phi[i] * __grad_z[j] * R;
						mat->_mat[0].coeffRef(I, J) += val;
						//mat->_mat[0].data().valuePtr()[map[I*N+J]] += val;
					}
				}
			}
			else
			{
				for (int i = 0; i < _nNode; i++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						int I = index[i];
						int J = index[j];
						double val = __grad_phi[i] * __grad_z[j] * R;
						//mat->_mat[0].coeffRef(I, J) += val;
						mat->_mat[0].data().valuePtr()[map[I * N + J]] += val;
					}
				}
			}
		}
		void D_phi_z(_mySparse* mat, long long* index, int N, long long* map, double sc)
		{
			if (map == 0)
			{
				for (int i = 0; i < _nNode; i++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						int I = index[i];
						int J = index[j];
						double val = _ref->__mat[i * _nNode + j] * sc;
						mat->_mat[0].coeffRef(I, J) += val;
						//mat->_mat[0].data().valuePtr()[map[I*N+J]] += val;
					}
				}
			}
			else
			{
				for (int i = 0; i < _nNode; i++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						int I = index[i];
						int J = index[j];
						double val = _ref->__mat[i * _nNode + j] * sc;
						//mat->_mat[0].coeffRef(I, J) += val;
						mat->_mat[0].data().valuePtr()[map[I * N + J]] += val;
					}
				}
			}
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
		double MB(int i, int k)
		{
			double _val3 = 0;
			for (int j = 0; j < 2; j++) {
				_val3 += _Sij[(i << 1) + j] * this->get_gi(j, k);
			}
			return _val3 * _ref->_refDv;
		}
		//membrane boundary term
		double MB3(int i, int I)
		{
			double _val3 = 0;
			for (int j = 0; j < 2; j++) {
				_val3 += _Sij[(i << 1) + j] * _ref->d1[j][I];
			}
			return _val3;
		}
		//membrane boundary term
		double MB2(int i, int k)
		{
			double _val3 = 0;
			for (int j = 0; j < 2; j++) {
				_val3 += _Sij[(i << 1) + j] * this->get_gi(j, k);
			}
			return _val3;
		}
		double S(int i, int j) {
			return _Sij[(i << 1) + j];
		}
		//bending boundary term
		double KB(int j, int k, int l, int m, double _la, double _mu)
		{
			double _val3 = 0;


			for (int g = 0; g < 2; g++)
			{
				for (int h = 0; h < 2; h++)
				{
					double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));
					double D = (_ref->d2[g * 2 + h][j] - Gammaijk[(g * 2 + h) * 2 + 0] * _ref->d1[0][j] - Gammaijk[(g * 2 + h) * 2 + 1] * _ref->d1[1][j]);
					_val3 += A * N[k] * (D);
				}
			}

			return _val3 * _ref->refDv * _ref->refDv;
		}
		//membrane boundary term
		double HB(int j, int k, int l, int m, double _la, double _mu)
		{
			double _val3 = 0;


			for (int g = 0; g < 2; g++)
			{
				for (int h = 0; h < 2; h++)
				{
					double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));
					//double D = (_ref->d2[g * 2 + h][j] - Gammaijk[(g * 2 + h) * 2 + 0] * _ref->d1[0][j] - Gammaijk[(g * 2 + h) * 2 + 1] * _ref->d1[1][j]);
					double D = (_ref->d1[g][j] * get_gi(h, k) + _ref->d1[h][j] * get_gi(g, k));
					_val3 += A * (D);
				}
			}

			return 0.5 * _val3;// * _ref->refDv * _ref->refDv;
		}
		//angle term
		double T(int i, int s, int m) {
			double val = 0;
			for (int l = 0; l < 2; l++)
			{

				val += _ref->d1[l][i] * N[s] * _ref->get__Gij(l, m);
			}
			return val * _ref->refDv;
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
			const static int ll[2]{ 0,1 };

			double _val4 = 0;

			for (const auto& l : ll)
			{
				for (const auto& m : ll)
				{
					for (const auto& g : ll)
					{
						for (const auto& h : ll)
						{

							double A = (_la * _ref->get__Gij(h, g) * _ref->get__Gij(m, l) + 2 * _mu * _ref->get__Gij(h, m) * _ref->get__Gij(g, l));

							double FF = (_ref->d1[g][j] * get_gi(h, k) + _ref->d1[h][j] * get_gi(g, k));
							double GG = (_ref->d1[l][i] * get_gi(m, k2) + _ref->d1[m][i] * get_gi(l, k2));
							_val4 += A * FF * GG;
						}
					}
				}
			}
			return _val4 * _ref->refDv * 0.25;
		}
		void H(_mySparse* M, int64_t* _index, double _la, double _mu, double sc)
		{
			const static int kk[3]{ 0,1,2 };
			const static int ll[2]{ 0,1 };
			//static std::map<_mySparse*, std::vector<Eigen::Triplet<double>>> dict;

			for (int i = 0; i < _nNode; i++)
			{
				int I = _index[i] * 3;
				for (int j = 0; j < _nNode; j++)
				{
					int J = _index[j] * 3;
					for (const auto& k : kk)
					{
						for (const auto& k2 : kk)
						{
							double _val4 = 0;
							_val4 = H(i, k, j, k2, _la, _mu);
							M->_mat[0].coeffRef(I + k, J + k2) += _val4 * sc;
						}
					}
				}
			}
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
			if (mode == "SLOPE")
			{
				std::memcpy(__mem->_Gi, Gi2, sizeof(double) * 6);
				std::memcpy(__mem->_Gij, Gij2, sizeof(double) * 4);
				std::memcpy(__mem->_gij, gij2, sizeof(double) * 4);
				__mem->_gi[2] = 0;
				__mem->_gi[5] = 0;

			}
			else {
				std::memcpy(__mem->_gij, gij, sizeof(double) * 4);
				std::memcpy(__mem->_Gi, Gi, sizeof(double) * 6);
				std::memcpy(__mem->_Gij, Gij, sizeof(double) * 4);
			}
			//std::memset(__mem->_gij, 0, sizeof(double) * 4);
			//std::memset(__mem->_Gij, 0, sizeof(double) * 4);
			std::memcpy(__mem->_bij, bij, sizeof(double) * 12);
			std::memcpy(__mem->_Gammaijk, Gammaijk, sizeof(double) * 8);

		}
	};

	public ref class memS_ref {
	public:
		double _lambda;
		double lambda;
		double getlambda()
		{
			return lambda;
		}
		void setlambda(double val)
		{
			lambda = val;
		}
		double get_lambda()
		{
			return _lambda;
		}
		void set_lambda(double val)
		{
			_lambda = val;
		}
		_memS_ref* __mem = 0;
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
		double _gi(int l, int s) {
			return __mem->get__gi(l, s);
		}
		double _Gi(int l, int s) {
			return __mem->get__Gi(l, s);
		}
		memS_ref() {
			__mem = new _memS_ref();
		}
		void dispose() {
			if (__mem != 0)
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
				__mem->set_buf_z(Z->_arr->__v.data(), nNode);
			}
			if (phi != nullptr)
			{
				__mem->set_buf_phi(phi->_arr->__v.data(), nNode);
			}
		}
		void update_xi_eta(int nNode, KingOfMonsters::myDoubleArray^ xi, KingOfMonsters::myDoubleArray^ eta) {
			if (xi != nullptr)
			{
				__mem->set_buf_xi(xi->_arr->__v.data(), nNode);
			}
			if (eta != nullptr)
			{
				__mem->set_buf_eta(eta->_arr->__v.data(), nNode);
			}
			
		}
		void update_mu(int nNode, KingOfMonsters::myDoubleArray^ mu) {
			if (mu != nullptr)
			{
				__mem->set_buf_mu(mu->_arr->__v.data(), nNode);
			}

		}
		void update3(int nNode, KingOfMonsters::myDoubleArray^ node, KingOfMonsters::myDoubleArray^ weights, KingOfMonsters::myDoubleArray^ def, bool ignoreZ) {
			if (node != nullptr) {
				__mem->set_node(node->_arr->__v.data(), nNode * 3);
			}
			if (def != nullptr)
			{
				__mem->set_def(def->_arr->__v.data(), nNode);
			}
		}
		void update3(int nNode, KingOfMonsters::myDoubleArray^ node, KingOfMonsters::myDoubleArray^ weights, KingOfMonsters::myDoubleArray^ def) {
			update3(nNode, node, weights, def, false);
		}
		void update(int nNode, int uDim, int vDim) {
			__mem->update(nNode, uDim, vDim);
		}
	};
	public ref class memS
	{
	public:
		_memS* __mem = 0;
	public:
		double x, y, z, Z, phi;
		double _x, _y, _z, _Z, _phi;
		double xi, eta,mu,nu;
		double refDv;
		double _refDv;
		double dv;
		double _dv;
	public:
		double orient(memS^ another) {
			double dot = another->__mem->N[0] * this->__mem->N[0] + another->__mem->N[1] * this->__mem->N[1] + another->__mem->N[2] * this->__mem->N[2];
			if (dot < 0)return -1;
			return 1;
		}
		void setRef(memS_ref^ _mem_) {
			__mem->_ref = _mem_->__mem;
			__mem->_ref->RAM = __mem->RAM;


		}
		double _Gi2(int i, int s)
		{
			return __mem->Gi2[i * 3 + s];
		}
		double _gi(int i, int s)
		{
			return __mem->gi[i * 3 + s];
		}
		double Gamma000()
		{
			return __mem->Gamma000();
		}
		double Gamma001()
		{
			return __mem->Gamma001();
		}
		double Gamma010()
		{
			return __mem->Gamma010();
		}
		double Gamma011()
		{
			return __mem->Gamma011();
		}
		double Gamma110()
		{
			return __mem->Gamma110();
		}
		double Gamma111()
		{
			return __mem->Gamma111();
		}
		double f2uu()
		{
			return __mem->f2uu();
		}
		double f2uv()
		{
			return __mem->f2uv();
		}
		double f2vv()
		{
			return __mem->f2vv();
		}
		double f1u()
		{
			return __mem->f1u();
		}
		double f1v()
		{
			return __mem->f1v();
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
		void G3(myDoubleArray^ vec, myIntArray^ index, double sc) {
			__mem->G3(vec->_arr, index->data(), sc);
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
		void MASS(mySparse^ mat, myIntArray^ index, double sc) {
			__mem->MASS(mat->dat, index->data(), sc);
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
		void U_phi_z(mySparse^ mat, myIntArray^ index, int N, myIntArray2^ map, double R)
		{
			if (map == nullptr)
				__mem->U_phi_z(mat->dat, index->_arr, N, 0, R);
			else
				__mem->U_phi_z(mat->dat, index->_arr, N, map->_arr, R);
		}
		void D_phi_z(mySparse^ mat, myIntArray^ index, int N, myIntArray2^ map, double sc)
		{
			if (map == nullptr)
				__mem->D_phi_z(mat->dat, index->_arr, N, 0, sc);
			else
				__mem->D_phi_z(mat->dat, index->_arr, N, map->_arr, sc);
		}

		double theta(bool accurate)
		{
			return __mem->theta(accurate);
		}
		void theta_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, bool accurate)
		{
			__mem->theta_phi(__mem->__grad, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void theta_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, bool accurate)
		{
			__mem->theta_Z(__mem->__grad, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		double curvature(double v1, double v2, bool accurate)
		{
			return __mem->curvature(v1, v2, accurate);
		}
		void  curvature_Z(myDoubleArray^ grad, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->curvature_Z(__mem->__grad, v1, v2, accurate);
			for (int i = 0; i < __mem->_nNode; i++)
			{
				grad->_arr->__v(index->data()[i]) += __mem->__grad[i] * c1 ;
			}
			//mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		
		double curvature1(double v1, double v2, bool accurate)
		{
			return __mem->curvature(v1, v2, accurate);
		}
		void  curvature1_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->curvature_xi(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void  curvature1_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->curvature_eta(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double curvature2(double v1, double v2, bool accurate)
		{
			return __mem->curvature(v1, v2, accurate);
		}
		void  curvature2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->curvature_xi(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void  curvature2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->curvature_eta(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double guideBC(double v1, double v2, bool accurate)
		{
			return __mem->guideBC(v1, v2, accurate,0);
		}
		void guideBC_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guideBC_xi(__mem->__grad, v1, v2, accurate, 0);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guideBC_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guideBC_eta(__mem->__grad, v1, v2, accurate, 0);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double guide2(double v1, double v2, bool accurate)
		{
			return __mem->guide2(v1, v2, accurate);
		}
		void guide2_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guide2_z(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		
		double guide3(double v1, double v2, bool accurate)
		{
			return __mem->guide3(v1, v2, accurate);
		}
		void guide3_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guide3_z(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		
		double guide(double v1, double v2, bool accurate,int mode)
		{
			return __mem->guide(v1, v2,  accurate,mode);
		}
		void guide_xi(mySparse^ mat, int ii, myIntArray^ index, double sc,double c1,double v1, double v2, bool accurate,int mode)
		{
			__mem->guide_xi(__mem->__grad,v1, v2, accurate,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
			memcpy(__mem->__grad_z_tmp, __mem->__grad, sizeof(double) * __mem->_nNode);
		}
		void guide_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate,int mode)
		{
			__mem->guide_eta(__mem->__grad,v1, v2, accurate,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
			memcpy(__mem->__grad_phi_tmp, __mem->__grad, sizeof(double) * __mem->_nNode);
		}
	
		void guide_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guide_Z(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		double guide4(double v1, double v2, bool accurate)
		{
			return __mem->guide4(v1, v2, accurate);
		}
		void guide4_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guide4_Z(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		double guide5(double v1, double v2, bool accurate)
		{
			return __mem->guide5(v1, v2, accurate);
		}
		void guide5_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,  bool accurate,bool remove, [Runtime::InteropServices::Out] double% lambda )
		{
			__mem->guide5_Z(__mem->__grad,  accurate);
			if (remove)
			{
				__mem->remove3(__mem->_nNode, __mem->__grad, __mem->__grad_phi_tmp, __mem->__grad_z_tmp);
			}
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide5_Z(myDoubleArray^ grad,myIntArray^ index, double c1, bool accurate,bool remove)
		{
			__mem->guide5_Z(__mem->__grad, accurate);
			double lambda = 1;
			if(remove)
			 __mem->remove3(__mem->_nNode, __mem->__grad, __mem->__grad_phi_tmp, __mem->__grad_z_tmp);
			for (int i = 0; i < __mem->_nNode; i++)
			{
				grad->_arr->__v(index->data()[i]) += __mem->__grad[i]*c1;
			}

			//mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		double guide6( bool accurate)
		{
			return __mem->guide6( accurate);
		}
		void guide6_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,  bool accurate, bool remove)
		{
			__mem->guide6_Z(__mem->__grad, accurate);
			
			if (remove)
				__mem->remove3(__mem->_nNode, __mem->__grad, __mem->__grad_phi_tmp, __mem->__grad_z_tmp);
			/*for (int i = 0; i < __mem->_nNode; i++)
			{
				grad->_arr->__v(index->data()[i]) += __mem->__grad[i] * c1;
			}*/

			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		
		}
		double guide7(bool accurate)
		{
			return __mem->guide7(accurate);
		}
		void guide7_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, bool accurate, bool remove)
		{
			__mem->guide7_Z(__mem->__grad, accurate);

			if (remove)
				__mem->remove3(__mem->_nNode, __mem->__grad, __mem->__grad_phi_tmp, __mem->__grad_z_tmp);
			/*for (int i = 0; i < __mem->_nNode; i++)
			{
				grad->_arr->__v(index->data()[i]) += __mem->__grad[i] * c1;
			}*/

			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);

		}
		double unit()
		{
			return __mem->unit();
		}
		void unit_mu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->unit_mu(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void unit_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->unit_nu(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);

		}
		double get_v1()
		{
			return __mem->get_v1();
		}
		double get_v2()
		{
			return __mem->get_v2();
		}

		double fair()
		{
			return __mem->fair();
		}
		void fair_mu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->fair_mu(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		
		double curl_free()
		{
			return __mem->curl_free();
		}
		void curl_free_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->curl_free_xi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void curl_free_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->curl_free_eta(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double vec(double V1, double V2)
		{
			return __mem->vec(V1, V2);
		}
		void vec_mu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,double V1,double V2)
		{
			__mem->vec_mu(__mem->__grad,V1,V2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void vec_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double V1, double V2)
		{
			__mem->vec_nu(__mem->__grad,V1, V2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double guide8(bool accurate,double v1,double v2,int mode)
		{
			return __mem->guide8(accurate,v1,v2,mode);
		}
	
		void guide8_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate,int mode)
		{
			__mem->guide8_xi(__mem->__grad, v1, v2, accurate,mode);

			//if (remove)
			//	__mem->remove3(__mem->_nNode, __mem->__grad, __mem->__grad_phi_tmp, __mem->__grad_z_tmp);
			/*for (int i = 0; i < __mem->_nNode; i++)
			{
				grad->_arr->__v(index->data()[i]) += __mem->__grad[i] * c1;
			}*/

			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);

		}
		void guide8_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate,int mode)
		{
			__mem->guide8_eta(__mem->__grad, v1, v2, accurate,mode);

			//if (remove)
			//	__mem->remove3(__mem->_nNode, __mem->__grad, __mem->__grad_phi_tmp, __mem->__grad_z_tmp);
			/*for (int i = 0; i < __mem->_nNode; i++)
			{
				grad->_arr->__v(index->data()[i]) += __mem->__grad[i] * c1;
			}*/

			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);

		}
		/*void guide_xi(double v1, double v2, bool accurate)
		{
			__mem->guide_xi( v1, v2, accurate);
			//mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_eta( double v1, double v2,  bool accurate)
		{
			__mem->guide_eta(v1, v2,  accurate);
			//mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}*/
	
		/*void guide_write(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			mat->dat->addrow(ii, index->_arr, __mem->__guide2_xieta, 0, sc, __mem->_nNode * 2, true, c1);
		}*/
		
		void mix_write(mySparse^ mat,  int ii, myIntArray^ index, myIntArray^ index2,double v1,double v2,double w1,double w2, double sc, double c1, bool accurate, bool remove)
		{
			__mem->mix_xi(__mem->__grad_xi , v1, v2, w1, w2, accurate);
			__mem->mix_eta(__mem->__grad_eta, v1, v2, w1, w2, accurate);
			//if(remove)
			//__mem->remove2(__mem->_nNode,__mem->__grad_xi, __mem->__grad_z_tmp, __mem->__grad_eta, __mem->__grad_phi_tmp, sc);
			mat->dat->addrow(ii, index->_arr, __mem->__grad_xi, 0, sc, __mem->_nNode, true, c1);
			mat->dat->addrow(ii, index2->_arr, __mem->__grad_eta, 0, sc, __mem->_nNode, false, c1);
		}
		void mix_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix_xi(__mem->__grad,v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void mix_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix_eta(__mem->__grad,v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
	
		double guide_supported(double v1, double v2, bool accurate,int mode)
		{
			return __mem->guide_supported(v1, v2, accurate,mode);
		}
		void guide_xi_supported(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate, int mode)
		{
			__mem->guide_supported_xi(__mem->__grad, v1, v2, accurate, mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_eta_supported(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate, int mode)
		{
			__mem->guide_supported_eta(__mem->__grad, v1, v2, accurate, mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		

		double guide_free(double v1, double v2, bool accurate,int mode)
		{
			return __mem->guide_free(v1, v2, accurate,mode);
		}
		void guide_xi_free(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate,int mode)
		{
			__mem->guide_free_xi(__mem->__grad, v1, v2, accurate,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_eta_free(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate, int mode)
		{
			__mem->guide_free_eta(__mem->__grad, v1, v2, accurate, mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		
		void guide_z_free(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate,int mode)
		{
			__mem->guide_free_z(__mem->__grad, v1, v2, accurate,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
	
		double testequal(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->testequal(v1, v2, w1, w2, accurate);
		}
		void testequal_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->testequal_xi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void testequal_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->testequal_eta(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double E22(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->E22(v1, v2, w1, w2, accurate);
		}
		void E22_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->E22_xi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void E22_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->E22_eta(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double setEij_A(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->setEij_A(v1, v2, w1, w2, accurate);
		}
		void setEij_A_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->setEij_A_Z(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void setEij_A_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->setEij_A_xi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad,0, sc, __mem->_nNode, true,c1);
		}
		void setEij_A_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->setEij_A_eta(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		double setEij_B(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->setEij_B(v1, v2, w1, w2, accurate);
		}
		void setEij_B_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->setEij_B_xi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void setEij_B_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->setEij_B_eta(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double mix_guide(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->mix_guide(v1, v2, w1, w2, accurate);
		}
		
		void mix_guide_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix_guide_phi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void mix_guide_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix_guide_Z(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		double mix(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->mix(v1, v2, w1, w2, accurate);
		}
		double e11(double cx,double cy,double L)
		{
			return __mem->e11(cx,cy,L);
		}
		void e11_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double cx, double cy, double L)
		{
			__mem->e11_xi(cx, cy, L,__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad,0, sc, __mem->_nNode,true, c1);
		}
		void e11_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double cx, double cy, double L)
		{
			__mem->e11_eta(cx, cy, L,__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		void mix_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix_phi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
			//memcpy(__mem->__grad_z_tmp , __mem->__grad, sizeof(double) * __mem->_nNode);
			for (int i = 0; i < __mem->_nNode;i++)
			{
				__mem->__grad_z_tmp[i] = __mem->__grad[i] * c1;
			}
		}
		void mix_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix_Z(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		double mix2(double v1, double v2, double w1, double w2, bool accurate)
		{
			return __mem->mix2(v1, v2, w1, w2, accurate);
		}

		void mix2_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix2_phi(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void mix2_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, bool accurate)
		{
			__mem->mix2_Z(__mem->__grad, v1, v2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
	
		double BCEQ(double _t1, double _t2, double s11,double a,double b,bool accurate)
		{
			return __mem->BCEQ(_t1, _t2,s11,a,b,accurate);
		}
		void BCEQ_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double _t1, double _t2,double  stt,double a,double b,bool accurate)
		{
			__mem->BCEQ_phi(__mem->__grad, _t1, _t2, stt,a,b, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void BCEQ_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double _t1, double _t2,double stt,double a,double b,bool accurate)
		{
			__mem->BCEQ_z(__mem->__grad, _t1,_t2, stt, a,b,accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		double mix_BC(double v1, double v2,  bool accurate)
		{
			return __mem->mix_BC(v1, v2, accurate);
		}
		void mix_BC_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->mix_BC_phi(__mem->__grad, v1, v2,accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void mix_BC_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->mix_BC_Z(__mem->__grad, v1, v2,accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void conjugate(double v1, double v2, [Runtime::InteropServices::Out] double% s1, [Runtime::InteropServices::Out] double% s2, bool accurate)
		{
			double _s1 = 0; double _s2 = 0;
			__mem->conjugate(v1, v2, &_s1, &_s2, accurate);
			s1 = _s1;
			s2 = _s2;
		}
		void transpose(double v1, double v2, [Runtime::InteropServices::Out] double% s1, [Runtime::InteropServices::Out] double% s2, bool accurate)
		{
			double _s1 = 0; double _s2 = 0;
			__mem->transpose(v1, v2, &_s1, &_s2, accurate);
			s1 = _s1;
			s2 = _s2;
		}
		double tt(double v1, double v2, bool accurate) {

			return __mem->tt(v1, v2, accurate);
		}
		double st(double v1, double v2, bool accurate)
		{

			return __mem->st(v1, v2,  accurate);
		}
		void st_z(myDoubleArray^ arr, double v1, double v2, bool accurate)
		{

			__mem->st_z(arr->_arr->__v.data(), v1, v2, accurate);

		}
		void st_phi(myDoubleArray^ arr, double v1, double v2, bool accurate)
		{

			__mem->st_phi(arr->_arr->__v.data(), v1, v2, accurate);

		}
		double st2(double v1, double v2, bool accurate)
		{

			return __mem->st2(v1, v2, accurate);
		}
		void st2_phi(myDoubleArray^ arr, double v1, double v2, bool accurate)
		{

			__mem->st2_phi(arr->_arr->__v.data(), v1, v2, accurate);

		}
		void tt_z(myDoubleArray^ arr, double v1, double v2, bool accurate)
		{
			__mem->tt_z(arr->_arr->__v.data(), v1, v2, accurate);
		}
		
		double tt_2(double v1, double v2, bool accurate) {

			return __mem->tt_2(v1, v2, accurate);
		}
		
		void tt2_phi(myDoubleArray^ arr, double v1, double v2, bool accurate)
		{
			__mem->tt2_phi(arr->_arr->__v.data(), v1, v2, accurate);
		}
		
		double CC(double v1, double v2, bool accurate)
		{
			return __mem->CC(v1, v2, accurate);
		}
		void CC_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->CC_phi(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void CC_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->CC_Z(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		double D1()
		{
			return __mem->D1();
		}
		double D2()
		{
			return __mem->D2();
		}
		void D_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double c2)
		{
			__mem->D1_phi(__mem->__grad);
			__mem->D2_phi(__mem->__grad2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, __mem->__grad2, sc, __mem->_nNode, c1, c2);
		}
		void D_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double c2)
		{
			__mem->D1_z(__mem->__grad);
			__mem->D2_z(__mem->__grad2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, __mem->__grad2, sc, __mem->_nNode, c1, c2);
		}

		void d0(mySparse^ mat, int ii, myIntArray^ index, double sc)
		{
			mat->dat->addrow(ii, index->_arr, __mem->_ref->d0, sc, __mem->_nNode);
		}
		void U_z(mySparse^ mat, int ii, myIntArray^ index, double sc)
		{

			mat->dat->addrow(ii, index->_arr, __mem->__grad_z, sc, __mem->_nNode);
			//memcpy(__mem->__grad_z_tmp, __mem->__grad_z, sizeof(double) * __mem->_nNode);

		}
		void U_mix(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double c2)
		{

			mat->dat->addrow(ii, index->_arr, __mem->__grad_phi, __mem->__grad_phi, sc, __mem->_nNode, c1, c2);

		}
		void U_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
		{
			//__U_z(__mem->__grad_z, true);
			__mem->__U_z(__mem->__grad_z, true);
			mat->dat->addrow(ii, index->_arr, __mem->__grad_z, sc, __mem->_nNode, coeff);
			//memcpy(__mem->__grad_z_tmp, __mem->__grad_z, sizeof(double) * __mem->_nNode);

		}
		void U_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff,double load,bool accurate)
		{
			__mem->__U_phi(__mem->__grad_phi, load,accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad_phi, sc, __mem->_nNode, coeff);
			memcpy(__mem->__grad_phi_tmp, __mem->__grad_phi, sizeof(double) * __mem->_nNode);

		}
		void U_z(array<double>^ X) {
			System::Runtime::InteropServices::Marshal::Copy((IntPtr)__mem->__grad_z, X, 0, __mem->_nNode);
			memcpy(__mem->__grad_z_tmp, __mem->__grad_z, sizeof(double) * __mem->_nNode);
		}
		void U_phi(array<double>^ x) {
			System::Runtime::InteropServices::Marshal::Copy((IntPtr)__mem->__grad_phi, x, 0, __mem->_nNode);
			memcpy(__mem->__grad_phi_tmp, __mem->__grad_phi, sizeof(double) * __mem->_nNode);
		}
		void U_z(myDoubleArray^ X) {
			memcpy(X->_arr->__v.data(), __mem->__grad_z, sizeof(double) * __mem->_nNode);
			memcpy(__mem->__grad_z_tmp, __mem->__grad_z, sizeof(double) * __mem->_nNode);
		}
		void U_phi(myDoubleArray^ x) {
			memcpy(x->_arr->__v.data(), __mem->__grad_phi, sizeof(double) * __mem->_nNode);
			memcpy(__mem->__grad_phi_tmp, __mem->__grad_phi, sizeof(double) * __mem->_nNode);
		}
		array<double>^ fM(double _la, double _mu)
		{
			double a, b, c, d;
			__mem->fM(_la, _mu, &a, &b, &c, &d);
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
		double shear_z(int uv)
		{
			return __mem->shear_z(uv);
		}
		double shear_phi(int uv)
		{
			return __mem->shear_phi(uv);
		}
		void shear(myDoubleArray^ arr, int uv)
		{
			__mem->shear(arr->_arr->__v.data(), uv);
		}
		double SLOPE_phi(int l) {
			return __mem->SLOPE_phi(l);
		}
		double SLOPE_z(int l) {
			return __mem->SLOPE_z(l);
		}
		double SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {
			return __mem->SLOPE_phi(dcdtstar0, dcdtstar1, sc);
		}
		double SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {
			return __mem->SLOPE_z(dcdtstar0, dcdtstar1, sc);
		}
		double SLOPE_mu(double dcdtstar0, double dcdtstar1) {
			return __mem->SLOPE_mu(dcdtstar0, dcdtstar1);
		}
		double SLOPE_xi(double dcdtstar0, double dcdtstar1) {
			return __mem->SLOPE_xi(dcdtstar0, dcdtstar1);
		}
		double SLOPE_eta(double dcdtstar0, double dcdtstar1) {
			return __mem->SLOPE_eta(dcdtstar0, dcdtstar1);
		}
		double SLOPE(int l, int i) {
			return __mem->SLOPE(l, i);
		}
		double constant_SLOPE_z(double dcdtstar0, double dcdtstar1, bool accurate) {
			return __mem->constant_SLOPE_z(dcdtstar0, dcdtstar1, accurate);
		}
		double constant_SLOPE_phi(double dcdtstar0, double dcdtstar1, bool accurate) {
			return __mem->constant_SLOPE_phi(dcdtstar0, dcdtstar1, accurate);
		}

		void constant_SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc, bool accurate)
		{
			__mem->constant_SLOPE(arr->_arr->__v.data(), dcdt1, dcdt2, sc, accurate);

		}
		void constant_SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc, int shift, bool accurate)
		{
			__mem->constant_SLOPE(arr->_arr->__v.data() + shift, dcdt1, dcdt2, sc, accurate);

		}
		void _SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2,double sc) {
			__mem->_SLOPE(arr->_arr->__v.data(), dcdt1, dcdt2,sc);
			//return __mem->SLOPE(l, i);
		}
		void _SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, int shift,double sc) {
			__mem->_SLOPE(arr->_arr->__v.data()+shift, dcdt1, dcdt2,sc);
			//return __mem->SLOPE(l, i);
		}
		void SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc) {
			__mem->SLOPE(arr->_arr->__v.data(), dcdt1, dcdt2, sc);
			//return __mem->SLOPE(l, i);
		}
		void SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc, int shift) {
			__mem->SLOPE(arr->_arr->__v.data() + shift, dcdt1, dcdt2, sc);
			//return __mem->SLOPE(l, i);
		}
		double eM(double _la, double _mu) {
			return __mem->eM(_la, _mu);
		}
		double eK(double _la, double _mu) {
			return __mem->eK(_la, _mu);
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
		double d2(int l, int i) {
			return __mem->_d2(l, i);
		}
		double gi(int i, int s) {
			return __mem->get_gi(i, s);
		}
		double Gi(int i, int s) {
			return __mem->get_Gi(i, s);
		}
		///error norm
		double error() {
			return __mem->error();
		}
		///stress function
		double F(int i, int j) {
			return __mem->F(i, j);
		}
		void F(mySparse^ mat, myIntArray^ index, double sc) {
			__mem->F(mat->dat, index->data(), sc);
		}
		//least squares
		double F2(int i, int j) {
			return __mem->F2(i, j);
		}
		System::String^ F2(mySparse^ mat, myIntArray^ index, double sc) {
			auto str = __mem->F2(mat->dat, index->data(), sc);
			return gcnew System::String(str.c_str());
		}
		System::String^ F2A(mySparse^ mat, myIntArray^ index, double sc) {
			auto str = __mem->F2A(mat->dat, index->data(), sc);
			return gcnew System::String(str.c_str());
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
		double Dc(double a, double b, int l) {
			return __mem->Dc(a, b, l);
		}
		double get_a(int l) {
			return __mem->__a(l);
		}
		double get_b(int l) {
			return __mem->__b(l);
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
		void F4_z(myDoubleArray^ arr, double dcdtstar0, double dcdtstar1, double sc)
		{
			__mem->F4_z(arr->_arr->__v.data(), dcdtstar0, dcdtstar1, sc);
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
		double K(int i, int s, int j, int k, double _la, double _mu) {
			return __mem->K(i, s, j, k, _la, _mu);
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
		double MB3(int l, int i)
		{
			return __mem->MB3(l, i);
		}
		///<summary>
		///bending boundary term
		///</summary>
		double KB(int j, int k, int l, int m, double _la, double _mu) {
			return __mem->KB(j, k, l, m, _la, _mu);
		}
		double HB(int j, int k, int l, int m, double _la, double _mu) {
			return __mem->HB(j, k, l, m, _la, _mu);
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
		void H(mySparse^ M, myIntArray^ index, double _la, double _mu, double sc)
		{
			__mem->H(M->dat, index->data(), _la, _mu, sc);
		}
		void K(mySparse^ M, myIntArray^ index, double _la, double _mu, double sc)
		{
			__mem->K(M->dat, index->data(), _la, _mu, sc);
		}
		void S(double val1, double val2, double val3) {
			__mem->set_Sij(val1, val2, val3);
		}
		void getmat(myIntArray^ index, workspace^ _dat) {
			__mem->getMat(index->data(), _dat->_dat);
		}
		void getmat_Galerkin(myIntArray^ index, workspace^ _dat, int i, double w) {
			__mem->getMat_Galerkin(index->data(), _dat->_dat, i, w);
		}
		void getmatG(myIntArray^ index, myIntArray^ index2, workspace^ _dat, double sc) {
			__mem->getMatG(index->data(), index2->data(), _dat->_dat, sc);
		}
		void getmat_d0(myIntArray^ index, workspace^ _dat) {
			__mem->getMat_d0(index->data(), _dat->_dat);
		}
		void getmat_slope(myIntArray^ index, workspace^ _dat, double dcdt1, double dcdt2, bool airy) {
			__mem->getMat_slope(index->data(), _dat->_dat, dcdt1, dcdt2, airy);
		}
		double detZ()
		{
			return __mem->detZ();
		}
		double detphi()
		{
			return __mem->detphi();
		}
		double _detZ()
		{
			return __mem->_detZ();
		}
		double _detZ_a()
		{
			return __mem->_detZ_a();
		}
		double _detZ_b()
		{
			return __mem->_detZ_b();
		}
		double _detZ_c()
		{
			return __mem->_detZ_c();
		}
		double _detphi()
		{
			return __mem->_detphi();
		}
		double bodyF(double load,bool accurate) {
			//double sc = 1.0/__mem->_ref->_refDv/ __mem->_ref->_refDv;
			if(accurate)
			{
				return __mem->sc * __mem->bodyF-load* __mem->dv/ __mem->_ref->_refDv;
			}
			else
			{
				return __mem->sc * __mem->bodyF - load;
			}
		}
		double BCF(int l) {
			//double sc = 1.0 / __mem->_ref->_refDv / __mem->_ref->_refDv;
			return  __mem->sc * __mem->BCF[l];
		}
		double BCF2(int l) {
			//double sc = 1.0 / __mem->_ref->_refDv / __mem->_ref->_refDv;
			return __mem->sc * __mem->BCF2[l];
		}
		double BCF6(int l) {
			return __mem->BCF6[l];
		}
		double dcdtstar(double x, double y, int i)
		{
			return x * __mem->get_gi(i, 0) + y * __mem->get_gi(i, 1);
		}
		double norm(double d1, double d2) {
			double x = __mem->get_Gi(0, 0) * d1 + __mem->get_Gi(1, 0) * d2;
			double y = __mem->get_Gi(0, 1) * d1 + __mem->get_Gi(1, 1) * d2;
			double z = __mem->get_Gi(0, 2) * d1 + __mem->get_Gi(1, 2) * d2;
			return sqrt(x * x + y * y + z * z);
		}
		void compute() {
			if (__mem->RAM == SAVE)
			{
				__mem->_ref->__mat = __mem->__mat;
				__mem->_ref->d0 = __mem->d0;

				__mem->_ref->d1 = __mem->d1;
				__mem->_ref->d2 = __mem->d2;
				__mem->_ref->d2_star = __mem->d2_star;
				__mem->_ref->B = __mem->B;
				__mem->_ref->tt0 = __mem->tt0;
				__mem->_ref->hh0 = __mem->hh0;
				__mem->_ref->tt1 = __mem->tt1;
				__mem->_ref->hh1 = __mem->hh1;
				__mem->_ref->tt2 = __mem->tt2;
				__mem->_ref->hh2 = __mem->hh2;

			}
			__mem->update2();
			x = __mem->x;
			y = __mem->y;
			z = __mem->z;
			Z = __mem->Z;
			phi = __mem->phi;
			dv = __mem->dv;
			_dv = __mem->_dv;
			_refDv = __mem->_dv;
			refDv = __mem->dv;
			xi = __mem->xi;
			eta = __mem->eta;
			mu = __mem->mu;
			//nu = __mem->nu;
		}
		void update_lo(array<double>^ lo) {
			__mem->set_lo(lo[0], lo[1]);
		}
		void update_elem_ref(int nNode, int uDim, int vDim, array<double, 3>^ M, array<int, 2>^ dd) {
			__mem->update_ref(nNode, uDim, vDim);
			if (__mem->RAM == SAVE)
			{
				__mem->_ref->M[0] = __mem->M[0];
				__mem->_ref->M[1] = __mem->M[1];
				__mem->_ref->dd = __mem->dd;

			}
			if (__mem->RAM == MAX && (!__mem->_ref->initialized))
			{
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
		}
		void update_elem(int nNode, int uDim, int vDim, array<double, 3>^ M, array<int, 2>^ dd) {
			__mem->update(nNode, uDim, vDim);

			if (__mem->RAM == SAVE)
			{
				int i = 0;
				for (int j = 0; j < uDim; j++) {
					for (int k = 0; k < uDim; k++) {
						__mem->set_M_save(i, j, k, M[i, j, k]);
					}
				}
				i = 1;
				for (int j = 0; j < vDim; j++) {
					for (int k = 0; k < vDim; k++) {
						__mem->set_M_save(i, j, k, M[i, j, k]);
					}
				}
				for (int i = 0; i < nNode; i++) {
					for (int j = 0; j < 2; j++) {
						__mem->set_dd_save(i, j, dd[i, j]);
					}
				}
			}
		}
		void memory(memS_ref^ _mem_) {
			__mem->memory(_mem_->__mem);
			_x = __mem->_ref->_x;
			_y = __mem->_ref->_y;
			_z = __mem->_ref->_z;
			Z = __mem->_ref->Z;
			//phi = __mem->phi;
		}
		memS(System::String^ RAM) {
			_RAM __RAM = SAVE;
			if (RAM == "SAVE")__RAM = SAVE;
			if (RAM == "MAX")__RAM = MAX;
			__mem = new _memS(std::string(""), __RAM);//MAX=cosume memory more but fast, SAVE=save memory and a little slow
		}
		void static func(int i) {

		}
		memS(System::String^ ultimate, System::String^ RAM) {

			_RAM __RAM = SAVE;
			if (RAM == "SAVE")__RAM = SAVE;
			if (RAM == "MAX")__RAM = MAX;

			if (ultimate == "U")
				__mem = new _memS("U", __RAM);
			else if (ultimate == "SLOPE")
				__mem = new _memS("SLOPE", __RAM);
			else if (ultimate == "SHELL")
				__mem = new _memS("SHELL", __RAM);
			else if (ultimate == "SENSITIVITY")
				__mem = new _memS("SENSITIVITY", __RAM);
			else
				__mem = new _memS("", __RAM);
		}
		void dispose() {
			if (__mem != 0)
				delete __mem;
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
