#pragma once
#include <cmath>
#include<vector>

using namespace System;
#include <cstring>

using std::vector;
using std::string;

#include "mySparseLibrary.h"
//#define EIGEN_DONT_PARALLELIZE
//#define EIGEN_MALLOC_ALREADY_ALIGNED  0
//double globalratio = 1;
namespace KingOfMonsters {
	
	public class _buffer {
	public:
		double mem[12000];
	};
	public ref class buffer {
	public:
		/*static void setratio(double r)
		{
			globalratio = r;
		}*/
		_buffer* _buf = 0;
		inline double at(int i)
		{
			return _buf->mem[i];
		}
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
		double oG11, oG12, oG22,og11,og12,og22,oX,oY,osc, orefDv;
		int _nNode;
		int dim[2]{ 0,0 };
		int _uDim, _vDim;
		double* node = 0;
		double* def = 0;
		double* buf_z = 0;
		double* buf_phi = 0;
		double* buf_xi = 0;
		double* buf_eta = 0;
		double* buf_nu = 0;
		double* buf_u = 0;
		double* buf_v = 0;
		double* buf_w = 0;
		double xi_1,xi_2, eta_1,eta_2;
		double _gi[6];
		double _Gi[6];
		double _ogi[6];
		double _oGi[6];
		double _gij[4];
		double _Gij[4];
		double _bij[12];
		double _Sij[4];
		double _Gammaijk[8];
		double oGammaijk[8];
		double _gammaijk[8];
		double* __dh[4]{ 0,0,0,0 };

		
		double* __mat = 0;
		double* __matF_xi = 0;
		double* __matF_eta = 0;
		double* __matF_phi = 0;
		double* __matD_xi = 0;
		double* __matD_eta = 0;
		double* __matD_phi = 0;
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
		double** d3 = 0;
		double** B = 0; //4
		double** tt0 = 0, ** hh0 = 0;//2
		double** tt1 = 0, ** hh1 = 0; //4
		double** tt2 = 0, ** hh2 = 0;//8
		double** tt3 = 0, ** hh3 = 0;//8

	public:
		bool initialized = false;
		double refDv = 0, _refDv = 0,load_dv=0;
		double _x = 0, _y = 0, _z = 0, __z = 0, Z = 0, _Z = 0, _u,_v=0,_w;
		double _xi = 0, _eta = 0,_nu=0;
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
			def = &buf[1000];
			buf_z = &buf[2000];
			buf_phi = &buf[3000];
			buf_xi = &buf[4000];
			buf_eta = &buf[5000];		
			buf_nu = &buf[6000];
			buf_u = &buf[7000];
			buf_v = &buf[8000];
			buf_w = &buf[9000];
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
		inline void set_buf_u(const int& i, const double& val) {
			buf_u[i] = val;
		}
		inline void set_buf_u(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_u, ptr, sizeof(double) * N);
		}
		inline void set_buf_v(const int& i, const double& val) {
			buf_v[i] = val;
		}
		inline void set_buf_v(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_v, ptr, sizeof(double) * N);
		}
		inline void set_buf_w(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_w, ptr, sizeof(double) * N);
		}
		inline void set_buf_nu(double* ptr, const int& N) {
			//buf_z[i] = val;
			memcpy(buf_nu, ptr, sizeof(double) * N);
		}
		inline void clear_buf_w(const int& N)
		{
			memset(buf_w, 0, sizeof(double) * N);
		}
		inline void clear_buf_nu(const int& N)
		{
			memset(buf_nu, 0, sizeof(double) * N);
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
		inline double& get__gammaijk(const int& i, const int& j, const int& k) {
			return _gammaijk[(((i << 1) + j) << 1) + k];
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
			d3 = 0;
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
			tt3 = 0;
			hh3 = 0;

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
				if (__mat != 0) { delete[] __mat; }
				__mat = 0;

				if (__matF_phi != 0) { delete[] __matF_phi; }
				__matF_phi = 0;
				if (__matF_xi != 0) { delete[] __matF_xi; }
				__matF_xi = 0;
				if (__matF_eta != 0) { delete[] __matF_eta; }
				__matF_eta = 0;


				if (__matD_phi != 0) { delete[] __matD_phi; }
				__matD_phi = 0;
				if (__matD_xi != 0) { delete[] __matD_xi; }
				__matD_xi = 0;
				if (__matD_eta != 0) { delete[] __matD_eta; }
				__matD_eta = 0;

				if (__dh[0] != 0)delete[] __dh[0];
				if (__dh[1] != 0)delete[] __dh[1];
				if (__dh[2] != 0)delete[] __dh[2];
				if (__dh[3] != 0)delete[] __dh[3];
				
				__dh[0] = 0;
				__dh[1] = 0;
				__dh[2] = 0;
				__dh[3] = 0;
			

				if (d0 != 0) {
					delete[] d0;
					d0 = 0;
				}
				if (d1 != 0)
				{
					delete[] d1[0];
					delete[] d1[1];					
					delete[] d1;
					d1 = 0;

				}
				if (d2 != 0)
				{
					delete[] d2[0];
					delete[] d2[1];
					delete[] d2[2];
					delete[] d2[3];
					
					delete[] d2;
					d2 = 0;
					delete[] d2_star[0];
					delete[] d2_star[1];
					delete[] d2_star[2];
					delete[] d2_star[3];
					
					delete[] d2_star;
					d2_star = 0;
				}
				if (d3 != 0)
				{
					delete[] d3[0];
					delete[] d3[1];
					delete[] d3[2];
					delete[] d3[3];
					delete[] d3[4];
					delete[] d3[5];
					delete[] d3[6];
					delete[] d3[7];
					delete[] d3;
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
					delete[] tt0[1];
					delete[] tt0;
					tt0 = 0;

					delete[] hh0[0];
					delete[] hh0[1];
					delete[] hh0;
					hh0 = 0;

					delete[] tt1[0];
					delete[] tt1[1];
					delete[] tt1[2];
					delete[] tt1[3];
					delete[] tt1;
					tt1 = 0;

					delete[] hh1[0];
					delete[] hh1[1];
					delete[] hh1[2];
					delete[] hh1[3];
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
					delete[] hh2;
					hh2 = 0;

					delete[] tt3[0];
					delete[] tt3[1];
					delete[] tt3[2];
					delete[] tt3[3];
					delete[] tt3[4];
					delete[] tt3[5];
					delete[] tt3[6];
					delete[] tt3[7];
					delete[] tt3[8];
					delete[] tt3[9];
					delete[] tt3[10];
					delete[] tt3[11];
					delete[] tt3[12];
					delete[] tt3[13];
					delete[] tt3[14];
					delete[] tt3[15];
					delete[] tt3;
					tt3 = 0;

					delete[] hh3[0];
					delete[] hh3[1];
					delete[] hh3[2];
					delete[] hh3[3];
					delete[] hh3[4];
					delete[] hh3[5];
					delete[] hh3[6];
					delete[] hh3[7];
					delete[] hh3[8];
					delete[] hh3[9];
					delete[] hh3[10];
					delete[] hh3[11];
					delete[] hh3[12];
					delete[] hh3[13];
					delete[] hh3[14];
					delete[] hh3[15];
					delete[] hh3;
					hh3 = 0;
				}
			}

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
					
					__dh[0] = new double[_nNode];
					__dh[1] = new double[_nNode];
					__dh[2] = new double[_nNode];
					__dh[3] = new double[_nNode];

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

					d3 = new double* [8];
					d3[0] = new double[nNode];
					d3[1] = new double[nNode];
					d3[2] = new double[nNode];
					d3[3] = new double[nNode];
					d3[4] = new double[nNode];
					d3[5] = new double[nNode];
					d3[6] = new double[nNode];
					d3[7] = new double[nNode];

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

					tt3 = new double* [16];
					tt3[0] = new double[uDim];
					tt3[1] = new double[vDim];
					tt3[2] = new double[uDim];
					tt3[3] = new double[vDim];
					tt3[4] = new double[uDim];
					tt3[5] = new double[vDim];
					tt3[6] = new double[uDim];
					tt3[7] = new double[vDim];
					tt3[8] = new double[uDim];
					tt3[9] = new double[vDim];
					tt3[10] = new double[uDim];
					tt3[11] = new double[vDim];
					tt3[12] = new double[uDim];
					tt3[13] = new double[vDim];
					tt3[14] = new double[uDim];
					tt3[15] = new double[vDim];
					hh3 = new double* [16];
					hh3[0] = new double[uDim];
					hh3[1] = new double[vDim];
					hh3[2] = new double[uDim];
					hh3[3] = new double[vDim];
					hh3[4] = new double[uDim];
					hh3[5] = new double[vDim];
					hh3[6] = new double[uDim];
					hh3[7] = new double[vDim];
					hh3[8] = new double[uDim];
					hh3[9] = new double[vDim];
					hh3[10] = new double[uDim];
					hh3[11] = new double[vDim];
					hh3[12] = new double[uDim];
					hh3[13] = new double[vDim];
					hh3[14] = new double[uDim];
					hh3[15] = new double[vDim];


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
		//double trEij = 0;
		//double treij = 0;
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
		double* __mix_phi = 0;
		double* __mix_z = 0;
		

		int _nNode;
	private:
		//double* __mat = 0;
		double* _F3 = 0;
		double _F3_phi[2];
		double _F3_z[2];

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
		double gammaijk[8];
		double Gammaijk2[8];
		double _ss[4];
		double _Sij[4];
		double eij[4];
		//double mij[4];
		//double fij[4];
		//double kij[4];
		double Eij[4];
		double gkij[8];
		//double hkij[8];
		double __Sij[4];
		double __hij[4];// , ___sij[4];
		//double v1, v2, s1, s2;
		
	public:
		double sc;
	public:
		///////shared memory//////
		double** M[2];
		int** dd;
		double* __mat = 0;
		double* __dsigma_11[3]{ 0,0,0 };
		double* __dsigma_12[3]{ 0,0,0 };
		double* __dsigma_22[3]{ 0,0,0 };
	
		//double* __dh[4]{ 0,0,0,0 };

		double* d0;
		double* d1[2];
		double* d2[4];
		double* d2_star[4];
		double* d3[8];
		const int star2[4]{ 1,-1,-1,1 };
		const int _star[4]{ 3,1,2,0 };
		double* B[4];
		double* tt0[2], * hh0[2], * tt1[4], * hh1[4], * tt2[8], * hh2[8], * tt3[16], * hh3[16];
		///////shared memory//////

		double* gradN[3]{ 0,0,0 };
		double* gradG = 0;
		const int ___ll[4]{ 0,3,6,9 };
	public:
		double x, y, z, Z, _Z, phi, eta, xi,nu,xi_1, xi_2, eta_1, eta_2,u,v,w;/// , chi;
		double dv, _dv, sigma_dv,sigma_trace;
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
		inline double& get_eij(const int& i, const int& j) {
			return eij[(i << 1) + j];
		}
		
		/*inline double& get_fij(const int& i, const int& j) {
			return fij[(i << 1) + j];
		}
		inline double& get_kij(const int& i, const int& j) {
			return kij[(i << 1) + j];
		}*/
		inline double& get_Eij(const int& i, const int& j) {
			return Eij[(i << 1) + j];
		}
		inline double& get__Sij(const int& i, const int& j) {
			return __Sij[((i << 1) + j)];
		}
		inline double& get__hij(const int& i, const int& j) {
			return __hij[((i << 1) + j)];
		}
		//inline double& get___sij(const int& i, const int& j) {
	//		return ___sij[((i << 1) + j)];
	//	}
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
		inline double& get__bij(const int& i, const int& j) {
			return bij[___ll[((i << 1) + j)] + 2];
		}
		inline double& get_Gammaijk(const int& i, const int& j, const int& k) {
			return Gammaijk[(((i << 1) + j) << 1) + k];
		}
		inline double& get_gammaijk(const int& i, const int& j, const int& k) {
			return gammaijk[(((i << 1) + j) << 1) + k];
		}
		inline double& get_Gammaijk2(const int& i, const int& j, const int& k) {
			return Gammaijk2[(((i << 1) + j) << 1) + k];
		}
		inline double& get_gkij(const int& k, const int& i, const int& j) {
			return gkij[(((i << 1) + j) << 1) + k];
		}
		//inline double& get_hkij(const int& k, const int& i, const int& j) {
		//	return hkij[(((i << 1) + j) << 1) + k];
		//}
	
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
		inline double& get_hh3(const int& i, const int& j, const int& k, const int& l, const int& s) {
			return _ref->hh3[(((((i << 1) + j) << 1) + k)<<1)+l][s];
		}
		inline double& get_tt3(const int& i, const int& j, const int& k, const int& l, const int& s) {
			return _ref->tt3[(((((i << 1) + j) << 1) + k)<<1)+l][s];
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
			double t = lo[j];//j:coordinate, k:which term in (x^3+x^2+x^1+1) 
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
		double __hh3(int s, int n, int m, int j, int k) {
			double t = lo[j];//j:coordinate, k:which term in (x^3+x^2+x^1+1)
			if (j != m && j != n&&j!=s)
			{
				return pow(t, (dim[j] - k - 1));
			}
			if ((j != m && j!=s && j == n) || (j == m && j != n && j!=s)|| (j == s && j != m && j != n))
			{
				if (k < dim[j] - 1)
				{
					return (dim[j] - k - 1) * pow(t, (dim[j] - k - 2));
				}
				else {
					return 0;
				}
			}
			if ((j == m && j == n && j!=s)|| (j == n && j == s && j != m)|| (j == m && j == s && j != n))
			{
				if (k < dim[j] - 2)
					return (dim[j] - k - 1) * (dim[j] - k - 2) * pow(t, (dim[j] - k - 3));
				else return  0;
			}
			if (j == m && j == n && j == s)
			{
				if (k < dim[j] - 3)
					return (dim[j] - k - 1) * (dim[j] - k - 2)* (dim[j] - k - 3) * pow(t, (dim[j] - k - 4));
				else return  0;

			}
			return 0;
		}
		double __tt3(int s, int n, int m, int j, int k) {
			double val = 0;
			for (int l = 0; l < dim[j]; l++)
			{
				val += get_hh3(s,n, m, j, l) * _ref->M[j][l][k];
			}
			return val;
		}
		double _shape(int k)
		{
			//shape function
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
			//first derivatives
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
			//second derivatives
			double D = 1.0;
			for (int j = 0; j < 2; j++)
			{
				D *= get_tt2(m, n, j, _ref->dd[k][j]);
			}
			return D;
		}
		double _H(int m, int n, int s, int k)
		{
			//third derivatives
			double H = 1.0;
			for (int j = 0; j < 2; j++)
			{
				H *= get_tt3(m, n, s,j, _ref->dd[k][j]);
			}
			return H;
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
			if (mode == "U") 
			{
				__grad = 0;
				__grad2 = 0;
				__grad_z = 0;
				__grad_phi = 0;
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
				//__dh[0] = 0;
				//__dh[1] = 0;
				//__dh[2] = 0;
				//__dh[3] = 0;
				__dsigma_11[0] = 0;
				__dsigma_11[1] = 0;
				__dsigma_11[2] = 0;
				__dsigma_12[0] = 0;
				__dsigma_12[1] = 0;
				__dsigma_12[2] = 0;
				__dsigma_22[0] = 0;
				__dsigma_22[1] = 0;
				__dsigma_22[2] = 0;
			
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
				tt3[0] = 0;
				tt3[1] = 0;
				tt3[2] = 0;
				tt3[3] = 0;
				tt3[4] = 0;
				tt3[5] = 0;
				tt3[6] = 0;
				tt3[7] = 0;
				tt3[8] = 0;
				tt3[9] = 0;
				tt3[10] = 0;
				tt3[11] = 0;
				tt3[12] = 0;
				tt3[13] = 0;
				tt3[14] = 0;
				tt3[15] = 0;
				hh3[0] = 0;
				hh3[1] = 0;
				hh3[2] = 0;
				hh3[3] = 0;
				hh3[4] = 0;
				hh3[5] = 0;
				hh3[6] = 0;
				hh3[7] = 0;
				hh3[8] = 0;
				hh3[9] = 0;
				hh3[10] = 0;
				hh3[11] = 0;
				hh3[12] = 0;
				hh3[13] = 0;
				hh3[14] = 0;
				hh3[15] = 0;
			}

		}
		~_memS() {
			del();
		}
		
		
	
		void update_optional()
		{
			double* pptr1 = &_ref->buf_xi[0];
			double* pptr2 = &_ref->d0[0];
			double val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			xi = val;


			pptr1 = &_ref->buf_eta[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			eta = val;


			pptr1 = &_ref->buf_nu[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			nu = val;


			pptr1 = &_ref->buf_xi[0];
			pptr2 = &_ref->d1[0][0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			xi_1 = val;

			pptr1 = &_ref->buf_xi[0];
			pptr2 = &_ref->d1[1][0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			xi_2 = val;

			pptr1 = &_ref->buf_eta[0];
			pptr2 = &_ref->d1[0][0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			eta_1 = val;

			pptr1 = &_ref->buf_eta[0];
			pptr2 = &_ref->d1[1][0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			eta_2 = val;

			double* xii = &_ref->buf_u[0];
			double* etai = &_ref->buf_v[0];
			//double* nui = &_ref->buf_nu[0];

			double* xij = &_ref->buf_u[0];
			double* etaj = &_ref->buf_v[0];
			//double* nuj = &_ref->buf_nu[0];

			double* g1i = _ref->d1[0];
			double* g2i = _ref->d1[1];
			double* g1j = _ref->d1[0];
			double* g2j = _ref->d1[1];
			memset(eij, 0, sizeof(double) * 4);
	
			memset(__dsigma_11[0], 0, sizeof(double) * _nNode);
			memset(__dsigma_11[1], 0, sizeof(double) * _nNode);

			memset(__dsigma_12[0], 0, sizeof(double) * _nNode);
			memset(__dsigma_12[1], 0, sizeof(double) * _nNode);

			memset(__dsigma_22[0], 0, sizeof(double) * _nNode);
			memset(__dsigma_22[1], 0, sizeof(double) * _nNode);


			memset(__dsigma_11[2], 0, sizeof(double) * _nNode);
			memset(__dsigma_12[2], 0, sizeof(double) * _nNode);
			memset(__dsigma_22[2], 0, sizeof(double) * _nNode);

			//memcpy(__dsigma_11[0], _ref->__dh[0], sizeof(double) * _nNode);
			//memcpy(__dsigma_12[0], _ref->__dh[1], sizeof(double) * _nNode);
			//memcpy(__dsigma_22[0], _ref->__dh[3], sizeof(double) * _nNode);

			//memcpy(__dsigma_11[2], _ref->__dh[0], sizeof(double) * _nNode);
			//memcpy(__dsigma_12[2], _ref->__dh[1], sizeof(double) * _nNode);
			//memcpy(__dsigma_22[2], _ref->__dh[3], sizeof(double) * _nNode);
			//memset(gkij, 0, sizeof(double) * 8);
			pptr1 = &_ref->buf_u[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			u = val;


			pptr1 = &_ref->buf_v[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			v = val;

			pptr1 = &_ref->buf_w[0];
			pptr2 = &_ref->d0[0];
			val = 0;
			for (int i = 0; i < _ref->_nNode; i++)
			{
				//val += _ref->d0[i] * _ref->buf_phi[i];
				val += (*pptr1) * (*pptr2);
				pptr1++;
				pptr2++;
			}
			w = val;
			for (int i = 0; i < _nNode; i++)
			{
				xij = &_ref->buf_u[0];
				etaj = &_ref->buf_v[0];
			//	nuj = &_ref->buf_nu[0];
				g1j = _ref->d1[0];
				g2j = _ref->d1[1];
				//eij[0] += _ref->__dh[0][i] * _ref->buf_nu[i];
				//eij[1] += _ref->__dh[1][i] * _ref->buf_nu[i];
				//eij[2] += _ref->__dh[2][i] * _ref->buf_nu[i];
				//eij[3] += _ref->__dh[3][i] * _ref->buf_nu[i];


				for (int j = 0; j < _ref->_nNode; j++)
				{
					__dsigma_11[0][i] += 2 * (*g1i * *g1j * *xij);
					__dsigma_11[1][i] += 2 * (*g1i * *g1j * *etaj);
					//__dsigma_11[2][i] += 2 * (*g1i * *g1j * *nuj);
					__dsigma_12[0][i] += (*g1i * *g2j * *xij) + (*g2i * *g1j * *xij);
					__dsigma_12[1][i] += (*g1i * *g2j * *etaj) + (*g2i * *g1j * *etaj);
					//__dsigma_12[2][i] += (*g1i * *g2j * *nuj) + (*g2i * *g1j * *nuj);
					__dsigma_22[0][i] += 2 * (*g2i * *g2j * *xij);
					__dsigma_22[1][i] += 2 * (*g2i * *g2j * *etaj);
					//__dsigma_22[2][i] += 2 * (*g2i * *g2j * *nuj);

					/*fij[0] += *g1i * *xii * *g1j * *xij;
					fij[0] += *g1i * *etai * *g1j * *etaj;
					fij[0] += *g1i * *nui * *g1j * *nuj;
					fij[1] += *g1i * *xii * *g2j * *xij;
					fij[1] += *g1i * *etai * *g2j * *etaj;
					fij[1] += *g1i * *nui * *g2j * *nuj;
					fij[3] += *g2i * *xii * *g2j * *xij;
					fij[3] += *g2i * *etai * *g2j * *etaj;
					fij[3] += *g2i * *nui * *g2j * *nuj;*/



					eij[0] += *g1i * *xii * *g1j * *xij;
					eij[0] += *g1i * *etai * *g1j * *etaj;
					//eij[0] += *g1i * *nui * *g1j * *nuj;
					eij[1] += *g1i * *xii * *g2j * *xij;
					eij[1] += *g1i * *etai * *g2j * *etaj;
					//eij[1] += *g1i * *nui * *g2j * *nuj;
					eij[3] += *g2i * *xii * *g2j * *xij;
					eij[3] += *g2i * *etai * *g2j * *etaj;
					//eij[3] += *g2i * *nui * *g2j * *nuj;

					g1j++;
					g2j++;
					xij++;
					etaj++;
					//nuj++;
				}
			
				g1i++;
				g2i++;
				xii++;
				etai++;
				//nui++;
			}
			eij[2] = eij[1];
		

	
			Eij[0] = eij[3] * sc;
			Eij[3] = eij[0] * sc;
			Eij[1] = -eij[1] * sc;
			Eij[2] = -eij[1] * sc;
		
			
		}
	
		void update2(string __mode) {
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
				ptr = _ref->hh3;
				for (auto const& s : ___ee)
				{
					for (auto const& n : ___ee)
					{
						for (auto const& m : ___ee)
						{
							for (auto const& j : ___ee)
							{
								for (int k = 0; k < dim[j]; k++) {
									*(*ptr + k) = __hh3(s,n, m, j, k);
									//hh2[(n * 2 + m) * 2 + j][k] = __hh2(n, m, j, k);
								}
								ptr++;
							}
						}
					}
				}
				ptr = _ref->tt3;
				for (auto const& s : ___ee)
				{
					for (auto const& n : ___ee)
					{
						for (auto const& m : ___ee)
						{
							for (auto const& j : ___ee)
							{
								for (int k = 0; k < dim[j]; k++) {
									*(*ptr + k) = __tt3(s, n, m, j, k);
									//tt2[(n * 2 + m) * 2 + j][k] = __tt2(n, m, j, k);
								}
								ptr++;
							}
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

				ptr = _ref->d3;
				for (auto const& i : ___ee) {
					for (auto const& ii : ___ee) {
						for (auto const& iii : ___ee) {
							ptr2 = *ptr;
							for (int j = 0; j < _nNode; j++) {
								//d2[i * 2 + ii][j] = _D(i, ii, j);
								//*ptr2 = _ref->buf_W[j] * _D(i, ii, j) / _ref->w;
								*ptr2 = _H(i, ii, iii, j);
								ptr2++;
							}
							ptr++;
						}
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


			if (__mode == "SIMPLE")return;



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

				sc = 1.0 / _det2(gij2);
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
				//	if (!_ref->initialized)
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

								val = 0;
								val += bij[eee + 0] * get_gi2(k, 0);
								val += bij[eee + 1] * get_gi2(k, 1);
								gammaijk[ccc] = val;
								ccc++;
							}
							eee += 3;
						}
					}
				}
			}
			else {

				//	if (!_ref->initialized)
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
							bij[e + 2] = fz;
						}
					}
					int ccc = 0;
					static const int fff[8]{ 0,2,1,3,4,6,5,7 };

					int eee = 0;
					for (auto const& i : ___ee) {
						for (auto const& j : ___ee) {
							for (auto const& k : ___ee) {
								double val = 0;

								val += bij[eee + 0] * get_Gi(k, 0);
								val += bij[eee + 1] * get_Gi(k, 1);
								//val += bij[eee + 2] * get_Gi(k, 2);
								Gammaijk[ccc] = val;
								val = 0;
								val += bij[eee + 0] * get_Gi2(k, 0);
								val += bij[eee + 1] * get_Gi2(k, 1);
								val += bij[eee + 2] * get_Gi2(k, 2);
								Gammaijk2[ccc] = val;
								val = 0;
								val += bij[eee + 0] * get_gi(k, 0);
								val += bij[eee + 1] * get_gi(k, 1);
								gammaijk[ccc] = val;
								ccc++;
							}
							eee += 3;
						}
					}
				}
			}

			if (__mode == "STANDARD")return;





			/*if (mode == "U")*/ {

				double* ptr4;




				double val = 0;
				double val2 = 0;
				double val3 = 0;
				double* pptr = 0;
				double* pptr1 = &_ref->buf_z[0];
				double* pptr2 = &_ref->d0[0];

				static const int sss[4]{ 0,1,2,3 };

				//if (!_ref->initialized || RAM == SAVE)
				{



					for (int i = 0; i < _nNode; i++)
					{
						_ref->__dh[0][i] = _ref->d2[0][i] - Gammaijk[0] * _ref->d1[0][i] - Gammaijk[1] * _ref->d1[1][i];
						_ref->__dh[1][i] = _ref->d2[1][i] - Gammaijk[2] * _ref->d1[0][i] - Gammaijk[3] * _ref->d1[1][i];
						_ref->__dh[2][i] = _ref->d2[2][i] - Gammaijk[4] * _ref->d1[0][i] - Gammaijk[5] * _ref->d1[1][i];
						_ref->__dh[3][i] = _ref->d2[3][i] - Gammaijk[6] * _ref->d1[0][i] - Gammaijk[7] * _ref->d1[1][i];

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


				{
					memset(__Sij, 0, sizeof(double) * 4);
					memset(__hij, 0, sizeof(double) * 4);
					double __h11 = 0, __h12 = 0, __h22 = 0;
					double* pptr1 = &_ref->buf_z[0];
					double* pptr2 = &_ref->buf_phi[0];
					double* pdh0 = _ref->__dh[0];
					double* pdh1 = _ref->__dh[1];
					double* pdh3 = _ref->__dh[3];
					double* qdh0 = _ref->__dh[0];
					double* qdh1 = _ref->__dh[1];
					double* qdh3 = _ref->__dh[3];
					for (int i = 0; i < _nNode; i++)
					{
						__Sij[0] += *pdh0 * *pptr1;
						__Sij[1] += *pdh1 * *pptr1;
						__Sij[3] += *pdh3 * *pptr1;
						__h11 += *qdh0 * *pptr2;
						__h12 += *qdh1 * *pptr2;
						__h22 += *qdh3 * *pptr2;
						pptr1++;
						pptr2++;
						pdh0++;
						pdh1++;
						pdh3++;
						qdh0++;
						qdh1++;
						qdh3++;
					}
					__Sij[2] = __Sij[1];
					__hij[0] = __h11;
					__hij[1] = __h12;
					__hij[2] = __h12;
					__hij[3] = __h22;
				}


			}
			_ref->initialized = true;
			if (mode == "SLOPE") {

				/*for (int i = 0; i < _nNode; i++) {
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
					*/
			}

			set_Sij(1, 0, 1);
			/*for (int i = 0; i < _nNode; i++) {
				for (auto const& l : ___ee) {
					_F3[l * _nNode + i] = __F3(l, i);
				}
			}
			for (auto const& l : ___ee) {
				_F3_phi[l] = __F3_phi(l);
			}
			for (auto const& l : ___ee) {
				_F3_z[l] = __F3_z(l);
			}*/


		}

		void del() {
			if (mode == "U") {
				if (__grad_z != 0)delete[] __grad_z;
				if (__grad_phi != 0)delete[] __grad_phi;
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
					delete[] __dsigma_11[0];
					delete[] __dsigma_11[1];
					delete[] __dsigma_11[2];
					delete[] __dsigma_12[0];
					delete[] __dsigma_12[1];
					delete[] __dsigma_12[2];
					delete[] __dsigma_22[0];
					delete[] __dsigma_22[1];
					delete[] __dsigma_22[2];
					//delete[] __dh[0];
					//delete[] __dh[1];
					//delete[] __dh[2];
					//delete[] __dh[3];
					__mat = 0;
					__dsigma_11[0] = 0;
					__dsigma_11[1] = 0;
					__dsigma_11[2] = 0;
					__dsigma_12[0] = 0;
					__dsigma_12[1] = 0;
					__dsigma_12[2] = 0;
					__dsigma_22[0] = 0;
					__dsigma_22[1] = 0;
					__dsigma_22[2] = 0;

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

					delete[] tt3[0];
					delete[] tt3[1];
					delete[] tt3[2];
					delete[] tt3[3];
					delete[] tt3[4];
					delete[] tt3[5];
					delete[] tt3[6];
					delete[] tt3[7];
					delete[] tt3[8];
					delete[] tt3[9];
					delete[] tt3[10];
					delete[] tt3[11];
					delete[] tt3[12];
					delete[] tt3[13];
					delete[] tt3[14];
					delete[] tt3[15];
					tt3[0] = 0;
					tt3[1] = 0;
					tt3[2] = 0;
					tt3[3] = 0;
					tt3[4] = 0;
					tt3[5] = 0;
					tt3[6] = 0;
					tt3[7] = 0;
					tt3[8] = 0;
					tt3[9] = 0;
					tt3[10] = 0;
					tt3[11] = 0;
					tt3[12] = 0;
					tt3[13] = 0;
					tt3[14] = 0;
					tt3[15] = 0;

					delete[] hh3[0];
					delete[] hh3[1];
					delete[] hh3[2];
					delete[] hh3[3];
					delete[] hh3[4];
					delete[] hh3[5];
					delete[] hh3[6];
					delete[] hh3[7];
					delete[] hh3[8];
					delete[] hh3[9];
					delete[] hh3[10];
					delete[] hh3[11];
					delete[] hh3[12];
					delete[] hh3[13];
					delete[] hh3[14];
					delete[] hh3[15];
					hh3[0] = 0;
					hh3[1] = 0;
					hh3[2] = 0;
					hh3[3] = 0;
					hh3[4] = 0;
					hh3[5] = 0;
					hh3[6] = 0;
					hh3[7] = 0;
					hh3[8] = 0;
					hh3[9] = 0;
					hh3[10] = 0;
					hh3[11] = 0;
					hh3[12] = 0;
					hh3[13] = 0;
					hh3[14] = 0;
					hh3[15] = 0;
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
					__dsigma_11[0] = new double[nNode];
					__dsigma_11[1] = new double[nNode];
					__dsigma_11[2] = new double[nNode];
					__dsigma_12[0] = new double[nNode];
					__dsigma_12[1] = new double[nNode];
					__dsigma_12[2] = new double[nNode];
					__dsigma_22[0] = new double[nNode];
					__dsigma_22[1] = new double[nNode];
					__dsigma_22[2] = new double[nNode];


					__grad = new double[nNode];
					__grad2 = new double[nNode];
					__grad_z = new double[nNode];
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
					//__dh[0] = new double[nNode];
					//__dh[1] = new double[nNode];
					//__dh[2] = new double[nNode];
					//__dh[3] = new double[nNode];
					d0 = new double[nNode];

					d1[0] = new double[nNode];
					d1[1] = new double[nNode];

					d2[0] = new double[nNode];
					d2[1] = new double[nNode];
					d2[2] = new double[nNode];
					d2[3] = new double[nNode];

					d3[0] = new double[nNode];
					d3[1] = new double[nNode];
					d3[2] = new double[nNode];
					d3[3] = new double[nNode];
					d3[4] = new double[nNode];
					d3[5] = new double[nNode];
					d3[6] = new double[nNode];
					d3[7] = new double[nNode];

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

					tt3[0] = new double[uDim];
					tt3[1] = new double[vDim];
					tt3[2] = new double[uDim];
					tt3[3] = new double[vDim];
					tt3[4] = new double[uDim];
					tt3[5] = new double[vDim];
					tt3[6] = new double[uDim];
					tt3[7] = new double[vDim];
					tt3[8] = new double[uDim];
					tt3[9] = new double[vDim];
					tt3[10] = new double[uDim];
					tt3[11] = new double[vDim];
					tt3[12] = new double[uDim];
					tt3[13] = new double[vDim];
					tt3[14] = new double[uDim];
					tt3[15] = new double[vDim];
					hh3[0] = new double[uDim];
					hh3[1] = new double[vDim];
					hh3[2] = new double[uDim];
					hh3[3] = new double[vDim];
					hh3[4] = new double[uDim];
					hh3[5] = new double[vDim];
					hh3[6] = new double[uDim];
					hh3[7] = new double[vDim];
					hh3[8] = new double[uDim];
					hh3[9] = new double[vDim];
					hh3[10] = new double[uDim];
					hh3[11] = new double[vDim];
					hh3[12] = new double[uDim];
					hh3[13] = new double[vDim];
					hh3[14] = new double[uDim];
					hh3[15] = new double[vDim];

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
		inline double ____SLOPE_phi(int l) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_phi[i] * _ref->get__Gij(k, l);
				}
			}
			return val;
		}
		inline double ____SLOPE_z(int l) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_z[i] * _ref->get__Gij(k, l);
				}
			}
			return val;

		}
		inline double ___SLOPE_phi(int l) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_phi[i] * this->get_Gij(k, l);
				}
			}
			return val;
		}
		inline double ___SLOPE_z(int l) {
			double val = 0;
			for (int k = 0; k < 2; k++)
			{
				for (int i = 0; i < _nNode; i++)
				{
					val += _ref->d1[k][i] * _ref->buf_z[i] * this->get_Gij(k, l);
				}
			}
			return val;

		}
		
		
		inline double ___SLOPE_phi3(double v1, double v2) {
		
			double t1 = v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
			double t2 = -v1 * _ref->get__gij(0, 0) - v2 * _ref->get__gij(1, 0);
			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length;
			t2 /= length;
			double val = 0;
			for (int i = 0; i < _nNode; i++)
			{
				val += _ref->d1[0][i] * _ref->buf_phi[i] * t1;
				val += _ref->d1[1][i] * _ref->buf_phi[i] * t2;
			}
			return val;
		}

		void ___SLOPE_phi3_phi(double* ptr, double v1, double v2) {
			
			double t1 = v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
			double t2 = -v1 * _ref->get__gij(0, 0) - v2 * _ref->get__gij(1, 0);
			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length;
			t2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				val = _ref->d1[0][i] * t1;
				val += _ref->d1[1][i] * t2;
				*ptr1 = val;
				ptr1++;
			}
		}

	
		
		inline double ___SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {

			return sc * (____SLOPE_phi(0) * dcdtstar0 + ____SLOPE_phi(1) * dcdtstar1);// val;
		}
		inline double ___SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {

			return sc * (____SLOPE_z(0) * dcdtstar0 + ____SLOPE_z(1) * dcdtstar1);// val;
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
		inline double _SLOPE(int l, int i) {
			return _ref->d1[0][i] * _ref->get__Gij(0, l) + _ref->d1[1][i] * _ref->get__Gij(1, l);

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
		inline void _SLOPE_phi_x(double* ptr, double dcdt1, double dcdt2, double sc) {
			double* ptr1 = ptr;
			double phi1 = 0, phi2 = 0;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				phi1 += _ref->d1[0][s] * _ref->buf_phi[s];
				phi2 += _ref->d1[1][s] * _ref->buf_phi[s];
			}

			for (int i = 0; i < _nNode; i++)
			{
				double _g11 = 2*_ref->d1[0][i] * _ref->get__gi(0, 0);
				double _g12 = _ref->d1[0][i] * _ref->get__gi(1, 0)+ _ref->d1[1][i] * _ref->get__gi(0, 0);
				double _g22 = 2 * _ref->d1[1][i] * _ref->get__gi(1, 0);
				double _g21 = _g12;
				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;
				*ptr1 = ((phi1 * _G11+ phi2 * _G21) * dcdt1 +
					(phi1 *_G12 + phi2 * _G22) * dcdt2) * sc;
				ptr1++;

			}

		}
		inline void _SLOPE_phi_y(double* ptr, double dcdt1, double dcdt2, double sc) {
			double* ptr1 = ptr;
			double phi1 = 0, phi2 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				phi1 += _ref->d1[0][s] * _ref->buf_phi[s];
				phi2 += _ref->d1[1][s] * _ref->buf_phi[s];
			}

			for (int i = 0; i < _nNode; i++)
			{
				double _g11 = 2 * _ref->d1[0][i] * _ref->get__gi(0, 1);
				double _g12 = _ref->d1[0][i] * _ref->get__gi(1, 1) + _ref->d1[1][i] * _ref->get__gi(0, 1);
				double _g22 = 2 * _ref->d1[1][i] * _ref->get__gi(1, 1);
				double _g21 = _g12;
				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;
				*ptr1 = ((phi1 * _G11 + phi2 * _G21) * dcdt1 +
					(phi1 * _G12 + phi2 * _G22) * dcdt2) * sc;
				ptr1++;

			}

		}
		

		inline double _SLOPE_symm_sigma(double w1, double w2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_xi[i];
				xi_v += _ref->d1[1][i] * _ref->buf_xi[i];
				eta_u += _ref->d1[0][i] * _ref->buf_eta[i];
				eta_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
		
			double xi = (xi_u * w1 + xi_v * w2);
			double eta = (eta_u * w1 + eta_v * w2);
			double xi0_u = _ref->xi_1;
			double xi0_v = _ref->xi_2;
			double eta0_u = _ref->eta_1;
			double eta0_v = _ref->eta_2;

			double xi0 = (xi0_u* w1+ xi0_v * w2);
			double eta0 = (eta0_u * w1 + eta0_v *w2);


			double val = xi * eta0 - eta * xi0;
			return val;
		}
		void _SLOPE_symm_sigma_xi(double *ptr,double w1, double w2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double* ptr1 = ptr;
			
			double xi0_u = _ref->xi_1;
			double xi0_v = _ref->xi_2;
			double eta0_u = _ref->eta_1;
			double eta0_v = _ref->eta_2;


			double xi0 = (xi0_u * w1 + xi0_v * w2);
			double eta0 = (eta0_u * w1 + eta0_v * w2);

			for (int s = 0; s< _nNode; s++)
			{
				double _xi_u = _ref->d1[0][s];
				double _xi_v = _ref->d1[1][s];
			

				double _xi = (_xi_u * w1 + _xi_v * w2);

				double val = _xi *eta0;
				*ptr1 = val;
				ptr1 ++ ;
			}
			
		}
		void _SLOPE_symm_sigma_eta(double* ptr,double w1, double w2) {

			double* ptr1 = ptr;
			
			double xi0_u = _ref->xi_1;
			double xi0_v = _ref->xi_2;
			double eta0_u = _ref->eta_1;
			double eta0_v = _ref->eta_2;


			double xi0 = (xi0_u * w1 + xi0_v * w2);
			double eta0 = (eta0_u * w1 + eta0_v * w2);

			for (int s = 0; s < _nNode; s++)
			{
				double _eta_u = _ref->d1[0][s];
				double _eta_v = _ref->d1[1][s];
			
				double _eta = _eta_u * w1 + _eta_v * w2;


				double val = -_eta * xi0;
				*ptr1 = val;
				ptr1++;
			}
		}


		inline double _SLOPE_symm_sigma2(double w1, double w2,double s1,double s2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_xi[i];
				xi_v += _ref->d1[1][i] * _ref->buf_xi[i];
				eta_u += _ref->d1[0][i] * _ref->buf_eta[i];
				eta_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
			double duu = xi_u * _ref->get__gi(0, 0) + eta_u * _ref->get__gi(0, 1);
			double duv = xi_u * _ref->get__gi(1, 1) + eta_u * _ref->get__gi(1, 1);
			double dvu = xi_v * _ref->get__gi(0, 0) + eta_v * _ref->get__gi(0, 1);
			double dvv = xi_v * _ref->get__gi(1, 1) + eta_v * _ref->get__gi(1, 1);


			double val = duu * w1 * s1 + duv * w1 * s2 + dvu * w2 * s1 + dvv * s2 * s2;
			return val;
		}
		void _SLOPE_symm_sigma2_u(double* ptr, double w1, double w2, double s1, double s2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_xi[i];
				xi_v += _ref->d1[1][i] * _ref->buf_xi[i];
				eta_u += _ref->d1[0][i] * _ref->buf_eta[i];
				eta_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
			double* ptr1 = ptr;


			for (int s = 0; s < _nNode; s++)
			{
			
				double duu = xi_u * _ref->d1[0][s];
				double duv = xi_u * _ref->d1[1][s];
				double dvu = xi_v * _ref->d1[0][s];
				double dvv = xi_v * _ref->d1[1][s];


				double val = duu * w1 * s1 + duv * w1 * s2 + dvu * w2 * s1 + dvv * s2 * s2;
				*ptr1 = val;
				ptr1++;
			}

		}
		void _SLOPE_symm_sigma2_v(double* ptr, double w1, double w2, double s1, double s2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_xi[i];
				xi_v += _ref->d1[1][i] * _ref->buf_xi[i];
				eta_u += _ref->d1[0][i] * _ref->buf_eta[i];
				eta_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
			double* ptr1 = ptr;


			for (int s = 0; s < _nNode; s++)
			{
			
				double duu = eta_u * _ref->d1[0][s];
				double duv = eta_u * _ref->d1[1][s];
				double dvu = eta_v * _ref->d1[0][s];
				double dvv = eta_v * _ref->d1[1][s];


				double val = duu * w1 * s1 + duv * w1 * s2 + dvu * w2 * s1 + dvv * s2 * s2;
				*ptr1 = val;
				ptr1++;
			}

		}
		void _SLOPE_symm_sigma2_xi(double* ptr, double w1, double w2,double s1,double s2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;

			double* ptr1 = ptr;

		
			for (int s = 0; s < _nNode; s++)
			{
				double xi_u = 0, xi_v = 0;
				double eta_u = 0, eta_v = 0;
				
					xi_u = _ref->d1[0][s];
					xi_v = _ref->d1[1][s] ;
					eta_u = 0;// _ref->d1[0][i] * _ref->buf_eta[i];
					eta_v = 0;// _ref->d1[1][i] * _ref->buf_eta[i];
					double duu = xi_u * _ref->get__gi(0, 0) + eta_u * _ref->get__gi(0, 1);
					double duv = xi_u * _ref->get__gi(1, 1) + eta_u * _ref->get__gi(1, 1);
					double dvu = xi_v * _ref->get__gi(0, 0) + eta_v * _ref->get__gi(0, 1);
					double dvv = xi_v * _ref->get__gi(1, 1) + eta_v * _ref->get__gi(1, 1);


				double val = duu * w1 * s1 + duv * w1 * s2 + dvu * w2 * s1 + dvv * s2 * s2;
				*ptr1 = val;
				ptr1++;
			}

		}
		void _SLOPE_symm_sigma2_eta(double* ptr, double w1, double w2,double s1,double s2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			double* ptr1 = ptr;

		

			for (int s = 0; s < _nNode; s++)
			{
				double xi_u = 0, xi_v = 0;
				double eta_u = 0, eta_v = 0;

				xi_u = 0;// _ref->d1[0][s];
				xi_v = 0;// _ref->d1[1][s];
				eta_u =  _ref->d1[0][s];
				eta_v = _ref->d1[1][s];

				double duu = xi_u * _ref->get__gi(0, 0) + eta_u * _ref->get__gi(0, 1);
				double duv = xi_u * _ref->get__gi(1, 1) + eta_u * _ref->get__gi(1, 1);
				double dvu = xi_v * _ref->get__gi(0, 0) + eta_v * _ref->get__gi(0, 1);
				double dvv = xi_v * _ref->get__gi(1, 1) + eta_v * _ref->get__gi(1, 1);


				double val = duu * w1 * s1 + duv * w1 * s2 + dvu * w2 * s1 + dvv * s2 * s2;
				*ptr1 = val;
				ptr1++;
			}
		}
		inline double cont2(_memS* other,double w1,double w2,double  w1_2,double  w2_2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(w1_2 * w1_2 * other->_ref->og11 + 2 * w1_2 * w2_2 * other->_ref->og12 + w2_2 * w2_2 * other->_ref->og22);
			w1_2 /= length;
			w2_2 /= length;
			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;

			xi_u = _ref->get__gi(0, 0);
			xi_v = _ref->get__gi(1, 0);
			eta_u = _ref->get__gi(0, 1);
			eta_v = _ref->get__gi(1, 1);

			double xi_u_other = 0, xi_v_other = 0;
			double eta_u_other = 0, eta_v_other = 0;
			xi_u_other = other->_ref->get__gi(0, 0);
			xi_v_other = other->_ref->get__gi(1, 0);
			eta_u_other = other->_ref->get__gi(0, 1);
			eta_v_other = other->_ref->get__gi(1, 1);

			double val = (xi_u * w1 + xi_v * w2) * (eta_u_other * w1_2 + eta_v_other * w2_2)
						- (eta_u * w1 + eta_v * w2) * (xi_u_other * w1_2 + xi_v_other * w2_2);

			return val;
		}
		void cont2_u(_memS* other, double* ptr,double *ptr2, double w1, double w2, double  w1_2, double  w2_2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(w1_2 * w1_2 * other->_ref->og11 + 2 * w1_2 * w2_2 * other->_ref->og12 + w2_2 * w2_2 * other->_ref->og22);
			w1_2 /= length;
			w2_2 /= length;

			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;
	

			xi_u = _ref->get__gi(0, 0);
			xi_v = _ref->get__gi(1, 0);
			eta_u = _ref->get__gi(0, 1);
			eta_v = _ref->get__gi(1, 1);

			double xi_u_other = 0, xi_v_other = 0;
			double eta_u_other = 0, eta_v_other = 0;
			xi_u_other = other->_ref->get__gi(0, 0);
			xi_v_other = other->_ref->get__gi(1, 0);
			eta_u_other = other->_ref->get__gi(0, 1);
			eta_v_other = other->_ref->get__gi(1, 1);
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi_u = _ref->d1[0][s];
				double _xi_v = _ref->d1[1][s];
				//double _eta_u = _ref->d1[0][s];
				//double _eta_v = _ref->d1[1][s];

				double val = (_xi_u * w1 + _xi_v * w2) * (eta_u_other * w1_2 + eta_v_other * w2_2);
					// -(eta_u * w1 + eta_v * w2) * (xi_u_other * w1 + xi_v_other * w2);
				*ptr1 = val;
				ptr1++;
			}
			ptr1 = ptr2;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi_u_other = other->_ref->d1[0][s];
				double _xi_v_other = other->_ref->d1[1][s];
				//double _eta_u = other->_ref->d1[0][s];
				//double _eta_v = other->_ref->d1[1][s];

				//double val = _xi_u * w1 + xi_v * w2) * (eta_u_other * w1 + eta_v_other * w2);
				double val= -(eta_u * w1 + eta_v * w2) * (_xi_u_other * w1_2 + _xi_v_other * w2_2);
				*ptr1 = val;
				ptr1++;
			}
		}
		void cont2_v(_memS* other, double* ptr,double *ptr2, double w1, double w2, double  w1_2, double  w2_2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(w1_2 * w1_2 * other->_ref->og11 + 2 * w1_2 * w2_2 * other->_ref->og12 + w2_2 * w2_2 * other->_ref->og22);
			w1_2 /= length;
			w2_2 /= length;

			double xi_u = 0, xi_v = 0;
			double eta_u = 0, eta_v = 0;
			
			xi_u = _ref->get__gi(0, 0);
			xi_v = _ref->get__gi(1, 0);
		    eta_u = _ref->get__gi(0, 1);
			eta_v = _ref->get__gi(1, 1);

			double xi_u_other = 0, xi_v_other = 0;
			double eta_u_other = 0, eta_v_other = 0;
			xi_u_other = other->_ref->get__gi(0, 0);
			xi_v_other = other->_ref->get__gi(1, 0);
			eta_u_other = other->_ref->get__gi(0, 1);
			eta_v_other = other->_ref->get__gi(1, 1);
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _xi_u = _ref->d1[0][s];
				//double _xi_v = _ref->d1[1][s];
				double _eta_u = _ref->d1[0][s];
				double _eta_v = _ref->d1[1][s];

				//double val = (_xi_u * w1 + _xi_v * w2) * (eta_u_other * w1 + eta_v_other * w2);
				double val= -(_eta_u * w1 + _eta_v * w2) * (xi_u_other * w1_2 + xi_v_other * w2_2);
				*ptr1 = val;
				ptr1++;
			}
			ptr1 = ptr2;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _xi_u_other = other->_ref->d1[0][s];
				//double _xi_v_other = other->_ref->d1[1][s];
				double _eta_u_other = other->_ref->d1[0][s];
				double _eta_v_other = other->_ref->d1[1][s];

				double val = (xi_u * w1 + xi_v * w2) * (_eta_u_other * w1_2 + _eta_v_other * w2_2);
				//double val = -(eta_u * w1 + eta_v * w2) * (_xi_u_other * w1 + _xi_v_other * w2);
				*ptr1 = val;
				ptr1++;
			}
		}
		inline double _SLOPE_symm2(double w1, double w2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			
			double x_u = 0, x_v = 0;
			double y_u = 0, y_v = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				x_u += _ref->d1[0][s] * _ref->buf_u[s];
				x_v += _ref->d1[1][s] * _ref->buf_u[s];
				y_u += _ref->d1[0][s] * _ref->buf_v[s];
				y_v += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double x0 = _ref->_ogi[0] * w1 + _ref->_ogi[3] * w2;
			double y0 = _ref->_ogi[1] * w1 + _ref->_ogi[4] * w2;

			double X = x_u * w1 + x_v * w2;
			double Y = y_u * w1 + y_v * w2;


			double val = X * y0 - Y * x0;
			return val;
		}
		void _SLOPE_symm2_u(double* ptr, double w1, double w2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			

			double* ptr1 = ptr;

	

			double x0 = _ref->_ogi[0] * w1 + _ref->_ogi[3] * w2;
			double y0 = _ref->_ogi[1] * w1 + _ref->_ogi[4] * w2;
			for (int s = 0; s < _nNode; s++)
			{
				double _xi_u = 0, _xi_v = 0;
				double _eta_u = 0, _eta_v = 0;

				_xi_u = _ref->d1[0][s];
				_xi_v = _ref->d1[1][s];
				_eta_u = 0;// _ref->d1[0][i] * _ref->buf_eta[i];
				_eta_v = 0;// _ref->d1[1][i] * _ref->buf_eta[i];


				double X = _xi_u * w1 + _xi_v * w2;
				double Y = _eta_u * w1 + _eta_v * w2;


				double val = X * y0 - Y * x0;
				ptr1++;
			}

		}
		void _SLOPE_symm2_v(double* ptr, double w1, double w2) {
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length; 
			
			double* ptr1 = ptr;
		
	

			double x0 = _ref->_ogi[0] * w1 + _ref->_ogi[3] * w2;
			double y0 = _ref->_ogi[1] * w1 + _ref->_ogi[4] * w2;
			for (int s = 0; s < _nNode; s++)
			{
				double _xi_u = 0, _xi_v = 0;
				double _eta_u = 0, _eta_v = 0;

				_xi_u = 0;
				_xi_v = 0;
				_eta_u = _ref->d1[0][s];
				_eta_v = _ref->d1[1][s];

				double X = _xi_u * w1 + _xi_v * w2;
				double Y = _eta_u * w1 + _eta_v * w2;


				double val = X * y0 - Y * x0;
				*ptr1 = val;
				ptr1++;
			}
		}

		inline double _SLOPE_xi(double dcdtstar0, double dcdtstar1) {
		

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
		inline double _SLOPE_eta(double dcdtstar0, double dcdtstar1) {

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
		inline double _SLOPE_u(double dcdtstar0, double dcdtstar1) {

			double xi_u = 0, xi_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_u[i];
				xi_v += _ref->d1[1][i] * _ref->buf_u[i];
			}
			double f_U = xi_u * _ref->get__Gij(0, 0) + xi_v * _ref->get__Gij(1, 0);
			double f_V = xi_u * _ref->get__Gij(0, 1) + xi_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
		inline double _SLOPE_v(double dcdtstar0, double dcdtstar1) {

			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				eta_u += _ref->d1[0][i] * _ref->buf_v[i];
				eta_v += _ref->d1[1][i] * _ref->buf_v[i];
			}
			double f_U = eta_u * _ref->get__Gij(0, 0) + eta_v * _ref->get__Gij(1, 0);
			double f_V = eta_u * _ref->get__Gij(0, 1) + eta_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
		inline double _SLOPE_nu(double dcdtstar0, double dcdtstar1) {

			double nu_u = 0, nu_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				nu_u += _ref->d1[0][i] * _ref->buf_nu[i];
				nu_v += _ref->d1[1][i] * _ref->buf_nu[i];
			}
			double f_U = nu_u * _ref->get__Gij(0, 0) + nu_v * _ref->get__Gij(1, 0);
			double f_V = nu_u * _ref->get__Gij(0, 1) + nu_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
		double _SLOPE3phiX() {
			
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			double val=f_u* _ref->get__Gi(0, 0)+ f_v * _ref->get__Gi(1, 0);
			return val;

		}
		void _SLOPE3phiX_phi(double *ptr) {

			double f_u = 0, f_v = 0;
			double *ptr1=ptr;
			for (int i = 0; i < _nNode; i++)
			{
				f_u = _ref->d1[0][i] ;
				f_v = _ref->d1[1][i] ;
			
			     double val = f_u * _ref->get__Gi(0, 0) + f_v * _ref->get__Gi(1, 0);
			
				*ptr1=val;
				ptr1++;
			}

		}
		double _SLOPE3phiY() {

			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			double val = f_u * _ref->get__Gi(0, 1) + f_v * _ref->get__Gi(1, 1);
			return val;

		}
		void _SLOPE3phiY_phi(double* ptr) {

			double f_u = 0, f_v = 0;
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				f_u = _ref->d1[0][i];
				f_v = _ref->d1[1][i];

				double val = f_u * _ref->get__Gi(0, 1) + f_v * _ref->get__Gi(1, 1);

				*ptr1 = val;
				ptr1++;
			}

		}


		double _SLOPE3zX() {

			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			double val = f_u * _ref->get__Gi(0, 0) + f_v * _ref->get__Gi(1, 0);
			return val;

		}
		void _SLOPE3zX_z(double* ptr) {

			double f_u = 0, f_v = 0;
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				f_u = _ref->d1[0][i];
				f_v = _ref->d1[1][i];

				double val = f_u * _ref->get__Gi(0, 0) + f_v * _ref->get__Gi(1, 0);

				*ptr1 = val;
				ptr1++;
			}

		}
		double _SLOPE3zY() {

			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			double val = f_u * _ref->get__Gi(0, 1) + f_v * _ref->get__Gi(1, 1);
			return val;

		}
		void _SLOPE3zY_z(double* ptr) {

			double f_u = 0, f_v = 0;
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				f_u = _ref->d1[0][i];
				f_v = _ref->d1[1][i];

				double val = f_u * _ref->get__Gi(0, 1) + f_v * _ref->get__Gi(1, 1);

				*ptr1 = val;
				ptr1++;
			}

		}

		double _SLOPE2_phi(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_z(double v1, double v2) {
			
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_xi(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_xi[i];
				f_v += _ref->d1[1][i] * _ref->buf_xi[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_eta(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_eta[i];
				f_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_nu(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_nu[i];
				f_v += _ref->d1[1][i] * _ref->buf_nu[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_u(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_u[i];
				f_v += _ref->d1[1][i] * _ref->buf_u[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_v(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_v[i];
				f_v += _ref->d1[1][i] * _ref->buf_v[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		double _SLOPE2_w(double v1, double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_w[i];
				f_v += _ref->d1[1][i] * _ref->buf_w[i];
			}
			//double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			//double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_u * v1 + f_v * v2);// val;

		}
		void _SLOPE2(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);

			v1 /= length;
			v2 /= length;
			double f_u = 0, f_v = 0;
			double val = 0;
			double* ptr1 = ptr;
			for (int i = 0; i < _nNode; i++)
			{
				double _f_u = _ref->d1[0][i];
				double _f_v = _ref->d1[1][i];
				val = (_f_u * v1 + _f_v * v2);
				*ptr1 = val;
				ptr1++;
			}
		}
		void _SLOPE2_phi_x(double* ptr, double v1, double v2)
		{

			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			//v1 /= length;
			//v2 /= length;
			double f_u = 0, f_v = 0;
			double val = 0;
			double* ptr1 = ptr;

			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			for (int i = 0; i < _nNode; i++)
			{
				//double _g11 = 2 * _ref->d1[0][i] * _ref->get__gi(0, 0);
				//double _g12 = _ref->d1[0][i] * _ref->get__gi(1, 0) + _ref->d1[1][i] * _ref->get__gi(0, 0);
				//double _g22 = 2 * _ref->d1[1][i] * _ref->get__gi(1, 0);
				//double _g21 = _g12;
				double _g11 = __dsigma_11[0][i];
				double _g12 = __dsigma_12[0][i];
				double _g22 = __dsigma_22[0][i];
				double _g21 = _g12;

				
				double _length = 0.5 / (length) * (v1 * v1 * _g11 + 2 * v1 * v2 * _g12 + v2 * v2 * _g22);
				double _v1 = -v1 / length/length * _length;
				double _v2 = -v2 / length / length * _length;

				val = f_u * _v1 + f_v * _v2;
				*ptr1 = val;

				ptr1++;
			}
		}
		void _SLOPE2_phi_y(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			//v1 /= length;
			//v2 /= length;
			double f_u = 0, f_v = 0;
			double val = 0;
			double* ptr1 = ptr;

			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			for (int i = 0; i < _nNode; i++)
			{
				//double _g11 = 2 * _ref->d1[0][i] * _ref->get__gi(0, 1);
				//double _g12 = _ref->d1[0][i] * _ref->get__gi(1, 1) + _ref->d1[1][i] * _ref->get__gi(0, 1);
				//double _g22 = 2 * _ref->d1[1][i] * _ref->get__gi(1, 1);
				//double _g21 = _g12;
				double _g11 = __dsigma_11[1][i];
				double _g12 = __dsigma_12[1][i];
				double _g22 = __dsigma_22[1][i];
				double _g21 = _g12;

				double _length = 0.5 / (length) * (v1 * v1 * _g11 + 2 * v1 * v2 * _g12 + v2 * v2 * _g22);
				double _v1 = -v1 / length / length * _length;
				double _v2 = -v2 / length / length * _length;
				val = f_u * _v1 + f_v * _v2;
				*ptr1 = val;

				ptr1++;
			}
		}
		
		void _SLOPE2_z_x(double* ptr, double v1, double v2)
		{


			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			//v1 /= length;
			//v2 /= length;
			double f_u = 0, f_v = 0;
			double val = 0;
			double* ptr1 = ptr;

			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			for (int i = 0; i < _nNode; i++)
			{
				//double _g11 = 2 * _ref->d1[0][i] * _ref->get__gi(0, 0);
				//double _g12 = _ref->d1[0][i] * _ref->get__gi(1, 0) + _ref->d1[1][i] * _ref->get__gi(0, 0);
				//double _g22 = 2 * _ref->d1[1][i] * _ref->get__gi(1, 0);
				//double _g21 = _g12;
				double _g11 = __dsigma_11[0][i];
				double _g12 = __dsigma_12[0][i];
				double _g22 = __dsigma_22[0][i];
				double _g21 = _g12;

				double _length = 0.5 / (length) * (v1 * v1 * _g11 + 2 * v1 * v2 * _g12 + v2 * v2 * _g22);
				double _v1 = -v1 / length / length * _length;
				double _v2 = -v2 / length / length * _length;
				val = f_u * _v1 + f_v * _v2;
				*ptr1 = val;

				ptr1++;
			}
		}
		void _SLOPE2_z_y(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			//v1 /= length;
			//v2 /= length;
			double f_u = 0, f_v = 0;
			double val = 0;
			double* ptr1 = ptr;

			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			for (int i = 0; i < _nNode; i++)
			{
				//double _g11 = 2 * _ref->d1[0][i] * _ref->get__gi(0, 1);
				//double _g12 = _ref->d1[0][i] * _ref->get__gi(1, 1) + _ref->d1[1][i] * _ref->get__gi(0, 1);
				//double _g22 = 2 * _ref->d1[1][i] * _ref->get__gi(1, 1);
				//double _g21 = _g12;
				double _g11 = __dsigma_11[1][i];
				double _g12 = __dsigma_12[1][i];
				double _g22 = __dsigma_22[1][i];
				double _g21 = _g12;

				double _length = 0.5 / (length) * (v1 * v1 * _g11 + 2 * v1 * v2 * _g12 + v2 * v2 * _g22);
				double _v1 = -v1 / length / length * _length;
				double _v2 = -v2 / length / length * _length;
				val = f_u * _v1 + f_v * _v2;
				*ptr1 = val;

				ptr1++;
			}
		}
		double _SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;

		}
		double _SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			double f_U = f_u * _ref->get__Gij(0, 0) + f_v * _ref->get__Gij(1, 0);
			double f_V = f_u * _ref->get__Gij(0, 1) + f_v * _ref->get__Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1)*sc;// val;

		}
		double SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_phi[i];
				f_v += _ref->d1[1][i] * _ref->buf_phi[i];
			}
			double f_U = f_u * this->get_Gij(0, 0) + f_v * this->get_Gij(1, 0);
			double f_V = f_u * this->get_Gij(0, 1) + f_v * this->get_Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1)*sc;// val;

		}
		double SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {
			double f_u = 0, f_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				f_u += _ref->d1[0][i] * _ref->buf_z[i];
				f_v += _ref->d1[1][i] * _ref->buf_z[i];
			}
			double f_U = f_u * this->get_Gij(0, 0) + f_v * this->get_Gij(1, 0);
			double f_V = f_u * this->get_Gij(0, 1) + f_v * this->get_Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1)*sc;// val;

		}
		inline double SLOPE_xi(double dcdtstar0, double dcdtstar1) {

			double xi_u = 0, xi_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				xi_u += _ref->d1[0][i] * _ref->buf_xi[i];
				xi_v += _ref->d1[1][i] * _ref->buf_xi[i];
			}
			double f_U = xi_u * get_Gij(0, 0) + xi_v * get_Gij(1, 0);
			double f_V = xi_u * get_Gij(0, 1) + xi_v * get_Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
		inline double SLOPE_eta(double dcdtstar0, double dcdtstar1) {

			double eta_u = 0, eta_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				eta_u += _ref->d1[0][i] * _ref->buf_eta[i];
				eta_v += _ref->d1[1][i] * _ref->buf_eta[i];
			}
			double f_U = eta_u * get_Gij(0, 0) + eta_v * get_Gij(1, 0);
			double f_V = eta_u * get_Gij(0, 1) + eta_v * get_Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
		inline double SLOPE_nu(double dcdtstar0, double dcdtstar1) {

			double nu_u = 0, nu_v = 0;
			for (int i = 0; i < _nNode; i++)
			{
				nu_u += _ref->d1[0][i] * _ref->buf_nu[i];
				nu_v += _ref->d1[1][i] * _ref->buf_nu[i];
			}
			double f_U = nu_u * get_Gij(0, 0) + nu_v * get_Gij(1, 0);
			double f_V = nu_u * get_Gij(0, 1) + nu_v * get_Gij(1, 1);

			return (f_U * dcdtstar0 + f_V * dcdtstar1);// val;
		}
	
		inline double Dc(double a, double b, int l) {
			double f1 = a * this->get_gi(0, 0) + b * this->get_gi(0, 1);
			double f2 = a * this->get_gi(1, 0) + b * this->get_gi(1, 1);
			return f1 * get_Gij(0, l) + f2 * get_Gij(1, l);
		}
		inline double _Dc(double a, double b, int l) {
			double f1 = a * _ref->get__gi(0, 0) + b * _ref->get__gi(0, 1);
			double f2 = a * _ref->get__gi(1, 0) + b * _ref->get__gi(1, 1);
			return f1 * _ref->get__Gij(0, l) + f2 * _ref->get__Gij(1, l);
		}
		inline double rot1(double v1, double v2)
		{
			return v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
		}
		inline double rot2(double v1, double v2)
		{
			return -v1 * _ref->get__gij(0, 0) - v2 * _ref->get__gij(1, 0);

		}
		inline double _Dc2(double a, double b, double v1,double v2) {
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;

			double f1 = a * _ref->get__gi(0,0) + b * _ref->get__gi(0, 1);
			double f2 = a * _ref->get__gi(1, 0) + b * _ref->get__gi(1, 1);
			return f1 * v1 + f2 * v2;
		}
		
		void _Dc2_x(double* ptr, double a, double b, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double* ptr1 = ptr;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double  _f1 = _ref->d1[0][s];
				double _f2 = _ref->d1[1][s];

				double val = a*(_f1 * v1 + _f2 * v2);
				*ptr1 = val;
				ptr1++;
				
			}
		}
		void _Dc2_y(double* ptr, double a, double b, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double* ptr1 = ptr;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double  _f1 = _ref->d1[0][s];
				double _f2 = _ref->d1[1][s];

				double val = b*(_f1 * v1 + _f2 * v2);
				*ptr1 = val;
				ptr1++;
			
			}
		}

		double ___a2(double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;

			return _ref->get__gi(0, 0) * v1 + _ref->get__gi(1, 0) * v2;
		}
		double ___b2(double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;

			return _ref->get__gi(0, 1) * v1 + _ref->get__gi(1, 1) * v2;
		}

		inline double __a(int l)
		{
			return this->get_gi(0, 0) * get_Gij(0, l) + this->get_gi(1, 0) * get_Gij(1, l);
		}
		inline double __b(int l)
		{
			return this->get_gi(0, 1) * get_Gij(0, l) + this->get_gi(1, 1) * get_Gij(1, l);
		}
		inline double ___a(int l)
		{
			return _ref->get__gi(0, 0) * _ref->get__Gij(0, l) + _ref->get__gi(1, 0) * _ref->get__Gij(1, l);
		}
		inline double ___b(int l)
		{
			return _ref->get__gi(0, 1) * _ref->get__Gij(0, l) + _ref->get__gi(1, 1) * _ref->get__Gij(1, l);
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

		void laplacian(_mySparse* mat, int64_t* index, double sc)
		{
			for (int i = 0; i < _nNode; i++)
			{
				int I = index[i];
				for (int j = 0; j < _nNode; j++)
				{
					int J = index[j];
					int e = i * _nNode + j;
					double val = 0.5 * (_ref->get__Gij(0,0) * _ref->B[0][e] + 2 * _ref->get__Gij(0, 1) * _ref->B[1][e] + _ref->get__Gij(1, 1) * _ref->B[3][e]);
					mat->adddat(I, J, val*sc);
				}
			}
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
			double* pptr1 = &(_ref->buf_z[0]);
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				double _h11 = _ref->__dh[0][i];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][i];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][i];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;


				a += pptr1[i] * _h11;
				c += pptr1[i] * _h22;
				b += pptr1[i] * _h12;
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
		double _tracenu()
		{

			/*double A11 = _ref->get__Gij(0, 0) * get_fij(0, 0) * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * get_fij(0, 1) * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * get_fij(1, 0) * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * get_fij(1, 1) * _ref->get__Gij(1, 0);
			double A12 = _ref->get__Gij(0, 0) * get_fij(0, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * get_fij(0, 1) * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * get_fij(1, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * get_fij(1, 1) * _ref->get__Gij(1, 1);
			double A22 = _ref->get__Gij(1, 0) * get_fij(0, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * get_fij(0, 1) * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * get_fij(1, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * get_fij(1, 1) * _ref->get__Gij(1, 1);
			double a22 = A11 / sc;
			double a11 = A22 / sc;
			double a12 = -A12 / sc;
			double a21 = a12;*/

			double k11 = get_eij(0, 0);// +a11;
			double k22 = get_eij(1, 1);// +a22;
			double k12 = get_eij(0, 1);// +a12;
			double k21 = k12;

			return k11 * _ref->get__Gij(0, 0) + 2 * k12 * _ref->get__Gij(0, 1) + k22 * _ref->get__Gij(1, 1);
		}
	
		void _tracenu_nu(double* ptr)
		{

			double* ptr1 = ptr;
			/*double A11 = _ref->get__Gij(0, 0) * get_fij(0, 0) * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * get_fij(0, 1) * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * get_fij(1, 0) * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * get_fij(1, 1) * _ref->get__Gij(1, 0);
			double A12 = _ref->get__Gij(0, 0) * get_fij(0, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * get_fij(0, 1) * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * get_fij(1, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * get_fij(1, 1) * _ref->get__Gij(1, 1);
			double A22 = _ref->get__Gij(1, 0) * get_fij(0, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * get_fij(0, 1) * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * get_fij(1, 0) * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * get_fij(1, 1) * _ref->get__Gij(1, 1);
			double a22 = A11 / sc;
			double a11 = A22 / sc;
			double a12 = -A12 / sc;
			double a21 = a12;*/



			double k11 = get_eij(0, 0);// +a11;
			double k22 = get_eij(1, 1);//+ a22;
			double k12 = get_eij(0, 1);// +a12;
			double k21 = k12;



			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e22 = __dsigma_22[2][s];

				double _e21 = _e12;
				double _k11 = _e11, _k22 = _e22, _k12 = _e12;
				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 * _e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				*ptr1 =_k11 * _ref->get__Gij(0, 0) + 2 * _k12 * _ref->get__Gij(0, 1) + _k22 * _ref->get__Gij(1, 1);

				ptr1++;
			}
		}


		/*double _tracexi()
		{

			
			return get_fij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_fij(0, 1) * _ref->get__Gij(0,1) + get_fij(1, 1) * _ref->get__Gij(1, 1);
		}*/
		/*void _tracexi_xi(double* ptr)
		{

			double* ptr1 = ptr;
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				double _f11 = _ref->__dh[0][i] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _f12 = _ref->__dh[1][i] ;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _f22 = _ref->__dh[3][i] ;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;


				a = _f11;
				c = _f22;
				b = _f12;
				*ptr1 = a * _ref->get__Gij(0, 0) + 2 *b* _ref->get__Gij(0, 1) + c * _ref->get__Gij(1, 1);

				ptr1++;
			}
		}

		double _detsigma()
		{

		

			double k11 = get_mij(0, 0);// +a11;
			double k22 = get_mij(1, 1);// +a22;
			double k12 = get_mij(0, 1);// +a12;
			double k21 = k12;

			return k11 * k22 - k12 * k12;
		}
		void _detsigma_xi(double* ptr)
		{

			double* ptr1 = ptr;
		

			double k11 = get_mij(0, 0);// +a11;
			double k22 = get_mij(1, 1);// +a22;
			double k12 = get_mij(0, 1);// +a12;
			double k21 = k12;




			for (int s = 0; s < _ref->_nNode; s ++)
			{
				double _k11 = __dsigma_11[0][s];
				double _k12 = __dsigma_12[0][s];
				double _k22 = __dsigma_22[0][s];
				double _k21 = _k12;
			



				*ptr1 = _k11 * k22 + k11 * _k22 - 2 * k12 * _k12;
				ptr1++;
			}
		}
		void _detsigma_eta( double* ptr)
		{

			double* ptr1 = ptr;


			double k11 = get_mij(0, 0);// +a11;
			double k22 = get_mij(1, 1);// +a22;
			double k12 = get_mij(0, 1);// +a12;
			double k21 = k12;




			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _k11 = __dsigma_11[1][s];
				double _k12 = __dsigma_12[1][s];
				double _k22 = __dsigma_22[1][s];
				double _k21 = _k12;




				*ptr1 = _k11 * k22 + k11 * _k22 - 2 * k12 * _k12;
				ptr1++;
			}
		}
		void _detsigma_nu(double* ptr)
		{

			double* ptr1 = ptr;


			double k11 = get_mij(0, 0);// +a11;
			double k22 = get_mij(1, 1);// +a22;
			double k12 = get_mij(0, 1);// +a12;
			double k21 = k12;




			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _k11 = __dsigma_11[2][s];
				double _k12 = __dsigma_12[2][s];
				double _k22 = __dsigma_22[2][s];
				double _k21 = _k12;




				*ptr1 = _k11 * k22 + k11 * _k22 - 2 * k12 * _k12;
				ptr1++;
			}
		}*/
		double _detphi()
		{

			double* pptr1 = &(_ref->buf_phi[0]);
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				double _h11 = _ref->__dh[0][i] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][i] ;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][i] ;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;


				a += pptr1[i] * _h11;
				c += pptr1[i] * _h22;
				b += pptr1[i] * _h12;
			}
			return a * c - b * b;
		}
		void _detphi_phi(double *ptr)
		{

			double* ptr1 = ptr;
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				double _h11 = _ref->__dh[0][i] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][i] ;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][i] ;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;


				a = _h11;
				c = _h22;
				b = _h12;
				*ptr1=a * get__hij(1,1)+c* get__hij(0, 0) - 2*b * get__hij(0, 1);
				ptr1 ++;
			}
		}
		void _detZ_z(double* ptr)
		{

			double* ptr1 = ptr;
			double a = 0;
			double b = 0;
			double c = 0;


			for (int i = 0; i < _ref->_nNode; i++)
			{
				double _h11 = _ref->__dh[0][i];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][i];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][i];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;


				a = _h11;
				c = _h22;
				b = _h12;
				*ptr1 = a * get__Sij(1, 1) + c * get__Sij(0, 0) - 2 * b * get__Sij(0, 1);
				ptr1++;
			}
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
				for (int64_t j = 0; j < _nNode; j++)
				{
					int64_t J = index[j];
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
				int64_t I = index[i];
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
		}

		double st(double t1, double t2, double n1, double n2)
		{

			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;

			double val = get__hij(0, 0) * n1 * t1 + get__hij(0, 1) * n1 * t2 + get__hij(1, 0) * n2 * t1 + get__hij(1, 1) * n2 * t2;

			return val;
		}


		void st_phi(double* ptr, double t1, double t2, double n1, double n2)
		{
			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				double _h11 = _ref->__dh[0][s];
				double _h12 = _ref->__dh[1][s];
				double _h22 = _ref->__dh[3][s];

				double _h21 = _h12;



				val = _h11 * n1 * t1 + _h12 * n1 * t2 + _h12 * n2 * t1 + _h22 * n2 * t2;


				*ptr1 = val;
				ptr1++;
			}
		}
		
		double fair3(double v1, double v2, double w1, double w2)
		{
			double length = v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22;
			v1 /= length;
			v2 /= length;
			
			length = w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22;
			w1 /= length;
			w2 /= length;

			double e11 = 0, e12 = 0, e22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += _ref->d1[0][s] * _ref->d1[0][t] * _ref->buf_u[s] * _ref->buf_u[t];
					e11 += _ref->d1[0][s] * _ref->d1[0][t] * _ref->buf_v[s] * _ref->buf_v[t];
					e12 += _ref->d1[0][s] * _ref->d1[1][t] * _ref->buf_u[s] * _ref->buf_u[t];
					e12 += _ref->d1[0][s] * _ref->d1[1][t] * _ref->buf_v[s] * _ref->buf_v[t];
					e22 += _ref->d1[1][s] * _ref->d1[1][t] * _ref->buf_u[s] * _ref->buf_u[t];
					e22 += _ref->d1[1][s] * _ref->d1[1][t] * _ref->buf_v[s] * _ref->buf_v[t];
				}
			}
			
			return e11*v1*w1+e12*v1*w2+e12*v2*w1+e22*v2*w2;
		}
		void  fair3_u(double *ptr,double v1, double v2, double w1, double w2)
		{
			double length = v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22;
			v1 /= length;
			v2 /= length;

			length = w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22;
			w1 /= length;
			w2 /= length;

		

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					e11 += 2*_ref->d1[0][s] * _ref->d1[0][t]  * _ref->buf_u[t];
					//e11 += 2 * _ref->d1[0][s] * _ref->d1[0][t] * _ref->buf_v[t];
					e12 += _ref->d1[0][s] * _ref->d1[1][t] * _ref->buf_u[t]+ _ref->d1[1][s] * _ref->d1[0][t] * _ref->buf_u[t];
					//e12 += _ref->d1[0][s] * _ref->d1[1][t]  * _ref->buf_v[t]+ _ref->d1[1][s] * _ref->d1[0][t] * _ref->buf_v[t];
					e22 += 2 * _ref->d1[1][s] * _ref->d1[1][t] * _ref->buf_u[t];
					//e22 += 2 * _ref->d1[1][s] * _ref->d1[1][t]  * _ref->buf_v[t];
				}
				double val=e11 * v1 * w1 + e12 * v1 * w2 + e12 * v2 * w1 + e22 * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void  fair3_v(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22;
			v1 /= length;
			v2 /= length;

			length = w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22;
			w1 /= length;
			w2 /= length;



			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double e11 = 0, e12 = 0, e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					//e11 += 2 * _ref->d1[0][s] * _ref->d1[0][t] * _ref->buf_u[t];
					e11 += 2 * _ref->d1[0][s] * _ref->d1[0][t] * _ref->buf_v[t];
					//e12 += _ref->d1[0][s] * _ref->d1[1][t] * _ref->buf_u[t] + _ref->d1[1][s] * _ref->d1[0][t] * _ref->buf_u[t];
					e12 += _ref->d1[0][s] * _ref->d1[1][t] * _ref->buf_v[t] + _ref->d1[1][s] * _ref->d1[0][t] * _ref->buf_v[t];
					//e22 += 2 * _ref->d1[1][s] * _ref->d1[1][t] * _ref->buf_u[t];
					e22 += 2 * _ref->d1[1][s] * _ref->d1[1][t] * _ref->buf_v[t];
				}
				double val = e11 * v1 * w1 + e12 * v1 * w2 + e12 * v2 * w1 + e22 * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double fair4(double v1, double v2)
		{
			double length = v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22;
			v1 /= length;
			v2 /= length;

			double gttt = this->get_gammaijk(0, 0, 0) * v1 * v1 * v1 + this->get_gammaijk(0, 0, 1) * v1 * v1 * v2 +
				2 * (this->get_gammaijk(0, 1, 0) * v1 * v2 * v1 + this->get_gammaijk(0, 1, 1) * v1 * v2 * v2) +
				this->get_gammaijk(1, 1, 0) * v2 * v2 * v1 + this->get_gammaijk(1, 1, 1) * v2 * v2 * v2;
			double u111 = 0,u112=0,u121=0,u122=0,u221=0,u222=0;
			double v111 = 0,v112 = 0,v121 = 0,v122 = 0, v221 = 0, v222 = 0;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u111 += _ref->d3[0][s] * _ref->node[s * 3 + 0];
				u112 += _ref->d3[1][s] * _ref->node[s * 3 + 0];
				u121 += _ref->d3[2][s] * _ref->node[s * 3 + 0];
				u122 += _ref->d3[3][s] * _ref->node[s * 3 + 0];
				u221 += _ref->d3[6][s] * _ref->node[s * 3 + 0];
				u222 += _ref->d3[7][s] * _ref->node[s * 3 + 0];
				v111 += _ref->d3[0][s] * _ref->node[s * 3 + 1];
				v112 += _ref->d3[1][s] * _ref->node[s * 3 + 1];
				v121 += _ref->d3[2][s] * _ref->node[s * 3 + 1];
				v122 += _ref->d3[3][s] * _ref->node[s * 3 + 1];
				v221 += _ref->d3[6][s] * _ref->node[s * 3 + 1];
				v222 += _ref->d3[7][s] * _ref->node[s * 3 + 1];
			}
			double u211 = u121, u212 = u122;
			double v211 = v121, v212 = v122;

			double ut = _ref->_ogi[0] * v1 + _ref->_ogi[3] * v2;
			double vt = _ref->_ogi[1] * v1 + _ref->_ogi[4] * v2;
			double Gtt = 1;// v1* v1* _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1);

			double omegatttt = (u111 * v1 * v1 * v1 + u112 * v1 * v1 * v2 + 2 * (u121 * v1 * v2 * v1 + u122 * v1 * v2 * v2) +
				u221 * v2 * v2 * v1 + u222 * v2 * v2 * v2) * ut +
				(v111 * v1 * v1 * v1 + v112 * v1 * v1 * v2 + 2 * (v121 * v1 * v2 * v1 + v122 * v1 * v2 * v2) +
					v221 * v2 * v2 * v1 + v222 * v2 * v2 * v2) * vt;


			return 2 * omegatttt -gttt * gttt * Gtt;
		}

		void fair4_u(double* ptr, double v1, double v2)
		{
			double length = v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22;
			v1 /= length;
			v2 /= length;

			double gttt = this->get_gammaijk(0, 0, 0) * v1 * v1 * v1 + this->get_gammaijk(0, 0, 1) * v1 * v1 * v2 +
				2 * (this->get_gammaijk(0, 1, 0) * v1 * v2 * v1 + this->get_gammaijk(0, 1, 1) * v1 * v2 * v2) +
				this->get_gammaijk(1, 1, 0) * v2 * v2 * v1 + this->get_gammaijk(1, 1, 1) * v2 * v2 * v2;
			double u111 = 0, u112 = 0, u121 = 0, u122 = 0, u221 = 0, u222 = 0;
			double v111 = 0, v112 = 0, v121 = 0, v122 = 0, v221 = 0, v222 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				u111 += _ref->d3[0][s] * _ref->node[s * 3 + 0];
				u112 += _ref->d3[1][s] * _ref->node[s * 3 + 0];
				u121 += _ref->d3[2][s] * _ref->node[s * 3 + 0];
				u122 += _ref->d3[3][s] * _ref->node[s * 3 + 0];
				u221 += _ref->d3[6][s] * _ref->node[s * 3 + 0];
				u222 += _ref->d3[7][s] * _ref->node[s * 3 + 0];
				v111 += _ref->d3[0][s] * _ref->node[s * 3 + 1];
				v112 += _ref->d3[1][s] * _ref->node[s * 3 + 1];
				v121 += _ref->d3[2][s] * _ref->node[s * 3 + 1];
				v122 += _ref->d3[3][s] * _ref->node[s * 3 + 1];
				v221 += _ref->d3[6][s] * _ref->node[s * 3 + 1];
				v222 += _ref->d3[7][s] * _ref->node[s * 3 + 1];
			}
			double u211 = u121, u212 = u122;
			double v211 = v121, v212 = v122;

			double ut = _ref->_ogi[0] * v1 + _ref->_ogi[3] * v2;
			double vt = _ref->_ogi[1] * v1 + _ref->_ogi[4] * v2; 
			
			double Gtt = 1;// v1* v1* _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1);
			
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _ut = _ref->d1[0][s] * v1 + _ref->d1[1][s] * v2;
				double _vt = 0;// _ref->get__gi(0, 1)* v1 + _ref->get__gi(1, 1) * v2;

				double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 0)+ _ref->d1[1][s] * _ref->get__gi(0, 0);
				double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);
				double _g21 = _g12;
				double _u111 = _ref->d3[0][s] ;
				double _u112 = _ref->d3[1][s] ;
				double _u121 = _ref->d3[2][s] ;
				double _u122 = _ref->d3[3][s] ;
				double _u221 = _ref->d3[6][s] ;
				double _u222 = _ref->d3[7][s] ;
				double _v111 = 0;//_ref->d3[0][s] * _ref->buf_v[s];
				double _v112 = 0;//_ref->d3[1][s] * _ref->buf_v[s];
				double _v121 = 0;//_ref->d3[2][s] * _ref->buf_v[s];
				double _v122 = 0;//_ref->d3[3][s] * _ref->buf_v[s];
				double _v221 = 0;//_ref->d3[6][s] * _ref->buf_v[s];
				double _v222 = 0;// _ref->d3[7][s] * _ref->buf_v[s];

				double _gamma111 = _ref->__dh[0][s] * _ref->get__gi(0, 0);
				double _gamma112 = _ref->__dh[0][s] * _ref->get__gi(1, 0);
				double _gamma121 = _ref->__dh[1][s] * _ref->get__gi(0, 0);
				double _gamma122 = _ref->__dh[1][s] * _ref->get__gi(1, 0);
				double _gamma221 = _ref->__dh[3][s] * _ref->get__gi(0, 0);
				double _gamma222 = _ref->__dh[3][s] * _ref->get__gi(1, 0);

				double _gttt = _gamma111 * v1 * v1 * v1 + _gamma112 * v1 * v1 * v2 +
					2 * (_gamma121 * v1 * v2 * v1 + _gamma122 * v1 * v2 * v2) +
					_gamma221 * v2 * v2 * v1 + _gamma222 * v2 * v2 * v2;


				//double _Gtt = v1 * v1 * _g11 + 2 * v1 * v2 * _g12 + v2 * v2 * _g22;

				double _omegatttt = (_u111 * v1 * v1 * v1 + _u112 * v1 * v1 * v2 + 2 * (_u121 * v1 * v2 * v1 + _u122 * v1 * v2 * v2) +
					_u221 * v2 * v2 * v1 + _u222 * v2 * v2 * v2) * ut +
					(_v111 * v1 * v1 * v1 + _v112 * v1 * v1 * v2 + 2 * (_v121 * v1 * v2 * v1 + _v122 * v1 * v2 * v2) +
						_v221 * v2 * v2 * v1 + _v222 * v2 * v2 * v2) * vt;
		
				//_omegatttt += (u111 * v1 * v1 * v1 + u112 * v1 * v1 * v2 + 2 * (u121 * v1 * v2 * v1 + u122 * v1 * v2 * v2) +
				//	u221 * v2 * v2 * v1 + u222 * v2 * v2 * v2) * _ut +
					//(v111 * v1 * v1 * v1 + v112 * v1 * v1 * v2 + 2 * (v121 * v1 * v2 * v1 + v122 * v1 * v2 * v2) +
						//v221 * v2 * v2 * v1 + v222 * v2 * v2 * v2) * _vt;

				double val = 2 * _omegatttt - 2 * _gttt * gttt * Gtt;// -2 * gttt * gttt * _Gtt;
				*ptr1 = val;
				ptr1++;
			}
		}
		void fair4_v(double* ptr, double v1, double v2)
		{
			double length = v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22;
			v1 /= length;
			v2 /= length;

			double gttt = this->get_gammaijk(0, 0, 0) * v1 * v1 * v1 + this->get_gammaijk(0, 0, 1) * v1 * v1 * v2 +
				2 * (this->get_gammaijk(0, 1, 0) * v1 * v2 * v1 + this->get_gammaijk(0, 1, 1) * v1 * v2 * v2) +
				this->get_gammaijk(1, 1, 0) * v2 * v2 * v1 + this->get_gammaijk(1, 1, 1) * v2 * v2 * v2;
			double u111 = 0, u112 = 0, u121 = 0, u122 = 0, u221 = 0, u222 = 0;
			double v111 = 0, v112 = 0, v121 = 0, v122 = 0, v221 = 0, v222 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				u111 += _ref->d3[0][s] * _ref->node[s * 3 + 0];
				u112 += _ref->d3[1][s] * _ref->node[s * 3 + 0];
				u121 += _ref->d3[2][s] * _ref->node[s * 3 + 0];
				u122 += _ref->d3[3][s] * _ref->node[s * 3 + 0];
				u221 += _ref->d3[6][s] * _ref->node[s * 3 + 0];
				u222 += _ref->d3[7][s] * _ref->node[s * 3 + 0];
				v111 += _ref->d3[0][s] * _ref->node[s * 3 + 1];
				v112 += _ref->d3[1][s] * _ref->node[s * 3 + 1];
				v121 += _ref->d3[2][s] * _ref->node[s * 3 + 1];
				v122 += _ref->d3[3][s] * _ref->node[s * 3 + 1];
				v221 += _ref->d3[6][s] * _ref->node[s * 3 + 1];
				v222 += _ref->d3[7][s] * _ref->node[s * 3 + 1];
			}
			double u211 = u121, u212 = u122;
			double v211 = v121, v212 = v122;

			double ut = _ref->_ogi[0] * v1 + _ref->_ogi[3] * v2;
			double vt = _ref->_ogi[1] * v1 + _ref->_ogi[4] * v2;
			double Gtt = 1;// v1* v1* _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);
				double _g21 = _g12;

				double _ut = 0;
				double _vt = _ref->d1[0][s] * v1 + _ref->d1[1][s] * v2;

				double _u111 = 0;// _ref->d3[0][s];
				double _u112 = 0;//_ref->d3[1][s];
				double _u121 = 0;//_ref->d3[2][s];
				double _u122 = 0;//_ref->d3[3][s];
				double _u221 = 0;//_ref->d3[6][s];
				double _u222 = 0;//_ref->d3[7][s];
				double _v111 = _ref->d3[0][s];// *_ref->buf_v[s];
				double _v112 = _ref->d3[1][s];// * _ref->buf_v[s];
				double _v121 = _ref->d3[2][s];// * _ref->buf_v[s];
				double _v122 = _ref->d3[3][s];// * _ref->buf_v[s];
				double _v221 = _ref->d3[6][s];// * _ref->buf_v[s];
				double _v222 = _ref->d3[7][s];// * _ref->buf_v[s];

				double _gamma111 = _ref->__dh[0][s] * _ref->get__gi(0, 1);
				double _gamma112 = _ref->__dh[0][s] * _ref->get__gi(1, 1);
				double _gamma121 = _ref->__dh[1][s] * _ref->get__gi(0, 1);
				double _gamma122 = _ref->__dh[1][s] * _ref->get__gi(1, 1);
				double _gamma221 = _ref->__dh[3][s] * _ref->get__gi(0, 1);
				double _gamma222 = _ref->__dh[3][s] * _ref->get__gi(1, 1);

				double _gttt = _gamma111 * v1 * v1 * v1 + _gamma112 * v1 * v1 * v2 +
					2 * (_gamma121 * v1 * v2 * v1 + _gamma122 * v1 * v2 * v2) +
					_gamma221 * v2 * v2 * v1 + _gamma222 * v2 * v2 * v2;

				//double _Gtt = v1 * v1 * _g11 + 2 * v1 * v2 * _g12 + v2 * v2 * _g22;

				double _omegatttt = (_u111 * v1 * v1 * v1 + _u112 * v1 * v1 * v2 + 2 * (_u121 * v1 * v2 * v1 + _u122 * v1 * v2 * v2) +
					_u221 * v2 * v2 * v1 + _u222 * v2 * v2 * v2) * ut +
					(_v111 * v1 * v1 * v1 + _v112 * v1 * v1 * v2 + 2 * (_v121 * v1 * v2 * v1 + _v122 * v1 * v2 * v2) +
						_v221 * v2 * v2 * v1 + _v222 * v2 * v2 * v2) * vt;
				//_omegatttt += (u111 * v1 * v1 * v1 + u112 * v1 * v1 * v2 + 2 * (u121 * v1 * v2 * v1 + u122 * v1 * v2 * v2) +
				//	u221 * v2 * v2 * v1 + u222 * v2 * v2 * v2)* _ut +
				//	(v111 * v1 * v1 * v1 + v112 * v1 * v1 * v2 + 2 * (v121 * v1 * v2 * v1 + v122 * v1 * v2 * v2) +
				//		v221 * v2 * v2 * v1 + v222 * v2 * v2 * v2) * _vt;
				

				double val = 2 * _omegatttt - 2 * _gttt * gttt * Gtt;// -2 * gttt * gttt * _Gtt;
				*ptr1 = val;
				ptr1++;
			}
		}
		void fair_u(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;

			double s1 = (v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1)) / _dv;
			double s2 = -(v1 * _ref->get__gij(0, 0) + v2 * _ref->get__gij(1, 0)) / _dv;


			double g111 = 0, g112 = 0, g121 = 0, g122 = 0, g221 = 0, g222 = 0;
			g111 = _ref->get__Gammaijk(0, 0, 0);
			g112 = _ref->get__Gammaijk(0, 0, 1);
			g121 = _ref->get__Gammaijk(0, 1, 0);
			g122 = _ref->get__Gammaijk(0, 1, 1);
			g221 = _ref->get__Gammaijk(1, 1, 0);
			g222 = _ref->get__Gammaijk(1, 1, 1);

			double A1 = (g111 * v1 * v1 + g121 * (v1 * v2 + v2 * v1) + g221 * v2 * v2);
			double A2 = (g112 * v1 * v1 + g122 * (v1 * v2 + v2 * v1) + g222 * v2 * v2);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _g112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _g121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _g122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _g221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _g222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);

				double _g11 = 2 * _ref->get__gi(0, 0) * _ref->d1[0][s];
				double _g12 = _ref->get__gi(0, 0) * _ref->d1[1][s] + _ref->get__gi(1, 0) * _ref->d1[0][s];
				double _g22 = 2 * _ref->get__gi(1, 0) * _ref->d1[1][s];


				double _A1 = (_g111 * v1 * v1 + _g121 * (v1 * v2 + v2 * v1) + _g221 * v2 * v2);
				double _A2 = (_g112 * v1 * v1 + _g122 * (v1 * v2 + v2 * v1) + _g222 * v2 * v2);
				
				double __dv = 0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0, 1) + _g22 * _ref->get__Gij(1, 1))*_dv;

				double _s1 = (v1 *_g12 + v2 * _g22) / _dv-s1/_dv*__dv;
				double _s2 = -(v1 * _g11 + v2 * _g12) / _dv - s2 / _dv * __dv;

			

				double val = 0;
				val = _A1 * s1 + _A2 * s2;
				val += A1 * _s1 + A2 * _s2;


				*ptr1 = val;
				ptr1++;
			}
		}
		
		void fair_v(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;

			double s1 = (v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1)) / _dv;
			double s2 = -(v1 * _ref->get__gij(0, 0) + v2 * _ref->get__gij(1, 0)) / _dv;


			double g111 = 0, g112 = 0, g121 = 0, g122 = 0, g221 = 0, g222 = 0;
			g111 = _ref->get__Gammaijk(0, 0, 0);
			g112 = _ref->get__Gammaijk(0, 0, 1);
			g121 = _ref->get__Gammaijk(0, 1, 0);
			g122 = _ref->get__Gammaijk(0, 1, 1);
			g221 = _ref->get__Gammaijk(1, 1, 0);
			g222 = _ref->get__Gammaijk(1, 1, 1);

			double A1 = (g111 * v1 * v1 + g121 * (v1 * v2 + v2 * v1) + g221 * v2 * v2);
			double A2 = (g112 * v1 * v1 + g122 * (v1 * v2 + v2 * v1) + g222 * v2 * v2);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _g112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _g121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _g122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _g221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _g222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);

				double _g11 = 2 * _ref->get__gi(0, 1) * _ref->d1[0][s];
				double _g12 = _ref->get__gi(0, 1) * _ref->d1[1][s] + _ref->get__gi(1, 1) * _ref->d1[0][s];
				double _g22 = 2 * _ref->get__gi(1, 1) * _ref->d1[1][s];


				double _A1 = (_g111 * v1 * v1 + _g121 * (v1 * v2 + v2 * v1) + _g221 * v2 * v2);
				double _A2 = (_g112 * v1 * v1 + _g122 * (v1 * v2 + v2 * v1) + _g222 * v2 * v2);

				double __dv = 0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0, 1) + _g22 * _ref->get__Gij(1, 1)) * _dv;

				double _s1 = (v1 * _g12 + v2 * _g22) / _dv - s1 / _dv * __dv;
				double _s2 = -(v1 * _g11 + v2 * _g12) / _dv - s2 / _dv * __dv;



				double val = 0;
				val = _A1 * s1 + _A2 * s2;
				val += A1 * _s1 + A2 * _s2;


				*ptr1 = val;
				ptr1++;
			}
		}
		double gamma(double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;

			double gtt = _ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2;

			double dttt= 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				dttt += (_ref->__dh[0][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v1 * v1 * v1;
				dttt += (_ref->__dh[0][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v1 * v1 * v1;
				dttt += (_ref->__dh[0][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v1 * v1 * v2;
				dttt += (_ref->__dh[0][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v1 * v1 * v2;
				dttt += (_ref->__dh[1][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v1 * v2 * v1;
				dttt += (_ref->__dh[1][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v1 * v2 * v1;
				dttt += (_ref->__dh[1][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v1 * v2 * v2;
				dttt += (_ref->__dh[1][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v1 * v2 * v2;
				dttt += (_ref->__dh[2][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v2 * v1 * v1;
				dttt += (_ref->__dh[2][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v2 * v1 * v1;
				dttt += (_ref->__dh[2][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v2 * v1 * v2;
				dttt += (_ref->__dh[2][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v2 * v1 * v2;
				dttt += (_ref->__dh[3][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v2 * v2 * v1;
				dttt += (_ref->__dh[3][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v2 * v2 * v1;
				dttt += (_ref->__dh[3][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v2 * v2 * v2;
				dttt += (_ref->__dh[3][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v2 * v2 * v2;
			}
			
			return dttt;
		}
		void gamma_u(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double* ptr1 = ptr;
			double utt = 0;
			double dttt = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				utt += (_ref->__dh[0][s] * _ref->buf_u[s]) * v1 * v1;
				utt += 2*(_ref->__dh[1][s] * _ref->buf_u[s]) * v1 * v2;
				utt += (_ref->__dh[3][s] * _ref->buf_u[s]) * v2 * v2;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double dttt = 0;
				dttt += (_ref->__dh[0][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v1 * v1 * v1;
				//dttt += (_ref->__dh[0][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v1 * v1 * v1;
				dttt += (_ref->__dh[0][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v1 * v1 * v2;
				//dttt += (_ref->__dh[0][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v1 * v1 * v2;
				dttt += (_ref->__dh[1][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v1 * v2 * v1;
				//dttt += (_ref->__dh[1][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v1 * v2 * v1;
				dttt += (_ref->__dh[1][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v1 * v2 * v2;
				//dttt += (_ref->__dh[1][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v1 * v2 * v2;
				dttt += (_ref->__dh[2][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v2 * v1 * v1;
				//dttt += (_ref->__dh[2][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v2 * v1 * v1;
				dttt += (_ref->__dh[2][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v2 * v1 * v2;
				//dttt += (_ref->__dh[2][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v2 * v1 * v2;
				dttt += (_ref->__dh[3][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v2 * v2 * v1;
				//dttt += (_ref->__dh[3][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v2 * v2 * v1;
				dttt += (_ref->__dh[3][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v2 * v2 * v2;
				//dttt += (_ref->__dh[3][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v2 * v2 * v2;
				dttt += utt * _ref->d1[0][s] * v1 + utt * _ref->d1[1][s] * v2;
				*ptr1 = dttt;
				ptr1++;
			}
		}
		void gamma_v(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double* ptr1 = ptr;
			double vtt = 0;
			double dttt = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				vtt += (_ref->__dh[0][s] * _ref->buf_v[s]) * v1 * v1;
				vtt += 2 * (_ref->__dh[1][s] * _ref->buf_v[s]) * v1 * v2;
				vtt += (_ref->__dh[3][s] * _ref->buf_v[s]) * v2 * v2;
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double dttt = 0;
				//dttt += (_ref->__dh[0][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v1 * v1 * v1;
				dttt += (_ref->__dh[0][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v1 * v1 * v1;
				//dttt += (_ref->__dh[0][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v1 * v1 * v2;
				dttt += (_ref->__dh[0][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v1 * v1 * v2;
				//dttt += (_ref->__dh[1][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v1 * v2 * v1;
				dttt += (_ref->__dh[1][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v1 * v2 * v1;
				//dttt += (_ref->__dh[1][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v1 * v2 * v2;
				dttt += (_ref->__dh[1][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v1 * v2 * v2;
				//dttt += (_ref->__dh[2][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v2 * v1 * v1;
				dttt += (_ref->__dh[2][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v2 * v1 * v1;
				//dttt += (_ref->__dh[2][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v2 * v1 * v2;
				dttt += (_ref->__dh[2][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v2 * v1 * v2;
				//dttt += (_ref->__dh[3][s] * _ref->buf_u[s]) * _ref->get__gi(0, 0) * v2 * v2 * v1;
				dttt += (_ref->__dh[3][s] * _ref->buf_v[s]) * _ref->get__gi(0, 1) * v2 * v2 * v1;
				//dttt += (_ref->__dh[3][s] * _ref->buf_u[s]) * _ref->get__gi(1, 0) * v2 * v2 * v2;
				dttt += (_ref->__dh[3][s] * _ref->buf_v[s]) * _ref->get__gi(1, 1) * v2 * v2 * v2;
				dttt += vtt * _ref->d1[0][s] * v1 + vtt * _ref->d1[1][s] * v2;
				*ptr1 = dttt;
				ptr1++;
			}
		}
		double length(double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double val = _ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2;
			return val;
		}
		void length_u(double *ptr,double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 0)+ _ref->d1[1][s] * _ref->get__gi(0, 0);
				double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);

				double val = _g11 * v1 * v1 + 2 *_g12* v1 * v2 + _g22 * v2 * v2;
				*ptr1=val;
				ptr1++;
			}
		}
		void length_v(double* ptr, double v1, double v2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

				double val = _g11 * v1 * v1 + 2 * _g12 * v1 * v2 + _g22 * v2 * v2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double fair2(double v1, double v2, double w1, double w2, double s1, double s2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			double g111 = 0, g112 = 0, g121 = 0, g122 = 0, g221 = 0, g222 = 0;
			g111 = _ref->get__Gammaijk(0, 0, 0);
			g112 = _ref->get__Gammaijk(0, 0, 1);
			g121 = _ref->get__Gammaijk(0, 1, 0);
			g122 = _ref->get__Gammaijk(0, 1, 1);
			g221 = _ref->get__Gammaijk(1, 1, 0);
			g222 = _ref->get__Gammaijk(1, 1, 1);

			double g211 = g121, g212 = g122;
			double val = 0;
			
			double A1 = (g111 * v1 * w1 + g121 * (v1 * w2 + v2 * w1) + g221 * v2 * w2);
			double A2 = (g112 * v1 * w1 + g122 * (v1 * w2 + v2 * w1) + g222 * v2 * w2);
			double B1 = (_ref->oGammaijk[0] * v1 * w1 + _ref->oGammaijk[2] * (v1 * w2 + v2 * w1) + _ref->oGammaijk[6] * v2 * w2);
			double B2 = (_ref->oGammaijk[1] * v1 * w1 + _ref->oGammaijk[3] * (v1 * w2 + v2 * w1) + _ref->oGammaijk[7] * v2 * w2);

			val = A1 * B2 - A2 * B1;

			return val;
		}
		void fair2_u(double* ptr, double v1, double v2, double w1, double w2, double s1, double s2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;

			double B1 = (_ref->oGammaijk[0] * v1 * w1 + _ref->oGammaijk[2] * (v1 * w2 + v2 * w1) + _ref->oGammaijk[6] * v2 * w2);
			double B2 = (_ref->oGammaijk[1] * v1 * w1 + _ref->oGammaijk[3] * (v1 * w2 + v2 * w1) + _ref->oGammaijk[7] * v2 * w2);


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double g111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double g112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double g121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double g122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double g221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double g222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);

				double A1 = (g111 * v1 * w1 + g121 * (v1 * w2 + v2 * w1) + g221 * v2 * w2);
				double A2 = (g112 * v1 * w1 + g122 * (v1 * w2 + v2 * w1) + g222 * v2 * w2);

				double val = 0;
				val = A1 * B2 - A2 * B1;

				*ptr1 = val;
				ptr1++;
			}
		}
		void fair2_v(double* ptr, double v1, double v2, double w1, double w2, double s1, double s2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;

			
			double B1 = (_ref->oGammaijk[0] * v1 * w1 + _ref->oGammaijk[2] * (v1 * w2 + v2 * w1) + _ref->oGammaijk[6] * v2 * w2);
			double B2 = (_ref->oGammaijk[1] * v1 * w1 + _ref->oGammaijk[3] * (v1 * w2 + v2 * w1) + _ref->oGammaijk[7] * v2 * w2);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double g111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double g112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double g121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double g122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double g221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double g222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);


				double g211 = g121, g212 = g122;
				double val = 0;
				double A1 = (g111 * v1 * w1 + g121 * (v1 * w2 + v2 * w1) + g221 * v2 * w2);
				double A2 = (g112 * v1 * w1 + g122 * (v1 * w2 + v2 * w1) + g222 * v2 * w2);

				val = A1 * B2 - A2 * B1;
				*ptr1 = val;
				ptr1++;
			}
		}
		double st2(double t1,double t2,double n1,double n2)
		{
			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;

			double val = get__Sij(0, 0) * n1 * t1 +  get__Sij(0, 1) * n1 * t2 + get__Sij(1, 0) * n2 * t1 + get__Sij(1, 1) * n2 * t2;
			
			return val;
		}


		void st2_z(double* ptr,double t1,double t2,double n1,double n2)
		{

			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;

				/*double _S11 = _ref->__dh[0][s];
				double _S12 = _ref->__dh[1][s];
				double _S22 = _ref->__dh[3][s];

				double _S21 = _S12;*/

				double _S11 = _ref->__dh[0][s];
				double _S12 = _ref->__dh[1][s];
				double _S22 = _ref->__dh[3][s];

				double _S21 = _S12;

				val = _S11 * n1 * t1 + _S12 * n1 * t2+ _S12 * n2 * t1 + _S22 * n2 *t2;
				

				*ptr1 = val;
				ptr1++;
			}
		}

		void st2_x(double* ptr, double t1, double t2, double n1, double n2)
		{

			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;

			

			double S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}


			double* ptr1 = ptr;
			double val = 0;
			



			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;



				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;
		

				val = _Suu * n1 * t1 + _Suv * n1 * t2 + _Suv * n2 * t1 + _Svv * n2 * t2;


				*ptr1 = val;
				ptr1++;
			}
		}

		void st2_y(double* ptr, double t1, double t2, double n1, double n2)
		{

			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;



			double S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}


			double* ptr1 = ptr;
			double val = 0;




			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val = 0;



				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;


				val = _Suu * n1 * t1 + _Suv * n1 * t2 + _Suv * n2 * t1 + _Svv * n2 * t2;


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
		double sx(double v1, double v2)
		{
			double sx = v1 * get_gi(0,0) + v2 * get_gi(1, 0);
			return sx;
		}
		double sy(double v1, double v2)
		{

			double sy = v1 * get_gi(0, 1) + v2 * get_gi(1, 1);
			return sy;
		}
		double symm_sigma(double sx,double sy,double __xi, double __eta)
		{
			double _x = xi - __xi;
			double _y = eta - __eta;

			return _x* sy - _y * sx;

		}
		void symm_sigma_xi(double* ptr, double sx, double sy)
		{
			//double _x = xi - __xi;
			//double _y = eta - __eta;
			//double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			//double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi = _ref->d0[s];
				val = _xi * sy;
				*ptr1 = val;
				ptr1++;
			}
		}
		void symm_sigma_eta(double* ptr,  double sx, double sy)
		{
			//double _x = xi - _ref->_xi;
			//double _y = eta - _ref->_eta;
			//double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			//double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _eta = _ref->d0[s];
				val = -_eta * sx;
				*ptr1 = val;
				ptr1++;
			}
		}


		double symm_gurtin(double sx, double sy, double __xi, double __eta)
		{
			double _x = phi - __xi;
			double _y = nu - __eta;

			return _x * sy - _y * sx;

		}
		void symm_gurtin_phi(double* ptr, double sx, double sy)
		{
			//double _x = xi - __xi;
			//double _y = eta - __eta;
			//double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			//double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi = _ref->d0[s];
				val = _xi * sy;
				*ptr1 = val;
				ptr1++;
			}
		}
		void symm_gurtin_nu(double* ptr, double sx, double sy)
		{
			//double _x = xi - _ref->_xi;
			//double _y = eta - _ref->_eta;
			//double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			//double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _eta = _ref->d0[s];
				val = -_eta * sx;
				*ptr1 = val;
				ptr1++;
			}
		}
		double symm(double sx, double sy, double __x, double __y)
		{
			double length = sqrt(sx * sx + sy * sy);
			sx /= length;
			sy /= length;
			double u = 0, v = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u += _ref->d0[s] * _ref->buf_u[s];
				v += _ref->d0[s] * _ref->buf_v[s];
			}

			double _x = u - __x;
			double _y = v - __y;

			return _x * sy - _y * sx;

		}
		void symm_x(double* ptr, double sx, double sy)
		{
			double length = sqrt(sx * sx + sy * sy);
			sx /= length;
			sy /= length;
			//double _x = xi - __xi;
			//double _y = eta - __eta;
			//double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			//double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xi = _ref->d0[s];
				val = _xi * sy;
				*ptr1 = val;
				ptr1++;
			}
		}
		void symm_y(double* ptr, double sx, double sy)
		{
			double length = sqrt(sx * sx + sy * sy);
			sx /= length;
			sy /= length;
			//double _x = xi - _ref->_xi;
			//double _y = eta - _ref->_eta;
			//double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			//double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _eta = _ref->d0[s];
				val = -_eta * sx;
				*ptr1 = val;
				ptr1++;
			}
		}

		double symm2(double v1, double v2)
		{
			double V1 = v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0);
			double V2 = v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1);

			double ex=0;// = v1 * _ref->get__Gi(0, 0) + v2 * _ref->get__Gi(1, 0);
			double ey = 0;// = v1 * _ref->get__Gi(0, 1) + v2 * _ref->get__Gi(1, 1);
			for (int i = 0; i < _ref->_nNode; i++)
			{
				ex += V1 * _ref->d1[0][i] * _ref->buf_xi[i] + V2 * _ref->d1[1][i] * _ref->buf_xi[i];
				ey += V1 * _ref->d1[0][i] * _ref->buf_eta[i] + V2 * _ref->d1[1][i] * _ref->buf_eta[i];
			}
			double sx = V1 * _ref->get__gi(0, 0) + V2 * _ref->get__gi(1, 0);
			double sy = V1 * _ref->get__gi(0, 1) + V2 * _ref->get__gi(1, 1);

			return ex * sy - ey * sx;

		}
		void symm2_xi(double* ptr, double v1, double v2)
		{
			double V1 = v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0);
			double V2 = v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1);
			double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _ex = V1 * _ref->d1[0][s]  + V2 * _ref->d1[1][s] ;
				val = _ex * sy;
				*ptr1 = val;
				ptr1++;
			}
		}
		void symm2_eta(double* ptr, double v1, double v2)
		{
			double V1 = v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0);
			double V2 = v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1);
			double sx = v1 * _ref->get__gi(0, 0) + v2 * _ref->get__gi(1, 0);
			double sy = v1 * _ref->get__gi(0, 1) + v2 * _ref->get__gi(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _ey = V1 * _ref->d1[0][s]  + V2 * _ref->d1[1][s];
				val = -_ey * sx;
				*ptr1 = val;
				ptr1++;
			}
		}

		double symm3(double v1, double v2)
		{
			double V1 = v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0);
			double V2 = v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1);

			double ez = 0;// = v1 * _ref->get__Gi(0, 0) + v2 * _ref->get__Gi(1, 0);
			for (int i = 0; i < _ref->_nNode; i++)
			{
				ez += V1 * _ref->d1[0][i] * _ref->buf_nu[i] + V2 * _ref->d1[1][i] * _ref->buf_nu[i];
			
			}

			return ez;

		}
		void symm3_nu(double* ptr, double v1, double v2)
		{
			double V1 = v1 * _ref->get__Gij(0, 0) + v2 * _ref->get__Gij(1, 0);
			double V2 = v1 * _ref->get__Gij(0, 1) + v2 * _ref->get__Gij(1, 1);
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _ez = V1 * _ref->d1[0][s] + V2 * _ref->d1[1][s];
				val = _ez;
				*ptr1 = val;
				ptr1++;
			}
		}
		double tt_2(double v1, double v2, bool accurate)
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
			double val = 0;
			val += v1 * S11 * v1;
			val += v1 * S12 * v2;
			val += v2 * S21 * v1;
			val += v2 * S22 * v2;
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
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
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

				S11 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				S12 = (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				S22 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);


				double S21 = S12;
				val += v1 * S11 * v1;
				val += v1 * S12 * v2;
				val += v2 * S21 * v1;
				val += v2 * S22 * v2;
				*ptr1 = val / norm;
				ptr1++;
			}
		}
		/*double __AREA_sigma()
		{
			//double tr= (get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1));
			return sqrt(get_eij(0, 0)* get_eij(1, 1)- get_eij(0, 1)* get_eij(0, 1));
		}
		double AREA_sigma(_myDoubleArray *v, int64_t* index1, int64_t* index2, int64_t* index3,double sc) {
			double _area = __AREA_sigma();
			
			double det = get_eij(0, 0) * get_eij(1, 1) - get_eij(0, 1) * get_eij(0, 1);
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				double ddet = get_eij(0, 0) * _e22 + _e11 * get_eij(1, 1) - 2 * _e12 * get_eij(0, 1);
			
				v->__v(index1[s]) = 0.5 / _area * ddet;
				
				
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				double ddet = get_eij(0, 0) * _e22 + _e11 * get_eij(1, 1) - 2 * _e12 * get_eij(0, 1);

				v->__v(index2[s]) = 0.5 / _area * ddet;

			}
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
				}
				double ddet = get_eij(0, 0) * _e22 + _e11 * get_eij(1, 1) - 2 * _e12 * get_eij(0, 1);

				v->__v(index3[s]) = 0.5 / _area * ddet;

			}
			return _area * sc;
		}

		double __Normal_sigma_X()
		{
			double g1x = 0, g1y = 0, g1z = 0;
			double g2x = 0, g2y = 0, g2z = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{				
				g1x += _ref->d1[0][s] * _ref->buf_xi[s];
				g1y += _ref->d1[0][s] * _ref->buf_eta[s];
				g1z += _ref->d1[0][s] * _ref->buf_nu[s];
				g2x += _ref->d1[1][s] * _ref->buf_xi[s];
				g2y += _ref->d1[1][s] * _ref->buf_eta[s];
				g2z += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double val = g1y * g2z - g1z * g2y;
			return val;
		}
		double __Normal_sigma_Y()
		{
			double g1x = 0, g1y = 0, g1z = 0;
			double g2x = 0, g2y = 0, g2z = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				g1x += _ref->d1[0][s] * _ref->buf_xi[s];
				g1y += _ref->d1[0][s] * _ref->buf_eta[s];
				g1z += _ref->d1[0][s] * _ref->buf_nu[s];
				g2x += _ref->d1[1][s] * _ref->buf_xi[s];
				g2y += _ref->d1[1][s] * _ref->buf_eta[s];
				g2z += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double val = g1z * g2x - g1x * g2z;
			return val;
		}
		double __Normal_sigma_Z()
		{
			double g1x = 0, g1y = 0, g1z = 0;
			double g2x = 0, g2y = 0, g2z = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				g1x += _ref->d1[0][s] * _ref->buf_xi[s];
				g1y += _ref->d1[0][s] * _ref->buf_eta[s];
				g1z += _ref->d1[0][s] * _ref->buf_nu[s];
				g2x += _ref->d1[1][s] * _ref->buf_xi[s];
				g2y += _ref->d1[1][s] * _ref->buf_eta[s];
				g2z += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double val = g1x* g2y - g1y * g2x;
			return val;
		}
		
		double Normal_sigma_X(_myDoubleArray* v, int64_t* index1, int64_t* index2, int64_t* index3, double sc) {
			double _NX = __Normal_sigma_X();
			double g1x = 0, g1y = 0, g1z = 0;
			double g2x = 0, g2y = 0, g2z = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				g1x += _ref->d1[0][s] * _ref->buf_xi[s];
				g1y += _ref->d1[0][s] * _ref->buf_eta[s];
				g1z += _ref->d1[0][s] * _ref->buf_nu[s];
				g2x += _ref->d1[1][s] * _ref->buf_xi[s];
				g2y += _ref->d1[1][s] * _ref->buf_eta[s];
				g2z += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g1x = 0, _g1y = 0, _g1z = 0;
				double _g2x = 0, _g2y = 0, _g2z = 0;

				_g1x = _ref->d1[0][s];
				_g1y = _ref->d1[0][s];
				_g1z = _ref->d1[0][s];
				_g2x = _ref->d1[1][s];
				_g2y = _ref->d1[1][s];
				_g2z = _ref->d1[1][s];
				val = _g1y * g2z - g1z * _g2y;
				v->__v(index2[s]) = val;
				val = g1y * _g2z - _g1z * g2y;
				v->__v(index3[s]) = val;
			}
			return _NX * sc;
		}
		double Normal_sigma_Y(_myDoubleArray* v, int64_t* index1, int64_t* index2, int64_t* index3, double sc) {
			double _NY = __Normal_sigma_Y();
			double g1x = 0, g1y = 0, g1z = 0;
			double g2x = 0, g2y = 0, g2z = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				g1x += _ref->d1[0][s] * _ref->buf_xi[s];
				g1y += _ref->d1[0][s] * _ref->buf_eta[s];
				g1z += _ref->d1[0][s] * _ref->buf_nu[s];
				g2x += _ref->d1[1][s] * _ref->buf_xi[s];
				g2y += _ref->d1[1][s] * _ref->buf_eta[s];
				g2z += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g1x = 0, _g1y = 0, _g1z = 0;
				double _g2x = 0, _g2y = 0, _g2z = 0;

				_g1x = _ref->d1[0][s];
				_g1y = _ref->d1[0][s];
				_g1z = _ref->d1[0][s];
				_g2x = _ref->d1[1][s];
				_g2y = _ref->d1[1][s];
				_g2z = _ref->d1[1][s];
				val = _g1z * g2x - g1x * _g2z;
				v->__v(index3[s]) = val;
				val = g1z * _g2x - _g1x * g2z;
				v->__v(index1[s]) = val;
			}
			return _NY * sc;
		}
		double Normal_sigma_Z(_myDoubleArray* v, int64_t* index1, int64_t* index2, int64_t* index3, double sc) {
			double _NZ = __Normal_sigma_Z();
			double g1x = 0, g1y = 0, g1z = 0;
			double g2x = 0, g2y = 0, g2z = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				g1x += _ref->d1[0][s] * _ref->buf_xi[s];
				g1y += _ref->d1[0][s] * _ref->buf_eta[s];
				g1z += _ref->d1[0][s] * _ref->buf_nu[s];
				g2x += _ref->d1[1][s] * _ref->buf_xi[s];
				g2y += _ref->d1[1][s] * _ref->buf_eta[s];
				g2z += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g1x = 0, _g1y = 0, _g1z = 0;
				double _g2x = 0, _g2y = 0, _g2z = 0;

				_g1x = _ref->d1[0][s];
				_g1y = _ref->d1[0][s];
				_g1z = _ref->d1[0][s];
				_g2x = _ref->d1[1][s];
				_g2y = _ref->d1[1][s];
				_g2z = _ref->d1[1][s];
				val = _g1x * g2y - g1y * _g2x;
				v->__v(index1[s]) = val;
				val = g1x * _g2y - _g1y * g2x;
				v->__v(index2[s]) = val;
			}
			return _NZ * sc;
		}*/
		double guideBC(double t1,double t2,bool accurate)
		{
			double val = 0;
			double q1 = 0, q2 = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (get_gij(0, 1) * t1 + get_gij(1, 1) * t2) / dv;
				q2 = (-get_gij(0, 0) * t1 - get_gij(0, 1) * t2) / dv;

			
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (_ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2) / _ref->_refDv;
				q2 = (-_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2) / _ref->_refDv;

			}
			val = get_eij(0, 0) * t1 * q1 + get_eij(0, 1) * t1 * q2 + get_eij(1, 0) * t2 * q1 + get_eij(1, 1) * t2 * q2;
			return val;
		}

		void guideBC_xi(double* ptr, double t1, double t2,bool accurate)
		{
			double val = 0;
			double q1 = 0, q2 = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (get_gij(0, 1) * t1 + get_gij(1, 1) * t2) / dv;
				q2 = (-get_gij(0, 0) * t1 - get_gij(0, 1) * t2) / dv;
			}
			else
			{
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (_ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2) / _ref->_refDv;
				q2 = (-_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2) / _ref->_refDv;
			}
			
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				double _e21 = _e12;

				val = _e11 * t1 * q1 + _e12 * t1 * q2 + _e21 * t2 * q1 +_e22 * t2 * q2;

				*ptr1 = val;
				ptr1++;
			}

		}
		void guideBC_eta(double* ptr, double t1, double t2,bool accurate)
		{
			double val = 0;
			double q1 = 0, q2 = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (get_gij(0, 1) * t1 + get_gij(1, 1) * t2) / dv;
				q2 = (-get_gij(0, 0) * t1 - get_gij(0, 1) * t2) / dv;
			}
			else
			{
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (_ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2) / _ref->_refDv;
				q2 = (-_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2) / _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				double _e21 = _e12;

				val = _e11 * t1 * q1 + _e12 * t1 * q2 + _e21 * t2 * q1 + _e22 * t2 * q2;

				*ptr1 = val;
				ptr1++;
			}

		}
		void guideBC_nu(double* ptr, double t1, double t2,bool accurate)
		{
			double val = 0;
			double q1 = 0, q2 = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (get_gij(0, 1) * t1 + get_gij(1, 1) * t2) / dv;
				q2 = (-get_gij(0, 0) * t1 - get_gij(0, 1) * t2) / dv;
			}
			else
			{
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
				q1 = (_ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2) / _ref->_refDv;
				q2 = (-_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2) / _ref->_refDv;
			}
	
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
				}
				double _e21 = _e12;

				val = _e11 * t1 * q1 + _e12 * t1 * q2 + _e21 * t2 * q1 + _e22 * t2 * q2;


				*ptr1 = val;
				ptr1++;
			}

		}
		double guide_trace( double n1, double n2, double w1, double w2,bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij2(0, 1) + t2 * get_gij2(1, 1);
				u2 = -t1 * get_gij2(0, 0) - t2 * get_gij2(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}
			
			val = (get_eij(0, 0) * t1 * t1 + get_eij(0, 1) * t1 * t2 + get_eij(1, 0) * t2 * t1 + get_eij(1, 1) * t2 * t2);
			val += (get_eij(0, 0) * u1 * u1 + get_eij(0, 1) * u1 * u2 + get_eij(1, 0) * u2 * u1 + get_eij(1, 1) * u2 * u2);
			return val;
		}

		void guide_trace_xi(double* ptr,  double n1, double n2, double w1, double w2,bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij2(0, 1) + t2 * get_gij2(1, 1);
				u2 = -t1 * get_gij2(0, 0) - t2 * get_gij2(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e22 = __dsigma_22[0][s];

				double _e21 = _e12;

				val = /*w1 **/ (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val += /*w2 **/ (_e11 * u1 * u1 + _e12 * u1 * u2 + _e21 * u2 * u1 + _e22 * u2 * u2);

				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_trace_eta(double* ptr, double n1, double n2, double w1, double w2,bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij2(0, 1) + t2 * get_gij2(1, 1);
				u2 = -t1 * get_gij2(0, 0) - t2 * get_gij2(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e22 = __dsigma_22[1][s];

				double _e21 = _e12;

				val = /*w1 **/ (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val += /*w2 **/ (_e11 * u1 * u1 + _e12 * u1 * u2 + _e21 * u2 * u1 + _e22 * u2 * u2);

				*ptr1 = val;
				ptr1++;
			}

		}

		void guide_trace_nu(double* ptr,  double n1, double n2, double w1, double w2,bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij2(0, 1) + t2 * get_gij2(1, 1);
				u2 = -t1 * get_gij2(0, 0) - t2 * get_gij2(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e22 = __dsigma_22[2][s];

				double _e21 = _e12;

				val = /*w1 **/ (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val += /*w2 **/ (_e11 * u1 * u1 + _e12 * u1 * u2 + _e21 * u2 * u1 + _e22 * u2 * u2);

				*ptr1 = val;
				ptr1++;
			}
			
		}
		double __area()
		{
	
			double val = _dv;
			return val;
		}
		void __area_u(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * _ref->d1[0][s] * get_gi2(0, 0);
				double _g12 = _ref->d1[0][s] * get_gi2(1, 0) + _ref->d1[1][s] * get_gi2(0, 0);
				double _g21 = _g12;
				double _g22 = 2 * _ref->d1[1][s] * get_gi2(1, 0);

				double val = 0.5*(_g11*get_Gij(0,0)+2* _g12 * get_Gij(0, 1)+ _g22 * get_Gij(1, 1))* _dv;
				*ptr1 = val;
				ptr1++;
			}
		}
		void __area_v(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * _ref->d1[0][s] * get_gi2(0, 1);
				double _g12 = _ref->d1[0][s] * get_gi2(1, 1) + _ref->d1[1][s] * get_gi2(0, 1);
				double _g21 = _g12;
				double _g22 = 2 * _ref->d1[1][s] * get_gi2(1, 1);

				double val = 0.5 * (_g11 * get_Gij(0, 0) + 2 * _g12 * get_Gij(0, 1) + _g22 * get_Gij(1, 1)) * _dv;
				*ptr1 = val;
				ptr1++;
			}
		}
		
		double det()
		{
			double val = 0;
			val = sqrt((get_eij(0, 0) * get_eij(1, 1)- get_eij(0, 1) * get_eij(0, 1)));
			return val;
		}

		void det_x(double *ptr)
		{
			double val = 0;
			double* ptr1 = ptr;
			double det = sqrt((get_eij(0, 0) * get_eij(1, 1) - get_eij(0, 1) * get_eij(0, 1)));
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e22 = __dsigma_22[0][s];

				double _e21 = _e12;

				double _det=(_e11 * get_eij(1, 1)+ _e22 * get_eij(0, 0) - 2*_e12* get_eij(0, 1));

				*ptr1 = 0.5/det*_det;
				ptr1++;
			}

		}
		void det_y(double* ptr)
		{
			double val = 0;
			double* ptr1 = ptr;
			double det = sqrt((get_eij(0, 0) * get_eij(1, 1) - get_eij(0, 1) * get_eij(0, 1)));
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e22 = __dsigma_22[1][s];

				double _e21 = _e12;

				double _det = (_e11 * get_eij(1, 1) + _e22 * get_eij(0, 0) - 2 * _e12 * get_eij(0, 1));

				*ptr1 = 0.5 / det * _det;
				ptr1++;
			}


		}
		double cont_xi(double w1,double w2)
		{
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				//yu += _ref->d1[0][s] * _ref->buf_eta[s];
				//yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xuu += _ref->__dh[0][s] * _ref->buf_xi[s];
				xuv += _ref->__dh[1][s] * _ref->buf_xi[s];
				xvu += _ref->__dh[2][s] * _ref->buf_xi[s];
				xvv += _ref->__dh[3][s] * _ref->buf_xi[s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];
			}
			double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
			sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
			//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

			sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
			//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

			sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
			//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

			sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;
			
			double val = sxuu * w1*w1 + 2 * sxuv * w1*w2 + sxvv * w2*w2;
			return val;
		}

		void  cont_xi_xi(double *ptr,double w1, double w2)
		{
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{



				double xuu = 0, xuv = 0, xvu = 0, xvv = 0;



				xuu += _ref->__dh[0][s];
				xuv += _ref->__dh[1][s];
				xvu += _ref->__dh[2][s];
				xvv += _ref->__dh[3][s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];

				double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
				sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
				//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

				sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
				//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

				sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
				//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

				sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

				double val = sxuu * w1 * w1 + 2 * sxuv * w1 * w2 + sxvv * w2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double cont_eta(double w1, double w2)
		{
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;

		

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xuu += _ref->__dh[0][s] * _ref->buf_xi[s];
				xuv += _ref->__dh[1][s] * _ref->buf_xi[s];
				xvu += _ref->__dh[2][s] * _ref->buf_xi[s];
				xvv += _ref->__dh[3][s] * _ref->buf_xi[s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];
			}
			double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
			sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
			//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

			sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
			//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

			sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
			//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

			sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

			double val = sxuu * w1 * w1 + 2 * sxuv * w1 * w2 + sxvv * w2 * w2;
			return val;
		}

		void  cont_eta_eta(double* ptr, double w1, double w2)
		{
			double length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{



				double xuu = 0, xuv = 0, xvu = 0, xvv = 0;



				xuu += _ref->__dh[0][s];
				xuv += _ref->__dh[1][s];
				xvu += _ref->__dh[2][s];
				xvv += _ref->__dh[3][s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];

				double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
				sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
				//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

				sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
				//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

				sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
				//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

				sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

				double val = sxuu * w1 * w1 + 2 * sxuv * w1 * w2 + sxvv * w2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}

		double symm_xi(double v1,double v2,double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				//yu += _ref->d1[0][s] * _ref->buf_eta[s];
				//yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xuu += _ref->__dh[0][s] * _ref->buf_xi[s];
				xuv += _ref->__dh[1][s] * _ref->buf_xi[s];
				xvu += _ref->__dh[2][s] * _ref->buf_xi[s];
				xvv += _ref->__dh[3][s] * _ref->buf_xi[s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];
			}
			double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
			sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
			//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

			sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
			//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

			sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
			//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

			sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

			double val = sxuu * v1 * w1 + sxuv * (v1 * w2+v2*w1) + sxvv * v2 * w2;
			return val;
		}

		void  symm_xi_xi(double* ptr, double v1,double v2,double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{



				double xuu = 0, xuv = 0, xvu = 0, xvv = 0;



				xuu += _ref->__dh[0][s];
				xuv += _ref->__dh[1][s];
				xvu += _ref->__dh[2][s];
				xvv += _ref->__dh[3][s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];

				double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
				sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
				//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

				sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
				//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

				sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
				//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

				sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

				double val = sxuu * v1 * w1 + sxuv * (v1 * w2+v2*w1) + sxvv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double symm_eta(double v1,double v2,double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;



			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xuu += _ref->__dh[0][s] * _ref->buf_xi[s];
				xuv += _ref->__dh[1][s] * _ref->buf_xi[s];
				xvu += _ref->__dh[2][s] * _ref->buf_xi[s];
				xvv += _ref->__dh[3][s] * _ref->buf_xi[s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];
			}
			double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
			sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
			//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

			sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
			//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

			sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
			//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

			sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

			double val = sxuu * v1 * w1 + sxuv * (v1 * w2 +v2*w1)+ sxvv * v2 * w2;
			return val;
		}

		void  symm_eta_eta(double* ptr, double v1,double v2,double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{



				double xuu = 0, xuv = 0, xvu = 0, xvv = 0;



				xuu += _ref->__dh[0][s];
				xuv += _ref->__dh[1][s];
				xvu += _ref->__dh[2][s];
				xvv += _ref->__dh[3][s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];

				double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
				sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
				//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

				sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
				//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

				sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
				//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

				sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;

				double val = sxuu * v1 * w1 + sxuv * (v1 * w2+v2*w1) + sxvv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}


		double harmonic_x()
		{

			double xu = 0, xv = 0;// , yu = 0, yv = 0;


			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xuu += (_ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s]) * _ref->buf_xi[s];
				xuv += (_ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s]) * _ref->buf_xi[s];
				xvu += (_ref->d2[2][s] - _ref->oGammaijk[4] * _ref->d1[0][s] - _ref->oGammaijk[5] * _ref->d1[1][s]) * _ref->buf_xi[s];
				xvv += (_ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s]) * _ref->buf_xi[s];

				//yuu += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s]  - _ref->_Gammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yuv += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s]  - _ref->_Gammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvu += (_ref->d2[2][s] - _ref->_Gammaijk[4] * _ref->d1[0][s]  - _ref->_Gammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				//yvv += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s]  - _ref->_Gammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];
			}
			double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
			sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
			//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

			sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
			//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

			sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
			//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

			sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;
			//syvv = yvv - _ref->_Gammaijk[6] * yu - _ref->_Gammaijk[7] * yv;
			double S11 = _ref->oG11;
			double S12 = _ref->oG12;
			double S22 = _ref->oG22;
			double S21 = _ref->oG12;
			double val = sxuu * S11 + 2 * sxuv * S12 + sxvv * S22;
			return val;
		}
		void harmonic_x_xi(double* ptr)
		{
			double S11 = _ref->oG11;
			double S12 = _ref->oG12;
			double S22 = _ref->oG22;
			double S21 = _ref->oG12;


			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;
			double xu = 0, xv = 0;// , yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xuu = (_ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s]);
				xuv = (_ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s]);
				xvu = (_ref->d2[2][s] - _ref->oGammaijk[4] * _ref->d1[0][s] - _ref->oGammaijk[5] * _ref->d1[1][s]);
				xvv = (_ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s]);

				double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
				sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
				//syuu = yuu - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

				sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
				//syuv = yuv - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

				sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
				//syvu = yvu - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

				sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;
				//syvv = yvv - _ref->_Gammaijk[6] * yu - _ref->_Gammaijk[7] * yv;

				double val = sxuu * S11 + 2 * sxuv * S12 + sxvv * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		double harmonic_y()
		{
			//double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double yuu = 0, yuv = 0, yvu = 0, yvv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				//xuu += (_ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s]) * _ref->buf_xi[s];
				//xuv += (_ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s]) * _ref->buf_xi[s];
				//xvu += (_ref->d2[2][s] - _ref->oGammaijk[4] * _ref->d1[0][s] - _ref->oGammaijk[5] * _ref->d1[1][s]) * _ref->buf_xi[s];
				//xvv += (_ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s]) * _ref->buf_xi[s];

				yuu += (_ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s]  - _ref->oGammaijk[1] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				yuv += (_ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s]  - _ref->oGammaijk[3] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				yvu += (_ref->d2[2][s] - _ref->oGammaijk[4] * _ref->d1[0][s]  - _ref->oGammaijk[5] * _ref->d1[1][s] ) * _ref->buf_eta[s];
				yvv += (_ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s]  - _ref->oGammaijk[7] * _ref->d1[1][s] ) * _ref->buf_eta[s];
			}
			double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
			//sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
			syuu = yuu;// -_ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

			//sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
			syuv = yuv;// - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;
			
			//sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
			syvu = yvu;// - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

			//sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;
			syvv = yvv;// -_ref->_Gammaijk[6] * yu - _ref->_Gammaijk[7] * yv;
			double S11 = _ref->oG11;
			double S12 = _ref->oG12;
			double S22 = _ref->oG22;
			double S21 = _ref->oG12;
			double val = syuu * S11 + 2 * syuv * S12 + syvv * S22;
			return val;
		}
		void harmonic_y_eta(double* ptr)
		{
			double S11 = _ref->oG11;
			double S12 = _ref->oG12;
			double S22 = _ref->oG22;
			double S21 = _ref->oG12;


			double yuu = 0, yuv = 0, yvu = 0, yvv = 0;
			double* ptr1 = ptr;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				yuu = (_ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s]);
				yuv = (_ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s]);
				yvu = (_ref->d2[2][s] - _ref->oGammaijk[4] * _ref->d1[0][s] - _ref->oGammaijk[5] * _ref->d1[1][s]);
				yvv = (_ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s]);

				double sxuu = 0, sxuv = 0, sxvu = 0, sxvv = 0, syuu = 0, syuv = 0, syvu = 0, syvv = 0;
				//sxuu = xuu;// -_ref->_Gammaijk[0] * xu - _ref->_Gammaijk[1] * xv;
				syuu = yuu;// - _ref->_Gammaijk[0] * yu - _ref->_Gammaijk[1] * yv;

				//sxuv = xuv;// -_ref->_Gammaijk[2] * xu - _ref->_Gammaijk[3] * xv;
				syuv = yuv;// - _ref->_Gammaijk[2] * yu - _ref->_Gammaijk[3] * yv;

				//sxvu = xvu;// -_ref->_Gammaijk[4] * xu - _ref->_Gammaijk[5] * xv;
				syvu = yvu;// - _ref->_Gammaijk[4] * yu - _ref->_Gammaijk[5] * yv;

				//sxvv = xvv;// -_ref->_Gammaijk[6] * xu - _ref->_Gammaijk[7] * xv;
				syvv = yvv;// -_ref->_Gammaijk[6] * yu - _ref->_Gammaijk[7] * yv;

				double val = syuu * S11 + 2 * syuv * S12 + syvv * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		double ortho3(double s1, double s2, double w1, double w2)
		{
			double length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double u1=0, u2=0, v1=0, v2 = 0;

			for (int s = 0; s < _ref->_nNode; s++) 
			{
				u1 += _ref->d1[0][s] * _ref->buf_xi[s];
				u2 += _ref->d1[1][s] * _ref->buf_xi[s];
				v1 += _ref->d1[0][s] * _ref->buf_eta[s];
				v2 += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double b11 = _ref->_ogi[0] * u1 + _ref->_ogi[1] * v1;
			double b12 = _ref->_ogi[0] * u2 + _ref->_ogi[1] * v2;
			double b21 = _ref->_ogi[3] * u1 + _ref->_ogi[4] * v1;
			double b22 = _ref->_ogi[3] * u2 + _ref->_ogi[4] * v2;
			double B11 = b11, B12 = 0.5 * (b12 + b21), B22 = b22;
			double val = B11 * s1 * w1 + B12 * (s1 * w2+s2 * w1) + B22 * s2 * w2;
			return val;
		}
		void ortho3_xi(double* ptr, double s1, double s2, double w1, double w2)
		{
			double length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * (w1 * w2) * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = _ref->_ogi[0] * _ref->d1[0][s];
				double _b12 = _ref->_ogi[0] * _ref->d1[1][s];
				double _b21 = _ref->_ogi[3] * _ref->d1[0][s];
				double _b22 = _ref->_ogi[3] * _ref->d1[1][s];
				double _B11 = _b11, _B12 = 0.5 * (_b12 + _b21), _B22 = _b22;
				double val = _B11 * s1 * w1 + _B12 * (s1 * w2 + s2 * w1) + _B22 * s2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void ortho3_eta(double* ptr, double s1, double s2, double w1, double w2)
		{
			double length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = _ref->_ogi[1] * _ref->d1[0][s];
				double _b12 = _ref->_ogi[1] * _ref->d1[1][s];
				double _b21 = _ref->_ogi[4] * _ref->d1[0][s];
				double _b22 = _ref->_ogi[4] * _ref->d1[1][s];

				double _B11 = _b11, _B12 = 0.5 * (_b12 + _b21), _B22 = _b22;
				double val = _B11 * s1 * w1 + _B12 * (s1 * w2 + s2 * w1) + _B22 * s2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double ortho4(double s1, double s2, double w1, double w2)
		{
			double length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
			double u1 = 0, u2 = 0, v1 = 0, v2 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_xi[s];
				u2 += _ref->d1[1][s] * _ref->buf_xi[s];
				v1 += _ref->d1[0][s] * _ref->buf_eta[s];
				v2 += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double b11 = _ref->_ogi[0] * u1 + _ref->_ogi[1] * v1;
			double b12 = _ref->_ogi[0] * u2 + _ref->_ogi[1] * v2;
			double b21 = _ref->_ogi[3] * u1 + _ref->_ogi[4] * v1;
			double b22 = _ref->_ogi[3] * u2 + _ref->_ogi[4] * v2;
			double B11 = b11, B12 = 0.5 * (b12 + b21), B22 = b22;
			double val = B11 * (w1 * w1) + 2 * B12 * (w1 * w2) + B22 * w2 * w2;
			val -= B11 * (s1 * s1) + 2 * B12 * (s1 * s2) + B22 * s2 * s2;

			return val;
		}
		void ortho4_xi(double* ptr, double s1, double s2, double w1, double w2)
		{
			double length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * (w1 * w2) * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = _ref->_ogi[0] * _ref->d1[0][s];
				double _b12 = _ref->_ogi[0] * _ref->d1[1][s];
				double _b21 = _ref->_ogi[3] * _ref->d1[0][s];
				double _b22 = _ref->_ogi[3] * _ref->d1[1][s];
				double _B11 = _b11, _B12 = 0.5 * (_b12 + _b21), _B22 = _b22;
				double val = _B11 * (w1 * w1) + 2 * _B12 * (w1 * w2) + _B22 * w2 * w2;
				val -= _B11 * (s1 * s1) + 2 * _B12 * (s1 * s2) + _B22 * s2 * s2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void ortho4_eta(double* ptr, double s1, double s2, double w1, double w2)
		{
			double length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
			s1 /= length;
			s2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = _ref->_ogi[1] * _ref->d1[0][s];
				double _b12 = _ref->_ogi[1] * _ref->d1[1][s];
				double _b21 = _ref->_ogi[4] * _ref->d1[0][s];
				double _b22 = _ref->_ogi[4] * _ref->d1[1][s];
				double _B11 = _b11, _B12 = 0.5 * (_b12 + _b21), _B22 = _b22;
				double val = _B11 * (w1 * w1) + 2 * _B12 * (w1 * w2) + _B22 * w2 * w2;
				val -= _B11 * (s1 * s1) + 2 * _B12 * (s1 * s2) + _B22 * s2 * s2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double ortho(double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double b11 = _ref->_ogi[0] * _ref->get__gi(0, 0) + _ref->_ogi[1] * _ref->get__gi(0, 1);
			double b12 = _ref->_ogi[0] * _ref->get__gi(1, 0) + _ref->_ogi[1] * _ref->get__gi(1, 1);
			double b21 = _ref->_ogi[3] * _ref->get__gi(0, 0) + _ref->_ogi[4] * _ref->get__gi(0, 1);
			double b22 = _ref->_ogi[3] * _ref->get__gi(1, 0) + _ref->_ogi[4] * _ref->get__gi(1, 1);

			double B11 = 2 * b11;
			double B12 = b12 + b21;
			double B21 = b12 + b21;
			double B22 = 2*b22;

			double val=B11 * v1 * w1 + B12 * v1*w2+B21*v2 * w1 + B22 * v2 * w2;
			return val;
		}
		void ortho_x(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * (w1  * w2) * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = _ref->_ogi[0] * _ref->d1[0][s];
				double _b12 = _ref->_ogi[0] * _ref->d1[1][s];
				double _b21 = _ref->_ogi[3] * _ref->d1[0][s];
				double _b22 = _ref->_ogi[3] * _ref->d1[1][s];

				double _B11 = 2 * _b11;
				double _B12 = _b12 + _b21;
				double _B21 = _b12 + _b21;
				double _B22 = 2 * _b22;

				double val = _B11 * v1 * w1 + _B12 * v1 * w2 + _B21 * v2 * w1 + _B22 * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void ortho_y(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = _ref->_ogi[1] * _ref->d1[0][s];
				double _b12 = _ref->_ogi[1] * _ref->d1[1][s];
				double _b21 = _ref->_ogi[4] * _ref->d1[0][s];
				double _b22 = _ref->_ogi[4] * _ref->d1[1][s];
				double _B11 = 2 * _b11;
				double _B12 = _b12 + _b21;
				double _B21 = _b12 + _b21;
				double _B22 = 2 * _b22;

				double val = _B11 * v1 * w1 + _B12 * v1 * w2 + _B21 * v2 * w1 + _B22 * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double ortho2(double v1, double v2, double w1, double w2)
		{		
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double b11 = 0, b12 = 0, b21 = 0, b22 = 0;
			double g1x = 0, g1y = 0, g2x = 0, g2y = 0;
			for (int t = 0; t < _ref->_nNode; t++)
			{
				g1x += _ref->d1[0][t] * _ref->node[t * 3 + 0];
				g1y += _ref->d1[0][t] * _ref->node[t * 3 + 1];
				g2x += _ref->d1[1][t] * _ref->node[t * 3 + 0];
				g2y += _ref->d1[1][t] * _ref->node[t * 3 + 1];

			}
			b11 = g1x *_ref->_ogi[0] + g1y * _ref->_ogi[1];
			b12 = g1x * _ref->_ogi[3] + g1y * _ref->_ogi[4];
			b21 = g2x * _ref->_ogi[0] + g2y * _ref->_ogi[1];
			b22 = g2x * _ref->_ogi[3] + g2y * _ref->_ogi[4];
			double B11 = 2 * b11;
			double B12 = b12 + b21;
			double B21 = B12;
			double B22 = 2 * b22;
	
			double val = B11 * v1 * w1 + B12 * v1 * w2 + B21 * v2 * w1 + B22 * v2 * w2;
			return val;
		}
		void ortho2_x(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;

			double val = 0;
			double* ptr1 = ptr;
		
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = 2*_ref->d1[0][s] * _ref->_ogi[0];// +g1y * _ref->_ogi[1];
				double _b12 = _ref->d1[0][s] * _ref->_ogi[3];
				double _b21 = _ref->d1[1][s] * _ref->_ogi[0];
				double _b22 = 2*_ref->d1[1][s] * _ref->_ogi[3];// + g2y * _ref->_ogi[3];
				double _B11 = 2 * _b11;
				double _B12 = _b12 + _b21;
				double _B21 = _B12;
				double _B22 = 2 * _b22;

				double val = _B11 * v1 * w1 + _B12 * v1 * w2 + _B21 * v2 * w1 + _B22 * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void ortho2_y(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->og11 + 2 * v1 * v2 * _ref->og12 + v2 * v2 * _ref->og22);
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->og11 + 2 * w1 * w2 * _ref->og12 + w2 * w2 * _ref->og22);
			w1 /= length;
			w2 /= length;
	

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _b11 = 2 * _ref->d1[0][s] * _ref->_ogi[1];// +g1y * _ref->_ogi[1];
				double _b12 = _ref->d1[0][s] * _ref->_ogi[4];
				double _b21 = _ref->d1[1][s] * _ref->_ogi[1];
				double _b22 = 2 * _ref->d1[1][s] * _ref->_ogi[4];// + g2y * _ref->_ogi[3];
				double _B11 = 2 * _b11;
				double _B12 = _b12 + _b21;
				double _B21 = _B12;
				double _B22 = 2 * _b22;

				double val = _B11 * v1 * w1 + _B12 * v1 * w2 + _B21 * v2 * w1 + _B22 * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		
		

		double curl_free()
		{
			double u1 = 0, u2 = 0, v1 = 0, v2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				u2 += _ref->d1[1][s] * _ref->buf_u[s];
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double scale = 1.0/sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);
			
			double b12 = u1 * _ref->_ogi[3] + v1 * _ref->_ogi[4];
			double b21 = u2 * _ref->_ogi[0] + v2 * _ref->_ogi[1];

			double val = (b12 - b21)*scale;
			return val;
		}
		void curl_free_x(double* ptr)
		{


			double* ptr1 = ptr;

			double scale = 1.0 / sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
				double _u1 = _ref->d1[0][s];// *_ref->buf_u[s];
				double _v1 = 0;// _ref->d1[0][s];// *_ref->buf_v[s];
				double _u2 = _ref->d1[1][s];// *_ref->buf_u[s];
				double _v2 = 0;//_ref->d1[1][s];// * _ref->buf_v[s];

				double _b12 = _u1 * _ref->_ogi[3] + _v1 * _ref->_ogi[4] ;
				double _b21 = _u2 * _ref->_ogi[0] + _v2 * _ref->_ogi[1] ;
				double val = (_b12 - _b21)*scale;
				*ptr1 = val;
				ptr1++;
			}
		}
		void curl_free_y(double* ptr)
		{

			double* ptr1 = ptr;

			double scale = 1.0 / sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				double _u1 = 0;// _ref->d1[0][s];// *_ref->buf_u[s];
				double _v1 = _ref->d1[0][s];// *_ref->buf_v[s];
				double _u2 = 0;// _ref->d1[1][s];// *_ref->buf_u[s];
				double _v2 = _ref->d1[1][s];// * _ref->buf_v[s];

				double _b12 = _u1 * _ref->_ogi[3] + _v1 * _ref->_ogi[4] ;
				double _b21 = _u2 * _ref->_ogi[0] + _v2 * _ref->_ogi[1] ;

				double val = (_b12 - _b21)*scale;
				*ptr1 = val;
				ptr1++;
			}
		}

		double curl_free2()
		{
			double u1 = 0, u2 = 0, v1 = 0, v2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_xi[s];
				v1 += _ref->d1[0][s] * _ref->buf_eta[s];
				u2 += _ref->d1[1][s] * _ref->buf_xi[s];
				v2 += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double scale = 1.0 / sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);

			double b12 = u1 * _ref->_ogi[3] + v1 * _ref->_ogi[4];
			double b21 = u2 * _ref->_ogi[0] + v2 * _ref->_ogi[1];

			double val = (b12 - b21) * scale;
			return val;
		}
		void curl_free2_xi(double* ptr)
		{


			double* ptr1 = ptr;

			double scale = 1.0 / sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
				double _u1 = _ref->d1[0][s];// *_ref->buf_u[s];
				double _v1 = 0;// _ref->d1[0][s];// *_ref->buf_v[s];
				double _u2 = _ref->d1[1][s];// *_ref->buf_u[s];
				double _v2 = 0;//_ref->d1[1][s];// * _ref->buf_v[s];

				double _b12 = _u1 * _ref->_ogi[3] + _v1 * _ref->_ogi[4];
				double _b21 = _u2 * _ref->_ogi[0] + _v2 * _ref->_ogi[1];
				double val = (_b12 - _b21) * scale;
				*ptr1 = val;
				ptr1++;
			}
		}
		void curl_free2_eta(double* ptr)
		{

			double* ptr1 = ptr;

			double scale = 1.0 / sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				double _u1 = 0;// _ref->d1[0][s];// *_ref->buf_u[s];
				double _v1 = _ref->d1[0][s];// *_ref->buf_v[s];
				double _u2 = 0;// _ref->d1[1][s];// *_ref->buf_u[s];
				double _v2 = _ref->d1[1][s];// * _ref->buf_v[s];

				double _b12 = _u1 * _ref->_ogi[3] + _v1 * _ref->_ogi[4];
				double _b21 = _u2 * _ref->_ogi[0] + _v2 * _ref->_ogi[1];

				double val = (_b12 - _b21) * scale;
				*ptr1 = val;
				ptr1++;
			}
		}



		double div_free(double v1,double v2)
		{
			double V1 = v1 * _ref->get__gij(0, 0) + v2 * _ref->get__gij(1, 0);
			double V2 = v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
			double length = sqrt(V1 * V1 * _ref->get__Gij(0, 0) + 2 * V1 * V2 * _ref->get__Gij(0, 1) + V2 * V2 * _ref->get__Gij(1, 1));
			V1 /= length;
			V2 /= length;
			double suuu = 0, svuv = 0;
			double suvu = 0, svvv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				suuu += _ref->d1[0][s] * _ref->buf_xi[s] +
					2 * (_ref->get__Gammaijk(0, 0, 0) * _ref->d0[s] * _ref->buf_xi[s] + _ref->get__Gammaijk(0, 0, 1) * _ref->d0[s] * _ref->buf_eta[s]);
				//suu, u + 2*(G_{ u,u }^ {1}s{1u} - G_{ u,u }^ {2}s{ 2u })

				svuv += _ref->d1[1][s] * _ref->buf_eta[s] +
					_ref->get__Gammaijk(0, 1, 0) * _ref->d0[s] * _ref->buf_eta[s] + _ref->get__Gammaijk(0, 1, 1) * _ref->d0[s] * _ref->buf_phi[s] +
					_ref->get__Gammaijk(1, 1, 0) * _ref->d0[s] * _ref->buf_xi[s] + _ref->get__Gammaijk(1, 1, 1) * _ref->d0[s] * _ref->buf_eta[s];
				//suv, v + G_{ u,v }^ {1}s{1v} + G_{ u,v }^ {2}s{ 2v }+ G_{ v,v }^ {1}s{u,1} + G_{ v,v }^ {2}s{ u,2 }

				suvu += _ref->d1[0][s] * _ref->buf_eta[s] +
					_ref->get__Gammaijk(1, 0, 0) * _ref->d0[s] * _ref->buf_xi[s] + _ref->get__Gammaijk(1, 0, 1) * _ref->d0[s] * _ref->buf_eta[s] +
					_ref->get__Gammaijk(0, 0, 0) * _ref->d0[s] * _ref->buf_eta[s] + _ref->get__Gammaijk(0, 0, 1) * _ref->d0[s] * _ref->buf_phi[s];
				//svu, u + G_{ v,u }^ {1}s{1u} + G_{ v,u }^ {2}s{ 2u }+ G_{ u,u }^ {1}s{v,1} + G_{ u,u }^ {2}s{ v,2 }

				svvv += _ref->d1[1][s] * _ref->buf_phi[s] +
					2 * (_ref->get__Gammaijk(1, 1, 0) * _ref->d0[s] * _ref->buf_eta[s] + _ref->get__Gammaijk(1, 1, 1) * _ref->d0[s] * _ref->buf_phi[s]);
				//svv, v + 2*(G_{ v,v }^ {1}s{1v} + G_{ v,v }^ {2}s{ 2v })

			}

		
			double Fu = suuu + svuv;
			double Fv = suvu + svvv;
			return Fu*V1+Fv*V2;
		}
		void div_free_xi(double* ptr, double v1, double v2)
		{

			double V1 = v1 * _ref->get__gij(0, 0) + v2 * _ref->get__gij(1, 0);
			double V2 = v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
			double length = sqrt(V1 * V1 * _ref->get__Gij(0, 0) + 2 * V1 * V2 * _ref->get__Gij(0, 1) + V2 * V2 * _ref->get__Gij(1, 1));
			V1 /= length;
			V2 /= length;

			double* ptr1 = ptr;

			
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double suuu = _ref->d1[0][s]  +
					2 * (_ref->get__Gammaijk(0, 0, 0) * _ref->d0[s]);
				//suu, u + 2*(G_{ u,u }^ {1}s{1u} - G_{ u,u }^ {2}s{ 2u })

				double svuv = 
				
					_ref->get__Gammaijk(1, 1, 0) * _ref->d0[s] ;
				//suv, v + G_{ u,v }^ {1}s{1v} + G_{ u,v }^ {2}s{ 2v }+ G_{ v,v }^ {1}s{u,1} + G_{ v,v }^ {2}s{ u,2 }

				double suvu =
					_ref->get__Gammaijk(1, 0, 0) * _ref->d0[s] ;
				//svu, u + G_{ v,u }^ {1}s{1u} + G_{ v,u }^ {2}s{ 2u }+ G_{ u,u }^ {1}s{v,1} + G_{ u,u }^ {2}s{ v,2 }

				double svvv = 0;

				double Fu = suuu + svuv;
				double Fv = suvu + svvv;
				double val=Fu * V1 + Fv * V2;

				*ptr1 = val;
				ptr1++;
			}
		}
		void div_free_eta(double* ptr, double v1, double v2)
		{
			double V1 = v1 * _ref->get__gij(0, 0) + v2 * _ref->get__gij(1, 0);
			double V2 = v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
			double length = sqrt(V1 * V1 * _ref->get__Gij(0, 0) + 2 * V1 * V2 * _ref->get__Gij(0, 1) + V2 * V2 * _ref->get__Gij(1, 1));
			V1 /= length;
			V2 /= length;
			double* ptr1 = ptr;

			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double suuu =
					2 * (_ref->get__Gammaijk(0, 0, 1) * _ref->d0[s]);
				//suu, u + 2*(G_{ u,u }^ {1}s{1u} - G_{ u,u }^ {2}s{ 2u })

				double svuv = _ref->d1[1][s]  +
					_ref->get__Gammaijk(0, 1, 0) * _ref->d0[s] +
					 _ref->get__Gammaijk(1, 1, 1) * _ref->d0[s];
				//suv, v + G_{ u,v }^ {1}s{1v} + G_{ u,v }^ {2}s{ 2v }+ G_{ v,v }^ {1}s{u,1} + G_{ v,v }^ {2}s{ u,2 }

				double suvu = _ref->d1[0][s]  +
					_ref->get__Gammaijk(1, 0, 1) * _ref->d0[s]  +
					_ref->get__Gammaijk(0, 0, 0) * _ref->d0[s] ;
				//svu, u + G_{ v,u }^ {1}s{1u} + G_{ v,u }^ {2}s{ 2u }+ G_{ u,u }^ {1}s{v,1} + G_{ u,u }^ {2}s{ v,2 }

				double svvv = 
					2 * (_ref->get__Gammaijk(1, 1, 0) * _ref->d0[s] );

				double Fu = suuu + svuv;
				double Fv = suvu + svvv;
				double val = Fu * V1 + Fv * V2;

				*ptr1 = val;
				ptr1++;
			}
		}
		void div_free_phi(double* ptr, double v1, double v2)
		{
			double V1 = v1 * _ref->get__gij(0, 0) + v2 * _ref->get__gij(1, 0);
			double V2 = v1 * _ref->get__gij(0, 1) + v2 * _ref->get__gij(1, 1);
			double length = sqrt(V1 * V1 * _ref->get__Gij(0, 0) + 2 * V1 * V2 * _ref->get__Gij(0, 1) + V2 * V2 * _ref->get__Gij(1, 1));
			V1 /= length;
			V2 /= length;
			double* ptr1 = ptr;

			double scale = 1.0 / sqrt(_ref->og11 * _ref->og22 - _ref->og12 * _ref->og12);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double suuu = 0;
				//suu, u + 2*(G_{ u,u }^ {1}s{1u} - G_{ u,u }^ {2}s{ 2u })

				double svuv =
					_ref->get__Gammaijk(0, 1, 1) * _ref->d0[s];
					
				//suv, v + G_{ u,v }^ {1}s{1v} + G_{ u,v }^ {2}s{ 2v }+ G_{ v,v }^ {1}s{u,1} + G_{ v,v }^ {2}s{ u,2 }

				double suvu = 
					
					 _ref->get__Gammaijk(0, 0, 1)* _ref->d0[s];
				//svu, u + G_{ v,u }^ {1}s{1u} + G_{ v,u }^ {2}s{ 2u }+ G_{ u,u }^ {1}s{v,1} + G_{ u,u }^ {2}s{ v,2 }

				double svvv = _ref->d1[1][s]  +
					2 * (_ref->get__Gammaijk(1, 1, 1) * _ref->d0[s]);

				double Fu = suuu + svuv;
				double Fv = suvu + svvv;
				double val = Fu * V1 + Fv * V2;

				*ptr1 = val;
				ptr1++;
			}
		}
		

		double harmonic_X()
		{
			


			double Xuu = 0, Xuv = 0, Xvu = 0, Xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				Xuu += _ref->__dh[0][s] * (_ref->buf_u[s]);
				Xuv += _ref->__dh[1][s] * (_ref->buf_u[s]);
				Xvu += _ref->__dh[2][s] * (_ref->buf_u[s]);
				Xvv += _ref->__dh[3][s] * (_ref->buf_u[s]);

			}
		
			double val = Xuu * _ref->get__Gij(0,0) + 2 * Xuv * _ref->get__Gij(0, 1) + Xvv * _ref->get__Gij(1, 1);
			return val;
		}


		void harmonic_X_u(double* ptr)
		{
			double Xuu = 0, Xuv = 0, Xvu = 0, Xvv = 0;
			
			double Xu = 0, Xv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Xu += _ref->d1[0][s] * _ref->buf_u[s];
				Xv += _ref->d1[1][s] * _ref->buf_u[s];
				Xuu += _ref->__dh[0][s] * (_ref->buf_u[s]);
				Xuv += _ref->__dh[1][s] * (_ref->buf_u[s]);
				Xvu += _ref->__dh[2][s] * (_ref->buf_u[s]);
				Xvv += _ref->__dh[3][s] * (_ref->buf_u[s]);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1,0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1,0);
				_g21 = _g12;

				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;

				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Xuu = 0, _Xuv = 0, _Xvu = 0, _Xvv = 0;
				_Xuu = _ref->__dh[0][s] - _Gamma111 * Xu - _Gamma112 * Xv;
				_Xuv = _ref->__dh[1][s] - _Gamma121 * Xu - _Gamma122 * Xv;
				_Xvv = _ref->__dh[3][s] - _Gamma221 * Xu - _Gamma222 * Xv;
			
				double val = _Xuu * _ref->get__Gij(0,0) + 2 * _Xuv * _ref->get__Gij(0, 1) + _Xvv * _ref->get__Gij(1, 1);
				val += Xuu * _G11 + 2 * Xuv * _G12 + Xvv * _G22;

				*ptr1 = val;
				ptr1++;
			}
		}

		void harmonic_X_v(double* ptr)
		{
			double Xuu = 0, Xuv = 0, Xvu = 0, Xvv = 0;

			double Xu = 0, Xv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Xu += _ref->d1[0][s] * _ref->buf_u[s];
				Xv += _ref->d1[1][s] * _ref->buf_u[s];
				Xuu += _ref->__dh[0][s] * (_ref->buf_u[s]);
				Xuv += _ref->__dh[1][s] * (_ref->buf_u[s]);
				Xvu += _ref->__dh[2][s] * (_ref->buf_u[s]);
				Xvv += _ref->__dh[3][s] * (_ref->buf_u[s]);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);
				_g21 = _g12;

				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;

				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Xuu = 0, _Xuv = 0, _Xvu = 0, _Xvv = 0;
				_Xuu = - _Gamma111 * Xu - _Gamma112 * Xv;
				_Xuv =  - _Gamma121 * Xu - _Gamma122 * Xv;
				_Xvv = - _Gamma221 * Xu - _Gamma222 * Xv;

				double val = _Xuu * _ref->get__Gij(0, 0) + 2 * _Xuv * _ref->get__Gij(0, 1) + _Xvv * _ref->get__Gij(1, 1);
				val += Xuu * _G11 + 2 * Xuv * _G12 + Xvv * _G22;

				*ptr1 = val;
				ptr1++;
			}
		}
		double harmonic_Y()
		{



			double Xuu = 0, Xuv = 0, Xvu = 0, Xvv = 0;
			//double yuu = 0, yuv = 0, yvu = 0, yvv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				Xuu += _ref->__dh[0][s] * (_ref->buf_v[s]);
				Xuv += _ref->__dh[1][s] * (_ref->buf_v[s]);
				Xvu += _ref->__dh[2][s] * (_ref->buf_v[s]);
				Xvv += _ref->__dh[3][s] * (_ref->buf_v[s]);

			}

			double val = Xuu * _ref->get__Gij(0, 0) + 2 * Xuv * _ref->get__Gij(0, 1) + Xvv * _ref->get__Gij(1, 1);
			return val;
		}


		void harmonic_Y_u(double* ptr)
		{
			double Xuu = 0, Xuv = 0, Xvu = 0, Xvv = 0;

			double Xu = 0, Xv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Xu += _ref->d1[0][s] * _ref->buf_v[s];
				Xv += _ref->d1[1][s] * _ref->buf_v[s];
				Xuu += _ref->__dh[0][s] * (_ref->buf_v[s]);
				Xuv += _ref->__dh[1][s] * (_ref->buf_v[s]);
				Xvu += _ref->__dh[2][s] * (_ref->buf_v[s]);
				Xvv += _ref->__dh[3][s] * (_ref->buf_v[s]);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);
				_g21 = _g12;

				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;

				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Xuu = 0, _Xuv = 0, _Xvu = 0, _Xvv = 0;
				_Xuu =  - _Gamma111 * Xu - _Gamma112 * Xv;
				_Xuv =  - _Gamma121 * Xu - _Gamma122 * Xv;
				_Xvv =  - _Gamma221 * Xu - _Gamma222 * Xv;

				double val = _Xuu * _ref->get__Gij(0, 0) + 2 * _Xuv * _ref->get__Gij(0, 1) + _Xvv * _ref->get__Gij(1, 1);
				val += Xuu * _G11 + 2 * Xuv * _G12 + Xvv * _G22;

				*ptr1 = val;
				ptr1++;
			}
		}

		void harmonic_Y_v(double* ptr)
		{
			double Xuu = 0, Xuv = 0, Xvu = 0, Xvv = 0;

			double Xu = 0, Xv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Xu += _ref->d1[0][s] * _ref->buf_v[s];
				Xv += _ref->d1[1][s] * _ref->buf_v[s];
				Xuu += _ref->__dh[0][s] * (_ref->buf_v[s]);
				Xuv += _ref->__dh[1][s] * (_ref->buf_v[s]);
				Xvu += _ref->__dh[2][s] * (_ref->buf_v[s]);
				Xvv += _ref->__dh[3][s] * (_ref->buf_v[s]);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);
				_g21 = _g12;

				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;

				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Xuu = 0, _Xuv = 0, _Xvu = 0, _Xvv = 0;
				_Xuu = _ref->__dh[0][s] -_Gamma111 * Xu - _Gamma112 * Xv;
				_Xuv = _ref->__dh[1][s] -_Gamma121 * Xu - _Gamma122 * Xv;
				_Xvv = _ref->__dh[3][s] -_Gamma221 * Xu - _Gamma222 * Xv;

				double val = _Xuu * _ref->get__Gij(0, 0) + 2 * _Xuv * _ref->get__Gij(0, 1) + _Xvv * _ref->get__Gij(1, 1);
				val += Xuu * _G11 + 2 * Xuv * _G12 + Xvv * _G22;

				*ptr1 = val;
				ptr1++;
			}
		}

		double rot_free()
		{


			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
			double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
			double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
			double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

			return xy - yx;
		}



		void rot_free_u(double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				xu = _ref->d1[0][s];
				xv = _ref->d1[1][s];
				//yu = _ref->d1[0][s];
				//yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xy - yx;

				*ptr1 = val;
				ptr1++;
			}
		}
		void rot_free_v(double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				//xu = _ref->d1[0][s];
				//xv = _ref->d1[1][s];
				yu = _ref->d1[0][s];
				yv = _ref->d1[1][s];
				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xy - yx;

				*ptr1 = val;
				ptr1++;
			}
		}
		double conformal_x()
		{


			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
			double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
			double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
			double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

			return xx - yy;
		}



		void conformal_x_xi(double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				xu = _ref->d1[0][s];
				xv = _ref->d1[1][s];
				//yu = _ref->d1[0][s];
				//yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xx - yy;

				*ptr1 = val;
				ptr1++;
			}
		}
		void conformal_x_eta( double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				//xu = _ref->d1[0][s];
				//xv = _ref->d1[1][s];
				yu = _ref->d1[0][s];
				yv = _ref->d1[1][s];
				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xx - yy;

				*ptr1 = val;
				ptr1++;
			}
		}
		double conformal_y()
		{


			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
			double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
			double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
			double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

			return xy + yx;
		}



		void conformal_y_xi( double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				xu = _ref->d1[0][s];
				xv = _ref->d1[1][s];
				//yu = _ref->d1[0][s];
				//yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xy + yx;

				*ptr1 = val;
				ptr1++;
			}
		}
		void conformal_y_eta( double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				//xu = _ref->d1[0][s];
				//xv = _ref->d1[1][s];
				yu = _ref->d1[0][s];
				yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xy + yx;

				*ptr1 = val;
				ptr1++;
			}
		}

		double conformal_X()
		{


			double xu = 0, xv = 0 , yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
			double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
			double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
			double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];
		
			return xx-yy;
		}



		void conformal_X_x( double* ptr)
		{
			
			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;
		

			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;

				
					xu = _ref->d1[0][s];
					xv = _ref->d1[1][s];
					//yu = _ref->d1[0][s];
					//yv = _ref->d1[1][s];
					double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
					double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
					double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
					double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val=xx - yy;

				*ptr1 = val;
				ptr1++;
			}
		}
		void conformal_X_y( double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				//xu = _ref->d1[0][s];
				//xv = _ref->d1[1][s];
				yu = _ref->d1[0][s];
				yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xx - yy;

				*ptr1 = val;
				ptr1++;
			}
		}
		double conformal_Y()
		{


			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
			double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
			double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
			double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

			return xy + yx;
		}



		void conformal_Y_x( double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				xu = _ref->d1[0][s];
				xv = _ref->d1[1][s];
				//yu = _ref->d1[0][s];
				//yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xy + yx;

				*ptr1 = val;
				ptr1++;
			}
		}
		void conformal_Y_y( double* ptr)
		{

			double xuu = 0, xuv = 0, xvu = 0, xvv = 0;
			double* ptr1 = ptr;


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;


				//xu = _ref->d1[0][s];
				//xv = _ref->d1[1][s];
				yu = _ref->d1[0][s];
				yv = _ref->d1[1][s];

				double xx = xu * _ref->_oGi[0] + xv * _ref->_oGi[3];
				double xy = xu * _ref->_oGi[1] + xv * _ref->_oGi[4];
				double yx = yu * _ref->_oGi[0] + yv * _ref->_oGi[3];
				double yy = yu * _ref->_oGi[1] + yv * _ref->_oGi[4];

				double val = xy + yx;

				*ptr1 = val;
				ptr1++;
			}
		}
		
		double guide_trace2()
		{
		
			double val = 0;
		

			val = get_eij(0,0) * _ref->oG11 + get_eij(0, 1) * _ref->oG12 + get_eij(1,0) * _ref->oG12 + get_eij(1, 1) * _ref->oG22;

		

			return val;
		}

		void guide_trace2_xi(double* ptr)
		{
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e22 = __dsigma_22[0][s];
					val = _e11 * _ref->oG11 + _e12 * _ref->oG12 + _e12 * _ref->oG12 + _e22 * _ref->oG22;


				
					*ptr1 = val;
				ptr1++;
			}

		}
		void guide_trace2_eta( double* ptr)
		{
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e22 = __dsigma_22[1][s];
				val = _e11 * _ref->oG11 + _e12 * _ref->oG12 + _e12 * _ref->oG12 + _e22 * _ref->oG22;



				*ptr1 = val;
				ptr1++;
			}

		}
	

		double guide_traceBC(double n1, double n2, double w1, double w2, bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij(0, 1) + t2 * get_gij(1, 1);
				u2 = -t1 * get_gij(0, 0) - t2 * get_gij(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			val = w1 * (get_eij(0, 0) * t1 * t1 + get_eij(0, 1) * t1 * t2 + get_eij(1, 0) * t2 * t1 + get_eij(1, 1) * t2 * t2);
			val += w2 * (get_eij(0, 0) * u1 * u1 + get_eij(0, 1) * u1 * u2 + get_eij(1, 0) * u2 * u1 + get_eij(1, 1) * u2 * u2);
			return val;
		}

		void guide_traceBC_xi(double* ptr, double n1, double n2, double w1, double w2, bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij(0, 1) + t2 * get_gij(1, 1);
				u2 = -t1 * get_gij(0, 0) - t2 * get_gij(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else
			 {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}
				double _e21 = _e12;

				val = w1 * (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val += w2 * (_e11 * u1 * u1 + _e12 * u1 * u2 + _e21 * u2 * u1 + _e22 * u2 * u2);

				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_traceBC_eta(double* ptr, double n1, double n2, double w1, double w2, bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij(0, 1) + t2 * get_gij(1, 1);
				u2 = -t1 * get_gij(0, 0) - t2 * get_gij(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}
				double _e21 = _e12;

				val = w1 * (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val += w2 * (_e11 * u1 * u1 + _e12 * u1 * u2 + _e21 * u2 * u1 + _e22 * u2 * u2);
				//val +=1./w * (_e11 * t1 * u1 + _e12 * t1 * u2 + _e21 * t2 * u1 + _e22 * t2 * u2);
				//val += 1./w * (_e11 * u1 * t1 + _e12 * u1 * t2 + _e21 * u2 * t1 + _e22 * u2 * t2);
				*ptr1 = val;
				ptr1++;
			}

		}

		void guide_traceBC_nu(double* ptr, double n1, double n2, double w1, double w2, bool accurate)
		{
			double val = 0;
			double t1 = n1, t2 = n2;
			double u1 = 0, u2 = 0;
			if (accurate) {
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * get_gij(0, 1) + t2 * get_gij(1, 1);
				u2 = -t1 * get_gij(0, 0) - t2 * get_gij(1, 0);

				u1 /= dv;
				u2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				u1 = t1 * _ref->get__gij(0, 1) + t2 * _ref->get__gij(1, 1);
				u2 = -t1 * _ref->get__gij(0, 0) - t2 * _ref->get__gij(1, 0);

				u1 /= _ref->_refDv;
				u2 /= _ref->_refDv;
			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
				}
				double _e21 = _e12;

				val = w1 * (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val += w2 * (_e11 * u1 * u1 + _e12 * u1 * u2 + _e21 * u2 * u1 + _e22 * u2 * u2);
				//val += 1./w * (_e11 * t1 * u1 + _e12 * t1 * u2 + _e21 * t2 * u1 + _e22 * t2 * u2);
				//val += 1./w * (_e11 * u1 * t1 + _e12 * u1 * t2 + _e21 * u2 * t1 + _e22 * u2 * t2);
				*ptr1 = val;
				ptr1++;
			}

		}
		double guideA(double t1, double t2, double s1, double s2)
		{
			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;

			double u1 = 0, u2 = 0, v1 = 0, v2 = 0 , w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				w1 += _ref->d1[0][s] * _ref->buf_w[s];
				u2 += _ref->d1[1][s] * _ref->buf_u[s];
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
			    w2 += _ref->d1[1][s] * _ref->buf_w[s];
			}
			double e11 = u1 * u1 + v1 * v1 +w1 * w1;
			double e12 = u1 * u2 + v1 * v2 +w1 * w2;
			double e22 = u2 * u2 + v2 * v2 +w2 * w2;

			double val=e11* t1* s1 + e12 * (t1 * s2 + t2 * s1) + e22 * t2 * s2;
			return val;
		}
		void guideA_u(double *ptr,double t1, double t2, double s1, double s2)
		{
			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;
			double u1 = 0, u2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];
				
				u2 += _ref->d1[1][s] * _ref->buf_u[s];
				
			}
			
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2*u1 * _ref->d1[0][s];
				double _e12 = u1 * _ref->d1[1][s]+u2 * _ref->d1[0][s];
				double _e22 = 2 * u2 * _ref->d1[1][s];

				double val=_e11* t1* s1 +_e12 * (t1 * s2 + t2 * s1) +_e22 * t2 * s2;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guideA_v(double* ptr, double t1, double t2, double s1, double s2)
		{
			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;

			double  v1 = 0, v2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
				
			}
			
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * v1 * _ref->d1[0][s];
				double _e12 = v1 * _ref->d1[1][s] + v2 * _ref->d1[0][s];
				double _e22 = 2 * v2 * _ref->d1[1][s];
				double val = _e11 * t1 * s1 + _e12 * (t1 * s2 + t2 * s1) + _e22 * t2 * s2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void guideA_w(double* ptr, double t1, double t2, double s1, double s2)
		{
			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;

			double w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				w1 += _ref->d1[0][s] * _ref->buf_w[s];

				w2 += _ref->d1[1][s] * _ref->buf_w[s];
			}
		
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * w1 * _ref->d1[0][s];
				double _e12 = w1 * _ref->d1[1][s] + w2 * _ref->d1[0][s];
				double _e22 = 2 * w2 * _ref->d1[1][s];
				double val = _e11 * t1 * s1 + _e12 * (t1 * s2 + t2 * s1) + _e22 * t2 * s2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double guideB()
		{
		

			double u1 = 0, u2 = 0, v1 = 0, v2 = 0 , w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				w1 += _ref->d1[0][s] * _ref->buf_w[s];
				u2 += _ref->d1[1][s] * _ref->buf_u[s];
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
				w2 += _ref->d1[1][s] * _ref->buf_w[s];
			}
			double e11 = u1 * u1 + v1 * v1 +w1 * w1;
			double e12 = u1 * u2 + v1 * v2 + w1 * w2;
			double e22 = u2 * u2 + v2 * v2 + w2 * w2;

			double val = (e11 * _ref->oG11 + 2*e12 * _ref->oG12+ e22 *_ref->oG22);
			
			return val;
		}
		void guideB_u(double* ptr)
		{
			

			double u1 = 0, u2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];
				
				u2 += _ref->d1[1][s] * _ref->buf_u[s];
			
			}
		
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * u1 * _ref->d1[0][s];
				double _e12 = u1 * _ref->d1[1][s] + u2 * _ref->d1[0][s];
				double _e22 = 2 * u2 * _ref->d1[1][s];

				double val = (_e11 * _ref->oG11 + 2 * _e12 * _ref->oG12 + _e22 * _ref->oG22);
				*ptr1 = val;
				ptr1++;
			}

		}
		void guideB_v(double* ptr)
		{
			

			double  v1 = 0, v2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{				
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * v1 * _ref->d1[0][s];
				double _e12 = v1 * _ref->d1[1][s] + v2 * _ref->d1[0][s];
				double _e22 = 2 * v2 * _ref->d1[1][s];
				double val = (_e11 * _ref->oG11 + 2 * _e12 * _ref->oG12 + _e22 * _ref->oG22);
				*ptr1 = val;
				ptr1++;
			}
		}
		void guideB_w(double* ptr)
		{
			

			double w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				w1 += _ref->d1[0][s] * _ref->buf_w[s];
				
				w2 += _ref->d1[1][s] * _ref->buf_w[s];
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * w1 * _ref->d1[0][s];
				double _e12 = w1 * _ref->d1[1][s] + w2 * _ref->d1[0][s];
				double _e22 = 2 * w2 * _ref->d1[1][s];
				double val = (_e11 * _ref->oG11 + 2 * _e12 * _ref->oG12 + _e22 * _ref->oG22);
				*ptr1 = val;
				ptr1++;
			}
		}
		double guideC(double t1,double t2,double s1,double s2,double q1,double q2)
		{

			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;

			double u1 = 0, u2 = 0, v1 = 0, v2 = 0 , w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				w1 += _ref->d1[0][s] * _ref->buf_w[s];
				u2 += _ref->d1[1][s] * _ref->buf_u[s];
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
				w2 += _ref->d1[1][s] * _ref->buf_w[s];
			}
			double e11 = u1 * u1 + v1 * v1 +w1 * w1;
			double e12 = u1 * u2 + v1 * v2 +w1 * w2;
			double e22 = u2 * u2 + v2 * v2 +w2 * w2;

			double val = q1 * (e11 * t1 * t1 + 2 * e12 * (t1 * t2 + t2 * t1) + e22 * (t2 * t2));
			val += q2 * (e11 * s1 * s1 + 2 * e12 * (s1 * s2 + s2 * s1) + e22 * (s2 * s2));

			return val;
		}
		void guideC_u(double* ptr, double t1, double t2, double s1, double s2, double q1, double q2)
		{
			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;


			double u1 = 0, u2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];

				u2 += _ref->d1[1][s] * _ref->buf_u[s];

			}

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * u1 * _ref->d1[0][s];
				double _e12 = u1 * _ref->d1[1][s] + u2 * _ref->d1[0][s];
				double _e22 = 2 * u2 * _ref->d1[1][s];

				double val = q1 * (_e11 * t1 * t1 + 2 * _e12 * (t1 * t2 + t2 * t1) + _e22 * (t2 * t2));
				val += q2 * (_e11 * s1 * s1 + 2 * _e12 * (s1 * s2 + s2 * s1) + _e22 * (s2 * s2));

				*ptr1 = val;
				ptr1++;
			}

		}
		void guideC_v(double* ptr, double t1, double t2, double s1, double s2, double q1, double q2)
		{
			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;


			double  v1 = 0, v2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				v1 += _ref->d1[0][s] * _ref->buf_v[s];
				v2 += _ref->d1[1][s] * _ref->buf_v[s];
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * v1 * _ref->d1[0][s];
				double _e12 = v1 * _ref->d1[1][s] + v2 * _ref->d1[0][s];
				double _e22 = 2 * v2 * _ref->d1[1][s];
				double val = q1 * (_e11 * t1 * t1 + 2 * _e12 * (t1 * t2 + t2 * t1) + _e22 * (t2 * t2));
				val += q2 * (_e11 * s1 * s1 + 2 * _e12 * (s1 * s2 + s2 * s1) + _e22 * (s2 * s2));
				*ptr1 = val;
				ptr1++;
			}
		}
		void guideC_w(double* ptr, double t1, double t2, double s1, double s2, double q1, double q2)
		{

			double length = sqrt(_ref->og11 * t1 * t1 + 2 * _ref->og12 * t1 * t2 + _ref->og22 * t2 * t2);
			t1 /= length;
			t2 /= length;
			length = sqrt(_ref->og11 * s1 * s1 + 2 * _ref->og12 * s1 * s2 + _ref->og22 * s2 * s2);
			s1 /= length;
			s2 /= length;

			double w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				w1 += _ref->d1[0][s] * _ref->buf_w[s];

				w2 += _ref->d1[1][s] * _ref->buf_w[s];
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 2 * w1 * _ref->d1[0][s];
				double _e12 = w1 * _ref->d1[1][s] + w2 * _ref->d1[0][s];
				double _e22 = 2 * w2 * _ref->d1[1][s];
				double val = q1 * (_e11 * t1 * t1 + 2 * _e12 * (t1 * t2 + t2 * t1) + _e22 * (t2 * t2));
				val += q2 * (_e11 * s1 * s1 + 2 * _e12 * (s1 * s2 + s2 * s1) + _e22 * (s2 * s2));
				*ptr1 = val;
				ptr1++;
			}
		}
		double guide(double t1,double t2,bool accurate)
		{
			double val = 0;
			double w1 = 0,  w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			val = (get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2);
			return val;
		}

		void guide_xi(double* ptr, double t1,double t2, bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e22 = __dsigma_22[0][s];
				/*for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}*/
				double _e21 = _e12;

				val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);
				//double dtr = _e11 * _ref->get__Gij(0, 0) + 2 * _e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);
				//val -= (get_eij(0, 0) * s1 * v1 + get_eij(1, 0) * s2 * v1 + get_eij(0, 1) * s1 * v2 + get_eij(1, 1) * s2 * v2) / treij / treij * dtr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_eta(double* ptr,double t1,double t2,bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				/*double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}*/
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e22 = __dsigma_22[1][s];
				double _e21 = _e12;

				val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);


				*ptr1 = val;
				ptr1++;
			}

		}

		void guide_nu(double* ptr,double t1,double t2,bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				/*double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
				}*/
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e22 = __dsigma_22[2][s];

				double _e21 = _e12;

				val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);


				*ptr1 = val;
				ptr1++;
			}

		}
	
		
		double guide1_2(double t1, double t2,double w1,double w2, bool accurate)
		{
			double val = 0;
			
			val = get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * (w1 * t2 + w2 * t1 )+ get_eij(1, 1) * w2 * t2;
			return val;
		}

		void guide1_2_xi(double* ptr, double t1, double t2,double w1,double w2,bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f22 = __dsigma_22[0][s];
				double _f21 = _f12;
				

				val = (_f11 * w1 * t1 + _f12 *( w1 * t2 + w2 * t1 )+ _f22 * w2 * t2);// / tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide1_2_eta( double* ptr, double t1, double t2, double w1, double w2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[1][s];
				double _f12 = __dsigma_12[1][s];
				double _f22 = __dsigma_22[1][s];
				double _f21 = _f12;


				val = (_f11 * w1 * t1 + _f12 * (w1 * t2 +  w2 * t1) + _f22 * w2 * t2);// / tr;
				*ptr1 = val;
				ptr1++;
			}

		}


		void guide1_2_nu(double* ptr, double t1, double t2, double w1,double w2,bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[2][s];
				double _f12 = __dsigma_12[2][s];
				double _f22 = __dsigma_22[2][s];
				double _f21 = _f12;


				val = (_f11 * w1 * t1 + _f12 * (w1 * t2 + w2 * t1) + _f22 * w2 * t2);// / tr;
				*ptr1 = val;
				ptr1++;
			}

		}

		
		double guide2(double t1, double t2, double c1, double c2, bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			val = c1 * (get_eij(0, 0) * t1 * t1 + get_eij(0, 1) * t1 * t2 + get_eij(1, 0) * t2 * t1 + get_eij(1, 1) * t2 * t2);
			val -= c2 * (get_eij(0, 0) * w1 * w1 + get_eij(0, 1) * w1 * w2 + get_eij(1, 0) * w2 * w1 + get_eij(1, 1) * w2 * w2);
			return val;
		}

		void guide2_xi(double* ptr, double t1, double t2, double c1, double c2, bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				/*for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_xi[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_xi[t];
				}*/
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e22 = __dsigma_22[0][s];

				double _e21 = _e12;

				val = c1 * (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val -= c2 * (_e11 * w1 * w1 + _e12 * w1 * w2 + _e21 * w2 * w1 + _e22 * w2 * w2);
				//double dtr = _e11 * _ref->get__Gij(0, 0) + 2 * _e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);
				//val -= (get_eij(0, 0) * s1 * v1 + get_eij(1, 0) * s2 * v1 + get_eij(0, 1) * s1 * v2 + get_eij(1, 1) * s2 * v2) / treij / treij * dtr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide2_eta(double* ptr, double t1, double t2, double c1, double c2, bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				/*double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_eta[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_eta[t];
				}*/
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e22 = __dsigma_22[1][s];

				double _e21 = _e12;

				val = c1 * (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val -= c2 * (_e11 * w1 * w1 + _e12 * w1 * w2 + _e21 * w2 * w1 + _e22 * w2 * w2);

				*ptr1 = val;
				ptr1++;
			}

		}

		void guide2_nu(double* ptr, double t1, double t2, double c1, double c2, bool accurate)
		{
			double val = 0;
			double w1 = 0, w2 = 0;
			if (accurate)
			{
				double length = get_gij2(0, 0) * t1 * t1 + 2 * get_gij2(0, 1) * t1 * t2 + get_gij2(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = get_gij2(0, 1) * t1 + get_gij2(1, 1) * t2;
				w2 = -get_gij2(0, 0) * t1 - get_gij2(0, 1) * t2;
				w1 /= dv;
				w2 /= dv;
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

				w1 = _ref->get__gij(0, 1) * t1 + _ref->get__gij(1, 1) * t2;
				w2 = -_ref->get__gij(0, 0) * t1 - _ref->get__gij(0, 1) * t2;
				w1 /= _ref->_refDv;
				w2 /= _ref->_refDv;
			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				/*double _e11 = 0, _e12 = 0, _e22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_e11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
					_e12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_nu[t];

					_e22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_nu[t];
				}*/
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e22 = __dsigma_22[2][s];

				double _e21 = _e12;

				val = c1 * (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);
				val -= c2 * (_e11 * w1 * w1 + _e12 * w1 * w2 + _e21 * w2 * w1 + _e22 * w2 * w2);

				*ptr1 = val;
				ptr1++;
			}

		}


		double guide2_2( double t1, double t2, double w1, double w2,double q1,double q2, bool accurate)
		{
			double val = 0;
			
			val = q1 * (get_eij(0, 0) * t1 * t1 + 2 * get_eij(0, 1) * t1 * t2 + get_eij(1, 1) * t2 * t2);
			val -= q2 * (get_eij(0, 0) * w1 * w1 + 2 * get_eij(0, 1) * w1 * w2 + get_eij(1, 1) * w2 * w2);

			return val;
		}

		void guide2_2_xi( double* ptr, double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f22 = __dsigma_22[0][s];
				double _f21 = _f12;
				
				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 *_e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = q1 * (_f11 * t1 * t1 + _f12 * t1 * t2 + _f21 * t2 * t1 + _f22 * t2 * t2);// / tr;
				val -= q2 * (_f11 * w1 * w1 + _f12 * w1 * w2 + _f21 * w2 * w1 + _f22 * w2 * w2);// / tr;
				//val += -(get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2) / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide2_2_eta( double* ptr, double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[1][s];
				double _f12 = __dsigma_12[1][s];
				double _f22 = __dsigma_22[1][s];
				double _f21 = _f12;

				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 *_e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = q1 * (_f11 * t1 * t1 + _f12 * t1 * t2 + _f21 * t2 * t1 + _f22 * t2 * t2);// / tr;
				val -= q2 * (_f11 * w1 * w1 + _f12 * w1 * w2 + _f21 * w2 * w1 + _f22 * w2 * w2);// / tr;
				//val += -(get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2) / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide2_2_nu( double* ptr, double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[2][s];
				double _f12 = __dsigma_12[2][s];
				double _f22 = __dsigma_22[2][s];
				double _f21 = _f12;

				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 *_e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = q1 * (_f11 * t1 * t1 + _f12 * t1 * t2 + _f21 * t2 * t1 + _f22 * t2 * t2);// / tr;
				val -= q2 * (_f11 * w1 * w1 + _f12 * w1 * w2 + _f21 * w2 * w1 + _f22 * w2 * w2);// / tr;
				//val += -(get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2) / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}

		}

		/*double traceA(double t1, double t2)
		{
			double val = 0;

			val =  (get_mij(0, 0) * t1 * t1 + 2 * get_mij(0, 1) * t1 * t2 + get_mij(1, 1) * t2 * t2);
		

			return val;
		}
		void traceA_xi( double* ptr, double t1, double t2)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e22 = __dsigma_22[0][s];

				double _e21 = _e12;


				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);// / tr;


				*ptr1 = val;
				ptr1++;
			}

		}
		void traceA_eta( double* ptr, double t1, double t2)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e22 = __dsigma_22[1][s];

				double _e21 = _e12;


				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);// / tr;


				*ptr1 = val;
				ptr1++;
			}

		}
		void traceA_nu( double* ptr, double t1, double t2)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e22 = __dsigma_22[2][s];

				double _e21 = _e12;


				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val =  (_e11 * t1 * t1 + _e12 * t1 * t2 + _e21 * t2 * t1 + _e22 * t2 * t2);// / tr;
				

				*ptr1 = val;
				ptr1++;
			}

		}*/

		double guide_supported(double t1, double t2,bool accurate)
		{
			double val = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

			}

			//tr = 1;
			val = (get_eij(0, 0) * t1 * t1 + 2 * get_eij(0, 1) * t1 * t2 + get_eij(1, 1) * t2 * t2);
			
			return val;
		}
		void guide_supported_xi(double* ptr, double t1, double t2,bool accurate)
		{
			double val = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[0][s];

				val =  (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);
				
				*ptr1 = val;
				ptr1++;
			}

		}
	


		void guide_supported_eta(double* ptr, double t1, double t2,bool accurate)
		{
			double val = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[1][s];


				val = (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);

				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_supported_nu(double* ptr, double t1, double t2,bool accurate)
		{
			double val = 0;
			if (accurate)
			{
				double length = get_gij(0, 0) * t1 * t1 + 2 * get_gij(0, 1) * t1 * t2 + get_gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);
			}
			else {
				double length = _ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2;
				t1 /= sqrt(length);
				t2 /= sqrt(length);

			}
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[2][s]; 

				val = (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);

				*ptr1 = val;
				ptr1++;
			}

		}
		/*double guide_supportedA()
		{
			
			return get_kij(0, 0);
		}
		double guide_supportedB()
		{
		
			return get_kij(1, 1);
		}
		double guide_supportedC()
		{

			return get_kij(0, 1);
		}*/
		/*void guide_supportedA_xi(double* ptr)
		{
			double val = 0;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f21 = _f12;
				double _f22 = __dsigma_22[0][s];
				double _A11 = _ref->get__Gij(0, 0) * _f11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _f12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _f21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _f22 * _ref->get__Gij(1, 0);
				double _A12 = _ref->get__Gij(0, 0) * _f11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _f12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _f21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _f22 * _ref->get__Gij(1, 1);
				double _A22 = _ref->get__Gij(1, 0) * _f11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _f12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _f21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _f22 * _ref->get__Gij(1, 1);

				val = _A11;

				*ptr1 = val;
				ptr1++;
			}

		}
	

		void guide_supportedB_xi( double* ptr)
		{
			double val = 0;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f21 = _f12;
				double _f22 = __dsigma_22[0][s];
				double _A11 = _ref->get__Gij(0, 0) * _f11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _f12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _f21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _f22 * _ref->get__Gij(1, 0);
				double _A12 = _ref->get__Gij(0, 0) * _f11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _f12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _f21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _f22 * _ref->get__Gij(1, 1);
				double _A22 = _ref->get__Gij(1, 0) * _f11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _f12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _f21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _f22 * _ref->get__Gij(1, 1);
				double val= _A22;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_supportedC_xi( double* ptr)
		{
			double val = 0;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f21 = _f12;
				double _f22 = __dsigma_22[0][s];
				double _A11 = _ref->get__Gij(0, 0) * _f11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _f12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _f21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _f22 * _ref->get__Gij(1, 0);
				double _A12 = _ref->get__Gij(0, 0) * _f11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _f12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _f21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _f22 * _ref->get__Gij(1, 1);
				double _A22 = _ref->get__Gij(1, 0) * _f11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _f12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _f21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _f22 * _ref->get__Gij(1, 1);
				double val = _A12;
				*ptr1 = val;
				ptr1++;
			}

		}*/
		double guide_supported2(double t1, double t2, bool accurate)
		{
			double val = 0;
			
			val = (get_eij(0, 0) * t1 * t1 + 2 * get_eij(0, 1) * t1 * t2 + get_eij(1, 1) * t2 * t2);

			return val;
		}
		void guide_supported2_xi(double* ptr, double t1, double t2, bool accurate)
		{
			double val = 0;
			
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[0][s];
				double _e12 = __dsigma_12[0][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[0][s];

				val = (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);

				*ptr1 = val;
				ptr1++;
			}

		}



		void guide_supported2_eta(double* ptr, double t1, double t2, bool accurate)
		{
			double val = 0;
			
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[1][s];


				val = (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);

				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_supported2_nu(double* ptr, double t1, double t2, bool accurate)
		{
			double val = 0;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[2][s];


				val = (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);

				*ptr1 = val;
				ptr1++;
			}

		}


		
		double guide_free2(double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;

			val = q1 * (get_eij(0, 0) * t1 * t1 + 2 * get_eij(0, 1) * t1 * t2 + get_eij(1, 1) * t2 * t2);
			val -= q2 * (get_eij(0, 0) * w1 * w1 + 2 * get_eij(0, 1) * w1 * w2 + get_eij(1, 1) * w2 * w2);

			return val;
		}

		void guide_free2_xi(double* ptr, double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f22 = __dsigma_22[0][s];
				double _f21 = _f12;

				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 *_e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = q1 * (_f11 * t1 * t1 + _f12 * t1 * t2 + _f21 * t2 * t1 + _f22 * t2 * t2);// / tr;
				val -= q2 * (_f11 * w1 * w1 + _f12 * w1 * w2 + _f21 * w2 * w1 + _f22 * w2 * w2);// / tr;
				//val += -(get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2) / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_free2_eta(double* ptr, double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[1][s];
				double _f12 = __dsigma_12[1][s];
				double _f22 = __dsigma_22[1][s];
				double _f21 = _f12;

				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 *_e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = q1 * (_f11 * t1 * t1 + _f12 * t1 * t2 + _f21 * t2 * t1 + _f22 * t2 * t2);// / tr;
				val -= q2 * (_f11 * w1 * w1 + _f12 * w1 * w2 + _f21 * w2 * w1 + _f22 * w2 * w2);// / tr;
				//val += -(get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2) / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void guide_free2_nu(double* ptr, double t1, double t2, double w1, double w2, double q1, double q2, bool accurate)
		{
			double val = 0;
			//double tr = get_eij(0, 0) * _ref->get__Gij(0, 0) + 2 * get_eij(0, 1) * _ref->get__Gij(0, 1) + get_eij(1, 1) * _ref->get__Gij(1, 1);

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//double _e11 = 0, _e12 = 0, _e22 = 0;
				double _f11 = __dsigma_11[2][s];
				double _f12 = __dsigma_12[2][s];
				double _f22 = __dsigma_22[2][s];
				double _f21 = _f12;

				//double _tr = _e11 * _ref->get__Gij(0, 0) + 2 *_e12 * _ref->get__Gij(0, 1) + _e22 * _ref->get__Gij(1, 1);

				//val = (_e11 * w1 * t1 + _e12 * w1 * t2 + _e21 * w2 * t1 + _e22 * w2 * t2);// / tr;
				val = q1 * (_f11 * t1 * t1 + _f12 * t1 * t2 + _f21 * t2 * t1 + _f22 * t2 * t2);// / tr;
				val -= q2 * (_f11 * w1 * w1 + _f12 * w1 * w2 + _f21 * w2 * w1 + _f22 * w2 * w2);// / tr;
				//val += -(get_eij(0, 0) * w1 * t1 + get_eij(0, 1) * w1 * t2 + get_eij(1, 0) * w2 * t1 + get_eij(1, 1) * w2 * t2) / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}
		}
	
			double dir()
			{
				double dx = 0;
					double dy = 0;

					for (int s = 0; s < _ref->_nNode; s++)
					{
						dx += _ref->d1[0][s] * _ref->buf_xi[s];
						dy += _ref->d1[0][s] * _ref->buf_eta[s];
					}
					return dx * _ref->get__gi(0, 1) - dy * _ref->get__gi(0, 0);

			}
			void dir_xi(double* ptr)
			{
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double val = _ref->d1[0][s] * _ref->get__gi(0, 1);
					*ptr1 = val;
					ptr1++;
				}
			}
			void dir_eta(double* ptr)
			{
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double val = -_ref->d1[0][s] * _ref->get__gi(0, 0);
					*ptr1 = val;
					ptr1++;
				}
			}
			double guide_free3( double t1, double t2, double s1, double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

			
				double huu = -globalratio * suu +get__hij(0, 0);
				double huv = -globalratio * suv + get__hij(0, 1);
				double hvu =  -globalratio * suv + get__hij(0, 1);
				double hvv =  -globalratio * svv + get__hij(1, 1);



				val =  (huu * s1 * t1 + huv * s1 * t2 + hvu * s2 * t1 + hvv * s2 * t2);
				


				return val;
			}	

			void guide_free3_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				double phi1 = 0, phi2 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					phi1 += _ref->d1[0][s] * _ref->buf_phi[s];
					phi2 += _ref->d1[1][s] * _ref->buf_phi[s];
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
	

					

					double _duu = xu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = xu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = xv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = xv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 0)+ _ref->d1[1][s] * _ref->get__gi(0, 0);
					double _g21 = _g12;
					double _g22 = 2*_ref->d1[1][s] * _ref->get__gi(1, 0);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0)-tr*_g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1)-tr*_g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0)-tr*_g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1)-tr*_g22);

					double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
					double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
					double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
					double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
					double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
					double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
					double _f11 = -_Gamma111 * phi1 - _Gamma112 * phi2;
					double _f12 = -_Gamma121 * phi1 - _Gamma122 * phi2;
					double _f22 = -_Gamma221 * phi1 - _Gamma222 * phi2;
					double _h11 = -globalratio * _suu + _f11;
					double _h12 = -globalratio * _suv + _f12;
					double _h21 = -globalratio * _suv + _f12;
					double _h22 = -globalratio * _svv + _f22;
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);

					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free3_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				double phi1 = 0, phi2 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					phi1 += _ref->d1[0][s] * _ref->buf_phi[s];
					phi2 += _ref->d1[1][s] * _ref->buf_phi[s];
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					;



					double _duu = yu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = yu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = yv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = yv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0,1);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);


					double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
					double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
					double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
					double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
					double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
					double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
					double _f11 = -_Gamma111 * phi1 - _Gamma112 * phi2;
					double _f12 = -_Gamma121 * phi1 - _Gamma122 * phi2;
					double _f22 = -_Gamma221 * phi1 - _Gamma222 * phi2;
					double _h11 = -globalratio * _suu + _f11;
					double _h12 = -globalratio * _suv + _f12;
					double _h21 = -globalratio * _suv + _f12;
					double _h22 = -globalratio * _svv + _f22;
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);

					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free3_phi(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _h11 = _ref->__dh[0][s];
					double _h12 = _ref->__dh[1][s];
					double _h21 = _ref->__dh[1][s];
					double _h22 = _ref->__dh[3][s];
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);

					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free3_xi( double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = _ref->d1[0][s];
					xv = _ref->d1[1][s];
					yu = 0;
					yv = 0;

					double duu = xu * _ref->get__gi(0, 0);// +yu * _ref->get__gi(0, 1);
					double duv = xu * _ref->get__gi(1, 0);// + yu * _ref->get__gi(1, 1);
					double dvu = xv * _ref->get__gi(0, 0);//+ yv * _ref->get__gi(0, 1);
					double dvv = xv * _ref->get__gi(1, 0);//+ yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

				
					double _h11 = -globalratio * suu;
					double _h12 = -globalratio * suv;
					double _h21 = -globalratio * suv;
					double _h22 = -globalratio * svv;
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);
					
					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free3_eta( double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
				xu = 0;
				xv = 0;
				yu = _ref->d1[0][s];
				yv = _ref->d1[1][s];

				double duu =  yu * _ref->get__gi(0, 1);
				double duv =  yu * _ref->get__gi(1, 1);
				double dvu =  yv * _ref->get__gi(0, 1);
				double dvv =  yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

				
				double _h11 = -globalratio * suu;
				double _h12 = -globalratio * suv;
				double _h21 = -globalratio * suv;
				double _h22 = -globalratio * svv;
					val =  (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);
					
					*ptr1 = val;
					ptr1++;
				}

			}
			
			double guide_free4(double t1, double t2,double s1,double s2, double globalratio,bool add,int mode)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));


				double huu = -globalratio * suu;
				double huv = -globalratio * suv;
				double hvu = -globalratio * suv;
				double hvv = -globalratio * svv;

				if (add)
				{
					huu += get__hij(0, 0);
					huv += get__hij(0, 1);
					hvu += get__hij(0, 1);
					hvv += get__hij(1, 1);
				}

				if (mode == 0)
				{
					val = (huu * t1 * t1 + huv * t1 * t2 + hvu * t2 * t1 + hvv * t2 * t2);
				}
				if (mode == 1)
				{
					val = (huu * t1 * s1 + huv * t1 * s2 + hvu * t2 * s1 + hvv * t2 * s2);
				}
				if (mode == 2)
				{
					val = (huu * s1 * s1 + huv * s1 * s2 + hvu * s2 * s1 + hvv * s2 * s2);
				}


				return val;
			}

			void guide_free4_phi(double* ptr, double t1, double t2, double s1, double s2, double globalratio,int mode)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _h11 = _ref->__dh[0][s];
					double _h12 = _ref->__dh[1][s];
					double _h21 = _ref->__dh[1][s];
					double _h22 = _ref->__dh[3][s];
					if (mode == 0)
					{
						val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					}
					if (mode == 1)
					{
						val = (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					}
					if (mode == 2)
					{
						val = (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);
					}


					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free4_xi(double* ptr, double t1, double t2, double s1, double s2, double globalratio,int mode)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = _ref->d1[0][s];
					xv = _ref->d1[1][s];
					yu = 0;
					yv = 0;

					double duu = xu * _ref->get__gi(0, 0);// +yu * _ref->get__gi(0, 1);
					double duv = xu * _ref->get__gi(1, 0);// + yu * _ref->get__gi(1, 1);
					double dvu = xv * _ref->get__gi(0, 0);//+ yv * _ref->get__gi(0, 1);
					double dvv = xv * _ref->get__gi(1, 0);//+ yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));


					double _h11 = -globalratio * suu;
					double _h12 = -globalratio * suv;
					double _h21 = -globalratio * suv;
					double _h22 = -globalratio * svv;
					if (mode == 0)
					{
						val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					}
					if (mode == 1)
					{
						val = (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					}
					if (mode == 2)
					{
						val = (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);
					}

					*ptr1 = val;
					ptr1++;
				}

			}

			void guide_free4_eta(double* ptr, double t1, double t2, double s1, double s2, double globalratio,int mode)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				 length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = 0;
					xv = 0;
					yu = _ref->d1[0][s];
					yv = _ref->d1[1][s];

					double duu =  yu * _ref->get__gi(0, 1);
					double duv =  yu * _ref->get__gi(1, 1);
					double dvu =  yv * _ref->get__gi(0, 1);
					double dvv =  yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));


					double _h11 = -globalratio * suu;
					double _h12 = -globalratio * suv;
					double _h21 = -globalratio * suv;
					double _h22 = -globalratio * svv;
					if (mode == 0)
					{
						val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					}
					if (mode == 1)
					{
						val = (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					}
					if (mode == 2)
					{
						val = (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);
					}

					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free4_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				 length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				double phi1 = 0, phi2 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					phi1 += _ref->d1[0][s] * _ref->buf_phi[s];
					phi2 += _ref->d1[1][s] * _ref->buf_phi[s];
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
		


					double _duu = xu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = xu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = xv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = xv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);

					double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
					double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
					double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
					double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
					double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
					double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
					double _f11 = -_Gamma111 * phi1 - _Gamma112 * phi2;
					double _f12 = -_Gamma121 * phi1 - _Gamma122 * phi2;
					double _f22 = -_Gamma221 * phi1 - _Gamma222 * phi2;

					double _h11 = -globalratio * _suu+_f11;
					double _h12 = -globalratio * _suv+_f12;
					double _h21 = -globalratio * _suv+_f12;
					double _h22 = -globalratio * _svv+_f22;
					val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					val -= (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);

					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free4_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				 length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				double phi1 = 0, phi2 = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					phi1 += _ref->d1[0][s] * _ref->buf_phi[s];
					phi2 += _ref->d1[1][s] * _ref->buf_phi[s];
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
				

					double _duu = yu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = yu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = yv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = yv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);

					double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
					double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
					double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
					double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
					double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
					double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
					double _f11 = -_Gamma111 * phi1 - _Gamma112 * phi2;
					double _f12 = -_Gamma121 * phi1 - _Gamma122 * phi2;
					double _f22 = -_Gamma221 * phi1 - _Gamma222 * phi2;

					double _h11 = -globalratio * _suu+_f11 ;
					double _h12 = -globalratio * _suv+_f12 ;
					double _h21 = -globalratio * _suv +_f12;
					double _h22 = -globalratio * _svv+_f22 ;
					val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					val -= (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);

					*ptr1 = val;
					ptr1++;
				}

			}
			double guide_free5(double t1, double t2, double s1, double s2, double globalratio, bool add, int mode)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				suv = 0.5 * (suv + svu);

				double huu = suu;
				double huv = suv;
				double hvu = suv;
				double hvv = svv;

				if (add)
				{
					huu += get__hij(0, 0);
					huv += get__hij(0, 1);
					hvu += get__hij(0, 1);
					hvv += get__hij(1, 1);
				}

				if (mode == 0)
				{
					val = (huu * t1 * t1 + huv * t1 * t2 + hvu * t2 * t1 + hvv * t2 * t2);
				}
				if (mode == 1)
				{
					val = (huu * t1 * s1 + huv * t1 * s2 + hvu * t2 * s1 + hvv * t2 * s2);
				}
				if (mode == 2)
				{
					val = (huu * s1 * s1 + huv * s1 * s2 + hvu * s2 * s1 + hvv * s2 * s2);
				}


				return val;
			}

			void guide_free5_phi(double* ptr, double t1, double t2, double s1, double s2, double globalratio, int mode)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _h11 = _ref->__dh[0][s];
					double _h12 = _ref->__dh[1][s];
					double _h21 = _ref->__dh[1][s];
					double _h22 = _ref->__dh[3][s];
					if (mode == 0)
					{
						val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					}
					if (mode == 1)
					{
						val = (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					}
					if (mode == 2)
					{
						val = (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);
					}


					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_free5_xi(double* ptr, double t1, double t2, double s1, double s2, double globalratio, int mode)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = _ref->d1[0][s];
					xv = _ref->d1[1][s];
					yu = 0;
					yv = 0;

					double _suu = xu * _ref->get__gi(0, 0);// +yu * _ref->get__gi(0, 1);
					double _suv = xu * _ref->get__gi(1, 0);// + yu * _ref->get__gi(1, 1);
					double _svu = xv * _ref->get__gi(0, 0);//+ yv * _ref->get__gi(0, 1);
					double _svv = xv * _ref->get__gi(1, 0);//+ yv * _ref->get__gi(1, 1);
					_suv = 0.5 * (_suv + _svu);
					double _h11 = _suu;
					double _h12 = _suv;
					double _h21 = _svu;
					double _h22 = _svv;
					if (mode == 0)
					{
						val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					}
					if (mode == 1)
					{
						val = (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					}
					if (mode == 2)
					{
						val = (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);
					}

					*ptr1 = val;
					ptr1++;
				}

			}

			void guide_free5_eta(double* ptr, double t1, double t2, double s1, double s2, double globalratio, int mode)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = 0;
					xv = 0;
					yu = _ref->d1[0][s];
					yv = _ref->d1[1][s];

					double _suu = xu * _ref->get__gi(0, 0);// +yu * _ref->get__gi(0, 1);
					double _suv = xu * _ref->get__gi(1, 0);// + yu * _ref->get__gi(1, 1);
					double _svu = xv * _ref->get__gi(0, 0);//+ yv * _ref->get__gi(0, 1);
					double _svv = xv * _ref->get__gi(1, 0);//+ yv * _ref->get__gi(1, 1);
					_suv = 0.5 * (_suv + _svu);
					double _h11 = _suu;
					double _h12 = _suv;
					double _h21 = _svu;
					double _h22 = _svv;
					if (mode == 0)
					{
						val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					}
					if (mode == 1)
					{
						val = (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					}
					if (mode == 2)
					{
						val = (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);
					}

					*ptr1 = val;
					ptr1++;
				}

			}
			


			double guide_Free3(double t1, double t2, double s1, double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));


				double huu = -suu;
				double huv = -suv;
				double hvu = -suv;
				double hvv = -svv;



				val = (huu * s1 * t1 + huv * s1 * t2 + hvu * s2 * t1 + hvv * s2 * t2);



				return val;
			}

			void guide_Free3_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _xu = _ref->d1[0][s];
					double _xv = _ref->d1[1][s];
					double _yu = 0;
					double _yv = 0;
					double _duu = _xu * _ref->get__gi(0, 0) + xu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = _xu * _ref->get__gi(1, 0) + xu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = _xv * _ref->get__gi(0, 0) + xv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = _xv * _ref->get__gi(1, 0) + xv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);

					
					double _h11 = -_suu;
					double _h12 = -_suv;
					double _h21 = -_suv;
					double _h22 = -_svv;
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);

					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_Free3_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _yu = _ref->d1[0][s];
					double _yv = _ref->d1[1][s];
				
					double _duu = _yu * _ref->get__gi(0, 1) + yu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = _yu * _ref->get__gi(1, 1) + yu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = _yv * _ref->get__gi(0, 1) + yv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = _yv * _ref->get__gi(1, 1) + yv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);


					double _h11 = -_suu;
					double _h12 = -_suv;
					double _h21 = -_suv;
					double _h22 = -_svv;
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);

					*ptr1 = val;
					ptr1++;
				}

			}
		
			double guide_Free4(double t1, double t2,double s1,double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));


				double huu = -suu;
				double huv = -suv;
				double hvu = -suv;
				double hvv = -svv;



				val = (huu * t1 * t1 + huv * t1 * t2 + hvu * t2 * t1 + hvv * t2 * t2);
				val -= (huu * t1 * s1 + huv * t1 * s2 + hvu * t2 * s1 + hvv * t2 * s2);



				return val;
			}

			void guide_Free4_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _xu = _ref->d1[0][s];
					double _xv = _ref->d1[1][s];
					double _yu = 0;
					double _yv = 0;
					double _duu = _xu * _ref->get__gi(0, 0) + xu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = _xu * _ref->get__gi(1, 0) + xu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = _xv * _ref->get__gi(0, 0) + xv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = _xv * _ref->get__gi(1, 0) + xv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);


					double _h11 = -_suu;
					double _h12 = -_suv;
					double _h21 = -_suv;
					double _h22 = -_svv;
					val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					val -= (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					*ptr1 = val;
					ptr1++;
				}

			}
			void guide_Free4_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _yu = _ref->d1[0][s];
					double _yv = _ref->d1[1][s];

					double _duu = _yu * _ref->get__gi(0, 1) + yu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
					double _duv = _yu * _ref->get__gi(1, 1) + yu * _ref->d1[1][s];// +yu * _ref->d1[1][s];
					double _dvu = _yv * _ref->get__gi(0, 1) + yv * _ref->d1[0][s];// +yv * _ref->d1[0][s];
					double _dvv = _yv * _ref->get__gi(1, 1) + yv * _ref->d1[1][s];// +yv * _ref->d1[1][s];

					double _g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
					double _g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
					double _g21 = _g12;
					double _g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;
					double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1)
						+ duu * _G11 + duv * _G12 + dvu * _G21 + dvv * _G22;

					double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0) - tr * _g11);
					double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1) - tr * _g12);
					double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0) - tr * _g21);
					double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1) - tr * _g22);


					double _h11 = -_suu;
					double _h12 = -_suv;
					double _h21 = -_suv;
					double _h22 = -_svv;
					val = (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					val -= (_h11 * t1 * s1 + _h12 * t1 * s2 + _h21 * t2 * s1 + _h22 * t2 * s2);
					*ptr1 = val;
					ptr1++;
				}

			}

			double guide_Free5(double t1, double t2, double s1, double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0 , wu = 0, wv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
					wu += _ref->d1[0][s] * _ref->buf_w[s];
					wv += _ref->d1[1][s] * _ref->buf_w[s];
				}
				double euu = xu * xu + yu * yu + wu * wu;
				double euv = xu * xv + yu * yv + wu * wv;
				double evv = xv * xv + yv * yv + wv * wv;
				double evu = euv;

				val = (euu * s1 * t1 + euv * s1 * t2 + evu * s2 * t1 + evv * s2 * t2);



				return val;
			}

			void guide_Free5_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * xu * _ref->d1[0][s];
					double _euv = xu * _ref->d1[1][s]+ xv * _ref->d1[0][s];
					double _evv = 2 * xv * _ref->d1[1][s];
					double _evu = _euv;
					double val= (_euu * s1 * t1 + _euv * s1 * t2 + _evu * s2 * t1 + _evv * s2 * t2);
					*ptr1 = val;
					ptr1++;

				}

			}
			void guide_Free5_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double  yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
				
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * yu * _ref->d1[0][s];
					double _euv = yu * _ref->d1[1][s] +yv * _ref->d1[0][s];
					double _evv = 2 * yv * _ref->d1[1][s];
					double _evu = _euv;
					double val = (_euu * s1 * t1 + _euv * s1 * t2 + _evu * s2 * t1 + _evv * s2 * t2);
					*ptr1 = val;
					ptr1++;

				}

			}
			void guide_Free5_w(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double  wu = 0, wv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					wu += _ref->d1[0][s] * _ref->buf_w[s];
					wv += _ref->d1[1][s] * _ref->buf_w[s];
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * wu * _ref->d1[0][s];
					double _euv = wu * _ref->d1[1][s] + wv * _ref->d1[0][s];
					double _evv = 2 * wv * _ref->d1[1][s];
					double _evu = _euv;
					double val = (_euu * s1 * t1 + _euv * s1 * t2 + _evu * s2 * t1 + _evv * s2 * t2);
					*ptr1 = val;
					ptr1++;

				}

			}

			double guide_Free6(double t1, double t2, double s1, double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0 , wu = 0, wv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
					wu += _ref->d1[0][s] * _ref->buf_w[s];
					wv += _ref->d1[1][s] * _ref->buf_w[s];
				}
				double euu = xu * xu + yu * yu+ wu * wu;
				double euv = xu * xv + yu * yv+ wu * wv;
				double evv = xv * xv + yv * yv+ wv * wv;
				double evu = euv;


				val = (euu * s1 * t1 + euv * s1 * t2 + evu * s2 * t1 + evv * s2 * t2);
				val -= (euu * t1 * t1 + euv * t1 * t2 + evu * t2 * t1 + evv * t2 * t2);



				return val;
			}

			void guide_Free6_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * xu * _ref->d1[0][s];
					double _euv = xu * _ref->d1[1][s] + xv * _ref->d1[0][s];
					double _evv = 2 * xv * _ref->d1[1][s];
					double _evu = _euv;
					double val = (_euu * s1 * t1 + _euv * s1 * t2 + _evu * s2 * t1 + _evv * s2 * t2);
					val -= (_euu * t1 * t1 + _euv * t1 * t2 + _evu * t2 * t1 + _evv * t2 * t2);
					*ptr1 = val;
					ptr1++;

				}

			}
			void guide_Free6_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
				
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * yu * _ref->d1[0][s];
					double _euv = yu * _ref->d1[1][s] + yv * _ref->d1[0][s];
					double _evv = 2 * yv * _ref->d1[1][s];
					double _evu = _euv;
					double val = (_euu * s1 * t1 + _euv * s1 * t2 + _evu * s2 * t1 + _evv * s2 * t2);
					val -= (_euu * t1 * t1 + _euv * t1 * t2 + _evu * t2 * t1 + _evv * t2 * t2);
					*ptr1 = val;
					ptr1++;

				}

			}
			void guide_Free6_w(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double wu = 0, wv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					wu += _ref->d1[0][s] * _ref->buf_w[s];
					wv += _ref->d1[1][s] * _ref->buf_w[s];
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * wu * _ref->d1[0][s];
					double _euv = wu * _ref->d1[1][s] + wv * _ref->d1[0][s];
					double _evv = 2 * wv * _ref->d1[1][s];
					double _evu = _euv;
					double val = (_euu * s1 * t1 + _euv * s1 * t2 + _evu * s2 * t1 + _evv * s2 * t2);
					val -= (_euu * t1 * t1 + _euv * t1 * t2 + _evu * t2 * t1 + _evv * t2 * t2);
					*ptr1 = val;
					ptr1++;

				}

			}

			double guide_Free7(double t1, double t2, double s1, double s2, double globalratio)
			{
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double val = 0;
				double xu = 0, xv = 0, yu = 0, yv = 0, wu = 0, wv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
					wu += _ref->d1[0][s] * _ref->buf_w[s];
					wv += _ref->d1[1][s] * _ref->buf_w[s];
				}
				double euu = xu * xu + yu * yu+ wu * wu;
				double euv = xu * xv + yu * yv+ wu * wv;
				double evv = xv * xv + yv * yv+ wv * wv;
				double evu = euv;

			
				val = (euu * t1 * t1 + euv * t1 * t2 + evu * t2 * t1 + evv * t2 * t2);
				val += (euu * s1 * s1 + euv * s1 * s2 + evu * s2 * s1 + evv * s2 * s2);
				val -= 2;


				return val;
			}

			void guide_Free7_u(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double xu = 0, xv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_u[s];
					xv += _ref->d1[1][s] * _ref->buf_u[s];
					
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * xu * _ref->d1[0][s];
					double _euv = xu * _ref->d1[1][s] + xv * _ref->d1[0][s];
					double _evv = 2 * xv * _ref->d1[1][s];
					double _evu = _euv;
					
					double val = (_euu * t1 * t1 + _euv * t1 * t2 + _evu * t2 * t1 + _evv * t2 * t2);
					val += (_euu * s1 * s1 + _euv * s1 * s2 + _evu * s2 * s1 + _evv * s2 * s2);
					*ptr1 = val;
					ptr1++;

				}

			}
			void guide_Free7_v(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					
					yu += _ref->d1[0][s] * _ref->buf_v[s];
					yv += _ref->d1[1][s] * _ref->buf_v[s];
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * yu * _ref->d1[0][s];
					double _euv = yu * _ref->d1[1][s] + yv * _ref->d1[0][s];
					double _evv = 2 * yv * _ref->d1[1][s];
					double _evu = _euv;
					
					double val = (_euu * t1 * t1 + _euv * t1 * t2 + _evu * t2 * t1 + _evv * t2 * t2);
					val += (_euu * s1 * s1 + _euv * s1 * s2 + _evu * s2 * s1 + _evv * s2 * s2);
					*ptr1 = val;
					ptr1++;

				}

			}
			void guide_Free7_w(double* ptr, double t1, double t2, double s1, double s2, double globalratio)
			{
				double val = 0;
				double length = sqrt(t1 * t1 * _ref->og11 + 2 * t1 * t2 * _ref->og12 + t2 * t2 * _ref->og22);
				t1 /= length;
				t2 /= length;
				length = sqrt(s1 * s1 * _ref->og11 + 2 * s1 * s2 * _ref->og12 + s2 * s2 * _ref->og22);
				s1 /= length;
				s2 /= length;
				double wu = 0, wv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					wu += _ref->d1[0][s] * _ref->buf_w[s];
					wv += _ref->d1[1][s] * _ref->buf_w[s];
				}
				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{

					double _euu = 2 * wu * _ref->d1[0][s];
					double _euv = wu * _ref->d1[1][s] + wv * _ref->d1[0][s];
					double _evv = 2 * wv * _ref->d1[1][s];
					double _evu = _euv;

					double val = (_euu * t1 * t1 + _euv * t1 * t2 + _evu * t2 * t1 + _evv * t2 * t2);
					val += (_euu * s1 * s1 + _euv * s1 * s2 + _evu * s2 * s1 + _evv * s2 * s2);
					*ptr1 = val;
					ptr1++;

				}

			}
			double __sc()
			{
				return sc;
			}
			double det_sigma()
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));
				double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
				double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
				double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

				double fuu = get__hij(0, 0) + Svv / sc;
				double fuv = get__hij(0, 1) - Suv / sc;
				double fvu = get__hij(1, 0) - Suv / sc;
				double fvv = get__hij(1, 1) + Suu / sc;


				return fuu * fvv - fuv * fvu;

			}
			double Suu_sigma()
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));
				double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
				double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
				double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

				double fuu = get__hij(0, 0) + Svv / sc;
				double fuv = get__hij(0, 1) - Suv / sc;
				double fvu = get__hij(1, 0) - Suv / sc;
				double fvv = get__hij(1, 1) + Suu / sc;


				return Suu;
			}
			double Suv_sigma()
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));
				double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
				double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
				double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

				double fuu = get__hij(0, 0) + Svv / sc;
				double fuv = get__hij(0, 1) - Suv / sc;
				double fvu = get__hij(1, 0) - Suv / sc;
				double fvv = get__hij(1, 1) + Suu / sc;


				return Suv;

			}
			double Svv_sigma()
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));
				double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
				double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
				double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

				double fuu = get__hij(0, 0) + Svv / sc;
				double fuv = get__hij(0, 1) - Suv / sc;
				double fvu = get__hij(1, 0) - Suv / sc;
				double fvv = get__hij(1, 1) + Suu / sc;


				return Svv;

			}
			double fuu_sigma()
			{

				double fuu = get__hij(0, 0);



				return fuu;
			}
			double fuv_sigma()
			{

				double fuv = get__hij(0, 1);



				return fuv;
			}
			double fvv_sigma()
			{

				double fvv = get__hij(1, 1);



				return fvv;
			}

			double parallel(double t1, double t2, bool add)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;


				
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

				double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
				double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
				double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

				double val=Suu* t1* t2 - t1 * t1 * Suv;



				return val;
			}


			void parallel_xi(double* ptr, double t1, double t2)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;


				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = _ref->d1[0][s];
					xv = _ref->d1[1][s];
					//yu = _ref->d1[0][s] ;
					//yv = _ref->d1[1][s] ;

					double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
					double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
					double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
					double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

					double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
					double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
					double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

					double val = Suu * t1 * t2 - t1 * t1 * Suv;

					*ptr1 = val;
					ptr1++;
				}

			}
			void parallel_eta(double* ptr, double t1, double t2)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;

			

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					//xu = _ref->d1[0][s];
					//xv = _ref->d1[1][s];
					yu = _ref->d1[0][s] ;
					yv = _ref->d1[1][s] ;

					double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
					double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
					double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
					double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

					double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
					double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
					double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

					double val = Suu * t1 * t2 - t1 * t1 * Suv;

					*ptr1 = val;
					ptr1++;
				}

			}


			double parallel2(double t1, double t2, bool add)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;


			
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double suu = (2 * duu - tr * _ref->get__gij(0, 0));
				double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

				double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
				double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
				double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);

				double val = Svv * t1 * t2 - t2 * t2 * Suv;



				return val;
			}


			void parallel2_xi(double* ptr, double t1, double t2)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;

				

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					xu = _ref->d1[0][s];
					xv = _ref->d1[1][s];
					//yu = _ref->d1[0][s] ;
					//yv = _ref->d1[1][s] ;

					double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
					double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
					double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
					double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

					double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
					double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
					double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);


					double val = Svv * t1 * t2 - t2 * t2 * Suv;

					*ptr1 = val;
					ptr1++;
				}

			}
			void parallel2_eta(double* ptr, double t1, double t2)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;

				

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double xu = 0, xv = 0, yu = 0, yv = 0;
					//xu = _ref->d1[0][s];
					//xv = _ref->d1[1][s];
					yu = _ref->d1[0][s];
					yv = _ref->d1[1][s];

					double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
					double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
					double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
					double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

					double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

					double suu = (2 * duu - tr * _ref->get__gij(0, 0));
					double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
					double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
					double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

					double Suu = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 0);
					double Suv = _ref->get__Gij(0, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * svv * _ref->get__Gij(1, 1);
					double Svv = _ref->get__Gij(1, 0) * suu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * suv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * svu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * svv * _ref->get__Gij(1, 1);


					double val = Svv * t1 * t2 - t2 * t2 * Suv;

					*ptr1 = val;
					ptr1++;
				}

			}
			
			
			double guide_freeA( double t1, double t2, double s1, double s2, double w1, double w2, bool accurate)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;
				length = sqrt(_ref->get__gij(0, 0) * s1 * s1 + 2 * _ref->get__gij(0, 1) * s1 * s2 + _ref->get__gij(1, 1) * s2 * s2);
				s1 /= length;
				s2 /= length;


				double val = 0;
				
				double suu = get__hij(0, 0);
				double suv = get__hij(0, 1);
				double svu = get__hij(1, 0);
				double svv = get__hij(1, 1);

				val = w1 * (suu * t1 * t1 + suv * t1 * t2 + svu * t2 * t1 + svv * t2 * t2);
				val -= w2 * (suu * s1 * s1 + suv * s1 * s2 + svu * s2 * s1 + svv * s2 * s2);

				return val;
			}
			void guide_freeA_phi( double* ptr, double t1, double t2, double s1, double s2, double w1, double w2, bool accurate)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;
				length = sqrt(_ref->get__gij(0, 0) * s1 * s1 + 2 * _ref->get__gij(0, 1) * s1 * s2 + _ref->get__gij(1, 1) * s2 * s2);
				s1 /= length;
				s2 /= length;
				double val = 0;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _h11 = _ref->__dh[0][s];
					double _h12 = _ref->__dh[1][s];
					double _h21 = _ref->__dh[2][s];
					double _h22 = _ref->__dh[3][s];


					val = w1 * (_h11 * t1 * t1 + _h12 * t1 * t2 + _h21 * t2 * t1 + _h22 * t2 * t2);
					val -= w2 * (_h11 * s1 * s1 + _h12 * s1 * s2 + _h21 * s2 * s1 + _h22 * s2 * s2);

					*ptr1 = val;
					ptr1++;
				}

			}
		
			double guide_freeB( double t1, double t2, double s1, double s2, double w1, double w2, bool accurate)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;
				length = sqrt(_ref->get__gij(0, 0) * s1 * s1 + 2 * _ref->get__gij(0, 1) * s1 * s2 + _ref->get__gij(1, 1) * s2 * s2);
				s1 /= length;
				s2 /= length;
			
				double suu = get__hij(0, 0);
				double suv = get__hij(0, 1);
				double svu = get__hij(1, 0);
				double svv = get__hij(1, 1);
				double val = 0;


				val = (suu * s1 * t1 + suv * s1 * t2 + svu * s2 * t1 + svv * s2 * t2);



				return val;
			}

			void guide_freeB_phi( double* ptr, double t1, double t2, double s1, double s2, double w1, double w2, bool accurate)
			{
				double length = sqrt(_ref->get__gij(0, 0) * t1 * t1 + 2 * _ref->get__gij(0, 1) * t1 * t2 + _ref->get__gij(1, 1) * t2 * t2);
				t1 /= length;
				t2 /= length;
				length = sqrt(_ref->get__gij(0, 0) * s1 * s1 + 2 * _ref->get__gij(0, 1) * s1 * s2 + _ref->get__gij(1, 1) * s2 * s2);
				s1 /= length;
				s2 /= length;
				double val = 0;

				double* ptr1 = ptr;
				for (int s = 0; s < _ref->_nNode; s++)
				{
					double _h11 = _ref->__dh[0][s];
					double _h12 = _ref->__dh[1][s];
					double _h21 = _ref->__dh[2][s];
					double _h22 = _ref->__dh[3][s];
					val = (_h11 * s1 * t1 + _h12 * s1 * t2 + _h21 * s2 * t1 + _h22 * s2 * t2);

					*ptr1 = val;
					ptr1++;
				}

			}
		
		/*void guide_free2_eta(double* ptr, double t1, double t2, double s1, double s2, double w1, double w2, bool accurate)
		{
			double val = 0;
			

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[1][s];
				double _e12 = __dsigma_12[1][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[1][s];



				val = w1 * (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);
				val -= w2 * (_e11 * s1 * s1 + 2 * _e12 * s1 * s2 + _e22 * s2 * s2);
				*ptr1 = val;
				ptr1++;
			}

		}*/
		/*void guide_free2_nu( double* ptr, double t1, double t2, double s1, double s2, double w1, double w2, bool accurate)
		{
			double val = 0;
		
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = __dsigma_11[2][s];
				double _e12 = __dsigma_12[2][s];
				double _e21 = _e12;
				double _e22 = __dsigma_22[2][s];


				val = w1 * (_e11 * t1 * t1 + 2 * _e12 * t1 * t2 + _e22 * t2 * t2);
				val -= w2 * (_e11 * s1 * s1 + 2 * _e12 * s1 * s2 + _e22 * s2 * s2);

				*ptr1 = val;
				ptr1++;
			}

		}*/

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
				val = (s11_phi * E11 * _S12 + s11_phi * E12 * _S22 + s12_phi * E21 * _S12 + s12_phi * E22 * _S22);// / sc;
				val -= (s21_phi * E11 * _S11 + s21_phi * E12 * _S21 + s22_phi * E21 * _S11 + s22_phi * E22 * _S21);// / sc;

				*ptr1 = val;
				ptr1++;
			}

		}

		double align_mix( double v1, double v2, double w1, double w2)
		{
			double val = 0;
			double length = get_gij2(0, 0) * v1 * v1 + 2 * get_gij2(0, 1) * v1 * v2 + get_gij2(1, 1) * v2 * v2;
			v1 /= sqrt(length);
			v2 /= sqrt(length);

			double s1 = (v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / dv;
			double s2 = (-v1 * get_gij2(0, 0) - v2 * get_gij2(0, 1)) / dv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			val = (get__hij(0, 0) * E11 * get__Sij(0, 1) + get__hij(0, 0) * E12 * get__Sij(1, 1) + get__hij(0, 1) * E21 * get__Sij(0, 1) + get__hij(0, 1) * E22 * get__Sij(1, 1)) * scale;
			val -= (get__hij(1, 0) * E11 * get__Sij(0, 0) + get__hij(1, 0) * E12 * get__Sij(1, 0) + get__hij(1, 1) * E21 * get__Sij(0, 0) + get__hij(1, 1) * E22 * get__Sij(1, 0)) * scale;

			return val;

		}


		void align_mix_z( double* ptr, double v1, double v2, double w1, double w2)
		{
			double val = 0;
			double length = get_gij2(0, 0) * v1 * v1 + 2 * get_gij2(0, 1) * v1 * v2 + get_gij2(1, 1) * v2 * v2;
			v1 /= sqrt(length);
			v2 /= sqrt(length);

			double s1 = (v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / dv;
			double s2 = (-v1 * get_gij2(0, 0) - v2 * get_gij2(0, 1)) / dv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				/*double _g11 = 0, _g12 = 0, _g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}*/
				//double _g21 = _g12;
				//double ddv = 0.5 * (_g11 * get_Gij2(0, 0) + 2 * _g12 * get_Gij2(0, 1) + _g22 * get_Gij2(1, 1)) * dv;
				//double _s1 = (v1 * _g12 + v2 * _g22) / dv;
				//double _s2 = (-v1 * _g11 - v2 * _g12) / dv;
				//_s1 += -(v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / dv / dv * ddv;
				//_s2 += -(-v1 * get_gij2(0, 0) - v2 * get_gij2(0, 1)) / dv / dv * ddv;
				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;
				val = (get__hij(0, 0) * E11 * _S12 + get__hij(0, 0) * E12 * _S22 + get__hij(0, 1) * E21 * _S12 + get__hij(0, 1) * E22 * _S22) * scale;
				val -= (get__hij(1, 0) * E11 * _S11 + get__hij(1, 0) * E12 * _S21 + get__hij(1, 1) * E21 * _S11 + get__hij(1, 1) * E22 * _S21) * scale;
				//val += (get__hij(0, 0) * (2 * w2 * _s1 * s1) * get__Sij(0, 1) + get__hij(0, 0) * (w2 * s1 * _s2 + w2 * _s1 * s2) * get__Sij(1, 1) + get__hij(0, 1) * (w2 * s1 * _s2 + w2 * _s1 * s2) * get__Sij(0, 1) + get__hij(0, 1) * (2 * w2 * s2 * _s2) * get__Sij(1, 1)) * scale;
				//val -= (get__hij(1, 0) * (2 * w2 * _s1 * s1) * get__Sij(0, 0) + get__hij(1, 0) * (w2 * s1 * _s2 + w2 * _s1 * s2) * get__Sij(1, 0) + get__hij(1, 1) * (w2 * s1 * _s2 + w2 * _s1 * s2) * get__Sij(0, 0) + get__hij(1, 1) * (2 * w2 * s2 * _s2) * get__Sij(1, 0)) * scale;

				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix_phi( double* ptr, double v1, double v2, double w1, double w2)
		{

			double val = 0;
			double length = get_gij2(0, 0) * v1 * v1 + 2 * get_gij2(0, 1) * v1 * v2 + get_gij2(1, 1) * v2 * v2;
			v1 /= sqrt(length);
			v2 /= sqrt(length);

			double s1 = (v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / dv;
			double s2 = (-v1 * get_gij2(0, 0) - v2 * get_gij2(0, 1)) / dv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][s] ;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][s] ;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;
				double _h21 = _h12;
				val = (_h11 * E11 * get__Sij(0, 1) + _h11 * E12 * get__Sij(1, 1) + _h12 * E21 * get__Sij(0, 1) + _h12 * E22 * get__Sij(1, 1)) * scale;
				val -= (_h21 * E11 * get__Sij(0, 0) + _h21 * E12 * get__Sij(1, 0) + _h22 * E21 * get__Sij(0, 0) + _h22 * E22 * get__Sij(1, 0)) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void __mix22_mat_phi(double v1, double v2, double s1, double s2, double w1, double w2)
		{
		

			if (_ref->__matD_phi == 0)
			{
				_ref->__matD_phi = new double[_nNode * _nNode];
				double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
				double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
				double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
				double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


				double scale = 1 / _ref->_refDv;
				


				for (int s = 0; s < _ref->_nNode; s++) {
					double huu = _ref->__dh[0][s];
					double huv = _ref->__dh[1][s];
					double hvu = _ref->__dh[2][s];
					double hvv = _ref->__dh[3][s];
					for (int t = 0; t < _ref->_nNode; t++) {

						double suu = _ref->__dh[0][t];
						double suv = _ref->__dh[1][t];
						double svu = _ref->__dh[2][t];
						double svv = _ref->__dh[3][t];

					
						_ref->__matD_phi[s * _nNode + t] = (suu*E11 * huv + suu*E12 * hvv + suv*E21 * huv + suv*E22 * hvv) * scale
						 - (svu*E11 * huu + svu*E12 * hvu + svv*E21 * huu + svv*E22 * hvu) * scale;
					}
				}
			}


		}

		void __mix22_mat_xi(double v1, double v2, double s1, double s2, double w1, double w2 )
		{
			if (_ref->__matD_xi == 0)
			{
				_ref->__matD_xi = new double[_nNode * _nNode];
				double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
				double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
				double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
				double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


				double scale = 1 / _ref->_refDv;


				for (int s = 0; s < _ref->_nNode; s++) {
					double huu = _ref->__dh[0][s];
					double huv = _ref->__dh[1][s];
					double hvu = _ref->__dh[2][s];
					double hvv = _ref->__dh[3][s];
					for (int t = 0; t < _ref->_nNode; t++) {

						double xu = 0, xv = 0, yu = 0, yv = 0;


						xu = _ref->d1[0][t];
						xv = _ref->d1[1][t];
						yu = 0;
						yv = 0;

						double duu = xu * _ref->get__gi(0, 0);
						double duv = xu * _ref->get__gi(1, 0);
						double dvu = xv * _ref->get__gi(0, 0);
						double dvv = xv * _ref->get__gi(1, 0);

						double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

						double suu = (2 * duu - tr * _ref->get__gij(0, 0));
						double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
						double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
						double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

						suu = -suu;
						suv = -suv;
						svu = suv;
						svv = -svv;

						_ref->__matD_xi[s * _nNode + t] = (suu * E11 * huv + suu * E12 * hvv + suv * E21 * huv + suv * E22 * hvv) * scale
							- (svu * E11 * huu + svu * E12 * hvu + svv * E21 * huu + svv * E22 * hvu) * scale;
					}
				}
			}
		}

		void __mix22_mat_eta(double v1, double v2, double s1, double s2, double w1, double w2)
		{
			if (_ref->__matD_eta == 0)
			{
				_ref->__matD_eta = new double[_nNode * _nNode];
				double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
				double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
				double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
				double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


				double scale = 1 / _ref->_refDv;


				for (int s = 0; s < _ref->_nNode; s++) {
					double huu = _ref->__dh[0][s];
					double huv = _ref->__dh[1][s];
					double hvu = _ref->__dh[2][s];
					double hvv = _ref->__dh[3][s];
					for (int t = 0; t < _ref->_nNode; t++) {

						double xu = 0, xv = 0, yu = 0, yv = 0;

						xu = 0;
						xv = 0;
						yu = _ref->d1[0][t];
						yv = _ref->d1[1][t];

						double duu = yu * _ref->get__gi(0, 1);
						double duv = yu * _ref->get__gi(1, 1);
						double dvu = yv * _ref->get__gi(0, 1);
						double dvv = yv * _ref->get__gi(1, 1);

						double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

						double suu = (2 * duu - tr * _ref->get__gij(0, 0));
						double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
						double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
						double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

						suu = -suu;
						suv = -suv;
						svu = suv;
						svv = -svv;
						_ref->__matD_eta[s * _nNode + t] = (suu * E11 * huv + suu * E12 * hvv + suv * E21 * huv + suv * E22 * hvv) * scale
							- (svu * E11 * huu + svu * E12 * hvu + svv * E21 * huu + svv * E22 * hvu) * scale;
					}
				}
			}
		}
		double align_mix22(double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool add)
		{
			double val = 0;
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}

			double fuu =  get__hij(0, 0);
			double fuv =  get__hij(0, 1);
			double fvv =  get__hij(1, 1);
			if (add)
			{
				
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double g11 = _ref->get__gij(0, 0);
				double g12 = _ref->get__gij(0, 1);
				double g22 = _ref->get__gij(1, 1);
				double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);

				double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;

				double huu = (2 * duu - tr * g11);
				double huv = ((duv + dvu) - tr * g12);
				double hvu = ((duv + dvu) - tr * g12);
				double hvv = (2 * dvv - tr * g22);


				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);
			}
			double fvu = fuv;
			
			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);

			
			double scale = 1 / _ref->_refDv;
			val = ((fuu) * E11 * Suv + (fuu) * E12 * Svv + (fuv) * E21 * Suv + (fuv) * E22 * Svv) * scale;
			val -= ((fvu) * E11 * Suu + (fvu) * E12 * Svu + (fvv) * E21 * Suu + (fvv) * E22 * Svu) * scale;

			return val;

		}
		
		void align_mix22_x( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio,bool add)
		{

			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}

			double fuu = get__hij(0, 0);
			double fuv = get__hij(0, 1);
			double fvv = get__hij(1, 1);

			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double g11 = _ref->get__gij(0, 0);
			double g12 = _ref->get__gij(0, 1);
			double g22 = _ref->get__gij(1, 1);
			double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);

			double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;
			if (add)
			{
			

				double huu = (2 * duu - tr * g11);
				double huv = ((duv + dvu) - tr * g12);
				double hvu = ((duv + dvu) - tr * g12);
				double hvv = (2 * dvv - tr * g22);

				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);
			}
			double fvu = fuv;

			double Suu = get__Sij(0,0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			double eu1 = 0, eu2 = 0, ev1 = 0, ev2 = 0;
			
	
			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;

			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}

			



			for (int s = 0; s < _ref->_nNode; s++)
			{
										
				double Xu = 0, Xv = 0, Yu = 0, Yv = 0;
				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;
				double _fuu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _fuv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _fvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);

				_g21 = _g12;
				if (add)
				{
					

					double _duu = xu * _ref->d1[0][s];
					double _duv = xu * _ref->d1[1][s];
					double _dvu = xv * _ref->d1[0][s];
					double _dvv = xv * _ref->d1[1][s];

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;

					double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;
					double _tr = duu * _G11 + duv * _G12 + dvu * _G12 + dvv * _G22 + _duu * _ref->get__Gij(0, 0) + _dvu * _ref->get__Gij(0, 1) + _duv * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);

					double _huu = (2 * _duu - tr * _g11 - _tr * g11);
					double _huv = (_duv + _dvu - tr * _g12 - _tr * g12);
					double _hvu = (_duv + _dvu - tr * _g12 - _tr * g12);
					double _hvv = (2 * _dvv - tr * _g22 - _tr * g22);
				
					_fuu += globalratio * (-_huu);///sc);
					_fuv += globalratio * (-_huv);// / sc);
					_fvv += globalratio * (-_hvv);///sc);
				}
				double _fvu = _fuv;

				double scale = 1 / _ref->_refDv;
				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2 * _ref->get__gij(0, 1) * _g12;
				double sqrtdv = sqrt(det);
				double _sqrtdv = 0.5  / sqrtdv * _det;
				double _scale = -1 / _ref->_refDv / _ref->_refDv * _sqrtdv;
				val = ((fuu)*E11 * _Suv + (fuu)*E12 * _Svv + (fuv)*E21 * _Suv + (fuv)*E22 * _Svv) * scale;
				val -= ((fvu)*E11 * _Suu + (fvu)*E12 * _Svu + (fvv)*E21 * _Suu + (fvv)*E22 * _Svu) * scale;
				val += ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				val += ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * _scale;
				val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * _scale;

				*ptr1 = val;
				ptr1++;
			}

		}

		void align_mix22_y(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio,bool add)
		{
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double g11 = _ref->get__gij(0, 0);
			double g12 = _ref->get__gij(0, 1);
			double g22 = _ref->get__gij(1, 1);
			double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);
			double G21 = G12;
			double fuu = get__hij(0, 0);
			double fuv = get__hij(0, 1);
			double fvv = get__hij(1, 1);
	

			double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;
			if (add)
			{
				double huu = (2 * duu - tr * g11);
				double huv = ((duv + dvu) - tr * g12);
				double hvu = ((duv + dvu) - tr * g12);
				double hvv = (2 * dvv - tr * g22);

				fuu += globalratio * (-huu/* / sc*/);
				fuv += globalratio * (-huv/* / sc*/);
				fvv += globalratio * (-hvv/* / sc*/);
			}
			double fvu = fuv;

			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			/*double eu1 = 0, eu2 = 0, ev1 = 0, ev2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				eu1 += _ref->d1[0][s] * _ref->buf_u[s];
				eu2 += _ref->d1[1][s] * _ref->buf_u[s];
				ev1 += _ref->d1[0][s] * _ref->buf_v[s];
				ev2 += _ref->d1[1][s] * _ref->buf_v[s];

			}
			double e11 = eu1 * eu1 + ev1 * ev1, e12 = eu1 * eu2 + ev1 * ev2, e22 = eu2 * eu2 + ev2 * ev2;
			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;*/
			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;

			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}

			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

				_g21 = _g12;

			

				double Xu = 0, Xv = 0, Yu = 0, Yv = 0;
				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;
				double _fuu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _fuv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _fvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				if (add)
				{
					double _duu = yu * _ref->d1[0][s];
					double _duv = yu * _ref->d1[1][s];
					double _dvu = yv * _ref->d1[0][s];
					double _dvv = yv * _ref->d1[1][s];

					double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
					double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
					double _G21 = _G12;

					double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;
					double _tr = duu * _G11 + duv * _G12 + dvu * _G12 + dvv * _G22 + _duu * _ref->get__Gij(0, 0) + _dvu * _ref->get__Gij(0, 1) + _duv * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);

					double _huu = (2 * _duu - tr * _g11 - _tr * g11);
					double _huv = (_duv + _dvu - tr * _g12 - _tr * g12);
					double _hvu = (_duv + _dvu - tr * _g12 - _tr * g12);
					double _hvv = (2 * _dvv - tr * _g22 - _tr * g22);
					//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
					//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
					//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


					_fuu += globalratio * (-_huu);// / sc);
					_fuv += globalratio * (-_huv);// / sc);
					
					_fvv += globalratio * (-_hvv);// / sc);
				}
				double _fvu = _fuv;

				double scale = 1 / _ref->_refDv;
				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2 * _ref->get__gij(0, 1) * _g12;
				double sqrtdv = sqrt(det);
				double _sqrtdv = 0.5 / sqrtdv * _det;
				double _scale = -1 / _ref->_refDv / _ref->_refDv * _sqrtdv;
				val = ((fuu)*E11 * _Suv + (fuu)*E12 * _Svv + (fuv)*E21 * _Suv + (fuv)*E22 * _Svv) * scale;
				val -= ((fvu)*E11 * _Suu + (fvu)*E12 * _Svu + (fvv)*E21 * _Suu + (fvv)*E22 * _Svu) * scale;
				val += ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				val += ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * _scale;
				val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * _scale;

				*ptr1 = val;
				ptr1++;
			}
		}
		void align_mix22_z( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio,bool add)
		{

			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double fuu = get__hij(0, 0);
			double fuv = get__hij(0, 1);
			double fvv = get__hij(1, 1);
		


			if (add)
			{
				

				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				double g11 = _ref->get__gij(0, 0);
				double g12 = _ref->get__gij(0, 1);
				double g22 = _ref->get__gij(1, 1);
				double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);
				double G21 = G12;

				double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;

				double huu = (2 * duu - tr * g11);
				double huv = ((duv + dvu) - tr * g12);
				double hvu = ((duv + dvu) - tr * g12);
				double hvv = (2 * dvv - tr * g22);
				fuu += globalratio * (-huu/* / sc*/);
				fuv += globalratio * (-huv/* / sc*/);
				fvv += globalratio * (-hvv/* / sc*/);
				
			}
			
			double fvu = fuv;
	

			double Suu = get__Sij(0,0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			
			
			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
				double _Suu = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _Suv = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _Svv = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _Svu = _Suv;


				double scale = 1 / _ref->_refDv;
				val = ((fuu)*E11 * _Suv + (fuu)*E12 * _Svv + (fuv)*E21 * _Suv + (fuv)*E22 * _Svv) * scale;
				val -= ((fvu)*E11 * _Suu + (fvu)*E12 * _Svu + (fvv)*E21 * _Suu + (fvv)*E22 * _Svu) * scale;
			
				*ptr1 = val;
				ptr1++;
			}

		}
	
		void align_mix22_phi( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{

			

		
			double Suu = get__Sij(0, 0);// -Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);// - Xvv * f1 - Yvv * f2;
			double Svu = Suv;


			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
				/*double eu1 = 0, eu2 = 0, ev1 = 0, ev2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				eu1 += _ref->d1[0][s] * _ref->buf_u[s];
				eu2 += _ref->d1[1][s] * _ref->buf_u[s];
				ev1 += _ref->d1[0][s] * _ref->buf_v[s];
				ev2 += _ref->d1[1][s] * _ref->buf_v[s];

			}
			double e11 = eu1 * eu1 + ev1 * ev1, e12 = eu1 * eu2 + ev1 * ev2, e22 = eu2 * eu2 + ev2 * ev2;

			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;*/

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double fuu = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double fuv = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double fvv = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double fvu = fuv;


				val = ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * scale;
				val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * scale;
				
				*ptr1 = val;
				ptr1++;
			}

		}

		
		void align_mix22_xi( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{


			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1   / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
			
				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = _ref->d1[0][s];
				_xv = _ref->d1[1][s];
				_yu = 0;
				_yv = 0;//

				double _duu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _duv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _dvu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _dvv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
			
				
				double _tr = _duu * _ref->get__Gij(0,0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);
	


				double huu = (2 * _duu - _tr * _ref->get__gij(0,0));
				double huv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvu = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				

				double _fuu = globalratio * (-huu) ;
				double _fuv = globalratio * (-huv) ;
				double _fvv = globalratio * (-hvv );
				double _fvu = _fuv;

			

				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix22_eta( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
		

			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = 0;
				_xv = 0;
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];

				double _duu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _duv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _dvu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _dvv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);


				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);



				double huu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double huv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvu = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				
				double _fuu = globalratio * (-huu);// / sc);
				double _fuv = globalratio * (-huv);// / sc);
				double _fvv = globalratio * (-hvv);// / sc);
				double _fvu = _fuv;



				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}

		double free_cond(double v1, double v2)
		{
			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 +  _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;
			double xu = 0, xv = 0, yu = 0, yv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}

		

			double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
			suv = 0.5 * (suv + svu);
			svu = suv;
			double tu = 0, tv = 0, wu = 0, wv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				tu += _ref->d1[0][s] * _ref->buf_phi[s];
				tv += _ref->d1[1][s] * _ref->buf_phi[s];
				wu += _ref->d1[0][s] * _ref->buf_nu[s];
				wv += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double duu = tu * _ref->get__gi(0, 0) + wu * _ref->get__gi(0, 1);
			double duv = tu * _ref->get__gi(1, 0) + wu * _ref->get__gi(1, 1);
			double dvu = tv * _ref->get__gi(0, 0) + wv * _ref->get__gi(0, 1);
			double dvv = tv * _ref->get__gi(1, 0) + wv * _ref->get__gi(1, 1);

			double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));


			suu += -huu;
			suv += -huv;
			svv += -hvv;
			double val = suu * v1 * v1 + 2 * suv * v1 * v2 + svv * v2 * v2;
			return val;


		}
		
		void free_cond_xi(double* ptr, double v1, double v2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 +  _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = _ref->d1[0][s];
				_xv = _ref->d1[1][s];
				_yu = 0;
				_yv = 0;//

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);

				double val = _suu * v1 * v1 + 2 * _suv * v1 * v2 + _svv * v2 * v2;
				*ptr1 = val;
				ptr1++;
			}
		}


		void free_cond_eta(double* ptr, double v1, double v2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = 0;
				_xv = 0;
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);

				double val = _suu * v1 * v1 + 2 * _suv * v1 * v2 + _svv * v2 * v2;
				*ptr1 = val;
				ptr1++;
			}
		}

		void free_cond_phi(double* ptr, double v1, double v2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

			

				
					double _tu = _ref->d1[0][s];
					double _tv = _ref->d1[1][s];
					double _wu = 0;
					double _wv = 0;
				
				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1));

				_suu = -_suu;
				_suv = -_suv;
				_svu = -_svu;
				_svv = -_svv;

				double val = _suu * v1 * v1 + 2 * _suv * v1 * v2 + _svv * v2 * v2;
				*ptr1 = val;
				ptr1++;
			}
		}
		void free_cond_nu(double* ptr, double v1, double v2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;

			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

			


				double _tu = 0;
				double _tv = 0;
				double _wu = _ref->d1[0][s];
				double _wv = _ref->d1[1][s];

				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1));

				_suu = -_suu;
				_suv = -_suv;
				_svu = -_svu;
				_svv = -_svv;

				double val = _suu * v1 * v1 + 2 * _suv * v1 * v2 + _svv * v2 * v2;
				*ptr1 = val;
				ptr1++;
			}
		}
		double free_cond2(double v1, double v2, double w1, double w2)
		{
			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;

			length = sqrt(_ref->get__gij(0, 0) * w1 * w1 + 2 * _ref->get__gij(0, 1) * w1 * w2 + _ref->get__gij(1, 1) * w2 * w2);
			w1 /= length; w2 /= length;

			double xu = 0, xv = 0, yu = 0, yv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}



			double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
			suv = 0.5 * (suv + svu);
			svu = suv;
			double tu = 0, tv = 0, wu = 0, wv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				tu += _ref->d1[0][s] * _ref->buf_phi[s];
				tv += _ref->d1[1][s] * _ref->buf_phi[s];
				wu += _ref->d1[0][s] * _ref->buf_nu[s];
				wv += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double duu = tu * _ref->get__gi(0, 0) + wu * _ref->get__gi(0, 1);
			double duv = tu * _ref->get__gi(1, 0) + wu * _ref->get__gi(1, 1);
			double dvu = tv * _ref->get__gi(0, 0) + wv * _ref->get__gi(0, 1);
			double dvv = tv * _ref->get__gi(1, 0) + wv * _ref->get__gi(1, 1);

			double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));


			suu += -huu;
			suv += -huv;
			svv += -hvv;
			double val = suu * v1 * w1 +  suv * (v1 * w2+v2*w1) + svv * v2 * w2;
			return val;


		}

		void free_cond2_xi(double* ptr, double v1, double v2, double w1, double w2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;
			length = sqrt(_ref->get__gij(0, 0) * w1 * w1 + 2 * _ref->get__gij(0, 1) * w1 * w2 + _ref->get__gij(1, 1) * w2 * w2);
			w1 /= length; w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = _ref->d1[0][s];
				_xv = _ref->d1[1][s];
				_yu = 0;
				_yv = 0;//

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);

				double val = _suu * v1 * w1 +  _suv * (v1 * w2+v2*w1) + _svv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}


		void free_cond2_eta(double* ptr, double v1, double v2, double w1, double w2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;
			length = sqrt(_ref->get__gij(0, 0) * w1 * w1 + 2 * _ref->get__gij(0, 1) * w1 * w2 + _ref->get__gij(1, 1) * w2 * w2);
			w1 /= length; w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = 0;
				_xv = 0;
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);

				double val = _suu * v1 * w1 +  _suv * (v1 * w2+v2*w1) + _svv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}

		void free_cond2_phi(double* ptr, double v1, double v2, double w1, double w2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;
			length = sqrt(_ref->get__gij(0, 0) * w1 * w1 + 2 * _ref->get__gij(0, 1) * w1 * w2 + _ref->get__gij(1, 1) * w2 * w2);
			w1 /= length; w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				


				double _tu = _ref->d1[0][s];
				double _tv = _ref->d1[1][s];
				double _wu = 0;
				double _wv = 0;

				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1));

				_suu = -_suu;
				_suv = -_suv;
				_svu = -_svu;
				_svv = -_svv;

				double val = _suu * v1 * w1 +_suv * (v1 * w2+v2*w1) + _svv * (v2 * w2);
				*ptr1 = val;
				ptr1++;
			}
		}
		void free_cond2_nu(double* ptr, double v1, double v2, double w1, double w2)
		{


			double length = sqrt(_ref->get__gij(0, 0) * v1 * v1 + 2 * _ref->get__gij(0, 1) * v1 * v2 + _ref->get__gij(1, 1) * v2 * v2);
			v1 /= length; v2 /= length;
			length = sqrt(_ref->get__gij(0, 0) * w1 * w1 + 2 * _ref->get__gij(0, 1) * w1 * w2 + _ref->get__gij(1, 1) * w2 * w2);
			w1 /= length; w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{



				double _tu = 0;
				double _tv = 0;
				double _wu = _ref->d1[0][s];
				double _wv = _ref->d1[1][s];

				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1));

				_suu = -_suu;
				_suv = -_suv;
				_svu = -_svu;
				_svv = -_svv;

				double val = _suu * v1 * w1 + _suv * (v1 * w2 + v2 * w1) + _svv * (v2 * w2);
				*ptr1 = val;
				ptr1++;
			}
		}



		double free_cond3(double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->get__gij(0, 0) + 2 * w1 * w2 * _ref->get__gij(0, 1) + w2 * w2 * _ref->get__gij(1, 1));
			w1 /= length;
			w2 /= length;
			double suu = 0, suv = 0, svv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				suu += _ref->d0[s] * _ref->buf_xi[s];
				suv += _ref->d0[s] * _ref->buf_eta[s];
				svv += _ref->d0[s] * _ref->buf_phi[s];
			}
			double Suu = svv / sc;
			double Suv = -suv / sc;
			double Svv = suu / sc;

			double val = Suu * v1 * w1 + Suv * (v1 * w2 + v2 * w1) + Svv * v2 * w2;
			return val;


		}

		void free_cond3_xi(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->get__gij(0, 0) + 2 * w1 * w2 * _ref->get__gij(0, 1) + w2 * w2 * _ref->get__gij(1, 1));
			w1 /= length;
			w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _suu = 0, _suv = 0, _svv = 0;
				
					_suu = _ref->d0[s];
					_suv = 0;
					_svv = 0;
					double Suu = _svv / sc;
					double Suv = -_suv / sc;
					double Svv = _suu / sc;

					double val = Suu * v1 * w1 + Suv * (v1 * w2 + v2 * w1) + Svv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}


		void free_cond3_eta(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->get__gij(0, 0) + 2 * w1 * w2 * _ref->get__gij(0, 1) + w2 * w2 * _ref->get__gij(1, 1));
			w1 /= length;
			w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _suu = 0, _suv = 0, _svv = 0;

				_suu = 0;
				_suv = _ref->d0[s];
				_svv = 0;
				double Suu = _svv / sc;
				double Suv = -_suv / sc;
				double Svv = _suu / sc;

				double val = Suu * v1 * w1 + Suv * (v1 * w2 + v2 * w1) + Svv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}

		void free_cond3_phi(double* ptr, double v1, double v2, double w1, double w2)
		{
			double length = sqrt(v1 * v1 * _ref->get__gij(0, 0) + 2 * v1 * v2 * _ref->get__gij(0, 1) + v2 * v2 * _ref->get__gij(1, 1));
			v1 /= length;
			v2 /= length;
			length = sqrt(w1 * w1 * _ref->get__gij(0, 0) + 2 * w1 * w2 * _ref->get__gij(0, 1) + w2 * w2 * _ref->get__gij(1, 1));
			w1 /= length;
			w2 /= length;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _suu = 0, _suv = 0, _svv = 0;

				_suu = 0;
				_suv = 0;
				_svv = _ref->d0[s];
				double Suu = _svv / sc;
				double Suv = -_suv / sc;
				double Svv = _suu / sc;

				double val = Suu * v1 * w1 + Suv * (v1 * w2 + v2 * w1) + Svv * v2 * w2;
				*ptr1 = val;
				ptr1++;
			}
		}
		
		double align_mix4(double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool add)
		{
			double val = 0;
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}

		
			
				double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				suv = 0.5 * (suv + svu);
				svu = suv;
				double tu = 0, tv = 0, wu = 0, wv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					tu += _ref->d1[0][s] * _ref->buf_phi[s];
					tv += _ref->d1[1][s] * _ref->buf_phi[s];
					wu += _ref->d1[0][s] * _ref->buf_nu[s];
					wv += _ref->d1[1][s] * _ref->buf_nu[s];
				}
				double duu = tu * _ref->get__gi(0, 0) + wu * _ref->get__gi(0, 1);
				double duv = tu * _ref->get__gi(1, 0) + wu * _ref->get__gi(1, 1);
				double dvu = tv * _ref->get__gi(0, 0) + wv * _ref->get__gi(0, 1);
				double dvv = tv * _ref->get__gi(1, 0) + wv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double huu = (2 * duu - tr * _ref->get__gij(0, 0));
				double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

				double Huu = (suu - huu);
				double Huv = (suv - huv) ;
				double Hvv = (svv - hvv) ;
				double Hvu = Huv;
			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			val = ((Huu)*E11 * Suv + (Huu)*E12 * Svv + (Huv)*E21 * Suv + (Huv)*E22 * Svv) * scale;
			val -= ((Hvu)*E11 * Suu + (Hvu)*E12 * Svu + (Hvv)*E21 * Suu + (Hvv)*E22 * Svu) * scale;

			return val;

		}

		
		void align_mix4_z(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool add)
		{

			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
		



				double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
				suv = 0.5 * (suv + svu);
				svu = suv;


				double tu = 0, tv = 0, wu = 0, wv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					tu += _ref->d1[0][s] * _ref->buf_phi[s];
					tv += _ref->d1[1][s] * _ref->buf_phi[s];
					wu += _ref->d1[0][s] * _ref->buf_nu[s];
					wv += _ref->d1[1][s] * _ref->buf_nu[s];
				}
				double duu = tu * _ref->get__gi(0, 0) + wu * _ref->get__gi(0, 1);
				double duv = tu * _ref->get__gi(1, 0) + wu * _ref->get__gi(1, 1);
				double dvu = tv * _ref->get__gi(0, 0) + wv * _ref->get__gi(0, 1);
				double dvv = tv * _ref->get__gi(1, 0) + wv * _ref->get__gi(1, 1);

				double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

				double huu = (2 * duu - tr * _ref->get__gij(0, 0));
				double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
				double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

				double Huu = (suu - huu);
				double Huv = (suv - huv);
				double Hvv = (svv - hvv);
				double Hvu = Huv;


		

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _Suu = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _Suv = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _Svv = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _Svu = _Suv;


				double scale = 1 / _ref->_refDv;
				val = ((Huu)*E11 * _Suv + (Huu)*E12 * _Svv + (Huv)*E21 * _Suv + (Huv)*E22 * _Svv) * scale;
				val -= ((Hvu)*E11 * _Suu + (Hvu)*E12 * _Svu + (Hvv)*E21 * _Suu + (Hvv)*E22 * _Svu) * scale;

				*ptr1 = val;
				ptr1++;
			}

		}



		void align_mix4_xi(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{


			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = _ref->d1[0][s];
				_xv = _ref->d1[1][s];
				_yu = 0;
				_yv = 0;//

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);



				double _fuu = globalratio * (_suu);
				double _fuv = globalratio * (_suv);
				double _fvv = globalratio * (_svv);
				double _fvu = _fuv;



				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix4_eta(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = 0;
				_xv = 0;
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];


				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);



				double _fuu = globalratio * (_suu);
				double _fuv = globalratio * (_suv);
				double _fvv = globalratio * (_svv);
				double _fvu = _fuv;



				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix4_phi(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _tu = _ref->d1[0][s];
				double _tv = _ref->d1[1][s];
				double _wu = 0;// _ref->d1[0][s] * _ref->buf_nu[s];
				double _wv = 0;//_ref->d1[1][s] * _ref->buf_nu[s];

				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _fuu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _fuv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _fvu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _fvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				_fuu = -_fuu;
				_fuv = -_fuv;
				_fvu = -_fvu;
				_fvv = -_fvv;


				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix4_nu(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _tu = 0;
				double _tv = 0;
				double _wu = _ref->d1[0][s];// _ref->d1[0][s] * _ref->buf_nu[s];
				double _wv = _ref->d1[1][s];//_ref->d1[1][s] * _ref->buf_nu[s];

				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _fuu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _fuv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _fvu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _fvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));


				_fuu = -_fuu;
				_fuv = -_fuv;
				_fvu = -_fvu;
				_fvv = -_fvv;
				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		double align_mix5(double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool add)
		{
			double val = 0;
			double suu = 0, suv = 0, svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				suu += _ref->d0[s] * _ref->buf_xi[s];
				suv += _ref->d0[s] * _ref->buf_eta[s];
				svv += _ref->d0[s] * _ref->buf_phi[s];
			}
			double Huu = svv/sc;
			double Huv = -suv/sc;
			double Hvv = suu/sc;
			double Hvu = Huv;


			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			val = ((Huu)*E11 * Suv + (Huu)*E12 * Svv + (Huv)*E21 * Suv + (Huv)*E22 * Svv) * scale;
			val -= ((Hvu)*E11 * Suu + (Hvu)*E12 * Svu + (Hvv)*E21 * Suu + (Hvv)*E22 * Svu) * scale;

			return val;

		}


		void align_mix5_z(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool add)
		{

			double suu = 0, suv = 0, svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				suu += _ref->d0[s] * _ref->buf_xi[s];
				suv += _ref->d0[s] * _ref->buf_eta[s];
				svv += _ref->d0[s] * _ref->buf_phi[s];
			}
			double Huu = svv / sc;
			double Huv = -suv / sc;
			double Hvv = suu / sc;
			double Hvu = Huv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _Suu = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _Suv = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _Svv = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _Svu = _Suv;


				double scale = 1 / _ref->_refDv;
				val = ((Huu)*E11 * _Suv + (Huu)*E12 * _Svv + (Huv)*E21 * _Suv + (Huv)*E22 * _Svv) * scale;
				val -= ((Hvu)*E11 * _Suu + (Hvu)*E12 * _Svu + (Hvv)*E21 * _Suu + (Hvv)*E22 * _Svu) * scale;

				*ptr1 = val;
				ptr1++;
			}

		}



		void align_mix5_xi(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{


			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				 double _suu = 0, _suv = 0, _svv = 0;

				
				 _suu = _ref->d0[s];
				 _suv = 0;//
				 _svv = 0;//

				 double _fuu = _svv / sc;
				 double _fuv = -_suv / sc;
				 double _fvv = _suu / sc;
				 double _fvu = _fuv;

				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix5_eta(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _suu = 0, _suv = 0, _svv = 0;

				_suu = 0;
				_suv = _ref->d0[s];//
				_svv = 0;//

				double _fuu = _svv / sc;
				double _fuv = -_suv / sc;
				double _fvv = _suu / sc;
				double _fvu = _fuv;

				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix5_phi(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2, double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);


			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _suu = 0, _suv = 0, _svv = 0;

				_suu = 0;
				_suv = 0;//
				_svv = _ref->d0[s];//

				double _fuu = _svv / sc;
				double _fuv = -_suv / sc;
				double _fvv = _suu / sc;
				double _fvu = _fuv;


				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}

		double align_mix32( double globalratio,bool add)
		{
			double val = 0;


			double fuu = +get__hij(0, 0);
			double fuv = +get__hij(0, 1);
			double fvv = +get__hij(1, 1);

			
			if (add)
			{
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double h11 = _ref->get__gij(0, 0);
			double h12 = _ref->get__gij(0, 1);
			double h22 = _ref->get__gij(1, 1);
			double H11 = _ref->get__Gij(0, 0), H12 = _ref->get__Gij(0, 1), H22 = _ref->get__Gij(1, 1);

			double tr = duu * H11 + duv * H12 + dvu * H12 + dvv * H22;

			double huu = (2 * duu - tr * h11);
			double huv = ((duv + dvu) - tr * h12);
			double hvu = ((duv + dvu) - tr * h12);
			double hvv = (2 * dvv - tr * h22);

			//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
			//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
			//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);

			}
			double fvu = fuv;


			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;
			
			double xu = 0, xv = 0, yu = 0, yv = 0, wu = 0, wv = 0;
			

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
				wu += _ref->d1[0][s] * _ref->buf_w[s];
				wv += _ref->d1[1][s] * _ref->buf_w[s];

			}
			double e11 = xu * xu + yu * yu + wu * wu, e12 = xu * xv + yu * yv + wu * wv, e22 = xv * xv + yv * yv + wv * wv;

			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;
			double scale = 1 / _ref->_refDv;
			val = ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * scale;
			val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * scale;

			return val;

		}
		void align_mix32_u(double* ptr, double globalratio,bool add)
		{
			double val = 0;


			double fuu = +get__hij(0, 0);
			double fuv = +get__hij(0, 1);
			double fvv = +get__hij(1, 1);
if(add)
{
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double h11 = _ref->get__gij(0, 0);
			double h12 = _ref->get__gij(0, 1);
			double h22 = _ref->get__gij(1, 1);
			double H11 = _ref->get__Gij(0, 0), H12 = _ref->get__Gij(0, 1), H22 = _ref->get__Gij(1, 1);

			double tr = duu * H11 + duv * H12 + dvu * H12 + dvv * H22;

			double huu = (2 * duu - tr * h11);
			double huv = ((duv + dvu) - tr * h12);
			double hvu = ((duv + dvu) - tr * h12);
			double hvv = (2 * dvv - tr * h22);

			//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
			//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
			//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


			
				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);

			}
			double fvu = fuv;

			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;
			double u1 = 0, u2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				u1 += _ref->d1[0][s] * _ref->buf_u[s];

				u2 += _ref->d1[1][s] * _ref->buf_u[s];

			}
			double scale = 1 / _ref->_refDv;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				
					_e11 += 2 * _ref->d1[0][s] * u1;
					_e12 += _ref->d1[0][s] *u2 + _ref->d1[1][s] * u1;
					_e22 += 2 * _ref->d1[1][s] *u2;
			
				double _E11 = _e22 * sc, _E22 = _e11 * sc, _E12 = -_e12 * sc;
				double _E21 = _E12;
				val = ((fuu)*_E11 * Suv + (fuu)*_E12 * Svv + (fuv)*_E21 * Suv + (fuv)*_E22 * Svv) * scale;
				val -= ((fvu)*_E11 * Suu + (fvu)*_E12 * Svu + (fvv)*_E21 * Suu + (fvv)*_E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}

		void  align_mix32_v(double* ptr, double globalratio,bool add)
		{
			double val = 0;
			double fuu = +get__hij(0, 0);
			double fuv = +get__hij(0, 1);
			double fvv = +get__hij(1, 1);

			if (add)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double h11 = _ref->get__gij(0, 0);
				double h12 = _ref->get__gij(0, 1);
				double h22 = _ref->get__gij(1, 1);
				double H11 = _ref->get__Gij(0, 0), H12 = _ref->get__Gij(0, 1), H22 = _ref->get__Gij(1, 1);

				double tr = duu * H11 + duv * H12 + dvu * H12 + dvv * H22;

				double huu = (2 * duu - tr * h11);
				double huv = ((duv + dvu) - tr * h12);
				double hvu = ((duv + dvu) - tr * h12);
				double hvv = (2 * dvv - tr * h22);

				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);

			}
			double fvu = fuv;

			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;
			double  v1 = 0, v2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				v1 += _ref->d1[0][s] * _ref->buf_v[s];

				v2 += _ref->d1[1][s] * _ref->buf_v[s];

			}
			double scale = 1 / _ref->_refDv;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				
					_e11 += 2 * _ref->d1[0][s] * v1;
					_e12 += _ref->d1[0][s] * v2 + _ref->d1[1][s] * v1;
					_e22 += 2 * _ref->d1[1][s] * v2;
				
				double _E11 = _e22 * sc, _E22 = _e11 * sc, _E12 = -_e12 * sc;
				double _E21 = _E12;
				val = ((fuu)*_E11 * Suv + (fuu)*_E12 * Svv + (fuv)*_E21 * Suv + (fuv)*_E22 * Svv) * scale;
				val -= ((fvu)*_E11 * Suu + (fvu)*_E12 * Svu + (fvv)*_E21 * Suu + (fvv)*_E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void  align_mix32_w(double* ptr,  double globalratio,bool add)
		{
			double val = 0;
			double fuu = +get__hij(0, 0);
			double fuv = +get__hij(0, 1);
			double fvv = +get__hij(1, 1);
			if(add)
			{
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double h11 = _ref->get__gij(0, 0);
			double h12 = _ref->get__gij(0, 1);
			double h22 = _ref->get__gij(1, 1);
			double H11 = _ref->get__Gij(0, 0), H12 = _ref->get__Gij(0, 1), H22 = _ref->get__Gij(1, 1);

			double tr = duu * H11 + duv * H12 + dvu * H12 + dvv * H22;

			double huu = (2 * duu - tr * h11);
			double huv = ((duv + dvu) - tr * h12);
			double hvu = ((duv + dvu) - tr * h12);
			double hvv = (2 * dvv - tr * h22);

			//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
			//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
			//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


		

			
				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);

			}
			double fvu = fuv;


			double Suu = get__Sij(0, 0);//-Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);//- Xvv * f1 - Yvv * f2;
			double Svu = Suv;
			double w1 = 0, w2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				w1 += _ref->d1[0][s] * _ref->buf_w[s];

				w2 += _ref->d1[1][s] * _ref->buf_w[s];

			}
			double scale = 1 / _ref->_refDv;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _e11 = 0, _e12 = 0, _e22 = 0;
				
					_e11 += 2 * _ref->d1[0][s] * w1;
					_e12 += _ref->d1[0][s] *w2 + _ref->d1[1][s] * w1;
					_e22 += 2 * _ref->d1[1][s] *w2;
			
				double _E11 = _e22 * sc, _E22 = _e11 * sc, _E12 = -_e12 * sc;
				double _E21 = _E12;
				val = ((fuu)*_E11 * Suv + (fuu)*_E12 * Svv + (fuv)*_E21 * Suv + (fuv)*_E22 * Svv) * scale;
				val -= ((fvu)*_E11 * Suu + (fvu)*_E12 * Svu + (fvv)*_E21 * Suu + (fvv)*_E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}

		void align_mix32_z(double* ptr, double globalratio,bool add)
		{
			double fuu = get__hij(0, 0);
			double fuv = get__hij(0, 1);
			double fvv = get__hij(1, 1);
			if (add)
			{
				double xu = 0, xv = 0, yu = 0, yv = 0;

				for (int s = 0; s < _ref->_nNode; s++)
				{
					xu += _ref->d1[0][s] * _ref->buf_xi[s];
					xv += _ref->d1[1][s] * _ref->buf_xi[s];
					yu += _ref->d1[0][s] * _ref->buf_eta[s];
					yv += _ref->d1[1][s] * _ref->buf_eta[s];
				}
				double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
				double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
				double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
				double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

				double g11 = _ref->get__gij(0, 0);
				double g12 = _ref->get__gij(0, 1);
				double g22 = _ref->get__gij(1, 1);
				double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);
				double G21 = G12;

				double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;

				double huu = (2 * duu - tr * g11);
				double huv = ((duv + dvu) - tr * g12);
				double hvu = ((duv + dvu) - tr * g12);
				double hvv = (2 * dvv - tr * g22);

				fuu += globalratio * (-huu /*/ sc*/);
				fuv += globalratio * (-huv /*/ sc*/);
				fvv += globalratio * (-hvv /*/ sc*/);

			}
			double fvu = fuv;



			
			double xu = 0, xv = 0, yu = 0, yv = 0, wu = 0, wv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
				wu += _ref->d1[0][s] * _ref->buf_w[s];
				wv += _ref->d1[1][s] * _ref->buf_w[s];

			}
			double e11 = xu * xu + yu * yu +wu * wu;
			double e12 = xu * xv + yu * yv + wu * wv;
			double e22 = xv * xv + yv * yv + wv * wv;

			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;
			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{








				double _Suu = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _Suv = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _Svv = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _Svu = _Suv;


				double scale = 1 / _ref->_refDv;
				val = ((fuu)*E11 * _Suv + (fuu)*E12 * _Svv + (fuv)*E21 * _Suv + (fuv)*E22 * _Svv) * scale;
				val -= ((fvu)*E11 * _Suu + (fvu)*E12 * _Svu + (fvv)*E21 * _Suu + (fvv)*E22 * _Svu) * scale;

				*ptr1 = val;
				ptr1++;
			}

		}

		void align_mix32_phi(double* ptr,  double globalratio)
		{




			double Suu = get__Sij(0, 0);// -Xuu * f1 - Yuu * f2;
			double Suv = get__Sij(0, 1);// - Xuv * f1 - Yuv * f2;
			double Svv = get__Sij(1, 1);// - Xvv * f1 - Yvv * f2;
			double Svu = Suv;


			double xu = 0, xv = 0, yu = 0, yv = 0, wu = 0, wv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
				wu += _ref->d1[0][s] * _ref->buf_w[s];
				wv += _ref->d1[1][s] * _ref->buf_w[s];

			}
			double e11 = xu * xu + yu * yu +wu * wu;
			double e12 = xu * xv + yu * yv + wu * wv;
			double e22 = xv * xv + yv * yv + wv * wv;

			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;
			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double fuu = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double fuv = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double fvv = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double fvu = fuv;


				val = ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * scale;
				val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * scale;

				*ptr1 = val;
				ptr1++;
			}

		}


		void align_mix32_xi(double* ptr,  double globalratio)
		{


			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;


			double xu = 0, xv = 0, yu = 0, yv = 0 , wu = 0, wv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
				wu += _ref->d1[0][s] * _ref->buf_w[s];
				wv += _ref->d1[1][s] * _ref->buf_w[s];

			}
			double e11 = xu * xu + yu * yu + wu * wu , e12 = xu * xv + yu * yv+ wu * wv, e22 = xv * xv + yv * yv + wv * wv;

			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;

			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = _ref->d1[0][s];
				_xv = _ref->d1[1][s];
				_yu = 0;
				_yv = 0;

				double _duu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _duv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _dvu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _dvv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);


				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);



				double huu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double huv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvu = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


				double _fuu = globalratio * (-huu);
				double _fuv = globalratio * (-huv);
				double _fvv = globalratio * (-hvv);
				double _fvu = _fuv;



				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix32_eta(double* ptr,  double globalratio)
		{



			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			
			double xu = 0, xv = 0, yu = 0, yv = 0 , wu = 0, wv = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_u[s];
				xv += _ref->d1[1][s] * _ref->buf_u[s];
				yu += _ref->d1[0][s] * _ref->buf_v[s];
				yv += _ref->d1[1][s] * _ref->buf_v[s];
				wu += _ref->d1[0][s] * _ref->buf_w[s];
				wv += _ref->d1[1][s] * _ref->buf_w[s];

			}
			double e11 = xu * xu + yu * yu  + wu * wu, e12 = xu * xv + yu * yv + wu * wv, e22 = xv * xv + yv * yv + wv * wv;

			double E11 = e22 * sc, E22 = e11 * sc, E12 = -e12 * sc, E21 = E12;

			double scale = 1 / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;


				_xu = 0;
				_xv = 0;
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];

				double _duu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _duv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _dvu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _dvv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);


				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);



				double huu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double huv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvu = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double hvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


				double _fuu = globalratio * (-huu);// / sc);
				double _fuv = globalratio * (-huv);// / sc);
				double _fvv = globalratio * (-hvv);// / sc);
				double _fvu = _fuv;



				val = ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}

		double align_metric(double v1, double v2, double s1, double s2, double w1, double w2)
		{
			double val = 0;
			
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);

			double e11 = _ref->og11;
			double e21 = _ref->og12;
			double e12 = _ref->og12;
			double e22 = _ref->og22;

			double scale = 1.0 / _ref->orefDv;
			val = (e11 * E11 * _ref->get__gij(0, 1) + e11 * E12 * _ref->get__gij(1, 1) + e12 * E21 * _ref->get__gij(0, 1) + e12 * E22 * _ref->get__gij(1, 1)) * scale;
			val -= (e21 * E11 * _ref->get__gij(0, 0) + e21 * E12 * _ref->get__gij(1, 0) + e22 * E21 * _ref->get__gij(0, 0) + e22 * E22 * _ref->get__gij(1, 0)) * scale;

			return val;

		}
		void align_metric_u(double *ptr,double v1, double v2, double s1, double s2, double w1, double w2)
		{
			double val = 0;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			double _g11 = 0, _g12 = 0, _g22 = 0, _g21 = 0;
			
			double e11 = _ref->og11;
			double e21 = _ref->og12;
			double e12 = _ref->og12;
			double e22 = _ref->og22;

			double scale = 1.0 / _ref->orefDv;
			double* ptr1 = ptr;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
			
				_g21 = _g12;
				val = (e11 * E11 * _g12 + e11 * E12 * _g22 + e12 * E21 * _g12 + e12 * E22 * _g22) * scale;
				val -= (e21 * E11 * _g11 + e21 * E12 * _g21 + e22 * E21 * _g11 + e22 * E22 * _g21) * scale;

				*ptr1 = val;
				ptr1++;
			}
		}
		void align_metric_v(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2)
		{
			double val = 0;

			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			double _g11 = 0, _g12 = 0, _g22 = 0, _g21 = 0;
			
				
			double e11 = _ref->og11;
			double e21 = _ref->og12;
			double e12 = _ref->og12;
			double e22 = _ref->og22;

			double scale = 1.0 / _ref->orefDv;
			double* ptr1 = ptr;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
		
			_g21 = _g12;
				val = (e11 * E11 * _g12 + e11 * E12 * _g22 + e12 * E21 * _g12 + e12 * E22 * _g22) * scale;
				val -= (e21 * E11 * _g11 + e21 * E12 * _g21 + e22 * E21 * _g11 + e22 * E22 * _g21) * scale;

				*ptr1 = val;
				ptr1++;
			}
		}
		double align_mix2(double v1, double v2, double s1,double s2,double w1, double w2)
		{
			double val = 0;
			/*double length = get_gij2(0, 0) * v1 * v1 + 2 * get_gij2(0, 1) * v1 * v2 + get_gij2(1, 1) * v2 * v2;
			v1 /= sqrt(length);
			v2 /= sqrt(length);
			length = get_gij2(0, 0) * s1 * s1 + 2 * get_gij2(0, 1) * s1 * s2 + get_gij2(1, 1) * s2 * s2;
			s1 /= sqrt(length);
			s2 /= sqrt(length);
			*/
			
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);

			/*double h11 = 0, h12 = 0, h22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				h11 += _ref->d2[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_phi[s];
				h12 += _ref->d2[1][s] * _ref->buf_phi[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_phi[s];
				h22 += _ref->d2[3][s] * _ref->buf_phi[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_phi[s];
			}
			double h21 = h12;*/
			double h11 = get__hij(0, 0), h12 = get__hij(0, 1), h22 = get__hij(1, 1);
			double h21 = h12;

			double scale = 1.0 / _ref->orefDv;
			val = (h11 * E11 * get__Sij(0, 1) +h11* E12 * get__Sij(1, 1) + h12 *  E21 * get__Sij(0, 1) + h12 *  E22 * get__Sij(1, 1))* scale;
			val -= (h21 * E11 * get__Sij(0, 0) + h21 *  E12 * get__Sij(1, 0) + h22 *  E21 * get__Sij(0, 0) + h22 *  E22 * get__Sij(1, 0)) * scale;

			return val;

		}

		
		void align_mix2_x( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2)
		{
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			double h11 = get__hij(0, 0), h12 = get__hij(0, 1), h22 = get__hij(1, 1);
			double h21 = h12;

			double scale = 1.0  /*/ trEij*/ / _ref->orefDv;
			double* ptr1 = ptr;
			double val = 0;
			double h1 = 0, h2 = 0;
			double S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				h1 += _ref->d1[0][s] * _ref->buf_phi[s];
				h2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{			
				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _h11 = (-_Gamma111 * h1 - _Gamma112 * h2);
				double _h12 = (-_Gamma121 * h1 - _Gamma122 * h2);
				double _h22 = (-_Gamma221 * h1 - _Gamma222 * h2);
				double _h21 = _h12;

				double _S11 = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _S12 = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _S22 = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _S21 = _S12;

				double _g11 = 2 * (_ref->d1[0][s] * _ref->get__gi(0, 0));
				double _g12 = (_ref->d1[0][s] * _ref->get__gi(1, 0)) + (_ref->d1[1][s] * _ref->get__gi(0, 0));
				double _g22 = 2 * (_ref->d1[1][s] * _ref->get__gi(1, 0));
				double _g21 = _g12;

				//_g11 = 2.0 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				//_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0,0);
				//_g22 = 2.0 * _ref->d1[1][s] * _ref->get__gi(1,0);

				//_g21 = _g12;

				//double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				//double dv = sqrt(det);
				//double _f = 0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0,1) + _g22 * _ref->get__Gij(1, 1));
				//double __dv =  _f* _dv;
				//double _scale = -1.0 / _dv / _dv * __dv;


				val = (_h11 * E11 * get__Sij(0,1) + _h11 * E12 * get__Sij(1, 1) + _h12 * E21 * get__Sij(0, 1) + _h12 * E22 * get__Sij(1, 1)) * scale;
				val -= (_h21 * E11 * get__Sij(0, 0) + _h21 * E12 * get__Sij(1, 0) + _h22 * E21 * get__Sij(0, 0) + _h22 * E22 * get__Sij(1, 0)) * scale;
				val += (h11 * E11 * _S12 + h11 * E12 * _S22 + h12 * E21 * _S12 + h12 * E22 * _S22) * scale;
				val -= (h21 * E11 * _S11 + h21 * E12 * _S21 + h22 * E21 * _S11 +h22 * E22 * _S21) * scale;
				//val += (h11 * E11 * get__Sij(0, 1) + h11 * E12 * get__Sij(1, 1) + h12 * E21 * get__Sij(0, 1) + h12 * E22 * get__Sij(1, 1)) * (_scale);
				//val -= (h21 * E11 * get__Sij(0, 0) + h21 * E12 * get__Sij(1, 0) + h22 * E21 * get__Sij(0, 0) + h22 * E22 * get__Sij(1, 0)) * (_scale);
				
				
				

				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix2_y( double* ptr, double v1, double v2, double s1, double s2, double w1, double w2)
		{
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			double h11 = get__hij(0, 0), h12 = get__hij(0, 1), h22 = get__hij(1, 1);
			double h21 = h12;


			double scale = 1.0  /*/ trEij*/ / _ref->orefDv;
			double* ptr1 = ptr;
			double val = 0;
			double h1 = 0, h2 = 0;
			double S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				h1 += _ref->d1[0][s] * _ref->buf_phi[s];
				h2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

			    double _h11 = (-_Gamma111 * h1 - _Gamma112 * h2);
				double _h12 = (-_Gamma121 * h1 - _Gamma122 * h2);
				double _h22 = (-_Gamma221 * h1 - _Gamma222 * h2);
				double _h21 = _h12;

				double _S11 = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _S12 = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _S22 = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _S21 = _S12;

				double _g11 = 2 * (_ref->d1[0][s] * _ref->get__gi(0, 1));
				double _g12 = (_ref->d1[0][s] * _ref->get__gi(1, 1)) + (_ref->d1[1][s] * _ref->get__gi(0, 1));
				double _g22 = 2 * (_ref->d1[1][s] * _ref->get__gi(1, 1));
				double _g21 = _g12;

				//double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				//double dv = sqrt(det);
				//double _f = 0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0, 1) + _g22 * _ref->get__Gij(1, 1));
				//double __dv = _f * _dv;
				//double _scale = -1.0 / _dv / _dv * __dv;


				val = (_h11 * E11 * get__Sij(0,1) + _h11 * E12 * get__Sij(1, 1) + _h12 * E21 * get__Sij(0, 1) + _h12 * E22 * get__Sij(1, 1)) * scale;
				val -= (_h21 * E11 * get__Sij(0, 0) + _h21 * E12 * get__Sij(1, 0) + _h22 * E21 * get__Sij(0, 0) + _h22 * E22 * get__Sij(1, 0)) * scale;
				val += (h11 * E11 * _S12 + h11 * E12 * _S22 + h12 * E21 * _S12 + h12 * E22 * _S22) * scale;
				val -= (h21 * E11 * _S11 + h21 * E12 * _S21 + h22 * E21 * _S11 + h22 * E22 * _S21) * scale;
				//val += (h11 * E11 * get__Sij(0, 1) + h11 * E12 * get__Sij(1, 1) + h12 * E21 * get__Sij(0, 1) + h12 * E22 * get__Sij(1, 1)) * (_scale);
				//val -= (h21 * E11 * get__Sij(0, 0) + h21 * E12 * get__Sij(1, 0) + h22 * E21 * get__Sij(0, 0) + h22 * E22 * get__Sij(1, 0)) * (_scale);




				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix2_z(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2)
		{

			/*double length = get_gij2(0, 0) * v1 * v1 + 2 * get_gij2(0, 1) * v1 * v2 + get_gij2(1, 1) * v2 * v2;
			v1 /= sqrt(length);
			v2 /= sqrt(length);
			length = get_gij2(0, 0) * s1 * s1 + 2 * get_gij2(0, 1) * s1 * s2 + get_gij2(1, 1) * s2 * s2;
			s1 /= sqrt(length);
			s2 /= sqrt(length);*/
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);
			double h11 = get__hij(0, 0), h12 = get__hij(0, 1), h22 = get__hij(1, 1);
			double h21 = h12;
			double scale = 1.0  /*/ trEij*/ / _ref->orefDv;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				/*double _g11 = 0, _g12 = 0, _g22 = 0;
				for (int t = 0; t < _ref->_nNode; t++)
				{
					_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

					_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
				}*/
				//double _g21 = _g12;
				//double ddv = 0.5 * (_g11 * get_Gij2(0, 0) + 2 * _g12 * get_Gij2(0, 1) + _g22 * get_Gij2(1, 1)) * dv;
				//double _s1 = (v1 * _g12 + v2 * _g22) / dv;
				//double _s2 = (-v1 * _g11 - v2 * _g12) / dv;
				//_s1 += -(v1 * get_gij2(0, 1) + v2 * get_gij2(1, 1)) / dv / dv * ddv;
				//_s2 += -(-v1 * get_gij2(0, 0) - v2 * get_gij2(0, 1)) / dv / dv * ddv;
				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;
				val = (h11  * E11 * _S12 + h11 *  E12 * _S22 + h12 * E21 * _S12 + h12 * E22 * _S22) * scale;
				val -= (h21 * E11 * _S11 + h21 *  E12 * _S21 + h22 * E21 * _S11 + h22 * E22 * _S21) * scale;
				
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_mix2_phi(double* ptr, double v1, double v2, double s1, double s2, double w1, double w2)
		{

			/*double length = get_gij2(0, 0) * v1 * v1 + 2 * get_gij2(0, 1) * v1 * v2 + get_gij2(1, 1) * v2 * v2;
			v1 /= sqrt(length);
			v2 /= sqrt(length);

			length = get_gij2(0, 0) * s1 * s1 + 2 * get_gij2(0, 1) * s1 * s2 + get_gij2(1, 1) * s2 * s2;
			s1 /= sqrt(length);
			s2 /= sqrt(length);*/
			double E11 = w1 * (v1 * v1) + w2 * (s1 * s1);
			double E12 = w1 * (v1 * v2) + w2 * (s1 * s2);
			double E21 = w1 * (v2 * v1) + w2 * (s2 * s1);
			double E22 = w1 * (v2 * v2) + w2 * (s2 * s2);

			double scale = 1.0  /*/ trEij*/ / _ref->orefDv;
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

			
				double _h11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _h12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _h22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _h21 = _h12;

				val = (_h11 * E11 * get__Sij(0, 1) + _h11 * E12 * get__Sij(1, 1) + _h12 * E21 * get__Sij(0, 1) + _h12 * E22 * get__Sij(1, 1)) * scale;
				val -= (_h21 * E11 * get__Sij(0, 0) + _h21 * E12 * get__Sij(1, 0) + _h22 * E21 * get__Sij(0, 0) + _h22 * E22 * get__Sij(1, 0)) * scale;
				*ptr1 = val;
				ptr1++;
			}

		}
		double singular()
		{
			double val = 0;

			double S11 = get__Sij(1, 1);
			double S12 = -get__Sij(0, 1);
			double S21 = -get__Sij(0, 1);
			double S22 = get__Sij(0, 0);

			val = S11 * get__hij(0, 0) + S12 * get__hij(1, 0) -
				S21 * get__hij(0, 1) - S22 * get__hij(1, 1);
			return val;
		}

		void singular_z( double* ptr)
		{
			double* grad = 0;
		

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = _ref->__dh[3][s];
				double _S12 = -_ref->__dh[1][s];
				double _S21 = -_ref->__dh[1][s];
				double _S22 = _ref->__dh[0][s];

				val = _S11 * get__hij(0, 0) + _S12 * get__hij(1, 0) -
					_S21 * get__hij(0, 1) - _S22 * get__hij(1, 1);
				*ptr1 = val;
				ptr1++;
			}
		}
		void singular_phi( double* ptr)
		{
			double* grad = 0;


			double* ptr1 = ptr;
			double val = 0;
			double S11 = get__Sij(1, 1);
			double S12 = -get__Sij(0, 1);
			double S21 = -get__Sij(0, 1);
			double S22 = get__Sij(0, 0);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _h11 = _ref->__dh[0][s];
				double _h12 = _ref->__dh[1][s];
				double _h21 = _ref->__dh[1][s];
				double _h22 = _ref->__dh[3][s];

				val = S11 * _h11 + S12 * _h21 -
					S21 * _h12 - S22 * _h22;
				*ptr1 = val;
				ptr1++;
			}
		}


		double singular2()
		{
			double val = 0;

			double S11 = get__Sij(1, 1);
			double S12 = -get__Sij(0, 1);
			double S21 = -get__Sij(0, 1);
			double S22 = get__Sij(0, 0);

			val = S11 * get__hij(0, 1) + S12 * get__hij(1, 1);
				
			return val;
		}

		void singular2_z( double* ptr)
		{
			double* grad = 0;


			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = _ref->__dh[3][s];
				double _S12 = -_ref->__dh[1][s];
				double _S21 = -_ref->__dh[1][s];
				double _S22 = _ref->__dh[0][s];

				val = _S11 * get__hij(0, 1) + _S12 * get__hij(1, 1);
				
				*ptr1 = val;
				ptr1++;
			}
		}
		void singular2_phi( double* ptr)
		{
			double* grad = 0;
			double S11 = get__Sij(1, 1);
			double S12 = -get__Sij(0, 1);
			double S21 = -get__Sij(0, 1);
			double S22 = get__Sij(0, 0);

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _h11 = _ref->__dh[0][s];
				double _h12 = _ref->__dh[1][s];
				double _h21 = _ref->__dh[1][s];
				double _h22 = _ref->__dh[3][s];

				val = S11 * _h12 + S12 * _h22;
				*ptr1 = val;
				ptr1++;
			}
		}

		double singular3()
		{
			double val = 0;

			double S11 = get__Sij(1, 1);
			double S12 = -get__Sij(0, 1);
			double S21 = -get__Sij(0, 1);
			double S22 = get__Sij(0, 0);

			val = S21 * get__hij(0, 0) + S22 * get__hij(1, 0);

			return val;
		}

		void singular3_z( double* ptr)
		{
			double* grad = 0;


			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = _ref->__dh[3][s];
				double _S12 = -_ref->__dh[1][s];
				double _S21 = -_ref->__dh[1][s];
				double _S22 = _ref->__dh[0][s];

				val = _S21 * get__hij(0, 0) + _S22 * get__hij(1, 0);

				*ptr1 = val;
				ptr1++;
			}
		}
		void singular3_phi( double* ptr)
		{
			double* grad = 0;


			double* ptr1 = ptr;
			double val = 0;
			double S11 = get__Sij(1, 1);
			double S12 = -get__Sij(0, 1);
			double S21 = -get__Sij(0, 1);
			double S22 = get__Sij(0, 0);
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _h11 = _ref->__dh[0][s];
				double _h12 = _ref->__dh[1][s];
				double _h21 = _ref->__dh[1][s];
				double _h22 = _ref->__dh[3][s];

				val = S21 * _h11 + S22 * _h21;
				*ptr1 = val;
				ptr1++;
			}
		}

		double singularA()
		{
			double val = 0;


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;
			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;
			val = E11 * h11 + E12 * h21;
			val -= E21 * h12 + E22 * h22;
			return val;

		}
		
		void singularA_phi( double* ptr)
		{


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _h12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _h22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _h21 = _h12;

				val = E11 * _h11 + E12 * _h21;
				val -= E21 * _h12 + E22 * _h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularA_xi( double* ptr)
		{


			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _E11 = __dsigma_22[0][s]*sc;
				double _E12 = -__dsigma_12[0][s] * sc;
				double _E22 = __dsigma_11[0][s] * sc;
				double _E21 = _E12;

				val = _E11 * h11 + _E12 * h21;
				val -= _E21 * h12 + _E22 * h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularA_eta( double* ptr)
		{


			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[1][s] * sc;
				double _E12 = -__dsigma_12[1][s] * sc;
				double _E22 = __dsigma_11[1][s] * sc;
				double _E21 = _E12;

				val = _E11 * h11 + _E12 * h21;
				val -= _E21 * h12 + _E22 * h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularA_nu( double* ptr)
		{
			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[2][s] * sc;
				double _E12 = -__dsigma_12[2][s] * sc;
				double _E22 = __dsigma_11[2][s] * sc;
				double _E21 = _E12;

				val = _E11 * h11 + _E12 * h21;
				val -= _E21 * h12 + _E22 * h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		
		double singularB()
		{
			double val = 0;


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;
			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;
			val = E21 * h11 + E22 * h21;
			
			return val;

		}

		void singularB_phi( double* ptr)
		{


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _h12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _h22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _h21 = _h12;
				val = E21 * _h11 + E22 * _h21;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularB_xi( double* ptr)
		{


			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[0][s] * sc;
				double _E12 = -__dsigma_12[0][s] * sc;
				double _E22 = __dsigma_11[0][s] * sc;
				double _E21 = _E12;

				val = _E21 * h11 + _E22 * h21;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularB_eta( double* ptr)
		{


			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[1][s] * sc;
				double _E12 = -__dsigma_12[1][s] * sc;
				double _E22 = __dsigma_11[1][s] * sc;
				double _E21 = _E12;

				val = _E21 * h11 + _E22 * h21;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularB_nu( double* ptr)
		{
			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[2][s] * sc;
				double _E12 = -__dsigma_12[2][s] * sc;
				double _E22 = __dsigma_11[2][s] * sc;
				double _E21 = _E12;

				val = _E21 * h11 + _E22 * h21;
				*ptr1 = val;
				ptr1++;
			}
		}
		double singularB2()
		{
			double val = 0;


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;
			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;
			val = E11 * h12 + E12 * h22;

			return val;

		}

		void singularB2_phi( double* ptr)
		{


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _h12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _h22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _h21 = _h12;
				val = E11 * _h12 + E12 * _h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularB2_xi( double* ptr)
		{


			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[0][s] * sc;
				double _E12 = -__dsigma_12[0][s] * sc;
				double _E22 = __dsigma_11[0][s] * sc;
				double _E21 = _E12;

				val = _E11 * h12 + _E12 * h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularB2_eta( double* ptr)
		{


			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[1][s] * sc;
				double _E12 = -__dsigma_12[1][s] * sc;
				double _E22 = __dsigma_11[1][s] * sc;
				double _E21 = _E12;

				val = _E11 * h12 + _E12 * h22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularB2_nu( double* ptr)
		{
			double h11 = get__hij(0, 0);
			double h22 = get__hij(1, 1);
			double h12 = get__hij(0, 1);
			double h21 = h12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[2][s] * sc;
				double _E12 = -__dsigma_12[2][s] * sc;
				double _E22 = __dsigma_11[2][s] * sc;
				double _E21 = _E12;

				val = _E11 * h12 + _E12 * h22;
				*ptr1 = val;
				ptr1++;
			}
		}

		double singularC()
		{
			double val = 0;


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;
			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;
			val = E11 * S11 + E12 * S21;
			val -= E21 * S12 + E22 * S22;
			return val;

		}

		void singularC_z( double* ptr)
		{


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;

				val = E11 * _S11 + E12 * _S21;
				val -= E21 * _S12 + E22 * _S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularC_xi( double* ptr)
		{


			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[0][s] * sc;
				double _E12 = -__dsigma_12[0][s] * sc;
				double _E22 = __dsigma_11[0][s] * sc;
				double _E21 = _E12;

				val = _E11 * S11 + _E12 * S21;
				val -= _E21 * S12 + _E22 * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularC_eta( double* ptr)
		{


			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[1][s] * sc;
				double _E12 = -__dsigma_12[1][s] * sc;
				double _E22 = __dsigma_11[1][s] * sc;
				double _E21 = _E12;

				val = _E11 * S11 + _E12 * S21;
				val -= _E21 * S12 + _E22 * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularC_nu( double* ptr)
		{
			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[2][s] * sc;
				double _E12 = -__dsigma_12[2][s] * sc;
				double _E22 = __dsigma_11[2][s] * sc;
				double _E21 = _E12;

				val = _E11 * S11 + _E12 * S21;
				val -= _E21 * S12 + _E22 * S22;
				*ptr1 = val;
				ptr1++;
			}
		}

		double singularD()
		{
			double val = 0;


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;
			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;
			val = E21 * S11 + E22 * S21;

			return val;

		}

		void singularD_z( double* ptr)
		{


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;
				val = E21 * _S11 + E22 * _S21;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularD_xi( double* ptr)
		{


			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[0][s] * sc;
				double _E12 = -__dsigma_12[0][s] * sc;
				double _E22 = __dsigma_11[0][s] * sc;
				double _E21 = _E12;

				val = _E21 * S11 + _E22 * S21;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularD_eta( double* ptr)
		{


			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[1][s] * sc;
				double _E12 = -__dsigma_12[1][s] * sc;
				double _E22 = __dsigma_11[1][s] * sc;
				double _E21 = _E12;

				val = _E21 * S11 + _E22 * S21;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularD_nu( double* ptr)
		{
			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[2][s] * sc;
				double _E12 = -__dsigma_12[2][s] * sc;
				double _E22 = __dsigma_11[2][s] * sc;
				double _E21 = _E12;

				val = _E21 * S11 + _E22 * S21;
				*ptr1 = val;
				ptr1++;
			}
		}

		double singularD2()
		{
			double val = 0;


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;
			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;
			val = E11 * S12 + E12 * S22;

			return val;

		}

		void singularD2_z( double* ptr)
		{


			double E11 = get_Eij(0, 0);
			double E22 = get_Eij(1, 1);
			double E12 = get_Eij(0, 1);
			double E21 = E12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;
				val = E11 * _S12 + E12 * _S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularD2_xi( double* ptr)
		{


			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[0][s] * sc;
				double _E12 = -__dsigma_12[0][s] * sc;
				double _E22 = __dsigma_11[0][s] * sc;
				double _E21 = _E12;

				val = _E11 * S12 + _E12 * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularD2_eta( double* ptr)
		{


			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[1][s] * sc;
				double _E12 = -__dsigma_12[1][s] * sc;
				double _E22 = __dsigma_11[1][s] * sc;
				double _E21 = _E12;

				val = _E11 * S12 + _E12 * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		void singularD2_nu( double* ptr)
		{
			double S11 = get__Sij(0, 0);
			double S22 = get__Sij(1, 1);
			double S12 = get__Sij(0, 1);
			double S21 = S12;

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _E11 = __dsigma_22[2][s] * sc;
				double _E12 = -__dsigma_12[2][s] * sc;
				double _E22 = __dsigma_11[2][s] * sc;
				double _E21 = _E12;

				val = _E11 * S12 + _E12 * S22;
				*ptr1 = val;
				ptr1++;
			}
		}
		double trace()
		{
			double val = get_Eij(0, 0) * _ref->og11 + 2 * get_Eij(0, 1) * _ref->og12 + get_Eij(1, 1) * _ref->og22;

			return val;
		}
		void trace_xi(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f22 = __dsigma_22[0][s];
				double _f21 = _f12;

				double _E11 = _f22 * sc;
				double _E22 = _f11 * sc;
				double _E12 = -_f12 * sc;
				double _E21 = -_f12 * sc;
				double val = _E11 * _ref->og11 + 2 * _E12 * _ref->og12 + _E22 * _ref->og22;

				*ptr1 = val;
				ptr1++;
			}

		}
		void trace_eta(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _f11 = __dsigma_11[1][s];
				double _f12 = __dsigma_12[1][s];
				double _f22 = __dsigma_22[1][s];
				double _f21 = _f12;

				double _E11 = _f22 * sc;
				double _E22 = _f11 * sc;
				double _E12 = -_f12 * sc;
				double _E21 = -_f12 * sc;
				double val = _E11 * _ref->og11 + 2 * _E12 * _ref->og12 + _E22 * _ref->og22;

				*ptr1 = val;
				ptr1++;
			}

		}
		void trace_nu(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _f11 = __dsigma_11[2][s];
				double _f12 = __dsigma_12[2][s];
				double _f22 = __dsigma_22[2][s];
				double _f21 = _f12;

				double _E11 = _f22 * sc;
				double _E22 = _f11 * sc;
				double _E12 = -_f12 * sc;
				double _E21 = -_f12 * sc;
				double val = _E11 * _ref->og11 + 2 * _E12 * _ref->og12 + _E22 * _ref->og22;

				*ptr1 = val;
				ptr1++;
			}

		}
		double align_sigma()
		{
			double val = 0;

	
		
			double scale = 1.0  /*/ trEij*/ / _ref->_refDv;
			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			val = (get__hij(0,0) * get_Eij(0,0) * get__Sij(0, 1) + get__hij(0, 0) * get_Eij(0, 1) * get__Sij(1, 1) + get__hij(0, 1) * get_Eij(1, 0) * get__Sij(0, 1) + get__hij(0, 1) *get_Eij(1, 1) * get__Sij(1, 1)) * scale;
			val -= (get__hij(1, 0) * get_Eij(0, 0)  * get__Sij(0,0) + get__hij(1, 0) * get_Eij(0, 1) * get__Sij(1, 0) + get__hij(1, 1) * get_Eij(1, 0) *  get__Sij(0, 0) + get__hij(1, 1) * get_Eij(1, 1) * get__Sij(1, 0)) * scale;
			//val /= tr;
			
			return val;

		}
		
		
		void align_sigma_z(double *ptr)
		{

			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			double S11 = 0;

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double* ptr1 = ptr;
			double val=0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;
				val = (get__hij(0, 0) * get_Eij(0,0) * _S12 + get__hij(0, 0) * get_Eij(0, 1) * _S22 + get__hij(0, 1) * get_Eij(1, 0) * _S12 + get__hij(0, 1) * get_Eij(1, 1) * _S22) * scale;
				val -= (get__hij(1, 0)  * get_Eij(0, 0) * _S11 + get__hij(1, 0) * get_Eij(0, 1)* _S21 + get__hij(1, 1) * get_Eij(1, 0) * _S11 + get__hij(1, 1) * get_Eij(1, 1) * _S21)* scale;
				//val /= tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_sigma_phi(double* ptr)
		{


			double S11 = 0;

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][s] ;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][s] ;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;
				double _h21 = _h12;
				val = (_h11 * get_Eij(0, 0) * get__Sij(0, 1) + _h11 * get_Eij(0, 1) *get__Sij(1, 1) + _h12 * get_Eij(1, 0) * get__Sij(0, 1) + _h12 * get_Eij(1, 1) * get__Sij(1, 1)) * scale;
				val -= (_h21 * get_Eij(0, 0) * get__Sij(0, 0) + _h21 * get_Eij(0, 1) *get__Sij(1, 0) + _h22 * get_Eij(1, 0)  * get__Sij(0, 0) + _h22 * get_Eij(1, 1)  * get__Sij(1, 0))* scale;
				//val /= tr;
				*ptr1 = val;
				ptr1++;
			}

		
}
		void align_sigma_x(double* ptr)
		{
			double fuu = get__hij(0, 0);
			double fuv = get__hij(0, 1);
			double fvu = get__hij(1, 0);
			double fvv = get__hij(1, 1);

			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = get_Eij(0, 0);
			double E12 = get_Eij(0, 1);
			double E21 = get_Eij(1, 0);
			double E22 = get_Eij(1, 1);

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			double* ptr1 = ptr;
			double val = 0;

			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}





			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);

				_g21 = _g12;


				double Xu = 0, Xv = 0, Yu = 0, Yv = 0;
				double _Gamma111 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(0, 0);
				double _Gamma112 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(1, 0);
				double _Gamma121 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(0, 0);
				double _Gamma122 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(1, 0);
				double _Gamma221 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(0, 0);
				double _Gamma222 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;

				double _fuu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _fuv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _fvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _fvu = _fuv;

				double scale = 1 / _ref->_refDv;
				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2 * _ref->get__gij(0, 1) * _g12;
				double sqrtdv = sqrt(det);
				double _sqrtdv = 0.5 / sqrtdv * _det;
				double _scale = -1 / _ref->_refDv / _ref->_refDv * _sqrtdv;
				//double _tr = get_Eij(0, 0) * _g11 + 2 * get_Eij(0, 1) * _g12 + get_Eij(1, 1) * _g22;

				val = ((fuu)*E11 * _Suv + (fuu)*E12 * _Svv + (fuv)*E21 * _Suv + (fuv)*E22 * _Svv) * scale;
				val -= ((fvu)*E11 * _Suu + (fvu)*E12 * _Svu + (fvv)*E21 * _Suu + (fvv)*E22 * _Svu) * scale;
				val += ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				val += ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * _scale;
				val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * _scale;
			
				double _sc = -1 / det / det * _det;

				double _E11 = get_eij(0,0) * _sc;
				double _E22 = get_eij(1, 1) * _sc;
				double _E12 = get_eij(0, 1) * _sc;
				double _E21 = get_eij(1, 0) * _sc;

				val += ((fuu)*_E11 * Suv + (fuu)*_E12 * Svv + (fuv)*_E21 * Suv + (fuv)*_E22 * Svv) * scale;
				val -= ((fvu)*_E11 * Suu + (fvu)*_E12 * Svu + (fvv)*_E21 * Suu + (fvv)*_E22 * Svu) * scale;

				*ptr1 = val;
				ptr1++;
			}

		}

		void align_sigma_y(double* ptr)
		{
			double fuu = get__hij(0, 0);
			double fuv = get__hij(0, 1);
			double fvu = get__hij(1, 0);
			double fvv = get__hij(1, 1);

			double Suu = get__Sij(0, 0);// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
			double Suv = get__Sij(0, 1);// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
			double Svv = get__Sij(1, 1);// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
			double Svu = Suv;

			double E11 = get_Eij(0, 0);
			double E12 = get_Eij(0, 1);
			double E21 = get_Eij(1, 0);
			double E22 = get_Eij(1, 1);

			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			double* ptr1 = ptr;
			double val = 0;

			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}





			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);

				_g21 = _g12;


				double Xu = 0, Xv = 0, Yu = 0, Yv = 0;
				double _Gamma111 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(0, 1);
				double _Gamma112 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(1, 1);
				double _Gamma121 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(0, 1);
				double _Gamma122 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(1, 1);
				double _Gamma221 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(0, 1);
				double _Gamma222 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;

				double _fuu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _fuv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _fvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _fvu = _fuv;

				double scale = 1 / _ref->_refDv;
				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2 * _ref->get__gij(0, 1) * _g12;
				double sqrtdv = sqrt(det);
				double _sqrtdv = 0.5 / sqrtdv * _det;
				double _scale = -1 / _ref->_refDv / _ref->_refDv * _sqrtdv;
				
				val = ((fuu)*E11 * _Suv + (fuu)*E12 * _Svv + (fuv)*E21 * _Suv + (fuv)*E22 * _Svv) * scale ;
				val -= ((fvu)*E11 * _Suu + (fvu)*E12 * _Svu + (fvv)*E21 * _Suu + (fvv)*E22 * _Svu) * scale;
				val += ((_fuu)*E11 * Suv + (_fuu)*E12 * Svv + (_fuv)*E21 * Suv + (_fuv)*E22 * Svv) * scale;
				val -= ((_fvu)*E11 * Suu + (_fvu)*E12 * Svu + (_fvv)*E21 * Suu + (_fvv)*E22 * Svu) * scale;
				val += ((fuu)*E11 * Suv + (fuu)*E12 * Svv + (fuv)*E21 * Suv + (fuv)*E22 * Svv) * _scale;
				val -= ((fvu)*E11 * Suu + (fvu)*E12 * Svu + (fvv)*E21 * Suu + (fvv)*E22 * Svu) * _scale;
				double _sc = -1 / det / det * _det;

				double _E11 = get_eij(0, 0) * _sc;
				double _E22 = get_eij(1, 1) * _sc;
				double _E12 = get_eij(0, 1) * _sc;
				double _E21 = get_eij(1, 0) * _sc;

				val += ((fuu)*_E11 * Suv + (fuu)*_E12 * Svv + (fuv)*_E21 * Suv + (fuv)*_E22 * Svv) * scale;
				val -= ((fvu)*_E11 * Suu + (fvu)*_E12 * Svu + (fvv)*_E21 * Suu + (fvv)*_E22 * Svu) * scale;

				*ptr1 = val;
				ptr1++;
			}
		}
		void align_sigma_xi( double* ptr)
		{			
			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _f11 = __dsigma_11[0][s];
				double _f12 = __dsigma_12[0][s];
				double _f22 = __dsigma_22[0][s];
				double _f21 = _f12;

				double _E11 = _f22*sc;
				double _E22 = _f11*sc;
				double _E12 = -_f12*sc;
				double _E21 = -_f12*sc;

				val = (get__hij(0,0) * _E11 * get__Sij(0, 1) + get__hij(0, 0) * _E12 * get__Sij(1, 1) + get__hij(0, 1) * _E21 * get__Sij(0, 1) + get__hij(0, 1) * _E22 * get__Sij(1, 1)) * scale;
				val -= (get__hij(1, 0) * _E11 * get__Sij(0, 0) + get__hij(1, 0) * _E12 * get__Sij(1, 0) + get__hij(1, 1) * _E21 * get__Sij(0, 0) + get__hij(1, 1) * _E22 * get__Sij(1, 0)) * scale;
				//val /= tr;
				//double _tr =_E11 * _ref->get__gij(0, 0) + 2 * _E12 * _ref->get__gij(0, 1) + _E22 * _ref->get__gij(1, 1);

				//val +=- (get__hij(0, 0) * get_Eij(0, 0) * get__Sij(0, 1) + get__hij(0, 0) * get_Eij(0, 1) * get__Sij(1, 1) + get__hij(0, 1) * get_Eij(1, 0) * get__Sij(0, 1) + get__hij(0, 1) * get_Eij(1, 1) * get__Sij(1, 1)) * scale/tr/tr*_tr;
				//val -= -(get__hij(1, 0) * get_Eij(0, 0) * get__Sij(0, 0) + get__hij(1, 0) * get_Eij(0, 1) * get__Sij(1, 0) + get__hij(1, 1) * get_Eij(1, 0) * get__Sij(0, 0) + get__hij(1, 1) * get_Eij(1, 1) * get__Sij(1, 0)) * scale/tr/tr*_tr;
				*ptr1 = val;
				ptr1++;
			}
		}
		void align_sigma_eta( double* ptr)
		{


			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _f11 = __dsigma_11[1][s];
				double _f12 = __dsigma_12[1][s];
				double _f22 = __dsigma_22[1][s];
				double _E11 = _f22 * sc;
				double _E22 = _f11 * sc;
				double _E12 = -_f12*sc;
				double _E21 = -_f12*sc;


			
				val = (get__hij(0, 0) * _E11 * get__Sij(0, 1) + get__hij(0, 0) * _E12 * get__Sij(1, 1) + get__hij(0, 1) * _E21 * get__Sij(0, 1) + get__hij(0, 1) * _E22 * get__Sij(1, 1)) * scale;
				val -= (get__hij(1, 0) * _E11 * get__Sij(0, 0) + get__hij(1, 0) * _E12 * get__Sij(1, 0) + get__hij(1, 1) * _E21 * get__Sij(0, 0) + get__hij(1, 1) * _E22 * get__Sij(1, 0)) * scale;
				//val /= tr;
				//double _tr = _E11 * _ref->get__gij(0, 0) + 2 * _E12 * _ref->get__gij(0, 1) + _E22 * _ref->get__gij(1, 1);

				//val += -(get__hij(0, 0) * get_Eij(0, 0) * get__Sij(0, 1) + get__hij(0, 0) * get_Eij(0, 1) * get__Sij(1, 1) + get__hij(0, 1) * get_Eij(1, 0) * get__Sij(0, 1) + get__hij(0, 1) * get_Eij(1, 1) * get__Sij(1, 1)) * scale / tr / tr * _tr;
				//val -= -(get__hij(1, 0) * get_Eij(0, 0) * get__Sij(0, 0) + get__hij(1, 0) * get_Eij(0, 1) * get__Sij(1, 0) + get__hij(1, 1) * get_Eij(1, 0) * get__Sij(0, 0) + get__hij(1, 1) * get_Eij(1, 1) * get__Sij(1, 0)) * scale / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}

		}
		void align_sigma_nu( double* ptr)
		{


			double scale = 1  /*/ trEij*/ / _ref->_refDv;
			double val = 0;
			double* ptr1 = ptr;
			//double tr = get_Eij(0, 0) * _ref->get__gij(0, 0) + 2 * get_Eij(0, 1) * _ref->get__gij(0, 1) + get_Eij(1, 1) * _ref->get__gij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _f11 = __dsigma_11[2][s];
				double _f12 = __dsigma_12[2][s];
				double _f22 = __dsigma_22[2][s];
				double _E11 = _f22 * sc;
				double _E22 = _f11 * sc;
				double _E12 = -_f12 * sc;
				double _E21 = -_f12 * sc;



				val = (get__hij(0, 0) * _E11 * get__Sij(0, 1) + get__hij(0, 0) * _E12 * get__Sij(1, 1) + get__hij(0, 1) * _E21 * get__Sij(0, 1) + get__hij(0, 1) * _E22 * get__Sij(1, 1)) * scale;
				val -= (get__hij(1, 0) * _E11 * get__Sij(0, 0) + get__hij(1, 0) * _E12 * get__Sij(1, 0) + get__hij(1, 1) * _E21 * get__Sij(0, 0) + get__hij(1, 1) * _E22 * get__Sij(1, 0)) * scale;
				//val /= tr;
				//double _tr = _E11 * _ref->get__gij(0, 0) + 2 * _E12 * _ref->get__gij(0, 1) + _E22 * _ref->get__gij(1, 1);

				//val += -(get__hij(0, 0) * get_Eij(0, 0) * get__Sij(0, 1) + get__hij(0, 0) * get_Eij(0, 1) * get__Sij(1, 1) + get__hij(0, 1) * get_Eij(1, 0) * get__Sij(0, 1) + get__hij(0, 1) * get_Eij(1, 1) * get__Sij(1, 1)) * scale / tr / tr * _tr;
				//val -= -(get__hij(1, 0) * get_Eij(0, 0) * get__Sij(0, 0) + get__hij(1, 0) * get_Eij(0, 1) * get__Sij(1, 0) + get__hij(1, 1) * get_Eij(1, 0) * get__Sij(0, 0) + get__hij(1, 1) * get_Eij(1, 1) * get__Sij(1, 0)) * scale / tr / tr * _tr;
				*ptr1 = val;
				ptr1++;
			}


		}

		void __bodyF2_mat_phi()
		{
			if (_ref->__matF_phi == 0)
			{
				_ref->__matF_phi = new double[_nNode * _nNode];

				for (int s = 0; s < _ref->_nNode; s++) {
					double huu = _ref->__dh[0][s];
					double huv = _ref->__dh[1][s];
					double hvu = _ref->__dh[2][s];
					double hvv = _ref->__dh[3][s];
					for (int t = 0; t < _ref->_nNode; t++) {

						double suu = _ref->__dh[0][t];
						double suv = _ref->__dh[1][t];
						double svu = _ref->__dh[2][t];
						double svv = _ref->__dh[3][t];

						double Suu = svv * sc;
						double Suv = -suv * sc;
						double Svu = -suv * sc;
						double Svv = suu * sc;
						_ref->__matF_phi[s * _nNode + t] = Suu * huu + 2 * Suv * huv + Svv * hvv;
					}
				}
			}


		}
	
		void __bodyF2_mat_xi()
		{
			if (_ref->__matF_xi == 0)
			{
				_ref->__matF_xi = new double[_nNode * _nNode];

				for (int s = 0; s < _ref->_nNode; s++) {
					double huu = _ref->__dh[0][s];
					double huv = _ref->__dh[1][s];
					double hvu = _ref->__dh[2][s];
					double hvv = _ref->__dh[3][s];
					for (int t = 0; t < _ref->_nNode; t++) {

						double xu = 0, xv = 0, yu = 0, yv = 0;


						xu = _ref->d1[0][t];
						xv = _ref->d1[1][t];
						yu = 0;
						yv = 0;

						double duu = xu * _ref->get__gi(0, 0);
						double duv = xu * _ref->get__gi(1, 0);
						double dvu = xv * _ref->get__gi(0, 0);
						double dvv = xv * _ref->get__gi(1, 0);

						double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

						double suu = (2 * duu - tr * _ref->get__gij(0, 0));
						double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
						double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
						double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

						double Suu = -svv * sc;
						double Suv = suv * sc;
						double Svu = suv * sc;
						double Svv = -suu * sc;
						_ref->__matF_xi[s * _nNode + t] = Suu * huu + 2 * Suv * huv + Svv * hvv;
					}
				}
			}
		}

		void __bodyF2_mat_eta()
		{
			if (_ref->__matF_eta == 0)
			{
				_ref->__matF_eta = new double[_nNode * _nNode];

				for (int s = 0; s < _ref->_nNode; s++) {
					double huu = _ref->__dh[0][s];
					double huv = _ref->__dh[1][s];
					double hvu = _ref->__dh[2][s];
					double hvv = _ref->__dh[3][s];
					for (int t = 0; t < _ref->_nNode; t++) {

						double xu = 0, xv = 0, yu = 0, yv = 0;

						xu = 0;
						xv = 0;
						yu = _ref->d1[0][t];
						yv = _ref->d1[1][t];

						double duu = yu * _ref->get__gi(0, 1);
						double duv = yu * _ref->get__gi(1, 1);
						double dvu = yv * _ref->get__gi(0, 1);
						double dvv = yv * _ref->get__gi(1, 1);

						double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

						double suu = (2 * duu - tr * _ref->get__gij(0, 0));
						double suv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
						double svu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
						double svv = (2 * dvv - tr * _ref->get__gij(1, 1));

						double Suu = -svv * sc;
						double Suv = suv * sc;
						double Svu = suv * sc;
						double Svv = -suu * sc;
						_ref->__matF_eta[s * _nNode + t] = Suu * huu + 2 * Suv * huv + Svv * hvv;
					}
				}
			}
		}
		double __bodyF5(double globalratio)
		{
			double val = 0;

			double suu = 0, suv = 0, svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				suu += _ref->d0[s] * _ref->buf_xi[s];
				suv += _ref->d0[s] * _ref->buf_eta[s];
				svv += _ref->d0[s] * _ref->buf_phi[s];

			}
			double Huu = suu;
			double Huv = suv;
			double Hvv = svv;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);

			val = Hvv * Svv + 2 * Huv * Suv + Huu * Suu;

			/*if (accurate_area)
				val -= load * dv / _ref->_refDv;
			else
				val -= load;*/
			return val;

		}
		double __bodyF5_rhs(double load, bool accurate_area)
		{
			double val = 0;

			if (accurate_area)
				val = load * dv / _ref->_refDv;
			else
				val = load;
			return val;

		}

		void __bodyF5_z(double* ptr, double load, bool accurate_area, double globalratio)
		{

			double* ptr1 = ptr;
			double val = 0;
			double suu = 0, suv = 0, svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				suu += _ref->d0[s] * _ref->buf_xi[s];
				suv += _ref->d0[s] * _ref->buf_nu[s];
				svv += _ref->d0[s] * _ref->buf_phi[s];

			}
			double Huu = suu;
			double Huv = suv;
			double Hvv = svv;

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];
				double _S12 = _ref->__dh[1][s];
				double _S22 = _ref->__dh[3][s];
				double _S21 = _S12;
				val = Hvv * _S22 + 2 * Huv * _S12 + Huu * _S11;
				if (accurate_area) {
					double _g11 = 0, _g12 = 0, _g22 = 0;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					}

					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					val += -load * ddv / _ref->_refDv;
				}
				*ptr1 = val;
				ptr1++;
			}

		}

		void __bodyF5_xi(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				double suu = 0, suv = 0, svv = 0;

				
					suu = _ref->d0[s] ;
					suv = 0;//
					svv = 0;//

				
					double _Huu = suu;
					double _Huv = suv;
					double _Hvv = svv;
				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF5_eta(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double suu = 0, suv = 0, svv = 0;


				suu = 0;
				suv = _ref->d0[s] ;//
				svv = 0;//


				double _Huu = suu;
				double _Huv = suv;
				double _Hvv = svv;

				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF5_phi(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double suu = 0, suv = 0, svv = 0;


				suu = 0;
				suv = 0;//
				svv = _ref->d0[s] ;//

				double _Huu = suu;
				double _Huv = suv;
				double _Hvv = svv;

				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
	
		double __bodyF4(double globalratio)
		{
			double val = 0;

			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
			suv = 0.5 * (suv + svu);
			svu = suv;

			double tu = 0, tv = 0, wu = 0, wv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				tu += _ref->d1[0][s] * _ref->buf_phi[s];
				tv += _ref->d1[1][s] * _ref->buf_phi[s];
				wu += _ref->d1[0][s] * _ref->buf_nu[s];
				wv += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double duu = tu * _ref->get__gi(0, 0) + wu * _ref->get__gi(0, 1);
			double duv = tu * _ref->get__gi(1, 0) + wu * _ref->get__gi(1, 1);
			double dvu = tv * _ref->get__gi(0, 0) + wv * _ref->get__gi(0, 1);
			double dvv = tv * _ref->get__gi(1, 0) + wv * _ref->get__gi(1, 1);

			double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

			double Huu = (svv - hvv) * sc;
			double Huv = -(suv - huv) * sc;
			double Hvv = (suu - huu) * sc;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);

			val = Hvv * Svv + 2 * Huv * Suv + Huu * Suu;

			/*if (accurate_area)
				val -= load * dv / _ref->_refDv;
			else
				val -= load;*/
			return val;

		}
		double __bodyF4_rhs(double load, bool accurate_area)
		{
			double val = 0;

			if (accurate_area)
				val = load * dv / _ref->_refDv;
			else
				val = load;
			return val;

		}
		
		void __bodyF4_z(double* ptr, double load, bool accurate_area, double globalratio)
		{

			double* ptr1 = ptr;
			double val = 0;
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double suu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double suv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double svu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double svv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);
			suv = 0.5 * (suv + svu);
			svu = suv;

			double tu = 0, tv = 0, wu = 0, wv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				tu += _ref->d1[0][s] * _ref->buf_phi[s];
				tv += _ref->d1[1][s] * _ref->buf_phi[s];
				wu += _ref->d1[0][s] * _ref->buf_nu[s];
				wv += _ref->d1[1][s] * _ref->buf_nu[s];
			}
			double duu = tu * _ref->get__gi(0, 0) + wu * _ref->get__gi(0, 1);
			double duv = tu * _ref->get__gi(1, 0) + wu * _ref->get__gi(1, 1);
			double dvu = tv * _ref->get__gi(0, 0) + wv * _ref->get__gi(0, 1);
			double dvv = tv * _ref->get__gi(1, 0) + wv * _ref->get__gi(1, 1);

			double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

			double Huu = (svv - hvv) * sc;
			double Huv = -(suv - huv) * sc;
			double Hvv = (suu - huu) * sc;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];
				double _S12 = _ref->__dh[1][s];
				double _S22 = _ref->__dh[3][s];
				double _S21 = _S12;
				val = Hvv * _S22 + 2 * Huv * _S12 + Huu * _S11;
				if (accurate_area) {
					double _g11 = 0, _g12 = 0, _g22 = 0;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					}

					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					val += -load * ddv / _ref->_refDv;
				}
				*ptr1 = val;
				ptr1++;
			}

		}

		void __bodyF4_xi(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;
				_xu = _ref->d1[0][s];
				_xv = _ref->d1[1][s];
				_yu = 0;
				_yv = 0;

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);
				_svu = _suv;

				double _Huu = (_svv ) * sc;
				double _Huv = -(_suv ) * sc;
				double _Hvv = (_suu ) * sc;

				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF4_eta(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;
				_xu = 0;// _ref->d1[0][s];
				_xv = 0;//_ref->d1[1][s];
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];

				double _suu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _suv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _svu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _svv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);
				_suv = 0.5 * (_suv + _svu);
				_svu = _suv;

				double _Huu = (_svv ) * sc;
				double _Huv = -(_suv ) * sc;
				double _Hvv = (_suu ) * sc;

				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF4_phi(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _tu = 0, _tv = 0, _wu = 0, _wv = 0;

				_tu = _ref->d1[0][s];
				_tv = _ref->d1[1][s];
				_wu = 0;// _ref->d1[0][s] * _ref->buf_nu[s];
				_wv = 0;//_ref->d1[1][s] * _ref->buf_nu[s];
				
				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				double _Huu = (-_svv)*sc;
				double _Huv = -(-_suv)*sc;
				double _Hvv = (-_suu)*sc;
				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF4_nu(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _tu = 0, _tv = 0, _wu = 0, _wv = 0;


				_tu = 0;// _ref->d1[0][s] * _ref->buf_phi[s];
				_tv = 0;//_ref->d1[1][s] * _ref->buf_phi[s];
				_wu = _ref->d1[0][s];
				_wv = _ref->d1[1][s];

				double _duu = _tu * _ref->get__gi(0, 0) + _wu * _ref->get__gi(0, 1);
				double _duv = _tu * _ref->get__gi(1, 0) + _wu * _ref->get__gi(1, 1);
				double _dvu = _tv * _ref->get__gi(0, 0) + _wv * _ref->get__gi(0, 1);
				double _dvv = _tv * _ref->get__gi(1, 0) + _wv * _ref->get__gi(1, 1);

				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _suu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _suv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				double _svu = ((_duv + _dvu) - _tr * _ref->get__gij(1, 0));
				double _svv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				double _Huu = (-_svv)*sc;
				double _Huv = -(-_suv)*sc;
				double _Hvv = (-_suu)*sc;
				val = (Suu * _Huu + 2 * Suv * _Huv + Svv * _Hvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		double __bodyF2( double globalratio)
		{
			double val = 0;

			double xu = 0, xv = 0, yu = 0, yv = 0;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			
			
			double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

			double huu = (2*duu - tr * _ref->get__gij(0, 0));
			double huv = ( (duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ( (duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2*dvv - tr * _ref->get__gij(1, 1));

			double Huu =  globalratio*(-hvv*sc)+get__hij(1, 1) * sc;
			double Huv = globalratio * (huv*sc)-get__hij(0, 1) * sc;
			double Hvv = globalratio * (-huu*sc)+get__hij(0, 0) * sc;


			double Suu =  get__Sij(0, 0);
			double Suv =  get__Sij(0, 1);
			double Svv =  get__Sij(1, 1);

			val = Hvv * Svv + 2 * Huv * Suv + Huu * Suu;

			/*if (accurate_area)
				val -= load * dv / _ref->_refDv;
			else
				val -= load;*/
			return val;

		}

		double __bodyF2_rhs(double load, bool accurate_area)
		{
			double val = 0;

			if (accurate_area)
				val = load * dv / _ref->_refDv;
			else
				val = load;
			return val;

		}
		void __bodyF2_phi(double* ptr,double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;


			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[3][s] * sc;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = -_ref->__dh[1][s] * sc;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[0][s] * sc;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;
				double _h21 = _h12;
				val = (_h11 * Suu + 2 * _h12 * Suv + _h22 * Svv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF2_z( double* ptr, double load, bool accurate_area, double globalratio)
		{

			double* ptr1 = ptr;
			double val = 0;
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);



			double tr = duu * _ref->get__Gij(0, 0) + duv * _ref->get__Gij(0, 1) + dvu * _ref->get__Gij(1, 0) + dvv * _ref->get__Gij(1, 1);

			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

			//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
			//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
			//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);
			double Huu = globalratio * (-hvv*sc)  +get__hij(1, 1) * sc;
			double Huv = globalratio * (huv*sc)  -get__hij(0, 1) * sc;
			double Hvv = globalratio * (-huu*sc)  +get__hij(0, 0) * sc;
			//huu = lambda * (Huu)  +get__hij(1, 1) * sc;
			//huv = lambda * (Huv)  -get__hij(0, 1) * sc;
			//hvv = lambda * (Hvv)  +get__hij(0, 0) * sc;

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];
				double _S12 = _ref->__dh[1][s];
				double _S22 = _ref->__dh[3][s];
				double _S21 = _S12;
				val = Hvv* _S22 + 2 * Huv * _S12 + Huu* _S11;
				if (accurate_area) {
					double _g11 = 0, _g12 = 0, _g22 = 0;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					}

					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					val += -load * ddv / _ref->_refDv;
				}
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF2_x(double* ptr, double load, bool accurate_area, double globalratio)
		{
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double g11 = _ref->get__gij(0, 0);
			double g12 = _ref->get__gij(0, 1);
			double g22 = _ref->get__gij(1, 1);
			double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);
			double G21 = G12;


			double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;

			
			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

			//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
			//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
			//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


			double Fuu = globalratio * (-hvv*sc)  + get__hij(1, 1)*sc;
			double Fuv = globalratio * (huv*sc)  - get__hij(0, 1)*sc;
			double Fvv = globalratio * (-huu*sc) + get__hij(0, 0)*sc;
			double Fvu = Fuv;


			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);
			double Svu = get__Sij(1, 0);
			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}

			
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;


				val = _Suu * Fuu + 2 * _Suv * Fuv + _Svv * Fvv;


				double _fuu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _fuv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _fvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _fvu = _fuv;
				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 0);
				_g21 = _g12;

				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2*_ref->get__gij(0, 1) * _g12;
				double _sc = -1 / det / det * _det;
				double _Fuu = globalratio * (-hvv * _sc) + get__hij(1, 1) * _sc;
				double _Fuv = globalratio * (huv * _sc) - get__hij(0, 1) * _sc;
				double _Fvv = globalratio * (-huu * _sc) + get__hij(0, 0) * _sc;

				_Fuu += _fvv * sc;
					_Fuv += -_fuv * sc;
					_Fvv += _fuu * sc;

				
			

				double _duu = xu * _ref->d1[0][s];// +yu * _ref->d1[0][s];
				double _duv = xu * _ref->d1[1][s];// + yu * _ref->d1[1][s];
				double _dvu = xv * _ref->d1[0][s];// + yv * _ref->d1[0][s];
				double _dvv = xv * _ref->d1[1][s];// + yv * _ref->d1[1][s];

				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;

				double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;
				double _tr = duu * _G11 + duv * _G12 + dvu * _G12 + dvv * _G22 + _duu * _ref->get__Gij(0, 0) + _dvu * _ref->get__Gij(0, 1) + _duv * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);

				double _huu = (2*_duu-tr * _g11 - _tr * g11);
				double _huv = (_duv+_dvu-tr * _g12 - _tr * g12);
				double _hvu = (_duv + _dvu -tr * _g12 - _tr * g12);
				double _hvv = (2*_dvv-tr * _g22 - _tr * g22);
				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);

				_Fuu += -globalratio * _hvv * sc;
				_Fuv += globalratio * _huv * sc;
				_Fvv += -globalratio * _huu * sc;
				val += Suu * _Fuu + 2 * Suv * _Fuv + Svv * _Fvv;
				if (accurate_area) {
				
					double _dv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					double __refDv = 0.5 * (_g11 * _ref->get__Gij(0, 0) + _g22 *_ref->get__Gij(1, 1) + 2 * _g12 * _ref->get__Gij(0, 1)) * _ref->_refDv;
					val += -load * _dv / _ref->_refDv;
					val += +load * dv / _ref->_refDv / _ref->_refDv * __refDv;
				}
				*ptr1 = val;
				ptr1++;
			}

		}

		void __bodyF2_y( double* ptr,double load,bool accurate_area, double globalratio)
		{
			double xu = 0, xv = 0, yu = 0, yv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				xu += _ref->d1[0][s] * _ref->buf_xi[s];
				xv += _ref->d1[1][s] * _ref->buf_xi[s];
				yu += _ref->d1[0][s] * _ref->buf_eta[s];
				yv += _ref->d1[1][s] * _ref->buf_eta[s];
			}
			double duu = xu * _ref->get__gi(0, 0) + yu * _ref->get__gi(0, 1);
			double duv = xu * _ref->get__gi(1, 0) + yu * _ref->get__gi(1, 1);
			double dvu = xv * _ref->get__gi(0, 0) + yv * _ref->get__gi(0, 1);
			double dvv = xv * _ref->get__gi(1, 0) + yv * _ref->get__gi(1, 1);

			double g11 = _ref->get__gij(0, 0);
			double g12 = _ref->get__gij(0, 1);
			double g22 = _ref->get__gij(1, 1);
			double G11 = _ref->get__Gij(0, 0), G12 = _ref->get__Gij(0, 1), G22 = _ref->get__Gij(1, 1);
			double G21 = G12;


			double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;


			double huu = (2 * duu - tr * _ref->get__gij(0, 0));
			double huv = ((duv + dvu) - tr * _ref->get__gij(0, 1));
			double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
			double hvv = (2 * dvv - tr * _ref->get__gij(1, 1));

			//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
			//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
			//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);


			double Fuu = globalratio * (-hvv*sc)+get__hij(1, 1) * sc;
			double Fuv = globalratio * (huv*sc)-get__hij(0, 1) * sc;
			double Fvv = globalratio * (-huu*sc)+get__hij(0, 0) * sc;
			double Fvu = Fuv;


			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);
			double Svu = get__Sij(1, 0);
			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}

			double scale = 1 / _ref->_refDv;
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double Xu = 0, Xv = 0, Yu = 0, Yv = 0;
				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _fuu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _fuv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _fvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _fvu = _fuv;
				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;
				_g11 = 2 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				_g22 = 2 * _ref->d1[1][s] * _ref->get__gi(1, 1);
				_g21 = _g12;


				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2 * _ref->get__gij(0, 1) * _g12;
				double _sc = -1 / det / det * _det;



				double _Fuu = globalratio * (-hvv * _sc) + get__hij(1, 1) * _sc;
				double _Fuv = globalratio * (huv * _sc) - get__hij(0, 1) * _sc;
				double _Fvv = globalratio * (-huu * _sc) + get__hij(0, 0) * _sc;

				_Fuu += _fvv * sc;
				_Fuv += -_fuv * sc;
				_Fvv += _fuu * sc;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;


				val = _Suu * Fuu + 2 * _Suv * Fuv + _Svv * Fvv;

			

		


				double _duu = yu * _ref->d1[0][s];
				double _duv = yu * _ref->d1[1][s];
				double _dvu = yv * _ref->d1[0][s];
				double _dvv = yv * _ref->d1[1][s];

				double _G11 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 0));
				double _G12 = -(_ref->get__Gij(0, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G22 = -(_ref->get__Gij(1, 0) * _g11 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * _g21 * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * _g12 * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * _g22 * _ref->get__Gij(1, 1));
				double _G21 = _G12;

				double tr = duu * G11 + duv * G12 + dvu * G12 + dvv * G22;
				double _tr = duu * _G11 + duv * _G12 + dvu * _G12 + dvv * _G22 + _duu * _ref->get__Gij(0, 0) + _dvu * _ref->get__Gij(0, 1) + _duv * _ref->get__Gij(0, 1) + _dvv * _ref->get__Gij(1, 1);

				double _huu = (2 * _duu - tr * _g11 - _tr * g11);
				double _huv = (_duv + _dvu - tr * _g12 - _tr * g12);
				
				double _hvv = (2 * _dvv - tr * _g22 - _tr * g22);
				double _hvu = _huv;
				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);

				_Fuu += -globalratio * _hvv*sc;
				_Fuv += globalratio * _huv*sc;
				_Fvv += -globalratio * _huu*sc;
				val += Suu * _Fuu + 2 * Suv * _Fuv + Svv * _Fvv;
				if (accurate_area) {
					
					double _dv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					double __refDv = 0.5 * (_g11 * _ref->get__Gij(0, 0) + _g22 * _ref->get__Gij(1, 1) + 2 * _g12 * _ref->get__Gij(0, 1)) * _ref->_refDv;
					val += -load * _dv / _ref->_refDv;
					val += +load * dv / _ref->_refDv / _ref->_refDv * __refDv;
				}

				*ptr1 = val;
				ptr1++;
			}
		}

		
		void __bodyF2_xi(double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;
		
			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;
				_xu = _ref->d1[0][s] ;
				_xv = _ref->d1[1][s] ;
				_yu = 0;
				_yv = 0;
				
				double _duu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _duv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _dvu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _dvv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);



				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _huu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _huv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				//double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double _hvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				double _hvu = _huv;
				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);
				double _fuu = -globalratio * (_hvv) * sc;
				double _fuv = globalratio * (_huv)*sc;
				double _fvv = -globalratio * (_huu) * sc;

				val =  (Suu * _fuu + 2 * Suv * _fuv + Svv * _fvv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF2_eta( double* ptr, double globalratio)
		{


			double val = 0;
			double* ptr1 = ptr;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);


			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _xu = 0, _xv = 0, _yu = 0, _yv = 0;
				_xu = 0;// _ref->d1[0][s];
				_xv = 0;//_ref->d1[1][s];
				_yu = _ref->d1[0][s];
				_yv = _ref->d1[1][s];

				double _duu = _xu * _ref->get__gi(0, 0) + _yu * _ref->get__gi(0, 1);
				double _duv = _xu * _ref->get__gi(1, 0) + _yu * _ref->get__gi(1, 1);
				double _dvu = _xv * _ref->get__gi(0, 0) + _yv * _ref->get__gi(0, 1);
				double _dvv = _xv * _ref->get__gi(1, 0) + _yv * _ref->get__gi(1, 1);



				double _tr = _duu * _ref->get__Gij(0, 0) + _duv * _ref->get__Gij(0, 1) + _dvu * _ref->get__Gij(1, 0) + _dvv * _ref->get__Gij(1, 1);

				double _huu = (2 * _duu - _tr * _ref->get__gij(0, 0));
				double _huv = ((_duv + _dvu) - _tr * _ref->get__gij(0, 1));
				//double hvu = ((duv + dvu) - tr * _ref->get__gij(1, 0));
				double _hvv = (2 * _dvv - _tr * _ref->get__gij(1, 1));
				double _hvu = _huv;
				//double Huu = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 0) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 0) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 0);
				//double Huv = _ref->get__Gij(0, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(0, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(0, 1) * hvv * _ref->get__Gij(1, 1);
				//double Hvv = _ref->get__Gij(1, 0) * huu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 0) * huv * _ref->get__Gij(1, 1) + _ref->get__Gij(1, 1) * hvu * _ref->get__Gij(0, 1) + _ref->get__Gij(1, 1) * hvv * _ref->get__Gij(1, 1);
				double _fuu = -globalratio * (_hvv)*sc;
				double _fuv = globalratio * (_huv)*sc;
				double _fvv = -globalratio * (_huu) * sc;

				val = (Suu * _fuu + 2 * Suv * _fuv + Svv * _fvv);
				*ptr1 = val;
				ptr1++;
			}

		}

		double __bodyU(double load, bool accurate_area)
		{
			double val = 0;

			double h11 = 0, h12 = 0, h22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				h11 += _ref->d2[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_phi[s];
				h12 += _ref->d2[1][s] * _ref->buf_phi[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_phi[s];
				h22 += _ref->d2[3][s] * _ref->buf_phi[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_phi[s];
			}

			double Huu = +h22 * _ref->osc;
			double Huv = -h12 * _ref->osc;
			double Hvv = +h11 * _ref->osc;

			double Suu = 0, Suv = 0, Svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Suu += _ref->d2[0][s] * _ref->buf_u[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_u[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_u[s];
				Suv += _ref->d2[1][s] * _ref->buf_u[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_u[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_u[s];
				Svv += _ref->d2[3][s] * _ref->buf_u[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_u[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_u[s];
			}

			val = Hvv * Svv + 2 * Huv * Suv + Huu * Suu;

			return val;

		}
		void __bodyU_phi(double* ptr)
		{

			double val = 0;
			double* ptr1 = ptr;

			double Suu = 0, Suv = 0, Svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Suu += _ref->d2[0][s] * _ref->buf_u[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_u[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_u[s];
				Suv += _ref->d2[1][s] * _ref->buf_u[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_u[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_u[s];
				Svv += _ref->d2[3][s] * _ref->buf_u[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_u[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_u[s];
			}

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s];
				double _h12 = _ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s];
				double _h22 = _ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s];
				double _H11 = _h22* _ref->osc, _H22 = _h11 * _ref->osc, _H12 = -_h12 * _ref->osc;

				val = (_H11 * Suu + 2 * _H12 * Suv + _H22 * Svv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyU_u(double* ptr, double load, bool accurate_area)
		{

			double* ptr1 = ptr;
			double val = 0;

			double h11 = 0, h12 = 0, h22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				h11 += _ref->d2[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_phi[s];
				h12 += _ref->d2[1][s] * _ref->buf_phi[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_phi[s];
				h22 += _ref->d2[3][s] * _ref->buf_phi[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_phi[s];
			}

			double Huu = +h22 * _ref->osc;
			double Huv = -h12 * _ref->osc;
			double Hvv = +h11 * _ref->osc;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = _ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s];
				double _S12 = _ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s];
				double _S22 = _ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s];
				val = Hvv * _S22 + 2 * Huv * _S12 + Huu * _S11;
				*ptr1 = val;
				ptr1++;
			}
		}
		double __bodyV(double load, bool accurate_area)
		{
			double val = 0;

			double h11 = 0, h12 = 0, h22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				h11 += _ref->d2[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_phi[s];
				h12 += _ref->d2[1][s] * _ref->buf_phi[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_phi[s];
				h22 += _ref->d2[3][s] * _ref->buf_phi[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_phi[s];
			}

			double Huu = +h22 * _ref->osc;
			double Huv = -h12 * _ref->osc;
			double Hvv = +h11 * _ref->osc;

			double Suu = 0, Suv = 0, Svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Suu += _ref->d2[0][s] * _ref->buf_v[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_v[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_v[s];
				Suv += _ref->d2[1][s] * _ref->buf_v[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_v[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_v[s];
				Svv += _ref->d2[3][s] * _ref->buf_v[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_v[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_v[s];
			}

			val = Hvv * Svv + 2 * Huv * Suv + Huu * Suu;

			return val;

		}
		void __bodyV_phi(double* ptr)
		{

			double val = 0;
			double* ptr1 = ptr;

			double Suu = 0, Suv = 0, Svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Suu += _ref->d2[0][s] * _ref->buf_v[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_v[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_v[s];
				Suv += _ref->d2[1][s] * _ref->buf_v[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_v[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_v[s];
				Svv += _ref->d2[3][s] * _ref->buf_v[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_v[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_v[s];
			}

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s];
				double _h12 = _ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s];
				double _h22 = _ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s];
				double _H11 = _h22 * _ref->osc, _H22 = _h11 * _ref->osc, _H12 = -_h12 * _ref->osc;

				val = (_H11 * Suu + 2 * _H12 * Suv + _H22 * Svv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyV_v(double* ptr, double load, bool accurate_area)
		{

			double* ptr1 = ptr;
			double val = 0;

			double h11 = 0, h12 = 0, h22 = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				h11 += _ref->d2[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_phi[s];
				h12 += _ref->d2[1][s] * _ref->buf_phi[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_phi[s];
				h22 += _ref->d2[3][s] * _ref->buf_phi[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_phi[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_phi[s];
			}

			double Huu = +h22 * _ref->osc;
			double Huv = -h12 * _ref->osc;
			double Hvv = +h11 * _ref->osc;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _S11 = _ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s];
				double _S12 = _ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s];
				double _S22 = _ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s];
				val = Hvv * _S22 + 2 * Huv * _S12 + Huu * _S11;
				*ptr1 = val;
				ptr1++;
			}
		}
		double __bodyF3(double load, bool accurate_area)
		{
			double val = 0;

		
			double h11 = get__hij(0, 0);
			double h12 = get__hij(0, 1);
			double h22 = get__hij(1, 1);

			double Huu = h22 *sc;
			double Huv =  - h12 * sc;
			double Hvv = h11 * sc;

			/*double Suu = 0, Suv = 0, Svv = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				Suu += _ref->d2[0][s] * _ref->buf_z[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_z[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_z[s];
				Suv += _ref->d2[1][s] * _ref->buf_z[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_z[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_z[s];
				Svv += _ref->d2[3][s] * _ref->buf_z[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_z[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_z[s];
			}*/

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);

			val = Hvv * Svv + 2 * Huv * Suv + Huu * Suu;

			if (accurate_area)
				val -= load * dv / _dv;
			else
				val -= load;
			return val;

		}
		void __bodyF3_phi(double* ptr)
		{


			double val = 0;
			double* ptr1 = ptr;

			

			/*for (int s = 0; s < _ref->_nNode; s++)
			{
				Suu += _ref->d2[0][s] * _ref->buf_z[s] - _ref->oGammaijk[0] * _ref->d1[0][s] * _ref->buf_z[s] - _ref->oGammaijk[1] * _ref->d1[1][s] * _ref->buf_z[s];
				Suv += _ref->d2[1][s] * _ref->buf_z[s] - _ref->oGammaijk[2] * _ref->d1[0][s] * _ref->buf_z[s] - _ref->oGammaijk[3] * _ref->d1[1][s] * _ref->buf_z[s];
				Svv += _ref->d2[3][s] * _ref->buf_z[s] - _ref->oGammaijk[6] * _ref->d1[0][s] * _ref->buf_z[s] - _ref->oGammaijk[7] * _ref->d1[1][s] * _ref->buf_z[s];
			}*/
			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);

			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;
				double _h21 = _h12;
				
				//double _h11 = _ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s];
				//double _h12 = _ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s];
				//double _h22 = _ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s];
				double _H11 = _h22* sc, _H22 = _h11 * sc, _H12 = -_h12 * sc;

				val = (_H11 * Suu + 2 * _H12 * Suv + _H22 * Svv);
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF3_z(double* ptr, double load, bool accurate_area)
		{

			double* ptr1 = ptr;
			double val = 0;
			

			double Huu = get__hij(1, 1) * sc;
			double Huv =  - get__hij(0, 1) * sc;
			double Hvv = get__hij(0, 0) * sc;
			
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _S11 = _ref->__dh[0][s];
				double _S12 = _ref->__dh[1][s];
				double _S22 = _ref->__dh[3][s];
				double _S21 = _S12;

				//double _S11 = _ref->d2[0][s] - _ref->oGammaijk[0] * _ref->d1[0][s] - _ref->oGammaijk[1] * _ref->d1[1][s];
				//double _S12 = _ref->d2[1][s] - _ref->oGammaijk[2] * _ref->d1[0][s] - _ref->oGammaijk[3] * _ref->d1[1][s];
				//double _S22 = _ref->d2[3][s] - _ref->oGammaijk[6] * _ref->d1[0][s] - _ref->oGammaijk[7] * _ref->d1[1][s];


				val = Hvv * _S22 + 2 * Huv * _S12 + Huu * _S11;
				if (accurate_area) {
					/*double _g11 = 0, _g12 = 0, _g22 = 0;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2.0 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g22 += 2.0 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					}*/
					double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

					_g11 = 2.0 * _ref->d1[0][s] * get_gi2(0, 2);
					_g12 = _ref->d1[0][s] * get_gi2(1, 2) + _ref->d1[1][s] * get_gi2(0, 2);
					_g22 = 2.0 * _ref->d1[1][s] * get_gi2(1, 2);
					_g21 = _g12;
					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					
					//double _ddv = 0.5 * (_g11 * this->get_Gij(0, 0) + _g22 * this->get_Gij(1, 1) + 2 * _g12 * this->get_Gij(0, 1)) * this->_dv;
					val += -load * ddv /  _dv;// / _dv;
					//val -= -load * dv / _dv / _dv * _ddv;
				}
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF3_x(double* ptr, double load, bool accurate_area)
		{



			double Huu = get__hij(1, 1) * _ref->osc;
			double Huv = -get__hij(0, 1) * _ref->osc;
			double Hvv = get__hij(0, 0) * _ref->osc;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);
			double Svu = get__Sij(1, 0);
			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
			    f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}
			

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 0);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 0);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 0);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 0);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 0);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;

				double _Huu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _Huv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _Hvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _Hvu = _Huv;



				val = _Suu * Huu + 2.0 * _Suv * Huv + _Svv * Hvv;


				
				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2.0 * _ref->d1[0][s] * _ref->get__gi(0, 0);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 0) + _ref->d1[1][s] * _ref->get__gi(0, 0);
				_g22 = 2.0 * _ref->d1[1][s] * _ref->get__gi(1, 0);
				_g21 = _g12;

				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2.0 * _ref->get__gij(0, 1) * _g12;
				double _sc =  -1.0 / det / det * _det;

				double _Fuu = get__hij(1, 1) * _sc;
				double _Fuv = - get__hij(0, 1) * _sc;
				double _Fvv =  get__hij(0, 0) * _sc;

				_Fuu += _Hvv * sc;
				_Fuv += -_Huv * sc;
				_Fvv += _Huu * sc;

				val += Suu * _Fuu + 2.0 * Suv * _Fuv + Svv * _Fvv;

				
				
				
				if (accurate_area) {

					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					
					double _ddv = 0.5 * (_g11 * this->get_Gij(0, 0) + _g22 * this->get_Gij(1, 1) + 2 * _g12 * this->get_Gij(0, 1)) * this->_dv;
					val += -load * ddv / _dv;
					val -= -load * dv / _dv / _dv * _ddv;
				
					//val += +load * dv / _ref->_refDv / _ref->_refDv * __refDv;
				}
				
				*ptr1 = val;
				ptr1++;
			}

		}

		void __bodyF3_y(double* ptr, double load, bool accurate_area)
		{
		

			double Huu = get__hij(1, 1) * _ref->osc;
			double Huv = -get__hij(0, 1) * _ref->osc;
			double Hvv = get__hij(0, 0) * _ref->osc;

			double Suu = get__Sij(0, 0);
			double Suv = get__Sij(0, 1);
			double Svv = get__Sij(1, 1);
			double Svu = get__Sij(1, 0);
			double f1 = 0, f2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				f1 += _ref->d1[0][s] * _ref->buf_phi[s];
				f2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}


			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{


				double _Gamma111 = _ref->__dh[0][s] * _ref->get__Gi(0, 1);
				double _Gamma112 = _ref->__dh[0][s] * _ref->get__Gi(1, 1);
				double _Gamma121 = _ref->__dh[1][s] * _ref->get__Gi(0, 1);
				double _Gamma122 = _ref->__dh[1][s] * _ref->get__Gi(1, 1);
				double _Gamma221 = _ref->__dh[3][s] * _ref->get__Gi(0, 1);
				double _Gamma222 = _ref->__dh[3][s] * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _Suu = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _Suv = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _Svv = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _Svu = _Suv;

				double _Huu = (-_Gamma111 * f1 - _Gamma112 * f2);
				double _Huv = (-_Gamma121 * f1 - _Gamma122 * f2);
				double _Hvv = (-_Gamma221 * f1 - _Gamma222 * f2);
				double _Hvu = _Huv;



				val = _Suu * Huu + 2.0 * _Suv * Huv + _Svv * Hvv;



				double _g11 = 0, _g12 = 0, _g21 = 0, _g22 = 0;

				_g11 = 2.0 * _ref->d1[0][s] * _ref->get__gi(0, 1);
				_g12 = _ref->d1[0][s] * _ref->get__gi(1, 1) + _ref->d1[1][s] * _ref->get__gi(0, 1);
				_g22 = 2.0 * _ref->d1[1][s] * _ref->get__gi(1, 1);
				_g21 = _g12;

				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _ref->get__gij(0, 0) * _g22 + _ref->get__gij(1, 1) * _g11 - 2.0 * _ref->get__gij(0, 1) * _g12;
				double _sc =  -1.0 / det / det * _det;

				double _Fuu = get__hij(1, 1) * _sc;
				double _Fuv = -get__hij(0, 1) * _sc;
				double _Fvv = get__hij(0, 0) * _sc;

				_Fuu += _Hvv *sc;
				_Fuv += -_Huv *sc;
				_Fvv += _Huu * sc;

				val += Suu * _Fuu + 2.0 * Suv * _Fuv + Svv * _Fvv;




				if (accurate_area) {

					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;

					double _ddv = 0.5 * (_g11 * this->get_Gij(0, 0) + _g22 * this->get_Gij(1, 1) + 2 * _g12 * this->get_Gij(0, 1)) * this->_dv;
					val += -load * ddv / _dv;
					val -= -load * dv / _dv / _dv * _ddv;

					//val += +load * dv / _ref->_refDv / _ref->_refDv * __refDv;
				}

				*ptr1 = val;
				ptr1++;
			}
		}


		


		double __bodyF( double load,bool accurate_area)
		{
			double val = 0;

		

			val = (get__hij(0, 0) * get__Sij(1, 1) - 2 * get__hij(0, 1)  * get__Sij(0, 1) + get__hij(1, 1) * get__Sij(0, 0)) * sc;
			
			if (accurate_area)
				val -= load * dv / _ref->_refDv;
			else
				val -= load;
			return val;

		}
		
		void __bodyF_x( double* ptr, double load, bool accurate_area)
		{

			double* ptr1 = ptr;
			double val = 0;
			double h1 = 0, h2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				h1 += _ref->d1[0][s] * _ref->buf_phi[s];
				h2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _Gamma111 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(0, 0);
				double _Gamma112 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(1, 0);
				double _Gamma121 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(0, 0);
				double _Gamma122 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(1, 0);
				double _Gamma221 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(0, 0);
				double _Gamma222 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(1, 0);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _h11 = (-_Gamma111 * h1 - _Gamma112 * h2) ;
				double _h12 = (-_Gamma121 * h1 - _Gamma122 * h2) ;
				double _h22 = (-_Gamma221 * h1 - _Gamma222 * h2) ;
				double _h21 = _h12;

				double _S11 = (-_Gamma111 * S1 - _Gamma112 * S2) ;
				double _S12 = (-_Gamma121 * S1 - _Gamma122 * S2) ;
				double _S22 = (-_Gamma221 * S1 - _Gamma222 * S2) ;
				double _S21 = _S12;


				//double _g11 = 2 * (_ref->d1[0][s] * _ref->get__gi(0, 0));
				//double _g12 = (_ref->d1[0][s] * _ref->get__gi(1, 0))+ (_ref->d1[1][s] * _ref->get__gi(0, 0));
				//double _g22 = 2 * (_ref->d1[1][s] * _ref->get__gi(1, 0));
				double _g11 = __dsigma_11[0][s];
				double _g12 = __dsigma_12[0][s];
				double _g22 = __dsigma_22[0][s];
				double _g21 = _g12;

				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _g11 * _ref->get__gij(1, 1)+ _g22 * _ref->get__gij(0, 0) - 2*_g12 * _ref->get__gij(0, 1);
				double _sc = -1 / det / det * _det;

				val = (get__hij(0, 0) * get__Sij(1, 1) - 2 * get__hij(0, 1) * get__Sij(0, 1) + get__hij(1, 1) * get__Sij(0, 0)) * _sc;
				val += (_h11 * get__Sij(1, 1) - 2 * _h12 * get__Sij(0, 1) + _h22 * get__Sij(0, 0)) * sc;
				val += (get__hij(0, 0) * _S22 - 2 * get__hij(0, 1) * _S12 + get__hij(1, 1) * _S11)*sc;
				if (accurate_area)
				{
					double _dv=
					val -= load * dv / _ref->_refDv * (0.5 * (_g11 * get_Gij2(0, 0) + 2 * _g12 * get_Gij2(0, 1) + _g22 * get_Gij2(1, 1)));
					val-=-load*dv/_ref->_refDv/_ref->_refDv* (0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0, 1) + _g22 * _ref->get__Gij(1, 1)));

				

				}else{
					val -= load*(0.5*(_g11*_ref->get__Gij(0,0)+ 2 * _g12 * _ref->get__Gij(0, 1)+ _g22 * _ref->get__Gij(1, 1)));
				}
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF_y( double* ptr, double load, bool accurate_area)
		{

			double* ptr1 = ptr;
			double val = 0;
			double h1 = 0, h2 = 0, S1 = 0, S2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				h1 += _ref->d1[0][s] * _ref->buf_phi[s];
				h2 += _ref->d1[1][s] * _ref->buf_phi[s];
				S1 += _ref->d1[0][s] * _ref->buf_z[s];
				S2 += _ref->d1[1][s] * _ref->buf_z[s];
			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _Gamma111 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(0, 1);
				double _Gamma112 = (_ref->d2[0][s] - (_ref->_Gammaijk[0] * _ref->d1[0][s] + _ref->_Gammaijk[1] * _ref->d1[1][s])) * _ref->get__Gi(1, 1);
				double _Gamma121 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(0, 1);
				double _Gamma122 = (_ref->d2[1][s] - (_ref->_Gammaijk[2] * _ref->d1[0][s] + _ref->_Gammaijk[3] * _ref->d1[1][s])) * _ref->get__Gi(1, 1);
				double _Gamma221 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(0, 1);
				double _Gamma222 = (_ref->d2[3][s] - (_ref->_Gammaijk[6] * _ref->d1[0][s] + _ref->_Gammaijk[7] * _ref->d1[1][s])) * _ref->get__Gi(1, 1);
				double _Gamma211 = _Gamma121;
				double _Gamma212 = _Gamma122;

				double _h11 = (-_Gamma111 * h1 - _Gamma112 * h2);
				double _h12 = (-_Gamma121 * h1 - _Gamma122 * h2);
				double _h22 = (-_Gamma221 * h1 - _Gamma222 * h2);
				double _h21 = _h12;

				double _S11 = (-_Gamma111 * S1 - _Gamma112 * S2);
				double _S12 = (-_Gamma121 * S1 - _Gamma122 * S2);
				double _S22 = (-_Gamma221 * S1 - _Gamma222 * S2);
				double _S21 = _S12;
				

				double _g11 = __dsigma_11[1][s];
				double _g12 = __dsigma_12[1][s];
				double _g22 = __dsigma_22[1][s];
				double _g21 = _g12; 
				
				double det = _ref->get__gij(0, 0) * _ref->get__gij(1, 1) - _ref->get__gij(0, 1) * _ref->get__gij(0, 1);
				double _det = _g11 * _ref->get__gij(1, 1) + _g22 * _ref->get__gij(0, 0) - 2 * _g12 * _ref->get__gij(0, 1);
				double _sc = -1 / det / det * _det;

				val = (get__hij(0, 0) * get__Sij(1, 1) - 2 * get__hij(0, 1) * get__Sij(0, 1) + get__hij(1, 1) * get__Sij(0, 0)) * _sc;
				val += (_h11 * get__Sij(1, 1) - 2 * _h12 * get__Sij(0, 1) + _h22 * get__Sij(0, 0)) * sc;
				val += (get__hij(0, 0) * _S22 - 2 * get__hij(0, 1) * _S12 + get__hij(1, 1) * _S11) * sc;
				if (accurate_area)
				{
					val -= load * dv / _ref->_refDv * (0.5 * (_g11 * get_Gij2(0, 0) + 2 * _g12 * get_Gij2(0, 1) + _g22 * get_Gij2(1, 1)));
					val -= -load * dv / _ref->_refDv / _ref->_refDv * (0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0, 1) + _g22 * _ref->get__Gij(1, 1)));
				}
				else {
					//val -= load * (0.5 * (_g11 * _ref->get__Gij(0, 0) + 2 * _g12 * _ref->get__Gij(0, 1) + _g22 * _ref->get__Gij(1, 1)));
				}
				*ptr1 = val;
				ptr1++;
			}

		}

		double mean()
		{
			double val = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s11 += _ref->__dh[0][s] * _ref->buf_nu[s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 += _ref->__dh[1][s] * _ref->buf_nu[s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s22 += _ref->__dh[3][s] * _ref->buf_nu[s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);


			}
			val = s11 * _ref->get__Gij(0, 0) + 2*s12 * _ref->get__Gij(0, 1) + s22 * _ref->get__Gij(1, 1);
			return val;
		}
		void mean_nu( double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
		
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _s11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _s12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _s22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _s21 = _s12;
				*ptr1 = _s11 * _ref->get__Gij(0, 0) + 2 * _s12 * _ref->get__Gij(0, 1) + _s22 * _ref->get__Gij(1, 1);

				ptr1++;
			}
		}
		double gauss()
		{
			double val = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s11 += _ref->__dh[0][s] * _ref->buf_nu[s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 += _ref->__dh[1][s] * _ref->buf_nu[s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s22 += _ref->__dh[3][s] * _ref->buf_nu[s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				

			}
			val = s11 * s22 - s12 * s12;
			return val;
		}
		void gauss_nu( double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			double s11 = 0, s12 = 0, s22 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				s11 += _ref->__dh[0][s] * _ref->buf_nu[s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				s12 += _ref->__dh[1][s] * _ref->buf_nu[s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				s22 += _ref->__dh[3][s] * _ref->buf_nu[s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);


			}
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _s11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _s12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _s22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _s21 = _s12;
				*ptr1 = _s11 * s22+s11*_s22 - 2*_s12 * s12;
				ptr1++;
			}
		}
		void __bodyF_z( double* ptr,double load,bool accurate_area)
		{

			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				
				double _S11 = _ref->__dh[0][s];// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				double _S12 = _ref->__dh[1][s];// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				double _S22 = _ref->__dh[3][s];// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);
				double _S21 = _S12;
				val = (get__hij(0, 0) * _S22 - 2 * get__hij(0, 1) * _S12 + get__hij(1, 1) * _S11) * sc;
				if (accurate_area) {
					double _g11 = 0, _g12 = 0, _g22 = 0;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * _ref->buf_z[t];

						_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * _ref->buf_z[t];
					}

					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					val += -load * ddv / _ref->_refDv;
				}
				*ptr1 = val;
				ptr1++;
			}

		}
		void __bodyF_phi( double* ptr)
		{

			
			double val = 0;
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{

				double _h11 = _ref->__dh[0][s] ;// (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s])* _ref->sc11;
				double _h12 = _ref->__dh[1][s] ;// (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s])* _ref->sc12;
				double _h22 = _ref->__dh[3][s] ;// (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->sc22;
				double _h21 = _h12;
				val = (_h11 * get__Sij(1, 1) - 2 * _h12 * get__Sij(0, 1) + _h22 * get__Sij(0, 0)) * sc;
				*ptr1 = val;
				ptr1++;
			}

		}
		/*void __bodyF_phi(double* ptr, bool accurate)
		{
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			double _s11 = 0;
			double _s12 = 0;
			double _s22 = 0;
			double _S11 = 0, _S12 = 0, _S22 = 0;
			
			//double _s21 = _s12;


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				//s22 = (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]);
				//s12 = -(_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]);
				//s11 = (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]);



				//double val = (_S11 * s11 + _S22 * s22 + 2 * _S12 * s12) * sc;
				double val = __grad_z[s];
				
					
				*ptr1 = val;
				
				ptr1++;
			}
		}
		void __bodyF_z(double *ptr,double load,bool accurate)
		{		
			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double val=__grad_phi[s] ;
				if (accurate)
				{
					double _g11 = 0, _g12 = 0, _g22 = 0;
					double* ptr1 = _ref->buf_z;
					for (int t = 0; t < _ref->_nNode; t++)
					{
						_g11 += 2 * (_ref->d1[0][s]) * (_ref->d1[0][t]) * *ptr1;

						_g12 += (_ref->d1[0][s]) * (_ref->d1[1][t]) * *ptr1;
						_g12 += (_ref->d1[1][s]) * (_ref->d1[0][t]) * *ptr1;

						_g22 += 2 * (_ref->d1[1][s]) * (_ref->d1[1][t]) * *ptr1;
						ptr1++;
					}
					double _g21 = _g12;
					//__mem->sc * __mem->bodyF-load* __mem->dv/ __mem->_ref->_refDv;
					
					double ddv = 0.5 * (_g11 * this->get_Gij2(0, 0) + _g22 * this->get_Gij2(1, 1) + 2 * _g12 * this->get_Gij2(0, 1)) * this->dv;
					val += -load*ddv / _ref->_refDv;
					*ptr1 = val; 
				}
				else {
					*ptr1 = val;
				}
				ptr1++;
			}
		}*/
		
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

		double BCEQ(double t1/*up*/, double t2/*up*/,double n1,double n2,double stt, double a,double b,bool accurate)
		{
			double val = 0;

			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;

			double d1 = 0, d2 = 0, D1 = 0, D2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				d1 += _ref->d1[0][s] * _ref->buf_phi[s];
				d2 += _ref->d1[1][s] * _ref->buf_phi[s];
				D1 += _ref->d1[0][s] * _ref->buf_z[s];
				D2 += _ref->d1[1][s] * _ref->buf_z[s];
				/*s11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_phi[s];
				s22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_phi[s];
				S11 += (_ref->d2[0][s] - _ref->_Gammaijk[0] * _ref->d1[0][s] - _ref->_Gammaijk[1] * _ref->d1[1][s]) * _ref->buf_z[s];
				S12 += (_ref->d2[1][s] - _ref->_Gammaijk[2] * _ref->d1[0][s] - _ref->_Gammaijk[3] * _ref->d1[1][s]) * _ref->buf_z[s];
				S22 += (_ref->d2[3][s] - _ref->_Gammaijk[6] * _ref->d1[0][s] - _ref->_Gammaijk[7] * _ref->d1[1][s]) * _ref->buf_z[s];*/
			}
			double S11 = get__Sij(0,0);
			double S12 = get__Sij(0, 1);
			double S22 = get__Sij(1, 1);
			double s11 = get__hij(0, 0);
			double s12 = get__hij(0, 1);
			double s22 = get__hij(1, 1);

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

				//double  tension=get_gammaijk(0, 0, 0)* t1* t1* n1 + get_gammaijk(0, 0, 1) * t1 * t1 * n2 + 2 * (get_gammaijk(0, 1, 0) * t1 * t2 * n1 + get_gammaijk(0, 1, 1) * t1 * t2 * n2 ) +
				//		get_gammaijk(1, 1, 0) * t2 * t2 * n1 + get_gammaijk(1, 1, 1) * t2 * t2 * n2;
				//val -= tension * (d1 * n1 + d2 * n2 - dc);
				val += (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2 - dc);//Gammaijk
			}
			else {
				val += (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * stt;

			}
			return val;
		}
		void BCEQ_z(double* ptr, double t1, double t2,double n1,double n2,double stt,double a,double b ,bool accurate)
		{


			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
		
			double gamma = 0;
			double eta = 0;
			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;
			double d1 = 0, d2 = 0, D1 = 0, D2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				D1 += _ref->d1[0][s] * _ref->buf_z[s];
				D2 += _ref->d1[1][s] * _ref->buf_z[s];
				S11 += _ref->__dh[0][s] * _ref->buf_z[s];
				S12 += _ref->__dh[1][s] * _ref->buf_z[s];
				S22 += _ref->__dh[3][s] * _ref->buf_z[s];
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
				s11 = _ref->__dh[0][s];
				s12 = _ref->__dh[1][s];
				s22 = _ref->__dh[3][s];
				double S21 = S12;
				double s21 = s12;
				//s11 = 0;
				{
					val += (s11 * t1 * t1 + 2 * s12 * t1 * t2 + s22 * t2 * t2) * (D1 * n1 + D2 * n2);
				}
				if (stt < -100)
				{
					//double  tension = get_gammaijk(0, 0, 0) * t1 * t1 * n1 + get_gammaijk(0, 0, 1) * t1 * t1 * n2 + 2 * (get_gammaijk(0, 1, 0) * t1 * t2 * n1 + get_gammaijk(0, 1, 1) * t1 * t2 * n2) +
					//	get_gammaijk(1, 1, 0) * t2 * t2 * n1 + get_gammaijk(1, 1, 1) * t2 * t2 * n2;
					//val -= tension * (d1 * n1 + d2 * n2 - dc);

					val += (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2);

				}
				else {
					//val -= (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (stt);
				}

				

				*ptr1 = val;
				ptr1 ++;
			}


		}
		void BCEQ_phi(double* ptr, double t1, double t2,double n1,double n2, double stt, double a, double b, bool accurate)
		{
			double S11 = 0;
			double S12 = 0;
			double S22 = 0;
			double s11 = 0;
			double s12 = 0;
			double s22 = 0;
			
			double gamma = 0;
			//double eta = 0;
			//double det = 0;
			double length = sqrt(t1 * t1 * _ref->get__gij(0, 0) + 2 * t1 * t2 * _ref->get__gij(0, 1) + t2 * t2 * _ref->get__gij(1, 1));
			t1 /= length; t2 /= length;
			length = sqrt(n1 * n1 * _ref->get__gij(0, 0) + 2 * n1 * n2 * _ref->get__gij(0, 1) + n2 * n2 * _ref->get__gij(1, 1));
			n1 /= length; n2 /= length;
			//det = eta * gamma;
			double d1 = 0, d2 = 0, D1 = 0, D2 = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				d1 += _ref->d1[0][s] * _ref->buf_phi[s];
				d2 += _ref->d1[1][s] * _ref->buf_phi[s];
				s11 += _ref->__dh[0][s] * _ref->buf_phi[s];
				s12 += _ref->__dh[1][s] * _ref->buf_phi[s];
				s22 += _ref->__dh[3][s] * _ref->buf_phi[s];
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
				S11 = _ref->__dh[0][s];
				S12 = _ref->__dh[1][s];
				S22 = _ref->__dh[3][s];
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
					val += (S11 * t1 * t1 + 2 * S12 * t1 * t2 + S22 * t2 * t2) * (d1 * n1 + d2 * n2-dc);


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
		void memory2(_memS_ref* __mem)
		{
			if (mode == "U")
			{

				__mem->oG11 = get_Gij(0, 0);
				__mem->oG12 = get_Gij(0, 1);
				__mem->oG22 = get_Gij(1, 1);
				__mem->og11 = get_gij(0, 0);
				__mem->og12 = get_gij(0, 1);
				__mem->og22 = get_gij(1, 1);
				memcpy(__mem->_ogi, gi, sizeof(double) * 6);
				memcpy(__mem->_oGi,Gi, sizeof(double) * 6);
			}
			else {
				__mem->oG11 = get_Gij2(0, 0);
				__mem->oG12 = get_Gij2(0, 1);
				__mem->oG22 = get_Gij2(1, 1);
				__mem->og11 = get_gij2(0, 0);
				__mem->og12 = get_gij2(0, 1);
				__mem->og22 = get_gij2(1, 1);
				memcpy(__mem->_ogi, gi2, sizeof(double) * 6);
				memcpy(__mem->_oGi, Gi2, sizeof(double) * 6);
			}
			memcpy(__mem->oGammaijk, Gammaijk, sizeof(double) * 8);
		    __mem->osc = sc;
			__mem->orefDv = _dv;
			__mem->oX = this->x;
			__mem->oY = this->y;
			__mem->xi_1 = xi_1;
			__mem->xi_2 = xi_2;
			__mem->eta_1 = eta_1;
			__mem->eta_2 = eta_2;

			

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
			__mem->_xi = xi;
			__mem->_eta = eta;
			__mem->_nu = nu;
			__mem->_u = u;
			__mem->_v = v;
			__mem->_w = w;
			__mem->xi_1 = xi_1;
			__mem->xi_2 = xi_2;
			__mem->eta_1 = eta_1;
			__mem->eta_2 = eta_2;

			std::memcpy(__mem->_Sij, _Sij, sizeof(double) * 4);
			
			if (mode == "SLOPE")
			{
				std::memcpy(__mem->_Gi, Gi2, sizeof(double) * 6);
				std::memcpy(__mem->_Gij, Gij2, sizeof(double) * 4);
				std::memcpy(__mem->_gij, gij2, sizeof(double) * 4);
				std::memcpy(__mem->_gi, gi2, sizeof(double) * 6);
				__mem->_gi[2] = 0;
				__mem->_gi[5] = 0;
				__mem->_Gi[2] = 0;
				__mem->_Gi[5] = 0;

			}
			else {
				std::memcpy(__mem->_gij, gij, sizeof(double) * 4);
				std::memcpy(__mem->_Gi, Gi, sizeof(double) * 6);
				std::memcpy(__mem->_Gij, Gij, sizeof(double) * 4);
				std::memcpy(__mem->_gi, gi, sizeof(double) * 6);
				__mem->_gi[2] = 0;
				__mem->_gi[5] = 0;
				__mem->_Gi[2] = 0;
				__mem->_Gi[5] = 0;
			}
			//std::memset(__mem->_gij, 0, sizeof(double) * 4);
			//std::memset(__mem->_Gij, 0, sizeof(double) * 4);
			std::memcpy(__mem->_bij, bij, sizeof(double) * 12);
			std::memcpy(__mem->_Gammaijk, Gammaijk, sizeof(double) * 8);
			std::memcpy(__mem->_gammaijk, gammaijk, sizeof(double) * 8);

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
		void update_xi_eta_nu(int nNode, KingOfMonsters::myDoubleArray^ xi, KingOfMonsters::myDoubleArray^ eta, KingOfMonsters::myDoubleArray^ nu) {
			if (xi != nullptr)
			{
				__mem->set_buf_xi(xi->_arr->__v.data(), nNode);
			}
			if (eta != nullptr)
			{
				__mem->set_buf_eta(eta->_arr->__v.data(), nNode);
			}
			if (nu != nullptr)
			{
				__mem->set_buf_nu(nu->_arr->__v.data(), nNode);
			}
			else {
				__mem->clear_buf_nu(nNode);
			}
		}
		void update_u_v_xi_eta(int nNode, KingOfMonsters::myDoubleArray^ u, KingOfMonsters::myDoubleArray^ v,KingOfMonsters::myDoubleArray^ xi, KingOfMonsters::myDoubleArray^ eta) {
			if (xi != nullptr)
			{
				__mem->set_buf_xi(xi->_arr->__v.data(), nNode);
			}
			if (eta != nullptr)
			{
				__mem->set_buf_eta(eta->_arr->__v.data(), nNode);
			}
			if (u != nullptr)
			{
				__mem->set_buf_u(u->_arr->__v.data(), nNode);
			}
			if (v != nullptr)
			{
				__mem->set_buf_v(v->_arr->__v.data(), nNode);
			}
		}
		void update_u_v_xi_eta(int nNode, KingOfMonsters::myDoubleArray^ u, KingOfMonsters::myDoubleArray^ v, KingOfMonsters::myDoubleArray^ xi, KingOfMonsters::myDoubleArray^ eta, KingOfMonsters::myDoubleArray^ nu) {
			if (xi != nullptr)
			{
				__mem->set_buf_xi(xi->_arr->__v.data(), nNode);
			}
			if (eta != nullptr)
			{
				__mem->set_buf_eta(eta->_arr->__v.data(), nNode);
			}
			if (nu != nullptr)
			{
				__mem->set_buf_nu(nu->_arr->__v.data(), nNode);
			}
			if (u != nullptr)
			{
				__mem->set_buf_u(u->_arr->__v.data(), nNode);
			}
			if (v != nullptr)
			{
				__mem->set_buf_v(v->_arr->__v.data(), nNode);
			}

		}
		void update_u_v_w_xi_eta(int nNode, KingOfMonsters::myDoubleArray^ u, KingOfMonsters::myDoubleArray^ v, KingOfMonsters::myDoubleArray^ w, KingOfMonsters::myDoubleArray^ xi, KingOfMonsters::myDoubleArray^ eta) {
			if (xi != nullptr)
			{
				__mem->set_buf_xi(xi->_arr->__v.data(), nNode);
			}
			if (eta != nullptr)
			{
				__mem->set_buf_eta(eta->_arr->__v.data(), nNode);
			}
			
			if (u != nullptr)
			{
				__mem->set_buf_u(u->_arr->__v.data(), nNode);
			}
			if (v != nullptr)
			{
				__mem->set_buf_v(v->_arr->__v.data(), nNode);
			}
			if (w != nullptr)
			{
				__mem->set_buf_w(w->_arr->__v.data(), nNode);
			}
			else {
				__mem->clear_buf_w(nNode);
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
		double x, y, z, Z, phi,u,v,w;
		double _x, _y, _z, _Z, _phi,_u,_v,_w;
		double xi, eta,mu,nu;
		double _xi, _eta, _nu;
		double refDv;
		double _refDv;
		
		double load_dv;
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

		
		
		double guideBC( double v1,double v2,bool accurate)
		{
			return __mem->guideBC(v1, v2,accurate);
		}
		void guideBC_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guideBC_xi(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guideBC_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guideBC_eta(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guideBC_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guideBC_nu(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double guideA(double v1, double v2, double w1, double w2)
		{
			return __mem->guideA(v1, v2,w1,w2);
		}

		void guideA_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->guideA_u(__mem->__grad, v1, v2,w1,w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guideA_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->guideA_v(__mem->__grad, v1, v2,w1,w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guideA_w(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->guideA_w(__mem->__grad, v1, v2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double guideB()
		{
			return __mem->guideB();
		}

		void guideB_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->guideB_u(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guideB_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->guideB_v(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guideB_w(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->guideB_w(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		double guideC(double t1, double t2, double s1, double s2, double q1, double q2)
		{
			return __mem->guideC(t1,t2,s1,s2,q1,q2);
		}

		void guideC_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, double s1, double s2, double q1, double q2)
		{
			__mem->guideC_u(__mem->__grad, t1, t2, s1, s2, q1, q2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guideC_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, double s1, double s2, double q1, double q2)
		{
			__mem->guideC_v(__mem->__grad, t1, t2, s1, s2, q1, q2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guideC_w(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, double s1, double s2, double q1, double q2)
		{
			__mem->guideC_w(__mem->__grad, t1, t2, s1, s2, q1, q2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double guide(double v1, double v2,  bool accurate)
		{
			return __mem->guide(v1, v2, accurate);
		}

		void guide_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guide_xi(__mem->__grad, v1, v2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2,  bool accurate)
		{
			__mem->guide_eta(__mem->__grad, v1, v2,  accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guide_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, bool accurate)
		{
			__mem->guide_nu(__mem->__grad, v1, v2,  accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}

		

		
		double guide_trace(double t1, double t2,double w1, double w2, bool accurate)
		{
			return __mem->guide_trace(t1, t2, w1, w2, accurate);
		}
		void guide_trace_xi( myDoubleArray^ vec, myIntArray^ index, double sc, double t1, double t2, double w1, double w2, bool accurate)
		{
			__mem->guide_trace_xi(__mem->__grad, t1, t2,  w1, w2, accurate);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		void guide_trace_eta(myDoubleArray^ vec, myIntArray^ index, double sc, double t1, double t2, double w1, double w2, bool accurate)
		{
			__mem->guide_trace_eta(__mem->__grad, t1, t2,  w1, w2, accurate);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		void guide_trace_nu( myDoubleArray^ vec, myIntArray^ index, double sc, double t1, double t2,double w1, double w2, bool accurate)
		{
			__mem->guide_trace_nu(__mem->__grad, t1, t2,  w1, w2, accurate);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		void guide_trace_xi( mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2,  double w1, double w2, bool accurate)
		{
			__mem->guide_trace_xi(__mem->__grad, t1, t2,  w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_trace_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2,  double w1, double w2, bool accurate)
		{
			__mem->guide_trace_eta( __mem->__grad, t1, t2,  w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guide_trace_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2,  double w1, double w2, bool accurate)
		{
			__mem->guide_trace_nu(__mem->__grad, t1, t2,  w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode,true, c1);
		}

		
		double guide2(double v1, double v2,  double c1, double c2, bool accurate)
		{
			return __mem->guide2(v1, v2, c1, c2, accurate);
		}

		void guide2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double c1, double c2, bool accurate)
		{
			__mem->guide2_xi(__mem->__grad, v1, v2,  c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c);
		}
		void guide2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2,  double c1, double c2, bool accurate)
		{
			__mem->guide2_eta(__mem->__grad, v1, v2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c);
		}
		void guide2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double c1, double c2, bool accurate)
		{
			__mem->guide2_nu(__mem->__grad, v1, v2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c);
		}

		double guide1_2(double v1, double v2, double s1, double s2, bool accurate)
		{
			return __mem->guide1_2(v1, v2, s1, s2, accurate);
		}

		void guide1_2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, bool accurate, bool add)
		{
			__mem->guide1_2_xi(__mem->__grad, v1, v2, s1, s2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void guide1_2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, bool accurate, bool add)
		{
			__mem->guide1_2_eta(__mem->__grad, v1, v2, s1, s2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void guide1_2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, bool accurate, bool add)
		{
			__mem->guide1_2_nu(__mem->__grad, v1, v2, s1, s2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}


		double guide2_2(double v1, double v2,double s1,double s2, double c1, double c2, bool accurate)
		{
			return __mem->guide2_2(v1, v2, s1,s2,c1, c2, accurate);
		}

		void guide2_2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, double c1, double c2, bool accurate,bool add)
		{
			__mem->guide2_2_xi(__mem->__grad, v1, v2, s1,s2,c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide2_2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, double c1, double c2, bool accurate,bool add)
		{
			__mem->guide2_2_eta(__mem->__grad, v1, v2,s1,s2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide2_2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, double c1, double c2,bool accurate, bool add)
		{
			__mem->guide2_2_nu(__mem->__grad, v1, v2,s1,s2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		double __area()
		{
			return __mem->__area();
		}
		void __area_u(myDoubleArray^ vec, myIntArray^ index, double sc)
		{
			__mem->__area_u(__mem->__grad);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		void __area_v(myDoubleArray^ vec, myIntArray^ index, double sc)
		{
			__mem->__area_v(__mem->__grad);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		
		void __area_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c)
		{
			__mem->__area_u(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c);
		}
		void __area_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c)
		{
			__mem->__area_v(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c);
		}
		
		double det()
		{
			return __mem->det();
		}
		double guide_free3(double v1, double v2, double s1, double s2,double globalratio)
		{
			return __mem->guide_free3(v1, v2, s1, s2, globalratio);
		}
		void guide_free3_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free3_phi(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}

		void guide_free3_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free3_xi(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_free3_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free3_eta(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_free3_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free3_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_free3_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free3_v(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		
		double cont_xi(double s1, double s2)
		{
			return __mem->cont_xi(s1,s2);
		}
		void cont_xi_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c,double s1, double s2, bool add, int shift)
		{
			__mem->cont_xi_xi(__mem->__grad,  s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}
		double cont_eta(double s1, double s2)
		{
			return __mem->cont_eta(s1, s2);
		}
		void cont_eta_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double s1, double s2, bool add, int shift)
		{
			__mem->cont_eta_eta(__mem->__grad, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}


		double symm_xi(double v1,double v2,double s1, double s2)
		{
			return __mem->symm_xi(v1,v2,s1, s2);
		}
		void symm_xi_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, int shift)
		{
			__mem->symm_xi_xi(__mem->__grad, v1, v2, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}
		double symm_eta(double v1, double v2, double s1, double s2)
		{
			return __mem->symm_eta(v1, v2, s1, s2);
		}
		void symm_eta_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, int shift)
		{
			__mem->symm_eta_eta(__mem->__grad, v1, v2, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}
		double guide_free4(double v1, double v2, double s1, double s2, double globalratio,bool add,int mode)
		{
			return __mem->guide_free4(v1, v2, s1, s2, globalratio,add,mode);
		}
		void guide_free4_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio,int mode,int shift)
		{
			__mem->guide_free4_phi(__mem->__grad, v1, v2, s1, s2, globalratio,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
		}
		void guide_free4_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio, int mode, int shift)
		{
			__mem->guide_free4_xi(__mem->__grad, v1, v2, s1, s2, globalratio,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
		}
		void guide_free4_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio, int mode, int shift)
		{
			__mem->guide_free4_eta(__mem->__grad, v1, v2, s1, s2, globalratio,mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
		}
		void guide_free4_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free4_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_free4_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_free4_v(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		double free_cond(double v1, double v2)
		{
			return __mem->free_cond(v1, v2);
		}
		

		void free_cond_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2,bool add,int shift)
		{

			__mem->free_cond_xi(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
		}

			void free_cond_eta(mySparse ^ mat, int ii, myIntArray ^ index, double sc, double c, double v1, double v2,bool add, int shift)
			{

				__mem->free_cond_eta(__mem->__grad, v1, v2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}
			void free_cond_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, bool add, int shift)
			{

				__mem->free_cond_phi(__mem->__grad, v1, v2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}
			void free_cond_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, bool add, int shift)
			{

				__mem->free_cond_nu(__mem->__grad, v1, v2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}
			double free_cond2(double v1, double v2,double w1,double w2)
			{
				return __mem->free_cond2(v1, v2,w1,w2);
			}
			
			void free_cond2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2,bool add, int shift)
			{

				__mem->free_cond2_xi(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}

			void free_cond2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2,bool add, int shift)
			{

				__mem->free_cond2_eta(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}
			void free_cond2_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2, bool add, int shift)
			{

				__mem->free_cond2_phi(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}
			void free_cond2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2, bool add, int shift)
			{

				__mem->free_cond2_nu(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, c);
			}


			double free_cond3(double v1, double v2, double w1, double w2)
			{
				return __mem->free_cond3(v1, v2, w1, w2);
			}

			void free_cond3_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2, bool add, int shift)
			{

				__mem->free_cond3_xi(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
			}

			void free_cond3_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2, bool add, int shift)
			{

				__mem->free_cond3_eta(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
			}
			void free_cond3_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double w1, double w2, bool add, int shift)
			{

				__mem->free_cond3_phi(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
			}
		
		double guide_free5(double v1, double v2, double s1, double s2, double globalratio, bool add, int mode)
		{
			return __mem->guide_free5(v1, v2, s1, s2, globalratio, add, mode);
		}
		void guide_free5_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio, int mode, int shift)
		{
			__mem->guide_free5_phi(__mem->__grad, v1, v2, s1, s2, globalratio, mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}
		void guide_free5_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio, int mode, int shift)
		{
			__mem->guide_free5_xi(__mem->__grad, v1, v2, s1, s2, globalratio, mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}
		void guide_free5_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio, int mode, int shift)
		{
			__mem->guide_free5_eta(__mem->__grad, v1, v2, s1, s2, globalratio, mode);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c);
		}
		double guide_Free3(double v1, double v2, double s1, double s2, double globalratio)
		{
			return __mem->guide_Free3(v1, v2, s1, s2, globalratio);
		}
		void guide_Free3_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free3_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free3_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free3_v(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}

		double guide_Free4(double v1, double v2,double s1,double s2,double globalratio)
		{
			return __mem->guide_Free4(v1, v2,s1,s2, globalratio);
		}
		void guide_Free4_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free4_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free4_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free4_v(__mem->__grad, v1, v2,s1,s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}

		double guide_Free5(double v1, double v2, double s1, double s2, double globalratio)
		{
			return __mem->guide_Free5(v1, v2, s1, s2, globalratio);
		}
		void guide_Free5_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free5_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free5_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free5_v(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free5_w(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free5_w(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		double guide_Free6(double v1, double v2, double s1, double s2, double globalratio)
		{
			return __mem->guide_Free6(v1, v2, s1, s2, globalratio);
		}
		void guide_Free6_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free6_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free6_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free6_v(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free6_w(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free6_w(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		double guide_Free7(double v1, double v2, double s1, double s2, double globalratio)
		{
			return __mem->guide_Free7(v1, v2, s1, s2, globalratio);
		}
		void guide_Free7_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free7_u(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_Free7_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
		{
			__mem->guide_Free7_v(__mem->__grad, v1, v2, s1, s2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
			void guide_Free7_w(mySparse ^ mat, int ii, myIntArray ^ index, double sc, double c, double v1, double v2, double s1, double s2, bool add, double globalratio)
			{
				__mem->guide_Free7_w(__mem->__grad, v1, v2, s1, s2, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
			}
		double guide_free2(double v1, double v2, double s1, double s2, double c1, double c2, bool accurate)
		{
			return __mem->guide_free2(v1, v2, s1, s2, c1, c2, accurate);
		}

		void guide_free2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, double c1, double c2, bool accurate, bool add)
		{
			__mem->guide_free2_xi(__mem->__grad, v1, v2, s1, s2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_free2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, double c1, double c2, bool accurate, bool add)
		{
			__mem->guide_free2_eta(__mem->__grad, v1, v2, s1, s2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}
		void guide_free2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c, double v1, double v2, double s1, double s2, double c1, double c2, bool accurate, bool add)
		{
			__mem->guide_free2_nu(__mem->__grad, v1, v2, s1, s2, c1, c2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c);
		}

		double guide_supported2(double t1, double t2, bool accurate)
		{
			return __mem->guide_supported2(t1, t2, accurate);
		}
		void guide_supported2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, bool accurate)
		{
			__mem->guide_supported2_xi(__mem->__grad, t1, t2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_supported2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, bool accurate)
		{
			__mem->guide_supported2_eta(__mem->__grad, t1, t2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guide_supported2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, bool accurate)
		{
			__mem->guide_supported2_nu(__mem->__grad, t1, t2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
	


		double guide_trace2()
		{
			return __mem->guide_trace2();
		}
		void guide_trace2_xi(myDoubleArray^ vec, myIntArray^ index, double sc)
		{
			__mem->guide_trace2_xi( __mem->__grad);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		void guide_trace2_eta(myDoubleArray^ vec, myIntArray^ index, double sc)
		{
			__mem->guide_trace2_eta( __mem->__grad);
			vec->_arr->plus_useindex(__mem->__grad, sc, __mem->_nNode, index->_arr);
		}
		
		void guide_trace2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,bool add)
		{
			__mem->guide_trace2_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void guide_trace2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,bool add)
		{
			__mem->guide_trace2_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode,add, c1);
		}
		

		double guide_traceBC(double n1, double n2, double w1, double w2, bool accurate)
		{
			return __mem->guide_traceBC(n1, n2, w1, w2, accurate);
		}

		void guide_traceBC_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double n1, double n2, double w1, double w2, bool accurate)
		{
			__mem->guide_traceBC_xi(__mem->__grad, n1, n2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_traceBC_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double n1, double n2, double w1, double w2, bool accurate)
		{
			__mem->guide_traceBC_eta(__mem->__grad, n1, n2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guide_traceBC_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double n1, double n2, double w1, double w2, bool accurate)
		{
			__mem->guide_traceBC_nu(__mem->__grad, n1, n2, w1, w2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		
		
		
		double sx(double v1, double v2) {
			return __mem->sx(v1, v2);
		}

		double sy(double v1, double v2) {
			return __mem->sy(v1, v2);
		}

		double symm_sigma(double sx,double sy,double __xi,double __eta)
		{
			return __mem->symm_sigma(sx,sy,__xi,__eta);
		}
		void symm_sigma_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double sx, double sy)
		{
			__mem->symm_sigma_xi(__mem->__grad,sx, sy);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void symm_sigma_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double sx, double sy)
		{
			__mem->symm_sigma_eta(__mem->__grad,  sx, sy);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);

		}
		double symm_gurtin(double sx, double sy, double __xi, double __eta)
		{
			return __mem->symm_gurtin(sx, sy, __xi, __eta);
		}
		void symm_gurtin_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double sx, double sy)
		{
			__mem->symm_gurtin_phi(__mem->__grad, sx, sy);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void symm_gurtin_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double sx, double sy)
		{
			__mem->symm_gurtin_nu(__mem->__grad, sx, sy);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);

		}
		double symm(double sx, double sy, double __x, double __y)
		{
			return __mem->symm(sx, sy, __x, __y);
		}
		void symm_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double sx, double sy)
		{
			__mem->symm_x(__mem->__grad, sx, sy);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void symm_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double sx, double sy)
		{
			__mem->symm_y(__mem->__grad, sx, sy);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);

		}

		double symm2(double v1, double v2)
		{
			return __mem->symm2(v1, v2);
		}
		void symm2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->symm2_xi(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void symm2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->symm2_eta(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double symm3(double v1, double v2)
		{
			return __mem->symm3(v1, v2);
		}
		void symm3_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->symm3_nu(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		
		
		
		double guide_supported(double t1, double t2,bool accurate)
		{
			return __mem->guide_supported(t1, t2,accurate);
		}
		void guide_supported_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, bool accurate)
		{
			__mem->guide_supported_xi(__mem->__grad, t1, t2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void guide_supported_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, bool accurate)
		{
			__mem->guide_supported_eta(__mem->__grad, t1, t2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void guide_supported_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, bool accurate)
		{
			__mem->guide_supported_nu(__mem->__grad, t1, t2, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		



		
		double dir()
		{
			return __mem->dir();
		}
		void dir_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->dir_xi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void dir_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->dir_eta(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		
		double parallel(double t1, double t2, bool add)
		{
			return __mem->parallel(t1, t2, add);
		}
		void parallel_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, int shift, double t1, double t2, bool add)
		{
			__mem->parallel_xi(__mem->__grad, t1, t2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c1);
		}
		void parallel_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, int shift, double t1, double t2, bool add)
		{
			__mem->parallel_eta(__mem->__grad, t1, t2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c1);
		}
		double parallel2(double t1, double t2, bool add)
		{
			return __mem->parallel2(t1, t2, add);
		}
		void parallel2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, int shift, double t1, double t2, bool add)
		{
			__mem->parallel2_xi(__mem->__grad, t1, t2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c1);
		}
		void parallel2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, int shift, double t1, double t2, bool add)
		{
			__mem->parallel2_eta(__mem->__grad, t1, t2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, c1);
		}
		
		double singular()
		{
			return __mem->singular();
		}
		
		void singular_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singular_z( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singular_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singular_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}

		double singular2()
		{
			return __mem->singular2();
		}

		void singular2_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singular2_z( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singular2_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singular2_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}


		double singular3()
		{
			return __mem->singular();
		}

		void singular3_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singular3_z( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singular3_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singular3_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}


		double singularA()
		{
			return __mem->singularA();
		}


		void singularA_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularA_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularA_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularA_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularA_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularA_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void singularA_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularA_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		
		double singularB()
		{
			return __mem->singularB();
		}


		void singularB_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularB_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularB_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void singularB_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double singularB2()
		{
			return __mem->singularB2();
		}


		void singularB2_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB2_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularB2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB2_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularB2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB2_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void singularB2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularB2_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}


		double singularC()
		{
			return __mem->singularC();
		}


		void singularC_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularC_z( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularC_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularC_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularC_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularC_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void singularC_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularC_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		double singularD()
		{
			return __mem->singularD();
		}


		void singularD_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD_z( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularD_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularD_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void singularD_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}


		double singularD2()
		{
			return __mem->singularD2();
		}


		void singularD2_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD2_z( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularD2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD2_xi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void singularD2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD2_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void singularD2_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->singularD2_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double trace()
		{
			return __mem->trace();
		}


		
		void trace_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->trace_xi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void trace_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->trace_eta(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void trace_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->trace_nu(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		double align_sigma()
		{  
			return __mem->align_sigma();
		}
		
		
		void align_sigma_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_z(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void align_sigma_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_phi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void align_sigma_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_x(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void align_sigma_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_y(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void align_sigma_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_xi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void align_sigma_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_eta( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
		void align_sigma_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1)
		{
			__mem->align_sigma_nu( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode,false, c1);
		}
		double align_mix(double v1, double v2,  double w1, double w2)
		{
			return __mem->align_mix( v1, v2,  w1, w2);
		}

		
		void align_mix_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->align_mix_z( __mem->__grad, v1, v2,  w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void align_mix_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2,  double w1, double w2)
		{
			__mem->align_mix_phi( __mem->__grad, v1, v2,  w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}

		void mix22_mat_xi(mySparse^ mat, myIntArray^ index1, myIntArray^ index2, double sc, double v1, double v2, double s1, double s2, double w1, double w2)
		{
			__mem->__mix22_mat_xi(v1,v2,s1,s2,w1,w2);
			for (int s = 0; s < __mem->_nNode; s++)
			{
				for (int t = 0; t < __mem->_nNode; t++)
				{
					mat->dat->adddat(index1->at(s), index2->at(t), __mem->_ref->__matD_xi[s * __mem->_nNode + t] * sc);
				}
			}
		}
		void mix22_mat_eta(mySparse^ mat, myIntArray^ index1, myIntArray^ index2, double sc,  double v1, double v2, double s1, double s2, double w1, double w2)
		{
			__mem->__mix22_mat_eta(v1, v2, s1, s2, w1, w2);
			for (int s = 0; s < __mem->_nNode; s++)
			{
				for (int t = 0; t < __mem->_nNode; t++)
				{
					mat->dat->adddat(index1->at(s), index2->at(t), __mem->_ref->__matD_eta[s * __mem->_nNode + t] * sc);
				}
			}
		}
		void mix22_mat_phi(mySparse^ mat, myIntArray^ index1, myIntArray^ index2, double sc, double v1, double v2, double s1, double s2, double w1, double w2)
		{
			__mem->__mix22_mat_phi(v1, v2, s1, s2, w1, w2);
			for (int s = 0; s < __mem->_nNode; s++)
			{
				for (int t = 0; t < __mem->_nNode; t++)
				{
					mat->dat->adddat(index1->at(s), index2->at(t), __mem->_ref->__matD_phi[s * __mem->_nNode + t] * sc);
				}
			}
		}


		double align_mix22(double v1, double v2, double s1, double s2, double w1, double w2,double globalratio,bool additional)
		{
			return __mem->align_mix22( v1, v2, s1, s2, w1, w2, globalratio,additional);
		}
		void align_mix22_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2,bool add, double globalratio, bool additional)
		{
			__mem->align_mix22_z( __mem->__grad, v1, v2, s1, s2, w1, w2, globalratio, additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		
		void align_mix22_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix22_phi( __mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix22_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix22_xi( __mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix22_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix22_eta( __mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix22_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio, bool additional)
		{
			__mem->align_mix22_x(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio, additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix22_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio, bool additional)
		{
			__mem->align_mix22_y(__mem->__grad, v1, v2, s1, s2, w1, w2,globalratio,additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		

		double align_mix4(double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool additional)
		{
			return __mem->align_mix4(v1, v2, s1, s2, w1, w2, globalratio, additional);
		}
		void align_mix4_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio, bool additional)
		{
			__mem->align_mix4_z(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio, additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}

	
		void align_mix4_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix4_xi(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix4_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix4_eta(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix4_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix4_phi(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix4_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix4_nu(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		double align_mix5(double v1, double v2, double s1, double s2, double w1, double w2, double globalratio, bool additional)
		{
			return __mem->align_mix5(v1, v2, s1, s2, w1, w2, globalratio, additional);
		}
		void align_mix5_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio, bool additional)
		{
			__mem->align_mix5_z(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio, additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}


		void align_mix5_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix5_xi(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix5_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix5_eta(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix5_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add, double globalratio)
		{
			__mem->align_mix5_phi(__mem->__grad, v1, v2, s1, s2, w1, w2, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		
		double align_mix32( double globalratio,bool additional)
		{
			return __mem->align_mix32(globalratio,additional);
		}
		void align_mix32_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,bool add, double globalratio, bool additional)
		{
			__mem->align_mix32_z(__mem->__grad, globalratio,additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}

		void align_mix32_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,  bool add, double globalratio)
		{
			__mem->align_mix32_phi(__mem->__grad,  globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix32_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,  bool add, double globalratio)
		{
			__mem->align_mix32_xi(__mem->__grad, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix32_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, bool add, double globalratio)
		{
			__mem->align_mix32_eta(__mem->__grad, globalratio);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
	
		void align_mix32_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, bool add, double globalratio, bool additional)
		{
			__mem->align_mix32_u(__mem->__grad,  globalratio, additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix32_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,  bool add, double globalratio, bool additional)
		{
			__mem->align_mix32_v(__mem->__grad, globalratio, additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix32_w(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1,  bool add, double globalratio, bool additional)
		{
			__mem->align_mix32_w(__mem->__grad, globalratio,additional);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		double align_metric(double v1, double v2, double s1, double s2, double w1, double w2)
		{
			return __mem->align_metric(v1, v2, s1, s2, w1, w2);
		}

		void align_metric_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add)
		{
			__mem->align_metric_u(__mem->__grad, v1, v2, s1, s2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_metric_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add)
		{
			__mem->align_metric_v(__mem->__grad, v1, v2, s1, s2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		double align_mix2(double v1, double v2, double s1, double s2, double w1, double w2)
		{
			return __mem->align_mix2( v1, v2, s1, s2, w1, w2);
		}

		void align_mix2_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2,bool add)
		{
			__mem->align_mix2_x( __mem->__grad, v1, v2, s1, s2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix2_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add)
		{
			__mem->align_mix2_y( __mem->__grad, v1, v2, s1, s2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix2_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add)
		{
			__mem->align_mix2_z( __mem->__grad, v1, v2, s1,s2,w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
		}
		void align_mix2_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double s1, double s2, double w1, double w2, bool add)
		{
			__mem->align_mix2_phi( __mem->__grad, v1, v2, s1,s2,w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, c1);
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
	
		double BCEQ(double t1, double t2,double s1,double s2, double s11,double a,double b,bool accurate)
		{
			return __mem->BCEQ(t1,t2,s1,s2,s11,a,b,accurate);
		}
		void BCEQ_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2,double s1,double s2,double  stt,double a,double b,bool accurate)
		{
			__mem->BCEQ_phi(__mem->__grad, t1, t2,s1,s2, stt,a,b, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void BCEQ_Z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2,double s1,double s2,double stt,double a,double b,bool accurate)
		{
			__mem->BCEQ_z(__mem->__grad, t1,t2,s1,s2, stt, a,b,accurate);
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
		double st(double t1, double t2, double n1, double n2)
		{

			return __mem->st(t1, t2, n1, n2);
		}
		void st_phi(myDoubleArray^ arr, double t1, double t2, double n1, double n2)
		{

			__mem->st_phi(arr->_arr->__v.data(), t1, t2, n1, n2);

		}
		double cont2(memS^ other, double w1, double w2,double w1_2,double w2_2) {
			return __mem->cont2(other->__mem, w1, w2,w1_2,w2_2);
		}
		void cont2_u(mySparse^ mat, memS^ other, int ii, myIntArray^ index, int shift, double sc, double c1, double w1, double w2, double w1_2, double w2_2,bool left,bool right)
		{
			__mem->cont2_u(other->__mem, __mem->__grad, other->__mem->__grad, w1, w2, w1_2, w2_2);
			if (left)
			{
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
				if (right)
				{
					mat->dat->addrow(ii, index->_arr, other->__mem->__grad - shift, shift, sc, __mem->_nNode, false, c1);
				}
			}else
			if (right)
			{
				mat->dat->addrow(ii, index->_arr, other->__mem->__grad - shift, shift, sc, __mem->_nNode, true, c1);
			}
		}
		void cont2_v(mySparse^ mat,memS^ other, int ii, myIntArray^ index, int shift, double sc, double c1, double w1, double w2, double w1_2, double w2_2, bool left, bool right)
		{
			__mem->cont2_v(other->__mem, __mem->__grad, other->__mem->__grad, w1, w2, w1_2, w2_2);
			if (left)
			{
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
			}
			if(right){
				mat->dat->addrow(ii, index->_arr, other->__mem->__grad - shift, shift, sc, __mem->_nNode, false, c1);
			}
		}
		double fair3(double v1, double v2, double w1, double w2)
		{

			return __mem->fair3(v1,v2,w1,w2);
		}

		void fair3_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->fair3_u(__mem->__grad,v1,v2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void fair3_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->fair3_v(__mem->__grad,v1,v2, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		double fair4(double v1, double v2, double w1, double w2)
		{

			return __mem->fair4(v1, v2);
		}

		void fair4_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->fair4_u(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, c1);
		}
		void fair4_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2)
		{
			__mem->fair4_v(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		double length(double v1, double v2)
		{
			return __mem->length(v1, v2);
		}
		void length_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->length_u(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void length_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->length_v(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}

		double gamma(double v1, double v2)
		{
			return __mem->gamma(v1, v2);
		}
		void gamma_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->gamma_u(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
		void gamma_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2)
		{
			__mem->gamma_v(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
		}
			double fair2(double v1, double v2, double w1, double w2, double s1, double s2)
			{

				return __mem->fair2(v1, v2, w1, w2, s1, s2);
			}
			void fair2_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, double s1, double s2)
			{
				__mem->fair2_u(__mem->__grad, v1, v2, w1, w2, s1, s2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
			}
			void fair2_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double v1, double v2, double w1, double w2, double s1, double s2)
			{
				__mem->fair2_v(__mem->__grad, v1, v2, w1, w2, s1, s2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, c1);
			}
		double st2(double t1,double t2,double n1,double n2)
		{

			return __mem->st2(t1, t2, n1, n2);
		}
		void st2_z(myDoubleArray^ arr, double t1, double t2, double n1, double n2)
		{

			__mem->st2_z(arr->_arr->__v.data(), t1, t2, n1, n2);

		}
		void st2_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, double n1, double n2)
		{
			__mem->st2_z(__mem->__grad, t1, t2, n1, n2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, c1);
		}
			void st2_x(mySparse ^ mat, int ii, myIntArray ^ index, double sc, double c1, double t1, double t2, double n1, double n2)
			{
				__mem->st2_x(__mem->__grad, t1, t2, n1, n2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad,0, sc, __mem->_nNode, false,c1);
			}
			void st2_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double t1, double t2, double n1, double n2)
			{
				__mem->st2_y(__mem->__grad, t1, t2, n1, n2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0,sc, __mem->_nNode,false,c1);
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

		void d0(mySparse^ mat, int ii, myIntArray^ index, double sc,double c1,bool add)
		{
			mat->dat->addrow(ii, index->_arr, __mem->_ref->d0,0, sc, __mem->_nNode,add,c1);
		}
		void d0(mySparse^ mat, int ii, myIntArray^ index,int shift, double sc, double c1, bool add)
		{
			mat->dat->addrow(ii, index->_arr, __mem->_ref->d0-shift,shift,  sc, __mem->_nNode, add, c1);
		}
		void U_z(mySparse^ mat, int ii, myIntArray^ index, double sc)
		{
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode);
		}
		void U_mix(mySparse^ mat, int ii, myIntArray^ index, double sc, double c1, double c2)
		{
			mat->dat->addrow(ii, index->_arr, __mem->__grad, __mem->__grad, sc, __mem->_nNode, c1, c2);
		}
		
			
			double gauss()
			{
				return __mem->gauss();
			}
			void gauss_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->gauss_nu( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, coeff);
			}
			double mean()
			{
				return __mem->mean();
			}
			void mean_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->mean_nu( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, coeff);
			}
			void bodyF2_mat_xi(mySparse^ mat, myIntArray^ index1,myIntArray^ index2,double sc)
			{
				__mem->__bodyF2_mat_xi();
				for (int s = 0; s < __mem->_nNode; s++)
				{
					for (int t = 0; t < __mem->_nNode; t++)
					{
						mat->dat->adddat(index1->at(s), index2->at(t), __mem->_ref->__matF_xi[s*__mem->_nNode+t]*sc);
					}
				}
			}
			void bodyF2_mat_eta(mySparse^ mat, myIntArray^ index1, myIntArray^ index2,double sc)
			{
				__mem->__bodyF2_mat_eta();
				for (int s = 0; s < __mem->_nNode; s++)
				{
					for (int t = 0; t < __mem->_nNode; t++)
					{
						mat->dat->adddat(index1->at(s), index2->at(t), __mem->_ref->__matF_eta[s * __mem->_nNode + t]*sc);
					}
				}
			}
			void bodyF2_mat_phi (mySparse^ mat, myIntArray^ index1, myIntArray^ index2, double sc)
			{
				__mem->__bodyF2_mat_phi();
				for (int s = 0; s < __mem->_nNode; s++)
				{
					for (int t = 0; t < __mem->_nNode; t++)
					{
						mat->dat->adddat(index1->at(s), index2->at(t), __mem->_ref->__matF_phi[s * __mem->_nNode + t] * sc);
					}
				}
			}
			double bodyF2( double globalratio)
			{
				return __mem->__bodyF2( globalratio);
			}
				double bodyF2_rhs(double load, bool accurate)
				{
					return __mem->__bodyF2_rhs(load, accurate);
				}
			void bodyF2_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate,bool add, double globalratio)
			{
				__mem->__bodyF2_z( __mem->__grad, load, accurate, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad,0, sc, __mem->_nNode,add, coeff);
			}
			
			void bodyF2_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF2_phi( __mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad,0, sc, __mem->_nNode, add,coeff);
			}
			void bodyF2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF2_xi( __mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0,sc, __mem->_nNode,add,coeff);
			}
			void bodyF2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF2_eta( __mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}

			void bodyF2_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff,double load, bool accurate, bool add, double globalratio)
			{
				__mem->__bodyF2_x( __mem->__grad,load,accurate, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF2_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate, bool add, double globalratio)
			{
				__mem->__bodyF2_y( __mem->__grad, load, accurate, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}

			double bodyF4(double globalratio)
			{
				return __mem->__bodyF4(globalratio);
			}
			double bodyF4_rhs(double load, bool accurate)
			{
				return __mem->__bodyF4_rhs(load, accurate);
			}
			void bodyF4_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate, bool add, double globalratio)
			{
				__mem->__bodyF4_z(__mem->__grad, load, accurate, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}

			
			void bodyF4_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF4_xi(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF4_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF4_eta(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF4_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF4_phi(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF4_nu(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF4_nu(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
	
			double bodyF5(double globalratio)
			{
				return __mem->__bodyF5(globalratio);
			}
			double bodyF5_rhs(double load, bool accurate)
			{
				return __mem->__bodyF5_rhs(load, accurate);
			}
			void bodyF5_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate, bool add, double globalratio)
			{
				__mem->__bodyF5_z(__mem->__grad, load, accurate, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}


			void bodyF5_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF5_xi(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF5_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF5_eta(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF5_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add, double globalratio)
			{
				__mem->__bodyF5_phi(__mem->__grad, globalratio);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			
			double bodyU(double load, bool accurate)
			{
				return __mem->__bodyU(load, accurate);
			}

			void bodyU_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add)
			{
				__mem->__bodyU_phi(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyU_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate, bool add)
			{
				__mem->__bodyU_u(__mem->__grad, load, accurate);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			double bodyV(double load, bool accurate)
			{
				return __mem->__bodyV(load, accurate);
			}
			void bodyV_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add)
			{
				__mem->__bodyV_phi(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyV_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate, bool add)
			{
				__mem->__bodyV_v(__mem->__grad, load, accurate);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}


			double bodyF3(double load, bool accurate)
			{
				return __mem->__bodyF3(load, accurate);
			}
			void bodyF3_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate,bool add)
			{
				__mem->__bodyF3_z(__mem->__grad, load, accurate);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0,sc, __mem->_nNode, add,coeff);
			}

			void bodyF3_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add)
			{
				__mem->__bodyF3_phi(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			
			void bodyF3_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate, bool add)
			{
				__mem->__bodyF3_x(__mem->__grad, load, accurate);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void bodyF3_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate,bool add)
			{
				__mem->__bodyF3_y(__mem->__grad, load, accurate);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			
			double ortho(double v1, double v2, double w1, double w2)
			{
				return __mem->ortho(v1, v2, w1, w2);
			}
			
			void ortho_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2,int shift)
			{
				__mem ->ortho_x(__mem->__grad,v1,v2,w1,w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift,shift, sc, __mem->_nNode, true, coeff);
			}
			void ortho_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2,int shift)
			{
				__mem->ortho_y(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, false, coeff);
			}
			double ortho2(double v1, double v2, double w1, double w2)
			{
				return __mem->ortho2(v1, v2, w1, w2);
			}

			
			void ortho2_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2, int shift)
			{
				__mem->ortho2_x(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, true, coeff);
			}
			void ortho2_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2, int shift)
			{
				__mem->ortho2_y(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, false, coeff);
			}

			double ortho3(double v1,double v2,double w1,double w2)
			{
				return __mem->ortho3(v1,v2,w1,w2);
			}


			void ortho3_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2,int shift)
			{
				__mem->ortho3_xi(__mem->__grad, v1, v2,w1,w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, true, coeff);
			}
			void ortho3_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2, int shift)
			{
				__mem->ortho3_eta(__mem->__grad, v1, v2,w1,w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, false, coeff);
			}
			
			double ortho4(double v1, double v2, double w1, double w2)
			{
				return __mem->ortho4(v1, v2, w1, w2);
			}

			void ortho4_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2, int shift)
			{
				__mem->ortho4_xi(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, true, coeff);
			}
			void ortho4_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, double w1, double w2, int shift)
			{
				__mem->ortho4_eta(__mem->__grad, v1, v2, w1, w2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, false, coeff);
			}
		
			double curl_free()
			{
				return __mem->curl_free();
			}

			void curl_free_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff,bool add)
			{
				__mem->curl_free_x(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void curl_free_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff,bool add)
			{
				__mem->curl_free_y(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			double curl_free2()
			{
				return __mem->curl_free2();
			}
			void curl_free2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add)
			{
				__mem->curl_free2_xi(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void curl_free2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, bool add)
			{
				__mem->curl_free2_eta(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}


			double div_free(double v1,double v2)
			{
				return __mem->div_free(v1,v2);
			}
			void div_free_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2,bool add)
			{
				__mem->div_free_xi(__mem->__grad, v1, v2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void div_free_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, bool add)
			{
				__mem->div_free_eta(__mem->__grad, v1, v2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			void div_free_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, bool add)
			{
				__mem->div_free_phi(__mem->__grad, v1, v2);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
			}
			double harmonic_x()
			{
				return __mem->harmonic_x();
			}
			void harmonic_x_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->harmonic_x_xi( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode,true, coeff);
			}
			
			double harmonic_y()
			{
				return __mem->harmonic_y();
			}
			
			void harmonic_y_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->harmonic_y_eta( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			
			
			double harmonic_X()
			{
				return __mem->harmonic_X();
			}
			
			
			void harmonic_X_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->harmonic_X_u( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void harmonic_X_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->harmonic_X_v(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}

			double harmonic_Y()
			{
				return __mem->harmonic_Y();
			}
			
			void harmonic_Y_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->harmonic_Y_u(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void harmonic_Y_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->harmonic_Y_v( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}
			
			double rot_free()
			{
				return __mem->rot_free();
			}
			void rot_free_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->rot_free_u(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void rot_free_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->rot_free_v(__mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}
			double conformal_x()
			{
				return __mem->conformal_x();
			}
			void conformal_x_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_x_xi( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void conformal_x_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_x_eta( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}
			double conformal_y()
			{
				return __mem->conformal_y();
			}
			void conformal_y_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_y_xi( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void conformal_y_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_y_eta( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}

			double conformal_X()
			{
				return __mem->conformal_X();
			}
			void conformal_X_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_X_x( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void conformal_X_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_X_y( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}
			double conformal_Y()
			{
				return __mem->conformal_Y();
			}
			void conformal_Y_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_Y_x( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
			}
			void conformal_Y_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
			{
				__mem->conformal_Y_y( __mem->__grad);
				mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
			}
			
		double bodyF(double load, bool accurate) 
		{			
				return __mem->__bodyF( load, accurate);
		}

		void bodyF_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate)
		{
			__mem->__bodyF_x( __mem->__grad, load, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
		}
		void bodyF_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double load, bool accurate)
		{
			__mem->__bodyF_y( __mem->__grad, load, accurate);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
		}
	    void bodyF_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff,  double load,bool accurate)
		{
			__mem->__bodyF_z(__mem->__grad,load, accurate);			
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, coeff);
		}
		void bodyF_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
		{
			__mem->__bodyF_phi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, coeff);
		}
		void _detZ_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
		{
			__mem->_detZ_z(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, coeff);
		}
		void _detphi_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff)
		{
			__mem->_detphi_phi( __mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, sc, __mem->_nNode, coeff);
		}
		double sc()
		{
			return __mem->__sc();
		}
		double det_sigma()
		{
			return __mem->det_sigma();
		}
		double Suu()
		{
			return __mem->Suu_sigma();
		}
		double Suv()
		{
			return __mem->Suv_sigma();
		}
		double Svv()
		{
			return __mem->Svv_sigma();
		}
		double fuu()
		{
			return __mem->fuu_sigma();
		}
		double fuv()
		{
			return __mem->fuv_sigma();
		}
		double fvv()
		{
			return __mem->fvv_sigma();
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
		double _SLOPE_phi(int l) {
			return __mem->____SLOPE_phi(l);
		}
		double _SLOPE_z(int l) {
			return __mem->____SLOPE_z(l);
		}
		double SLOPE_phi(int l) {
			return __mem->___SLOPE_phi(l);
		}
		double SLOPE_z(int l) {
			return __mem->___SLOPE_z(l);
		}
		double SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {
			return __mem->SLOPE_phi(dcdtstar0, dcdtstar1, sc);
		}
		double SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {
			return __mem->SLOPE_z(dcdtstar0, dcdtstar1, sc);
		}
		double _SLOPE_phi(double dcdtstar0, double dcdtstar1, double sc) {
			return __mem->___SLOPE_phi(dcdtstar0, dcdtstar1, sc);
		}
		double _SLOPE_z(double dcdtstar0, double dcdtstar1, double sc) {
			return __mem->___SLOPE_z(dcdtstar0, dcdtstar1, sc);
		}
		
		double _SLOPE_xi(double dcdtstar0, double dcdtstar1) {
			return __mem->_SLOPE_xi(dcdtstar0, dcdtstar1);
		}
		double _SLOPE_eta(double dcdtstar0, double dcdtstar1) {
			return __mem->_SLOPE_eta(dcdtstar0, dcdtstar1);
		}
		double _SLOPE_nu(double dcdtstar0, double dcdtstar1) {
			return __mem->_SLOPE_nu(dcdtstar0, dcdtstar1);
		}
		double _SLOPE_u(double dcdtstar0, double dcdtstar1) {
			return __mem->_SLOPE_u(dcdtstar0, dcdtstar1);
		}
		double _SLOPE_v(double dcdtstar0, double dcdtstar1) {
			return __mem->_SLOPE_v(dcdtstar0, dcdtstar1);
		}
		double SLOPE_xi(double dcdtstar0, double dcdtstar1) {
			return __mem->SLOPE_xi(dcdtstar0, dcdtstar1);
		}
		double SLOPE_eta(double dcdtstar0, double dcdtstar1) {
			return __mem->SLOPE_eta(dcdtstar0, dcdtstar1);
		}
		double SLOPE_nu(double dcdtstar0, double dcdtstar1) {
			return __mem->SLOPE_nu(dcdtstar0, dcdtstar1);
		}
		double _SLOPE3phiX()
		{
			return __mem->_SLOPE3phiX();
		}
		void _SLOPE3phiX_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, int shift, bool add)
		{

			__mem->_SLOPE3phiX_phi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, coeff);
		}
		double _SLOPE3phiY()
		{
			return __mem->_SLOPE3phiY();
		}
		void _SLOPE3phiY_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, int shift, bool add)
		{

			__mem->_SLOPE3phiY_phi(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, coeff);
		}
		double _SLOPE3zX()
		{
			return __mem->_SLOPE3zX();
		}
		void _SLOPE3zX_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, int shift, bool add)
		{

			__mem->_SLOPE3zX_z(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, coeff);
		}
		double _SLOPE3zY()
		{
			return __mem->_SLOPE3zY();
		}
		void _SLOPE3zY_z(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff,int shift,bool add)
		{
			
				__mem->_SLOPE3zY_z(__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, coeff);
		}
		double _SLOPE2_phi(double v1, double v2)
		{
			return __mem->_SLOPE2_phi(v1, v2);
		}
		double _SLOPE2_z(double v1, double v2)
		{
			return __mem->_SLOPE2_z(v1, v2);
		}
		double _SLOPE2_xi(double v1, double v2)
		{
			return __mem->_SLOPE2_xi(v1, v2);
		}
		double _SLOPE2_eta(double v1, double v2)
		{
			return __mem->_SLOPE2_eta(v1, v2);
		}
		double _SLOPE2_nu(double v1, double v2)
		{
			return __mem->_SLOPE2_nu(v1, v2);
		}
		double _SLOPE2_u(double v1, double v2)
		{
			return __mem->_SLOPE2_u(v1, v2);
		}
		double _SLOPE2_v(double v1, double v2)
		{
			return __mem->_SLOPE2_v(v1, v2);
		}
			double _SLOPE2_w(double v1, double v2)
			{
				return __mem->_SLOPE2_w(v1, v2);
			}
		void _SLOPE2(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2,int shift, bool add)
		{
			__mem->_SLOPE2(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift,shift, sc, __mem->_nNode, add, coeff);
		}
		void _SLOPE2(myDoubleArray^ arr, double v1, double v2,int shift)
		{
			__mem->_SLOPE2(arr->_arr->__v.data()+shift, v1, v2);
			
		}

		void _SLOPE2_phi_x(myDoubleArray^ arr, double v1, double v2, int shift)
		{
			__mem->_SLOPE2_phi_x(arr->_arr->__v.data(), v1, v2);
	
		}
		void _SLOPE2_phi_y(myDoubleArray^ arr, double v1, double v2, int shift)
		{
			__mem->_SLOPE2_phi_y(arr->_arr->__v.data(), v1, v2);
		}

		double _SLOPE_phi3(double v1, double v2)
		{
			return __mem->___SLOPE_phi3(v1, v2);
		}
		void _SLOPE_phi3_phi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, int shift, bool add)
		{
			__mem->___SLOPE_phi3_phi(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_nNode, add, coeff);
		}
		
		void _SLOPE2_phi_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, int shift,bool add)
		{
			__mem->_SLOPE2_phi_x(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, coeff);
		}
		void _SLOPE2_phi_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, int shift,bool add)
		{
			__mem->_SLOPE2_phi_y(__mem->__grad , v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, coeff);
		}
		void _SLOPE2_z_x(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, int shift,bool add)
		{
			__mem->_SLOPE2_z_x(__mem->__grad, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, coeff);
		}
		void _SLOPE2_z_y(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double v1, double v2, int shift, bool add)
		{
			__mem->_SLOPE2_z_y(__mem->__grad , v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_nNode, add, coeff);
		}
		double SLOPE(int l, int i) {
			return __mem->SLOPE(l, i);
		}
		double _SLOPE(int l, int i) {
			return __mem->_SLOPE(l, i);
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
		double  _SLOPE_symm_sigma(double w1, double w2) {
			return __mem->_SLOPE_symm_sigma(w1, w2);
		}
		void _SLOPE_symm_sigma_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2) {
			__mem->_SLOPE_symm_sigma_xi(__mem->__grad, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
		}


		void _SLOPE_symm_sigma_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2) {
			__mem->_SLOPE_symm_sigma_eta(__mem->__grad, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
		}
		double  _SLOPE_symm_sigma2(double w1, double w2,double s1,double s2) {
			return __mem->_SLOPE_symm_sigma2(w1, w2, s1, s2);
		}
		void _SLOPE_symm_sigma2_xi(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2, double s1, double s2, bool add) {
			__mem->_SLOPE_symm_sigma2_xi(__mem->__grad, w1, w2, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
		}


		void _SLOPE_symm_sigma2_eta(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2, double s1, double s2, bool add) {
			__mem->_SLOPE_symm_sigma2_eta(__mem->__grad, w1, w2, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
		}

		void _SLOPE_symm_sigma2_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2, double s1, double s2,bool add) {
			__mem->_SLOPE_symm_sigma2_u(__mem->__grad, w1, w2, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode,add, coeff);
		}


		void _SLOPE_symm_sigma2_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2, double s1, double s2,bool add) {
			__mem->_SLOPE_symm_sigma2_v(__mem->__grad, w1, w2, s1, s2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
		}

		double  _SLOPE_symm2(double w1, double w2) {
			return __mem->_SLOPE_symm2(w1, w2);
		}
		void _SLOPE_symm2_u(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2) {
			__mem->_SLOPE_symm2_u(__mem->__grad, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, true, coeff);
		}


		void _SLOPE_symm2_v(mySparse^ mat, int ii, myIntArray^ index, double sc, double coeff, double w1, double w2) {
			__mem->_SLOPE_symm2_v(__mem->__grad, w1, w2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, false, coeff);
		}
		void _SLOPE_phi_x(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc) {
			__mem->_SLOPE_phi_x(arr->_arr->__v.data(), dcdt1, dcdt2, sc);
		}
		void _SLOPE_phi_y(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc) {
			__mem->_SLOPE_phi_y(arr->_arr->__v.data(), dcdt1, dcdt2, sc);
		}

		void _SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2,double sc) {
			__mem->_SLOPE(arr->_arr->__v.data(), dcdt1, dcdt2,sc);
			//return __mem->SLOPE(l, i);
		}
		void _SLOPE(myDoubleArray^ arr, double dcdt1, double dcdt2, double sc, int shift) {
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
		void laplacian(mySparse^ mat, myIntArray^ index, double sc)
		{
			mat->dat->_nt = 1;
			__mem->laplacian(mat->dat, index->data(), sc);
		}
		void d0(myDoubleArray^ arr)
		{
			memcpy(arr->_arr->__v.data(), this->__mem->_ref->d0, sizeof(double) * __mem->_ref->_nNode);

		}
			void d0(myDoubleArray ^ arr,int shift)
			{
				memcpy(arr->_arr->__v.data()+shift, this->__mem->_ref->d0, sizeof(double) * __mem->_ref->_nNode);
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
		double _Dc(double a, double b, int l) {
			return __mem->_Dc(a, b, l);
		}
		double get_a(int l) {
			return __mem->__a(l);
		}
		double get_b(int l) {
			return __mem->__b(l);
		}
		double get__a(int l) {
			return __mem->___a(l);
		}
		double get__b(int l) {
			return __mem->___b(l);
		}
		double rot1(double v1, double v2)
		{
			return __mem->rot1(v1, v2);
		}
		double rot2(double v1, double v2)
		{
			return __mem->rot2(v1, v2);
		}
		double _Dc2(double a, double b, double v1,double v2) {
			return __mem->_Dc2(a, b, v1,v2);
		}
		void _Dc2_x(mySparse^ mat, int ii, myIntArray^ index, double a, double b, double v1, double v2, double sc, double coeff, bool add)
		{
			__mem->_Dc2_x(__mem->__grad, a, b, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
		}
		void _Dc2_y(mySparse^ mat, int ii, myIntArray^ index, double a, double b, double v1, double v2, double sc, double coeff, bool add)
		{
			__mem->_Dc2_y(__mem->__grad, a, b, v1, v2);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_nNode, add, coeff);
		}

		double get__a2(double v1, double v2) {
			return __mem->___a2(v1,v2);
		}
		double get__b2(double v1, double v2) {
			return __mem->___b2(v1,v2);
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
		double get_sgima_dv()
		{
			return __mem->sigma_dv;
		}
		double get_sigma_trace()
		{
			return __mem->sigma_trace;
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
		
		void bodyF_freeze()
		{
			this->load_dv = __mem->_ref->load_dv;
			
		}
		void bodyF_update()
		{
			__mem->_ref->load_dv = __mem->dv;
			this->load_dv = __mem->_ref->load_dv;
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
			compute("FULL");
		}
		double compute(String ^simple) {		
			if (__mem->RAM == SAVE)
			{
				__mem->_ref->__mat = __mem->__mat;
				__mem->_ref->d0 = __mem->d0;
				
				__mem->_ref->d1 = __mem->d1;
				__mem->_ref->d2 = __mem->d2;
				__mem->_ref->d3 = __mem->d3;
				__mem->_ref->d2_star = __mem->d2_star;
				__mem->_ref->B = __mem->B;
				__mem->_ref->tt0 = __mem->tt0;
				__mem->_ref->hh0 = __mem->hh0;
				__mem->_ref->tt1 = __mem->tt1;
				__mem->_ref->hh1 = __mem->hh1;
				__mem->_ref->tt2 = __mem->tt2;
				__mem->_ref->hh2 = __mem->hh2;
				__mem->_ref->tt3 = __mem->tt3;
				__mem->_ref->hh3 = __mem->hh3;

			}
			std::chrono::high_resolution_clock::time_point begin = high_resolution_clock::now();

			
			
				if(simple=="SIMPLE")
					__mem->update2("SIMPLE");
				else if (simple == "FULL")
					__mem->update2( "FULL");
				else
					__mem->update2( "STANDARD");

		
			// 3. 現在日時を再度取得
			high_resolution_clock::time_point end = high_resolution_clock::now();

			// 経過時間を取得
			microseconds elapsed_time = duration_cast<microseconds>(end - begin);
			

			
			

			x = __mem->x;
			y = __mem->y;
			z = __mem->z;
			Z = __mem->Z;
			u = __mem->u;
			v = __mem->v;
			w = __mem->w;
			phi = __mem->phi;
			dv = __mem->dv;
			_dv = __mem->_dv;
			
			_refDv = __mem->_dv;
			refDv = __mem->dv;
			
			//nu = __mem->nu;
			return elapsed_time.count();
		}
		double __orefDv()
		{
			return __mem->_ref->orefDv;
		}
		void update_optional()
		{
			__mem->update_optional();
			u = __mem->u;
			v = __mem->v;
			w = __mem->w;
			xi = __mem->xi;
			eta = __mem->eta;
			nu = __mem->nu;
			
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
		void memory2(memS_ref^ mem)
		{
			__mem-> memory2(mem->__mem);
		}
		void memory(memS_ref^ _mem_) {
			__mem->memory(_mem_->__mem);
			_x = __mem->_ref->_x;
			_y = __mem->_ref->_y;
			_z = __mem->_ref->_z;
			_u = __mem->_ref->_u;
			_v = __mem->_ref->_v;
			_w = __mem->_ref->_w;
			Z = __mem->_ref->Z;
			_xi = __mem->_ref->_xi;
			_eta = __mem->_ref->_eta;
			_nu = __mem->_ref->_nu;;

			//phi = __mem->phi;
		}
		memS(System::String^ RAM) {
			_RAM __RAM = SAVE;
			if (RAM == "SAVE")__RAM = SAVE;
			if (RAM == "MAX")__RAM = MAX;
			__mem = new _memS(std::string(""), __RAM);//MAX=cosume memory more but fast, SAVE=save memory and a little slower
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
