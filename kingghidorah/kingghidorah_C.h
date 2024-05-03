#pragma once
#include <cmath>

#include<vector>
using namespace System;
#include <cstring>
using std::vector;

#include "mySparseLibrary.h"
#include "kingghidorah_S.h"
namespace KingOfMonsters {
	public class _memC_ref {
	public:
		int _nNode;
		
		double* node=0;
		double* def = 0;
		double* buf_z = 0;
		double* buf_phi = 0;
		//double* buf_b = 0;
		//double* buf_D = 0;
		double* buf_W = 0;

		double _gi[3];
		double _Gi[3];
		double _gij[1];
		double _Gij[1];
		double _bij[3];
		double _Sij[1];
		double _Gammaijk[1];
		const double _up[3]{ 0,0,1 };
	public:
		double refDv,_refDv;
		double _x, _y, _z,__z,w;
		inline void set_z(double z) {
			this->_z = z;
		}
		inline void set__z(double z) {
			this->__z = z;
		}
		inline void set_buffer(double* buf)
		{
			node = buf;
			def = &buf[2000];
			buf_z = &buf[4000];
			buf_phi = &buf[6000];
			//buf_b = &buf[8000];
			//buf_D = &buf[10000];
			buf_W = &buf[8000];
		}
		inline void set_node(int i, int s, double val) {
			node[i * 3 + s] = val;
		}
		inline void set_buf_z(int i, double val) {
			buf_z[i] = val;
		}
		inline void set_buf_phi(int i, double val) {
			buf_phi[i] = val;
		}
		inline void set_buf_W(int i, double val) {
			buf_W[i] = val;
		}
		inline void set_def(int i, int s, double val) {
			def[i * 3 + s] = val;
		}
		inline double get_node(int i, int s) {
			return node[i * 3 + s];
		}
		inline double get__gt(int s) {
			return _gi[s];
		}
		inline double get__Gt(int s) {
			return _Gi[s];
		}
		inline double get__gtt() {
			return _gij[0];
		}
		inline double get__Gtt() {
			return _Gij[0];
		}
		inline double get__btt(int s) {
			return _bij[s];
		}
		inline double get__Gammattt() {
			return _Gammaijk[0];
		}
		inline double get_up(int i) {
			return _up[i];
		}
	public:
		_memC_ref() {
			
			node = 0;
			def = 0;
			_nNode = 0;

			__z = -10000000;
		}
		~_memC_ref() {
			del();
		}
		
		void del(){
			_nNode = 0;

			__z = -10000000;
		}
		void update(int nNode, int dim) {
			if (_nNode != nNode) {
				_nNode = nNode;
			}
		}
	};
	public class _memC {
	public:
		_memC_ref* _ref = 0;
	private:

		double** M = 0;
		double lo;
		int* dd = 0;

		int _nNode;
		int _dim;
		
		double gi[3];
		double gi2[3];
		double Gi[3];
		double Gi2[3];
		double gij[1];
		double gij2[1];
		double Gij[1];
		double Gij2[1];
		double bij[3];
		double bij2[3];
		double Gammaijk[1];
		double Gammaijk2[1];
		double _ss[4];
		double _Sij[1];
		double N[3], H[3];
		double* d0 = 0;
		double* d1 = 0;
		double* d2 = 0;
		double* B = 0;
		double* tt0 = 0, * hh0 = 0, * tt1 = 0, * hh1 = 0, * tt2 = 0, * hh2 = 0;
		double* gradN[3]{ 0,0,0 };
		double* gradG = 0;
	public:
		double x, y, z;
		double dv, _dv;
		double* __grad = 0;
	public:
		inline void set_lo(double L1) {
			lo = L1;
		}
		inline void set_M(int j, int k, double val) {
			M[j][k] = val;
		}
		inline void set_dd(int i, int val) {
			dd[i] = val;
		}
		inline void set_L(double _L) {

			_Sij[0] = _L / this->_dv / this->_dv;

		}

		inline double get_gt(int s) {
			return gi[s];
		}
		inline double get_gt2(int s) {
			return gi2[s];
		}
		inline double get_Gt(int s) {
			return Gi[s];
		}
		inline double get_Gt2(int s) {
			return Gi2[s];
		}
		inline double get_gtt() {
			return gij[0];
		}
		inline double get_gtt2() {
			return gij2[0];
		}
		inline double get_Gtt() {
			return Gij[0];
		}
		inline double get_Gtt2() {
			return Gij2[0];
		}
		inline double get_btt(int s) {
			return bij[s];
		}
		inline double get_Gammattt() {
			return Gammaijk[0];
		}

		inline double get_tt0(int k) {
			return tt0[k];
		}
		inline double get_hh0(int k) {
			return hh0[k];
		}
		inline double get_tt1(int k) {
			return tt1[k];
		}
		inline double get_hh1(int k) {
			return hh1[k];
		}
		inline double get_tt2(int k) {
			return tt2[k];
		}
		inline double get_hh2(int k) {
			return hh2[k];
		}

	public:
		double __hh0(int k) {
			double t = lo;
			return pow(t, (_dim - k - 1));
		}
		double __tt0(int k) {
			double val = 0;
			double t = lo;
			for (int l = 0; l < _dim; l++)
			{
				val += get_hh0(l) * M[l][k];
			}
			return val;
		}
		double __hh1(int k) {
			double t = lo;
			if (k < _dim - 1)
			{
				return (_dim - k - 1) * pow(t, (_dim - k - 2));
			}
			else {
				return 0;
			}
		}
		double __tt1(int k) {
			double val = 0;
			double t = lo;

			for (int l = 0; l < _dim; l++)
			{

				val += get_hh1(l) * M[l][k];
			}
			return val;
		}
		double __hh2(int k) {
			double t = lo;
			if (k < _dim - 2)
				return (_dim - k - 1) * (_dim - k - 2) * pow(t, (_dim - k - 3));
			else return  0;
		}
		double __tt2(int k) {
			double val = 0;
			for (int l = 0; l < _dim; l++)
			{
				val += get_hh2(l) * M[l][k];
			}
			return val;
		}
		double _shape(int k)
		{
			double shape = 1.0;
			shape *= get_tt0(dd[k]);
			//return shape[k];
			return shape;
		}
		double _C(int k)
		{
			double C = 1.0;
			C *= get_tt1(dd[k]);
			//return C[m, k];
			return C;
		}
		double _D(int k)
		{
			double D = 1.0;

			D *= get_tt2(dd[k]);

			return D;
		}
		double _B(int u, int v) {
			//対称化
			double val = d1[u] * d1[v];
			val += d1[u] * d1[v];
			return val;
		}
	public:
		_memC() {

			M = 0;
			dd = 0;

			_nNode = 0;
			_dim = -1;

			gradN[0] = 0;
			gradN[1] = 0;
			gradN[2] = 0;
			gradG = 0;
			__grad = 0;
			d0 = 0;
			d1 = 0;
			d2 = 0;

			B = 0;
			tt0 = 0;

			hh0 = 0;

			tt1 = 0;

			hh1 = 0;

			tt2 = 0;

			hh2 = 0;


		}
		~_memC() {
			del();
		}
		void update2() {

			double* ptr;
			ptr = hh0;
			for (int k = 0; k < _dim; k++)
			{
				*ptr = __hh0(k);
				ptr++;
				//hh0[k] = __hh0(k);
			}

			ptr = tt0;
			for (int k = 0; k < _dim; k++) {
				*ptr = __tt0(k);
				ptr++;
				//tt0[k] = __tt0(k);
			}
			ptr = hh1;
			for (int k = 0; k < _dim; k++) {
				*ptr = __hh1(k);
				ptr++;
				//hh1[k] = __hh1(k);
			}

			for (int k = 0; k < _dim; k++) {
				tt1[k] = __tt1(k);
			}

			for (int k = 0; k < _dim; k++) {
				hh2[k] = __hh2(k);
			}

			for (int k = 0; k < _dim; k++) {
				tt2[k] = __tt2(k);
			}

			ptr = d0;
			for (int j = 0; j < _nNode; j++) {
				*ptr = _shape(j);
				ptr++;
			}
			_ref->w = 0;
			double W = 0;
			for (int j = 0; j < _nNode; j++) {
				W += d0[j] * _ref->buf_W[j];
			}
			_ref->w = W;

			for (int j = 0; j < _nNode; j++) {
				//d0[j]*=_ref->buf_W[j]/_ref->w;
			}

			ptr = d1;
			for (int j = 0; j < _nNode; j++) {
				//*ptr = _ref->buf_W[j] * _C(j) / _ref->w;
				*ptr = _C(j);
				ptr++;
			}

			ptr = d2;
			for (int j = 0; j < _nNode; j++) {
				//*ptr = _ref->buf_W[j] * _ref->w * _D(j) / _ref->w;
				*ptr = _D(j);
				ptr++;
			}
			ptr = B;
			for (int j = 0; j < _nNode; j++) {
				for (int jj = 0; jj < _nNode; jj++) {
					//*ptr = _ref->buf_W[j] * _ref->buf_W[jj] * _B(j, jj) / _ref->w / _ref->w;
					*ptr = _B(j, jj);
					ptr++;
				}
			}

			double X = 0, Y = 0, Z = 0;
			ptr = d0;
			double* ptr2;
			ptr2 = _ref->node;
			for (int i = 0; i < _nNode; i++)
			{
				X += *ptr * *ptr2;
				ptr2++;
				Y += *ptr * *ptr2;
				ptr2++;
				Z += *ptr * *ptr2;
				ptr2++;
				ptr++;
			}
			x = X;
			y = Y;
			z = Z;
			//covariant base vectors

			double fx = 0, fy = 0, fz = 0;
			ptr = d1;
			ptr2 = _ref->node;
			for (int i = 0; i < _nNode; i++)
			{
				fx += *ptr * *ptr2;
				ptr2++;
				fy += *ptr * *ptr2;
				ptr2++;
				//fz += *ptr * *ptr2;
				fz += *ptr * _ref->buf_z[i];
				ptr2++;
				ptr++;
			}
			gi[0] = fx;
			gi[1] = fy;
			gi[2] = fz;

			gij2[0] = gi[0] * gi[0] + gi[1] * gi[1];

			gi2[0] = fx;
			gi2[1] = fy;
			gi2[2] = 0;


			_dv = sqrt(_det1(gij2));
			_inv1(gij2, Gij2);
			double Fx = 0, Fy = 0, Fz = 0;
			Fx += gi2[0] * Gij2[0];
			Fy += gi2[1] * Gij2[0];
			Gi2[0] = Fx;
			Gi2[1] = Fy;
			Gi2[2] = 0;


			gij[0] = gi[0] * gi[0] + gi[1] * gi[1] + gi[2] * gi[2];


			_inv1(gij, Gij);
			dv = sqrt(_det1(gij));

			//contravatiant base vectors
			Fx = 0; Fy =0 ; Fz = 0;
			Fx += get_gt(0) * Gij[0];
			Fy += get_gt(1) * Gij[0];
			Fz += get_gt(2) * Gij[0];
			Gi[0] = Fx;
			Gi[1] = Fy;
			Gi[2] = Fz;


			fx = 0, fy = 0, fz = 0;
			ptr = d2;
			ptr2 = _ref->node;
			for (int i = 0; i < _nNode; i++)
			{
				fx += *ptr * *ptr2;
				ptr2++;
				fy += *ptr * *ptr2;
				ptr2++;
				//fz += *ptr * *ptr2;
				fz += *ptr * _ref->buf_z[i];
				ptr2++;
				ptr++;
			}
			bij[0] = fx;
			bij[1] = fy;
			bij[2] = fz;
			bij2[0] = fx;
			bij2[1] = fy;
			bij2[2] = 0;
			ptr = bij;
			Gammaijk[0] = bij[0] * get_Gt(0) + bij[1] * get_Gt(1) + bij[2] * get_Gt(2);
			Gammaijk2[0] = bij2[0] * get_Gt2(0) + bij2[1] * get_Gt2(1);
			double gx = get_gt(0);
			double gy = get_gt(1);
			double gz = get_gt(2);
			double upx = _ref->get_up(0);
			double upy = _ref->get_up(1);
			double upz = _ref->get_up(2);

			double Hx = gy * upz - gz * upy;
			double Hy = gz * upx - gx * upz;
			double Hz = gx * upy - gy * upx;

			//double Hx = gy;
			//double Hy = -gx;
			//double Hz = 0;
			double Nx = -(gy * Hz - gz * Hy);
			double Ny = -(gz * Hx - gx * Hz);
			double Nz = -(gx * Hy - gy * Hx);
			double  norm = sqrt(Nx * Nx + Ny * Ny + Nz * Nz);
			Nz = Nz / norm;
			Nx = Nx / norm;
			Ny = Ny / norm;
			norm = sqrt(Hx * Hx + Hy * Hy + Hz * Hz);
			Hz = Hz / norm;
			Hx = Hx / norm;
			Hy = Hy / norm;

			N[0] = Nx;
			N[1] = Ny;
			N[2] = Nz;
			H[0] = Hx;
			H[1] = Hy;
			H[2] = Hz;
		}
		void del() {
			if (__grad != 0)
			{
				delete[]__grad;
				__grad = 0;
			}
			if (M != 0) {
				for (int i = 0; i < _dim; i++) {
					delete[] M[i];
				}
				delete[] M;
				M = 0;
			}

			if (dd != 0)
			{
				delete[] dd;
				dd = 0;
			}

			if (d0 != 0) {
				delete[] d0;
				d0 = 0;
			}
			if (d1 != 0)
			{
				delete[] d1;
				d1 = 0;
			}

			if (d2 != 0)
			{
				delete[] d2;
				d2 = 0;
			}

			if (B != 0) {
				delete[] B;
				B = 0;
			}
			if (gradN[0] != 0) {
				delete[] gradN[0];
				delete[] gradN[1];
				delete[] gradN[2];
				gradN[0] = 0;
				gradN[1] = 0;
				gradN[2] = 0;

			}
			if (gradG != 0)delete[] gradG;
			if (tt0 != 0) {
				delete[] tt0;
				tt0 = 0;
				delete[] tt0;
				tt0 = 0;
				delete[] hh0;
				hh0 = 0;
				delete tt1;
				tt1 = 0;
				delete[] hh1;
				hh1 = 0;
				delete[] tt2;
				tt2 = 0;

				delete[] hh2;
				hh2 = 0;
			}
		}
		void update(int nNode, int Dim) {
			if (_nNode != nNode || _dim != Dim) {

				del();
				_nNode = nNode;
				_dim = Dim;

				d0 = new double[nNode];

				d1 = new double[nNode];

				d2 = new double[nNode];

				B = new double[nNode * nNode];

				tt0 = new double[Dim];

				hh0 = new double[Dim];

				tt1 = new double[Dim];

				hh1 = new double[Dim];

				tt2 = new double[Dim];

				hh2 = new double[Dim];

				gradN[0] = new double[3 * _nNode];
				gradN[1] = new double[3 * _nNode];
				gradN[2] = new double[3 * _nNode];
				gradG = new double[8 * _nNode];
				__grad = new double[_nNode];
				M = new double* [Dim];
				for (int i = 0; i < Dim; i++) {
					M[i] = new double[Dim];
				}
				dd = new int[nNode];
			}
		}
	private:
		void _inv1(double* From, double* to)
		{

			to[0] = 1 / From[0];

		}
		double _det1(double* m)
		{
			return m[0];
		}
	public:
		double __length()
		{
			return _dv;
		}
		void __length_u(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * (d1[s] * get_gt2(0));
				

				val = 0.5 * (_g11 * get_Gtt2()) * _dv;
				*ptr1 = val;
				ptr1++;
			}

		}
		void __length_v(double* ptr)
		{
			double* ptr1 = ptr;
			double val = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _g11 = 2 * (d1[s] * get_gt2(1));


				val = 0.5 * (_g11 * get_Gtt2()) * _dv;
				*ptr1 = val;
				ptr1++;
			}

		}
		double gammattt()
		{

			double val = 0;

			for (int s = 0; s < _ref->_nNode; s++)
			{
				val += d2[s] * _ref->node[s * 2 + 0] * get_gt2(0) +
					d2[s] * _ref->node[s * 2 + 1] * get_gt2(1);
			}
			double scale =1.0/(sqrt(get_gtt2()) *get_gtt2());
			return val*scale;
		}
		void gammattt_u(double* ptr)
		{
			double sx = 0;
			double sy = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				sx+= d2[s] * _ref->node[s * 2 + 0];
				sy += d2[s] * _ref->node[s * 2 + 1];
			}
			double* ptr1 = ptr;
			double scale = 1.0 / (sqrt(get_gtt2()) * get_gtt2());
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _scale = -1.5 /(sqrt(get_gtt2()* get_gtt2() *get_gtt2())) * (2*d1[s]*get_gt2(0));
				double val = d2[s] * get_gt2(0)*scale
					+ sx * d1[s]*scale;
				val += (sx * get_gt2(0) + sy * get_gt2(1)) * _scale;
				*ptr1 = val;
				ptr1++;
			}
		}
		void gammattt_v(double* ptr)
		{

			double sx = 0;
			double sy = 0;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				sx += d2[s] * _ref->node[s * 2 + 0];
				sy += d2[s] * _ref->node[s * 2 + 1];
			}
			double* ptr1 = ptr;
			double scale = 1.0 / (sqrt(get_gtt2()) * get_gtt2());
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _scale = -1.5 / (sqrt(get_gtt2() * get_gtt2() * get_gtt2())) * (2 * d1[s] * get_gt2(1));
				double val = d2[s] * get_gt2(1) * scale
					+ sy * d1[s] * scale;
				val += (sx * get_gt2(0)+sy*get_gt2(1)) * _scale;
				*ptr1 = val;
				ptr1++;
			}
		}
		double angle(_memC* other)
		{
			double vx = this->get_gt2(0);
			double vy = this->get_gt2(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt2(0);
			double wy = other->get_gt2(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;
			
			double val = (vx * wx + vy * wy+vz*wz) /(lv*lw);
			return val;
		}
		/*void angle_z1(_memC* other, double* ptr)
		{
			double vx = this->get_gt(0);
			double vy = this->get_gt(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt(0);
			double wy = other->get_gt(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;
			


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _vx = 0;
				double _vy = 0;
				double _vz = d1[s];
				double _gtt = 2 * d1[s] * this->get_gt2(2);
				double _gamma = 0.5 * _gtt * get_Gtt2() * sqrt(get_gtt2());
				double val = (_vx * wx + _vy * wy +_vz*wz) / (lv * lw);
				val += -(vx * wx + vy * wy+vz*wz) / (lv * lv * lw) * _gamma;
				*ptr1 = val;
				ptr1++;
			}
		}*/
		void angle_u1(_memC* other,double *ptr)
		{
			double vx = this->get_gt2(0);
			double vy = this->get_gt2(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt2(0);
			double wy = other->get_gt2(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;


			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _vx = d1[s];
				double _vy = 0;
				double _vz = 0;
				double _gtt = 2 * d1[s] * this->get_gt2(0);
				double _gamma = 0.5 * _gtt * get_Gtt2() * sqrt(get_gtt2());
				double val = (_vx * wx + _vy * wy + _vz * wz) / (lv * lw);
				val += -(vx * wx + vy * wy + vz * wz) / (lv * lv * lw) * _gamma;
				*ptr1 = val;
				ptr1++;
			}
		}
		void angle_v1(_memC* other, double* ptr)
		{
			double vx = this->get_gt2(0);
			double vy = this->get_gt2(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt2(0);
			double wy = other->get_gt2(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _vx = 0;
				double _vy = d1[s];
				double _vz = 0;
				double _gtt = 2 * d1[s] * this->get_gt2(1);
				double _gamma = 0.5 * _gtt * get_Gtt2() * sqrt(get_gtt2());
				double val = (_vx * wx + _vy * wy + _vz * wz) / (lv * lw);
				val += -(vx * wx + vy * wy + vz * wz) / (lv * lv * lw) * _gamma;
				*ptr1 = val;
				ptr1++;
			}
		}
		/*void angle_z2(_memC* other, double* ptr)
		{
			double vx = this->get_gt(0);
			double vy = this->get_gt(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt(0);
			double wy = other->get_gt(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _wx = 0;
				double _wy = 0;
				double _wz = other->d1[s];
				double _gtt = 2 * other->d1[s] * other->get_gt(2);
				double _gamma = 0.5 * _gtt * other->get_Gtt() * sqrt(other->get_gtt());
				double val = (vx * _wx + vy * _wy + vz * _wz) / (lv * lw);
				val += -(vx * wx + vy * wy + vz * wz) / (lv * lw * lw) * _gamma;
				*ptr1 = val;
				ptr1++;
			}
		}*/
		void angle_u2(_memC* other, double* ptr)
		{
			double vx = this->get_gt2(0);
			double vy = this->get_gt2(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt2(0);
			double wy = other->get_gt2(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _wx = other->d1[s];
				double _wy = 0;
				double _wz = 0;
				double _gtt = 2 * other->d1[s] * other->get_gt2(0);
				double _gamma = 0.5 * _gtt * other->get_Gtt2() * sqrt(other->get_gtt2());
				double val = (vx * _wx + vy * _wy + vz * _wz) / (lv * lw);
				val += -(vx * wx + vy * wy + vz * wz) / (lv * lw * lw) * _gamma;
				*ptr1 = val;
				ptr1++;
			}
		}


		void angle_v2(_memC* other, double* ptr)
		{
			double vx = this->get_gt2(0);
			double vy = this->get_gt2(1);
			double vz = 0;// this->get_gt2(2);
			double wx = other->get_gt2(0);
			double wy = other->get_gt2(1);
			double wz = 0;// other->get_gt2(2);
			double lv = _dv;
			double lw = other->_dv;

			double* ptr1 = ptr;
			for (int s = 0; s < _ref->_nNode; s++)
			{
				double _wx = 0;
				double _wy = other->d1[s];
				double _wz = 0;
				double _gtt = 2 * other->d1[s] * other->get_gt2(1);
				double _gamma = 0.5 * _gtt * other->get_Gtt2() * sqrt(other->get_gtt2());
				double val = (vx * _wx + vy * _wy + vz * _wz) / (lv * lw);
				val += -(vx * wx + vy * wy + vz * wz) / (lv * lw * lw) * _gamma;
				*ptr1 = val;
				ptr1++;
			}
		}
		

		//stress function
		double F(int i, int j) {
			int e = i * _nNode + j;
			return 0.5 * (_Sij[0] * B[e]) * _ref->refDv;
		}
		//stress function L2
		double F2(int i, int j) {
			return
				(_Sij[0] * (d2[i] - this->_ref->get__Gammattt() * d1[i]))
				* (_Sij[0] * (d2[j] - this->_ref->get__Gammattt() * d1[j])) * _ref->refDv;

		}
		//stress function L2 half
		double F4(int i) {
			return
				(_Sij[0] * (d2[i] - Gammaijk[0] * d1[i]));

		}
		double K(int i)
		{
			return(d2[i] - Gammaijk[0] * d1[i]);
		}
		void K(double*ptr,double sc)
		{
			double* ptr1 = ptr;
			double* ptr2 = d2;
			double* ptr3 = d1;

			for (int i = 0; i < _nNode; i++)
			{
				*ptr1 = (*ptr2 - Gammaijk[0] * *ptr3)*sc;
				ptr1++;
				ptr2++;
				ptr3++;
			}
			//return(d2[i] - Gammaijk[0] * d1[i]);
		}
		double K_z() {
			double val = 0;
			for (int I = 0; I < _nNode; I++) {
				val += (d2[I] - Gammaijk[0] * d1[I]) * _ref->buf_z[I];
			}
			return val;
		}
		//boundary term
		double MB2(int i)
		{
			return d1[i] * _Sij[0];
		}
		//error term
		double error() {
			double bij00 = this->get_btt(2) - this->get_Gammattt() * this->get_gt(2);
			
			return _ref->_Sij[0] * bij00;

		}
		//body force term
		double G(int i) {
			return d0[i] * _ref->refDv;
		}
		double G4(int i) {
			return d0[i] * _ref->_refDv;
		}
		//body force term accurate element area
		double G2(int i) {
			return d0[i] * dv;
		}
		double _d0(int i) {
			return d0[i];
		}
		double G3(int i) {
			return (_Sij[0] * (d2[i] - this->_ref->get__Gammattt() * d1[i])) * _ref->refDv;// / _ref->refDv / _ref->refDv;
		}
		//Hessian
		double __B(int i, int j) {
			double val = 0;

			val += 0.5 * d0[i] * d1[j] * this->get_gt(2) * this->dv * this->get_Gtt();
			val += 0.5 * d0[i] * d1[j] * this->get_gt(2) * this->dv * this->get_Gtt();
			return val;
		}
		//linear spring term
		double R(int i) {
			return -this->d0[i] * _ref->_z;
		}
		double R2(int i) {
			return -this->d0[i] * pow((_ref->_z - _ref->__z), 3) * _ref->_refDv;
		}
		double R3(int i) {
			return this->d0[i] * _ref->_refDv;
		}
		double MASS(int i, int j) {
			return this->d0[i] * this->d0[j] * _ref->_refDv;
		}
		double area() {
			return _ref->_refDv;
		}
		double A() {
			return this->dv;
		}
		double dA(int i) {
			double da = 0;
			da += 0.5 * get_Gtt() * (get_gt(2) * d1[i] + get_gt(2) * d1[i]) * dv;
			return da;
		}

		double fM()
		{
			//vector<double> ret;
			double S[1]; //covariant
			//double _S[0];//covariant-contravariant
			S[0] = 0;

			double val = 0;



			double A =  _ref->get__Gtt() * _ref->get__Gtt();


			double D = 0;
			//double E = 0;
			for (int s = 0; s < 3; s++)
			{
				D += _ref->get__gt(s) * (get_gt(s) - _ref->get__gt(s));
				D += _ref->get__gt(s) * (get_gt(s) - _ref->get__gt(s));
				//E += _ref->get__gi(k, s) * (get_gi(l, s) - _ref->get__gi(l, s));
				//E += _ref->get__gi(l, s) * (get_gi(k, s) - _ref->get__gi(k, s));
			}
			val += A * D;

			S[0] = val;

			return S[0];
		}
		double eM() {
			double membrane = 0;
			double A = 0.2 * _ref->get__Gtt() * _ref->get__Gtt() + _ref->get__Gtt() * _ref->get__Gtt();

			double D = 0;
			double E = 0;
			for (int s = 0; s < 3; s++)
			{
				D += _ref->get__gt(s) * (get_gt(s) - _ref->get__gt(s));
				D += _ref->get__gt(s) * (get_gt(s) - _ref->get__gt(s));
				E += _ref->get__gt(s) * (get_gt(s) - _ref->get__gt(s));
				E += _ref->get__gt(s) * (get_gt(s) - _ref->get__gt(s));
			}
			membrane += 0.25 * A * D * E * _ref->refDv;
			return membrane;
		}

		double eK() {
			double bending = 0;
			double __bij2 = get_btt(0) * N[0] + get_btt(1) * N[1] + get_btt(2) * N[2];
			double ___bij2 = _ref->get__btt(0) * N[0] + _ref->get__btt(1) * N[1] + _ref->get__btt(2) * N[2];
			double __bij = get_btt(0) * N[0] + get_btt(1) * N[1] + get_btt(2) * N[2];
			double ___bij = _ref->get__btt(0) * N[0] + _ref->get__btt(1) * N[1] + _ref->get__btt(2) * N[2];

			double A = 0.2 * _ref->get__Gtt() * _ref->get__Gtt() + 0.8 * _ref->get__Gtt() * _ref->get__Gtt();

			double D = 0;
			double E = 0;
			for (int s = 0; s < 3; s++)
			{
				D += N[s] * ((get_btt(s) - _ref->get__btt(s)) - _ref->get__Gammattt() * (get_gt(s) - _ref->get__gt(s)));
				E += N[s] * ((get_btt(s) - _ref->get__btt(s)) - _ref->get__Gammattt() * (get_gt(s) - _ref->get__gt(s)));
				//E += 0.5 * N[s] * ((tup.bij[l, k, s] - tup._bij[l, l, s]) - tup._Gammaijk[l, k, 0] * (tup.gi[0, s] - tup._gi[0, s]) - tup._Gammaijk[l, m, 1] * (tup.gi[1, s] - tup._gi[1, s]));
			}
			bending += A * D * E * _ref->refDv;

			return bending;
		}
		void computeGrads() {
			for (int i = 0; i < _nNode; i++) {
				for (int s = 0; s < 3; s++) {
					double valx = 0;
					double valy = 0;
					double valz = 0;
					valx += -d1[i] * (N[0]) * get_Gt(s);
					valy += -d1[i] * (N[1]) * get_Gt(s);
					valz += -d1[i] * (N[2]) * get_Gt(s);

					gradN[s][i * 3 + 0] = valx;
					gradN[s][i * 3 + 1] = valy;
					gradN[s][i * 3 + 2] = valz;
				}
			}
			//gradient of Gammaijk
			for (int A = 0; A < _nNode; A++) {
				double val = 0;
				val = d2[A] * get_Gt(2);
				val += get_btt(2) * get_Gtt() * d1[A];
				double AA = -get_Gtt() * (d1[A] * get_gt(2) + d1[A] * get_gt(2)) * get_Gtt();
				for (int oo = 0; oo < 3; oo++) {
					val += get_btt(oo) * get_gt(oo) * AA;
				}
				gradG[A] = val;
			}
		}
		//gradient of the bending matrix
		double L(int alpha)
		{

			double ddv = 0;

			ddv += 0.5 * get_Gtt() * (get_gt(2) * d1[alpha] + get_gt(2) * d1[alpha]);


			double coeff = dv;// * tup.area * tup.K;
			double membrane = 0;
			double bending = 0;




			double ik = 0;
			double jl = 0;
			double ij = 0;
			double kl = 0;
			ik += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			jl += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			ij += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			kl += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			double ijkl = 0.2 * (ij * get_Gtt() + kl * get_Gtt()) + ik * get_Gtt() + get_Gtt() * jl;

			double ff = get_Gtt();

			double ijkl2 = 0.2 * get_Gtt() * get_Gtt() + get_Gtt() * get_Gtt();

			for (int I = 0; I < _nNode; I++) {

				for (int J = 0; J < _nNode; J++) {

					//if(fx[I] ||fx[J])continue;
					double _E1 = d2[I] - get_Gammattt() * d1[I] - get_Gammattt() * d1[I];
					double _E2 = d2[J] - get_Gammattt() * d1[J] - get_Gammattt() * d1[J];


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


			ijkl = 0.2 * get_Gtt() * get_Gtt() + get_Gtt() * get_Gtt();

			for (int I = 0; I < _nNode; I++) {

				for (int J = 0; J < _nNode; J++) {

					//if(fx[I] ||fx[J])continue;
					double _E1 = (d2[J] - get_Gammattt() * d1[J] - get_Gammattt() * d1[J]);
					double _E2 = (d2[I] - get_Gammattt() * d1[I] - get_Gammattt() * d1[I]);
					double A1 = d2[I];
					double A2 = d2[J];



					double _B11 = 0, _B21 = 0;
					double _B12 = 0, _B22 = 0;

					_B11 += -gradG[alpha] * (d1[I]);

					_B21 += -get_Gammattt() * d1[I];
					_B12 += -gradG[alpha] * (d1[J]);
					_B22 += -get_Gammattt() * d1[J];


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

			return bending;
		}

		//gradient of the membrane stiffness matrix
		double _M(int alpha)
		{

			double ddv = 0;
			for (int n = 0; n < 2; n++) {
				for (int m = 0; m < 2; m++) {
					ddv += 0.5 * get_Gtt() * (get_gt(2) * d1[alpha] + get_gt(2) * d1[alpha]);
				}
			}


			double coeff = dv;// * tup.area * tup.K;
			double membrane = 0;
			double bending = 0;


			double ik = 0;
			double jl = 0;
			double ij = 0;
			double kl = 0;

			ik += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			jl += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			ij += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();
			kl += -get_Gtt() * (d1[alpha] * get_gt(2) + d1[alpha] * get_gt(2)) * get_Gtt();

			double ijkl = 0.2 * (ij * get_Gtt() + kl * get_Gtt()) + ik * get_Gtt() + get_Gtt() * jl;



			double ijkl2 = 0.2 * get_Gtt() * get_Gtt() + get_Gtt() * get_Gtt();

			for (int I = 0; I < _nNode; I++) {

				for (int J = 0; J < _nNode; J++) {

					for (int _i = 0; _i < 3; _i++) {
						for (int _j = 0; _j < 3; _j++) {

							double FF = 0.5 * (d1[I] * get_gt(_i) + d1[I] * get_gt(_i));
							double GG = 0.5 * (d1[J] * get_gt(_j) + d1[J] * get_gt(_j));

							double BB = _ref->def[I * 3 + _i] * _ref->def[J * 3 + _j] * FF * GG * coeff;  //hh term
							double CC = ijkl + ijkl2 * ddv;

							membrane += BB * CC;  //hh term
						}
					}
				}
			}


			ijkl = 0.2 * get_Gtt() * get_Gtt() + get_Gtt() * get_Gtt();

			for (int I = 0; I < _nNode; I++) {

				for (int J = 0; J < _nNode; J++) {

					//if(fx[I] ||fx[J])continue;

					double _F1 = d1[alpha] * d1[I] + d1[alpha] * d1[I];//=zero when x or y
					double _F2 = d1[alpha] * d1[J] + d1[alpha] * d1[J];//=zero when x or y
					for (int mm = 0; mm < 3; mm++) {
						for (int nn = 0; nn < 3; nn++) {
							double F1 = 0.25 * _F1 * (get_gt(mm) * d1[J] + get_gt(mm) * d1[J]);
							double F2 = 0.25 * _F2 * (get_gt(nn) * d1[I] + get_gt(nn) * d1[I]);
							if (nn == 2)membrane += _ref->def[I * 3 + 2] * _ref->def[J * 3 + mm] * ijkl * (F1)*coeff;
							if (mm == 2)membrane += _ref->def[I * 3 + nn] * _ref->def[J * 3 + 2] * ijkl * (F2)*coeff;
						}
					}

				}
			}

			return membrane;
		}

		//bending term strong axis
		double KN(int i, int k2, int j, int k)
		{
			double _val3 = 0;
			double _val4 = 0;

			double A = _ref->get__Gtt() * _ref->get__Gtt();
			double D = (d2[j] - Gammaijk[0] * d1[j]);
			double E = (d2[i] - Gammaijk[0] * d1[i]);
			_val3 += A * N[k] * N[k2] * (D) * (E);

			return _val3 * _ref->refDv;
		}
		void KN(_mySparse* M,  int64_t* _index, double sc)
		{
			for (int i = 0; i < _nNode; i++)
			{
				int I = _index[i] * 3;
				for (int s = 0; s < 3; s++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						int J = _index[j] * 3;
						for (int ss = 0; ss < 3; ss++)
						{

							double _val4 = KN(i, s, j, ss);
							M->_mat[0].coeffRef(I + s, J + ss) += _val4 * sc;
						}
					}
				}
			}

		}
		//bending term weak axis
		double KH(int i, int k2, int j, int k)
		{
			double _val3 = 0;
			double _val4 = 0;

			double A = _ref->get__Gtt() * _ref->get__Gtt();
			double D = (d2[j] - Gammaijk[0] * d1[j]);
			double E = (d2[i] - Gammaijk[0] * d1[i]);
			_val3 += A * H[k] * H[k2] * (D) * (E);

			return _val3 * _ref->refDv;
		}
		void KH(_mySparse* M, int64_t* _index, double sc)
		{
			for (int i = 0; i < _nNode; i++)
			{
				int I = _index[i] * 3;
				for (int s = 0; s < 3; s++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						int J = _index[j] * 3;
						for (int ss = 0; ss < 3; ss++)
						{

							double _val4 = KH(i, s, j, ss);
							M->_mat[0].coeffRef(I + s, J + ss) += _val4 * sc;
						}
					}
				}
			}

		}
		//bending boundary term
		double KB(int j, int k,int s)
		{
			double _val3 = 0;

			for (int g = 0; g < 2; g++)
			{
				for (int h = 0; h < 2; h++)
				{
					double A = (0.0 * _ref->get__Gtt() * _ref->get__Gtt() + _ref->get__Gtt() * _ref->get__Gtt());
					double D = (d2[j] - Gammaijk[0] * d1[j]);
					_val3 += A * N[k] * (D)*H[s];
					_val3 += A * H[k] * (D)*N[s];
				}
			}

			return _val3;// *_ref->refDv;
		}
		//angle term
		double T(int i, int s, int s2) {
			double val = 0;
			val += d1[i] * N[s] * _ref->get__Gtt() * H[s2];
			val += d1[i] * H[s] * _ref->get__Gtt() * N[s2];

			return val * _ref->refDv;
		}
		/*//bending boundary term
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
		//angle term
		*/



		//membrane term
		double _H(int i, int k2, int j, int k)
		{
			double _val4 = 0;
			double A = _ref->get__Gtt() * _ref->get__Gtt();

			double FF = (d1[j] * get_gt(k) + d1[j] * get_gt(k));
			double GG = (d1[i] * get_gt(k2) + d1[i] * get_gt(k2));
			_val4 += A * FF * GG;
			return _val4 * _ref->refDv * 0.25;
		}
		void _H(_mySparse* M, int64_t* _index, double sc)
		{
			for (int i = 0; i < _nNode; i++)
			{
				int64_t I = _index[i] * 3;
				for (int s = 0; s < 3; s++)
				{
					for (int j = 0; j < _nNode; j++)
					{
						long long J = _index[j] * 3;
						for (int ss = 0; ss < 3; ss++)
						{

							double _val4 = _H(i, s, j, ss);
							M->_mat[0].coeffRef(I + s, J + ss) += _val4 * sc;
						}
					}
				}
			}
		}
		void memory(_memC_ref* __mem) {
			if (__mem->__z < -1000) {
				__mem->__z = z;
			}
			__mem->_x = x;
			__mem->_y = y;
			__mem->_z = z;

			__mem->refDv = dv;
			__mem->_refDv = _dv;
			std::memcpy(__mem->_gi, gi2, sizeof(double) * 3);
			std::memcpy(__mem->_Gi, Gi2, sizeof(double) * 3);
			std::memcpy(__mem->_gij, gij2, sizeof(double) * 1);
			std::memcpy(__mem->_Gij, Gij2, sizeof(double) * 1);
			std::memcpy(__mem->_bij, bij2, sizeof(double) * 3);
			std::memcpy(__mem->_Gammaijk, Gammaijk2, sizeof(double) * 1);
			std::memcpy(__mem->_Sij, _Sij, sizeof(double) * 1);
		}
	};
	public ref class memC_ref {
	public:
		_memC_ref* __mem;
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

		memC_ref() {
			__mem = new _memC_ref();
		}
		void dispose() {
			delete __mem;
		}
		~memC_ref() {
			dispose();
		}
		!memC_ref() {
			dispose();
		}
		void set_z(double z)
		{
			__mem->set_z(z);
		}
		void update_z_phi(int nNode, KingOfMonsters::myDoubleArray ^ Z, KingOfMonsters::myDoubleArray^ phi) {
			if (Z != nullptr)
			{
				for (int i = 0; i < nNode; i++) {
					__mem->set_buf_z(i, (Z->_arr->__v)(i));
				}
			}
			if (phi != nullptr)
			{
				for (int i = 0; i < nNode; i++) {
					int e = i;
					__mem->set_buf_phi(i, (phi->_arr->__v)(e));
				}
			}
		}
		void update3(int nNode, KingOfMonsters::myDoubleArray^ node, KingOfMonsters::myDoubleArray^ weights, array<double>^ def,bool ignorez) {
			if (node != nullptr) {
				for (int i = 0; i < nNode; i++) {
					int e = i * 3;
					__mem->set_node(i, 0, (node->_arr->__v)(e + 0));
					__mem->set_node(i, 1, (node->_arr->__v)(e + 1));
					if (!ignorez) {
						__mem->set_node(i, 2, (node->_arr->__v)(e + 2));
					}
					else {
						__mem->set_node(i, 2, 0);
					}
					
				}
			}
			if (weights != nullptr)
			{
				for (int i = 0; i < nNode; i++) {
					__mem->set_buf_W(i, (weights->_arr->__v)(i));
				}
			}
			if (def != nullptr)
			{
				for (int i = 0; i < nNode; i++) {
					int e = i * 3;
					__mem->set_def(i, 0, def[e + 0]);
					__mem->set_def(i, 1, def[e + 1]);
					if (!ignorez) {
						__mem->set_def(i, 2, def[e + 2]);
					}
					else {
						__mem->set_def(i, 2, 0);
					}
				}
			}
		}
		void update3(int nNode, KingOfMonsters::myDoubleArray^ node, KingOfMonsters::myDoubleArray^ weights,array<double>^ def) {
			update3(nNode, node, weights, def, false);
		}
		void update(int nNode, int Dim) {
			__mem->update(nNode, Dim);
		}
		
	};
	public ref class memC {
	public:
		_memC* __mem=0;
	public:
		int _nNode;
		int dim;
		double x, y, z;
		double _x, _y, _z;
		double refDv;
		double _refDv;
		double dv;
	public:
		array<double> ^gi() {
			array<double> ^ret = gcnew array<double>(3);
			ret[0] = this->__mem->get_gt(0);
			ret[1] = this->__mem->get_gt(1);
			ret[2] = this->__mem->get_gt(2);
			return ret;
		}
		void setRef(memC_ref^ _mem_) {
			__mem->_ref = _mem_->__mem;
		}
		void computeGrads() {
			__mem->computeGrads();
		}
		double B(int i, int j) {
			return __mem->_B(i, j);
		}
		double G2(int i) {
			return __mem->G2(i);
		}
		double G3(int i) {
			return __mem->G3(i);
		}
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
		double area() {
			return __mem->area();
		}
		double fM()
		{
			return __mem->fM();
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
		double eM() {
			return __mem->eM();
		}
		double eK() {
			return __mem->eK();
		}
		double L(int a) {
			return __mem->L(a);
		}
		double M(int a) {
			return __mem->_M(a);
		}
		double d0(int i) {
			return __mem->_d0(i);
		}
		double K(int i) {
			return __mem->K(i);
		}
		void K(myDoubleArray ^arr,double sc,int shift) {
			__mem->K(arr->_arr->__v.data()+shift,sc);
		}
		double K_z() {
			return __mem->K_z();
		}
		double F(int i, int j) {
			return __mem->F(i, j);
		}
		double F2(int i, int j) {
			return __mem->F2(i, j);
		}
		double F4(int i) {
			return __mem->F4(i);
		}
		double error() {
			return __mem->error();
		}
		double G(int i) {
			return __mem->G(i);
		}
		double MB2(int i) {
			return __mem->MB2(i);
		}
		//bending around the strong axis
		double KN(int i, int s, int j, int k) {
			return __mem->KN(i, s, j, k);
		}
		void KN(mySparse^ M, myIntArray^ index, double sc)
		{
			__mem->KN(M->dat, index->data(), sc);
		}

		//bending around the weak axis
		double KH(int i, int s, int j, int k) {
			return __mem->KH(i, s, j, k);
		}
		void KH(mySparse^ M, myIntArray^ index, double sc)
		{
			__mem->KH(M->dat, index->data(), sc);
		}

		//bending boundary term -- to be corrected
		double KB(int j, int k,int s2) {
			return __mem->KB(j, k,s2);
		}
		double T(int i, int s, int s2) {
			return __mem->T(i, s, s2);
		}
		//axial force
		double H(int i, int s, int j, int k) {
			return __mem->_H(i, s, j, k);
		}
		void H(mySparse^ M, myIntArray^ index, double sc)
		{
			__mem->_H(M->dat, index->data(), sc);
		}
		void L(double val1) {
			__mem->set_L(val1);
		}
		double length()
		{
			return __mem->__length();
		}
		void length_u(mySparse^ mat, myIntArray^ index, int ii, double sc, double coeff)
		{
			this->__mem->__length_u(this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad , 0, sc, __mem->_ref->_nNode, true, coeff);
		}
		void length_v(mySparse^ mat, myIntArray^ index, int ii, double sc, double coeff)
		{
			this->__mem->__length_v(this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_ref->_nNode, false, coeff);
		}
		double gammattt()
		{
			return this->__mem->gammattt();
		}		
		void gammattt_u(mySparse^ mat, int ii, myIntArray^ index,  double sc, double coeff, Int64 shift)
		{
			this->__mem->gammattt_u( this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad , 0,sc, __mem->_ref->_nNode, true, coeff);
		}
		void gammattt_v(mySparse^ mat, int ii, myIntArray^ index,  double sc, double coeff, Int64 shift)
		{
			this->__mem->gammattt_v(this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad, 0, sc, __mem->_ref->_nNode, false, coeff);
		}
		double angle(memC^ other)
		{
			return this->__mem->angle(other->__mem);
		}
		/*void angle_z1(mySparse^ mat, memC^ other, myIntArray^ index, int ii, double sc, double coeff, Int64 shift)
		{
			this->__mem->angle_z1(other->__mem, this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_ref->_nNode, true, coeff);
		}*/
		void angle_u1(mySparse^ mat, memC^ other, myIntArray^ index, int ii, double sc, double coeff,Int64 shift)
		{
			this->__mem->angle_u1(other->__mem, this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad-shift, shift, sc, __mem->_ref->_nNode, true, coeff);
		}
		void angle_v1(mySparse^ mat, memC^ other, myIntArray^ index, int ii, double sc, double coeff, Int64 shift)
		{
			this->__mem->angle_v1(other->__mem, this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_ref->_nNode, false, coeff);
		}
		/*void angle_z2(mySparse^ mat, memC^ other, myIntArray^ index, int ii, double sc, double coeff, Int64 shift)
		{
			this->__mem->angle_z2(other->__mem, this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_ref->_nNode, false, coeff);
		}*/
		void angle_u2(mySparse^ mat, memC^ other, myIntArray^ index, int ii, double sc, double coeff, Int64 shift)
		{
			this->__mem->angle_u2(other->__mem, this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_ref->_nNode, false, coeff);
		}
		void angle_v2(mySparse^ mat, memC^ other, myIntArray^ index, int ii, double sc, double coeff, Int64 shift)
		{
			this->__mem->angle_v2(other->__mem, this->__mem->__grad);
			mat->dat->addrow(ii, index->_arr, __mem->__grad - shift, shift, sc, __mem->_ref->_nNode, false, coeff);
		}
		void compute() {
			
			__mem->update2();
			x = __mem->x;
			y = __mem->y;
			z = __mem->z;
			dv = __mem->dv;
			_refDv = __mem->_dv;
			refDv = __mem->dv;
		}
		void update_lo(double lo) {
			__mem->set_lo(lo);
		}
		void update_elem(int nNode, int Dim,array<double, 2>^ M, array<int>^ dd) {
			__mem->update(nNode, Dim);

			int i = 0;
			for (int j = 0; j < Dim; j++) {
				for (int k = 0; k < Dim; k++) {
					__mem->set_M(j, k, M[j, k]);
				}

				for (int i = 0; i < nNode; i++) {
					__mem->set_dd(i, dd[i]);
				}
			}
		}
		void memory(memC_ref^ _mem_) {
			__mem->memory(_mem_->__mem);
			_x = __mem->_ref->_x;
			_y = __mem->_ref->_y;
			_z = __mem->_ref->_z;
		}
		memC() {
			__mem = new _memC();
		}
		void dispose() {
			if(__mem!=0)
			delete __mem;
			__mem = 0;
		}
		~memC() {
			dispose();
		}
		!memC() {
			dispose();
		}
	};

}
