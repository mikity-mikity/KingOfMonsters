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
#define EIGEN_USE_LAPACK
#define EIGEN_USE_MKL_ALL
int64_t previdentiyN = 0;
//std::vector<cudaStream_t> streams;
Eigen::MatrixXd I;
#ifndef _CPU
#define STREAMCOUNT 16
bool __cuinit = false;
#endif

//static std::vector<Eigen::SparseMatrix<double>> e;
//std::vector<std::vector<cusparseHandle_t>> sp_handle;
std::map<std::tuple<KingOfMonsters::_mySparse*,int64_t>, KingOfMonsters::spgemm> dict;
std::map<KingOfMonsters::_mySparse*, Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> dict2;
std::map< Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>*, std::vector<int64_t>>  map;

void KingOfMonsters::cuda::disable()
{
#ifndef _CPU
	__cuinit = false;
#endif
}
double KingOfMonsters::_helper::VarPro(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd *_r1, Eigen::VectorXd* _r2, double dt, int tt)
{
	Eigen::VectorXd X0 = *zz;
	Eigen::VectorXd x0 = *phi;

	int _mt = 0;
	int n = __U->cols();
	int m = _mats1.size();
	_mt = omp_get_max_threads();
	int _mt2 = 0;
#pragma omp parallel
	{
#pragma omp single
		_mt2 = omp_get_num_threads();
	}
	if (_mt2 > _mt)_mt = _mt2;
	Eigen::setNbThreads(_mt);
	//int m = __W->rows();

	int nr1 = _r1->size();
	int nr2 = _r2->size();

	Eigen::VectorXd r1 = /*__W->transpose() * */*_r1;
	Eigen::VectorXd r2 = /*__W->transpose() * */*_r2;

	Eigen::VectorXd b1(m), b2(m);
	b1.setZero();
	b2.setZero();

	Eigen::VectorXd ___r1(m), ___r2(m);

	___r1 = r1; ___r2 = r2;
	Eigen::MatrixXd Jx(m, n);
	Eigen::MatrixXd JX(m, n);
	Jx.setZero();
	JX.setZero();
	Eigen::VectorXd rhsx(m);
	Eigen::VectorXd rhsX(m);
	Eigen::MatrixXd __I(n, n);
	__I.setIdentity();
	double residual = 0;
	double residual2 = 0;
	//if (tt % 2 == 0)
	//{
#pragma omp parallel for
		for (int ii = 0; ii < _mt; ii++)
		{
			int S = ii * m / _mt;
			int E = (ii + 1) * m / _mt;
			for (int i = S; i < E; i++)
			{

				auto M1 = _mats1[i];
				auto M2 = _mats2[i];
				auto M3 = _mats3[i];//sparse
				if (M1->rows() == n && M1->cols() == n)
				{
					Jx.row(i) = (X0.transpose() * *M1);
					JX.row(i) = x0.transpose() * *M2;
					b1(i) = X0.transpose() * *M1 * x0;
					b2(i) = x0.transpose() * *M2 * X0;
					double a = b1(i);
					double b = b2(i);
					double c = r1(i);
					double d = r2(i);

					residual += b1(i);
					residual2 +=___r1(i);
				}
				else if (M1->rows() == 1 && M1->cols() == n)
				{
					Jx.row(i) = *M1;
					b1(i) = (*M1 * x0)(0, 0);
				}
				else  if (M2->rows() == 1 && M2->cols() == n) {
					JX.row(i) = *M2;
					b2(i) = (*M2 * X0)(0, 0);
				}
			}
		}

		rhsx = (b1 - ___r1).transpose() * coeff->asDiagonal() * Jx;
		rhsX = (b2 - ___r2).transpose() * coeff->asDiagonal() * JX;
		double norm1 = rhsx.norm();
		double norm2 = rhsX.norm();

		Eigen::FullPivLU<Eigen::MatrixXd> qr;
		Eigen::MatrixXd ee = Jx.transpose() * coeff->asDiagonal() * Jx;
		qr.compute(ee);
		Eigen::VectorXd dx = -qr.solve(rhsx);
		x0 += dx * 1;
	//}
	//else {

		Jx.setZero();
		JX.setZero();
		rhsx.setZero();
		rhsX.setZero();
		b1.setZero();
		b2.setZero();
#pragma omp parallel for
		for (int ii = 0; ii < _mt; ii++)
		{
			int S = ii * m / _mt;
			int E = (ii + 1) * m / _mt;
			for (int i = S; i < E; i++)
			{

				auto M1 = _mats1[i];
				auto M2 = _mats2[i];
				auto M3 = _mats3[i];//sparse
				if (M1->rows() == n && M1->cols() == n)
				{
					Jx.row(i) = (X0.transpose() * *M1);
					JX.row(i) = x0.transpose() * *M2;
					b1(i) = X0.transpose() * *M1 * x0;
					b2(i) = x0.transpose() * *M2 * X0;
				}
				else if (M1->rows() == 1 && M1->cols() == n)
				{
					Jx.row(i) = *M1;
					b1(i) = (*M1 * x0)(0, 0);
				}
				else  if (M2->rows() == 1 && M2->cols() == n) {
					JX.row(i) = *M2;
					b2(i) = (*M2 * X0)(0, 0);
				}
			}
		}
		//b2 += (*_mats3[n + 1] * X0);
		//b1 += (*_mats3[n] * x0);

		//JX += *_mats2[n];
		//Jx += *_mats1[n];

		rhsx = (b1 - ___r1).transpose() * coeff->asDiagonal() * Jx;
		double norm = rhsx.norm();
		rhsX = (b2 - ___r2).transpose() * coeff->asDiagonal() * JX;
		ee = Jx.transpose() * coeff->asDiagonal() * Jx;
		
		qr.compute(ee);
		//eei = ee.inverse();
		//Eigen::MatrixXd ffff=I-(ee* eei); 
		//double ffffff = ffff.norm();
		Eigen::MatrixXd ff = JX.transpose() * coeff->asDiagonal() * Jx * qr.inverse();
		Eigen::MatrixXd gg = JX.transpose() * coeff->asDiagonal() * JX;
		//ff.setZero();
		Eigen::MatrixXd kk = ff * Jx.transpose() * coeff->asDiagonal() * JX;
		Eigen::MatrixXd mm = gg - kk;

		qr.compute(mm);

		Eigen::VectorXd ll = rhsX - ff * rhsx;
		Eigen::VectorXd dX = -qr.solve(ll);
		X0 += dX * dt;
		//}
	//}
	*zz = X0;
	*phi = x0;
	return (residual-residual2)/residual;
}

double KingOfMonsters::_helper::ALT(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt)
{
	Eigen::VectorXd X0 = *zz;
	Eigen::VectorXd x0 =  *phi;

	int _mt = 0;
	int n = __U->cols();
	int m = _mats1.size();
	_mt = omp_get_max_threads();
	int _mt2 = 0;
#pragma omp parallel
	{
#pragma omp single
		_mt2 = omp_get_num_threads();
	}
	if (_mt2 > _mt)_mt = _mt2;
	Eigen::setNbThreads(_mt);
	//int m = __W->rows();

	int nr1 = _r1->size();
	int nr2 = _r2->size();

	Eigen::VectorXd r1 = /*__W->transpose() * */*_r1;
	Eigen::VectorXd r2 = /*__W->transpose() * */*_r2;

	Eigen::VectorXd b1(m), b2(m);
	b1.setZero();
	b2.setZero();

	Eigen::VectorXd ___r1(m), ___r2(m);

	___r1 = r1; ___r2 = r2;
	Eigen::MatrixXd Jx(m, n);
	Eigen::MatrixXd JX(m, n);
	Jx.setZero();
	JX.setZero();
	Eigen::VectorXd rhsx(m);
	Eigen::VectorXd rhsX(m);

#pragma omp parallel for
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * m / _mt;
		int E = (ii + 1) * m / _mt;
		for (int i = S; i < E; i++)
		{

			auto M1 = _mats1[i];
			auto M2 = _mats2[i];
			auto M3 = _mats3[i];//sparse
			if (M1->rows() == n && M1->cols() == n)
			{
				Jx.row(i) = (*M2* X0).transpose();
				JX.row(i) = x0.transpose() * *M2;
				b1(i) = x0.transpose() * *M2 * X0;
				b2(i) = x0.transpose() * *M2 * X0;
			}
			else if (M1->rows() == 1 && M1->cols() == n)
			{
				Jx.row(i) = *M1;
				b1(i) = ( *M1 *x0)(0, 0);
			}
			else  if (M2->rows() == 1 && M2->cols() == n) {
				JX.row(i) = *M2;
				b2(i) = (*M2 * X0)(0, 0);
			}
		}
	}
	//b2 += (*_mats3[n + 1] * X0);
	//b1 += (*_mats3[n] * x0);

	//JX += *_mats2[n];
	//Jx += *_mats1[n];

	rhsx = (b1 - ___r1).transpose() * coeff->asDiagonal() * Jx;
	rhsX = (b2 - ___r2).transpose() * coeff->asDiagonal() * JX;

	//Console::WriteLine("Jx=" + Jx.sum().ToString());
	//Console::WriteLine("JX=" + JX.sum().ToString());
	//Console::WriteLine("residual normx=" + rhsx.norm().ToString());
	//Console::WriteLine("residual normX=" + rhsX.norm().ToString());
	Eigen::FullPivLU<Eigen::MatrixXd> qr;
	Eigen::MatrixXd I(n, n);
	I.setIdentity();
	Eigen::MatrixXd ee = Jx.transpose() * coeff->asDiagonal() * Jx;
	qr.compute(ee);
	//Eigen::MatrixXd eei = ee.inverse();
	Eigen::VectorXd dx = -qr.solve(rhsx);
	//Eigen::VectorXd dX = -(JX.transpose() * JX + I * 0.000000000001).inverse() * rhsX;
	//X0 += dX * dt;
	x0 += dx * dt;


	Jx.setZero();
	JX.setZero();
	rhsx.setZero();
	rhsX.setZero();
	b1.setZero();
	b2.setZero();
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * m / _mt;
		int E = (ii + 1) * m / _mt;
		for (int i = S; i < E; i++)
		{

			auto M1 = _mats1[i];
			auto M2 = _mats2[i];
			auto M3 = _mats3[i];//sparse
			if (M1->rows() == n && M1->cols() == n)
			{
				Jx.row(i) = (*M2 * X0).transpose();
				JX.row(i) = x0.transpose() * *M2;

				b1(i) = x0.transpose() * *M2 * X0;
				b2(i) = x0.transpose() * *M2 * X0;
			}
			else if (M1->rows() == 1 && M1->cols() == n)
			{
				Jx.row(i) = *M1;
				b1(i) = (*M1 *x0)(0, 0);
			}
			else  if (M2->rows() == 1 && M2->cols() == n) {
				JX.row(i) = *M2;
				b2(i) = (*M2 * X0)(0, 0);
			}
		}
	}

	rhsx = (b1 - ___r1).transpose() * coeff->asDiagonal() * Jx;
	double norm = rhsx.norm();
	rhsX = (b2 - ___r2).transpose() * coeff->asDiagonal() * JX;
	qr.compute(JX.transpose() * coeff->asDiagonal() * JX);

	Eigen::VectorXd dX = -qr.solve(rhsX);
	X0 += dX * dt;
	//}
	*zz =X0;
	*phi = x0;
	return rhsX.norm();
}
double KingOfMonsters::_helper::Simple(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt)
{
	Eigen::VectorXd X0 = *zz;

	int _mt = 0;
	int n = X0.size();
	int m = _mats1.size();
	_mt = omp_get_max_threads();
	int _mt2 = 0;
#pragma omp parallel
	{
#pragma omp single
		_mt2 = omp_get_num_threads();
	}
	if (_mt2 > _mt)_mt = _mt2;
	Eigen::setNbThreads(_mt);
	//int m = __W->rows();

	int nr1 = _r2->size();

	Eigen::VectorXd r2 = /*__W->transpose() * */*_r2;

	Eigen::VectorXd b2(m);
	b2.setZero();

	Eigen::VectorXd ___r2(m);

	___r2 = r2;
	Eigen::MatrixXd JX(m, n);
	JX.setZero();
	Eigen::VectorXd rhsX(m);

#pragma omp parallel for
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * m / _mt;
		int E = (ii + 1) * m / _mt;
		for (int i = S; i < E; i++)
		{

			auto M1 = _mats1[i];
			auto M2 = _mats2[i];
			auto M3 = _mats3[i];//sparse
			if (M1->rows() == n && M1->cols() == n)
			{
				JX.row(i) = X0.transpose() * *M2;
				b2(i) = X0.transpose() * *M2 * X0;
			}
			if (M2->rows() == 1 && M2->cols() == n) {
				JX.row(i) = *M2;
				b2(i) = (*M2 * X0)(0, 0);
			}
		}
	}
	//b2 += (*_mats3[n + 1] * X0);
	//b1 += (*_mats3[n] * x0);

	//JX += *_mats2[n];
	//Jx += *_mats1[n];

	rhsX = (b2 - ___r2).transpose() * coeff->asDiagonal() * JX;

	//Console::WriteLine("Jx=" + Jx.sum().ToString());
	//Console::WriteLine("JX=" + JX.sum().ToString());
	//Console::WriteLine("residual normx=" + rhsx.norm().ToString());
	//Console::WriteLine("residual normX=" + rhsX.norm().ToString());
	Eigen::FullPivLU<Eigen::MatrixXd> qr;
	Eigen::MatrixXd I(n, n);
	I.setIdentity();
	Eigen::MatrixXd ee = JX.transpose() * coeff->asDiagonal() * JX;
	qr.compute(ee+I*0.000001);
	//Eigen::MatrixXd eei = ee.inverse();
	Eigen::VectorXd dX = -qr.solve(rhsX);
	//Eigen::VectorXd dX = -(JX.transpose() * JX + I * 0.000000000001).inverse() * rhsX;
	//X0 += dX * dt;
	X0 += dX * dt;


	*zz = X0;
	*phi = X0;
	return rhsX.norm();
}

double KingOfMonsters::_helper::GN(Eigen::VectorXd* coeff, Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt,KingOfMonsters::cuda *cuda)
{
	Eigen::VectorXd X0 = *zz;
	Eigen::VectorXd x0 = *phi;

	int _mt = 0;
	int n = __U->cols();
	int m = _mats1.size();
	_mt = omp_get_max_threads();
	int _mt2 = 0;
#pragma omp parallel
	{
#pragma omp single
		_mt2 = omp_get_num_threads();
	}
	if (_mt2 > _mt)_mt = _mt2;
	Eigen::setNbThreads(_mt);
	//int m = __W->rows();

	int nr1 = _r1->size();
	int nr2 = _r2->size();

	Eigen::VectorXd r1 = /*__W->transpose() * */*_r1;
	Eigen::VectorXd r2 = /*__W->transpose() * */*_r2;

	Eigen::VectorXd b1(m), b2(m);
	b1.setZero();
	b2.setZero();

	Eigen::VectorXd ___r1(m), ___r2(m);

	___r1 = r1; ___r2 = r2;
	Eigen::MatrixXd Jx(m, n);
	Eigen::MatrixXd JX(m, n);
	Jx.setZero();
	JX.setZero();
	Eigen::VectorXd rhsx(m);
	Eigen::VectorXd rhsX(m);

#pragma omp parallel for
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * m / _mt;
		int E = (ii + 1) * m / _mt;
		for (int i = S; i < E; i++)
		{

			auto M1 = _mats1[i];
			auto M2 = _mats2[i];
			auto M3 = _mats3[i];//sparse
			if (M1->rows() == n && M1->cols() == n)
			{
				Jx.row(i) = (*M2 * X0).transpose();
				JX.row(i) = x0.transpose() * *M2;
				b1(i) = x0.transpose() * *M2 * X0;
				b2(i) = x0.transpose() * *M2 * X0;
			}
			else if (M1->rows() == 1 && M1->cols() == n)
			{
				Jx.row(i) = *M1;
				b1(i) = (*M1 * x0)(0, 0);
			}
			else  if (M2->rows() == 1 && M2->cols() == n) {
				JX.row(i) = *M2;
				b2(i) = (*M2 * X0)(0, 0);
			}
		}
	}
	//b2 += (*_mats3[n + 1] * X0);
	//b1 += (*_mats3[n] * x0);

	//JX += *_mats2[n];
	//Jx += *_mats1[n];

	rhsx = (b1 - ___r1).transpose() * coeff->asDiagonal() * Jx;
	rhsX = (b2 - ___r2).transpose() * coeff->asDiagonal() * JX;

	//Eigen::FullPivLU<Eigen::MatrixXd> qr;
	Eigen::MatrixXd I(2*n, 2*n);
	I.setIdentity();
	
	Eigen::MatrixXd SYS(n * 2, n * 2);
	SYS.topLeftCorner(n, n) = Jx.transpose() * Jx;
	SYS.bottomRightCorner(n, n) = JX.transpose() * JX;
	SYS.topRightCorner(n, n) = JX.transpose() * Jx;
	SYS.bottomLeftCorner(n, n) = Jx.transpose() * JX;
	Eigen::VectorXd rhs(2 * n);
	rhs.topRows(n) = rhsx;
	rhs.topRows(n) = rhsX;
	SYS = SYS + I * 0.000000001;
	//qr.compute(SYS+I*0.0000001);
	//Eigen::MatrixXd eei = ee.inverse();
	//Eigen::VectorXd dxX = -qr.solve(rhs);
	Eigen::VectorXd ret;
	KingOfMonsters::_mySparse* mm = new KingOfMonsters::_mySparse();
	mm->_dmat = SYS;
	mm->_solve0_gpu(cuda,& rhs, &ret, cuda->fastest());
	//auto dx = dxX.topRows(n);
	//auto dX = dxX.bottomRows(n);
	auto dx = ret.topRows(n);
	auto dX = ret.bottomRows(n);

	x0 += dx * 1;
	X0 += dX * 1;


	*zz = X0;
	*phi = x0;
	delete mm;
	return rhsX.norm();
}

void KingOfMonsters::_helper::write(Eigen::VectorXd* coeff, Eigen::VectorXd* phi0, Eigen::VectorXd* zz0,Eigen::VectorXd* phi, Eigen::VectorXd* zz, Eigen::MatrixXd* __U, Eigen::MatrixXd* __V, Eigen::MatrixXd* __W, std::vector<Eigen::SparseMatrix<double>*> _mats1, std::vector<Eigen::SparseMatrix<double>*> _mats2, std::vector<Eigen::SparseMatrix<double>*> _mats3, Eigen::VectorXd* _r1, Eigen::VectorXd* _r2, double dt, int tt)
{
	Eigen::VectorXd X0 = *zz;
	Eigen::VectorXd x0 = *phi;

	int _mt = 1;
	int n = __U->cols();
	int m = _mats1.size();

	//int m = __W->rows();

	int nr1 = _r1->size();
	int nr2 = _r2->size();

	Eigen::VectorXd r1 = *_r1;
	Eigen::VectorXd r2 = *_r2;

	Eigen::VectorXd b1(m), b2(m);
	b1.setZero();
	b2.setZero();

	Eigen::VectorXd ___r1(m), ___r2(m);

	___r1 = r1; ___r2 = r2;
	Eigen::MatrixXd Jx(m, n);
	Eigen::MatrixXd JX(m, n);
	Jx.setZero();
	JX.setZero();
	Eigen::VectorXd rhsx(m);
	Eigen::VectorXd rhsX(m);
	std::ofstream Aijk;  //i:x, j:X
	std::ofstream Bjm;   //j:X
	std::ofstream Din;   //i:x
	
	std::ofstream _x0;
	std::ofstream _X0;
	std::ofstream xsol;
	std::ofstream Xsol;

	std::ofstream rhok;  //i:x, j:X
	std::ofstream cm;   //j:X
	std::ofstream en;   //i:x

	Aijk.open("example01_Aijk.txt", std::ios::out);
	Bjm.open("example01_Bjm.txt", std::ios::out);
	Din.open("example01_Din.txt", std::ios::out);
	rhok.open("example01_rhok.txt", std::ios::out);
    cm.open("example01_cm.txt", std::ios::out);
	en.open("example01_en.txt", std::ios::out);
	_x0.open("example01_phi0.txt", std::ios::out);
	_X0.open("example01_z0.txt", std::ios::out);
	xsol.open("example01_phisol.txt", std::ios::out);
	Xsol.open("example01_zsol.txt", std::ios::out);

	int kk = 0;
	int nn = 0;
	int mm = 0;
	for (int ii = 0; ii < _mt; ii++)
	{
		int S = ii * m / _mt;
		int E = (ii + 1) * m / _mt;
		for (int i = S; i < E; i++)
		{

			auto M1 = _mats1[i];
			auto M2 = _mats2[i];
			auto M3 = _mats3[i];//sparse
			if (M1->rows() == n && M1->cols() == n)
			{
				Jx.row(i) = (*M2 * X0).transpose();
				JX.row(i) = x0.transpose() * *M2;
				for (int k = 0; k < M2->outerSize(); ++k) {
					for (Eigen::SparseMatrix<double>::InnerIterator it(*M2, k); it; ++it) {
						Aijk << it.row() + 1 << " , " << it.col() + 1 << " , " << kk + 1 << " , " << it.value() << std::endl;
					}
				}
				Aijk << n << " , " << n << " , " << kk+1 << " , " << 0 << std::endl;
				b1(i) = x0.transpose() * *M2 * X0;
				b2(i) = x0.transpose() * *M2 * X0;
				rhok << kk+1 << " , " << ___r1(i) << std::endl;

				kk++;

			}
			else if (M1->rows() == 1 && M1->cols() == n)
			{
				//Cin
				Jx.row(i) = *M1;
				for (int k = 0; k < M1->outerSize(); ++k) {
					for (Eigen::SparseMatrix<double>::InnerIterator it(*M1, k); it; ++it) {
						Din << it.col() + 1 << " , " << nn + 1 << " , " << it.value() << std::endl;
					}
				}
				Din  << n << " , " << nn+1 << " , " << 0 << std::endl;
				b1(i) = (*M1 * x0)(0, 0);
				en << nn + 1 << " , " << ___r1(i) << std::endl;
				nn++;
			}
			else  if (M2->rows() == 1 && M2->cols() == n) {
				//Bjm
				JX.row(i) = *M2;
				for (int k = 0; k < M2->outerSize(); ++k) {
					for (Eigen::SparseMatrix<double>::InnerIterator it(*M2, k); it; ++it) {
						Bjm << it.col() + 1 << " , " << mm + 1 << " , " << it.value() << std::endl;
					}
				}
				Bjm << n << " , " << mm+1 << " , " << 0 << std::endl;
				b2(i) = (*M2 * X0)(0, 0);
				cm << mm + 1 << " , " << ___r2(i) << std::endl;
				mm++;
			}
		}
	}

	for (int i = 0; i < n; i++)
	{
		_x0 << i + 1 << " , " << (*phi0)(i) << std::endl;
		_X0 << i + 1 << " , " << (*zz0)(i) << std::endl;
		xsol << i + 1 << " , " << (*phi)(i) << std::endl;
		Xsol << i + 1 << " , " << (*zz)(i) << std::endl;
	}
	Aijk.close();
	Bjm.close();
	Din.close();
	rhok.close();
	cm.close();
	en.close();
	_x0.close();
	_X0.close();
	xsol.close();
	Xsol.close();
	
}

KingOfMonsters::cuda::cuda(int64_t N) {
#ifndef _CPU
	dict.clear();
	dict2.clear();
	map.clear();
	Eigen::initParallel();
	//Eigen::setNbThreads(omp_get_max_threads());
	//e.shrink_to_fit();
	//e.clear();
	I.resize(0, 0);
	omp_set_dynamic(false);
	//omp_set_num_threads(16);
	//prevT_A = 0;
	//prevN = 0;
	//prevwn = 0;
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
	//solver_handleSp.resize(_count);
	//cusparse_handle.resize(_count);
	_streams.resize(_count);
	for (int ii = 0; ii < _count; ii++)
	{
		cudaSetDevice(ii);
		solver_handle[ii].resize(STREAMCOUNT);
		//cusparse_handle[ii].resize(STREAMCOUNT);
		//solver_handleSp[ii].resize(STREAMCOUNT);
		_streams[ii].resize(STREAMCOUNT);
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

	//for (int i = 0; i < MAXDEVICE; i++)
	//{
	//	cublas_handle[i] = 0;
	//}
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		for (int j = 0; j < STREAMCOUNT; j++)
		{
			//cusolverSpCreate(&solver_handleSp[ii][j]);
			//cusparseCreate(&cusparse_handle[ii][j]);
		}
	}
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		cusolverStatus_t status;
		for (int j = 0; j < STREAMCOUNT; j++)
			status = cusolverDnCreate(&solver_handle[ii][j]);

		//auto status2 = cublasCreate(&cublas_handle[ii]);

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
		//if (status2 == cublasStatus_t::CUBLAS_STATUS_SUCCESS)
		{
			initialized = true;
			failed = false;
		}
		//else {
		//	initialized = false;
		//	failed = true;
			//cublas_handle[ii] = 0;
		//	return;
		//}
	}
	//cusparseCreate(&sp_handle);
	if (!initialized || failed)return;


	//__mgM2 = 0;
	//__mgrhs2 = 0;

	for (int i = 0; i < MAXDEVICE; i++)
	{
		__mgM[i] = 0;
		__mgrhs[i] = 0;
		__mgC[i] = 0;
		__work[i] = 0;
		__work2[i] = 0;
		work_size[i] = 0;
		work_size2[i] = 0;
		__info[i] = 0;
	}
	_fastest = 0;

	int64_t _N = 1000;
	for (int ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);		
		auto err = cudaMalloc(&__mgM[ii], sizeof(double) * _N * _N);
		err = cudaMalloc(&__mgrhs[ii], sizeof(double) * _N * (5));
		//cudaMalloc(&__mgC[ii], sizeof(double) * _N * (10 + (_N / count())));
		cudaMalloc(&__info[ii], sizeof(int64_t) * 10);
	}

	//cudaMallocHost(&__mgM2, sizeof(double) * _N * _N);
	//cudaMallocHost(&__mgrhs2, sizeof(double) * _N * _N);

	//for (int64_t ii = 0; ii < count(); ii++)
	{
		//cudaSetDevice(ii);
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
			for (int64_t i = 0; i < _N; i++)
			{
				for (int64_t j = i - 5; j < i + 5; j++)
				{
					if (j >= 0 && j < _N)
						m.adddat(i, j, i);
				}
			}
			Eigen::VectorXd rhs(_N);
			for (int64_t i = 0; i < _N; i++)rhs[i] = i;
			m.ofDat();
			m.clearcoeff();
			m._ofAtA(&m);
			Eigen::VectorXd ret(_N);
			m._solveLU_gpu(this, &rhs, &ret, ii);
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<microseconds>(stop - start);
			speed[ii] = duration.count();
		}
	}
	_fastest=std::distance(speed.begin(), std::min_element(speed.begin(), speed.end()));
	//_fastest = 0;
	for (int64_t ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		cudaFree(__mgM[ii]);
		cudaFree(__mgrhs[ii]);
		cudaFree(__info[ii]);
		//cudaFree(__mgC[ii]);
	}
	//cudaFreeHost(__mgM2);
	//cudaFreeHost(__mgrhs2);
	for (int64_t ii = 0; ii < count(); ii++)
	{
		cudaSetDevice(ii);
		
		cudaMalloc(&__mgM[ii], sizeof(double) * N * N);
		
		cudaMalloc(&__mgrhs[ii], sizeof(double) * N*(2));
		cudaMalloc(&__mgC[ii], sizeof(double)* N*(10));// *(1 + (N) / count()));
		cudaMalloc(&__info[ii], sizeof(int64_t) * 10);
	}




#else
initialized = true;
failed = false;
return;
#endif
	
}


int* KingOfMonsters::cuda::devicelist()
{
	return _deviceList;
}
/*
double** KingOfMonsters::cuda::array_d_A()
{
	return this->_array_d_A;
}
double** KingOfMonsters::cuda::array_d_B()
{
	return this->_array_d_B;
}
double** KingOfMonsters::cuda::array_d_work()
{
	return this->_array_d_work;
}
bool KingOfMonsters::cuda::canpeeraccess(int64_t i, int64_t j)
{
	int64_t canAccessPeer = 0;
	cudaDeviceCanAccessPeer(&canAccessPeer, i, j);
	return canAccessPeer;
}*/
/*double* KingOfMonsters::cuda::L()
{
	return this->_L;
}*/
bool KingOfMonsters::cuda::canpeer()
{
	return _canpeer;
}
cudaStream_t& KingOfMonsters::cuda::__streams(int64_t i, int64_t j)
{
	return _streams[i][j];
}
void KingOfMonsters::cuda::dispose() {
#ifndef _CPU
	if (valid())
	{
		/*if (prevwn != 0)
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
		}*/
		for (int64_t ii = 0; ii < _count; ii++)
		{
			cudaSetDevice(ii);
			for (int64_t kk = 0; kk < STREAMCOUNT; kk++)
				cudaStreamDestroy(_streams[ii][kk]);
		}
		previdentiyN = 0;
		/*if (_L != 0)
			delete[] _L;
		_L = 0;*/
		/*if (prevN != 0)
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
		for (int64_t i = 0; i < MAXDEVICE; i++)
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
		__mgrhs2 = 0;*/
		for (int64_t i = 0; i < count(); i++)
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
			if (__work2[i] != 0)
			{
				cudaFree(__work2[i]);
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
			__work2[i] = 0;
			work_size[i] = 0;
			work_size2[i] = 0;
			for (int64_t j = 0; j < STREAMCOUNT; j++)
				if (solver_handle[i][j] != 0)
				{
					cusolverDnDestroy(solver_handle[i][j]);
					solver_handle[i][j] = 0;
				}
			for (int64_t j = 0; j < STREAMCOUNT; j++)
			{
				/*if (cusparse_handle[i][j] != 0)
				{
					cusparseDestroy(cusparse_handle[i][j]);
					cusparse_handle[i][j] = 0;
				}*/
				//if (solver_handleSp[i][j] != 0)
				//{
				//	cusolverSpDestroy(solver_handleSp[i][j]);
			//		solver_handleSp[i][j] = 0;
			//	}
			}
			//if (cublas_handle != 0)
			//{
			//	cublasDestroy(cublas_handle[i]);
			//}
			//cublas_handle[i] = 0;
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
#else
initialized = false;
#endif
}
/*
cusolverMgHandle_t KingOfMonsters::cuda::mgsolver() {
	return mg_solver;
}*/
#ifndef _CPU
int* KingOfMonsters::cuda::info(int i) {
	return __info[i];
}
#endif
KingOfMonsters::cuda::~cuda() {

	dispose();
}
#ifndef _CPU
double* KingOfMonsters::cuda::work_M(int i) {
	return __mgM[i];
}
double* KingOfMonsters::cuda::work_rhs(int i) {
	return __mgrhs[i];
}
double* KingOfMonsters::cuda::work_C(int i) {
	return __mgC[i];
}
double* KingOfMonsters::cuda::work(int64_t N, int device) {
	if (N < work_size[device])
	{
		//do nothing
		return __work[device];
	}
	else {
		if (__work[device] != 0)
		{
			cudaSetDevice(device);
			cudaFree(__work[device]);
		}
		int64_t _N = (int64_t)N * 1.3;
		cudaMalloc(&__work[device], _N * sizeof(double));
		work_size[device] = _N;
		return __work[device];
	}
}
double* KingOfMonsters::cuda::work2(int64_t N, int device) {
	if (N < work_size2[device])
	{
		//do nothing
		return __work2[device];
	}
	else {
		if (__work2[device] != 0)
		{
			cudaSetDevice(device);
			cudaFree(__work2[device]);
		}
		int64_t _N = (int64_t)N * 1.3;
		cudaMalloc(&__work2[device], _N * sizeof(double));
		work_size2[device] = _N;
		return __work2[device];
	}
}
double* KingOfMonsters::cuda::work(int64_t N, int device, cudaStream_t stream) {
	if (N < work_size[device])
	{
		return __work[device];
	}
	else {
		if (__work[device] != 0)
		{
			cudaSetDevice(device);
			cudaFree(__work[device]);
		}
		int64_t _N = (int64_t)N * 1.3;
		cudaMalloc(&__work[device], _N * sizeof(double));
		work_size[device] = _N;
		return __work[device];
	}
}
#endif
bool KingOfMonsters::cuda::valid() {
	//return true;
	return initialized && (!failed);
}
string  KingOfMonsters::cuda::device_name() {
#ifndef _CPU
	if (valid()) {
		char name[100];

		if (count() == 0) {
			failed = true;
			return "device not found";
		}
		std::stringstream ss;
		for (int64_t i = 0; i < count(); i++)
		{
			CUdevice dev;
			cuDeviceGet(&dev, i);
			cuDeviceGetName(name, 100, dev);
			ss << name << "::" << "SPEED<<" << speed[i] << std::endl;
		}
		if (count() > 1)
		{
			ss << "multiple GPUs available!" << std::endl;
			/*for (int64_t i = 0; i < count(); i++)
			{
				for (int64_t j = i + 1; j < count(); j++)
				{
					ss << "peer access" << "(" << i << "," << j << ")" << canpeeraccess(i, j) << std::endl;
				}
			}*/

		}
		return ss.str();
	}
	return "invalid";
#else
	return "CPU";
#endif
}
cusolverDnHandle_t& KingOfMonsters::cuda::solver(int64_t ii, int64_t kk) {
	return solver_handle[ii][kk];
}
//cusolverSpHandle_t& KingOfMonsters::cuda::solverSp(int64_t ii, int64_t kk) {
//	return solver_handleSp[ii][kk];
//}
//cublasHandle_t& KingOfMonsters::cuda::blas(int64_t ii) {
//	return cublas_handle[ii];
//}
int& KingOfMonsters::cuda::count() {
	return _count;
}
int& KingOfMonsters::cuda::fastest() {
#ifndef _CPU
	return _fastest;
#else
	_fastest = -1;
	return _fastest;
#endif
}
KingOfMonsters::_mySparse::_mySparse()
{
	//Eigen::initParallel();
	dat.reserve(1000);
	coeff.reserve(1000);
	_coeff.reserve(1000);
	_mat.reserve(1000);
	_mt = omp_get_max_threads();
	/*int _mt2 = 0;
#pragma omp parallel
	{
#pragma omp single
		_mt = omp_get_num_threads();
	}
	if (_mt2 > _mt)_mt = _mt2;*/
	Eigen::setNbThreads(_mt);
	omp_set_num_threads(_mt);
	//prevmat.resize(1, 1);
	//e = new Eigen::SparseMatrix<double>[200];
	//e2 = new Eigen::SparseMatrix<double>[200];
	//_smat.resize(1);
	//_smat = new Eigen::SparseMatrix<double>();
}
KingOfMonsters::_mySparse::~_mySparse()
{
	
	
}

// TODO: This is an example of a library function

KingOfMonsters::_myPermutation::_myPermutation(int64_t* ptr, int64_t N)
{
	//Eigen::VectorXi indices(N);
	Eigen::Vector<int64_t,Eigen::Dynamic> indices(N);
	for (int64_t i = 0; i < N; i++)
	{
		indices[i] = *ptr;
		ptr++;
	}
	perm.indices() = indices;
}

std::string KingOfMonsters::_mySparse::_testopenmp()
{
	std::stringstream ss;
	int64_t mt = omp_get_max_threads();
	ss << "num threads:" << mt << std::endl;
#pragma omp parallel for
	for (int64_t i = 0; i < 100; i++)
	{
		int64_t ct = omp_get_thread_num();
#pragma omp critical
		{
			ss << "current threads:" << ct << std::endl;
		}
		for (int64_t t = 0; t < ct * 100; t++)
		{
			auto f = new double[10];
			delete[] f;
		}
	}
	return ss.str();
}

void KingOfMonsters::_mySparse::freeze(bool _do) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, this->rows(), this->cols());
	//__r = this->rows(); __c = this->cols();
	_dmat.setZero(this->rows(), this->cols());
	if (_do)
		for (int64_t ii = 0; ii < _nt; ii++)
		{
			if(coeff[ii].size()!=0)
			_dmat += coeff[ii].asDiagonal() * this->_mat[ii];
		}
	else
		for (int64_t ii = 0; ii < _nt; ii++)
		{
			if(this->_mat[ii].rows()!=0)
			_dmat += this->_mat[ii];
		}
}
double KingOfMonsters::_mySparse::L2Norm(Eigen::VectorXd* a, Eigen::VectorXd* b) {
	return (*a).transpose() * this->_mat[0] * (*b);
}
Eigen::VectorXd KingOfMonsters::_mySparse::Vector(double* ptr1, int64_t N1) {
	auto a = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(ptr1, N1);
	return this->_mat[0] * a;
}
Eigen::VectorXd KingOfMonsters::_mySparse::Vector(Eigen::VectorXd* vec) {
	return this->_mat[0] * *vec;
}


void KingOfMonsters::_mySparse::Vector(Eigen::VectorXd* vec, Eigen::VectorXd* ret) {
	*ret = this->_mat[0] * *vec;
}
void KingOfMonsters::_mySparse::plus(_mySparse* m, double sc, bool dense, bool sparse) {
	if (sparse)
		this->_mat[0] += m->_mat[0] * sc;
	if (dense)
	{
		//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
		//this->_mat[0] = this->_mat[0] + m->_mat[0] * sc;
		_dmat += m->_mat[0] * sc;
	}
}
void KingOfMonsters::_mySparse::plus(Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> *m) {
	this->_mat[0] += *m;
}
void KingOfMonsters::_mySparse::setmat(Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> &mat, int64_t ii) {
	this->_mat[ii] = mat;
}
void KingOfMonsters::_mySparse::setmat(const Eigen::MatrixXd& mat) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, mat.rows(), mat.cols());
	//__r = mat.rows();
	//__c = mat.cols();
	_dmat.resize(mat.rows(), mat.cols());
	_dmat = mat;
}
double KingOfMonsters::_mySparse::at(int64_t i, int64_t ii) {
	return this->_mat[ii].data().value(i);
}
double KingOfMonsters::_mySparse::__at(int64_t i, int64_t j) {
	return this->_mat[0].coeffRef(i,j);
}
int64_t KingOfMonsters::_mySparse::num_elem(int64_t ii)
{
	return this->_mat[ii].data().size();
}
double KingOfMonsters::_mySparse::_at(int64_t i) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	return _dmat.data()[i];
}
double KingOfMonsters::_mySparse::_at(int64_t i, int64_t j) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	return _dmat(i, j);
}
int64_t KingOfMonsters::_mySparse::cols() {
	return _mat[0].cols();
}
void KingOfMonsters::_mySparse::_resize(int64_t n, int64_t m) {

	//Eigen::Map<Eigen::MatrixXd> map1(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> map2(___dmat, n, m);
	_dmat.conservativeResize(n, m);
	//__r = n;
	//__c = m;
}
void KingOfMonsters::_mySparse::setmiddlecolum(Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>& f, int64_t start, int64_t end) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	_dmat.middleCols(start, end - start) = f;
}
void KingOfMonsters::_mySparse::permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int64_t>& perm)
{
	for (int64_t ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii] * perm.transpose();
}

void KingOfMonsters::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int64_t>& perm, bool sparse, bool dense)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> _tmp(___dmat+__r*__c, __r, __c);
	int64_t nn = _dmat.rows();// __r;
	int64_t numthreads = 0;
	//_mt=omp_get_max_threads();
	int64_t S = nn / _mt / 2;
	auto pt = perm.transpose();
	if (sparse||dense)
		if (_mat.size() >= 1)
		{
			//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
			{
				_mat[0] = perm * (_mat[0]) * pt;
			}
		}
	if (false)
	{
		//_tmp.noalias() = _dmat * pt;
		//_dmat.noalias() = perm * _tmp;
		//_dmat.applyOnTheLeft(perm);
		//_dmat.applyOnTheRight(perm.transpose());
		//prrm.transpose().applyThisOnTheRight(_dmat);
		//perm.applyThisOnTheLeft(_dmat);

#pragma omp parallel for
		for (int64_t i = 0; i < nn; i += S)
		{
			int64_t start = i;
			int64_t end = i + S;
			if (end > nn)end = nn;
			_dmat.middleRows(start, end - start).applyOnTheRight(pt);// = _dmat.middleRows(start, end - start) * pt;
			//_tmp.middleRows(start, end - start).noalias()= _dmat.middleRows(start, end - start) * pt;
			//perm.transpose().applyThisOnTheRight(_dmat.middleRows(start, end - start));
		}

#pragma omp parallel for
		for (int64_t i = 0; i < nn; i += S)
		{
			int64_t start = i;
			int64_t end = i + S;
			if (end > nn)end = nn;
			//_dmat.middleCols(start, end - start).noalias() = perm * _dmat.middleCols(start, end - start);
			_dmat.middleCols(start, end - start).applyOnTheLeft(perm);// = perm * _dmat.middleCols(start, end - start);
			//perm.applyThisOnTheLeft(_dmat.middleCols(start, end - start));
		}
	}

}

void KingOfMonsters::_mySparse::_permuteCols(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int64_t>& perm, bool sparse, bool dense)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> _tmp(___dmat+__r*__c, __r, __c);
	int64_t nn = _dmat.rows();// __r;
	int64_t numthreads = 0;
	//_mt=omp_get_max_threads();
	int64_t S = nn / _mt / ((int64_t)2);
	auto pt = perm.transpose();
	if (sparse || dense)
		if (_mat.size() >= 1)
		{
			//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
			{
				_mat[0] = (_mat[0]) * pt;
			}
		}
	if (false)
	{
		//_tmp.noalias() = _dmat * pt;
		//_dmat.noalias() = perm * _tmp;
		//_dmat.applyOnTheLeft(perm);
		//_dmat.applyOnTheRight(perm.transpose());
		//prrm.transpose().applyThisOnTheRight(_dmat);
		//perm.applyThisOnTheLeft(_dmat);

#pragma omp parallel for
		for (int64_t i = 0; i < nn; i += S)
		{
			int64_t start = i;
			int64_t end = i + S;
			if (end > nn)end = nn;
			_dmat.middleRows(start, end - start).applyOnTheRight(pt);// = _dmat.middleRows(start, end - start) * pt;
			//_tmp.middleRows(start, end - start).noalias()= _dmat.middleRows(start, end - start) * pt;
			//perm.transpose().applyThisOnTheRight(_dmat.middleRows(start, end - start));
		}

#pragma omp parallel for
		for (int64_t i = 0; i < nn; i += S)
		{
			int64_t start = i;
			int64_t end = i + S;
			if (end > nn)end = nn;
			//_dmat.middleCols(start, end - start).noalias() = perm * _dmat.middleCols(start, end - start);
			_dmat.middleCols(start, end - start).applyOnTheLeft(perm);// = perm * _dmat.middleCols(start, end - start);
			//perm.applyThisOnTheLeft(_dmat.middleCols(start, end - start));
		}
	}

}
void KingOfMonsters::_mySparse::_permuteRows(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm, bool sparse, bool dense)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> _tmp(___dmat+__r*__c, __r, __c);
	int64_t nn = _dmat.cols();// __r;
	int64_t numthreads = 0;
	//_mt=omp_get_max_threads();
	int64_t S = nn / _mt / ((int64_t)2);
	//auto pt = perm.transpose();
	if (sparse || dense)
		if (_mat.size() >= 1)
		{
			//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
			{
				_mat[0] = perm * (_mat[0]);// *pt;
			}
		}
/*	if (false)
	{
		//_tmp.noalias() = _dmat * pt;
		//_dmat.noalias() = perm * _tmp;
		//_dmat.applyOnTheLeft(perm);
		//_dmat.applyOnTheRight(perm.transpose());
		//prrm.transpose().applyThisOnTheRight(_dmat);
		//perm.applyThisOnTheLeft(_dmat);

#pragma omp parallel for
		for (int64_t i = 0; i < nn; i += S)
		{
			int64_t start = i;
			int64_t end = i + S;
			if (end > nn)end = nn;
			_dmat.middleRows(start, end - start).applyOnTheRight(pt);// = _dmat.middleRows(start, end - start) * pt;
			//_tmp.middleRows(start, end - start).noalias()= _dmat.middleRows(start, end - start) * pt;
			//perm.transpose().applyThisOnTheRight(_dmat.middleRows(start, end - start));
		}

#pragma omp parallel for
		for (int64_t i = 0; i < nn; i += S)
		{
			int64_t start = i;
			int64_t end = i + S;
			if (end > nn)end = nn;
			//_dmat.middleCols(start, end - start).noalias() = perm * _dmat.middleCols(start, end - start);
			_dmat.middleCols(start, end - start).applyOnTheLeft(perm);// = perm * _dmat.middleCols(start, end - start);
			//perm.applyThisOnTheLeft(_dmat.middleCols(start, end - start));
		}
	}
	*/
}
void KingOfMonsters::_mySparse::shrink(int64_t M)
{
	for (int64_t ii = 0; ii < _nt; ii++)
		_mat[ii] = _mat[ii].leftCols(M);
}
void KingOfMonsters::_mySparse::_shrink(int64_t M, bool sparse, bool dense)
{
	if (true)
	{
		if (_mat.size() >= 1)
		{
			//if ((_mat[0].rows() == _dmat.rows()) && (_mat[0].cols() == _dmat.cols()))
			{
				_mat[0] = _mat[0].block(0, 0, M, M);
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
	/*int64_t count = 0;
	int64_t count2 = 0;

	for (int64_t i = 0; i < M; i++)
	{
		if (_mat[0].coeff(i, i) < 0)count++;
		if(dense)if (_dmat.coeff(i, i) < 0)count2++;

	}*/
}
void KingOfMonsters::_mySparse::_shrinkCols(int64_t M, bool sparse, bool dense)
{
	if (true)
	{
		if (_mat.size() >= 1)
		{
			//if ((_mat[0].rows() == _dmat.rows()) && (_mat[0].cols() == _dmat.cols()))
			{
				_mat[0] = _mat[0].leftCols(M);
			}
		}
	}
	if (dense)
	{
		//__r = M;
		//__c = M;
		//Eigen::Map<Eigen::MatrixXd> map1(___dmat, __r, __c);
		_dmat.resize(_dmat.rows(), M);
		_dmat = _mat[0];
		//_resize(M, M);
		//_dmat.conservativeResize(M, M);
	}
	/*int64_t count = 0;
	int64_t count2 = 0;

	for (int64_t i = 0; i < M; i++)
	{
		if (_mat[0].coeff(i, i) < 0)count++;
		if(dense)if (_dmat.coeff(i, i) < 0)count2++;

	}*/
}
void KingOfMonsters::_mySparse::_permute(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int64_t>& perm, Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic,int64_t>& perm2)
{
	
	auto pt = perm2.transpose();

	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0] = perm * (_mat[0]) * pt;
		}
	}
}
void KingOfMonsters::_mySparse::__permuteRows(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm)
{

	

	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0] = perm * (_mat[0]) ;
		}
	}
}
void KingOfMonsters::_mySparse::__permuteCols(Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int64_t>& perm)
{

	auto pt = perm.transpose();

	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0] =  (_mat[0]) * pt;
		}
	}
}
void KingOfMonsters::_mySparse::_shrink(int64_t M, int64_t N,bool sparse,bool dense)
{
	if (_mat.size() >= 1)
	{
		//if (_mat[0].rows() == _dmat.rows() && _mat[0].cols() == _dmat.cols())
		{
			_mat[0].conservativeResize(M, N);
		}
	}
	if (dense)
	{
		
			//__r = M;
			//__c = M;
			//Eigen::Map<Eigen::MatrixXd> map1(___dmat, __r, __c);
			_dmat.resize(M,N);
			_dmat = _mat[0];
			//_resize(M, M);
			//_dmat.conservativeResize(M, M);
	
	}
	//_dmat.conservativeResize(M, N);//;// _dmat.topLeftCorner(M, N);// f;
}

int64_t KingOfMonsters::_mySparse::rows() {
	int64_t _ret = 0;
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		_ret += _mat[ii].rows();
	}
	return _ret;
}
void  KingOfMonsters::_mySparse::scale(int i, double sc)
{
	if (_mat.size() > 0)
	{
		_mat[0].row(i) *= sc;
	}
}
void  KingOfMonsters::_mySparse::scale(double sc)
{
	_mat[0] *= sc;
}
void  KingOfMonsters::_mySparse::_scale(double sc)
{
	_dmat *= sc;
}
int64_t KingOfMonsters::_mySparse::_rows() {
	return _dmat.rows();// __r;
}
int64_t KingOfMonsters::_mySparse::_cols() {
	return _dmat.cols();// __c;
}
int64_t KingOfMonsters::_mySparse::__rows() {
	int64_t _ret = 0;
	for (int64_t ii = 0; ii < this->_coeff.size(); ii++)
	{
		_ret += _coeff[ii].size();
	}
	return _ret;
}
void KingOfMonsters::_mySparse::Clear() {
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
	//for (int64_t ii = 0; ii < _nt; ii++)
	{
		//this->coeff[ii].setZero();
		//this->_mat[ii].setZero();
	}
}

void KingOfMonsters::_mySparse::init(int64_t n, int64_t m)
{
	this->_nt = 1;

	resize(n, m);
}
int64_t KingOfMonsters::_mySparse::resize(int64_t n, int64_t m) {
	this->_nt = this->dat.size();
	if (_nt == 0 || _nt == 1)
	{
		_nt = 1;
		_coeff.resize(_nt);
		dat.resize(_nt);
		_mat.resize(_nt);
		coeff.resize(_nt);
		_mat[0].resize(n, m);
		_mat[0].setZero();
		_mat[0].makeCompressed();
		//_mat[0].reserve(n * m / 10);
	}

	return this->_nt;
}
void KingOfMonsters::_mySparse::reserve(int64_t n) {
	_mat[0].reserve(n);
}
void KingOfMonsters::_mySparse::addemptyrow(int64_t ii) {
	eigen_assert(_coeff[0].size() == ii );
	_coeff[0].push_back(1);
}
void KingOfMonsters::_mySparse::addrow(int64_t ii, int64_t* ptr, double* data, double sc, int64_t N)
{
	addrow(ii, ptr, data, 0, sc, N, true, 1.0);
}
void KingOfMonsters::_mySparse::addrow(int64_t ii, int64_t* ptr, double* data, double sc, int64_t N,double coeff)
{
	addrow(ii, ptr, data, 0, sc, N, true, coeff);
}
void KingOfMonsters::_mySparse::addrow(int64_t ii, int64_t* ptr, double* data, double* data2,double sc, int64_t N, double c1,double c2)
{
	//data += shift;
	//data += shift;
	//ptr += shift;
	if (dat.size() != 1)dat.resize(1);
	if (_coeff.size() != 1)_coeff.resize(1);
	dat[0].reserve(dat[0].size() + N);
	for (int64_t i = 0; i < N; i++)
	{
		dat[0].push_back(Eigen::Triplet<double>(ii, *ptr, (*data) * c1+(*data2)*c2));
		ptr++;
		data++;
		data2++;
	}
	//if (add)
	{
		eigen_assert(ii == _coeff[0].size());
		_coeff[0].push_back(sc);
	}
}

void KingOfMonsters::_mySparse::addrow(int64_t ii, int64_t* ptr, double* data, int64_t shift, double sc, int64_t N, bool add,double __coeff)
{
	data += shift;
	ptr += shift;
	if (dat.size() != 1)dat.resize(1);
	if (_coeff.size() != 1)_coeff.resize(1);
	dat[0].reserve(dat[0].size() + N);
	for (int64_t i = 0; i < N; i++)
	{
		dat[0].push_back(Eigen::Triplet<double>(ii, *ptr, (*data)*__coeff));
		ptr++;
		data++;
	}
	if (add)
	{
		eigen_assert(ii == _coeff[0].size());
		_coeff[0].push_back(sc);
	}
}

void KingOfMonsters::_mySparse::adddat(int64_t ii, int64_t j, double value)
{
	dat[0].push_back(Eigen::Triplet<double>(ii, j, value));
}
void KingOfMonsters::_mySparse::addcoeff(double sc) {
	_coeff[0].push_back(sc);
}
Eigen::VectorXd KingOfMonsters::_mySparse::get_coeff(int64_t ii) {
	return this->coeff[ii];
}
void  KingOfMonsters::_mySparse::copycoefffrom(KingOfMonsters::_mySparse* mat)
{
	this->_nt = mat->_nt;
	this->coeff.resize(_nt);
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		this->coeff[ii] = Eigen::VectorXd(mat->get_coeff(ii));
	}
}

void KingOfMonsters::_mySparse::begin_construct()
{
	_dat_count = 0;
}
double KingOfMonsters::_mySparse::mulboth(Eigen::VectorXd& u,Eigen::VectorXd& v)
{
	return (u.transpose() * this->_mat[0]) * v;
}
	
void KingOfMonsters::_mySparse::mulright(Eigen::VectorXd& v, Eigen::VectorXd& ret,double sc)
{
	ret+=this->_mat[0] * v*sc;
}
void KingOfMonsters::_mySparse::setzero(int row)
{
	this->_mat[0].prune([row] (int i, int j, float val) {return i != row; });
}
void KingOfMonsters::_mySparse::mulleft(Eigen::VectorXd& v, Eigen::VectorXd& ret, double sc)
{
	ret+=(v.transpose() * this->_mat[0]).transpose()*sc;
}
void KingOfMonsters::_mySparse::end_construct(int64_t cc)
{
	_nt = _dat_count;
	//_mat.clear();
	//_mat.shrink_to_fit();
	_mat.resize(_nt);
	coeff.resize(_nt);
#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		_mat[ii].resize(_coeff[ii].size(), cc);
	}
}

void KingOfMonsters::_mySparse::addmat(_mySparse* mat)
{
	_dat_count++;
	if (dat.size() < _dat_count)dat.resize(_dat_count);
	if (_coeff.size() < _dat_count)_coeff.resize(_dat_count);
	dat[_dat_count - 1] = mat->dat[0];
	_coeff[_dat_count - 1] = (mat->_coeff[0]);
}
void KingOfMonsters::_mySparse::OfDuplicate(_mySparse* mat)
{
	this->_nt = mat->_nt;
	this->_mat.resize(_nt);
	this->_coeff.resize(_nt);
	this->coeff.resize(_nt);
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		this->_mat[ii].resize(mat->_mat[ii].rows(), mat->_mat[ii].cols());
		this->_mat[ii].reserve(mat->_mat[ii].data().allocatedSize());
		this->_mat[ii] = mat->_mat[ii];// M;
		
		this->coeff[ii].resize(mat->coeff[ii].size());
		this->_coeff[ii].resize(mat->coeff[ii].size());

		this->coeff[ii] =  mat->coeff[ii];
		this->_coeff[ii] = mat->_coeff[ii];
	}
}
void KingOfMonsters::_mySparse::_OfDuplicate(_mySparse* mat)
{
	if (this->_mat.size() == 0)
		this->_mat.resize(1);
	this->_mat[0] = mat->_mat[0];
	//this->_dmat = this->_mat[0];// mat->_dmat;
}
void KingOfMonsters::_mySparse::ofDat()
{
	if (_mat.size() != _nt)_mat.resize(_nt);
#pragma omp parallel for schedule(dynamic,2)
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		_mat[ii].setZero();
		//_mat[ii].reserve(dat[ii].size()*2);
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

std::map<Eigen::SparseMatrix<double, 0, int64_t>*, std::map<std::tuple<int64_t, int64_t>,int>> map2;
void KingOfMonsters::_mySparse::makePattern()
{

	std::map<std::tuple<int64_t, int64_t>,int> __map2;
	for (int64_t k = 0; k < _mat[0].outerSize(); ++k) {
		for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(_mat[0], k); it; ++it) {
			int64_t S = _mat[0].outerIndexPtr()[k];
			int64_t E = _mat[0].outerIndexPtr()[k + 1];

			for (int64_t tt = S; tt < E; tt++)
			{
				if (_mat[0].innerIndexPtr()[tt] == it.row())
					__map2[std::tuple<int64_t,int64_t>(it.row(), it.col())] = tt;
			}
		}
	}
	map2[&this->_mat[0]] = __map2;
}
void KingOfMonsters::_mySparse::computeQR() {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	int count = 0;
	int nnz = 0;
	for (auto f : _mat)
	{
		count += f.rows();
		nnz += f.nonZeros();
	}

	Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t> M(count, _mat[0].cols());
	Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> N(count, _mat[0].cols());
	M.reserve(nnz);
	int offset = 0;
	for (auto f : _mat)
	{
		M.middleRows(offset, f.rows()) = f;
		offset += f.rows();
	}
	N = M;


	Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::RowMajor, int64_t>, Eigen::COLAMDOrdering<int64_t>> qr;
	qr.setPivotThreshold(0.00000000001);
	qr.compute(N);

}
void KingOfMonsters::_mySparse::freezecoeff() {
//#pragma omp parallel for schedule(dynamic,1)
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		coeff[ii] = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(_coeff[ii].data(), _coeff[ii].size());
	}
}
int64_t KingOfMonsters::_mySparse::numBlocks()
{
	return this->dat.size();
}


std::string KingOfMonsters::_mySparse::ofAtA( _mySparse* A, bool sparse)
{	
	static std::vector<std::vector<int64_t>> index;
	static std::map<_mySparse*,std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>>> ___e;
	static std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> e2;
	auto ss = std::stringstream();
	auto now = high_resolution_clock::now();
	int64_t nn = this->cols();
	int64_t mm = this->cols();
	//int64_t _mt = omp_get_max_threads();
	//omp_set_num_threads(_mt);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	int64_t __mt = _mt;// std::min(_nt / 10, _mt * 10);
	std::vector < Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>>* e;
	if (___e.contains(this))e = &___e[this]; else { std::vector < Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> _e; ___e[this] = _e; e = &___e[this]; }
	if (e->size() < __mt)
	{
		e->resize(__mt);
		e2.resize(__mt);
	}
	//Eigen::initParallel();
	//Eigen::setNbThreads(1);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	ss << _nt << ":nt:" << std::endl;
	ss << _mt << ":_mt:" << std::endl;
	//Eigen::initParallel();
	int64_t job = 0;
	int64_t sss = _nt / __mt / 8;
	int64_t __nt2 = _nt;
	if (sss == 0)sss = 1;
	Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>* prevmat;
	if (dict2.contains(this)) {
		prevmat = &dict2[this];
	}
	else {
		Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> _prevmat(nn, nn);
		dict2[this] = _prevmat;
		prevmat = &dict2[this];
		prevmat->resize(nn, nn);
		prevmat->reserve(nn * nn / 10);
	}
	//prevmat->makeCompressed();
	std::vector<int64_t>* _map = 0;
	if (map.contains(prevmat))
	{
		_map = &map[prevmat];
	}
#pragma omp parallel for
	for (int64_t i = 0; i < __mt; i++) {
		if (_map==0||(*e)[i].nonZeros()!=prevmat->nonZeros())
		{
			(*e)[i].resize(nn, mm);
			(*e)[i] = *prevmat;
			(*e)[i].makeCompressed();
		}
		memset((*e)[i].valuePtr(), 0, sizeof(double) * prevmat->nonZeros());
		e2[i].resize(nn, mm);
		e2[i].reserve(nn * mm / 100);
	}
	index.resize(__mt);
#pragma omp parallel for
	for (int64_t _ii = 0; _ii < __mt; _ii++)
	{
		int64_t S = 0;
		int64_t E = 0;
		//int64_t K = 0;
		//auto _e = e[_ii];
		for (int64_t tt = 0; tt < 4000; tt++)
		{
#pragma omp critical
			{
				S = job;
				E = job + sss;
				if (E > _nt)E = _nt;
				job = E;
			}
			if (S >= _nt)break;
			for (int64_t ii = S; ii < E; ii++)
			{
				
				//(*e)[_ii] += this->_mat[ii].transpose() * coeff[ii].asDiagonal() * this->_mat[ii];
//#pragma omp critical
				
					if (_map == 0)
					{
						e2[_ii] = this->_mat[ii].transpose() * coeff[ii].asDiagonal() * this->_mat[ii];
						(*e)[_ii] += e2[_ii];
					}
					else {
						auto cc = coeff[ii].asDiagonal();
						
						e2[_ii] = this->_mat[ii].transpose() * coeff[ii].asDiagonal() * this->_mat[ii];
						if (e2[_ii].nonZeros() > 0)
						{
							//int64_t* ptr = &index[_ii][0];
							for (int64_t k = 0; k < mm; ++k) {
								for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(e2[_ii], k); it; ++it) {
									//e[0].coeffRef(it.row(), it.col()) += it.value();
									*((*e)[_ii].valuePtr() + (*_map)[it.row() * mm + it.col()]) += it.value();
									//ptr++;
								}
							}
						}
						//e[_ii].makeCompressed();
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
	//Eigen::setNbThreads(_mt);

	for (int64_t tt = 0; tt < 4000; tt++)
	{
#pragma omp parallel for schedule(static,1)
		for (int64_t i = 0; i < __mt; i += 2)
		{
			if (i + 1 < __mt) {
				if (_map == 0 || (*e)[i].nonZeros() != (*e)[i + 1].nonZeros())
				{
					(*e)[i] += (*e)[i + 1];
				}
				else {
					Eigen::Map<Eigen::VectorXd> map1((*e)[i].valuePtr(), (*e)[i].nonZeros());
					Eigen::Map<Eigen::VectorXd> map2((*e)[i + 1].valuePtr(), (*e)[i].nonZeros());
					map1 += map2;
				}
			}
		}
		int64_t _ct = 0;
#pragma omp parallel for ordered schedule(static,1)
		for (int64_t i = 0; i < __mt; i += 2)
		{
#pragma omp ordered
			if (_map == 0 || (*e)[i].nonZeros() != (*e)[i / 2].nonZeros())
			{
				(*e)[i / 2] = (*e)[i];
			}
			else
			{
				memcpy((*e)[i / 2].valuePtr(), (*e)[i].valuePtr(), sizeof(double) * (*e)[i].nonZeros());
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
		this->_mat[0] = (*e)[0];
		//for (int64_t i = 1; i < __mt; i++) {
		//	this->_mat[0] += e[i];
		//}
		this->_mat[0].makeCompressed();
		if (_map == 0 || prevmat->nonZeros() != this->_mat[0].nonZeros())
		{
			*prevmat = this->_mat[0];
			prevmat->makeCompressed();
			dict2[this] = *prevmat;
			//build map
			std::vector<int64_t> __map;
			__map.resize(nn*nn);
			for (int64_t k = 0; k < prevmat->outerSize(); ++k) {
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(*prevmat, k); it; ++it) {
					int64_t S = prevmat->outerIndexPtr()[k];
					int64_t E = prevmat->outerIndexPtr()[k + 1];

					for (int64_t tt = S; tt < E; tt++)
					{
						if (prevmat->innerIndexPtr()[tt] == it.row())
							__map[it.row()*mm+it.col()] = tt;
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
		
		//__r = nn;
		//__c = nn;
		
		//this->_dmat = x;
	}
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	//this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	ss << "sum:" << (*e)[0].sum() << std::endl;
	return ss.str();
}
int KingOfMonsters::_mySparse::find_location(int64_t i, int64_t j)
{
	return map2[&this->_mat[0]][std::tuple<int64_t, int64_t>(i, j)];
}
void KingOfMonsters::_mySparse::add_usemap(int64_t i, int64_t j, double val)
{
	*(_mat[0].valuePtr() + (map2[&this->_mat[0]][std::tuple<int64_t, int64_t>(i, j)])) += val;

}
void KingOfMonsters::_mySparse::set_usemap(int64_t i, int64_t j, double val)
{
	*(_mat[0].valuePtr() + (map2[&this->_mat[0]][std::tuple<int64_t, int64_t>(i, j)])) = val;
}
void KingOfMonsters::_mySparse::add_uselocation(int64_t location, double val)
{
	*(_mat[0].valuePtr() + (location)) += val;

}
void KingOfMonsters::_mySparse::set_uselocation(int64_t location, double val)
{
	*(_mat[0].valuePtr() + (location)) += val;

}


std::string KingOfMonsters::_mySparse::ofAtA_gpu(cuda* _cuda, _mySparse* A, bool sparse)
{
#ifndef _CPU
	auto prevmat = &dict2[this];
	auto __index = map[prevmat];
	int64_t nnz = prevmat->nonZeros();
	int64_t* cols = new int64_t[nnz];
	int64_t* rows = new int64_t[nnz];
	
	for (int64_t i = 0; i < _nt; i++)
	{
	}
		return "";
#else
	return "";
#endif

}


void KingOfMonsters::_mySparse::_freeze() {
	//this->_dmat = this->_mat[0];
}
std::string KingOfMonsters::_mySparse::_ofAtA(_mySparse* A)
{

	//__r = A->cols();
	//__c = A->cols();
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, A->cols(), A->cols());
	std::stringstream ss;
#ifdef _DEBUG
	ss << _nt << std::endl;
	for (int64_t ii = 0; ii < _nt; ii++)
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
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{

			auto _ret = (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
			_dmat += _ret;
		}
	}
	return ss.str();
}
std::string KingOfMonsters::_mySparse::_ofAtA_sparse(_mySparse* A)
{

	//__r = A->cols();
	//__c = A->cols();
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, A->cols(), A->cols());
	std::stringstream ss;
#ifdef _DEBUG
	ss << _nt << std::endl;
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{
			ss << ii << "th" << std::endl;
			ss << A->_mat[ii].rows() << "," << A->_mat[ii].cols() << "," << coeff[ii].size() << std::endl;

		}
	}
#endif

	//this->_dmat.resize(A->cols(), A->cols());
	_mat[0].resize(A->cols(), A->cols());
	_mat[0].setZero();
	//_tmp.setZero(__r, __c);
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{

			auto _ret = (A->_mat[ii].transpose() * coeff[ii].asDiagonal() * A->_mat[ii]);
			_mat[0] += _ret;
		}
	}
	return ss.str();
}
std::string KingOfMonsters::_mySparse::info()
{
	return "viennacl has been abandoned";
}


void KingOfMonsters::_mySparse::_ofAtB(_mySparse* B, _mySparse* C)
{
	this->freeze(true);
	B->freeze(false);

	int64_t nn = this->_cols();
	int64_t mm = B->_cols();
	
	//C->__r = nn;
	//C->__c = mm;
	C->_dmat.resize(nn, mm);
	//C->_tmp.resize(nn, mm);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, nn, mm);
	//Eigen::Map<Eigen::MatrixXd> C_dmat(C->___dmat, nn, mm);
	//C_dmat.setZero();

	//int64_t mt = omp_get_max_threads();

	int64_t ss = mm / _mt;
	auto left = _dmat.transpose();
	auto right = B->_mat[0];

#pragma omp parallel for
	for (int64_t ii = 0; ii < mm; ii += ss)
	{
		int64_t S = ii;
		int64_t E = ii + ss;
		if (E >= mm)E = mm;
		C->_dmat.middleCols(S, E - S) = left * right.middleCols(S, E - S);
	}
	if (C->_mat.size() == 0)C->_mat.resize(1);
	C->_mat[0] = C->_dmat.sparseView(1.0, 0.0000000000001);
	C->_mat[0].makeCompressed();
}
/*void KingOfMonsters::_mySparse::_ofBtAB_qr(_mySparse* B, Eigen::VectorXd* b, _mySparse* C, Eigen::VectorXd* ret)
{
	this->join();
	B->join();

	static Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> _a;
	static Eigen::SparseMatrix<double,Eigen::ColMajor> q;
	static Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t > tmp;
	this->_mat[0].makeCompressed();
	_a = this->_mat[0];
	Eigen::SparseQR<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>, Eigen::COLAMDOrdering<int64_t>> qr;
	qr.setPivotThreshold(0.00000000001);
	qr.compute(_a);
	//q = qr.m_Q;
	//tmp = q * B->_mat[0];
	//C->_dmat = tmp.transpose() * tmp;
}*/


void KingOfMonsters::_mySparse::_ofBtAB(_mySparse* B,/* Eigen::VectorXd* b, */_mySparse* C/*, Eigen::VectorXd* ret*/)
{
	Eigen::MatrixXd D;

	int64_t nn = B->_mat[0].cols();
	int64_t kk = _dmat.cols();// __c;

	C->_dmat.setZero(nn, nn);

	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];

	D.resize(nn, kk);
	int64_t ss = kk / _mt / 2;

#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < kk; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= kk)_E = kk;
		D.middleCols(_S, _E - _S).noalias() = left * mid.middleCols(_S, _E - _S);
	}
	//D.noalias() = left * mid;
	ss = nn / _mt / 2;
#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < nn; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= nn)_E = nn;
		C->_dmat.middleCols(_S, _E - _S).noalias() = D * right.middleCols(_S, _E - _S);
	}


	//*ret = D.transpose() * *b;
}
void KingOfMonsters::_mySparse::_ofBtAB(_mySparse* B, _mySparse* B2,/* Eigen::VectorXd* b, */_mySparse* C/*, Eigen::VectorXd* ret*/)
{
	Eigen::MatrixXd D;
	Eigen::MatrixXd D2;

	int64_t nn = B->_mat[0].cols();
	int64_t mm = B2->_mat[0].cols();
	int64_t kk = _dmat.cols();

	C->_dmat.setZero(nn+mm, nn+mm);

	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];
	auto left2 = B2->_mat[0].transpose();
	auto right2 = B2->_mat[0];

	D.resize(nn, kk);
	D2.resize(mm, kk);
	int64_t ss = kk / _mt / 2;

//#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < kk; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= kk)_E = kk;
		D.middleCols(_S, _E - _S).noalias() = left * mid.middleCols(_S, _E - _S);
	}
//#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < kk; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= kk)_E = kk;
		D2.middleCols(_S, _E - _S).noalias() = left2 * mid.middleCols(_S, _E - _S);
	}
	//D.noalias() = left * mid;
//#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < nn; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= nn)_E = nn;
		C->_dmat.middleCols(_S, _E - _S).topRows(nn).noalias() = D * right.middleCols(_S, _E - _S);
	}
	ss = mm / _mt / 2;
//#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < mm; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= mm)_E = mm;
		C->_dmat.middleCols(_S+nn, _E - _S).bottomRows(mm).noalias() = D2 * right2.middleCols(_S, _E - _S);
	}
	ss = nn / _mt / 2;
//#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < nn; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= nn)_E = nn;
		C->_dmat.middleCols(_S, _E - _S).bottomRows(mm).noalias() = D2 * right.middleCols(_S, _E - _S);
	}
	ss = mm / _mt / 2;
//#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < mm; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= mm)_E = mm;
		C->_dmat.middleCols(_S+nn, _E - _S).topRows(nn).noalias() = D * right2.middleCols(_S, _E - _S);
	}


	//*ret = D.transpose() * *b;
}
void KingOfMonsters::_mySparse::_ofCBtAB(_mySparse* B, _mySparse* C, _mySparse* D)
{
	Eigen::MatrixXd E;

	int64_t nn = B->_mat[0].cols();
	int64_t kk = _dmat.cols();// __c;



	//C->__r = nn;
	//C->__c = nn;
	//C->_dmat.resize(nn, mm);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> C_dmat(C->___dmat, nn, nn);


	//C->_dmat.resize(nn, nn);
	D->_dmat.setZero(nn, nn);
	//C->_tmp.setZero(nn,nn);

	//int64_t mt = omp_get_max_threads();

	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];

	E.resize(nn, kk);
	int64_t ss = kk / _mt / 2;

#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < kk; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= kk)_E = kk;
		E.middleCols(_S, _E - _S).noalias() = left * mid.middleCols(_S, _E - _S);
	}
	//D.noalias() = left * mid;
	ss = nn / _mt / 2;
#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < nn; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= nn)_E = nn;
		D->_dmat.middleCols(_S, _E - _S).noalias() = E * right.middleCols(_S, _E - _S);
	}
	D->_dmat = C->_dmat * D->_dmat;

	//Eigen::EigenSolver<Eigen::MatrixXd> eigen;
	//eigen.compute(D->_dmat);

	//singularvalues->__v.resize(D->_dmat.rows());
	//singularvalues->__v = eigen.eigenvalues().real();
}
void KingOfMonsters::_mySparse::_ofCBtAB2(_mySparse* B, _mySparse* C, _mySparse* D, _mySparse* E)//this:invZ, C:invW, B:Jx, E:JX, D:result
{
	Eigen::MatrixXd F;

	int64_t nn = B->_mat[0].rows();
	
	/*Eigen::SparseMatrix<double> I(nn, nn);
	I.setIdentity();

	D->_dmat=B->_mat[0].transpose()* C->_dmat* B->_mat[0] * E->_mat[0].transpose()* this->_dmat* E->_mat[0];
	*/
	int64_t mm = B->_mat[0].cols();
	int64_t kk = E->_mat[0].cols();
	//int64_t kk = _dmat.cols();// __c;




	D->_dmat.setZero(mm+kk, mm+kk);
	Eigen::SparseMatrix<double> I(nn, nn);
	I.setIdentity();
	I *= 0.0000000000001;
 	Eigen::MatrixXd ii = (I+B->_mat[0] * B->_mat[0].transpose()+ E->_mat[0] * E->_mat[0].transpose()).toDense().inverse();

	D->_dmat.topLeftCorner(mm, mm) = B->_mat[0].transpose() * ii * B->_mat[0];
	D->_dmat.bottomRightCorner(kk, kk) = E->_mat[0].transpose() * ii * E->_mat[0];

	D->_dmat.topRightCorner(mm, kk) = B->_mat[0].transpose() * ii * E->_mat[0];
	D->_dmat.bottomLeftCorner(kk, mm) = E->_mat[0].transpose() * ii * B->_mat[0];
	Eigen::SparseMatrix<double> I2(mm+kk, mm+kk);
	I2.setIdentity();
	D->_dmat = I2 - D->_dmat;
	

	//D->_dmat = (B->_mat[0] * this->_dmat* B->_mat[0].transpose()*E->_mat[0] * C->_dmat * E->_mat[0].transpose());

/*	auto left = B->_mat[0].transpose();
	auto mid = _dmat;
	auto right = B->_mat[0];

	F.resize(nn, kk);
	int64_t ss = kk / _mt / 2;

#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < kk; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= kk)_E = kk;
		F.middleCols(_S, _E - _S).noalias() = left * mid.middleCols(_S, _E - _S);
	}
	//D.noalias() = left * mid;
	ss = nn / _mt / 2;
#pragma omp parallel for schedule(dynamic,4)
	for (int64_t ii = 0; ii < nn; ii += ss)
	{
		int64_t _S = ii;
		int64_t _E = ii + ss;
		if (_E >= nn)_E = nn;
		D->_dmat.middleCols(_S, _E - _S).noalias() = F * right.middleCols(_S, _E - _S);
	}
	D->_dmat = C->_dmat * D->_dmat;
	*/
	//Eigen::EigenSolver<Eigen::MatrixXd> eigen;
	//eigen.compute(D->_dmat);

	//singularvalues->__v.resize(D->_dmat.rows());
	//singularvalues->__v = eigen.eigenvalues().real();
}

void KingOfMonsters::_mySparse::ofAtB(_mySparse* B, bool sparse)
{
	static std::vector<std::vector<int64_t>> index;
	static std::map<_mySparse*, std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>>> ___e;
	static std::vector<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> e2;
	auto ss = std::stringstream();
	auto now = high_resolution_clock::now();
	int64_t nn = this->cols();
	int64_t mm = B->cols();
	int64_t kk = this->rows();
	eigen_assert(kk == B->rows());
	//int64_t _mt = omp_get_max_threads();
	//omp_set_num_threads(_mt);
	auto end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	int64_t __mt = _mt;// std::min(_nt / 10, _mt * 10);
	std::vector < Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>>* e;
	if (___e.contains(this))e = &___e[this]; else { std::vector < Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> _e; ___e[this] = _e; e = &___e[this]; }
	if (e->size() < __mt)
	{
		e->resize(__mt);
		e2.resize(__mt);
	}
	//Eigen::initParallel();
	//Eigen::setNbThreads(1);
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	ss << _nt << ":nt:" << std::endl;
	ss << _mt << ":_mt:" << std::endl;
	//Eigen::initParallel();
	int64_t job = 0;
	int64_t sss = _nt / __mt / 8;
	int64_t __nt2 = _nt;
	if (sss == 0)sss = 1;
	Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>* prevmat;
	if (dict2.contains(this)) {
		prevmat = &dict2[this];
	}
	else {
		Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> _prevmat(nn, mm);
		dict2[this] = _prevmat;
		prevmat = &dict2[this];
		prevmat->resize(nn, mm);
		prevmat->reserve(nn * mm / 10);
	}
	//prevmat->makeCompressed();
	std::vector<int64_t>* _map = 0;
	if (map.contains(prevmat))
	{
		_map = &map[prevmat];
	}
#pragma omp parallel for
	for (int64_t i = 0; i < __mt; i++) {
		if (_map == 0 || (*e)[i].nonZeros() != prevmat->nonZeros())
		{
			(*e)[i].resize(nn, mm);
			(*e)[i] = *prevmat;
			(*e)[i].makeCompressed();
		}
		memset((*e)[i].valuePtr(), 0, sizeof(double) * prevmat->nonZeros());
		e2[i].resize(nn, mm);
		e2[i].reserve(nn * mm / 100);
	}
	index.resize(__mt);
#pragma omp parallel for
	for (int64_t _ii = 0; _ii < __mt; _ii++)
	{
		int64_t S = 0;
		int64_t E = 0;
		//int64_t K = 0;
		//auto _e = e[_ii];
		for (int64_t tt = 0; tt < 4000; tt++)
		{
#pragma omp critical
			{
				S = job;
				E = job + sss;
				if (E > _nt)E = _nt;
				job = E;
			}
			if (S >= _nt)break;
			for (int64_t ii = S; ii < E; ii++)
			{
		
				//(*e)[_ii] += this->_mat[ii].transpose() * coeff[ii].asDiagonal() * this->_mat[ii];
//#pragma omp critical
				{

					e2[_ii] = this->_mat[ii].transpose() * coeff[ii].asDiagonal() * B->_mat[ii];
	

					//#pragma omp critical
					{
						if (_map == 0)
						{
							(*e)[_ii] += e2[_ii];
						}
						else {
							if (e2[_ii].nonZeros() > 0)
							{
								//int64_t* ptr = &index[_ii][0];
								for (int64_t k = 0; k < mm; ++k) {
									for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(e2[_ii], k); it; ++it) {
										//e[0].coeffRef(it.row(), it.col()) += it.value();
										*((*e)[_ii].valuePtr() + (*_map)[it.row() * mm + it.col()]) += it.value();
										//ptr++;
									}
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
	//Eigen::setNbThreads(_mt);

	for (int64_t tt = 0; tt < 4000; tt++)
	{
#pragma omp parallel for schedule(static,1)
		for (int64_t i = 0; i < __mt; i += 2)
		{
			if (i + 1 < __mt) {
				if (_map == 0 || (*e)[i].nonZeros() != (*e)[i + 1].nonZeros())
				{
					(*e)[i] += (*e)[i + 1];
				}
				else {
					Eigen::Map<Eigen::VectorXd> map1((*e)[i].valuePtr(), (*e)[i].nonZeros());
					Eigen::Map<Eigen::VectorXd> map2((*e)[i + 1].valuePtr(), (*e)[i].nonZeros());
					map1 += map2;
				}
			}
		}
		int64_t _ct = 0;
#pragma omp parallel for ordered schedule(static,1)
		for (int64_t i = 0; i < __mt; i += 2)
		{
#pragma omp ordered
			if (_map == 0 || (*e)[i].nonZeros() != (*e)[i / 2].nonZeros())
			{
				(*e)[i / 2] = (*e)[i];
			}
			else
			{
				memcpy((*e)[i / 2].valuePtr(), (*e)[i].valuePtr(), sizeof(double) * (*e)[i].nonZeros());
			}
#pragma omp atomic
			_ct++;
		}
		__mt = _ct;
		if (__mt == 1)break;
	}
	if (true) {
		if (this->_mat.size() == 0)this->_mat.resize(1);
		this->_mat[0].resize(nn, mm);
		this->_mat[0].reserve(nn * mm / 20);
		this->_mat[0] = (*e)[0];
		this->_mat[0].makeCompressed();
		if (_map == 0 || prevmat->nonZeros() != this->_mat[0].nonZeros())
		{
			*prevmat = this->_mat[0];
			prevmat->makeCompressed();
			dict2[this] = *prevmat;
			//build map
			std::vector<int64_t> __map;
			__map.resize(nn * mm);
			for (int64_t k = 0; k < prevmat->outerSize(); ++k) {
				for (Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>::InnerIterator it(*prevmat, k); it; ++it) {
					int64_t S = prevmat->outerIndexPtr()[k];
					int64_t E = prevmat->outerIndexPtr()[k + 1];

					for (int64_t tt = S; tt < E; tt++)
					{
						if (prevmat->innerIndexPtr()[tt] == it.row())
							__map[it.row() * mm + it.col()] = tt;
					}
				}
			}
			map[prevmat] = __map;
		}
	}

	double fff = this->_mat[0].sum();

	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	if (!sparse) {

	}
	end = high_resolution_clock::now();
	duration = duration_cast<milliseconds>(now - end);
	ss << duration.count() << "ms" << std::endl;
	now = high_resolution_clock::now();
	//this->_mat[0] = this->_dmat.sparseView(1.0, 0.0000000000001);
	ss << "sum:" << (*e)[0].sum() << std::endl;
	return;// ss.str();
}
void KingOfMonsters::_mySparse::_plus(int64_t i, int64_t j, double val)
{
	this->_mat[0].coeffRef(i, j) += val;
}
void KingOfMonsters::_mySparse::Atb(double* ptr, double* ptr2, double sc,int64_t N, Eigen::VectorXd* c)
{
	static std::map<int64_t, std::vector<Eigen::VectorXd>> __dict;
	c->resize(this->cols());
	c->setZero();
	int64_t nn = this->cols();// _mat[0].cols();
	int64_t offset = 0;
	std::vector<Eigen::VectorXd>* mm = 0;
	if (__dict.contains(nn))
	{
		mm = &__dict[nn];
	}
	else {
		std::vector<Eigen::VectorXd> _mm(_mt);
		for (int64_t i = 0; i < _mt; i++)
		{
			_mm[i].resize(nn);
		}
		__dict[nn] = _mm;
		mm= &__dict[nn];
	}
	int64_t job = 0;
	int64_t ss = _nt / _mt / 4;
	if (ss == 0)ss = 1;
#pragma omp parallel for
	for (int64_t ii = 0; ii < _mt; ii++)
	{

		int64_t S = 0;
		int64_t E = 0;
		//int64_t cpu =  omp_get_thread_num();
		auto tmp = (*mm)[ii];
		tmp.setZero();
		for (int64_t kk = 0; kk < 200; kk++)
		{
#pragma omp critical
			{
				S = job;
				E = S + ss;
				if (E > _nt)E = _nt;
				job = E;
			}
			if (S >= _nt)break;
			int64_t offset = 0;
			for (int64_t _t = 0; _t < S; _t++)
			{
				offset += coeff[_t].rows();
			}
			for (int64_t i = S; i < E; i++)
			{
				//Eigen::VectorXd tmp(this->cols());
				int64_t ee = coeff[i].rows();
		

				Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
				if (ptr2 != 0)
				{
					Eigen::Map<Eigen::VectorXd> b2(ptr2 + offset, ee);
					if (this->_mat[i].rows() > 0 && this->_mat[i].cols() > 0)
					{
						tmp += _mat[i].transpose() * coeff[i].asDiagonal() * (b + sc * b2);

					}
				}
				else {
					if (this->_mat[i].rows() > 0 && this->_mat[i].cols() > 0)
					{
						tmp += _mat[i].transpose() * coeff[i].asDiagonal() * (b);

					}
				}
				offset += ee;
			}
		}
#pragma omp critical
		{
			*c += tmp;// _mat[ii].transpose()* coeff[ii].asDiagonal()* (b + sc * b2);
		}
		//double sum = (*c).sum();
		//offset += ee;
	}
	//return ret;
 }
Eigen::VectorXd KingOfMonsters::_mySparse::Atb(double* ptr, int64_t N)
{
	Eigen::VectorXd ret(this->cols());
	ret.setZero();
	int64_t offset = 0;
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		int64_t ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{
			ret += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		}
		offset += ee;
	}
	return ret;
}
void KingOfMonsters::_mySparse::Atb(double* ptr, int64_t N,Eigen::VectorXd* c)
{
	//Eigen::VectorXd ret(this->cols());
	c->resize(this->cols());
	c->setZero();
	//ret.setZero();
	int64_t offset = 0;
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		int64_t ee = coeff[ii].rows();
		Eigen::Map<Eigen::VectorXd> b(ptr + offset, ee);
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)
		{
			*c += _mat[ii].transpose() * coeff[ii].asDiagonal() * b;
		}
		offset += ee;
	}
	//return ret;
}
Eigen::VectorXd KingOfMonsters::_mySparse::_Atb(double* ptr, int64_t N)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::Map<Eigen::VectorXd> b(ptr, N);
	return _dmat.transpose() * b;
}
void KingOfMonsters::_mySparse::merge()
{
	this->ofDat();
}
void KingOfMonsters::_mySparse::computeLU()
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	lu.compute(_dmat);
}
void KingOfMonsters::_mySparse::computeLLT(Eigen::LLT<Eigen::MatrixXd>* _LLT)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	_LLT->compute(_dmat);
}
int64_t KingOfMonsters::_mySparse::nonzeros() {
	int64_t _ret = 0;
	for (int64_t ii = 0; ii < _nt; ii++)
	{
		if (this->_mat[ii].rows() > 0 && this->_mat[ii].cols() > 0)

		_ret += _mat[ii].nonZeros();
	}
	return _ret;
}
void KingOfMonsters::_mySparse::solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	lu.compute(_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_dmat.cols());
	//Eigen::VectorXd x(_dmat.cols());
	//x.setZero();
	*ret = lu.solve(*rhs);
}
void KingOfMonsters::_mySparse::LSsolve(Eigen::VectorXd* rhs, Eigen::VectorXd* ret,double salt,int mode) {
#ifdef _CPU
	//Eigen::LLT<Eigen::MatrixXd> lu;
	//MKL_Set_Num_Threads(16);
	//MKL_Set_Dynamic(false);
	Eigen::PardisoLU< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	lu.pardisoParameterArray()[59] = 1;
	//Eigen::BiCGSTAB< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	//Eigen::MatrixXd m(this->_mat[0].rows(), this->_mat[0].cols());
	//m = this->_mat[0];
	//lu.setTolerance(DBL_EPSILON);
	//lu.setMaxIterations(this->_mat[0].rows() * 10);
	Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> id;
	if (mode==0)//minimum norm
	{
		id.resize(_mat[0].rows(), _mat[0].rows());
		id.setIdentity();
		lu.compute(this->_mat[0] * this->_mat[0].transpose()+ id* salt);
		ret->resize(this->_mat[0].cols());
		*ret = this->_mat[0].transpose() * lu.solve(*rhs);
		return;
	}
	if (mode==1)//least squares
	{
		id.resize(_mat[0].cols(), _mat[0].cols());
		id.setIdentity();
		lu.compute(this->_mat[0].transpose() * this->_mat[0] + id * salt);
		ret->resize(this->_mat[0].rows());
		Eigen::VectorXd vv = this->_mat[0].transpose() * *rhs;
		*ret = lu.solve(vv);
		return;
	}
#endif
}

void KingOfMonsters::_mySparse::Project(Eigen::VectorXd* rhs, Eigen::VectorXd* ret, double salt) {
	//Eigen::LLT<Eigen::MatrixXd> lu;
	//MKL_Set_Num_Threads(16);
	//MKL_Set_Dynamic(false);
	#ifdef _CPU
	Eigen::PardisoLU< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	lu.pardisoParameterArray()[59] = 1;
	//Eigen::MatrixXd m(this->_mat[0].rows(), this->_mat[0].cols());
	//m = this->_mat[0];
	Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> id;
	if (this->_mat[0].rows() <= this->_mat[0].cols())
	{
		id.resize(_mat[0].rows(), _mat[0].rows());
		lu.compute(this->_mat[0] * this->_mat[0].transpose() + id * salt);
		ret->resize(this->_mat[0].cols());
		Eigen::VectorXd v(this->_mat[0] * *rhs);
		*ret = *rhs-this->_mat[0].transpose() * lu.solve(v);
		return;
}
	if (this->_mat[0].rows() > this->_mat[0].cols())
	{
		id.resize(_mat[0].cols(), _mat[0].cols());
		lu.compute(this->_mat[0].transpose() * this->_mat[0] + id * salt) ;
		ret->resize(this->_mat[0].rows());
		Eigen::VectorXd v(((*rhs).transpose() * this->_mat[0]).transpose());
		*ret = lu.solve(v).transpose() * this->_mat[0].transpose();
		return;
	}
#endif
}


std::string KingOfMonsters::_mySparse::_solve0_gpu(KingOfMonsters::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int64_t device) {
	std::stringstream ss;
#ifdef _CPU
	Eigen::FullPivLU<Eigen::MatrixXd> lu(this->_dmat);
	*ret = lu.solve(*rhs);
	ss << "cpu_mode";
	return ss.str();
#else
	if (device < 0)
	{
		Eigen::FullPivLU<Eigen::MatrixXd> lu(this->_dmat);
		*ret = lu.solve(*rhs);
		ss << "cpu_mode";
		return ss.str();
	}
	else {
		//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
		Eigen::VectorXd x = *ret;
		//this->_freeze();
		int64_t N = rhs->rows();
		if (!cuda->valid())
		{
			x(0) = 10;
			return "";
		}
		cudaSetDevice(device);
		auto solver = cuda->solver(device, 0);
		//auto blas = cuda->blas(device);
		//auto stream=streams[device];
		//cudaStreamCreate(&stream);
		//cusolverDnSetStream(solver, stream);
		//Eigen::Map<Eigen::VectorXd> b(rhs, N);
		x.resize(N);
		//x.setZero();


		double* gpu_rhs = cuda->work_rhs(device);
		double* gpu_matrix = cuda->work_M(device);
		auto err = cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
		ss << "," << err;
		err = cudaMemcpy(gpu_rhs, rhs->data(), N * sizeof(double), cudaMemcpyHostToDevice);
		ss << "," << err;

		size_t work_size = 0;
		size_t work_size_host = 0;

		int* devInfo_on_gpu = cuda->info(device);
		//cudaMalloc(&devInfo_on_gpu, sizeof(int64_t));

		// --- CUDA CHOLESKY initialization
		cusolverDnParams_t params;
		cusolverDnCreateParams(&params);
		double* work_host;
		cusolverDnSetAdvOptions(params, cusolverDnFunction_t::CUSOLVERDN_POTRF, cusolverAlgMode_t::CUSOLVER_ALG_0);
		auto err2 = cusolverDnXpotrf_bufferSize(solver,params, CUBLAS_FILL_MODE_LOWER, N, cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F ,&work_size, &work_size_host);
		ss << "," << err2;
		if (work_size_host == 0)
		{
			work_host = 0;
		}
		else {
			work_host = (double*)malloc(sizeof(double) * work_size_host);
		}
		// --- CUDA POTRF execution
		double* work = 0;
		work = cuda->work(work_size, device);
		//cudaMalloc(&work, sizeof(double) * work_size);
		err = cudaMemset(work, 0, sizeof(double) * work_size);
		ss << "," << err;
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		ss << "," << err;
		//auto now = std::chrono::high_resolution_clock::now();
		err2 = cusolverDnXpotrf (solver,params, CUBLAS_FILL_MODE_LOWER, N,cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F, work, work_size, work_host,work_size_host,devInfo_on_gpu);

		ss << "," << err2;
		int devInfo_on_cpu;
		cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
		ss << "devInfo," << devInfo_on_cpu;

		//int64_t devInfo_on_cpu = 0;
		//cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int64_t), cudaMemcpyDeviceToHost);


		//if (0 != devInfo_on_cpu) {
		//	x(0) = devInfo_on_cpu;
		//	return;
		//}


		err2 = cusolverDnXpotrs(solver,params, CUBLAS_FILL_MODE_LOWER, N, 1, cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F, gpu_rhs, N, devInfo_on_gpu);
		ss << "," << err2;

		// auto end = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
		//std::cout << "Dn:" << duration.count() << "ms" << std::endl;
		devInfo_on_cpu;
		cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
		ss << "devInfo," << devInfo_on_cpu;
		//if (devInfo_on_cpu != 0) {
		//	x(0) = 24;
		//	return;
		//}

		cudaMemcpy(ret->data(), gpu_rhs, sizeof(double) * N, cudaMemcpyDeviceToHost);

		//cudaFree(work);

		//cudaFree(devInfo_on_gpu);
		cusolverDnDestroyParams(params);
		if (work_size_host == 0)
		{
		}
		else {
			free(work_host);
		}
		cudaDeviceSynchronize();
		ss << "size()" << ret->size();
		ss << "norm" << ret->norm();
		//cudaStreamDestroy(stream);
		return ss.str();
}
#endif
}
double KingOfMonsters::_mySparse::computeeigen(_mySparse* i1, _mySparse* i2, _mySparse* f1, int64_t N)
{
	auto xX = i1->_dmat * f1->_mat[0];
	auto Xx = i2->_dmat * f1->_mat[0].transpose();
	auto xXXx = xX * Xx;
	Eigen::MatrixXd mat;
	mat = xXXx;

	int64_t D = xXXx.cols();
	static Eigen::VectorXd vec(D);

	vec(0) = 1;
	vec.normalize();

	for (int64_t i = 0; i < N; i++)
	{
		vec = xXXx * vec;
	}
	vec.normalize();
	vec = xXXx * vec;
	return vec.norm();
	//return eigen.eigenvalues();

}
void  KingOfMonsters::_mySparse::project(_mySparse* i1/*JxxJ*/, _mySparse* i2/*JXXJ*/, _mySparse* f1/*JxJX*/,_myDoubleArray* v1,_myDoubleArray* v2, _myDoubleArray* _ret1, _myDoubleArray* _ret2)
{
	Eigen::MatrixXd JxxJ = i1->_mat[0];
	Eigen::MatrixXd JXXJ = i2->_mat[0];
	auto JxXJ = f1->_mat[0];
	JxxJ += 0.000001 * Eigen::MatrixXd::Identity(v2->__v.size(), v2->__v.size());
	JXXJ += 0.000001 * Eigen::MatrixXd::Identity(v1->__v.size(), v1->__v.size());
	Eigen::MatrixXd xX = JxxJ.inverse() * JxXJ;
	Eigen::MatrixXd Xx = JXXJ.inverse() * JxXJ.transpose();

	Eigen::MatrixXd xXXx = xX * Xx;
	Eigen::MatrixXd XxxX = Xx * xX;

	auto I = Eigen::MatrixXd::Identity(v1->__v.size(), v1->__v.size());
	double lambda = 0;
	Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr;
	for (int i = 0; i < 6; i++)
	{
		qr.compute((XxxX - I).transpose() * (XxxX - I) + lambda * I);
		double norm1 = v1->__v.transpose() * qr.solve(qr.solve(v1->__v));
		double norm2 = v1->__v.transpose() * qr.solve(v1->__v);
		lambda = norm2 / norm1;
	}
	qr.compute((XxxX - I).transpose() * (XxxX - I) + lambda * I);
	Eigen::VectorXd ret = qr.solve(v1->__v);// +(xX).transpose() * v2->__v);
	//ret.normalize();
	ret *= lambda;// lambda;
	Eigen::VectorXd ret2 = xX*ret;
	
	_ret1->__v = ret;
	_ret2->__v = ret2;
	
}
/*std::string KingOfMonsters::_mySparse::_QR_cpu(Eigen::MatrixXd* Q, Eigen::MatrixXd* R) {

	std::stringstream ss;

	//dmat is a column-major matrix
	//CUDA also uses a column-major format
	//this->_freeze();
	int64_t M = this->_dmat.rows();
	int64_t N = this->_dmat.cols();


	Q->resize(M, N);
	R->resize(N, N);
	Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr;
	qr.compute(this->_dmat);
	Eigen::MatrixXd thinQ = qr.householderQ().setLength(qr.nonzeroPivots());
	*Q = qr.householderQ();
	*R = qr.matrixR();

	
	return ss.str();
}*/
/*int KingOfMonsters::_mySparse::_QR_cpu(Eigen::MatrixXd* Q, Eigen::MatrixXd* R) {

	int64_t M = this->_dmat.rows();
	int64_t N = this->_dmat.cols();
	Eigen::FullPivHouseholderQR<Eigen::MatrixXd> qr(this->_dmat);

	*Q = qr.matrixQ();
	return qr.rank();
}*/
std::string KingOfMonsters::_mySparse::_solveLU_dense_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret)
{
	this->_mat[0].makeCompressed();
	//Eigen::PartialPivLU<Eigen::MatrixXd> lu(this->_dmat);
	//Eigen::FullPivHouseholderQR< Eigen::MatrixXd> lu;
	Eigen::SparseQR< Eigen::SparseMatrix<double,0,int64_t>,Eigen::COLAMDOrdering<int64_t>> lu;

	//Eigen::JacobiSVD<Eigen::MatrixXd> lu;
	lu.setPivotThreshold(0.0000000001);
	//lu.setThreshold(0.000000000001);
	//lu.compute(this->_dmat);
	lu.compute(this->_mat[0]);
	*ret=lu.solve(*rhs);
	//if (lu.info() == Eigen::ComputationInfo::Success)
	{

		//*ret = lu.solve(*rhs);
		return "success";
	}
	/*
	else {
		std::stringstream ss;
		ss << lu.info();
		return ss.str();
			;
	}*/
}
void KingOfMonsters::_mySparse::turnDense()
{
	this->_dmat = this->_mat[0];
}
std::string KingOfMonsters::_mySparse::_solveLU_sparse_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret)
{
	//MKL_Set_Num_Threads(16);
	//MKL_Set_Dynamic(false);
#ifdef _CPU
	_mt = MKL_Get_Max_Threads();
	if(_mt>1)
	MKL_Set_Num_Threads(_mt - 1);

	Eigen::setNbThreads(_mt);
	
	this->_mat[0].makeCompressed();
	//Eigen::SparseLU< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	//Eigen::SparseQR< Eigen::SparseMatrix<double, 0, int64_t>, Eigen::COLAMDOrdering<int64_t>>lu;
	//Eigen::BiCGSTAB< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	
	MKL_Set_Dynamic(true);
	
	Eigen::PardisoLU < Eigen::SparseMatrix<double, 0, int64_t>> lu;
	//lu.pardisoParameterArray()[59] = 1;
	//pardiso.compute(this->_mat[0]);


	//lu.setPivotThreshold(0.0000000001);l;ll;
	//lu.setMaxIterations(rhs->size() * 0.5);
	lu.compute(this->_mat[0]);
	if (lu.info() == Eigen::ComputationInfo::Success)
	{
		if (rhs->size() > this->_mat[0].rows())
		{
			Eigen::MatrixXd _rhs(this->_mat[0].rows(), rhs->size() / this->_mat[0].rows());
			for (int i = 0; i < rhs->size(); i += this->_mat[0].rows())
			{
				_rhs.col(i/this->_mat[0].rows())= rhs->segment(i, this->_mat[0].rows());
			}
			Eigen::MatrixXd _ret = lu.solve(_rhs);
			ret->resize(rhs->size());
			ret->setZero();
			for (int i = 0; i < rhs->size(); i += this->_mat[0].rows())
			{
				ret->segment(i, this->_mat[0].rows())=_ret.col(i / this->_mat[0].rows());
			}
		}
		else {
			*ret = lu.solve(*rhs);
		}
		return "success";
	}
	else {
		std::stringstream ss;
		ss << lu.info();
		return ss.str();
		;
	}
#endif
	return "_cpu only";
}
std::string KingOfMonsters::_mySparse::_solveCG_sparse_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret)
{
	//MKL_Set_Num_Threads(16);
	//MKL_Set_Dynamic(false);
	this->_mat[0].makeCompressed();
	//Eigen::SparseLU< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	//Eigen::SparseQR< Eigen::SparseMatrix<double, 0, int64_t>, Eigen::COLAMDOrdering<int64_t>>lu;
	Eigen::BiCGSTAB< Eigen::SparseMatrix<double, 0, int64_t>> lu;
	//MKL_Set_Num_Threads(16);
	//MKL_Set_Dynamic(false);
	//Eigen::PardisoLU < Eigen::SparseMatrix<double, 0, int64_t>> lu;
	//lu.pardisoParameterArray()[59] = 1;
	//pardiso.compute(this->_mat[0]);


	//lu.setPivotThreshold(0.0000000001);l;ll;
	lu.setMaxIterations(rhs->size() * 0.2);
	lu.compute(this->_mat[0]);
	if (lu.info() == Eigen::ComputationInfo::Success)
	{

		*ret = lu.solve(*rhs);
		return "success";
	}
	else {
		std::stringstream ss;
		ss << lu.info();
		return ss.str();
		;
	}

}
std::string KingOfMonsters::_mySparse::_solveLU_gpu(KingOfMonsters::cuda* cuda, Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int64_t device) {
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	std::stringstream ss;

#ifdef _CPU
	Eigen::FullPivLU<Eigen::MatrixXd> lu(this->_dmat);
	*ret = lu.solve(*rhs);
	ss << "cpu_mode";
	return ss.str();
#else
	if (device < 0)
	{
		Eigen::FullPivLU<Eigen::MatrixXd> lu(this->_dmat);
		*ret = lu.solve(*rhs);
		ss << "cpu_mode";
		return ss.str();
	}
	else
	{
		Eigen::VectorXd x = *ret;
		//this->_freeze();
		int64_t N = rhs->rows();
		if (!cuda->valid())
		{
			x(0) = 10;
			return "";
		}
		cudaSetDevice(device);
		auto solver = cuda->solver(device, 0);
		//auto blas = cuda->blas(device);
		//auto stream=streams[device];
		//cudaStreamCreate(&stream);
		//cusolverDnSetStream(solver, stream);
		//Eigen::Map<Eigen::VectorXd> b(rhs, N);
		x.resize(N);
		//x.setZero();
		cusolverDnParams_t params;
		cusolverDnCreateParams(&params);
		cusolverDnSetAdvOptions(params, cusolverDnFunction_t::CUSOLVERDN_GETRF, cusolverAlgMode_t::CUSOLVER_ALG_0);
		double* gpu_rhs = cuda->work_rhs(device);
		double* gpu_matrix = cuda->work_M(device);
		auto err = cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
		ss << "," << err;
		err = cudaMemcpy(gpu_rhs, rhs->data(), N * sizeof(double), cudaMemcpyHostToDevice);
		ss << "," << err;

		size_t work_size = 0;
		size_t work_size_host = 0;
		int* devInfo_on_gpu = cuda->info(device);
		//cudaMalloc(&devInfo_on_gpu, sizeof(int64_t));
		double* work_host;
		// --- CUDA CHOLESKY initialization
		auto err2 = cusolverDnXgetrf_bufferSize(solver, params, N, N ,cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F, &work_size,&work_size_host);
		ss << "," << err2;
		if (work_size_host == 0)
		{
			work_host = 0;
		}
		else {
			work_host = (double*)malloc(sizeof(double) * work_size_host);
		}
		// --- CUDA POTRF execution
		double* work = 0;
		work = cuda->work(work_size, device);
		//cudaMalloc(&work, sizeof(double) * work_size);
		err = cudaMemset(work, 0, sizeof(double) * work_size);
		ss << "," << err;
		cudaDeviceSynchronize();
		err = cudaGetLastError();
		ss << "," << err;
		//auto now = std::chrono::high_resolution_clock::now();
		int64_t *ipiv;
		cudaMalloc(&ipiv, sizeof(int64_t)*N);
		err2 = cusolverDnXgetrf(solver, params,N, N, cudaDataType::CUDA_R_64F, gpu_matrix ,N,ipiv,  cudaDataType::CUDA_R_64F, work, work_size,work_host,work_size_host, devInfo_on_gpu);
		ss << "," << err2;
		int devInfo_on_cpu;
		cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
		ss << "devInfo," << devInfo_on_cpu;

		//int64_t devInfo_on_cpu = 0;
		//cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int64_t), cudaMemcpyDeviceToHost);


		//if (0 != devInfo_on_cpu) {
		//	x(0) = devInfo_on_cpu;
		//	return;
		//}


		err2 = cusolverDnXgetrs(solver, params,cublasOperation_t::CUBLAS_OP_N, N, 1, cudaDataType::CUDA_R_64F, gpu_matrix, N,ipiv, cudaDataType::CUDA_R_64F, gpu_rhs, N, devInfo_on_gpu);
		ss << "," << err2;

		// auto end = std::chrono::high_resolution_clock::now();
		//auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - now);
		//std::cout << "Dn:" << duration.count() << "ms" << std::endl;
		cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);
		ss << "devInfo," << devInfo_on_cpu;
		//if (devInfo_on_cpu != 0) {
		//	x(0) = 24;
		//	return;
		//}

		cudaMemcpy(ret->data(), gpu_rhs, sizeof(double) * N, cudaMemcpyDeviceToHost);

		//cudaFree(work);

		//cudaFree(devInfo_on_gpu);
		
		if (work_size_host == 0) {}
		else {

			free(work_host);

		}
		cudaFree(ipiv);
		cusolverDnDestroyParams(params);
		cudaDeviceSynchronize();
		ss << "size()" << ret->size();
		ss << "norm" << ret->norm();
		//cudaStreamDestroy(stream);
		return ss.str();
	}
#endif
}
/*Eigen::MatrixXd KingOfMonsters::_mySparse::_solve0(_myLLT* LLT, _mySparse* mat)
{
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//this function assumes that LLT decomposition has been done already
	Eigen::MatrixXd ret(this->_dmat.cols(), mat->_dmat.cols());
	int64_t nn = mat->cols();
	int64_t numthreads = _mt;
	int64_t S = nn / numthreads / 2;

	//ret = LLT->LLT->solve(mat->_dmat);

#pragma omp parallel for schedule(dynamic,1)
	for (int64_t i = 0; i < nn; i += S)
	{
		int64_t start = i;
		int64_t end = i + S;
		if (end > nn)end = nn;
		ret.middleRows(start, end - start) = LLT->LLT->solve(_dmat.middleCols(start, end - start)).transpose();

	}
	return ret.transpose();
}*/


int64_t KingOfMonsters::_mySparse::_solveI(_mySparse* ret)
{
	//ret->_dmat = this->_dmat.inverse();
	//return 0;
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> llt;
	llt.compute(_mat[0]);
	auto I = Eigen::SparseMatrix<double>(_mat[0].cols(), _mat[0].cols());
	I.setIdentity();
	ret->_dmat = llt.solve(I);
	return 0;
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);	
	//Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::ColMajor> llt;
	//Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>> llt;
	//Eigen::FullPivLU<
}
int64_t KingOfMonsters::_mySparse::_solveI_dense(_mySparse* ret)
{
	//ret->_dmat = this->_dmat.inverse();
	//return 0;
	ret->_dmat = this->_dmat.inverse();
	return 0;
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);	
	//Eigen::SimplicialLLT<Eigen::SparseMatrix<double>, Eigen::ColMajor> llt;
	//Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>> llt;
	//Eigen::FullPivLU<
}

void initidentiy(KingOfMonsters::cuda* cuda, int64_t N,bool mp) {
	double _a = 0;
	double _b = 1;
#ifndef _CPU
	//if (I.cols() != N || I.rows() != N)
	{
		//I.setIdentity(N, N);
#pragma omp parallel for
		for (int64_t ii = 0; ii < cuda->count(); ii++)
		{
			int64_t S = N * ii / cuda->count();
			int64_t E = N * (ii+1) / cuda->count();
			if (!mp)S = 0;
			if (!mp)E = N;

			cudaSetDevice(ii);

			double* gpu_rhs = cuda->work_C(ii);
			auto stream = cuda->__streams(ii, 0);
			//cudaMemcpyAsync(gpu_rhs, I.data(), sizeof(double) * N * N, cudaMemcpyHostToDevice, stream);
			cudaMemsetAsync(gpu_rhs, 0, N * (E-S) * sizeof(double), stream);
			for (int64_t i = S; i < E; i++)
			{
				cudaMemcpyAsync(gpu_rhs + (i + (i-S) * N), &_b, 1 * sizeof(double), cudaMemcpyHostToDevice, stream);
			}
		}
	}
#endif
}
std::string KingOfMonsters::_mySparse::_solveI_gpu_single(KingOfMonsters::cuda* cuda, _mySparse* ret)
{	std::stringstream sss;

#ifdef _CPU
MKL_Set_Num_Threads(16);
MKL_Set_Dynamic(false);
ret->_dmat = this->_dmat.inverse();
	sss << "cpu_mode";
	return sss.str();
#else
	this->_freeze();
	int64_t N = _dmat.cols();// __c;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

	//ret->__r = N;
	//ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int64_t nn = cuda->count();
	initidentiy(cuda, N,false);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret->_dmat.setZero(N, N);
	//ret->_tmp.setZero(N, N);




	auto solver = cuda->solver(cuda->fastest(), 0);

	cudaSetDevice(cuda->fastest());
	//auto stream = streams[cuda->fastest()];
	double* gpu_matrix = cuda->work_M(cuda->fastest());
	cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	int* devInfo_on_gpu = cuda->info(cuda->fastest());

	size_t work_size = 0;
	size_t work_size_host = 0;
	cusolverDnParams_t params;
	cusolverDnCreateParams(&params);
	cusolverDnSetAdvOptions(params, cusolverDnFunction_t::CUSOLVERDN_POTRF, cusolverAlgMode_t::CUSOLVER_ALG_0);
	// --- CUDA CHOLESKY initialization
	auto err2 = cusolverDnXpotrf_bufferSize(solver, params,CUBLAS_FILL_MODE_LOWER , N, cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F, &work_size, &work_size_host);
	//ss << "," << err2;
	// --- CUDA POTRF execution
	double* work = cuda->work(work_size, cuda->fastest());
	double* work_host;
	if (work_size_host == 0)
	{
		work_host = 0;
	}
	else {
		work_host = (double*)malloc(sizeof(double) * work_size_host);
	}
	
	cudaMemset(work, 0, work_size * sizeof(double));
	//cusolverDnSetStream(solver, stream);
	//int64_t* ipiv;
	//cudaMalloc(&ipiv, sizeof(int64_t) * N);
	auto err=cusolverDnXpotrf(solver, params, CUBLAS_FILL_MODE_LOWER, N, cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F, work, work_size, work_host, work_size_host, devInfo_on_gpu);


	int devInfo_on_cpu = 0;

	double* gpu_rhs = cuda->work_C(cuda->fastest());
	err2 = cusolverDnXpotrs(solver, params, CUBLAS_FILL_MODE_LOWER, N, 1,cudaDataType::CUDA_R_64F, gpu_matrix, N,  cudaDataType::CUDA_R_64F, gpu_rhs, N, devInfo_on_gpu);
	cudaMemcpyAsync(ret->_dmat.data() , gpu_rhs, N * N * sizeof(double), cudaMemcpyDeviceToHost);
	
	
	//int64_t _S = 0;
	//int64_t _E = 0;


	bool exit = false;

	
	//cudaFree(ipiv);
	if (work_size_host == 0)
	{	
	}
	else {
		free(work_host);
	}
	cusolverDnDestroyParams(params);
	cudaDeviceSynchronize();
	return sss.str();
#endif
}
std::string KingOfMonsters::_mySparse::_solveI_gpu_omp(KingOfMonsters::cuda* cuda, _mySparse* ret)
{
	std::stringstream sss;

#ifdef _CPU
	ret->_dmat = this->_dmat.inverse();
	sss << "cpu_mode";
	return sss.str();
#else
	this->_freeze();
	int64_t N = this->_dmat.cols();

	//ret->__r = N;
	//ret->__c = N;
	//ret->_dmat.resize(N, N);

	if (!cuda->valid())return "";
	int64_t nn = cuda->count();
	initidentiy(cuda, N,true);
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
		for (int64_t ii = 0; ii < nn; ii++)
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
		for (int64_t ss = 0; ss < STREAMCOUNT; ss++)
		{
			cudaSetDevice(cuda->fastest());
			int64_t S = ss * N / STREAMCOUNT;
			int64_t E = (ss + 1) * N / STREAMCOUNT;
			cudaStream_t _streamX = cuda->__streams(cuda->fastest(), ss);
			cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_matrix + S * N, sizeof(double) * N * (E - S), cudaMemcpyDeviceToHost, _streamX);
			cudaStreamSynchronize(_streamX);
			//cudaDeviceSynchronize();

#pragma omp parallel for
			for (int64_t ii = 0; ii < nn; ii++)
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

	int64_t job = 0;
	int64_t ss = N / cuda->count() / STREAMCOUNT / 4;
	if (ss == 0) ss = 1;
	for (int64_t i = 0; i < nn; i++)
	{
		cudaSetDevice(i);
		cudaDeviceSynchronize();
	}
	std::vector<int64_t>count(nn);
#pragma omp parallel for
	for (int64_t ii = 0; ii < nn; ii++)
	{
		int64_t S = N * ii / nn;
		int64_t E = N * (ii+1) / nn;

		auto solver = cuda->solver(ii, 0);

		cudaSetDevice(ii);
		//auto stream=streams[ii];
		//cusolverDnSetStream(solver, stream);

		double* gpu_matrix = cuda->work_M(ii);
		double* gpu_rhs = cuda->work_C(ii);
		int* devInfo_on_gpu = cuda->info(ii);

	    int devInfo_on_cpu = 0;

		bool exit = false;
		//while (true)
		//{
//#pragma omp parallel for
			//for (int64_t kk = 0; kk < STREAMCOUNT; kk++)
			//{
				//int64_t S = 0;
				//int64_t E = 0;

//#pragma omp critical
//				{
//					S = job;
//					E = S + ss;
//					if (E > N)E = N;
//					job = E;
//				}
//				if (S >= N) {
//					exit = true;
//					break;
//				}
//				count[ii]++;
				//int64_t S = _S + (_E - _S) * kk / 4;
				//int64_t E = _S + (_E - _S) * (kk+1) / 4;
				cudaStream_t _stream = cuda->__streams(ii, 0);
				auto _solver = cuda->solver(ii, 0);
				cusolverDnSetStream(_solver, _stream);

				cusolverDnDpotrs(_solver, CUBLAS_FILL_MODE_LOWER, N, E - S, gpu_matrix, N, gpu_rhs, N, devInfo_on_gpu);
				cudaMemcpyAsync(ret->_dmat.data() + S * N, gpu_rhs, (E - S) * N * sizeof(double), cudaMemcpyDeviceToHost, _stream);
			//}
			//if (exit)break;
		//}

		//cusolverDnDestroy(_solver);
		//cusolverDnSetStream(solver,streams[ii]);
	}
	cudaDeviceSynchronize();
	sss << "jobs taken: ";
	for (int64_t i = 0; i < nn; i++)
	{
		sss << i << "-" << count[i] << "  ";
	}
	sss << std::endl;
	return sss.str();
#endif
}

std::string KingOfMonsters::_mySparse::_solveI_gpu_sparse(KingOfMonsters::cuda* cuda, _mySparse* ret)
{
	ret->_dmat = _dmat.inverse();
	return "";

}
/*
std::string KingOfMonsters::_mySparse::_solveI_gpu(KingOfMonsters::cuda* cuda, _mySparse* ret)
{
	std::stringstream ss;

#ifdef _CPU
	ret->_dmat = this->_dmat.inverse();
	ss << "cpu_mode";
	return ss.str();
#else
	//this->_freeze();
	int64_t N = _dmat.cols();// __c;


	int64_t count = 0;
	for (int64_t i = 0; i < N; i++)
	{
		if (this->_mat[0].coeff(i, i) < 0)count++;
	}
	ss << "negative diagonal count=" << count << ",";


	//ret->__r = N;
	//ret->__c = N;
	//ret->_dmat.resize(N, N);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, N, N);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, N);
	ret->_dmat.setZero(N, N);
	//ret->_tmp.setZero(N, N);
	//initidentiy(cuda, N);
	if (!cuda->valid())return ss.str();

	int64_t ii = cuda->fastest();
	auto solver = cuda->solver(ii, 0);
	//cudaStream_t stream = streams[ii];
	cusolverDnSetStream(solver, cuda->__streams(cuda->fastest(), 0));
	cudaSetDevice(ii);
	double* m_gpu = cuda->work_M(ii);
	cusolverDnParams_t params;
	cusolverDnCreateParams(&params);
	auto err=cudaMemcpyAsync(m_gpu, _dmat.data(), sizeof(double) * N * N, cudaMemcpyHostToDevice, cuda->__streams(cuda->fastest(), 0));
	ss << "," << err;
	double* work;
	//int work_size;
	size_t work_size = 0;
	size_t work_size_host = 0;

	//auto err2=cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, &work_size1);
	auto err2 = cusolverDnXpotrf_bufferSize(solver, params, CUBLAS_FILL_MODE_LOWER, N, cudaDataType::CUDA_R_64F, gpu_matrix, N, cudaDataType::CUDA_R_64F, &work_size, &work_size_host);
	ss << "," << err2;

	double* work_host;
	work_host = (double*)malloc(sizeof(double) * work_size_host);
	cusolverDnSetAdvOptions(params, cusolverDnFunction_t::CUSOLVERDN_POTRF, cusolverAlgMode_t::CUSOLVER_ALG_0);

	//work_size = std::max(work_size1, work_size2);
	work = cuda->work(work_size, ii, cuda->__streams(cuda->fastest(), 0));
	//cudaMalloc(&work, work_size1 * sizeof(double));
	//cudaMallocAsync(&work, sizeof(double) * work_size, stream);
	cudaDeviceSynchronize();
	err=cudaMemset(work, 0, sizeof(double) * work_size1);
	ss << "," << err;
	cudaDeviceSynchronize();
	err = cudaGetLastError();
	ss << "," << err;
	ss << "," << work_size << "," << work_size1 << "," << work_size2;
	int* devInfo = cuda->info(ii);

	err2=cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work, work_size1, devInfo);
	ss << "," << err2;
	double* work2 = cuda->work_rhs(ii);
	//cudaMalloc(&work2, sizeof(double) * work_size2);
	err=cudaMemset(work2, 0, sizeof(double) * work_size2);
	ss << "," << err;
	err2=cusolverDnDpotri(solver, CUBLAS_FILL_MODE_LOWER, N, m_gpu, N, work2, work_size2, devInfo);
	ss << "," << err2;
	//cudaFree(work2);
	cusolverDnDestroyParams(params);
	err=cudaMemcpy(ret->_dmat.data(), m_gpu, sizeof(double) * N * N, cudaMemcpyDeviceToHost);
	ss << "," << err2;
	cudaDeviceSynchronize();
	//cudaFree(work);
	ss << "sum" << ret->_dmat.sum();
	ret->_dmat.triangularView<Eigen::Upper>() = ret->_dmat.triangularView<Eigen::Lower>().transpose();
	ss << "sum2" << ret->_dmat.sum();
	return ss.str();
#endif
}
*/

void KingOfMonsters::_mySparse::_solve0_gpu(KingOfMonsters::cuda* cuda, _mySparse* mat, _mySparse* ret)
{
#ifdef _CPU
	//MKL_Set_Num_Threads(16);
	//MKL_Set_Dynamic(false);
	Eigen::FullPivLU<Eigen::MatrixXd> lu(this->_dmat);
	ret->_dmat=lu.solve(mat->_dmat);
#else
	int64_t nn = mat->_dmat.cols();// __c;
	int64_t N = this->_dmat.cols();// __c;
	//Eigen::MatrixXd x(N,nn);
	
	//ret->__c = nn;
	//ret->__r = N;
	if (!cuda->valid())return;
	auto solver = cuda->solver(cuda->fastest(), 0);
	//auto blas = cuda->blas(cuda->fastest());
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> ret_dmat(ret->___dmat, N, nn);
	//Eigen::Map<Eigen::MatrixXd> mat_dmat(mat->___dmat, mat->__r, mat->__c);
	ret->_dmat.setZero(N, nn);
	//ret->_tmp.setZero(N, N);
	int64_t job = 0;

	cudaSetDevice(cuda->fastest());
	double* gpu_matrix = cuda->work_M(cuda->fastest());
	cudaMemcpy(gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
	int* devInfo_on_gpu;
	cudaMalloc(&devInfo_on_gpu, sizeof(int64_t));
	double* work;
	int work_size = 0;

	// --- CUDA CHOLESKY initialization
	cusolverDnDpotrf_bufferSize(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, &work_size);

	// --- CUDA POTRF execution	
	work = cuda->work(work_size, cuda->fastest());
	//cudaMalloc(&work, work_size * sizeof(double));
	auto now = std::chrono::high_resolution_clock::now();
	cusolverDnDpotrf(solver, CUBLAS_FILL_MODE_LOWER, N, gpu_matrix, N, work, work_size, devInfo_on_gpu);
	int  devInfo_on_cpu = 0;
	cudaMemcpy(&devInfo_on_cpu, devInfo_on_gpu, sizeof(int), cudaMemcpyDeviceToHost);

	if (0 != devInfo_on_cpu) {
		ret->_dmat(0, 0) = 2;
		return;
	}
	cudaMemcpy(_dmat.data(), gpu_matrix, N * N * sizeof(double), cudaMemcpyDeviceToHost);
	//cudaFree(work);
	cudaFree(devInfo_on_gpu);
	bool exit = false;

	int64_t S = nn / cuda->count() / 5;

#pragma omp parallel for
	for (int64_t i = 0; i < cuda->count(); i++) {
		cudaSetDevice(i);
		auto solver = cuda->solver(i, 0);
		//auto blas = cuda->blas(i);
		double* _gpu_matrix = cuda->work_M(i);
		double* gpu_rhs = cuda->work_rhs(i);
		cudaMemcpy(_gpu_matrix, _dmat.data(), N * N * sizeof(double), cudaMemcpyHostToDevice);
		cudaMemcpy(gpu_rhs, mat->_dmat.data(), N * nn * sizeof(double), cudaMemcpyHostToDevice);
		int* _devInfo_on_gpu = 0;
		int _devInfo_on_cpu = 0;
		cudaMalloc(&_devInfo_on_gpu, sizeof(int));
		while (true)
		{
			int64_t nextjob = -1;
#pragma omp critical
			{
				nextjob = job;
				job += S;
			}
			if (nextjob >= nn)break;
			int64_t start = nextjob;
			int64_t end = nextjob + S;
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
#endif
}


std::string KingOfMonsters::_mySparse::_solve0_lu_cpu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret, int ordering) {
	Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor, int64_t>,Eigen::COLAMDOrdering<int64_t>> lu;
	//Eigen::SimplicialLDLT< Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>, Eigen::Lower, Eigen::COLAMDOrdering<int64_t>> lu;
	lu.compute(_mat[0]);

	if (lu.info()== Eigen::ComputationInfo::Success) {

		ret->conservativeResize(_mat[0].cols());
		ret->setZero();
		*ret = lu.solve(*rhs);
		//solver->solve(rhs->data(), ret->data());
		return "SUCCESS";
	}
	else {

		return ("FAILED ");

	}

}
void KingOfMonsters::_mySparse::solve0_lu(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
	//Eigen::SparseLU<Eigen::SparseMatrix<double,Eigen::ColMajor>> lu;
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	lu.compute(this->_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_mat[0].cols());
	ret->setZero();
	//Eigen::VectorXd x(_mat[0].rows());
	//x.setZero();
	*ret = lu.solve(*rhs);
	//return x;
}
void KingOfMonsters::_mySparse::_solve0_lu_cg(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
	//Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>> lu;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>, Eigen::Lower | Eigen::Upper> cg;
	cg.compute(_mat[0]);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->conservativeResize(_mat[0].cols());
	ret->setZero();
	//Eigen::VectorXd x(_mat[0].rows());
	//x.setZero();
	*ret = cg.solve(*rhs);
	//return x;
}
void KingOfMonsters::_mySparse::__solve0(Eigen::VectorXd* rhs, Eigen::VectorXd* ret) {
	Eigen::LLT<Eigen::MatrixXd> LLT;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

	LLT.compute(_dmat);
	//Eigen::Map<Eigen::VectorXd> b(rhs, N);
	ret->resize(_dmat.cols());
	ret->setZero();
	*ret = LLT.solve(*rhs);
}
Eigen::MatrixXd KingOfMonsters::_mySparse::inv() {
	Eigen::PartialPivLU<Eigen::MatrixXd> lu;
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	lu.compute(_dmat);
	Eigen::MatrixXd I(_dmat.rows(), _dmat.cols());
	I.setIdentity();
	return lu.solve(I);
}

Eigen::MatrixXd KingOfMonsters::_mySparse::solve0(_mySparse* rhs)
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


void KingOfMonsters::_mySparse::minus(_mySparse* m) {
	this->_freeze();
	//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);
	//Eigen::Map<Eigen::MatrixXd> m_dmat(m->___dmat, m->__r, m->__c);

	_dmat = _dmat - m->_dmat;
	//_mat[0] = _dmat.sparseView(1.0, 0.00000000001);
}
void KingOfMonsters::_mySparse::clearcoeff() {
	for (int64_t ii = 0; ii < _nt; ii++)
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
Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t> id;
void KingOfMonsters::_mySparse::addsmallidentity(double salt, bool sparse, bool dense,int m) {

	if (dense)
	{
		id.resize(_dmat.rows(), _dmat.cols());
		id.setZero();
		for (int i = 0; i < m; i++)
			id.coeffRef(i, i) = 1;
		//id.setIdentity();
		//Eigen::Map<Eigen::MatrixXd> _dmat(___dmat, __r, __c);

		_dmat += (id * salt);
	}
	if (sparse)

	{
		if (this->_mat.size() >= 1)
		{
			id.resize(this->_mat[0].rows(), this->_mat[0].cols());
			id.setZero();
			for (int i = 0; i < m; i++)
				id.coeffRef(i, i) = 1;
			this->_mat[0] += (id * salt);
		}
	}
}
void KingOfMonsters::_mySparse::addsmallidentity(double salt, bool sparse, bool dense) {

	if (dense)
	{
		id.resize(_dmat.rows(), _dmat.cols());
		
		
		id.setIdentity();


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
