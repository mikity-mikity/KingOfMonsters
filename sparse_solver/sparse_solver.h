#pragma once

#include <chrono>
#include <vector>
#include <map>
#include <string>
using namespace std::chrono;
using std::vector;
using std::string;
#include "_sparse_solver.h"
using namespace System;

namespace sparse_solver {
	public ref class sparse_solver abstract sealed
	{
	public:
		static array<double>^ findSearchDirection(int n, array<int>^ columnptr, array<int>^ rowindices, array<double>^ values, array<array<double>^>^ vecs,int vecs_count,array<double>^rhs,double salt)
		{
			Eigen::SparseMatrix<double> mat;
			std::vector<Eigen::VectorXd> _vecs(vecs_count);
			for (int i = 0; i < vecs_count; i++)
			{
				Eigen::VectorXd v(n);
				for (int j = 0; j < n; j++)
				{
					double vv = vecs[i][j];
					v(j) = vv;
				}
				_vecs[i] = v;
			}
			std::vector<Eigen::Triplet<double, int64_t>> tripletList;
			tripletList.reserve(values->Length);
			for (int i = 0; i < n; i++)
			{
				int start = columnptr[i];
				int end = columnptr[i + 1];
				for (int j = start; j < end; j++)
				{
					int row = rowindices[j];
					double val = values[j];
					tripletList.push_back(Eigen::Triplet<double, int64_t>(row, i, val));
					
				}
			}

			mat.resize(n, n);
			mat.setFromTriplets(tripletList.begin(), tripletList.end());
			mat.makeCompressed();
			Eigen::VectorXd _rhs(n);
			for (int i = 0; i < n; i++)
			{
				double v = rhs[i];
				_rhs(i) = v;
			}
			Eigen::VectorXd r=_sparse_solver::findSearchDirection(mat,_vecs,_rhs,salt);
			array<double>^ ret=gcnew array<double>(r.size());
			for (int i = 0; i < r.size(); i++)
			{
				ret[i] = r(i);
			}
			return ret;
		}
		static int solve_CHOLESKY(int n, array<int>^ columnptr, array<int>^ rowindices, array<double>^ values, array<double>^ rhs, array<double>^ ret,double salt)
		{
			

			Eigen::SparseMatrix<double> mat;
			Eigen::VectorXd b(n);
			Eigen::VectorXd _r(n);
			for (int i = 0; i < n; i++)
			{
				double v = rhs[i];
				b(i) = v;
			}
			std::vector<Eigen::Triplet<double, int64_t>> tripletList;
			tripletList.reserve(values->Length);
			for (int i = 0; i < n; i++)
			{
				int start = columnptr[i];
				int end = columnptr[i + 1];
				for (int j = start; j < end; j++)
				{
					int row = rowindices[j];
					double val = values[j];
					tripletList.push_back(Eigen::Triplet<double, int64_t>(row, i, val));
						//	if (row != i)
					//	tripletList.push_back(Eigen::Triplet<double, int64_t>(i, row, val));
				}
			}
			for (int i = 0; i < n; i++)
			{
				tripletList.push_back(Eigen::Triplet<double, int64_t>(i, i, salt));
				
			}
			mat.resize(n, n);
			mat.setFromTriplets(tripletList.begin(), tripletList.end());
			mat.makeCompressed();

			_r=_sparse_solver::solve_CHOLECKY(mat, b,salt);
			for (int i = 0; i < n; i++)
			{
				ret[i]=_r(i);
			}

			return 1;
		}
		static void genEigen2(double a, double b, double d, double A, double B, double D,
			[Runtime::InteropServices::Out]double% aa, [Runtime::InteropServices::Out]double% bb,
			[Runtime::InteropServices::Out]double% cc, [Runtime::InteropServices::Out]double% dd,
			[Runtime::InteropServices::Out]double% l1, [Runtime::InteropServices::Out]double% l1i,
			[Runtime::InteropServices::Out]double% l2, [Runtime::InteropServices::Out]double% l2i)
		{
			// 1. 行列の初期化（カンマ初期化子を使うと簡潔です）
			Eigen::MatrixXd ma(2, 2);
			ma << a, b,
				b, d;

			Eigen::MatrixXd mb(2, 2);
			mb << A, B,
				B, D;

			// 2. ネイティブの計算結果を受け取るための変数を用意
			// double% (マネージ参照) のアドレスは直接取れないため、ローカル変数を使います
			double nl1, nl2, nl1i, nl2i;

			// 3. ネイティブ関数の呼び出し
			// 返り値が複素数行列を想定している場合、Matrix2cd を使用します
			Eigen::MatrixXcd r = 
				_sparse_solver::genEigen2(ma, mb, &nl1, &nl2, &nl1i, &nl2i);

			// 4. マネージ参照（引数）へ結果を戻す
			aa = r.real()(0, 0);
			bb = r.real()(0, 1);
			cc = r.real()(1, 0);
			dd = r.real()(1, 1);

			l1 = nl1;
			l2 = nl2;
			l1i = nl1i;
			l2i = nl2i;
		}
	};
}
