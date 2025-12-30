#pragma once
#include "eigen-5.0.0/Eigen/PardisoSupport"
#include "eigen-5.0.0/Eigen/Sparse"
#include "eigen-5.0.0/Eigen/Dense"
#include "eigen-5.0.0/Eigen/SparseQR"
#include "eigen-5.0.0/Eigen/SparseLU"
#include "eigen-5.0.0/Eigen/SparseCholesky"	
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int64_t
#define EIGEN_DONT_ALIGN_STATICALLY
#define EIGEN_MAX_ALIGN_BYTES 0
#define EIGEN_DONT_VECTORIZE
#include <chrono>
#include <vector>
#include <map>
#include <string>
using namespace std::chrono;
using std::vector;
using std::string;

using namespace System;

namespace sparse_solver {
	public ref class sparse_solver abstract sealed
	{
	public:
		static int solve_CHOLESKY(int n, array<int>^ columnptr, array<int>^ rowindices, array<double>^ values, array<double>^ rhs, array<double>^ ret)
		{
			int rv = 0;
			Eigen::PardisoLDLT< Eigen::SparseMatrix<double,0, int64_t>> chol;
			Eigen::SparseMatrix<double> mat;
			Eigen::VectorXd b(n);
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
			
			mat.resize(n, n);
			mat.setFromTriplets(tripletList.begin(), tripletList.end());
			mat.makeCompressed();
			
			chol.compute(mat);
			if (chol.info() != Eigen::Success) {
				rv = chol.info();
				return rv;
			}
			Eigen::VectorXd _r(n);
			
			_r= chol.solve(b);
			for (int i = 0; i < n; i++)
			{
				ret[i] = _r(i);
			}
			if (chol.info() != Eigen::Success) {
				rv = chol.info();
				return rv+10;
			}

			return chol.info();
		}
	};
}
