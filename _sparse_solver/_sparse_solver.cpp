// _sparse_solver.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"

#include "_sparse_solver.h"
namespace _sparse_solver {
	Eigen::MatrixXcd genEigen2(Eigen::MatrixXd ma, Eigen::MatrixXd mb, double* l1, double* l2, double* l1i, double* l2i)
	{
		Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solve(ma, mb, true);
		*l1 = solve.eigenvalues()(0).real();
		*l2 = solve.eigenvalues()(1).real();
		*l1i = solve.eigenvalues()(0).imag();
		*l2i = solve.eigenvalues()(1).imag();

		return solve.eigenvectors().real();
	}

	Eigen::VectorXd solve_CHOLECKY(Eigen::SparseMatrix<double, 0, int64_t> mat, Eigen::VectorXd rhs)
	{
		/*auto _mt = omp_get_max_threads();
		int _mt2 = 0;
#pragma omp parallel
		{
#pragma omp single
			_mt2 = omp_get_num_threads();
		}
		if (_mt2 > _mt)_mt = _mt2;
		Eigen::setNbThreads(_mt - 1);
		omp_set_num_threads(_mt - 1);*/
		Eigen::PardisoLDLT< Eigen::SparseMatrix<double, 0, int64_t>> chol;
		
		int n = mat.rows();
		chol.compute(mat);
		int rv = 0;
		if (chol.info() != Eigen::Success) {
			rv = chol.info();
			return Eigen::VectorXd(n);
		}
		Eigen::VectorXd _r(n);

		_r = chol.solve(rhs);
		return _r;
	}
}
