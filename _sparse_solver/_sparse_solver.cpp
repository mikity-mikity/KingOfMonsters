// _sparse_solver.cpp : Defines the functions for the static library.
//

#include "pch.h"
#include "framework.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstdlib>
#include <string>
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
      
  

}
