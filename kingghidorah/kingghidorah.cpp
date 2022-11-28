//#include "pch.h"
#include"mySparseLibrary.h"
#include "kingghidorah_S.h"
#include "kingghidorah_C.h"

void main()
{
	Eigen::SparseMatrix<double, 1, int> mat(5, 5);
	mat.coeffRef(0, 0) = 1;
	mat.coeffRef(1, 1) = 1;
	mat.coeffRef(2, 2) = 1;
	mat.coeffRef(3, 3) = 1;
	mat.coeffRef(4, 4) = 1;
	Eigen::VectorXd rhs(5);
	rhs(0) = 0;
	rhs(1) = 1;
	rhs(2) = 0;
	rhs(3) = 0;
	rhs(4) = 0;
	Eigen::PardisoLU<Eigen::SparseMatrix<double, 0, int>> lu;
	lu.compute(mat);
	Eigen::VectorXd ret=lu.solve(rhs);
	std::cout << ret;
}