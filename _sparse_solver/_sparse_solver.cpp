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
    Eigen::MatrixXcd genEigen2(
        Eigen::MatrixXd ma,
        Eigen::MatrixXd mb,
        double* l1, double* l2,
        double* l1i, double* l2i)
    {
        const int n = static_cast<int>(ma.rows());

        auto fallback = [&]() -> Eigen::MatrixXcd
            {
                *l1 = 1.0;
                *l2 = 1.0;
                *l1i = 0.0;
                *l2i = 0.0;

                return Eigen::MatrixXcd::Identity(n, n);
            };

        if (n != 2 || ma.cols() != 2 ||
            mb.rows() != 2 || mb.cols() != 2 ||
            !ma.allFinite() || !mb.allFinite())
        {
            return fallback();
        }

        try
        {
            Eigen::GeneralizedEigenSolver<Eigen::MatrixXd> solve;
            solve.compute(ma, mb, true);

            if (solve.info() != Eigen::Success)
                return fallback();

            auto eigenvalues = solve.eigenvalues();
            auto eigenvectors = solve.eigenvectors();

            if (!eigenvalues.allFinite() ||
                !eigenvectors.allFinite())
                return fallback();

            const double eps = 1e-10;

            if (std::abs(eigenvalues(0).imag()) > eps ||
                std::abs(eigenvalues(1).imag()) > eps)
                return fallback();

            *l1 = eigenvalues(0).real();
            *l2 = eigenvalues(1).real();
            *l1i = eigenvalues(0).imag();
            *l2i = eigenvalues(1).imag();

            return eigenvectors;
        }
        catch (...)
        {
            return fallback();
        }
    }
      
  

}
