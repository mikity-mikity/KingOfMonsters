#pragma once

#include <vector>
#include <cmath>
#include <cstdint>
#include <algorithm>
#include <limits>
#include "_sparse_solver.h"

using namespace System;

namespace sparse_solver {

    public ref class sparse_solver abstract sealed
    {
    public:

        static array<double>^ findSearchDirection(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,
            array<array<double>^>^ vecs,
            array<double>^ rhs,
            double epsilon,
            [System::Runtime::InteropServices::Out]
            array<double>^% parameters)
        {
            using SparseMatrix64 = Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;
            using Triplet64 = Eigen::Triplet<double, int64_t>;

            if (n <= 0) throw gcnew ArgumentException("n must be positive");
            if (columnptr == nullptr) throw gcnew ArgumentNullException("columnptr");
            if (rowindices == nullptr) throw gcnew ArgumentNullException("rowindices");
            if (values == nullptr) throw gcnew ArgumentNullException("values");
            if (vecs == nullptr) throw gcnew ArgumentNullException("vecs");
            if (rhs == nullptr) throw gcnew ArgumentNullException("rhs");
            if (!std::isfinite(epsilon) || epsilon < 0.0)
                throw gcnew ArgumentException("epsilon must be finite and non-negative");
            if (rhs->Length != n) throw gcnew ArgumentException("rhs has incompatible size");
            if (columnptr->Length != n + 1)
                throw gcnew ArgumentException("columnptr must have length n + 1");
            if (rowindices->Length != values->Length)
                throw gcnew ArgumentException("rowindices and values have incompatible sizes");
            if (columnptr[0] != 0 || columnptr[n] != values->Length)
                throw gcnew ArgumentException("columnptr has invalid start or end");

            const Eigen::Index reducedSize = static_cast<Eigen::Index>(vecs->Length);
            if (reducedSize == 0)
            {
                parameters = gcnew array<double>(0);
                return gcnew array<double>(n);
            }

            // Keep the original basis: no orthonormalization.
            Eigen::MatrixXd P(n, reducedSize);
            for (Eigen::Index j = 0; j < reducedSize; ++j)
            {
                array<double>^ src = vecs[static_cast<int>(j)];
                if (src == nullptr) throw gcnew ArgumentException("vecs contains a null vector");
                if (src->Length != n)
                    throw gcnew ArgumentException("vecs contains a vector with incompatible size");
                for (int i = 0; i < n; ++i)
                {
                    double value = src[i];
                    if (!std::isfinite(value))
                        throw gcnew ArgumentException("vecs contains NaN or infinity");
                    P(i, j) = value;
                }
            }
            if (!P.allFinite())
                throw gcnew InvalidOperationException("Basis matrix contains NaN or infinity");

            SparseMatrix64 matA(n, n);
            std::vector<Triplet64> tripletsA;
            tripletsA.reserve(static_cast<std::size_t>(values->Length));
            for (int column = 0; column < n; ++column)
            {
                const int start = columnptr[column];
                const int end = columnptr[column + 1];
                if (start < 0 || end < start || end > values->Length)
                    throw gcnew ArgumentException("columnptr contains an invalid range");
                for (int index = start; index < end; ++index)
                {
                    const int row = rowindices[index];
                    const double value = values[index];
                    if (row < 0 || row >= n)
                        throw gcnew ArgumentException("rowindices contains an invalid row index");
                    if (!std::isfinite(value))
                        throw gcnew ArgumentException("values contains NaN or infinity");
                    tripletsA.emplace_back(static_cast<int64_t>(row),
                        static_cast<int64_t>(column), value);
                }
            }
            matA.setFromTriplets(tripletsA.begin(), tripletsA.end());
            matA.makeCompressed();

            Eigen::VectorXd b(n);
            for (int i = 0; i < n; ++i)
            {
                const double value = rhs[i];
                if (!std::isfinite(value))
                    throw gcnew ArgumentException("rhs contains NaN or infinity");
                b(i) = value;
            }

            Eigen::MatrixXd AP(n, reducedSize);
            AP.noalias() = matA * P;
            if (!AP.allFinite())
                throw gcnew InvalidOperationException("A * P contains NaN or infinity");

            Eigen::MatrixXd H(reducedSize, reducedSize);
            Eigen::VectorXd g(reducedSize);
            H.noalias() = P.transpose() * AP;
            g.noalias() = P.transpose() * b;
            if (!H.allFinite())
                throw gcnew InvalidOperationException("Reduced matrix contains NaN or infinity");
            if (!g.allFinite())
                throw gcnew InvalidOperationException("Reduced rhs contains NaN or infinity");

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(
                H, Eigen::ComputeThinU | Eigen::ComputeThinV);
            const Eigen::VectorXd singularValues = svd.singularValues();
            if (!singularValues.allFinite())
                throw gcnew InvalidOperationException("Reduced SVD returned NaN or infinity");

            Eigen::VectorXd alpha = Eigen::VectorXd::Zero(reducedSize);
            if (singularValues.size() > 0)
            {
                const double sigmaMax = singularValues(0);
                if (!std::isfinite(sigmaMax) || sigmaMax < 0.0)
                    throw gcnew InvalidOperationException("Largest singular value is invalid");
                if (sigmaMax > 0.0)
                {
                    const double threshold = epsilon * sigmaMax;
                    Eigen::VectorXd transformedRhs = svd.matrixU().transpose() * g;
                    for (Eigen::Index i = 0; i < singularValues.size(); ++i)
                    {
                        if (singularValues(i) > threshold)
                            transformedRhs(i) /= singularValues(i);
                        else
                            transformedRhs(i) = 0.0;
                    }
                    alpha.noalias() = svd.matrixV() * transformedRhs;
                }
            }
            if (!alpha.allFinite())
                throw gcnew InvalidOperationException("Reduced solution contains NaN or infinity");

            Eigen::VectorXd direction(n);
            direction.noalias() = P * alpha;
            if (!direction.allFinite())
                throw gcnew InvalidOperationException("Full-space direction contains NaN or infinity");

            array<double>^ ret = gcnew array<double>(n);
            parameters = gcnew array<double>(static_cast<int>(reducedSize));
            for (int i = 0; i < n; ++i) ret[i] = direction(i);
            for (Eigen::Index j = 0; j < reducedSize; ++j)
                parameters[static_cast<int>(j)] = alpha(j);
            return ret;
        }

        static int countEigenvaluesBelow(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,
            double threshold)
        {
            using SparseMatrix64 = Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;
            using Triplet64 = Eigen::Triplet<double, int64_t>;

            std::vector<Triplet64> tripletList;
            tripletList.reserve(values->Length + n);
            for (int64_t i = 0; i < n; i++)
            {
                int start = columnptr[i];
                int end = columnptr[i + 1];
                for (int j = start; j < end; j++)
                {
                    int64_t row = rowindices[j];
                    double val = values[j];
                    tripletList.push_back(Triplet64(row, i, val));
                }
            }
            // B = A - threshold * I
            for (int i = 0; i < n; i++)
                tripletList.push_back(Triplet64(i, i, -threshold));

            SparseMatrix64 mat(n, n);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
            mat.makeCompressed();

            Eigen::PardisoLDLT<SparseMatrix64> solver;
            solver.compute(mat);
            if (solver.info() != Eigen::Success) return -1;
            auto& iparm = solver.pardisoParameterArray();
            return iparm[22];
        }

        static int solve_CHOLESKY(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,
            array<double>^ rhs,
            array<double>^ ret,
            double epsilon,
            int numthreads)
        {
            using SparseMatrix64 = Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;
            using Triplet64 = Eigen::Triplet<double, int64_t>;
            SparseMatrix64 mat;
            Eigen::VectorXd b(n);
            Eigen::VectorXd _r(n);

            for (int i = 0; i < n; i++) b(i) = rhs[i];

            std::vector<Triplet64> tripletList;
            tripletList.reserve(values->Length + n);
            double matrixScale = 0.0;
            for (int i = 0; i < n; i++)
            {
                int start = columnptr[i];
                int end = columnptr[i + 1];
                for (int j = start; j < end; j++)
                {
                    int row = rowindices[j];
                    double val = values[j];
                    tripletList.push_back(Triplet64(row, i, val));
                    if (row == i)
                    {
                        double a = std::abs(val);
                        if (a > matrixScale) matrixScale = a;
                    }
                }
            }
            if (matrixScale == 0.0) matrixScale = 1.0;
            double salt = epsilon * matrixScale;
            for (int i = 0; i < n; i++)
                tripletList.push_back(Triplet64(i, i, salt));

            mat.resize(n, n);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
            mat.makeCompressed();

            Eigen::PardisoLU<SparseMatrix64> solver;
            solver.compute(mat);
            for (int i = 0; i < 50; i++)
            {
                if (solver.info() != Eigen::Success)
                {
                    double _salt = salt * std::pow(10.0, i + 1);
                    for (int j = 0; j < n; j++) mat.coeffRef(j, j) += _salt;
                    solver.compute(mat);
                }
                else break;
            }
            _r = solver.solve(b);
            for (int i = 0; i < n; i++) ret[i] = _r(i);
            return 1;
        }

        static int solve_CHOLESKY_MULTI(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,
            array<array<double>^>^ rhs,
            array<array<double>^>^ ret,
            double epsilon,
            int numthreads)
        {
            using SparseMatrix64 = Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;
            using Triplet64 = Eigen::Triplet<double, int64_t>;

            int nrhs = rhs->Length;
            if (ret->Length != nrhs)
                throw gcnew ArgumentException("rhs and ret must have same number of vectors.");
            for (int k = 0; k < nrhs; k++)
            {
                if (rhs[k] == nullptr || rhs[k]->Length != n)
                    throw gcnew ArgumentException("Each rhs vector must have length n.");
                if (ret[k] == nullptr || ret[k]->Length != n)
                    throw gcnew ArgumentException("Each ret vector must have length n.");
            }

            SparseMatrix64 mat;
            std::vector<Triplet64> tripletList;
            tripletList.reserve(values->Length + n);
            double matrixScale = 0.0;
            for (int i = 0; i < n; i++)
            {
                int start = columnptr[i];
                int end = columnptr[i + 1];
                for (int j = start; j < end; j++)
                {
                    int row = rowindices[j];
                    double val = values[j];
                    tripletList.push_back(Triplet64(static_cast<int64_t>(row),
                        static_cast<int64_t>(i), val));
                    if (row == i)
                    {
                        double a = std::abs(val);
                        if (a > matrixScale) matrixScale = a;
                    }
                }
            }
            if (matrixScale == 0.0) matrixScale = 1.0;
            double salt = epsilon * matrixScale;
            for (int i = 0; i < n; i++)
                tripletList.push_back(Triplet64(static_cast<int64_t>(i),
                    static_cast<int64_t>(i), salt));

            mat.resize(n, n);
            mat.setFromTriplets(tripletList.begin(), tripletList.end());
            mat.makeCompressed();

            Eigen::PardisoLU<SparseMatrix64> solver;
            solver.compute(mat);
            for (int i = 0; i < 50; i++)
            {
                if (solver.info() == Eigen::Success) break;
                double _salt = salt * std::pow(10.0, i + 1);
                for (int j = 0; j < n; j++) mat.coeffRef(j, j) += _salt;
                mat.makeCompressed();
                solver.compute(mat);
            }
            if (solver.info() != Eigen::Success) return 0;

            Eigen::VectorXd b(n);
            for (int k = 0; k < nrhs; k++)
            {
                for (int i = 0; i < n; i++) b(i) = rhs[k][i];
                Eigen::VectorXd x = solver.solve(b);
                if (solver.info() != Eigen::Success) return 0;
                for (int i = 0; i < n; i++) ret[k][i] = x(i);
            }
            return 1;
        }

        static void genEigen2(
            double a, double b, double d, double A, double B, double D,
            [Runtime::InteropServices::Out] double% aa,
            [Runtime::InteropServices::Out] double% bb,
            [Runtime::InteropServices::Out] double% cc,
            [Runtime::InteropServices::Out] double% dd,
            [Runtime::InteropServices::Out] double% l1,
            [Runtime::InteropServices::Out] double% l1i,
            [Runtime::InteropServices::Out] double% l2,
            [Runtime::InteropServices::Out] double% l2i)
        {
            Eigen::MatrixXd ma(2, 2);
            ma << a, b,
                b, d;
            Eigen::MatrixXd mb(2, 2);
            mb << A, B,
                B, D;

            double nl1, nl2, nl1i, nl2i;
            Eigen::MatrixXcd r =
                _sparse_solver::genEigen2(ma, mb, &nl1, &nl2, &nl1i, &nl2i);

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

} // namespace sparse_solver
