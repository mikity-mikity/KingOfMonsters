#pragma once

#include <chrono>
#include <vector>
#include <map>
#include <random>
#include <string>
using namespace std::chrono;
using std::vector;
using std::string;
#include "_sparse_solver.h"
using namespace System;
using namespace System;
using namespace System::Collections::Generic;

namespace sparse_solver {
   
	public ref class sparse_solver abstract sealed
	{
	public:
        static array<double>^ removePAndAPComponents(
            int n,

            // Primary matrix A (CSC)
            array<int>^ columnptrA,
            array<int>^ rowindicesA,
            array<double>^ valuesA,

            // P = [vecs[0] vecs[1] ...]
            array<array<double>^>^ vecs,

            // Vector from which the P and AP components are removed
            array<double>^ rhsb,

            // Relative threshold for determining rank([P, AP]).
            // For example: 1.0e-12.
            double tsvdThreshold)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;

            using Triplet64 =
                Eigen::Triplet<double, int64_t>;


            // =====================================================
            // Validation
            // =====================================================

            if (n <= 0)
                throw gcnew ArgumentException("n must be positive");

            if (columnptrA == nullptr)
                throw gcnew ArgumentNullException("columnptrA");

            if (rowindicesA == nullptr)
                throw gcnew ArgumentNullException("rowindicesA");

            if (valuesA == nullptr)
                throw gcnew ArgumentNullException("valuesA");

            if (vecs == nullptr)
                throw gcnew ArgumentNullException("vecs");

            if (rhsb == nullptr)
                throw gcnew ArgumentNullException("rhsb");

            if (rhsb->Length != n)
                throw gcnew ArgumentException("rhsb has incompatible size");

            if (columnptrA->Length != n + 1 ||
                columnptrA[0] != 0 ||
                columnptrA[n] != valuesA->Length ||
                rowindicesA->Length != valuesA->Length)
            {
                throw gcnew ArgumentException("Invalid CSC arrays for A");
            }

            if (!std::isfinite(tsvdThreshold) ||
                tsvdThreshold < 0.0)
            {
                throw gcnew ArgumentException(
                    "tsvdThreshold must be finite and nonnegative");
            }


            // =====================================================
            // Build sparse primary matrix A.
            // =====================================================

            SparseMatrix64 matA(n, n);

            std::vector<Triplet64> tripletsA;
            tripletsA.reserve(static_cast<std::size_t>(valuesA->Length));

            for (int column = 0; column < n; ++column)
            {
                const int start = columnptrA[column];
                const int end = columnptrA[column + 1];

                if (start < 0 || end < start || end > valuesA->Length)
                    throw gcnew ArgumentException(
                        "columnptrA contains an invalid range");

                for (int index = start; index < end; ++index)
                {
                    const int row = rowindicesA[index];
                    const double value = valuesA[index];

                    if (row < 0 || row >= n)
                        throw gcnew ArgumentException(
                            "rowindicesA contains an invalid row index");

                    if (!std::isfinite(value))
                        throw gcnew ArgumentException(
                            "valuesA contains NaN or infinity");

                    tripletsA.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }

            matA.setFromTriplets(tripletsA.begin(), tripletsA.end());
            matA.makeCompressed();


            // =====================================================
            // Build P directly from vecs and compute AP = A P.
            // =====================================================

            const Eigen::Index pCount =
                static_cast<Eigen::Index>(vecs->Length);

            Eigen::MatrixXd P(n, pCount);

            for (Eigen::Index j = 0; j < pCount; ++j)
            {
                array<double>^ source = vecs[static_cast<int>(j)];

                if (source == nullptr)
                    throw gcnew ArgumentException("vecs contains a null vector");

                if (source->Length != n)
                    throw gcnew ArgumentException(
                        "vecs contains a vector with incompatible size");

                for (int i = 0; i < n; ++i)
                {
                    const double value = source[i];

                    if (!std::isfinite(value))
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");

                    P(i, j) = value;
                }
            }

            Eigen::MatrixXd AP(n, pCount);
            AP.noalias() = matA * P;

            if (!AP.allFinite())
                throw gcnew InvalidOperationException(
                    "A * P contains NaN or infinity");


            // =====================================================
            // Build the forbidden subspace
            //
            //     span(C) = span(P, AP).
            //
            // Each nonzero column is normalized first.  Therefore the
            // rank decision does not depend on the arbitrary magnitudes
            // of individual vectors in P or AP.
            // =====================================================

            Eigen::MatrixXd C(n, 2 * pCount);

            for (Eigen::Index j = 0; j < pCount; ++j)
            {
                C.col(j) = P.col(j);
                C.col(pCount + j) = AP.col(j);
            }

            for (Eigen::Index j = 0; j < C.cols(); ++j)
            {
                const double columnNorm = C.col(j).norm();

                if (!std::isfinite(columnNorm))
                    throw gcnew InvalidOperationException(
                        "P or A * P has an invalid column norm");

                if (columnNorm > 0.0)
                    C.col(j) /= columnNorm;
            }


            // =====================================================
            // Construct an orthonormal basis Q for span(P, AP).
            //
            // Only this temporary constraint basis is orthonormalized.
            // The original vecs are not modified.
            // =====================================================

            Eigen::MatrixXd Q(n, 0);

            if (C.cols() > 0)
            {
                Eigen::JacobiSVD<Eigen::MatrixXd> svd(
                    C,
                    Eigen::ComputeThinU);

                const Eigen::VectorXd singularValues =
                    svd.singularValues();

                Eigen::Index rank = 0;

                if (singularValues.size() > 0 && singularValues(0) > 0.0)
                {
                    const double cutoff =
                        tsvdThreshold * singularValues(0);

                    for (Eigen::Index i = 0;
                        i < singularValues.size();
                        ++i)
                    {
                        if (singularValues(i) > cutoff)
                            ++rank;
                    }
                }

                Q = svd.matrixU().leftCols(rank);
            }


            // =====================================================
            // Copy rhsb and remove its components in span(P, AP):
            //
            //     direction = (I - Q Q^T) rhsb.
            //
            // Thus, up to floating-point error,
            //
            //     P^T direction    = 0,
            //     (AP)^T direction = 0.
            // =====================================================

            Eigen::VectorXd direction(n);

            for (int i = 0; i < n; ++i)
            {
                const double value = rhsb[i];

                if (!std::isfinite(value))
                    throw gcnew ArgumentException(
                        "rhsb contains NaN or infinity");

                direction(i) = value;
            }

            if (Q.cols() > 0)
            {
                Eigen::VectorXd coefficients = Q.transpose() * direction;
                direction -= Q * coefficients;

                // Reproject once to remove accumulated roundoff error.
                coefficients.noalias() = Q.transpose() * direction;
                direction -= Q * coefficients;
            }

            if (!direction.allFinite())
                throw gcnew InvalidOperationException(
                    "Projected direction contains NaN or infinity");


            // =====================================================
            // Return the projected vector.
            // =====================================================

            array<double>^ ret = gcnew array<double>(n);

            for (int i = 0; i < n; ++i)
                ret[i] = direction(i);

            return ret;
        }


        static array<double>^ findSearchDirection(
            int n,

            // =====================================================
            // Matrix A
            // CSC format
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Original basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // NO orthonormalization is performed.
            // =====================================================
            array<array<double>^>^ vecs,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,

            // =====================================================
            // Relative TSVD threshold for
            //
            //     P^T A P
            // =====================================================
            double epsilon,

            // =====================================================
            // Output:
            //
            // coefficients alpha with respect to ORIGINAL basis P
            //
            //     direction = P alpha
            // =====================================================
            array<double>^% parameters)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (columnptr == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "columnptr");
            }

            if (rowindices == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rowindices");
            }

            if (values == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "values");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (rhs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs");
            }

            if (!std::isfinite(epsilon) ||
                epsilon < 0.0)
            {
                throw gcnew ArgumentException(
                    "epsilon must be finite and non-negative");
            }

            if (rhs->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs has incompatible size");
            }

            if (columnptr->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");
            }

            if (rowindices->Length != values->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");
            }

            if (columnptr[0] != 0 ||
                columnptr[n] != values->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr has invalid start or end");
            }


            // =====================================================
            // Number of basis vectors
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    vecs->Length);


            // =====================================================
            // Empty basis
            // =====================================================

            if (reducedSize == 0)
            {
                parameters =
                    gcnew array<double>(0);

                return gcnew array<double>(n);
            }


            // =====================================================
            // Build original dense basis
            //
            //     P = [p0 p1 ...]
            //
            // NO QR / orthonormalization.
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                array<double>^ src =
                    vecs[
                        static_cast<int>(j)];

                if (src == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a null vector");
                }

                if (src->Length != n)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a vector with incompatible size");
                }

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double value =
                        src[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");
                    }

                    P(i, j) =
                        value;
                }
            }


            if (!P.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Basis matrix contains NaN or infinity");
            }


            // =====================================================
            // Build sparse matrix A
            // =====================================================

            SparseMatrix64 matA(
                n,
                n);

            std::vector<Triplet64> tripletsA;

            tripletsA.reserve(
                static_cast<std::size_t>(
                    values->Length));


            for (int column = 0;
                column < n;
                ++column)
            {
                const int start =
                    columnptr[column];

                const int end =
                    columnptr[column + 1];

                if (start < 0 ||
                    end < start ||
                    end > values->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr contains an invalid range");
                }


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices[index];

                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "rowindices contains an invalid row index");
                    }


                    const double value =
                        values[index];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "values contains NaN or infinity");
                    }


                    tripletsA.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }


            matA.setFromTriplets(
                tripletsA.begin(),
                tripletsA.end());

            matA.makeCompressed();


            // =====================================================
            // RHS
            // =====================================================

            Eigen::VectorXd b(
                n);

            for (int i = 0;
                i < n;
                ++i)
            {
                const double value =
                    rhs[i];

                if (!std::isfinite(value))
                {
                    throw gcnew ArgumentException(
                        "rhs contains NaN or infinity");
                }

                b(i) =
                    value;
            }


            // =====================================================
            // Compute
            //
            //     AP = A P
            // =====================================================

            Eigen::MatrixXd AP(
                n,
                reducedSize);

            AP.noalias() =
                matA * P;


            if (!AP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "A * P contains NaN or infinity");
            }


            // =====================================================
            // Reduced primary system
            //
            //     H = P^T A P
            //
            //     g = P^T rhs
            //
            // Solve
            //
            //     H alpha = g
            //
            // =====================================================

            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);


            H.noalias() =
                P.transpose() * AP;

            g.noalias() =
                P.transpose() * b;


            if (!H.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix contains NaN or infinity");
            }

            if (!g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced rhs contains NaN or infinity");
            }


            // =====================================================
            // TSVD solve
            //
            //     (P^T A P) alpha = P^T rhs
            //
            // No basis orthonormalization is performed.
            //
            // If P contains linearly dependent or nearly dependent
            // vectors, TSVD handles the rank deficiency here.
            // =====================================================

            Eigen::JacobiSVD<
                Eigen::MatrixXd>
                svd(
                    H,
                    Eigen::ComputeThinU |
                    Eigen::ComputeThinV);


            const Eigen::VectorXd singularValues =
                svd.singularValues();


            if (!singularValues.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced SVD returned NaN or infinity");
            }


            Eigen::VectorXd alpha(
                reducedSize);

            alpha.setZero();


            if (singularValues.size() > 0)
            {
                const double sigmaMax =
                    singularValues(0);


                if (!std::isfinite(sigmaMax) ||
                    sigmaMax < 0.0)
                {
                    throw gcnew InvalidOperationException(
                        "Largest singular value is invalid");
                }


                if (sigmaMax > 0.0)
                {
                    const double threshold =
                        epsilon* sigmaMax;


                    Eigen::VectorXd transformedRhs =
                        svd.matrixU().transpose() *
                        g;


                    for (Eigen::Index i = 0;
                        i < singularValues.size();
                        ++i)
                    {
                        const double sigma =
                            singularValues(i);


                        if (sigma > threshold)
                        {
                            transformedRhs(i) /=
                                sigma;
                        }
                        else
                        {
                            transformedRhs(i) =
                                0.0;
                        }
                    }


                    alpha.noalias() =
                        svd.matrixV() *
                        transformedRhs;
                }
            }


            if (!alpha.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Full-space direction
            //
            //     direction = P alpha
            //
            // Since P is the ORIGINAL basis,
            // alpha is directly the coefficient vector that we want.
            // =====================================================

            Eigen::VectorXd direction(
                n);

            direction.noalias() =
                P * alpha;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(n);

            parameters =
                gcnew array<double>(
                    static_cast<int>(
                        reducedSize));


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    direction(i);
            }


            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                parameters[
                    static_cast<int>(j)] =
                    alpha(j);
            }


            return ret;
        }
        static array<double>^ findSearchDirectionFixed(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,
            array<array<double>^>^ vecs,
            array<double>^ rhs,
            double epsilon,
            int fixedIndex,
            double fixedValue,
            array<double>^% parameters)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<double, Eigen::ColMajor, int64_t>;
            using Triplet64 =
                Eigen::Triplet<double, int64_t>;

            if (n <= 0)
                throw gcnew ArgumentException("n must be positive");
            if (columnptr == nullptr)
                throw gcnew ArgumentNullException("columnptr");
            if (rowindices == nullptr)
                throw gcnew ArgumentNullException("rowindices");
            if (values == nullptr)
                throw gcnew ArgumentNullException("values");
            if (vecs == nullptr)
                throw gcnew ArgumentNullException("vecs");
            if (rhs == nullptr)
                throw gcnew ArgumentNullException("rhs");
            if (!std::isfinite(epsilon) || epsilon < 0.0)
                throw gcnew ArgumentException("epsilon must be finite and non-negative");
            if (!std::isfinite(fixedValue))
                throw gcnew ArgumentException("fixedValue must be finite");
            if (fixedIndex < 0 || fixedIndex >= vecs->Length)
                throw gcnew ArgumentOutOfRangeException("fixedIndex");
            if (rhs->Length != n)
                throw gcnew ArgumentException("rhs has incompatible size");
            if (columnptr->Length != n + 1)
                throw gcnew ArgumentException("columnptr must have length n + 1");
            if (rowindices->Length != values->Length)
                throw gcnew ArgumentException("rowindices and values have incompatible sizes");
            if (columnptr[0] != 0 || columnptr[n] != values->Length)
                throw gcnew ArgumentException("columnptr has invalid start or end");

            const Eigen::Index basisCount =
                static_cast<Eigen::Index>(vecs->Length);
            const Eigen::Index freeCount = basisCount - 1;

            Eigen::MatrixXd P(n, basisCount);
            for (Eigen::Index j = 0; j < basisCount; ++j)
            {
                array<double>^ src = vecs[static_cast<int>(j)];
                if (src == nullptr || src->Length != n)
                    throw gcnew ArgumentException("vecs contains an invalid vector");

                for (int i = 0; i < n; ++i)
                {
                    if (!std::isfinite(src[i]))
                        throw gcnew ArgumentException("vecs contains NaN or infinity");
                    P(i, j) = src[i];
                }
            }

            Eigen::VectorXd b(n);
            for (int i = 0; i < n; ++i)
            {
                if (!std::isfinite(rhs[i]))
                    throw gcnew ArgumentException("rhs contains NaN or infinity");
                b(i) = rhs[i];
            }

            SparseMatrix64 matA(n, n);
            std::vector<Triplet64> triplets;
            triplets.reserve(static_cast<std::size_t>(values->Length));

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
                        throw gcnew ArgumentException("rowindices contains an invalid index");
                    if (!std::isfinite(value))
                        throw gcnew ArgumentException("values contains NaN or infinity");

                    triplets.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }

            matA.setFromTriplets(triplets.begin(), triplets.end());
            matA.makeCompressed();

            Eigen::MatrixXd AP = matA * P;
            if (!AP.allFinite())
                throw gcnew InvalidOperationException("A * P contains NaN or infinity");

            // 固定係数の寄与を右辺から差し引く。
            Eigen::VectorXd adjustedRhs =
                b - fixedValue * AP.col(fixedIndex);

            Eigen::MatrixXd Pfree(n, freeCount);
            Eigen::MatrixXd APfree(n, freeCount);

            for (Eigen::Index j = 0, k = 0; j < basisCount; ++j)
            {
                if (j == fixedIndex)
                    continue;

                Pfree.col(k) = P.col(j);
                APfree.col(k) = AP.col(j);
                ++k;
            }

            Eigen::VectorXd alpha =
                Eigen::VectorXd::Zero(basisCount);
            alpha(fixedIndex) = fixedValue;

            if (freeCount > 0)
            {
                Eigen::MatrixXd H = Pfree.transpose() * APfree;
                Eigen::VectorXd g = Pfree.transpose() * adjustedRhs;

                if (!H.allFinite() || !g.allFinite())
                    throw gcnew InvalidOperationException(
                        "Reduced system contains NaN or infinity");

                Eigen::JacobiSVD<Eigen::MatrixXd> svd(
                    H, Eigen::ComputeThinU | Eigen::ComputeThinV);

                const Eigen::VectorXd singularValues = svd.singularValues();
                if (!singularValues.allFinite())
                    throw gcnew InvalidOperationException(
                        "Reduced SVD returned NaN or infinity");

                Eigen::VectorXd alphaFree =
                    Eigen::VectorXd::Zero(freeCount);

                if (singularValues.size() > 0 && singularValues(0) > 0.0)
                {
                    const double threshold =
                        epsilon * singularValues(0);

                    Eigen::VectorXd transformedRhs =
                        svd.matrixU().transpose() * g;

                    for (Eigen::Index i = 0; i < singularValues.size(); ++i)
                    {
                        if (singularValues(i) > threshold)
                            transformedRhs(i) /= singularValues(i);
                        else
                            transformedRhs(i) = 0.0;
                    }

                    alphaFree = svd.matrixV() * transformedRhs;
                }

                for (Eigen::Index j = 0, k = 0; j < basisCount; ++j)
                {
                    if (j != fixedIndex)
                        alpha(j) = alphaFree(k++);
                }
            }

            Eigen::VectorXd direction = P * alpha;
            if (!alpha.allFinite() || !direction.allFinite())
                throw gcnew InvalidOperationException(
                    "Solution contains NaN or infinity");

            parameters = gcnew array<double>(static_cast<int>(basisCount));
            array<double>^ ret = gcnew array<double>(n);

            for (Eigen::Index j = 0; j < basisCount; ++j)
                parameters[static_cast<int>(j)] = alpha(j);

            for (int i = 0; i < n; ++i)
                ret[i] = direction(i);

            return ret;
        }
        static array<double>^ findSearchDirection2(
            int n,

            // =====================================================
            // Matrix A
            // CSC format
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Original basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // NO orthonormalization is performed.
            // =====================================================
            array<array<double>^>^ vecs,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,

            // =====================================================
            // Reduced-space LM damping for
            //
            //     P^T A P + lambda P^T P
            // =====================================================
            double lambda,

            // =====================================================
            // Output:
            //
            // coefficients alpha with respect to ORIGINAL basis P
            //
            //     direction = P alpha
            // =====================================================
            array<double>^% parameters)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (columnptr == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "columnptr");
            }

            if (rowindices == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rowindices");
            }

            if (values == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "values");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (rhs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs");
            }

            if (!std::isfinite(lambda) ||
                lambda <= 0.0)
            {
                throw gcnew ArgumentException(
                    "lambda must be finite and positive");
            }

            if (rhs->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs has incompatible size");
            }

            if (columnptr->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");
            }

            if (rowindices->Length != values->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");
            }

            if (columnptr[0] != 0 ||
                columnptr[n] != values->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr has invalid start or end");
            }


            // =====================================================
            // Number of basis vectors
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    vecs->Length);


            // =====================================================
            // Empty basis
            // =====================================================

            if (reducedSize == 0)
            {
                parameters =
                    gcnew array<double>(0);

                return gcnew array<double>(n);
            }


            // =====================================================
            // Build original dense basis
            //
            //     P = [p0 p1 ...]
            //
            // NO QR / orthonormalization.
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                array<double>^ src =
                    vecs[
                        static_cast<int>(j)];

                if (src == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a null vector");
                }

                if (src->Length != n)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a vector with incompatible size");
                }

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double value =
                        src[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");
                    }

                    P(i, j) =
                        value;
                }
            }


            if (!P.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Basis matrix contains NaN or infinity");
            }


            // =====================================================
            // Build sparse matrix A
            // =====================================================

            SparseMatrix64 matA(
                n,
                n);

            std::vector<Triplet64> tripletsA;

            tripletsA.reserve(
                static_cast<std::size_t>(
                    values->Length));


            for (int column = 0;
                column < n;
                ++column)
            {
                const int start =
                    columnptr[column];

                const int end =
                    columnptr[column + 1];

                if (start < 0 ||
                    end < start ||
                    end > values->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr contains an invalid range");
                }


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices[index];

                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "rowindices contains an invalid row index");
                    }


                    const double value =
                        values[index];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "values contains NaN or infinity");
                    }


                    tripletsA.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }


            matA.setFromTriplets(
                tripletsA.begin(),
                tripletsA.end());

            matA.makeCompressed();


            // =====================================================
            // RHS
            // =====================================================

            Eigen::VectorXd b(
                n);

            for (int i = 0;
                i < n;
                ++i)
            {
                const double value =
                    rhs[i];

                if (!std::isfinite(value))
                {
                    throw gcnew ArgumentException(
                        "rhs contains NaN or infinity");
                }

                b(i) =
                    value;
            }


            // =====================================================
            // Compute
            //
            //     AP = A P
            // =====================================================

            Eigen::MatrixXd AP(
                n,
                reducedSize);

            AP.noalias() =
                matA * P;


            if (!AP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "A * P contains NaN or infinity");
            }


            // =====================================================
            // Reduced primary system
            //
            //     H = P^T A P + lambda P^T P
            //
            //     g = P^T rhs
            //
            // Solve
            //
            //     H alpha = g
            //
            // =====================================================

            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);


            H.noalias() =
                P.transpose() * AP;

            H.noalias() +=
                lambda * (P.transpose() * P);

            // Remove very small asymmetry caused by floating-point
            // accumulation. Eigen::LLT expects a symmetric matrix.
            H =
                0.5 * (H + H.transpose());

            g.noalias() =
                P.transpose() * b;


            if (!H.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix contains NaN or infinity");
            }

            if (!g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced rhs contains NaN or infinity");
            }


            // =====================================================
            // Cholesky solve
            //
            //     (P^T A P + lambda P^T P) alpha
            //         = P^T rhs
            //
            // No basis orthonormalization is performed.
            //
            // If P has full column rank and lambda > 0,
            // lambda P^T P makes the reduced matrix positive
            // definite even when A is rank deficient.
            // =====================================================

            Eigen::LLT<Eigen::MatrixXd> cholesky;

            cholesky.compute(H);

            if (cholesky.info() != Eigen::Success)
            {
                throw gcnew InvalidOperationException(
                    "Cholesky factorization failed. "
                    "The basis P may be linearly dependent, "
                    "or lambda may be too small.");
            }


            Eigen::VectorXd alpha =
                cholesky.solve(g);


            if (cholesky.info() != Eigen::Success)
            {
                throw gcnew InvalidOperationException(
                    "Cholesky solve failed");
            }


            if (!alpha.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Full-space direction
            //
            //     direction = P alpha
            //
            // Since P is the ORIGINAL basis,
            // alpha is directly the coefficient vector that we want.
            // =====================================================

            Eigen::VectorXd direction(
                n);

            direction.noalias() =
                P * alpha;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(n);

            parameters =
                gcnew array<double>(
                    static_cast<int>(
                        reducedSize));


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    direction(i);
            }


            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                parameters[
                    static_cast<int>(j)] =
                    alpha(j);
            }


            return ret;
        }
        static array<double>^ findSearchDirection3(
            int n,

            // =====================================================
            // Matrix A
            // CSC format
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Original basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // NO orthonormalization is performed.
            // =====================================================
            array<array<double>^>^ vecs,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,

            // =====================================================
            // Radius of the coefficient-space ball constraint
            //
            //     ||alpha|| <= trustRadius
            // =====================================================
            double trustRadius,

            // =====================================================
            // Output:
            //
            // coefficients alpha with respect to ORIGINAL basis P
            //
            //     direction = P alpha
            // =====================================================
            array<double>^% parameters,

            // LM multiplier selected by the trust-region solve.
            // Zero means that the unconstrained step was accepted.
            double% lambdaUsed)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (columnptr == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "columnptr");
            }

            if (rowindices == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rowindices");
            }

            if (values == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "values");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (rhs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs");
            }

            if (!std::isfinite(trustRadius) ||
                trustRadius <= 0.0)
            {
                throw gcnew ArgumentException(
                    "trustRadius must be finite and positive");
            }

            if (rhs->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs has incompatible size");
            }

            if (columnptr->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");
            }

            if (rowindices->Length != values->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");
            }

            if (columnptr[0] != 0 ||
                columnptr[n] != values->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr has invalid start or end");
            }


            // =====================================================
            // Number of basis vectors
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    vecs->Length);


            // =====================================================
            // Empty basis
            // =====================================================

            if (reducedSize == 0)
            {
                parameters =
                    gcnew array<double>(0);

                lambdaUsed =
                    0.0;

                return gcnew array<double>(n);
            }


            // =====================================================
            // Build original dense basis
            //
            //     P = [p0 p1 ...]
            //
            // NO QR / orthonormalization.
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                array<double>^ src =
                    vecs[
                        static_cast<int>(j)];

                if (src == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a null vector");
                }

                if (src->Length != n)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a vector with incompatible size");
                }

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double value =
                        src[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");
                    }

                    P(i, j) =
                        value;
                }
            }


            if (!P.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Basis matrix contains NaN or infinity");
            }


            // =====================================================
            // Build sparse matrix A
            // =====================================================

            SparseMatrix64 matA(
                n,
                n);

            std::vector<Triplet64> tripletsA;

            tripletsA.reserve(
                static_cast<std::size_t>(
                    values->Length));


            for (int column = 0;
                column < n;
                ++column)
            {
                const int start =
                    columnptr[column];

                const int end =
                    columnptr[column + 1];

                if (start < 0 ||
                    end < start ||
                    end > values->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr contains an invalid range");
                }


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices[index];

                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "rowindices contains an invalid row index");
                    }


                    const double value =
                        values[index];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "values contains NaN or infinity");
                    }


                    tripletsA.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }


            matA.setFromTriplets(
                tripletsA.begin(),
                tripletsA.end());

            matA.makeCompressed();


            // =====================================================
            // RHS
            // =====================================================

            Eigen::VectorXd b(
                n);

            for (int i = 0;
                i < n;
                ++i)
            {
                const double value =
                    rhs[i];

                if (!std::isfinite(value))
                {
                    throw gcnew ArgumentException(
                        "rhs contains NaN or infinity");
                }

                b(i) =
                    value;
            }


            // =====================================================
            // Compute
            //
            //     AP = A P
            // =====================================================

            Eigen::MatrixXd AP(
                n,
                reducedSize);

            AP.noalias() =
                matA * P;


            if (!AP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "A * P contains NaN or infinity");
            }


            // =====================================================
            // Reduced primary system
            //
            //     H0 = P^T A P
            //
            //     g = P^T rhs
            //
            // Solve
            //
            //     (H0 + lambda I) alpha = g
            //
            // =====================================================

            Eigen::MatrixXd H0(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);


            H0.noalias() =
                P.transpose() * AP;

            // Remove very small asymmetry caused by floating-point
            // accumulation. Eigen::LLT expects a symmetric matrix.
            H0 =
                0.5 * (H0 + H0.transpose());

            g.noalias() =
                P.transpose() * b;


            if (!H0.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix contains NaN or infinity");
            }

            if (!g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced rhs contains NaN or infinity");
            }


            // =====================================================
            // Coefficient-space ball trust-region solve
            //
            //     minimize 0.5 alpha^T H0 alpha - g^T alpha
            //
            //     subject to
            //
            //         alpha^T alpha <= trustRadius^2
            //
            // Its KKT equation is
            //
            //     (H0 + lambda I) alpha = g,
            //
            // with lambda >= 0.  If the unconstrained Cholesky
            // solution exists and is inside the ball, lambda = 0.
            // Otherwise lambda is bracketed and found by bisection.
            // =====================================================

            const double radiusSquared =
                trustRadius * trustRadius;

            const double radiusTolerance =
                1.0e-10;

            Eigen::VectorXd alpha(
                reducedSize);

            alpha.setZero();

            bool acceptedUnconstrained =
                false;

            Eigen::LLT<Eigen::MatrixXd> cholesky;

            cholesky.compute(H0);

            if (cholesky.info() == Eigen::Success)
            {
                Eigen::VectorXd candidate =
                    cholesky.solve(g);

                if (cholesky.info() == Eigen::Success &&
                    candidate.allFinite())
                {
                    const double normSquared =
                        candidate.squaredNorm();

                    if (std::isfinite(normSquared) &&
                        normSquared <=
                        radiusSquared *
                        (1.0 + radiusTolerance))
                    {
                        alpha =
                            candidate;

                        lambdaUsed =
                            0.0;

                        acceptedUnconstrained =
                            true;
                    }
                }
            }


            if (!acceptedUnconstrained)
            {
                // Obtain a scale only for the initial bracketing value.
                // The final lambda is determined by the requested radius.
                double hScale =
                    0.0;

                for (Eigen::Index i = 0;
                    i < reducedSize;
                    ++i)
                {
                    hScale =
                        std::max(
                            hScale,
                            std::abs(H0(i, i)));

                }

                double lambdaLow =
                    0.0;

                double lambdaHigh =
                    1.0e-12;

                if (hScale > 0.0)
                {
                    lambdaHigh =
                        1.0e-12 * hScale;
                }

                if (!std::isfinite(lambdaHigh) ||
                    lambdaHigh <= 0.0)
                {
                    lambdaHigh =
                        1.0e-12;
                }

                bool bracketed =
                    false;

                Eigen::VectorXd candidate(
                    reducedSize);

                for (int iteration = 0;
                    iteration < 100;
                    ++iteration)
                {
                    Eigen::MatrixXd H =
                        H0;

                    H.diagonal().array() +=
                        lambdaHigh;

                    H =
                        0.5 * (H + H.transpose());

                    cholesky.compute(H);

                    if (cholesky.info() == Eigen::Success)
                    {
                        candidate =
                            cholesky.solve(g);

                        if (cholesky.info() == Eigen::Success &&
                            candidate.allFinite())
                        {
                            const double normSquared =
                                candidate.squaredNorm();

                            if (std::isfinite(normSquared) &&
                                normSquared <= radiusSquared)
                            {
                                alpha =
                                    candidate;

                                bracketed =
                                    true;

                                break;
                            }
                        }
                    }

                    lambdaLow =
                        lambdaHigh;

                    lambdaHigh *=
                        10.0;

                    if (!std::isfinite(lambdaHigh))
                    {
                        break;
                    }
                }

                if (!bracketed)
                {
                    throw gcnew InvalidOperationException(
                        "Failed to bracket the trust-region multiplier");
                }


                // For lambda > 0, ||alpha(lambda)|| decreases
                // monotonically.  Bisect until the boundary is reached.
                for (int iteration = 0;
                    iteration < 80;
                    ++iteration)
                {
                    const double lambdaMiddle =
                        0.5 * (lambdaLow + lambdaHigh);

                    Eigen::MatrixXd H =
                        H0;

                    H.diagonal().array() +=
                        lambdaMiddle;

                    H =
                        0.5 * (H + H.transpose());

                    cholesky.compute(H);

                    if (cholesky.info() != Eigen::Success)
                    {
                        lambdaLow =
                            lambdaMiddle;

                        continue;
                    }

                    candidate =
                        cholesky.solve(g);

                    if (cholesky.info() != Eigen::Success ||
                        !candidate.allFinite())
                    {
                        lambdaLow =
                            lambdaMiddle;

                        continue;
                    }

                    const double normSquared =
                        candidate.squaredNorm();

                    if (!std::isfinite(normSquared) ||
                        normSquared > radiusSquared)
                    {
                        lambdaLow =
                            lambdaMiddle;
                    }
                    else
                    {
                        lambdaHigh =
                            lambdaMiddle;

                        alpha =
                            candidate;
                    }

                    const double relativeRadiusError =
                        std::abs(normSquared - radiusSquared) /
                        radiusSquared;

                    if (relativeRadiusError <=
                        radiusTolerance)
                    {
                        break;
                    }
                }

                lambdaUsed =
                    lambdaHigh;
            }


            if (!alpha.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Full-space direction
            //
            //     direction = P alpha
            //
            // Since P is the ORIGINAL basis,
            // alpha is directly the coefficient vector that we want.
            // =====================================================

            Eigen::VectorXd direction(
                n);

            direction.noalias() =
                P * alpha;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(n);

            parameters =
                gcnew array<double>(
                    static_cast<int>(
                        reducedSize));


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    direction(i);
            }


            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                parameters[
                    static_cast<int>(j)] =
                    alpha(j);
            }


            return ret;
        }
        static array<double>^ findSearchDirection3(
            int n,

            // =====================================================
            // Matrix A
            // CSC format
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Original basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // No orthonormalization is performed.
            // =====================================================
            array<array<double>^>^ vecs,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,

            // =====================================================
            // Independent trust radius for each basis coefficient
            //
            //     |alpha_i| <= trustRadii[i]
            //
            // =====================================================
            array<double>^ trustRadii,

            // =====================================================
            // Output:
            // coefficients alpha with respect to ORIGINAL basis P
            //
            //     direction = P alpha
            // =====================================================
            array<double>^% parameters,

            // =====================================================
            // Maximum diagonal regularization actually used
            // =====================================================
            double% lambdaUsed)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // Input validation
            // =====================================================

            if (n <= 0)
                throw gcnew ArgumentException(
                    "n must be positive");

            if (columnptr == nullptr)
                throw gcnew ArgumentNullException(
                    "columnptr");

            if (rowindices == nullptr)
                throw gcnew ArgumentNullException(
                    "rowindices");

            if (values == nullptr)
                throw gcnew ArgumentNullException(
                    "values");

            if (vecs == nullptr)
                throw gcnew ArgumentNullException(
                    "vecs");

            if (rhs == nullptr)
                throw gcnew ArgumentNullException(
                    "rhs");

            if (trustRadii == nullptr)
                throw gcnew ArgumentNullException(
                    "trustRadii");


            if (rhs->Length != n)
                throw gcnew ArgumentException(
                    "rhs has incompatible size");

            if (columnptr->Length != n + 1)
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");

            if (rowindices->Length != values->Length)
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");

            if (columnptr[0] != 0 ||
                columnptr[n] != values->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr has invalid start or end");
            }


            const Eigen::Index basisCount =
                static_cast<Eigen::Index>(
                    vecs->Length);


            if (basisCount <= 0)
                throw gcnew ArgumentException(
                    "vecs must contain at least one basis vector");


            if (trustRadii->Length !=
                static_cast<int>(basisCount))
            {
                throw gcnew ArgumentException(
                    "trustRadii must have the same length as vecs");
            }


            for (int j = 0;
                j < trustRadii->Length;
                ++j)
            {
                if (!std::isfinite(trustRadii[j]) ||
                    trustRadii[j] <= 0.0)
                {
                    throw gcnew ArgumentException(
                        "Every trust radius must be finite and positive");
                }
            }


            // =====================================================
            // Build basis matrix P
            //
            //     P = [p0 p1 ...]
            // =====================================================

            Eigen::MatrixXd P(
                n,
                basisCount);


            for (Eigen::Index j = 0;
                j < basisCount;
                ++j)
            {
                array<double>^ src =
                    vecs[
                        static_cast<int>(j)];


                if (src == nullptr ||
                    src->Length != n)
                {
                    throw gcnew ArgumentException(
                        "vecs contains an invalid vector");
                }


                for (int i = 0;
                    i < n;
                    ++i)
                {
                    if (!std::isfinite(src[i]))
                    {
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");
                    }


                    P(i, j) =
                        src[i];
                }
            }


            // =====================================================
            // Build RHS
            // =====================================================

            Eigen::VectorXd b(n);


            for (int i = 0;
                i < n;
                ++i)
            {
                if (!std::isfinite(rhs[i]))
                {
                    throw gcnew ArgumentException(
                        "rhs contains NaN or infinity");
                }


                b(i) =
                    rhs[i];
            }


            // =====================================================
            // Build sparse matrix A
            // =====================================================

            SparseMatrix64 matA(
                n,
                n);


            std::vector<Triplet64>
                tripletList;


            tripletList.reserve(
                static_cast<std::size_t>(
                    values->Length));


            for (int column = 0;
                column < n;
                ++column)
            {
                const int start =
                    columnptr[column];

                const int end =
                    columnptr[column + 1];


                if (start < 0 ||
                    end < start ||
                    end > values->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr contains an invalid range");
                }


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices[index];

                    const double value =
                        values[index];


                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "rowindices contains an invalid row index");
                    }


                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "values contains NaN or infinity");
                    }


                    tripletList.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }


            matA.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            matA.makeCompressed();


            // =====================================================
            // Compute A P
            // =====================================================

            Eigen::MatrixXd AP =
                matA * P;


            if (!AP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "A * P contains NaN or infinity");
            }


            // =====================================================
            // Reduced system
            //
            //     H = P^T A P
            //
            //     g = P^T rhs
            //
            //
            // Solve
            //
            //     min 1/2 alpha^T H alpha - g^T alpha
            //
            // subject to
            //
            //     -trustRadii[i]
            //         <= alpha_i
            //         <= trustRadii[i]
            //
            // =====================================================

            Eigen::MatrixXd H =
                P.transpose() * AP;


            // enforce symmetry
            H =
                0.5 *
                (H + H.transpose());


            Eigen::VectorXd g =
                P.transpose() * b;


            if (!H.allFinite() ||
                !g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced system contains NaN or infinity");
            }


            // =====================================================
            // Initial solution
            // =====================================================

            Eigen::VectorXd alpha =
                Eigen::VectorXd::Zero(
                    basisCount);


            // =====================================================
            // Active set status
            //
            //  0 : free
            //
            // +1 : alpha_i = +trustRadii[i]
            //
            // -1 : alpha_i = -trustRadii[i]
            // =====================================================

            std::vector<int> status(
                static_cast<std::size_t>(
                    basisCount),
                0);


            lambdaUsed =
                0.0;


            // =====================================================
            // Active-set parameters
            // =====================================================

            const int maxIterations =
                20 *
                static_cast<int>(
                    basisCount)
                + 100;


            const double boundTolerance =
                1.0e-10;


            const double kktTolerance =
                1.0e-10;


            bool converged =
                false;


            // =====================================================
            // Active-set loop
            // =====================================================

            for (int outer = 0;
                outer < maxIterations;
                ++outer)
            {
                // =================================================
                // Collect free and active variables
                // =================================================

                std::vector<int>
                    freeIndices;


                std::vector<int>
                    activeIndices;


                freeIndices.reserve(
                    static_cast<std::size_t>(
                        basisCount));


                activeIndices.reserve(
                    static_cast<std::size_t>(
                        basisCount));


                for (int i = 0;
                    i < static_cast<int>(
                        basisCount);
                    ++i)
                {
                    if (status[
                        static_cast<std::size_t>(i)]
                        == 0)
                    {
                        freeIndices.push_back(i);
                    }
                    else
                    {
                        activeIndices.push_back(i);
                    }
                }


                // =================================================
                // Solve free variables
                //
                //
                // Partition:
                //
                //     H_FF alpha_F
                //       + H_FA alpha_A
                //       = g_F
                //
                //
                // therefore
                //
                //     H_FF alpha_F
                //       =
                //     g_F - H_FA alpha_A
                //
                // =================================================

                if (!freeIndices.empty())
                {
                    const Eigen::Index freeCount =
                        static_cast<Eigen::Index>(
                            freeIndices.size());


                    Eigen::MatrixXd Hff(
                        freeCount,
                        freeCount);


                    Eigen::VectorXd gf(
                        freeCount);


                    // =============================================
                    // Construct Hff and reduced RHS
                    // =============================================

                    for (Eigen::Index ii = 0;
                        ii < freeCount;
                        ++ii)
                    {
                        const int i =
                            freeIndices[
                                static_cast<std::size_t>(
                                    ii)];


                        double rhsReduced =
                            g(i);


                        for (int j :
                        activeIndices)
                        {
                            rhsReduced -=
                                H(i, j) *
                                alpha(j);
                        }


                        gf(ii) =
                            rhsReduced;


                        for (Eigen::Index jj = 0;
                            jj < freeCount;
                            ++jj)
                        {
                            const int j =
                                freeIndices[
                                    static_cast<std::size_t>(
                                        jj)];


                            Hff(ii, jj) =
                                H(i, j);
                        }
                    }


                    // enforce symmetry
                    Hff =
                        0.5 *
                        (Hff +
                            Hff.transpose());


                    if (!Hff.allFinite() ||
                        !gf.allFinite())
                    {
                        throw gcnew
                            InvalidOperationException(
                                "Free reduced system contains NaN or infinity");
                    }


                    // =============================================
                    // Solve
                    //
                    //     Hff x = gf
                    //
                    //
                    // Try LLT first.
                    //
                    // If Hff is numerically singular,
                    // add the smallest diagonal damping needed.
                    // =============================================

                    Eigen::LLT<
                        Eigen::MatrixXd>
                        cholesky;


                    Eigen::VectorXd solution(
                        freeCount);


                    bool solved =
                        false;


                    double hScale =
                        0.0;


                    for (Eigen::Index i = 0;
                        i < freeCount;
                        ++i)
                    {
                        hScale =
                            std::max(
                                hScale,
                                std::abs(
                                    Hff(i, i)));
                    }


                    double lambda =
                        0.0;


                    for (int attempt = 0;
                        attempt < 30;
                        ++attempt)
                    {
                        Eigen::MatrixXd Hsolve =
                            Hff;


                        if (lambda > 0.0)
                        {
                            Hsolve
                                .diagonal()
                                .array()
                                += lambda;
                        }


                        cholesky.compute(
                            Hsolve);


                        if (cholesky.info() ==
                            Eigen::Success)
                        {
                            solution =
                                cholesky.solve(
                                    gf);


                            if (cholesky.info() ==
                                Eigen::Success &&
                                solution.allFinite())
                            {
                                solved =
                                    true;

                                if (lambda > lambdaUsed)
                                    lambdaUsed = lambda;
//                                    std::max(
  //                                      lambdaUsed,
    //                                    lambda);


                                break;
                            }
                        }


                        // =========================================
                        // Increase damping only if necessary
                        // =========================================

                        if (lambda == 0.0)
                        {
                            lambda =
                                std::max(
                                    1.0e-14,
                                    1.0e-12 *
                                    std::max(
                                        1.0,
                                        hScale));
                        }
                        else
                        {
                            lambda *=
                                10.0;
                        }
                    }


                    if (!solved)
                    {
                        throw gcnew
                            InvalidOperationException(
                                "Failed to solve reduced active-set system");
                    }


                    // =============================================
                    // Insert solution into alpha
                    // =============================================

                    for (Eigen::Index ii = 0;
                        ii < freeCount;
                        ++ii)
                    {
                        const int index =
                            freeIndices[
                                static_cast<std::size_t>(
                                    ii)];


                        alpha(index) =
                            solution(ii);
                    }
                }


                // =================================================
                // Find trust-radius violations
                //
                //
                // We activate ONE bound at a time.
                //
                // Choose the coefficient with the largest
                // relative violation.
                // =================================================

                int violatingIndex =
                    -1;


                int violatingSide =
                    0;


                double largestViolation =
                    0.0;


                for (int i = 0;
                    i < static_cast<int>(
                        basisCount);
                    ++i)
                {
                    if (status[
                        static_cast<std::size_t>(i)]
                        != 0)
                    {
                        continue;
                    }


                    const double radius =
                        trustRadii[i];


                    // =============================================
                    // upper bound violation
                    // =============================================

                    if (alpha(i) >
                        radius *
                        (1.0 +
                            boundTolerance))
                    {
                        const double violation =
                            alpha(i) /
                            radius
                            - 1.0;


                        if (violation >
                            largestViolation)
                        {
                            largestViolation =
                                violation;


                            violatingIndex =
                                i;


                            violatingSide =
                                +1;
                        }
                    }


                    // =============================================
                    // lower bound violation
                    // =============================================

                    if (alpha(i) <
                        -radius *
                        (1.0 +
                            boundTolerance))
                    {
                        const double violation =
                            -alpha(i) /
                            radius
                            - 1.0;


                        if (violation >
                            largestViolation)
                        {
                            largestViolation =
                                violation;


                            violatingIndex =
                                i;


                            violatingSide =
                                -1;
                        }
                    }
                }


                // =================================================
                // Activate violated bound
                // =================================================

                if (violatingIndex >= 0)
                {
                    const double radius =
                        trustRadii[
                            violatingIndex];


                    if (violatingSide > 0)
                    {
                        alpha(
                            violatingIndex)
                            =
                            radius;


                        status[
                            static_cast<std::size_t>(
                                violatingIndex)]
                            =
                            +1;
                    }
                    else
                    {
                        alpha(
                            violatingIndex)
                            =
                            -radius;


                        status[
                            static_cast<std::size_t>(
                                violatingIndex)]
                            =
                            -1;
                    }


                    continue;
                }


                // =================================================
                // Current alpha is feasible.
                //
                // Check KKT conditions.
                //
                //
                // Objective:
                //
                //     f(alpha)
                //       =
                //     1/2 alpha^T H alpha
                //       - g^T alpha
                //
                //
                // Gradient:
                //
                //     grad = H alpha - g
                //
                //
                // At upper bound:
                //
                //     grad_i <= 0
                //
                //
                // At lower bound:
                //
                //     grad_i >= 0
                //
                // =================================================

                Eigen::VectorXd gradient =
                    H * alpha
                    - g;


                if (!gradient.allFinite())
                {
                    throw gcnew
                        InvalidOperationException(
                            "Gradient contains NaN or infinity");
                }


                int releaseIndex =
                    -1;


                double largestKktViolation =
                    0.0;


                for (int i = 0;
                    i < static_cast<int>(
                        basisCount);
                    ++i)
                {
                    const int currentStatus =
                        status[
                            static_cast<std::size_t>(i)];


                    // =============================================
                    // Upper bound
                    //
                    // Need:
                    //
                    //     gradient(i) <= 0
                    // =============================================

                    if (currentStatus == +1)
                    {
                        const double violation =
                            gradient(i);


                        if (violation >
                            kktTolerance &&
                            violation >
                            largestKktViolation)
                        {
                            largestKktViolation =
                                violation;


                            releaseIndex =
                                i;
                        }
                    }


                    // =============================================
                    // Lower bound
                    //
                    // Need:
                    //
                    //     gradient(i) >= 0
                    // =============================================

                    else if (
                        currentStatus == -1)
                    {
                        const double violation =
                            -gradient(i);


                        if (violation >
                            kktTolerance &&
                            violation >
                            largestKktViolation)
                        {
                            largestKktViolation =
                                violation;


                            releaseIndex =
                                i;
                        }
                    }
                }


                // =================================================
                // Release an incorrectly active bound
                // =================================================

                if (releaseIndex >= 0)
                {
                    status[
                        static_cast<std::size_t>(
                            releaseIndex)]
                        =
                        0;


                    continue;
                }


                // =================================================
                // Feasible and KKT satisfied
                // =================================================

                converged =
                    true;


                break;
            }


            // =====================================================
            // Active-set convergence check
            // =====================================================

            if (!converged)
            {
                throw gcnew
                    InvalidOperationException(
                        "Box trust-region active-set solver did not converge");
            }


            // =====================================================
            // Numerical safety clamp
            //
            // This only removes tiny floating point violations.
            // =====================================================

            for (int i = 0;
                i < static_cast<int>(
                    basisCount);
                ++i)
            {
                const double radius =
                    trustRadii[i];


                if (alpha(i) >
                    radius)
                {
                    alpha(i) =
                        radius;
                }


                if (alpha(i) <
                    -radius)
                {
                    alpha(i) =
                        -radius;
                }
            }


            // =====================================================
            // Full-space direction
            //
            //     direction = P alpha
            // =====================================================

            Eigen::VectorXd direction =
                P * alpha;


            if (!alpha.allFinite() ||
                !direction.allFinite())
            {
                throw gcnew
                    InvalidOperationException(
                        "Solution contains NaN or infinity");
            }


            // =====================================================
            // Output coefficients
            // =====================================================

            parameters =
                gcnew array<double>(
                    static_cast<int>(
                        basisCount));


            for (Eigen::Index j = 0;
                j < basisCount;
                ++j)
            {
                parameters[
                    static_cast<int>(j)]
                    =
                    alpha(j);
            }


            // =====================================================
            // Output full-space direction
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(
                    n);


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    direction(i);
            }


            return ret;
        }
        static array<double>^ findSearchDirection4(
            int n,

            // =====================================================
            // Matrix A
            // CSC format
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Matrix B (secondary)
            // CSC format
            // =====================================================
            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            // =====================================================
            // Original basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // NO orthonormalization is performed.
            // =====================================================
            array<array<double>^>^ vecs,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,

            // Secondary RHS
            array<double>^ rhs2,

            // =====================================================
            // Mixing weight
            //
            //     M = A + mu B
            //     b = rhs + mu rhs2
            // =====================================================
            double mu,

            // =====================================================
            // Radius of the coefficient-space ball constraint
            //
            //     ||alpha|| <= trustRadius
            // =====================================================
            double trustRadius,

            // Index of the full-problem GN direction in vecs.
            int fullDirectionIndex,

            // Relative full-direction anchor weight.
            //
            //     gamma = anchorRelativeWeight * max_i |H0_ii|
            //
            // The added objective is
            //
            //     0.5 * gamma * ||alpha - e_full||^2.
            //
            // A typical initial value is 1.0e-6.
            double anchorRelativeWeight,

            // Relative eigenvalue threshold for the unconstrained
            // reduced TSVD solve.  A typical value is 1.0e-10.
            double tsvdThreshold,

            // =====================================================
            // Output:
            //
            // coefficients alpha with respect to ORIGINAL basis P
            //
            //     direction = P alpha
            // =====================================================
            array<double>^% parameters,

            // LM multiplier selected by the trust-region solve.
            // Zero means that the unconstrained step was accepted.
            double% lambdaUsed)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (columnptr == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "columnptr");
            }

            if (rowindices == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rowindices");
            }

            if (values == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "values");
            }

            if (columnptr2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "columnptr2");
            }

            if (rowindices2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rowindices2");
            }

            if (values2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "values2");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (rhs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs");
            }

            if (rhs2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs2");
            }

            if (!std::isfinite(mu) ||
                mu < 0.0)
            {
                throw gcnew ArgumentException(
                    "mu must be finite and nonnegative");
            }

            if (!std::isfinite(trustRadius) ||
                trustRadius <= 0.0)
            {
                throw gcnew ArgumentException(
                    "trustRadius must be finite and positive");
            }

            if (!std::isfinite(anchorRelativeWeight) ||
                anchorRelativeWeight < 0.0)
            {
                throw gcnew ArgumentException(
                    "anchorRelativeWeight must be finite and nonnegative");
            }

            if (!std::isfinite(tsvdThreshold) ||
                tsvdThreshold <= 0.0)
            {
                throw gcnew ArgumentException(
                    "tsvdThreshold must be finite and positive");
            }

            if (rhs->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs has incompatible size");
            }

            if (rhs2->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs2 has incompatible size");
            }

            if (columnptr->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");
            }

            if (rowindices->Length != values->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");
            }

            if (columnptr2->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr2 must have length n + 1");
            }

            if (rowindices2->Length != values2->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices2 and values2 have incompatible sizes");
            }

            if (columnptr[0] != 0 ||
                columnptr[n] != values->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr has invalid start or end");
            }

            if (columnptr2[0] != 0 ||
                columnptr2[n] != values2->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr2 has invalid start or end");
            }


            // =====================================================
            // Number of basis vectors
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    vecs->Length);

            // =====================================================
            // Empty basis
            // =====================================================

            if (reducedSize == 0)
            {
                parameters =
                    gcnew array<double>(0);

                lambdaUsed =
                    0.0;

                return gcnew array<double>(n);
            }

            if (fullDirectionIndex < 0 ||
                fullDirectionIndex >= vecs->Length)
            {
                throw gcnew ArgumentOutOfRangeException(
                    "fullDirectionIndex");
            }


            // =====================================================
            // Build original dense basis
            //
            //     P = [p0 p1 ...]
            //
            // NO QR / orthonormalization.
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                array<double>^ src =
                    vecs[
                        static_cast<int>(j)];

                if (src == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a null vector");
                }

                if (src->Length != n)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a vector with incompatible size");
                }

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double value =
                        src[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");
                    }

                    P(i, j) =
                        value;
                }
            }


            if (!P.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Basis matrix contains NaN or infinity");
            }


            // =====================================================
            // Build mixed sparse matrix
            //
            //     M = A + mu B
            // =====================================================

            SparseMatrix64 matM(
                n,
                n);

            std::vector<Triplet64> tripletsM;

            tripletsM.reserve(
                static_cast<std::size_t>(
                    values->Length + values2->Length));


            for (int column = 0;
                column < n;
                ++column)
            {
                const int start =
                    columnptr[column];

                const int end =
                    columnptr[column + 1];

                if (start < 0 ||
                    end < start ||
                    end > values->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr contains an invalid range");
                }


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices[index];

                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "rowindices contains an invalid row index");
                    }


                    const double value =
                        values[index];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "values contains NaN or infinity");
                    }


                    tripletsM.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }


            for (int column = 0;
                column < n;
                ++column)
            {
                const int start =
                    columnptr2[column];

                const int end =
                    columnptr2[column + 1];

                if (start < 0 ||
                    end < start ||
                    end > values2->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 contains an invalid range");
                }


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices2[index];

                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "rowindices2 contains an invalid row index");
                    }


                    const double value =
                        values2[index];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "values2 contains NaN or infinity");
                    }


                    tripletsM.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        mu * value);
                }
            }


            // Duplicate entries are summed here, so overlapping
            // nonzeros of A and B are combined correctly.
            matM.setFromTriplets(
                tripletsM.begin(),
                tripletsM.end());

            matM.makeCompressed();


            // =====================================================
            // Mixed RHS
            //
            //     b = rhs + mu rhs2
            // =====================================================

            Eigen::VectorXd b(
                n);

            for (int i = 0;
                i < n;
                ++i)
            {
                const double value1 =
                    rhs[i];

                const double value2 =
                    rhs2[i];

                if (!std::isfinite(value1))
                {
                    throw gcnew ArgumentException(
                        "rhs contains NaN or infinity");
                }

                if (!std::isfinite(value2))
                {
                    throw gcnew ArgumentException(
                        "rhs2 contains NaN or infinity");
                }

                b(i) =
                    value1 + mu * value2;
            }


            // =====================================================
            // Compute
            //
            //     MP = (A + mu B) P
            // =====================================================

            Eigen::MatrixXd MP(
                n,
                reducedSize);

            MP.noalias() =
                matM * P;


            if (!MP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "(A + mu B) * P contains NaN or infinity");
            }


            // =====================================================
            // Reduced mixed system
            //
            //     H0 = P^T (A + mu B) P
            //
            //     g = P^T (rhs + mu rhs2)
            //
            // Full-direction anchor
            //
            //     gamma = anchorRelativeWeight * max_i |H0_ii|
            //
            //     H = H0 + gamma I
            //     gAnchor = g + gamma e_full
            //
            // =====================================================

            Eigen::MatrixXd H0(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);


            H0.noalias() =
                P.transpose() * MP;

            // Remove very small asymmetry caused by floating-point
            // accumulation. Eigen::LLT expects a symmetric matrix.
            H0 =
                0.5 * (H0 + H0.transpose());

            g.noalias() =
                P.transpose() * b;


            if (!H0.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix contains NaN or infinity");
            }

            if (!g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced rhs contains NaN or infinity");
            }


            double hScale =
                0.0;

            for (Eigen::Index i = 0;
                i < reducedSize;
                ++i)
            {
                hScale =
                    std::max(
                        hScale,
                        std::abs(H0(i, i)));
            }

            if (!std::isfinite(hScale))
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix scale is not finite");
            }

            const double gamma =
                anchorRelativeWeight * hScale;

            Eigen::MatrixXd HAnchor =
                H0;

            HAnchor.diagonal().array() +=
                gamma;

            HAnchor =
                0.5 * (HAnchor + HAnchor.transpose());

            Eigen::VectorXd gAnchor =
                g;

            gAnchor(
                static_cast<Eigen::Index>(
                    fullDirectionIndex)) +=
                gamma;

            if (!HAnchor.allFinite() ||
                !gAnchor.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Anchored reduced system contains NaN or infinity");
            }


            // =====================================================
            // Coefficient-space ball trust-region solve
            //
            //     minimize
            //
            //         0.5 alpha^T H0 alpha - g^T alpha
            //       + 0.5 gamma ||alpha - e_full||^2
            //
            //     subject to
            //
            //         alpha^T alpha <= trustRadius^2
            //
            // Its KKT equation is
            //
            //     (H0 + (gamma + lambda) I) alpha
            //         = g + gamma e_full,
            //
            // with lambda >= 0.  The unconstrained anchored system is
            // solved first by a symmetric TSVD.  If that solution is
            // inside the ball, lambda = 0.  Otherwise lambda is
            // bracketed and found by bisection.
            // =====================================================

            const double radiusSquared =
                trustRadius * trustRadius;

            const double radiusTolerance =
                1.0e-10;

            Eigen::VectorXd alpha(
                reducedSize);

            alpha.setZero();

            bool acceptedUnconstrained =
                false;

            // =====================================================
            // Unconstrained anchored solution by symmetric TSVD
            // =====================================================

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>
                eigenSolver;

            eigenSolver.compute(HAnchor);

            if (eigenSolver.info() != Eigen::Success)
            {
                throw gcnew InvalidOperationException(
                    "Eigenvalue decomposition of the reduced system failed");
            }

            const Eigen::VectorXd eigenvalues =
                eigenSolver.eigenvalues();

            const Eigen::MatrixXd eigenvectors =
                eigenSolver.eigenvectors();

            double maximumEigenvalue =
                0.0;

            for (Eigen::Index i = 0;
                i < reducedSize;
                ++i)
            {
                maximumEigenvalue =
                    std::max(
                        maximumEigenvalue,
                        std::abs(eigenvalues(i)));
            }

            if (!std::isfinite(maximumEigenvalue))
            {
                throw gcnew InvalidOperationException(
                    "Reduced eigenvalues contain NaN or infinity");
            }

            Eigen::VectorXd candidate(
                reducedSize);

            candidate.setZero();

            if (maximumEigenvalue > 0.0)
            {
                const double cutoff =
                    tsvdThreshold * maximumEigenvalue;

                const Eigen::VectorXd projectedRhs =
                    eigenvectors.transpose() * gAnchor;

                for (Eigen::Index i = 0;
                    i < reducedSize;
                    ++i)
                {
                    const double eigenvalue =
                        eigenvalues(i);

                    if (eigenvalue > cutoff)
                    {
                        candidate +=
                            eigenvectors.col(i) *
                            (projectedRhs(i) / eigenvalue);
                    }
                }
            }

            if (candidate.allFinite())
            {
                const double normSquared =
                    candidate.squaredNorm();

                if (std::isfinite(normSquared) &&
                    normSquared <=
                    radiusSquared *
                    (1.0 + radiusTolerance))
                {
                    alpha =
                        candidate;

                    lambdaUsed =
                        0.0;

                    acceptedUnconstrained =
                        true;
                }
            }

            Eigen::LLT<Eigen::MatrixXd> cholesky;


            if (!acceptedUnconstrained)
            {
                double lambdaLow =
                    0.0;

                double lambdaHigh =
                    1.0e-12;

                if (hScale > 0.0)
                {
                    lambdaHigh =
                        1.0e-12 * hScale;
                }

                if (!std::isfinite(lambdaHigh) ||
                    lambdaHigh <= 0.0)
                {
                    lambdaHigh =
                        1.0e-12;
                }

                bool bracketed =
                    false;

                Eigen::VectorXd candidate(
                    reducedSize);

                for (int iteration = 0;
                    iteration < 100;
                    ++iteration)
                {
                    Eigen::MatrixXd H =
                        HAnchor;

                    H.diagonal().array() +=
                        lambdaHigh;

                    H =
                        0.5 * (H + H.transpose());

                    cholesky.compute(H);

                    if (cholesky.info() == Eigen::Success)
                    {
                        candidate =
                            cholesky.solve(gAnchor);

                        if (cholesky.info() == Eigen::Success &&
                            candidate.allFinite())
                        {
                            const double normSquared =
                                candidate.squaredNorm();

                            if (std::isfinite(normSquared) &&
                                normSquared <= radiusSquared)
                            {
                                alpha =
                                    candidate;

                                bracketed =
                                    true;

                                break;
                            }
                        }
                    }

                    lambdaLow =
                        lambdaHigh;

                    lambdaHigh *=
                        10.0;

                    if (!std::isfinite(lambdaHigh))
                    {
                        break;
                    }
                }

                if (!bracketed)
                {
                    throw gcnew InvalidOperationException(
                        "Failed to bracket the trust-region multiplier");
                }


                // For lambda > 0, ||alpha(lambda)|| decreases
                // monotonically.  Bisect until the boundary is reached.
                for (int iteration = 0;
                    iteration < 80;
                    ++iteration)
                {
                    const double lambdaMiddle =
                        0.5 * (lambdaLow + lambdaHigh);

                    Eigen::MatrixXd H =
                        HAnchor;

                    H.diagonal().array() +=
                        lambdaMiddle;

                    H =
                        0.5 * (H + H.transpose());

                    cholesky.compute(H);

                    if (cholesky.info() != Eigen::Success)
                    {
                        lambdaLow =
                            lambdaMiddle;

                        continue;
                    }

                    candidate =
                        cholesky.solve(gAnchor);

                    if (cholesky.info() != Eigen::Success ||
                        !candidate.allFinite())
                    {
                        lambdaLow =
                            lambdaMiddle;

                        continue;
                    }

                    const double normSquared =
                        candidate.squaredNorm();

                    if (!std::isfinite(normSquared) ||
                        normSquared > radiusSquared)
                    {
                        lambdaLow =
                            lambdaMiddle;
                    }
                    else
                    {
                        lambdaHigh =
                            lambdaMiddle;

                        alpha =
                            candidate;
                    }

                    const double relativeRadiusError =
                        std::abs(normSquared - radiusSquared) /
                        radiusSquared;

                    if (relativeRadiusError <=
                        radiusTolerance)
                    {
                        break;
                    }
                }

                lambdaUsed =
                    lambdaHigh;
            }


            if (!alpha.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Full-space direction
            //
            //     direction = P alpha
            //
            // Since P is the ORIGINAL basis,
            // alpha is directly the coefficient vector that we want.
            // =====================================================

            Eigen::VectorXd direction(
                n);

            direction.noalias() =
                P * alpha;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(n);

            parameters =
                gcnew array<double>(
                    static_cast<int>(
                        reducedSize));


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    direction(i);
            }


            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                parameters[
                    static_cast<int>(j)] =
                    alpha(j);
            }


            return ret;
        }

        static array<array<double>^>^ GenerateRandomVectors(
            int n,
            int N,
            int randomSeed)
        {
            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (N < 0)
            {
                throw gcnew ArgumentException(
                    "N must be non-negative");
            }

            array<array<double>^>^ ret =
                gcnew array<array<double>^>(N);

            std::mt19937_64 generator(
                static_cast<
                std::mt19937_64::result_type>(
                    randomSeed));

            std::uniform_int_distribution<int>
                distribution(0, 1);

            // ±1/sqrt(n)なので、各ベクトルのノルムは
            // 丸め誤差を除いて1になる。
            const double scale =
                1.0 /
                std::sqrt(
                    static_cast<double>(n));

            for (int k = 0;
                k < N;
                ++k)
            {
                array<double>^ v =
                    gcnew array<double>(n);

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    if (distribution(generator) == 0)
                    {
                        v[i] = -scale;
                    }
                    else
                    {
                        v[i] = scale;
                    }
                }

                ret[k] = v;
            }

            return ret;
        }
        static int countEigenvaluesBelow(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,
            double threshold)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;

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
                    tripletList.push_back(Triplet64(
                        row,
                        i,
                        val));
                }
            }

            // B = A - threshold * I
            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Triplet64(
                        i,
                        i,
                        -threshold));
            }

            SparseMatrix64 mat(n, n);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();

            // IMPORTANT:
            // symmetric indefinite factorization
            Eigen::PardisoLDLT<SparseMatrix64> solver;

            solver.compute(mat);

            if (solver.info() != Eigen::Success)
                return -1;

            auto& iparm = solver.pardisoParameterArray();

            // number of negative eigenvalues of A - threshold I
            //
            // = number of eigenvalues of A below threshold
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
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;
            SparseMatrix64 mat;
            Eigen::VectorXd b(n);
            Eigen::VectorXd _r(n);

            // =====================================================
            // RHS
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                b(i) = rhs[i];
            }


            // =====================================================
            // Build matrix A
            //
            // Also find
            //
            //     matrixScale = max_i |A_ii|
            //
            // so that epsilon is relative.
            // =====================================================

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

                    tripletList.push_back(
                        Triplet64(
                            row,
                            i,
                            val));

                    if (row == i)
                    {
                        double a = std::abs(val);

                        if (a > matrixScale)
                            matrixScale = a;
                    }
                }
            }


            // =====================================================
            // Relative regularization
            //
            //     A + epsilon * matrixScale * I
            //
            // epsilon is therefore dimensionless.
            // =====================================================

            if (matrixScale == 0.0)
                matrixScale = 1.0;

            double salt = epsilon *matrixScale;

            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Triplet64(
                        i,
                        i,
                        salt));
            }


            // =====================================================
            // Solve
            // =====================================================

            mat.resize(n, n);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();

            
            // =====================================================
            Eigen::PardisoLU<SparseMatrix64> solver;

            solver.compute(mat);
            for (int i = 0; i < 50; i++)
            {
                if (solver.info() != Eigen::Success) {
                    double _salt = salt * std::pow(10.0, i + 1);

                    for (int j = 0; j < n; j++)
                    {
                        mat.coeffRef(j, j) += _salt;
                    }
                    solver.compute(mat);
                }
                else {
                    break;
                }

            }
			_r = solver.solve(b);
            // =====================================================
            // Return
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                ret[i] = _r(i);
            }

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
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;

            // =====================================================
            // Check
            // =====================================================

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


            // =====================================================
            // Build matrix A
            // =====================================================

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

                    tripletList.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(i),
                            val));

                    if (row == i)
                    {
                        double a = std::abs(val);

                        if (a > matrixScale)
                            matrixScale = a;
                    }
                }
            }


            // =====================================================
            // Relative regularization
            //
            // A + epsilon * matrixScale * I
            // =====================================================

            if (matrixScale == 0.0)
                matrixScale = 1.0;

            double salt = epsilon * matrixScale;

            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Triplet64(
                        static_cast<int64_t>(i),
                        static_cast<int64_t>(i),
                        salt));
            }


            // =====================================================
            // Assemble sparse matrix
            // =====================================================

            mat.resize(n, n);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();


            // =====================================================
            // Factorize ONCE
            // =====================================================

            Eigen::PardisoLU<SparseMatrix64> solver;

            solver.compute(mat);

            // =====================================================
            // If factorization fails, progressively add diagonal salt
            // =====================================================

            for (int i = 0; i < 50; i++)
            {
                if (solver.info() == Eigen::Success)
                    break;

                double _salt = salt * std::pow(10.0, i + 1);

                for (int j = 0; j < n; j++)
                {
                    mat.coeffRef(j, j) += _salt;
                }

                mat.makeCompressed();

                solver.compute(mat);
            }


            if (solver.info() != Eigen::Success)
                return 0;


            // =====================================================
            // Reuse factorization for all RHS
            // =====================================================

            Eigen::VectorXd b(n);

            for (int k = 0; k < nrhs; k++)
            {
                // ---------------------------------------------
                // Copy RHS
                // ---------------------------------------------

                for (int i = 0; i < n; i++)
                {
                    b(i) = rhs[k][i];
                }

                // ---------------------------------------------
                // Solve using SAME factorization
                // ---------------------------------------------

                Eigen::VectorXd x = solver.solve(b);

                if (solver.info() != Eigen::Success)
                    return 0;

                // ---------------------------------------------
                // Return
                // ---------------------------------------------

                for (int i = 0; i < n; i++)
                {
                    ret[k][i] = x(i);
                }
            }

            return 1;
        }
        static int solve_LEAST_SQUARES(
            int m,                              // rows of J
            int n,                              // cols of J

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            array<double>^ rhs,                 // size m
            array<double>^ ret,                 // size n

            double epsilon,
            int numthreads)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;

            // ============================================================
            // Build J
            // ============================================================

            std::vector<Triplet64> tripletList;
            tripletList.reserve(values->Length);

            for (int col = 0; col < n; col++)
            {
                int start = columnptr[col];
                int end = columnptr[col + 1];

                for (int k = start; k < end; k++)
                {
                    int row = rowindices[k];
                    double val = values[k];
                    tripletList.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(col),
                            val));
                }
            }

            SparseMatrix64 J(m, n);

            J.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            J.makeCompressed();


            // ============================================================
            // b
            // ============================================================

            Eigen::VectorXd b(m);

            for (int i = 0; i < m; i++)
                b(i) = rhs[i];


            Eigen::VectorXd x(n);


            // ============================================================
            // CASE 1:
            //
            // m >= n
            //
            // Least squares:
            //
            //   min ||J x - b||^2
            //
            // solve
            //
            //   (J^T J + lambda I)x = J^T b
            //
            // ============================================================

            if (m >= n)
            {
                SparseMatrix64 JTJ =
                    J.transpose() * J;

                JTJ.makeCompressed();

                Eigen::VectorXd JTb =
                    J.transpose() * b;


                // --------------------------------------------------------
                // Relative LM scale
                //
                // max diag(J^T J)
                // --------------------------------------------------------

                double matrixScale = 0.0;

                for (int i = 0; i < n; i++)
                {
                    double a =
                        std::abs(JTJ.coeff(i, i));

                    if (a > matrixScale)
                        matrixScale = a;
                }

                if (matrixScale == 0.0)
                    matrixScale = 1.0;


                double salt =
                    epsilon * matrixScale;


                // ========================================================
                // Try increasing LM regularization if factorization fails
                // ========================================================

                Eigen::PardisoLLT<SparseMatrix64> solver;

                bool success = false;

                for (int attempt = 0; attempt < 50; attempt++)
                {
                    double currentSalt =
                        salt * std::pow(10.0, attempt);

                    SparseMatrix64 A = JTJ;

                    for (int i = 0; i < n; i++)
                    {
                        A.coeffRef(i, i) += currentSalt;
                    }

                    A.makeCompressed();

                    solver.compute(A);

                    if (solver.info() == Eigen::Success)
                    {
                        x = solver.solve(JTb);

                        if (solver.info() == Eigen::Success)
                        {
                            success = true;
                            break;
                        }
                    }
                }

                if (!success)
                    return 0;
            }


            // ============================================================
            // CASE 2:
            //
            // m < n
            //
            // Minimum-norm / damped minimum-norm solution:
            //
            //   (J J^T + lambda I)y = b
            //
            //   x = J^T y
            //
            //
            // lambda -> 0 gives the Moore-Penrose
            // minimum-norm solution.
            //
            // ============================================================

            else
            {
                SparseMatrix64 JJT =
                    J * J.transpose();

                JJT.makeCompressed();


                // --------------------------------------------------------
                // Relative LM scale
                //
                // max diag(J J^T)
                // --------------------------------------------------------

                double matrixScale = 0.0;

                for (int i = 0; i < m; i++)
                {
                    double a =
                        std::abs(JJT.coeff(i, i));

                    if (a > matrixScale)
                        matrixScale = a;
                }

                if (matrixScale == 0.0)
                    matrixScale = 1.0;


                double salt =
                    epsilon * matrixScale;


                Eigen::PardisoLLT<SparseMatrix64> solver;

                Eigen::VectorXd y(m);

                bool success = false;

                for (int attempt = 0; attempt < 50; attempt++)
                {
                    double currentSalt =
                        salt * std::pow(10.0, attempt);

                    SparseMatrix64 A = JJT;

                    for (int i = 0; i < m; i++)
                    {
                        A.coeffRef(i, i) += currentSalt;
                    }

                    A.makeCompressed();

                    solver.compute(A);

                    if (solver.info() == Eigen::Success)
                    {
                        y = solver.solve(b);

                        if (solver.info() == Eigen::Success)
                        {
                            x = J.transpose() * y;

                            success = true;
                            break;
                        }
                    }
                }

                if (!success)
                    return 0;
            }


            // ============================================================
            // Return
            // ============================================================

            for (int i = 0; i < n; i++)
            {
                ret[i] = x(i);
            }

            return 1;
        }
        static int solve_CHOLESKY2(
            int n,

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            array<double>^ rhs,
            array<double>^ rhs2,

            array<double>^ ret,

            double matrixCoeff1,   // a : A
            double matrixCoeff2,   // b : B
            double rhsCoeff1,      // c : rhs
            double rhsCoeff2,      // d : rhs2

            double epsilon,
            int numthreads)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // System
            //
            //   (a A + b B + epsilon * matrixScale * I) x
            //
            //       = c rhs + d rhs2
            //
            // Primary / Secondary are NOT automatically rescaled.
            //
            // =====================================================


            // =====================================================
            // RHS
            // =====================================================

            Eigen::VectorXd b(n);

            for (int i = 0; i < n; i++)
            {
                b(i) =
                    rhsCoeff1 * rhs[i]
                    +
                    rhsCoeff2 * rhs2[i];
            }


            // =====================================================
            // Build matrix
            // =====================================================

            std::vector<Triplet64> tripletList;

            tripletList.reserve(
                values->Length
                + values2->Length
                + n);


            double matrixScale = 0.0;


            // =====================================================
            // Primary : matrixCoeff1 * A
            // =====================================================

            for (int col = 0; col < n; col++)
            {
                int start = columnptr[col];
                int end = columnptr[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices[j];

                    double val =
                        matrixCoeff1 * values[j];

                    tripletList.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(col),
                            val));

                    if (row == col)
                    {
                        matrixScale =
                            std::max(
                                matrixScale,
                                std::abs(val));
                    }
                }
            }


            // =====================================================
            // Secondary : matrixCoeff2 * B
            // =====================================================

            for (int col = 0; col < n; col++)
            {
                int start = columnptr2[col];
                int end = columnptr2[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices2[j];

                    double val =
                        matrixCoeff2 * values2[j];

                    tripletList.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(col),
                            val));

                    if (row == col)
                    {
                        matrixScale =
                            std::max(
                                matrixScale,
                                std::abs(val));
                    }
                }
            }


            // =====================================================
            // Relative epsilon
            //
            //   salt = epsilon * max |diagonal(aA + bB)|
            //
            // =====================================================

            if (!(matrixScale > 0.0))
                matrixScale = 1.0;

            const double salt =
                epsilon;// *matrixScale;


            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Triplet64(
                        static_cast<int64_t>(i),
                        static_cast<int64_t>(i),
                        salt));
            }


            // =====================================================
            // Sparse matrix
            // =====================================================

            SparseMatrix64 mat(n, n);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();


            // =====================================================
            // Solve
            // =====================================================
			Eigen::PardisoLLT<SparseMatrix64> solver;
            solver.compute(mat);
			for (int i = 0; i < 40; i++)
			{
                if (solver.info() != Eigen::Success) {
                    double _salt = salt*std::pow(10.0, i + 1);
					
					for (int j = 0; j < n; j++)
					{
						mat.coeffRef(j, j) += _salt;
					}
                    solver.compute(mat);
                }
                else {
                    break;
                }

            }
            Eigen::VectorXd solution = solver.solve(b);


            // =====================================================
            // Return
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                ret[i] = solution(i);
            }

            return 1;
        }
        static array<double>^ solve_dense(
            array<double, 2>^ A,
            array<double, 1>^ b)
        {
            int n0 = A->GetLength(0);
            int n1 = A->GetLength(1);

            if (n0 != n1)
                throw gcnew ArgumentException("A must be square.");

            if (b->Length != n0)
                throw gcnew ArgumentException("b size must match A.");

            int n = n0;

            // =====================================================
            // Copy managed arrays -> Eigen
            // =====================================================

            Eigen::MatrixXd mat(n, n);
            Eigen::VectorXd rhs(n);

            for (int i = 0; i < n; i++)
            {
                rhs(i) = b[i];

                for (int j = 0; j < n; j++)
                {
                    mat(i, j) = A[i, j];
                }
            }

            // =====================================================
            // Solve
            //
            //     A x = b
            // =====================================================

            Eigen::PartialPivLU<Eigen::MatrixXd> solver(mat);

            Eigen::VectorXd x = solver.solve(rhs);

            // =====================================================
            // Copy Eigen -> managed array
            // =====================================================

            array<double>^ ret = gcnew array<double>(n);

            for (int i = 0; i < n; i++)
            {
                ret[i] = x(i);
            }

            return ret;
        }        static int solve_NULLSPACE(
            int n,

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,       // A

            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,      // B

            array<double>^ rhs2,

            array<double>^ ret,

            double mu,
            double epsilon,
            int numIterations,
            int numthreads)
        {
            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet64 =
                Eigen::Triplet<
                double,
                int64_t>;


            // ============================================================
            // We repeatedly apply
            //
            //   x0      = (A + mu B + epsilon I)^(-1) mu rhs2
            //
            //   x_{k+1} = (A + mu B + epsilon I)^(-1) mu B x_k
            //
            // No normalization.
            // ============================================================


            // ============================================================
            // Build A + mu B + epsilon I
            // ============================================================

            std::vector<Triplet64> tripletList;

            tripletList.reserve(
                values->Length
                + values2->Length
                + n);


            // A
            for (int col = 0; col < n; col++)
            {
                int start = columnptr[col];
                int end = columnptr[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices[j];
                    double val = values[j];
                    tripletList.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(col),
                             val));
                }
            }


            // mu B
            for (int col = 0; col < n; col++)
            {
                int start = columnptr2[col];
                int end = columnptr2[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices2[j];

                    tripletList.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(col),
                            mu * values2[j]));
                }
            }


            // epsilon I
            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Triplet64(
                        static_cast<int64_t>(i),
                        static_cast<int64_t>(i),
                        epsilon));
            }


            SparseMatrix64 mat(n, n);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();


            // ============================================================
            // B itself
            // ============================================================

            std::vector<Triplet64> tripletListB;

            tripletListB.reserve(values2->Length);

            for (int col = 0; col < n; col++)
            {
                int start = columnptr2[col];
                int end = columnptr2[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices2[j];
					double val = values2[j];
                    tripletListB.push_back(
                        Triplet64(
                            static_cast<int64_t>(row),
                            static_cast<int64_t>(col),
                            val));
                }
            }


            SparseMatrix64 matB(n, n);

            matB.setFromTriplets(
                tripletListB.begin(),
                tripletListB.end());

            matB.makeCompressed();


            // ============================================================
            // Factorization
            // ============================================================

            Eigen::PardisoLLT<SparseMatrix64> solver;

            solver.compute(mat);


            // Same kind of fallback as solve_CHOLESKY2
            double addedSalt = epsilon;

            for (int k = 0; k < 10; k++)
            {
                if (solver.info() == Eigen::Success)
                    break;

                addedSalt *= 10.0;

                for (int i = 0; i < n; i++)
                {
                    mat.coeffRef(i, i) += addedSalt;
                }

                mat.makeCompressed();

                solver.compute(mat);
            }


            if (solver.info() != Eigen::Success)
                return 0;


            // ============================================================
            // Initial solve
            //
            // x = M^{-1} mu rhs2
            // ============================================================

            Eigen::VectorXd b(n);

            for (int i = 0; i < n; i++)
            {
                b(i) = mu * rhs2[i];
            }


            Eigen::VectorXd x = solver.solve(b);

            if (solver.info() != Eigen::Success)
                return 0;


            // ============================================================
            // Null-space filtering iteration
            //
            // x <- M^{-1} mu B x
            // ============================================================

            for (int iter = 0; iter < numIterations; iter++)
            {
                b.noalias() = mu * (matB * x);

                x = solver.solve(b);

                if (solver.info() != Eigen::Success)
                    return 0;
            }


            // ============================================================
            // Return
            // ============================================================

            for (int i = 0; i < n; i++)
            {
                ret[i] = x(i);
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
