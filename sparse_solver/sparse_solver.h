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
        static array<double>^ findSearchDirection6(
            int n,

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            // Original basis:
            //
            //     P0 = [p0 p1 ...]
            //
            // q is appended internally:
            //
            //     P = [P0 q]
            //
            array<array<double>^>^ vecs,

            // Only the LAST vector is used:
            //
            //     q = vecs2[vecs2->Length - 1]
            //
            array<array<double>^>^ vecs2,

            array<double>^ rhs,
            array<double>^ rhs2,

            // M = a A + b B
            // R = c rhs + d rhs2
            double a,
            double b,
            double c,
            double d,

            // Relative TSVD threshold
            double epsilon,

            // Soft constraint weight:
            //
            // rho / 2 *
            //
            //     ( q^T direction / q^T q - 1 )^2
            //
            // qDotWeight = 0:
            //     ordinary reduced GN
            //
            double qDotWeight,

            // Reduced coefficients.
            //
            // Last entry = coefficient of q.
            //
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

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (vecs2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs2");
            }

            if (vecs2->Length == 0)
            {
                throw gcnew ArgumentException(
                    "vecs2 must contain at least one vector");
            }

            if (!std::isfinite(a) ||
                !std::isfinite(b) ||
                !std::isfinite(c) ||
                !std::isfinite(d))
            {
                throw gcnew ArgumentException(
                    "a, b, c, and d must be finite");
            }

            if (!std::isfinite(epsilon) ||
                epsilon < 0.0)
            {
                throw gcnew ArgumentException(
                    "epsilon must be finite and non-negative");
            }

            if (!std::isfinite(qDotWeight) ||
                qDotWeight < 0.0)
            {
                throw gcnew ArgumentException(
                    "qDotWeight must be finite and non-negative");
            }

            parameters = nullptr;


            // =====================================================
            // Available inputs
            // =====================================================

            const bool hasA =
                columnptr != nullptr &&
                rowindices != nullptr &&
                values != nullptr;

            const bool hasB =
                columnptr2 != nullptr &&
                rowindices2 != nullptr &&
                values2 != nullptr;

            const bool hasRhs =
                rhs != nullptr;

            const bool hasRhs2 =
                rhs2 != nullptr;


            const double effectiveA =
                hasA ? a : 0.0;

            const double effectiveB =
                hasB ? b : 0.0;

            const double effectiveC =
                hasRhs ? c : 0.0;

            const double effectiveD =
                hasRhs2 ? d : 0.0;


            // =====================================================
            // Basis dimensions
            //
            // P = [P0 q]
            // =====================================================

            const Eigen::Index originalSize =
                static_cast<Eigen::Index>(
                    vecs->Length);

            const Eigen::Index reducedSize =
                originalSize + 1;

            const Eigen::Index qIndex =
                originalSize;


            // =====================================================
            // Read q
            // =====================================================

            array<double>^ qSrc =
                vecs2[
                    vecs2->Length - 1];

            if (qSrc == nullptr)
            {
                throw gcnew ArgumentException(
                    "last vector of vecs2 is null");
            }

            if (qSrc->Length != n)
            {
                throw gcnew ArgumentException(
                    "last vector of vecs2 has incompatible size");
            }


            Eigen::VectorXd q(n);

            for (int i = 0; i < n; ++i)
            {
                const double value =
                    qSrc[i];

                if (!std::isfinite(value))
                {
                    throw gcnew ArgumentException(
                        "last vector of vecs2 contains NaN or infinity");
                }

                q(i) =
                    value;
            }


            const double qNorm2 =
                q.squaredNorm();

            if (!std::isfinite(qNorm2))
            {
                throw gcnew ArgumentException(
                    "q has invalid squared norm");
            }

            if (qNorm2 <= 0.0)
            {
                throw gcnew ArgumentException(
                    "q must be nonzero");
            }


            // =====================================================
            // Build augmented basis
            //
            //     P = [P0 q]
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);


            for (Eigen::Index j = 0;
                j < originalSize;
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


                for (int i = 0; i < n; ++i)
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


            // q is the LAST basis vector

            P.col(qIndex) =
                q;


            // =====================================================
            // Sparse matrices
            // =====================================================

            SparseMatrix64 matA(
                n,
                n);

            SparseMatrix64 matB(
                n,
                n);

            SparseMatrix64 matM(
                n,
                n);


            std::vector<Triplet64> tripletsA;
            std::vector<Triplet64> tripletsB;


            const int nnzA =
                hasA
                ? values->Length
                : 0;

            const int nnzB =
                hasB
                ? values2->Length
                : 0;


            tripletsA.reserve(
                static_cast<std::size_t>(
                    nnzA));

            tripletsB.reserve(
                static_cast<std::size_t>(
                    nnzB));


            // =====================================================
            // Build A
            // =====================================================

            if (effectiveA != 0.0)
            {
                if (columnptr->Length != n + 1)
                {
                    throw gcnew ArgumentException(
                        "columnptr must have length n + 1");
                }

                if (rowindices->Length !=
                    values->Length)
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
                            static_cast<int64_t>(
                                row),
                            static_cast<int64_t>(
                                column),
                            value);
                    }
                }


                matA.setFromTriplets(
                    tripletsA.begin(),
                    tripletsA.end());

                matA.makeCompressed();
            }
            else
            {
                matA.setZero();
            }


            // =====================================================
            // Build B
            // =====================================================

            if (effectiveB != 0.0)
            {
                if (columnptr2->Length != n + 1)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 must have length n + 1");
                }

                if (rowindices2->Length !=
                    values2->Length)
                {
                    throw gcnew ArgumentException(
                        "rowindices2 and values2 have incompatible sizes");
                }

                if (columnptr2[0] != 0 ||
                    columnptr2[n] != values2->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 has invalid start or end");
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


                        tripletsB.emplace_back(
                            static_cast<int64_t>(
                                row),
                            static_cast<int64_t>(
                                column),
                            value);
                    }
                }


                matB.setFromTriplets(
                    tripletsB.begin(),
                    tripletsB.end());

                matB.makeCompressed();
            }
            else
            {
                matB.setZero();
            }


            // =====================================================
            // M = a A + b B
            // =====================================================

            if (effectiveA != 0.0 &&
                effectiveB != 0.0)
            {
                matM =
                    effectiveA * matA +
                    effectiveB * matB;
            }
            else if (effectiveA != 0.0)
            {
                matM =
                    effectiveA * matA;
            }
            else if (effectiveB != 0.0)
            {
                matM =
                    effectiveB * matB;
            }
            else
            {
                matM.setZero();
            }

            matM.makeCompressed();


            // =====================================================
            // R = c rhs + d rhs2
            // =====================================================

            Eigen::VectorXd combinedRhs(n);

            combinedRhs.setZero();


            if (effectiveC != 0.0)
            {
                if (rhs->Length != n)
                {
                    throw gcnew ArgumentException(
                        "rhs has incompatible size");
                }


                for (int i = 0; i < n; ++i)
                {
                    const double value =
                        rhs[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "rhs contains NaN or infinity");
                    }

                    combinedRhs(i) +=
                        effectiveC *
                        value;
                }
            }


            if (effectiveD != 0.0)
            {
                if (rhs2->Length != n)
                {
                    throw gcnew ArgumentException(
                        "rhs2 has incompatible size");
                }


                for (int i = 0; i < n; ++i)
                {
                    const double value =
                        rhs2[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "rhs2 contains NaN or infinity");
                    }

                    combinedRhs(i) +=
                        effectiveD *
                        value;
                }
            }


            // =====================================================
            // MP = M P
            // =====================================================

            Eigen::MatrixXd MP(
                n,
                reducedSize);


            if (effectiveA == 0.0 &&
                effectiveB == 0.0)
            {
                MP.setZero();
            }
            else
            {
                MP.noalias() =
                    matM *
                    P;
            }


            if (!MP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "M * P contains NaN or infinity");
            }


            // =====================================================
            // Reduced GN
            //
            // direction = P x
            //
            // Base:
            //
            //     H0 = P^T M P
            //     g0 = P^T R
            //
            //
            // Soft constraint:
            //
            //     rho / 2 *
            //
            //       ( q^T direction / q^T q - 1 )^2
            //
            //
            // Define:
            //
            //     u = P^T q / (q^T q)
            //
            // Then:
            //
            //     H = H0 + rho u u^T
            //
            //     g = g0 + rho u
            //
            // =====================================================

            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);


            // Base primary / mixed GN

            H.noalias() =
                P.transpose() *
                MP;

            g.noalias() =
                P.transpose() *
                combinedRhs;


            // =====================================================
            // q projection
            //
            //     rawProjectedQ = P^T q
            //
            //     projectedQ =
            //         P^T q / (q^T q)
            // =====================================================

            Eigen::VectorXd rawProjectedQ(
                reducedSize);

            rawProjectedQ.noalias() =
                P.transpose() *
                q;


            Eigen::VectorXd projectedQ(
                reducedSize);

            projectedQ.noalias() =
                rawProjectedQ /
                qNorm2;


            // =====================================================
            // Soft:
            //
            //     q^T direction ~= q^T q
            //
            // normalized as
            //
            //     q^T direction / q^T q ~= 1
            // =====================================================

            if (qDotWeight != 0.0)
            {
                H.noalias() +=
                    qDotWeight *
                    projectedQ *
                    projectedQ.transpose();

                g.noalias() +=
                    qDotWeight *
                    projectedQ;
            }


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
            // TSVD
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


            Eigen::VectorXd solution(
                reducedSize);

            Eigen::VectorXd transformedRhs(
                reducedSize);

            solution.setZero();


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
                        epsilon *
                        sigmaMax;


                    transformedRhs.noalias() =
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


                    solution.noalias() =
                        svd.matrixV() *
                        transformedRhs;
                }
            }


            if (!solution.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Recover full direction
            //
            //     direction = P solution
            //
            // Since
            //
            //     P = [P0 q]
            //
            // direction =
            //
            //     P0 alpha + beta q
            //
            // =====================================================

            Eigen::VectorXd direction(n);

            direction.noalias() =
                P *
                solution;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return
            //
            // parameters[last]
            //     = coefficient beta of q
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(n);


            parameters =
                gcnew array<double>(
                    static_cast<int>(
                        reducedSize));


            for (int i = 0; i < n; ++i)
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
                    solution(j);
            }


            return ret;
        }
        static array<double>^ findSearchDirection5(
            int n,

            // =====================================================
            // Matrix A
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Matrix B
            // =====================================================
            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            // =====================================================
            // Basis
            //
            // Each basis vector:
            //
            //     [ vecs[i]  ]
            //     [ vecs2[i] ]
            //
            // size = 2N
            // =====================================================
            array<array<double>^>^ vecs,
            array<array<double>^>^ vecs2,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,
            array<double>^ rhs2,

            // =====================================================
            // CHOLESKY3 epsilon
            // =====================================================
            double epsilon)
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
            // TSVD threshold
            //
            // Separate from CHOLESKY3 epsilon.
            // =====================================================

            const double tsvdEpsilon =
                1.0E-12;


            const Eigen::Index fullSize =
                static_cast<Eigen::Index>(
                    2 * n);


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (vecs2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs2");
            }

            if (vecs->Length != vecs2->Length)
            {
                throw gcnew ArgumentException(
                    "vecs and vecs2 must contain the same number of vectors");
            }

            if (rhs == nullptr ||
                rhs2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs and rhs2 must not be null");
            }

            if (rhs->Length != n ||
                rhs2->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs and rhs2 must have length n");
            }

            if (!std::isfinite(epsilon) ||
                epsilon < 0.0)
            {
                throw gcnew ArgumentException(
                    "epsilon must be finite and non-negative");
            }


            // =====================================================
            // Number of basis vectors
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    vecs->Length);

            if (reducedSize == 0)
            {
                return gcnew array<double>(
                    2 * n);
            }


            // =====================================================
            // Build 2N-dimensional basis
            //
            // P =
            //
            // [ vecs[0]   vecs[1]   ... ]
            // [ vecs2[0]  vecs2[1]  ... ]
            //
            // size = 2N x reducedSize
            //
            // No QR / orthonormalization.
            // =====================================================

            Eigen::MatrixXd P(
                fullSize,
                reducedSize);

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                array<double>^ upper =
                    vecs[
                        static_cast<int>(j)];

                array<double>^ lower =
                    vecs2[
                        static_cast<int>(j)];


                if (upper == nullptr ||
                    lower == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs or vecs2 contains a null vector");
                }

                if (upper->Length != n ||
                    lower->Length != n)
                {
                    throw gcnew ArgumentException(
                        "All basis vectors must have length n");
                }


                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double v1 =
                        upper[i];

                    const double v2 =
                        lower[i];


                    if (!std::isfinite(v1) ||
                        !std::isfinite(v2))
                    {
                        throw gcnew ArgumentException(
                            "Basis contains NaN or infinity");
                    }


                    P(i, j) =
                        v1;

                    P(n + i, j) =
                        v2;
                }
            }


            // =====================================================
            // Full 2N x 2N matrix
            //
            // M =
            //
            // [ I       A          ]
            // [ A      -epsilon B  ]
            //
            // =====================================================

            SparseMatrix64 M(
                fullSize,
                fullSize);

            std::vector<Triplet64> triplets;

            triplets.reserve(
                static_cast<std::size_t>(
                    n +
                    2 * values->Length +
                    values2->Length));


            // =====================================================
            // Upper-left: I
            // =====================================================

            for (int i = 0;
                i < n;
                ++i)
            {
                triplets.emplace_back(
                    static_cast<int64_t>(i),
                    static_cast<int64_t>(i),
                    1.0);
            }


            // =====================================================
            // A blocks
            //
            // Upper-right: A
            // Lower-left : A
            // =====================================================

            if (columnptr == nullptr ||
                rowindices == nullptr ||
                values == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "Matrix A is null");
            }

            if (columnptr->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");
            }

            if (rowindices->Length !=
                values->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");
            }


            for (int col = 0;
                col < n;
                ++col)
            {
                const int start =
                    columnptr[col];

                const int end =
                    columnptr[col + 1];


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
                            "Matrix A contains invalid row index");
                    }

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "Matrix A contains NaN or infinity");
                    }


                    // ---------------------------------------------
                    // Upper-right A
                    // ---------------------------------------------

                    triplets.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(n + col),
                        value);


                    // ---------------------------------------------
                    // Lower-left A
                    // ---------------------------------------------

                    triplets.emplace_back(
                        static_cast<int64_t>(n + row),
                        static_cast<int64_t>(col),
                        value);
                }
            }


            // =====================================================
            // Lower-right: -epsilon B
            // =====================================================

            if (columnptr2 == nullptr ||
                rowindices2 == nullptr ||
                values2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "Matrix B is null");
            }

            if (columnptr2->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr2 must have length n + 1");
            }

            if (rowindices2->Length !=
                values2->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices2 and values2 have incompatible sizes");
            }


            for (int col = 0;
                col < n;
                ++col)
            {
                const int start =
                    columnptr2[col];

                const int end =
                    columnptr2[col + 1];


                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices2[index];

                    const double value =
                        values2[index];


                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "Matrix B contains invalid row index");
                    }

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "Matrix B contains NaN or infinity");
                    }


                    triplets.emplace_back(
                        static_cast<int64_t>(n + row),
                        static_cast<int64_t>(n + col),
                        -epsilon * value);
                }
            }


            // =====================================================
            // Build full matrix
            // =====================================================

            M.setFromTriplets(
                triplets.begin(),
                triplets.end());

            M.makeCompressed();


            // =====================================================
            // Full RHS
            //
            // R =
            //
            // [ rhs           ]
            // [ -epsilon rhs2 ]
            //
            // =====================================================

            Eigen::VectorXd R(
                fullSize);


            for (int i = 0;
                i < n;
                ++i)
            {
                const double r1 =
                    rhs[i];

                const double r2 =
                    rhs2[i];


                if (!std::isfinite(r1) ||
                    !std::isfinite(r2))
                {
                    throw gcnew ArgumentException(
                        "RHS contains NaN or infinity");
                }


                R(i) =
                    r1;

                R(n + i) =
                    -epsilon * r2;
            }


            // =====================================================
            // Reduced problem
            //
            // MP = M P
            //
            // H = P^T M P
            //
            // g = P^T R
            //
            // =====================================================

            Eigen::MatrixXd MP(
                fullSize,
                reducedSize);

            MP.noalias() =
                M * P;


            if (!MP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "M * P contains NaN or infinity");
            }


            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            H.noalias() =
                P.transpose() *
                MP;


            Eigen::VectorXd g(
                reducedSize);

            g.noalias() =
                P.transpose() *
                R;


            if (!H.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix contains NaN or infinity");
            }

            if (!g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced RHS contains NaN or infinity");
            }


            // =====================================================
            // TSVD
            //
            // H c = g
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


            Eigen::VectorXd coefficients(
                reducedSize);

            coefficients.setZero();


            if (singularValues.size() > 0)
            {
                const double sigmaMax =
                    singularValues(0);


                if (sigmaMax > 0.0)
                {
                    const double threshold =
                        tsvdEpsilon *
                        sigmaMax;


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


                    coefficients.noalias() =
                        svd.matrixV() *
                        transformedRhs;
                }
            }


            if (!coefficients.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Full 2N-dimensional result
            //
            // direction = P c
            //
            // Do NOT split lambda/x.
            // =====================================================

            Eigen::VectorXd direction(
                fullSize);

            direction.noalias() =
                P * coefficients;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return one 2N-dimensional vector
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(
                    2 * n);


            for (int i = 0;
                i < 2 * n;
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
            // =====================================================
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Matrix B
            // =====================================================
            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            // =====================================================
            // Basis
            //
            // vecs[i]  : upper N components
            // vecs2[i] : lower N components
            //
            // Each basis vector is
            //
            //       [ vecs[i]  ]
            //       [ vecs2[i] ]
            //
            // =====================================================
            array<array<double>^>^ vecs,
            array<array<double>^>^ vecs2,

            // =====================================================
            // RHS
            // =====================================================
            array<double>^ rhs,
            array<double>^ rhs2,
            array<double>^ rhs3,

            // =====================================================
            // Parameters
            // =====================================================
            double mu,
            double lambda)
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
            // Settings
            // =====================================================

            const double epsilon = 1.0E-12;

            const Eigen::Index fullSize =
                static_cast<Eigen::Index>(2 * n);


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (vecs2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs2");
            }

            if (vecs->Length != vecs2->Length)
            {
                throw gcnew ArgumentException(
                    "vecs and vecs2 must contain the same number of vectors");
            }

            if (rhs == nullptr ||
                rhs2 == nullptr ||
                rhs3 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "rhs, rhs2, and rhs3 must not be null");
            }

            if (rhs->Length != n ||
                rhs2->Length != n ||
                rhs3->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs, rhs2, and rhs3 must have length n");
            }

            if (!std::isfinite(mu) ||
                !std::isfinite(lambda))
            {
                throw gcnew ArgumentException(
                    "mu and lambda must be finite");
            }

            if (mu < 0.0)
            {
                throw gcnew ArgumentException(
                    "mu must be non-negative");
            }

            if (lambda < 0.0)
            {
                throw gcnew ArgumentException(
                    "lambda must be non-negative");
            }


            // =====================================================
            // Number of basis vectors
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    vecs->Length);

            if (reducedSize == 0)
            {
                return gcnew array<double>(2 * n);
            }


            // =====================================================
            // Dense basis matrix
            //
            // P =
            //
            // [ vecs[0]   vecs[1]   ... ]
            // [ vecs2[0]  vecs2[1]  ... ]
            //
            // size = 2N x reducedSize
            //
            // No QR / orthonormalization.
            // =====================================================

            Eigen::MatrixXd P(
                fullSize,
                reducedSize);

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                array<double>^ upper =
                    vecs[static_cast<int>(j)];

                array<double>^ lower =
                    vecs2[static_cast<int>(j)];

                if (upper == nullptr ||
                    lower == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs or vecs2 contains a null vector");
                }

                if (upper->Length != n ||
                    lower->Length != n)
                {
                    throw gcnew ArgumentException(
                        "All basis vectors must have length n");
                }

                for (int i = 0; i < n; ++i)
                {
                    const double v1 =
                        upper[i];

                    const double v2 =
                        lower[i];

                    if (!std::isfinite(v1) ||
                        !std::isfinite(v2))
                    {
                        throw gcnew ArgumentException(
                            "Basis contains NaN or infinity");
                    }

                    P(i, j) =
                        v1;

                    P(n + i, j) =
                        v2;
                }
            }


            // =====================================================
            // Build full 2N x 2N sparse matrix
            //
            // M =
            //
            // [ A + lambda I       -lambda I      ]
            // [ -lambda I       mu B + lambda I   ]
            //
            // =====================================================

            SparseMatrix64 M(
                fullSize,
                fullSize);

            std::vector<Triplet64> triplets;

            triplets.reserve(
                static_cast<std::size_t>(
                    values->Length +
                    values2->Length +
                    4 * n));


            // -----------------------------------------------------
            // Upper-left block: A
            // -----------------------------------------------------

            if (columnptr == nullptr ||
                rowindices == nullptr ||
                values == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "Matrix A is null");
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

            for (int col = 0;
                col < n;
                ++col)
            {
                const int start =
                    columnptr[col];

                const int end =
                    columnptr[col + 1];

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
                            "Matrix A contains invalid row index");
                    }

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "Matrix A contains NaN or infinity");
                    }

                    triplets.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(col),
                        value);
                }
            }


            // -----------------------------------------------------
            // Lower-right block: mu B
            // -----------------------------------------------------

            if (columnptr2 == nullptr ||
                rowindices2 == nullptr ||
                values2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "Matrix B is null");
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

            for (int col = 0;
                col < n;
                ++col)
            {
                const int start =
                    columnptr2[col];

                const int end =
                    columnptr2[col + 1];

                for (int index = start;
                    index < end;
                    ++index)
                {
                    const int row =
                        rowindices2[index];

                    const double value =
                        values2[index];

                    if (row < 0 ||
                        row >= n)
                    {
                        throw gcnew ArgumentException(
                            "Matrix B contains invalid row index");
                    }

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "Matrix B contains NaN or infinity");
                    }

                    triplets.emplace_back(
                        static_cast<int64_t>(n + row),
                        static_cast<int64_t>(n + col),
                        mu * value);
                }
            }


            // -----------------------------------------------------
            // lambda ||x-y||^2
            //
            // [ +lambda I   -lambda I ]
            // [ -lambda I   +lambda I ]
            // -----------------------------------------------------

            for (int i = 0;
                i < n;
                ++i)
            {
                // upper-left
                triplets.emplace_back(
                    static_cast<int64_t>(i),
                    static_cast<int64_t>(i),
                    lambda);

                // lower-right
                triplets.emplace_back(
                    static_cast<int64_t>(n + i),
                    static_cast<int64_t>(n + i),
                    lambda);

                // upper-right
                triplets.emplace_back(
                    static_cast<int64_t>(i),
                    static_cast<int64_t>(n + i),
                    -lambda);

                // lower-left
                triplets.emplace_back(
                    static_cast<int64_t>(n + i),
                    static_cast<int64_t>(i),
                    -lambda);
            }


            M.setFromTriplets(
                triplets.begin(),
                triplets.end());

            M.makeCompressed();


            // =====================================================
            // Full RHS
            //
            // R =
            //
            // [ rhs       + lambda rhs3 ]
            // [ mu rhs2   - lambda rhs3 ]
            //
            // Here rhs3 = current x-y.
            //
            // This sign convention assumes that the returned
            // direction is later SUBTRACTED from the variables:
            //
            //     x_new = x_old - alpha * direction
            //
            // =====================================================

            Eigen::VectorXd R(fullSize);

            for (int i = 0;
                i < n;
                ++i)
            {
                const double r1 =
                    rhs[i];

                const double r2 =
                    rhs2[i];

                const double r3 =
                    rhs3[i];

                if (!std::isfinite(r1) ||
                    !std::isfinite(r2) ||
                    !std::isfinite(r3))
                {
                    throw gcnew ArgumentException(
                        "RHS contains NaN or infinity");
                }

                R(i) =
                    r1 +
                    lambda * r3;

                R(n + i) =
                    mu * r2 -
                    lambda * r3;
            }


            // =====================================================
            // Reduced matrix
            //
            // MP = M P
            //
            // H = P^T M P
            //
            // g = P^T R
            // =====================================================

            Eigen::MatrixXd MP(
                fullSize,
                reducedSize);

            MP.noalias() =
                M * P;

            if (!MP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "M * P contains NaN or infinity");
            }


            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            H.noalias() =
                P.transpose() *
                MP;


            Eigen::VectorXd g(
                reducedSize);

            g.noalias() =
                P.transpose() *
                R;


            if (!H.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced matrix contains NaN or infinity");
            }

            if (!g.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced RHS contains NaN or infinity");
            }


            // =====================================================
            // TSVD
            //
            // H c = g
            // =====================================================

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(
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


            Eigen::VectorXd coefficients(
                reducedSize);

            coefficients.setZero();


            if (singularValues.size() > 0)
            {
                const double sigmaMax =
                    singularValues(0);

                if (sigmaMax > 0.0)
                {
                    const double threshold =
                        epsilon *
                        sigmaMax;

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

                    coefficients.noalias() =
                        svd.matrixV() *
                        transformedRhs;
                }
            }


            if (!coefficients.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Full 2N-dimensional direction
            //
            // direction = P coefficients
            // =====================================================

            Eigen::VectorXd direction(
                fullSize);

            direction.noalias() =
                P * coefficients;


            if (!direction.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }


            // =====================================================
            // Return one 2N vector.
            //
            // No splitting into x/y vectors.
            // =====================================================

            array<double>^ ret =
                gcnew array<double>(2 * n);

            for (int i = 0;
                i < 2 * n;
                ++i)
            {
                ret[i] =
                    direction(i);
            }

            return ret;
        }
        static array<array<double>^>^ findSearchDirection3(
            int n,
            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            int n2,
            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            array<array<double>^>^ vecs,
            array<array<double>^>^ nearNullVecs,

            array<double>^ rhs,
            array<double>^ rhs2,
            double mu,
            double tau,
            double lambda,
            double epsilon,
            bool removeNearNullFromPrimary,
        int method)
        {
            array<double>  ^ projectedBasisNorms =
                gcnew array<double>(0);

            if (n <= 0 ||
                n2 <= 0)
            {
                throw gcnew ArgumentException(
                    "n and n2 must be positive");
            }

            if (n != n2)
            {
                throw gcnew ArgumentException(
                    "mat and mat2 must have the same size");
            }

            if (columnptr == nullptr ||
                rowindices == nullptr ||
                values == nullptr ||
                columnptr2 == nullptr ||
                rowindices2 == nullptr ||
                values2 == nullptr ||
                vecs == nullptr ||
                nearNullVecs == nullptr ||
                rhs == nullptr ||
                rhs2 == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "Input array is null");
            }

            if (columnptr->Length != n + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr must have length n + 1");
            }

            if (columnptr2->Length != n2 + 1)
            {
                throw gcnew ArgumentException(
                    "columnptr2 must have length n2 + 1");
            }

            if (rowindices->Length !=
                values->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices and values have incompatible sizes");
            }

            if (rowindices2->Length !=
                values2->Length)
            {
                throw gcnew ArgumentException(
                    "rowindices2 and values2 have incompatible sizes");
            }

            if (rhs->Length != n)
            {
                throw gcnew ArgumentException(
                    "rhs has incompatible size");
            }

            if (rhs2->Length != n2)
            {
                throw gcnew ArgumentException(
                    "rhs2 has incompatible size");
            }

            if (!std::isfinite(epsilon) ||
                epsilon < 0.0)
            {
                throw gcnew ArgumentException(
                    "epsilon must be finite and non-negative");
            }

            if (!std::isfinite(mu) ||
                mu < 0.0)
            {
                throw gcnew ArgumentException(
                    "mu must be finite and non-negative");
            }

            if (columnptr[0] != 0 ||
                columnptr[n] != values->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr has invalid start or end");
            }

            if (columnptr2[0] != 0 ||
                columnptr2[n2] != values2->Length)
            {
                throw gcnew ArgumentException(
                    "columnptr2 has invalid start or end");
            }

            using SparseMatrix64 =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            // =================================================
            // Primary basis
            // =================================================

            std::vector<Eigen::VectorXd>
                nativeVecs(vecs->Length);

            for (int i = 0;
                i < vecs->Length;
                ++i)
            {
                if (vecs[i] == nullptr)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a null vector");
                }

                if (vecs[i]->Length != n)
                {
                    throw gcnew ArgumentException(
                        "vecs contains a vector with incompatible size");
                }

                Eigen::VectorXd v(n);

                for (int j = 0;
                    j < n;
                    ++j)
                {
                    const double value =
                        vecs[i][j];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "vecs contains NaN or infinity");
                    }

                    v(j) =
                        value;
                }

                nativeVecs[i] =
                    std::move(v);
            }

            // =================================================
            // Near-null basis
            // =================================================

            std::vector<Eigen::VectorXd>
                nativeNearNullVecs(
                    nearNullVecs->Length);

            for (int i = 0;
                i < nearNullVecs->Length;
                ++i)
            {
                if (nearNullVecs[i] == nullptr)
                {
                    throw gcnew ArgumentException(
                        "nearNullVecs contains a null vector");
                }

                if (nearNullVecs[i]->Length != n)
                {
                    throw gcnew ArgumentException(
                        "nearNullVecs contains a vector with incompatible size");
                }

                Eigen::VectorXd v(n);

                for (int j = 0;
                    j < n;
                    ++j)
                {
                    const double value =
                        nearNullVecs[i][j];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "nearNullVecs contains NaN or infinity");
                    }

                    v(j) =
                        value;
                }

                nativeNearNullVecs[i] =
                    std::move(v);
            }

            // =================================================
            // First sparse matrix
            // =================================================

            SparseMatrix64 mat(n, n);

            std::vector<
                Eigen::Triplet<double, int64_t>>
                triplets;

            triplets.reserve(
                values->Length);

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

                    triplets.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }

            mat.setFromTriplets(
                triplets.begin(),
                triplets.end());

            mat.makeCompressed();

            // =================================================
            // Second sparse matrix
            // =================================================

            SparseMatrix64 mat2(n2, n2);

            std::vector<
                Eigen::Triplet<double, int64_t>>
                triplets2;

            triplets2.reserve(
                values2->Length );

          

            for (int column = 0;
                column < n2;
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
                        row >= n2)
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

                    triplets2.emplace_back(
                        static_cast<int64_t>(row),
                        static_cast<int64_t>(column),
                        value);
                }
            }

            mat2.setFromTriplets(
                triplets2.begin(),
                triplets2.end());
                
       

            mat2.makeCompressed();

            // =================================================
            // First RHS
            // =================================================

            Eigen::VectorXd nativeRhs(n);

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

                nativeRhs(i) =
                    value;
            }

            // =================================================
            // Second RHS
            // =================================================

            Eigen::VectorXd nativeRhs2(n2);

            for (int i = 0;
                i < n2;
                ++i)
            {
                const double value =
                    rhs2[i];

                if (!std::isfinite(value))
                {
                    throw gcnew ArgumentException(
                        "rhs2 contains NaN or infinity");
                }

                nativeRhs2(i) =
                    value;
            }

            // =================================================
            // Rebuild the system used by the secondary solve:
            //
            //     A_secondary <- A_first + mu * A_secondary
            //     b_secondary <- b_first + mu * b_secondary
            //
            // The first system itself remains unchanged.
            // =================================================

           
            mat2 *= mu;
            mat2 += mat;
            
            mat2.makeCompressed();

           

            nativeRhs2 =
            // nativeRhs + 
             mu * nativeRhs2;

            // =================================================
            // Call the native Eigen function
            // =================================================

            std::vector<double>
                nativeProjectedBasisNorms;
            std::array<Eigen::VectorXd, 4> result;
            if (method == 0) {
                result =
                    _sparse_solver::findSearchDirection3(
                        mat,
                        mat2,
                        nativeVecs,
                        nativeNearNullVecs,
                        nativeRhs,
                        nativeRhs2,
                        nativeProjectedBasisNorms,
                        mu,
                        tau,
                        epsilon,
                        false
                    );
            }
            else if(method == 1){
                 result =
                    _sparse_solver::findSearchDirection4(
                        mat,
                        mat2,
                        nativeVecs,
                        nativeNearNullVecs,
                        nativeRhs,
                        nativeRhs2,
                        nativeProjectedBasisNorms,
                        mu,
                        tau,
                        epsilon,
                        false
                    );
            }
            else if (method==2) {
                result =
                    _sparse_solver::findSearchDirection5(
                        mat,
                        mat2,
                        nativeVecs,
                        nativeNearNullVecs,
                        nativeRhs,
                        nativeRhs2,
                        nativeProjectedBasisNorms,
                        mu,
                        tau,
                        epsilon,
                        false
                    );
            }
            else if (method == 3) {
                result =
                    _sparse_solver::findSearchDirection6(
                        mat,
                        mat2,
                        nativeVecs,
                        nativeNearNullVecs,
                        nativeRhs,
                        nativeRhs2,
                        nativeProjectedBasisNorms,
                        mu,
                        tau,
                        epsilon,
                        false
                    );
			}
			else if (method==4){
				result =
                    _sparse_solver::findSearchDirection7(
                        mat,
                        mat2,
                        nativeVecs,
                        nativeNearNullVecs,
                        nativeRhs,
                        nativeRhs2,
                        nativeProjectedBasisNorms,
                        mu,
                        tau,
                        epsilon,
                        false
                    );
			}
else if (method == 5) {
    result =
        _sparse_solver::findSearchDirection8(
            mat,
            mat2,
            nativeVecs,
            nativeNearNullVecs,
            nativeRhs,
            nativeRhs2,
            nativeProjectedBasisNorms,
            mu,
            tau,
            epsilon,
            false
        );
        }
else if(method==6){
            result =
                _sparse_solver::findSearchDirection9(
                    mat,
                    mat2,
                    nativeVecs,
                    nativeNearNullVecs,
                    nativeRhs,
                    nativeRhs2,
                    nativeProjectedBasisNorms,
                    lambda,
                    tau,
                    epsilon,
                    false
                );
                }
else if(method==7) {
  
        result =
            _sparse_solver::findSearchDirection10(
                mat,
                mat2,
                nativeVecs,
                nativeNearNullVecs,
                nativeRhs,
                nativeRhs2,
                nativeProjectedBasisNorms,
                lambda,
                tau,
                epsilon,
                false
            );
 
                }
else if(method==8) {
                    result =
                        _sparse_solver::findSearchDirection11(
                            mat,
                            mat2,
                            nativeVecs,
                            nativeNearNullVecs,
                            nativeRhs,
                            nativeRhs2,
                            nativeProjectedBasisNorms,
                            lambda,
                            tau,
                            epsilon,
                            false
                        );

                        }
else if(method==9) {
                            result =
                                _sparse_solver::findSearchDirection12(
                                    mat,
                                    mat2,
                                    nativeVecs,
                                    nativeNearNullVecs,
                                    nativeRhs,
                                    nativeRhs2,
                                    nativeProjectedBasisNorms,
                                    lambda,
                                    tau,
                                    epsilon,
                                    false
                                );

                                }
else {
    result =
        _sparse_solver::findSearchDirection13(
            mat,
            mat2,
            nativeVecs,
            nativeNearNullVecs,
            nativeRhs,
            nativeRhs2,
            nativeProjectedBasisNorms,
            lambda,
            tau,
            epsilon,
            false
        );

        }
                        
                

            if (result.size() != 4)
            {
                throw gcnew InvalidOperationException(
                    "findSearchDirection3 must return two vectors");
            }

            


          
  
          

            // =================================================
            // Convert projected gradient norm
            // =================================================

            projectedBasisNorms =
                gcnew array<double>(3);

            projectedBasisNorms[0] =
                nativeProjectedBasisNorms[0];
            projectedBasisNorms[1] =
                nativeProjectedBasisNorms[1];
            projectedBasisNorms[2] =
                nativeProjectedBasisNorms[2];
            // =================================================
            // Convert search direction
            // =================================================

            array<array<double>^>^ ret =
                gcnew array<array<double>^>(4);

            for (int i = 0;
                i < 4;
                ++i)
            {
                ret[i] =
					gcnew array<double>(result[i].size());
				for (int j = 0;
					j < result[i].size();
					++j)
				{
					ret[i][j] =
						result[i](j);
				}
            }

            return ret;
        }
        

        static array<double>^ findSearchDirection(
            int n,

            // =====================================================
            // Matrix A
            //
            // CSC format
            //
            // If unavailable, coefficient a is forced to zero.
            // =====================================================

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Matrix B
            //
            // CSC format
            //
            // If unavailable, coefficient b is forced to zero.
            // =====================================================

            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            // =====================================================
            // Basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // Used exactly as supplied.
            //
            // No QR or orthonormalization is performed.
            // Linear dependence is handled by the TSVD of
            // the reduced matrix P^T M P.
            // =====================================================

            array<array<double>^>^ vecs,

            // =====================================================
            // RHS vectors
            // =====================================================

            array<double>^ rhs,
            array<double>^ rhs2,

            // =====================================================
            // Mixing coefficients
            //
            //     M = a A + b B
            //
            //     R = c rhs + d rhs2
            // =====================================================

            double a,
            double b,
            double c,
            double d,

            // =====================================================
            // Relative TSVD threshold
            // =====================================================

            double epsilon,

            // =====================================================
            // Previous reduced coefficients
            //
            // If nullptr, temporal regularization is disabled.
            // No normalization or orthogonalization is applied.
            // =====================================================

            array<double>^ previousParameters,

            // =====================================================
            // Temporal regularization weight
            //
            // Adds
            //
            //   previousParameterWeight / 2
            //       * ||x - previousParameters||^2
            //
            // to the reduced quadratic objective.
            // =====================================================

            double previousParameterWeight,

            // =====================================================
            // Output reduced coefficients used for this direction
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
            // Validate common parameters
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (!std::isfinite(a) ||
                !std::isfinite(b) ||
                !std::isfinite(c) ||
                !std::isfinite(d))
            {
                throw gcnew ArgumentException(
                    "a, b, c, and d must be finite");
            }

            if (!std::isfinite(epsilon) ||
                epsilon < 0.0)
            {
                throw gcnew ArgumentException(
                    "epsilon must be finite and non-negative");
            }

            if (!std::isfinite(previousParameterWeight) ||
                previousParameterWeight < 0.0)
            {
                throw gcnew ArgumentException(
                    "previousParameterWeight must be finite and non-negative");
            }

            parameters = nullptr;


            // =====================================================
            // Determine available inputs
            // =====================================================

            const bool hasA =
                columnptr != nullptr &&
                rowindices != nullptr &&
                values != nullptr;

            const bool hasB =
                columnptr2 != nullptr &&
                rowindices2 != nullptr &&
                values2 != nullptr;

            const bool hasRhs =
                rhs != nullptr;

            const bool hasRhs2 =
                rhs2 != nullptr;


            const double effectiveA =
                hasA ? a : 0.0;

            const double effectiveB =
                hasB ? b : 0.0;

            const double effectiveC =
                hasRhs ? c : 0.0;

            const double effectiveD =
                hasRhs2 ? d : 0.0;


            // =====================================================
            // Convert managed basis vectors to Eigen vectors
            // =====================================================

            std::vector<Eigen::VectorXd> inputBasis;

            inputBasis.reserve(
                static_cast<std::size_t>(
                    vecs->Length));

            for (int j = 0;
                j < vecs->Length;
                ++j)
            {
                array<double>^ src =
                    vecs[j];

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

                inputBasis.emplace_back(n);

                Eigen::VectorXd& vector =
                    inputBasis.back();

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

                    vector(i) =
                        value;
                }
            }


            // =====================================================
            // Keep the original basis exactly as supplied
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    inputBasis.size());

            if (reducedSize == 0)
            {
                parameters = gcnew array<double>(0);
                return gcnew array<double>(n);
            }


            // =====================================================
            // Previous coefficient vector
            //
            // The coefficient coordinates are meaningful across
            // iterations because P is used exactly as supplied:
            // NO QR, NO normalization, NO orthogonalization.
            // =====================================================

            Eigen::VectorXd previous =
                Eigen::VectorXd::Zero(reducedSize);

            const bool hasPrevious =
                previousParameters != nullptr;

            if (hasPrevious)
            {
                if (previousParameters->Length !=
                    static_cast<int>(reducedSize))
                {
                    throw gcnew ArgumentException(
                        "previousParameters has incompatible size");
                }

                for (Eigen::Index j = 0;
                    j < reducedSize;
                    ++j)
                {
                    const double value =
                        previousParameters[
                            static_cast<int>(j)];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "previousParameters contains NaN or infinity");
                    }

                    previous(j) =
                        value;
                }
            }


            // =====================================================
            // Number of nonzeros
            // =====================================================

            const int nnzA =
                hasA
                ? values->Length
                : 0;

            const int nnzB =
                hasB
                ? values2->Length
                : 0;


            // =====================================================
            // Local working storage
            //
            // Every invocation owns all mutable buffers. Therefore,
            // concurrent calls do not share Eigen matrices, vectors,
            // sparse matrices, or triplet arrays.
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);

            Eigen::MatrixXd MP(
                n,
                reducedSize);

            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);

            Eigen::VectorXd combinedRhs(
                n);

            Eigen::VectorXd solution(
                reducedSize);

            Eigen::VectorXd transformedRhs(
                reducedSize);

            SparseMatrix64 matA(
                n,
                n);

            SparseMatrix64 matB(
                n,
                n);

            SparseMatrix64 matM(
                n,
                n);

            std::vector<Triplet64> tripletsA;
            std::vector<Triplet64> tripletsB;

            tripletsA.reserve(
                static_cast<std::size_t>(
                    nnzA));

            tripletsB.reserve(
                static_cast<std::size_t>(
                    nnzB));


            // =====================================================
            // Build dense basis matrix without QR
            //
            //     P = [p0 p1 ...]
            //
            // The original order, scale, zero columns, and linear
            // dependence are preserved. Rank deficiency is handled
            // later by the TSVD of H = P^T M P.
            // =====================================================

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                P.col(j) =
                    inputBasis[
                        static_cast<std::size_t>(j)];
            }


            // =====================================================
            // Build matrix A
            //
            // Triplet storage is local to this invocation.
            // =====================================================

            if (effectiveA != 0.0)
            {
                if (columnptr->Length != n + 1)
                {
                    throw gcnew ArgumentException(
                        "columnptr must have length n + 1");
                }

                if (rowindices->Length !=
                    values->Length)
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

                tripletsA.clear();

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
                            static_cast<int64_t>(
                                row),

                            static_cast<int64_t>(
                                column),

                            value);
                    }
                }

                // -------------------------------------------------
                // Rebuild numerical sparse matrix.
                //
                // Triplet buffer itself is reused.
                // -------------------------------------------------

                matA.setZero();

                matA.setFromTriplets(
                    tripletsA.begin(),
                    tripletsA.end());

                matA.makeCompressed();
            }
            else
            {
                matA.setZero();
            }


            // =====================================================
            // Build matrix B
            // =====================================================

            if (effectiveB != 0.0)
            {
                if (columnptr2->Length != n + 1)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 must have length n + 1");
                }

                if (rowindices2->Length !=
                    values2->Length)
                {
                    throw gcnew ArgumentException(
                        "rowindices2 and values2 have incompatible sizes");
                }

                if (columnptr2[0] != 0 ||
                    columnptr2[n] != values2->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 has invalid start or end");
                }

                tripletsB.clear();

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

                        tripletsB.emplace_back(
                            static_cast<int64_t>(
                                row),

                            static_cast<int64_t>(
                                column),

                            value);
                    }
                }

                matB.setZero();

                matB.setFromTriplets(
                    tripletsB.begin(),
                    tripletsB.end());

                matB.makeCompressed();
            }
            else
            {
                matB.setZero();
            }


            // =====================================================
            // Build mixed matrix
            //
            //     M = a A + b B
            // =====================================================

            if (effectiveA != 0.0 &&
                effectiveB != 0.0)
            {
                matM =
                    effectiveA * matA +
                    effectiveB * matB;
            }
            else if (effectiveA != 0.0)
            {
                matM =
                    effectiveA * matA;
            }
            else if (effectiveB != 0.0)
            {
                matM =
                    effectiveB * matB;
            }
            else
            {
                matM.setZero();
            }

            matM.makeCompressed();


            // =====================================================
            // Build mixed RHS
            //
            //     R = c rhs + d rhs2
            //
            // combinedRhs is local to this invocation.
            // =====================================================

            combinedRhs.setZero();

            if (effectiveC != 0.0)
            {
                if (rhs->Length != n)
                {
                    throw gcnew ArgumentException(
                        "rhs has incompatible size");
                }

                double* dst =
                    combinedRhs.data();

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

                    dst[i] +=
                        effectiveC *
                        value;
                }
            }

            if (effectiveD != 0.0)
            {
                if (rhs2->Length != n)
                {
                    throw gcnew ArgumentException(
                        "rhs2 has incompatible size");
                }

                double* dst =
                    combinedRhs.data();

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double value =
                        rhs2[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "rhs2 contains NaN or infinity");
                    }

                    dst[i] +=
                        effectiveD *
                        value;
                }
            }


            // =====================================================
            // FAST sparse-dense multiplication
            //
            //     MP = M P
            //
            // MP is local to this invocation.
            // =====================================================

            if (effectiveA == 0.0 &&
                effectiveB == 0.0)
            {
                MP.setZero();
            }
            else
            {
                MP.noalias() =
                    matM * P;
            }


            if (!MP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "M * P contains NaN or infinity");
            }


            // =====================================================
            // FAST reduced matrix
            //
            //     H = P^T M P
            //
            // This replaces reducedSize^2 explicit VectorXd::dot()
            // operations with one dense matrix multiplication.
            // =====================================================

            H.noalias() =
                P.transpose() *
                MP;


            // =====================================================
            // Reduced RHS
            //
            //     g = P^T R
            // =====================================================

            g.noalias() =
                P.transpose() *
                combinedRhs;


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
            // Weak temporal regularization in coefficient space
            //
            // Original reduced objective:
            //
            //     1/2 x^T H x - g^T x
            //
            // Add:
            //
            //     w/2 ||x - x_prev||^2
            //
            // Therefore:
            //
            //     (H + w I) x = g + w x_prev
            //
            // If previousParameters == nullptr, no temporal
            // regularization is applied even when w > 0.
            // =====================================================

            if (hasPrevious &&
                previousParameterWeight > 0.0)
            {
                H.diagonal().array() +=
                    previousParameterWeight;

                g.noalias() +=
                    previousParameterWeight *
                    previous;
            }


            // =====================================================
            // TSVD
            //
            //     H = U S V^T
            //
            // Solve
            //
            //     H x = g
            //
            // with
            //
            //     sigma_i > epsilon * sigma_max
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


            solution.setZero();


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
                        epsilon *
                        sigmaMax;


                    transformedRhs.noalias() =
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


                    solution.noalias() =
                        svd.matrixV() *
                        transformedRhs;
                }
            }


            if (!solution.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Return the full-space direction
            //
            //     direction = P solution
            // =====================================================

            combinedRhs.noalias() =
                P * solution;


            if (!combinedRhs.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }

            array<double>^ ret =
                gcnew array<double>(n);

            parameters =
                gcnew array<double>(
                    static_cast<int>(reducedSize));


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    combinedRhs(i);
            }

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                parameters[
                    static_cast<int>(j)] =
                    solution(j);
            }


            return ret;
        }
        static array<double>^ findSearchDirection2(
            int n,

            // =====================================================
            // Matrix A
            //
            // CSC format
            //
            // If unavailable, coefficient a is forced to zero.
            // =====================================================

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            // =====================================================
            // Matrix B
            //
            // CSC format
            //
            // If unavailable, coefficient b is forced to zero.
            // =====================================================

            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            // =====================================================
            // Basis vectors
            //
            //     P = [p0 p1 ...]
            //
            // Used exactly as supplied.
            //
            // No QR or orthonormalization is performed.
            // Linear dependence is handled by the TSVD of
            // the reduced matrix P^T M P.
            // =====================================================

            array<array<double>^>^ vecs,

            // =====================================================
            // RHS vectors
            // =====================================================

            array<double>^ rhs,
            array<double>^ rhs2,

            // =====================================================
            // Mixing coefficients
            //
            //     M = a A + b B
            //
            //     R = c rhs + d rhs2
            // =====================================================

            double a,
            double b,
            double c,
            double d,

            // =====================================================
            // Relative TSVD threshold
            // =====================================================

            double epsilon,

            // =====================================================
            // Output reduced coefficients used for this direction
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
            // Validate common parameters
            // =====================================================

            if (n <= 0)
            {
                throw gcnew ArgumentException(
                    "n must be positive");
            }

            if (vecs == nullptr)
            {
                throw gcnew ArgumentNullException(
                    "vecs");
            }

            if (!std::isfinite(a) ||
                !std::isfinite(b) ||
                !std::isfinite(c) ||
                !std::isfinite(d))
            {
                throw gcnew ArgumentException(
                    "a, b, c, and d must be finite");
            }

            if (!std::isfinite(epsilon) ||
                epsilon < 0.0)
            {
                throw gcnew ArgumentException(
                    "epsilon must be finite and non-negative");
            }

            parameters = nullptr;


            // =====================================================
            // Determine available inputs
            // =====================================================

            const bool hasA =
                columnptr != nullptr &&
                rowindices != nullptr &&
                values != nullptr;

            const bool hasB =
                columnptr2 != nullptr &&
                rowindices2 != nullptr &&
                values2 != nullptr;

            const bool hasRhs =
                rhs != nullptr;

            const bool hasRhs2 =
                rhs2 != nullptr;


            const double effectiveA =
                hasA ? a : 0.0;

            const double effectiveB =
                hasB ? b : 0.0;

            const double effectiveC =
                hasRhs ? c : 0.0;

            const double effectiveD =
                hasRhs2 ? d : 0.0;


            // =====================================================
            // Convert managed basis vectors to Eigen vectors
            // =====================================================

            std::vector<Eigen::VectorXd> inputBasis;

            inputBasis.reserve(
                static_cast<std::size_t>(
                    vecs->Length));

            for (int j = 0;
                j < vecs->Length;
                ++j)
            {
                array<double>^ src =
                    vecs[j];

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

                inputBasis.emplace_back(n);

                Eigen::VectorXd& vector =
                    inputBasis.back();

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

                    vector(i) =
                        value;
                }
            }


            // =====================================================
            // Keep the original basis exactly as supplied
            // =====================================================

            const Eigen::Index reducedSize =
                static_cast<Eigen::Index>(
                    inputBasis.size());

            if (reducedSize == 0)
            {
                parameters = gcnew array<double>(0);
                return gcnew array<double>(n);
            }


            // =====================================================
            // Number of nonzeros
            // =====================================================

            const int nnzA =
                hasA
                ? values->Length
                : 0;

            const int nnzB =
                hasB
                ? values2->Length
                : 0;


            // =====================================================
            // Local working storage
            //
            // Every invocation owns all mutable buffers. Therefore,
            // concurrent calls do not share Eigen matrices, vectors,
            // sparse matrices, or triplet arrays.
            // =====================================================

            Eigen::MatrixXd P(
                n,
                reducedSize);

            Eigen::MatrixXd MP(
                n,
                reducedSize);

            Eigen::MatrixXd H(
                reducedSize,
                reducedSize);

            Eigen::VectorXd g(
                reducedSize);

            Eigen::VectorXd combinedRhs(
                n);

            Eigen::VectorXd solution(
                reducedSize);

            Eigen::VectorXd transformedRhs(
                reducedSize);

            SparseMatrix64 matA(
                n,
                n);

            SparseMatrix64 matB(
                n,
                n);

            SparseMatrix64 matM(
                n,
                n);

            std::vector<Triplet64> tripletsA;
            std::vector<Triplet64> tripletsB;

            tripletsA.reserve(
                static_cast<std::size_t>(
                    nnzA));

            tripletsB.reserve(
                static_cast<std::size_t>(
                    nnzB));


            // =====================================================
            // Build dense basis matrix without QR
            //
            //     P = [p0 p1 ...]
            //
            // The original order, scale, zero columns, and linear
            // dependence are preserved. Rank deficiency is handled
            // later by the TSVD of H = P^T M P.
            // =====================================================

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                P.col(j) =
                    inputBasis[
                        static_cast<std::size_t>(j)];
            }


            // =====================================================
            // Build matrix A
            //
            // Triplet storage is local to this invocation.
            // =====================================================

            if (effectiveA != 0.0)
            {
                if (columnptr->Length != n + 1)
                {
                    throw gcnew ArgumentException(
                        "columnptr must have length n + 1");
                }

                if (rowindices->Length !=
                    values->Length)
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

                tripletsA.clear();

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
                            static_cast<int64_t>(
                                row),

                            static_cast<int64_t>(
                                column),

                            value);
                    }
                }

                // -------------------------------------------------
                // Rebuild numerical sparse matrix.
                //
                // Triplet buffer itself is reused.
                // -------------------------------------------------

                matA.setZero();

                matA.setFromTriplets(
                    tripletsA.begin(),
                    tripletsA.end());

                matA.makeCompressed();
            }
            else
            {
                matA.setZero();
            }


            // =====================================================
            // Build matrix B
            // =====================================================

            if (effectiveB != 0.0)
            {
                if (columnptr2->Length != n + 1)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 must have length n + 1");
                }

                if (rowindices2->Length !=
                    values2->Length)
                {
                    throw gcnew ArgumentException(
                        "rowindices2 and values2 have incompatible sizes");
                }

                if (columnptr2[0] != 0 ||
                    columnptr2[n] != values2->Length)
                {
                    throw gcnew ArgumentException(
                        "columnptr2 has invalid start or end");
                }

                tripletsB.clear();

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

                        tripletsB.emplace_back(
                            static_cast<int64_t>(
                                row),

                            static_cast<int64_t>(
                                column),

                            value);
                    }
                }

                matB.setZero();

                matB.setFromTriplets(
                    tripletsB.begin(),
                    tripletsB.end());

                matB.makeCompressed();
            }
            else
            {
                matB.setZero();
            }


            // =====================================================
            // Build mixed matrix
            //
            //     M = a A + b B
            // =====================================================

            if (effectiveA != 0.0 &&
                effectiveB != 0.0)
            {
                matM =
                    effectiveA * matA +
                    effectiveB * matB;
            }
            else if (effectiveA != 0.0)
            {
                matM =
                    effectiveA * matA;
            }
            else if (effectiveB != 0.0)
            {
                matM =
                    effectiveB * matB;
            }
            else
            {
                matM.setZero();
            }

            matM.makeCompressed();


            // =====================================================
            // Build mixed RHS
            //
            //     R = c rhs + d rhs2
            //
            // combinedRhs is local to this invocation.
            // =====================================================

            combinedRhs.setZero();

            if (effectiveC != 0.0)
            {
                if (rhs->Length != n)
                {
                    throw gcnew ArgumentException(
                        "rhs has incompatible size");
                }

                double* dst =
                    combinedRhs.data();

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

                    dst[i] +=
                        effectiveC *
                        value;
                }
            }

            if (effectiveD != 0.0)
            {
                if (rhs2->Length != n)
                {
                    throw gcnew ArgumentException(
                        "rhs2 has incompatible size");
                }

                double* dst =
                    combinedRhs.data();

                for (int i = 0;
                    i < n;
                    ++i)
                {
                    const double value =
                        rhs2[i];

                    if (!std::isfinite(value))
                    {
                        throw gcnew ArgumentException(
                            "rhs2 contains NaN or infinity");
                    }

                    dst[i] +=
                        effectiveD *
                        value;
                }
            }


            // =====================================================
            // FAST sparse-dense multiplication
            //
            //     MP = M P
            //
            // MP is local to this invocation.
            // =====================================================

            if (effectiveA == 0.0 &&
                effectiveB == 0.0)
            {
                MP.setZero();
            }
            else
            {
                MP.noalias() =
                    matM * P;
            }


            if (!MP.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "M * P contains NaN or infinity");
            }


            // =====================================================
            // FAST reduced matrix
            //
            //     H = P^T M P
            //
            // This replaces reducedSize^2 explicit VectorXd::dot()
            // operations with one dense matrix multiplication.
            // =====================================================

            H.noalias() =
                P.transpose() *
                MP;


            // =====================================================
            // Reduced RHS
            //
            //     g = P^T R
            // =====================================================

            g.noalias() =
                P.transpose() *
                combinedRhs;


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
            // TSVD
            //
            //     H = U S V^T
            //
            // Solve
            //
            //     H x = g
            //
            // with
            //
            //     sigma_i > epsilon * sigma_max
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


            solution.setZero();


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
                        epsilon *
                        sigmaMax;


                    transformedRhs.noalias() =
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


                    solution.noalias() =
                        svd.matrixV() *
                        transformedRhs;
                }
            }


            if (!solution.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Reduced solution contains NaN or infinity");
            }


            // =====================================================
            // Return the full-space direction
            //
            //     direction = P solution
            // =====================================================

            combinedRhs.noalias() =
                P * solution;


            if (!combinedRhs.allFinite())
            {
                throw gcnew InvalidOperationException(
                    "Full-space direction contains NaN or infinity");
            }

            array<double>^ ret =
                gcnew array<double>(n);

            parameters =
                gcnew array<double>(
                    static_cast<int>(reducedSize));


            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    combinedRhs(i);
            }

            for (Eigen::Index j = 0;
                j < reducedSize;
                ++j)
            {
                parameters[
                    static_cast<int>(j)] =
                    solution(j);
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
        
       
		static int solve_CHOLESKY(int n, array<int>^ columnptr, array<int>^ rowindices, array<double>^ values, array<double>^ rhs, array<double>^ ret,double scale,double salt,int numthreads)
		{
			

			Eigen::SparseMatrix<double> mat;
			Eigen::VectorXd b(n);
			Eigen::VectorXd _r(n);
			for (int i = 0; i < n; i++)
			{
				double v = rhs[i];
				b(i) = v*scale;
			}
			std::vector<Eigen::Triplet<double, int64_t>> tripletList;
			tripletList.reserve(values->Length+n);
			for (int i = 0; i < n; i++)
			{
				int start = columnptr[i];
				int end = columnptr[i + 1];
				for (int j = start; j < end; j++)
				{
					int row = rowindices[j];
					double val = values[j];
					tripletList.push_back(Eigen::Triplet<double, int64_t>(row, i, val*scale));
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

			_r=_sparse_solver::solve_CHOLECKY(mat, b,numthreads);
			for (int i = 0; i < n; i++)
			{
				ret[i]=_r(i);
			}

			return 1;
		}
        static int solve_CHOLESKY3(
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
            double epsilon,
            int numthreads)
        {
            const int totalSize = 2 * n;

            Eigen::SparseMatrix<double, 0, int64_t> mat;
            Eigen::VectorXd b(totalSize);
            Eigen::VectorXd solution(totalSize);


            // =====================================================
            // RHS
            //
            // [ rhs            ]
            // [ -epsilon rhs2  ]
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                b(i) =
                    rhs[i];

                b(n + i) =
                    -epsilon * rhs2[i];
            }


            // =====================================================
            // Matrix
            //
            // [ I      A          ]
            // [ A    -epsilon B   ]
            // =====================================================

            std::vector<Eigen::Triplet<double, int64_t>> tripletList;

            tripletList.reserve(
                n +
                2 * values->Length +
                values2->Length);


            // =====================================================
            // Upper-left: I
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        i,
                        i,
                        1.0));
            }


            // =====================================================
            // A blocks
            //
            // Upper-right: A
            // Lower-left : A
            // =====================================================

            for (int col = 0; col < n; col++)
            {
                int start =
                    columnptr[col];

                int end =
                    columnptr[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row =
                        rowindices[j];

                    double val =
                        values[j];


                    // Upper-right: A
                    tripletList.push_back(
                        Eigen::Triplet<double, int64_t>(
                            row,
                            n + col,
                            val));


                    // Lower-left: A
                    tripletList.push_back(
                        Eigen::Triplet<double, int64_t>(
                            n + row,
                            col,
                            val));
                }
            }


            // =====================================================
            // Lower-right: -epsilon B
            // =====================================================

            for (int col = 0; col < n; col++)
            {
                int start =
                    columnptr2[col];

                int end =
                    columnptr2[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row =
                        rowindices2[j];

                    double val =
                        values2[j];

                    tripletList.push_back(
                        Eigen::Triplet<double, int64_t>(
                            n + row,
                            n + col,
                            -epsilon * val));
                }
            }


            // =====================================================
            // Build matrix
            // =====================================================

            mat.resize(
                totalSize,
                totalSize);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();


            // =====================================================
            // Solve using PardisoLU
            //
            // Matrix is symmetric indefinite.
            // =====================================================

            Eigen::PardisoLU<
                Eigen::SparseMatrix<double,0,int64_t>>
                solver;

            solver.compute(mat);

            if (solver.info() != Eigen::Success)
            {
                throw gcnew InvalidOperationException(
                    "PardisoLU factorization failed");
            }

            solution =
                solver.solve(b);

            if (solver.info() != Eigen::Success)
            {
                throw gcnew InvalidOperationException(
                    "PardisoLU solve failed");
            }


            // =====================================================
            // solution =
            //
            // [ lambda ]
            // [ x      ]
            //
            // Return lambda and x
            // =====================================================

            for (int i = 0; i < totalSize; i++)
            {
                ret[i] = solution(i);
            }

            // kept for compatibility
            (void)numthreads;

            return 1;
        }
        static int solve_CHOLESKY4(
            int n,

            array<int>^ columnptr,
            array<int>^ rowindices,
            array<double>^ values,

            array<int>^ columnptr2,
            array<int>^ rowindices2,
            array<double>^ values2,

            array<double>^ rhs,
            array<double>^ rhs2,

            // =====================================================
            // q vector
            //
            // If nullptr:
            //     ordinary solve
            //
            // If non-null and qDotWeight > 0:
            //
            //     q^T direction ~= q^T q
            //
            // is weakly imposed.
            // =====================================================

            array<double>^ orthogonalVector,

            array<double>^ ret,

            // =====================================================
            // Base system
            //
            // M =
            //     matrixCoeff1 A
            //   + matrixCoeff2 B
            //   + epsilon I
            //
            // b =
            //     rhsCoeff1 rhs
            //   + rhsCoeff2 rhs2
            // =====================================================

            double matrixCoeff1,
            double matrixCoeff2,
            double rhsCoeff1,
            double rhsCoeff2,

            // =====================================================
            // Soft constraint weight
            //
            // rho / 2 *
            //
            //   (q^T direction / q^T q - 1)^2
            //
            // qDotWeight = 0:
            //     ordinary solve
            // =====================================================

            double qDotWeight,

            double epsilon,
            int numthreads)
        {
            using SparseMatrix =
                Eigen::SparseMatrix<
                double,
                Eigen::ColMajor,
                int64_t>;

            using Triplet =
                Eigen::Triplet<
                double,
                int64_t>;


            // =====================================================
            // Validate
            // =====================================================

            if (n <= 0)
            {
                return 0;
            }

            if (columnptr == nullptr ||
                rowindices == nullptr ||
                values == nullptr)
            {
                return 0;
            }

            if (columnptr2 == nullptr ||
                rowindices2 == nullptr ||
                values2 == nullptr)
            {
                return 0;
            }

            if (rhs == nullptr ||
                rhs2 == nullptr ||
                ret == nullptr)
            {
                return 0;
            }

            if (rhs->Length != n ||
                rhs2->Length != n ||
                ret->Length != n)
            {
                return 0;
            }

            if (!std::isfinite(matrixCoeff1) ||
                !std::isfinite(matrixCoeff2) ||
                !std::isfinite(rhsCoeff1) ||
                !std::isfinite(rhsCoeff2) ||
                !std::isfinite(qDotWeight) ||
                !std::isfinite(epsilon))
            {
                return 0;
            }

            if (qDotWeight < 0.0)
            {
                return 0;
            }


            // =====================================================
            // Build base RHS
            //
            //     b =
            //         rhsCoeff1 rhs
            //       + rhsCoeff2 rhs2
            // =====================================================

            Eigen::VectorXd b(n);

            for (int i = 0; i < n; ++i)
            {
                b(i) =
                    rhsCoeff1 * rhs[i]
                    + rhsCoeff2 * rhs2[i];
            }


            if (!b.allFinite())
            {
                return 0;
            }


            // =====================================================
            // Build sparse base matrix
            //
            //     M =
            //         matrixCoeff1 A
            //       + matrixCoeff2 B
            //       + epsilon I
            //
            // IMPORTANT:
            //
            // No q q^T term is added here.
            // =====================================================

            SparseMatrix mat(n, n);

            std::vector<Triplet> tripletList;

            tripletList.reserve(
                values->Length
                + values2->Length
                + n);


            // =====================================================
            // matrixCoeff1 A
            // =====================================================

            if (matrixCoeff1 != 0.0)
            {
                for (int col = 0;
                    col < n;
                    ++col)
                {
                    const int start =
                        columnptr[col];

                    const int end =
                        columnptr[col + 1];

                    for (int j = start;
                        j < end;
                        ++j)
                    {
                        const int row =
                            rowindices[j];

                        const double val =
                            matrixCoeff1 *
                            values[j];

                        tripletList.push_back(
                            Triplet(
                                row,
                                col,
                                val));
                    }
                }
            }


            // =====================================================
            // matrixCoeff2 B
            // =====================================================

            if (matrixCoeff2 != 0.0)
            {
                for (int col = 0;
                    col < n;
                    ++col)
                {
                    const int start =
                        columnptr2[col];

                    const int end =
                        columnptr2[col + 1];

                    for (int j = start;
                        j < end;
                        ++j)
                    {
                        const int row =
                            rowindices2[j];

                        const double val =
                            matrixCoeff2 *
                            values2[j];

                        tripletList.push_back(
                            Triplet(
                                row,
                                col,
                                val));
                    }
                }
            }


            // =====================================================
            // epsilon I
            // =====================================================

            if (epsilon != 0.0)
            {
                for (int i = 0;
                    i < n;
                    ++i)
                {
                    tripletList.push_back(
                        Triplet(
                            i,
                            i,
                            epsilon));
                }
            }


            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();


            // =====================================================
            // Factorize M only ONCE
            //
            // This is still the original sparse matrix.
            //
            // No dense rank-1 matrix is constructed.
            // =====================================================

            mkl_set_num_threads(
                numthreads);

            Eigen::PardisoLLT<
                SparseMatrix> solver;


            solver.analyzePattern(
                mat);

            if (solver.info() !=
                Eigen::Success)
            {
                return 0;
            }


            solver.factorize(
                mat);

            if (solver.info() !=
                Eigen::Success)
            {
                return 0;
            }


            // =====================================================
            // No q constraint
            //
            // or
            //
            // qDotWeight = 0
            //
            // Ordinary:
            //
            //     M d = b
            // =====================================================

            if (orthogonalVector == nullptr ||
                qDotWeight == 0.0)
            {
                Eigen::VectorXd solution =
                    solver.solve(
                        b);

                if (solver.info() !=
                    Eigen::Success)
                {
                    return 0;
                }

                if (!solution.allFinite())
                {
                    return 0;
                }


                for (int i = 0;
                    i < n;
                    ++i)
                {
                    ret[i] =
                        solution(i);
                }


                return 1;
            }


            // =====================================================
            // Read q
            // =====================================================

            if (orthogonalVector->Length != n)
            {
                return 0;
            }


            Eigen::VectorXd q(n);

            for (int i = 0;
                i < n;
                ++i)
            {
                const double value =
                    orthogonalVector[i];

                if (!std::isfinite(value))
                {
                    return 0;
                }

                q(i) =
                    value;
            }


            const double qNorm2 =
                q.squaredNorm();


            if (!std::isfinite(qNorm2) ||
                qNorm2 <= 0.0)
            {
                return 0;
            }


            // =====================================================
            // Normalized constraint vector
            //
            //     u = q / (q^T q)
            //
            //
            // Then:
            //
            //     u^T d ~= 1
            //
            // means exactly
            //
            //     q^T d ~= q^T q
            //
            //
            // Objective:
            //
            //     1/2 d^T M d
            //     - b^T d
            //
            //   + rho/2 (u^T d - 1)^2
            //
            //
            // Normal equation:
            //
            //     (M + rho u u^T)d
            //         = b + rho u
            //
            //
            // We DO NOT form u u^T.
            // =====================================================

            Eigen::VectorXd u(n);

            u.noalias() =
                q /
                qNorm2;


            // =====================================================
            // Solve two RHS with ONE factorization
            //
            //     x = M^-1 b
            //
            //     z = M^-1 u
            //
            // We can solve them simultaneously.
            // =====================================================

            Eigen::MatrixXd RHS(
                n,
                2);

            RHS.col(0) =
                b;

            RHS.col(1) =
                u;


            Eigen::MatrixXd X =
                solver.solve(
                    RHS);


            if (solver.info() !=
                Eigen::Success)
            {
                return 0;
            }


            if (!X.allFinite())
            {
                return 0;
            }


            const Eigen::VectorXd x =
                X.col(0);

            const Eigen::VectorXd z =
                X.col(1);


            // =====================================================
            // Sherman-Morrison
            //
            // Need to solve:
            //
            //     (M + rho u u^T)d
            //         = b + rho u
            //
            //
            // with
            //
            //     x = M^-1 b
            //     z = M^-1 u
            //
            //
            // Exact result:
            //
            //               rho (1 - u^T x)
            //     d = x + --------------------- z
            //               1 + rho u^T z
            //
            // =====================================================

            const double uTx =
                u.dot(x);

            const double uTz =
                u.dot(z);


            const double denominator =
                1.0 +
                qDotWeight *
                uTz;


            if (!std::isfinite(denominator) ||
                denominator == 0.0)
            {
                return 0;
            }


            const double correctionCoeff =
                qDotWeight *
                (1.0 - uTx) /
                denominator;


            if (!std::isfinite(correctionCoeff))
            {
                return 0;
            }


            Eigen::VectorXd solution(n);

            solution.noalias() =
                x +
                correctionCoeff *
                z;


            if (!solution.allFinite())
            {
                return 0;
            }


            // =====================================================
            // Return
            // =====================================================

            for (int i = 0;
                i < n;
                ++i)
            {
                ret[i] =
                    solution(i);
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
            array<double>^ rhs3,
            array<double>^ ret,
            double mu,
            double lambda,
            double epsilon,
            int numthreads)
        {
            const int totalSize = 2 * n;

            Eigen::SparseMatrix<double> mat;
            Eigen::VectorXd b(totalSize);
            Eigen::VectorXd solution(totalSize);


            // =====================================================
            // RHS
            //
            // [ rhs       + lambda rhs3 ]
            // [ mu rhs2   - lambda rhs3 ]
            //
            // rhs3 = current x - y
            //
            // update:
            //
            //     variable_new
            //         = variable_old - step * direction
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                b(i) =
                    rhs[i]
                    + lambda * rhs3[i];

                b(n + i) =
                    mu * rhs2[i]
                    - lambda * rhs3[i];
            }


            // =====================================================
            // Matrix
            //
            // [ A + (lambda+epsilon) I      -lambda I ]
            // [ -lambda I      mu B + (lambda+epsilon) I ]
            // =====================================================

            std::vector<Eigen::Triplet<double, int64_t>> tripletList;

            tripletList.reserve(
                values->Length +
                values2->Length +
                6 * n);


            // =====================================================
            // Upper-left: A
            // =====================================================

            for (int col = 0; col < n; col++)
            {
                int start = columnptr[col];
                int end = columnptr[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices[j];
                    double val = values[j];

                    tripletList.push_back(
                        Eigen::Triplet<double, int64_t>(
                            row,
                            col,
                            val));
                }
            }


            // =====================================================
            // Lower-right: mu B
            // =====================================================

            for (int col = 0; col < n; col++)
            {
                int start = columnptr2[col];
                int end = columnptr2[col + 1];

                for (int j = start; j < end; j++)
                {
                    int row = rowindices2[j];
                    double val = values2[j];

                    tripletList.push_back(
                        Eigen::Triplet<double, int64_t>(
                            n + row,
                            n + col,
                            mu * val));
                }
            }


            // =====================================================
            // lambda ||x-y||^2 coupling
            //
            // lambda *
            //
            // [ +I  -I ]
            // [ -I  +I ]
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                // Upper-left: +lambda I
                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        i,
                        i,
                        lambda));

                // Lower-right: +lambda I
                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        n + i,
                        n + i,
                        lambda));

                // Upper-right: -lambda I
                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        i,
                        n + i,
                        -lambda));

                // Lower-left: -lambda I
                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        n + i,
                        i,
                        -lambda));
            }


            // =====================================================
            // Numerical regularization
            //
            // epsilon I on both diagonal blocks.
            // Does NOT couple x and y.
            // =====================================================

            for (int i = 0; i < n; i++)
            {
                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        i,
                        i,
                        epsilon));

                tripletList.push_back(
                    Eigen::Triplet<double, int64_t>(
                        n + i,
                        n + i,
                        epsilon));
            }


            // =====================================================
            // Build matrix
            // =====================================================

            mat.resize(
                totalSize,
                totalSize);

            mat.setFromTriplets(
                tripletList.begin(),
                tripletList.end());

            mat.makeCompressed();


            // =====================================================
            // Solve
            // =====================================================

            solution =
                _sparse_solver::solve_CHOLECKY(
                    mat,
                    b,
                    numthreads);


            // =====================================================
            // Return [dx; dy]
            // ret must have length 2*n
            // =====================================================

            for (int i = 0; i < totalSize; i++)
            {
                ret[i] = solution(i);
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
