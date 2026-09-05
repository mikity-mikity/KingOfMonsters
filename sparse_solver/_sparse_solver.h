#pragma once
#include "eigen-3.4.0/Eigen/PardisoSupport"
#include "eigen-3.4.0/Eigen/Sparse"
#include "eigen-3.4.0/Eigen/SVD"
#include "eigen-3.4.0/Eigen/Dense"
#include "eigen-3.4.0/Eigen/SparseQR"
#include "eigen-3.4.0/Eigen/SparseLU"
#include "eigen-3.4.0/Eigen/SparseCholesky"	

#include <limits>
#include <cmath>
#include <algorithm>

#define EIGEN_NO_DEBUG
#define EIGEN_NO_STATIC_ASSERT
#define EIGEN_USE_LAPACK
#define EIGEN_USE_MKL_ALL
//#define EIGEN_DONT_ALIGN_STATICALLY
//#define EIGEN_MAX_ALIGN_BYTES 0
//#define EIGEN_DONT_VECTORIZE
#define EIGEN_DEFAULT_DENSE_INDEX_TYPE int64_t
using SpMat = Eigen::SparseMatrix<double, 0, int64_t>;
using Vec = Eigen::VectorXd;

namespace _sparse_solver
{
 
	Eigen::MatrixXcd genEigen2(Eigen::MatrixXd ma, Eigen::MatrixXd mb, double* l1, double* l2, double* l1i, double* l2i);
    std::vector<Eigen::VectorXd>
        buildOrthonormalBasis(
            const std::vector<Eigen::VectorXd>& inputVecs,
            Eigen::Index fullSize,
            double relativeThreshold);

    std::array<Eigen::VectorXd, 4> findSearchDirection3(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection4(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection5(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection6(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection7(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection8(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );

    std::array<Eigen::VectorXd, 4> findSearchDirection9(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection10(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection11(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
    std::array<Eigen::VectorXd, 4> findSearchDirection12(
        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat,

        const Eigen::SparseMatrix<
        double,
        Eigen::ColMajor,
        int64_t>& mat2,

        const std::vector<Eigen::VectorXd>& vecs,
        const std::vector<Eigen::VectorXd>& nearNullVecs,

        const Eigen::VectorXd& rhs,
        const Eigen::VectorXd& rhs2,
        std::vector<double>& projectedBasisNorms,
        double secondaryWeight,
        double couplingWeight,
        double epsilon,
        bool removeNearNullFromPrimary
    );
        std::array<Eigen::VectorXd, 4> findSearchDirection13(
            const Eigen::SparseMatrix<
            double,
            Eigen::ColMajor,
            int64_t>& mat,

            const Eigen::SparseMatrix<
            double,
            Eigen::ColMajor,
            int64_t>& mat2,

            const std::vector<Eigen::VectorXd>& vecs,
            const std::vector<Eigen::VectorXd>& nearNullVecs,

            const Eigen::VectorXd& rhs,
            const Eigen::VectorXd& rhs2,
            std::vector<double>& projectedBasisNorms,
            double secondaryWeight,
            double couplingWeight,
            double epsilon,
            bool removeNearNullFromPrimary
        );
    Eigen::VectorXd solve_CHOLECKY(Eigen::SparseMatrix<double, 0, int64_t> mat, Eigen::VectorXd rhs,int numthreads);
	Eigen::VectorXd solve_projection(Eigen::SparseMatrix<double, 0, int64_t> mat, Eigen::SparseMatrix<double, 0, int64_t>, Eigen::VectorXd rhs);
	Vec apply_near_null_filter(const SpMat& A, const Vec& v, double mu);
}
