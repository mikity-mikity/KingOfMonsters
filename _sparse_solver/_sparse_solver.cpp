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
    Eigen::VectorXd solveReducedNonnegativeActiveSet(
        const std::vector<Eigen::VectorXd>& vecs,
        const Eigen::MatrixXd& mat,
        const Eigen::VectorXd& rhs,
        int maxIter = 25,
        double tol = 1e-10,
        double salt = 1e-10)
    {
        const int n = static_cast<int>(vecs.size());
        Eigen::MatrixXd M(n, n);
        Eigen::VectorXd rhs2(n);

        // mat * vecs[j] を先に計算して使い回す
        std::vector<Eigen::VectorXd> mvecs(n);
        for (int j = 0; j < n; ++j) {
            mvecs[j] = mat * vecs[j];
        }

        for (int i = 0; i < n; ++i) {
            rhs2(i) = vecs[i].dot(rhs);
            for (int j = 0; j < n; ++j) {
                M(i, j) = vecs[i].dot(mvecs[j]);
            }
        }

        M = 0.5 * (M + M.transpose());
        M.diagonal().array() += salt;

        std::vector<char> active(n, 0);
        Eigen::VectorXd alpha = Eigen::VectorXd::Zero(n);

        for (int iter = 0; iter < maxIter; ++iter) {
            std::vector<int> freeIds;
            freeIds.reserve(n);
            for (int i = 0; i < n; ++i) {
                if (!active[i]) freeIds.push_back(i);
            }

            if (freeIds.empty()) {
                alpha.setZero();
                break;
            }

            const int nf = static_cast<int>(freeIds.size());
            Eigen::MatrixXd Mff(nf, nf);
            Eigen::VectorXd bf(nf);

            for (int r = 0; r < nf; ++r) {
                const int rr = freeIds[r];
                bf(r) = rhs2(rr);
                for (int c = 0; c < nf; ++c) {
                    Mff(r, c) = M(rr, freeIds[c]);
                }
            }

            // FullPivLU より軽い
            Eigen::VectorXd alphaF;
            Eigen::LLT<Eigen::MatrixXd> llt(Mff);
            if (llt.info() == Eigen::Success) {
                alphaF = llt.solve(bf);
            }
            else {
                Eigen::LDLT<Eigen::MatrixXd> ldlt(Mff);
                alphaF = ldlt.solve(bf);
            }

            Eigen::VectorXd newAlpha = Eigen::VectorXd::Zero(n);
            for (int k = 0; k < nf; ++k) {
                newAlpha(freeIds[k]) = alphaF(k);
            }

            bool addedToActive = false;
            for (int i = 0; i < n; ++i) {
                if (!active[i] && newAlpha(i) < -tol) {
                    active[i] = 1;
                    addedToActive = true;
                }
            }
            if (addedToActive) continue;

            alpha = newAlpha;

            Eigen::VectorXd grad = M * alpha - rhs2;

            int releaseIdx = -1;
            double mostNegative = -tol;
            for (int i = 0; i < n; ++i) {
                if (active[i] && grad(i) < mostNegative) {
                    mostNegative = grad(i);
                    releaseIdx = i;
                }
            }

            if (releaseIdx >= 0) {
                active[releaseIdx] = 0;
                continue;
            }

            break;
        }

        alpha = alpha.cwiseMax(0.0);
        return alpha;
    }
	Eigen::VectorXd findSearchDirection(Eigen::SparseMatrix<double, 0, int64_t> mat, std::vector<Eigen::VectorXd> vecs,Eigen::VectorXd rhs,double salt)
	{
		


        Eigen::VectorXd sol = solveReducedNonnegativeActiveSet(vecs, mat, rhs, 40, 1e-10,salt);
        return sol;

        /*int n = vecs.size();
        Eigen::MatrixXd m(n, n);
		for (int i = 0; i < n; i++)
        {
			for (int j = 0; j < n; j++)
            {
                m(i, j) = vecs[i].dot(mat * vecs[j]);
            }
        }
		Eigen::FullPivLU <Eigen::MatrixXd> lu(m);
        Eigen::VectorXd rhs2(n);
        for (int i = 0; i < n; i++)
        {
            rhs2(i) = vecs[i].dot(rhs);
		}
        Eigen::VectorXd sol = lu.solve(rhs2);
		return sol;*/

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
		for (int ss = 0; ss < 10; ss++)
		{
			chol.compute(mat);
			if (chol.info() == Eigen::Success) {
				break;
			}
			else
				{
                Eigen::VectorXd _r(n);
                _r.setZero();
                return _r;
				}
		}
		Eigen::VectorXd _r(n);

		_r = chol.solve(rhs);
		return _r;
	}
  

}
