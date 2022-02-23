#ifndef __CUSPARSE_CHOLESKY_SOLVER_H__
#define __CUSPARSE_CHOLESKY_SOLVER_H__

#include <memory>
#include <string>
template <typename T>
class CuSparseCholeskySolver
{
public:

	enum Info
	{
		SUCCESS,
		NUMERICAL_ISSUE
	};

	using Ptr = std::unique_ptr<CuSparseCholeskySolver>;

	static Ptr create(int size = 0);
	virtual void resize(int size) = 0;
	virtual std::string analyze(int nnz, double* val, const int* csrRowPtr, const int* csrColInd) = 0;
	virtual std::string analyze_cpu(int nnz) = 0;

	virtual std::string factorize(const T* A, const int* csrRowPtr, const int* csrColInd, double* rhs, double* ret, int ordering) = 0;
	virtual std::string factorize_cpu(const T* A, const int* csrRowPtr, const int* csrColInd, double* rhs, double* ret, int ordering) = 0;
	virtual void solve(const T* b, T* x) = 0;
	virtual void setPermutaion(int size, const int* P) = 0;
	virtual Info info() const = 0;

	virtual ~CuSparseCholeskySolver();
};

#endif // !__CUSPARSE_CHOLESKY_SOLVER_H__
