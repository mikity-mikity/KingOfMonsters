#include "cusparse_cholesky_solver.h"

#include <vector>
#include <algorithm>
#include <numeric>

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>

#include "device_buffer.h"
#include "cusparse_wrapper.h"

struct CusparseHandle
{
	CusparseHandle() { init(); }
	~CusparseHandle() { destroy(); }
	void init() { cusparseCreate(&handle); }
	void destroy() { cusparseDestroy(handle); }
	operator cusparseHandle_t() const { return handle; }
	CusparseHandle(const CusparseHandle&) = delete;
	CusparseHandle& operator=(const CusparseHandle&) = delete;
	cusparseHandle_t handle;
};

struct CusolverHandle
{
	CusolverHandle() { init(); }
	~CusolverHandle() { destroy(); }
	void init() { cusolverSpCreate(&handle); }
	void destroy() { cusolverSpDestroy(handle); }
	operator cusolverSpHandle_t() const { return handle; }
	CusolverHandle(const CusolverHandle&) = delete;
	CusolverHandle& operator=(const CusolverHandle&) = delete;
	cusolverSpHandle_t handle;
};

struct CusparseMatDescriptor
{
	CusparseMatDescriptor() { init(); }
	~CusparseMatDescriptor() { destroy(); }

	void init()
	{
		cusparseCreateMatDescr(&desc);
		cusparseSetMatType(desc, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(desc, CUSPARSE_INDEX_BASE_ZERO);
		cusparseSetMatDiagType(desc, CUSPARSE_DIAG_TYPE_NON_UNIT);		
	}

	void destroy() { cusparseDestroyMatDescr(desc); }
	operator cusparseMatDescr_t() const { return desc; }
	CusparseMatDescriptor(const CusparseMatDescriptor&) = delete;
	CusparseMatDescriptor& operator=(const CusparseMatDescriptor&) = delete;
	cusparseMatDescr_t desc;
};

template <typename T>
class SparseSquareMatrixCSR
{

public:

	SparseSquareMatrixCSR() : size_(0), nnz_(0) {}

	void resize(int size)
	{
		size_ = size;
		rowPtr_.allocate(size + 1);
	}

	void resizeNonZeros(int nnz)
	{
		nnz_ = nnz;
		values_.allocate(nnz);
		colInd_.allocate(nnz);
	}

	void upload(const T* values = nullptr, const int* rowPtr = nullptr, const int* colInd = nullptr)
	{
		if (values)
			values_.upload(values);
		if (rowPtr)
			rowPtr_.upload(rowPtr);
		if (colInd)
			colInd_.upload(colInd);
	}

	T* val() { return values_.data; }
	int* rowPtr() { return rowPtr_.data; }
	int* colInd() { return colInd_.data; }

	const T* val() const { return values_.data; }
	const int* rowPtr() const { return rowPtr_.data; }
	const int* colInd() const { return colInd_.data; }

	int size() const { return size_; }
	int nnz() const { return nnz_; }

	cusparseMatDescr_t desc() const { return desc_; }

public:

	DeviceBuffer<T> values_;
	DeviceBuffer<int> rowPtr_;
	DeviceBuffer<int> colInd_;
	int size_, nnz_;
	CusparseMatDescriptor desc_;
};

template <typename T>
class Sparsecholesky
{

public:

	void init(cusolverSpHandle_t handle)
	{
		handle_ = handle;

		// create info
		cusolverSpCreateCsrluInfoHost(&info_);
	}

	std::string allocateBuffer(const SparseSquareMatrixCSR<T>& A, const int* csrRowPtr, const int* csrColInd)
	{
		size_t internalData, workSpace;
		cusolverSpDcsrluBufferInfoHost(handle_, A.size(), A.nnz(), A.desc(),
			(double*) A.val(), csrRowPtr, csrColInd, info_, &internalData, &workSpace);
		//std::cout << "cusparse workspace size=" << workSpace << std::endl;
		std::string res=buffer_.allocate(workSpace);
		return res;
	}

	bool hasZeroPivot(int* position = nullptr) const
	{
		const T tol = static_cast<T>(1e-14);
		int singularity = -1;
		cusolverSpXcsrluZeroPivot(handle_, info_, tol, &singularity);
		if (position)
			*position = singularity;
		return singularity >= 0;
	}

	std::string analyze(const SparseSquareMatrixCSR<T>& A, const int* csrRowPtr, const int* csrColInd)
	{
		//auto err=cusolverSpXcsrluAnalysisHost(handle_, A.size(), A.nnz(), A.desc(), csrRowPtr, csrColInd, info_);
		//std::string res=allocateBuffer(A, csrRowPtr, csrColInd );
		//return res;
		return "p";
	}

	std::string factorize(SparseSquareMatrixCSR<T>& A, double* val,  const int* csrRowPtr, const int* csrColInd,double* rhs, double* ret,int ordering)
	{
		/*cusolverSpXcsrluFactor(handle_, A.size(), A.nnz(), A.desc(),
			val, csrRowPtr, csrColInd, info_, buffer_.data);
		int zeropivot = -1;
		if (!hasZeroPivot(&zeropivot))return "SUCCESS";
		*/
		int singularity = -1;
		auto err = cusolverSpDcsrlsvchol(

			handle_, A.size(), A.nnz(), A.desc(),
			val, csrRowPtr, csrColInd,
			rhs,
			0.0000000001,
			ordering,
			ret,
			&singularity);
		if (err == 2)
		{
			return "ALLOCERR";
		}
		if (singularity >= 0)
		{
			return "PIVOT";// "found pivot = " + std::to_string(singularity) + "err " + std::to_string(err) + "nnz=" + std::to_string(A.nnz());
		}
		{
			return "SUCCESS";
		}
	}
	std::string factorize_cpu(SparseSquareMatrixCSR<T>& A, double* val, const int* csrRowPtr, const int* csrColInd, double* rhs, double* ret, int ordering)
	{
		/*cusolverSpXcsrluFactor(handle_, A.size(), A.nnz(), A.desc(),
			val, csrRowPtr, csrColInd, info_, buffer_.data);
		int zeropivot = -1;
		if (!hasZeroPivot(&zeropivot))return "SUCCESS";
		*/
		int singularity = -1;
		auto err = cusolverSpDcsrlsvcholHost(

			handle_, A.size(), A.nnz(), A.desc(),
			val, csrRowPtr, csrColInd,
			rhs,
			0.0000000001,
			ordering,
			ret,
			&singularity);
		if (err == 2)
		{
			return "ALLOCERR";
		}
		if (singularity >= 0)
		{
			return "PIVOT";// "found pivot = " + std::to_string(singularity) + "err " + std::to_string(err) + "nnz=" + std::to_string(A.nnz());
		}
		{
			return "SUCCESS";
		}
	}
	void solve(int size, const T* b, T* x)
	{
		cusolverSpXcsrluSolve(handle_, size, b, x, info_, (void*)buffer_.data);
	}

	void destroy()
	{
		cusolverSpDestroyCsrluInfoHost(info_);
	}

	~Sparsecholesky() { destroy(); }

private:

	cusolverSpHandle_t handle_;
	csrluInfoHost_t info_;
	DeviceBuffer<unsigned char> buffer_;
};

class Twist
{

public:

	Twist()
	{
		cusolverSpCreate(&handle_);
		cusparseCreateMatDescr(&desc_);
		cusparseSetMatType(desc_, CUSPARSE_MATRIX_TYPE_GENERAL);
		cusparseSetMatIndexBase(desc_, CUSPARSE_INDEX_BASE_ZERO);
	}

	~Twist()
	{
		cusolverSpDestroy(handle_);
		cusparseDestroyMatDescr(desc_);
	}

	void allocateBuffer(int n, int nnz, const int* csrRowPtr, const int* csrColInd, const int* P)
	{
		size_t bufSize;
		cusolverSpXcsrperm_bufferSizeHost(handle_, n, n, nnz, desc_,
			csrRowPtr, csrColInd, P, P, &bufSize);
		buffer_.resize(bufSize);
	}

	void compute(int n, int nnz, const int* csrRowPtr, const int* csrColInd, std::vector<int>& P,
		std::vector<int>& permRowPtr, std::vector<int>& permColInd, std::vector<int>& map)
	{
		const int* p = P.data();

		allocateBuffer(n, nnz, csrRowPtr, csrColInd, p);

		permRowPtr.resize(n + 1);
		permColInd.resize(nnz);
		map.resize(nnz);

		std::copy(csrRowPtr, csrRowPtr + (n + 1), std::begin(permRowPtr));
		std::copy(csrColInd, csrColInd + nnz, std::begin(permColInd));
		std::iota(std::begin(map), std::end(map), 0);

		cusolverSpXcsrpermHost(handle_, n, n, nnz, desc_,
			permRowPtr.data(), permColInd.data(), p, p, map.data(), buffer_.data());
	}

	void operator()(int n, int nnz, const int* csrRowPtr, const int* csrColInd, std::vector<int>& P,
		std::vector<int>& permRowPtr, std::vector<int>& permColInd, std::vector<int>& map)
	{
		compute(n, nnz, csrRowPtr, csrColInd, P, permRowPtr, permColInd, map);
	}

private:

	cusolverSpHandle_t handle_;
	cusparseMatDescr_t desc_;
	std::vector<unsigned char> buffer_;
};

template <typename T>
class CuSparseCholeskySolverImpl : public CuSparseCholeskySolver<T>
{

public:

	using Info = typename CuSparseCholeskySolver<T>::Info;

	CuSparseCholeskySolverImpl(int size)
	{
		init();

		if (size > 0)
			resize(size);
	}

	void init()
	{
		cholesky.init(cusolver);
		doOrdering = false;
		information = Info::SUCCESS;
	}

	void resize(int size) override
	{
		Acsr.resize(size);
		d_b.allocate(size);
		d_x.allocate(size);
		d_y.allocate(size);
	}

	void setPermutaion(int size, const int* P) override
	{
		h_P.resize(size);
		h_PT.resize(size);

		for (int i = 0; i < size; i++)
		{
			h_P[i] = P[i];
			h_PT[P[i]] = i;
		}

		d_P.assign(size, h_P.data());
		d_PT.assign(size, h_PT.data());
		d_z.allocate(size);

		doOrdering = true;
	}

	std::string analyze(int nnz, double* vals, const int* csrRowPtr, const int* csrColInd) override
	{
		// if permutation is set, apply permutation to A (P*A*PT)
		//if (doOrdering)
		//{
		//	twist(Acsr.size(), nnz, csrRowPtr, csrColInd, h_P, permRowPtr, permColInd, h_map);
		//	d_map.assign(nnz, h_map.data());
		//	csrRowPtr = permRowPtr.data();
		//	csrColInd = permColInd.data();
		//}

		// copy input data to device memory
		Acsr.resizeNonZeros(nnz);
		Acsr.upload((T*)vals, csrRowPtr, csrColInd);

		return cholesky.analyze(Acsr, csrRowPtr, csrColInd);
	}
	std::string analyze_cpu(int nnz) override
	{
		Acsr.resizeNonZeros(nnz);
		return "";
	}

	std::string factorize(const T* A,const int* csrRowPtr, const int* csrColInd,double *rhs,double *ret,int ordering) override
	{
		double a1[3];
		d_b.upload((T*)rhs);

		/*if (doOrdering)
		{
			d_values.assign(Acsr.nnz(), A);
			permute(Acsr.nnz(), d_values.data, Acsr.val(), d_map.data);
		}
		else*/
		{
			//Acsr.upload(A);
			
			CUDA_CHECK(cudaMemcpy(&a1, Acsr.values_.data, sizeof(T) * 3, cudaMemcpyDeviceToHost));

		}

		// M = L * LT
		std::string fff = cholesky.factorize(Acsr, (double*) Acsr.val(), Acsr.rowPtr(), Acsr.colInd(),(double*)d_b.data,(double*)d_x.data,ordering);

		if (fff== "SUCCESS")
		{
			information = Info::SUCCESS;
			d_x.download((T*)ret);
			return "FACTORIZE SUCCESS" +fff;
		}
		else {
			information = Info::NUMERICAL_ISSUE;
			return "FACTORIZE FAILED"+fff+","+std::to_string(a1[0])+","+std::to_string(a1[1])+ "," + std::to_string(a1[2]);
		}
	}
	std::string factorize_cpu(const T* A, const int* csrRowPtr, const int* csrColInd, double* rhs, double* ret, int ordering) override
	{
		double a1[3];
		//d_b.upload((T*)rhs);

		/*if (doOrdering)
		{
			d_values.assign(Acsr.nnz(), A);
			permute(Acsr.nnz(), d_values.data, Acsr.val(), d_map.data);
		}
		else*/
		{
			//Acsr.upload(A);

			//CUDA_CHECK(cudaMemcpy(&a1, Acsr.values_.data, sizeof(T) * 3, cudaMemcpyDeviceToHost));

		}

		// M = L * LT
		std::string fff = cholesky.factorize_cpu(Acsr, (double*)A, csrRowPtr, csrColInd, (double*)rhs, (double*)ret, ordering);

		if (fff == "SUCCESS")
		{
			information = Info::SUCCESS;
			//d_x.download((T*)ret);
			return "FACTORIZE SUCCESS" + fff;
		}
		else {
			information = Info::NUMERICAL_ISSUE;
			return "FACTORIZE FAILED" + fff + "," + std::to_string(a1[0]) + "," + std::to_string(a1[1]) + "," + std::to_string(a1[2]);
		}
	}
	void solve(const T* b, T* x) override
	{
		d_b.upload(b);

		if (doOrdering)
		{
			// y = P * b
			permute(Acsr.size(), d_b.data, d_y.data, d_P.data);

			// solve A * z = y
			cholesky.solve(Acsr.size(), d_y.data, d_z.data);

			// x = PT * z
			permute(Acsr.size(), d_z.data, d_x.data, d_PT.data);
		}
		else
		{
			// solve A * x = b
			cholesky.solve(Acsr.size(), d_b.data, d_x.data);
		}

		d_x.download(x);
	}

	Info info() const override
	{
		return information;
	}

	void permute(int size, const T* src, T* dst, const int* P)
	{
		//cusparseXgthr(cusparse, size, src, dst, P, CUSPARSE_INDEX_BASE_ZERO);
	}

	void destroy()
	{
	}

	~CuSparseCholeskySolverImpl()
	{
		destroy();
	}


	public:
		SparseSquareMatrixCSR<T> Acsr;

private:
	DeviceBuffer<T> d_b;
	DeviceBuffer<T> d_x;
	DeviceBuffer<T> d_y;
	DeviceBuffer<T> d_z;
	DeviceBuffer<T> d_values;
	DeviceBuffer<int> d_P, d_PT, d_map;

	CusparseHandle cusparse;
	CusolverHandle cusolver;

	Sparsecholesky<T> cholesky;

	Twist twist;
	std::vector<int> h_P, h_PT, h_map, permRowPtr, permColInd;
	bool doOrdering;

	Info information;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename T>
typename CuSparseCholeskySolver<T>::Ptr CuSparseCholeskySolver<T>::create(int size)
{
	return std::make_unique<CuSparseCholeskySolverImpl<T>>(size);
}

template<typename T>
CuSparseCholeskySolver<T>::~CuSparseCholeskySolver()
{
}

template class CuSparseCholeskySolver<double>;
template class CuSparseCholeskySolver<float>;
