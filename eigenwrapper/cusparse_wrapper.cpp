#include "cusparse_wrapper.h"

/* Description: Gather of non-zero elements from dense vector y into
   sparse vector x. */
//cusparseStatus_t CUSPARSEAPI cusparseXgthr(cusparseHandle_t handle,
//	int nnz,/
// const float *y,
//	float *xVal,
//	const int *xInd,
//	cusparseIndexBase_t idxBase)
//{
//	return
//		cusparseSgthr(handle,
//			nnz,
//			y,
//			xVal,
//			xInd,
//			idxBase);
//}

/*
 * Low level API for GPU luesky
 *
 */
//cusparseStatus_t CUSPARSEAPI cusparseXgthr(cusparseHandle_t handle,
//	int nnz,
//	const double *y,
//	double *xVal,
//	const int *xInd,
//	cusparseIndexBase_t idxBase)
//{
	//cusparseDnVecDescr_t vecY;
	//cusparseSpVecDescr_t vecX;
	//CHECK_CUSPARSE(cusparseCreateSpVec(&vecX, size, nnz, dX_indices, dX_values,
//		CUSPARSE_INDEX_32I,
//		CUSPARSE_INDEX_BASE_ZERO, CUDA_R_32F));
		// Create dense vector y
	//CHECK_CUSPARSE(cusparseCreateDnVec(&vecY, size, dY, CUDA_R_32F));
	//*/
	//return CUSPARSE_STATUS_SUCCESS;
		//cusparseGather(handle, vecY, vecX);
		/*cusparseDgthr(handle,
			nnz,
			y,
			xVal,
			xInd,
			idxBase);*/
//}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluBufferInfo(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const float *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	size_t *internalDataInBytes,
	size_t *workspaceInBytes)
{
	return
		cusolverSpScsrluBufferInfoHost(
			handle,
			n,
			nnzA,
			descrA,
			csrValA,
			csrRowPtrA,
			csrColIndA,
			info,
			internalDataInBytes,
			workspaceInBytes);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluBufferInfo(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const double *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	size_t *internalDataInBytes,
	size_t *workspaceInBytes)
{
	return
		cusolverSpDcsrluBufferInfoHost(
			handle,
			n,
			nnzA,
			descrA,
			csrValA,
			csrRowPtrA,
			csrColIndA,
			info,
			internalDataInBytes,
			workspaceInBytes);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluFactor(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const float *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	void *pBuffer)
{
	return
		cusolverSpScsrluFactorHost(handle,
			n,
			nnzA,
			descrA,
			csrValA,
			csrRowPtrA,
			csrColIndA,
			info,
			0.00000001,
			pBuffer);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluFactor(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const double *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	void *pBuffer)
{
	return
		cusolverSpDcsrluFactorHost(handle,
			n,
			nnzA,
			descrA,
			csrValA,
			csrRowPtrA,
			csrColIndA,
			info,
			0.0000000001,
			pBuffer);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluZeroPivot(
	cusolverSpHandle_t handle,
	csrluInfoHost_t info,
	float tol,
	int *position)
{
	return
		cusolverSpScsrluZeroPivotHost(
			handle,
			info,
			tol,
			position);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrqrZeroPivot(
	cusolverSpHandle_t handle,
	csrqrInfoHost_t info,
	double tol,
	int *position)
{
	return
		cusolverSpDcsrqrZeroPivotHost(
			handle,
			info,
			tol,
			position);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrqrSolve    (
	cusolverSpHandle_t handle,
	int n,
	int m,
	float *b,
	float *x,
	csrqrInfoHost_t info,
	void *pBuffer)
{
	return
		cusolverSpScsrqrSolveHost(
			handle,
			n,
			m,
			b,
			x,
			info,
			pBuffer);
}

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrqrSolve(
	cusolverSpHandle_t handle,
	int n,
	int m,
	double *b,
	double *x,
	csrqrInfoHost_t info,
	void *pBuffer)
{
	return
		cusolverSpDcsrqrSolveHost(
			handle,
			n,
			m,
			b,
			x,
			info,
			pBuffer);
}
