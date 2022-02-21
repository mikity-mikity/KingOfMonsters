#ifndef __CUSPARSE_WRAPPER_H__
#define __CUSPARSE_WRAPPER_H__

#include <cusparse.h>
#include <cusolverSp.h>
#include <cusolverSp_LOWLEVEL_PREVIEW.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
/* Description: Gather of non-zero elements from dense vector y into
   sparse vector x. */
/*cusparseStatus_t CUSPARSEAPI cusparseXgthr(cusparseHandle_t handle,
	int nnz,
	const float *y,
	float *xVal,
	const int *xInd,
	cusparseIndexBase_t idxBase);

cusparseStatus_t CUSPARSEAPI cusparseXgthr(cusparseHandle_t handle,
	int nnz,
	const double *y,
	double *xVal,
	const int *xInd,
	cusparseIndexBase_t idxBase);
*/
/*
 * Low level API for GPU luesky
 *
 */
cusolverStatus_t CUSOLVERAPI cusolverSpScsrluBufferInfoHost(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const float *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	size_t *internalDataInBytes,
	size_t *workspaceInBytes);

cusolverStatus_t CUSOLVERAPI cusolverSpDcsrluBufferInfoHost(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const double *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	size_t *internalDataInBytes,
	size_t *workspaceInBytes);

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluFactor(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const float *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	void *pBuffer);

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluFactor(
	cusolverSpHandle_t handle,
	int n,
	int nnzA,
	const cusparseMatDescr_t descrA,
	const double *csrValA,
	const int *csrRowPtrA,
	const int *csrColIndA,
	csrluInfoHost_t info,
	void *pBuffer);

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluZeroPivot(
	cusolverSpHandle_t handle,
	csrluInfoHost_t info,
	float tol,
	int *position);

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluZeroPivot(
	cusolverSpHandle_t handle,
	csrluInfoHost_t info,
	double tol,
	int *position);

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluSolve(
	cusolverSpHandle_t handle,
	int n,
	const float *b,
	float *x,
	csrluInfoHost_t info,
	void *pBuffer);

cusolverStatus_t CUSOLVERAPI cusolverSpXcsrluSolve(
	cusolverSpHandle_t handle,
	int n,
	const double *b,
	double *x,
	csrluInfoHost_t info,
	void *pBuffer);

#endif // !__CUSPARSE_WRAPPER_H__
