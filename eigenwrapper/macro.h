#ifndef __MACRO_H__
#define __MACRO_H__

#include <cstdio>
# ifndef _CPU
#define CUDA_CHECK(err) \
do {\
	if (err != cudaSuccess) { \
		printf("[CUDA Error] %s (code: %d) at %s:%d\n", cudaGetErrorString(err), err, __FILE__, __LINE__); \
	} \
} while (0)

# endif // !__MACRO_H__
#endif