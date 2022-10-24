#ifndef _CUDA_UTIL_H_
#define _CUDA_UTIL_H_

#include <cstdio>
#include <cstdlib>             /* EXIT_FAILURE */

#define CUDA_SAFE_CALL(call) {                                            \
	hipError_t err = call;                                                    \
	if(hipSuccess != err) {                                                \
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
				 __FILE__, __LINE__, hipGetErrorString( err) );             \
		exit(EXIT_FAILURE);                                                  \
	} }

#define CUDA_SAFE_CALL_SYNC(call) {                                       \
	CUDA_SAFE_CALL_NO_SYNC(call);                                            \
	hipError_t err |= hipDeviceSynchronize();                                \
	if(hipSuccess != err) {                                                \
		fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
				__FILE__, __LINE__, hipGetErrorString( err) );              \
		exit(EXIT_FAILURE);                                                  \
	} }

#define CUDA_CHECK_ERROR(errorMessage) {                                    \
    hipError_t err = hipGetLastError();                                    \
    if( hipSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, hipGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    } }

#endif
