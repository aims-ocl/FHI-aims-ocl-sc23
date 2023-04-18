#ifndef CUDA_ERROR_H
#define CUDA_ERROR_H

#include <cublas_v2.h>
#include <cuda_runtime.h>

extern void cudaError(
      cudaError_t err,
      const char *file,
      int line);

extern void cublasError(
      cublasStatus_t err,
      const char *file,
      int line);

extern void checkError(
      const char *file,
      int line);

extern void promptError(
      const char* condition,
      const char* file,
      int line);

extern void promptError(
      const char* condition,
      const char* file,
      const char* info,
      int line);

extern void promptErrorVar(
      const char* condition,
      const int val,
      const char* file,
      int line);

extern void promptErrorVars(
      const char* condition,
      const int val1,
      const int val2,
      const char* file,
      int line);

#define HANDLE_CUBLAS( err ) (cublasError( err, __FILE__, __LINE__ ))
#define AIMS_EXIT() \
           printf("\n***ERROR in %s, line %i !\n",__FILE__,__LINE__); \
           exit( EXIT_FAILURE );
#define PROMPT_ERR(expr) \
           printf("\n***ERROR in %s, line %i !\n",__FILE__,__LINE__);\
           printf("%s\n\n",#expr); \
           exit( EXIT_FAILURE );
#define HANDLE_CUDA( err ) (cudaError( err, __FILE__, __LINE__ ))
#define CHECK_FOR_ERROR() (checkError( __FILE__, __LINE__ ))

// CHECK should be sunset by CHECK_INFO whenever possible.  For low-level
// functionality like aimsSetVector, the output of CHECK is essentially
// meaningless.  There needs to be an easy way to "tag" where and why the code
// is calling this function, which CHECK_INFO provides.
#define CHECK(expr) \
         if (!(expr)) {promptError(#expr,__FILE__,__LINE__);}
#define CHECK_INFO(expr,info) \
         if (!(expr)) {promptError(#expr,__FILE__,info,__LINE__);}
#define CHECK_VAR(expr,val) \
         if (!(expr)) {promptErrorVar(#expr,val,__FILE__,__LINE__);}
#define CHECK_VARS(expr,val1,val2) \
         if (!(expr)) {promptErrorVars(#expr,val1,val2,__FILE__,__LINE__);}

#endif /*CUDA_ERROR_H*/
