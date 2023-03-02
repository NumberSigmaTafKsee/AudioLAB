// todo: table <=> array, table <=> vector, table <=> matrix etc
// csv read/write


%module vf
%{
#include "viperfish.h"
%}
%include "viperfish.h"


//%include "stdint.i"
typedef signed char int8_t;
typedef unsigned char uint8_t;
typedef short int16_t;
typedef unsigned short uint16_t;
typedef int int32_t;
typedef unsigned int uint32_t;
typedef long int int64_t;
typedef unsigned long int uint64_t;


%include "std_vector.i"

%template(float_vector) std::vector<float>;
%template(double_vector) std::vector<double>;
%template(complex_vector) std::vector<cuComplex>;
%template(double_complex_vector) std::vector<cuDoubleComplex>;

/*
%include "cuda_runtime.h"
%include "cublas_v2.h"
*/

typedef enum cudaDataType_t
{
    CUDA_R_16F  =  2, /* real as a half */
    CUDA_C_16F  =  6, /* complex as a pair of half numbers */
    CUDA_R_16BF = 14, /* real as a nv_bfloat16 */
    CUDA_C_16BF = 15, /* complex as a pair of nv_bfloat16 numbers */
    CUDA_R_32F  =  0, /* real as a float */
    CUDA_C_32F  =  4, /* complex as a pair of float numbers */
    CUDA_R_64F  =  1, /* real as a double */
    CUDA_C_64F  =  5, /* complex as a pair of double numbers */
    CUDA_R_4I   = 16, /* real as a signed 4-bit int */
    CUDA_C_4I   = 17, /* complex as a pair of signed 4-bit int numbers */
    CUDA_R_4U   = 18, /* real as a unsigned 4-bit int */
    CUDA_C_4U   = 19, /* complex as a pair of unsigned 4-bit int numbers */
    CUDA_R_8I   =  3, /* real as a signed 8-bit int */
    CUDA_C_8I   =  7, /* complex as a pair of signed 8-bit int numbers */
    CUDA_R_8U   =  8, /* real as a unsigned 8-bit int */
    CUDA_C_8U   =  9, /* complex as a pair of unsigned 8-bit int numbers */
    CUDA_R_16I  = 20, /* real as a signed 16-bit int */
    CUDA_C_16I  = 21, /* complex as a pair of signed 16-bit int numbers */
    CUDA_R_16U  = 22, /* real as a unsigned 16-bit int */
    CUDA_C_16U  = 23, /* complex as a pair of unsigned 16-bit int numbers */
    CUDA_R_32I  = 10, /* real as a signed 32-bit int */
    CUDA_C_32I  = 11, /* complex as a pair of signed 32-bit int numbers */
    CUDA_R_32U  = 12, /* real as a unsigned 32-bit int */
    CUDA_C_32U  = 13, /* complex as a pair of unsigned 32-bit int numbers */
    CUDA_R_64I  = 24, /* real as a signed 64-bit int */
    CUDA_C_64I  = 25, /* complex as a pair of signed 64-bit int numbers */
    CUDA_R_64U  = 26, /* real as a unsigned 64-bit int */
    CUDA_C_64U  = 27  /* complex as a pair of unsigned 64-bit int numbers */
} cudaDataType;

typedef enum libraryPropertyType_t
{
    MAJOR_VERSION,
    MINOR_VERSION,
    PATCH_LEVEL
} libraryPropertyType;


typedef enum cudaDataType_t cudaDataType_t;
typedef enum libraryPropertyType_t libraryPropertyType_t;

/* CUBLAS status type returns */
typedef enum {
CUBLAS_STATUS_SUCCESS = 0,
CUBLAS_STATUS_NOT_INITIALIZED = 1,
CUBLAS_STATUS_ALLOC_FAILED = 3,
CUBLAS_STATUS_INVALID_VALUE = 7,
CUBLAS_STATUS_ARCH_MISMATCH = 8,
CUBLAS_STATUS_MAPPING_ERROR = 11,
CUBLAS_STATUS_EXECUTION_FAILED = 13,
CUBLAS_STATUS_INTERNAL_ERROR = 14,
CUBLAS_STATUS_NOT_SUPPORTED = 15,
CUBLAS_STATUS_LICENSE_ERROR = 16
} cublasStatus_t;

typedef enum { CUBLAS_FILL_MODE_LOWER = 0, CUBLAS_FILL_MODE_UPPER = 1, CUBLAS_FILL_MODE_FULL = 2 } cublasFillMode_t;
typedef enum { CUBLAS_DIAG_NON_UNIT = 0, CUBLAS_DIAG_UNIT = 1 } cublasDiagType_t;
typedef enum { CUBLAS_SIDE_LEFT = 0, CUBLAS_SIDE_RIGHT = 1 } cublasSideMode_t;

typedef enum {
CUBLAS_OP_N = 0,
CUBLAS_OP_T = 1,
CUBLAS_OP_C = 2,
CUBLAS_OP_HERMITAN = 2, /* synonym if CUBLAS_OP_C */
CUBLAS_OP_CONJG = 3     /* conjugate, placeholder - not supported in the current release */
} cublasOperation_t;

typedef enum { CUBLAS_POINTER_MODE_HOST = 0, CUBLAS_POINTER_MODE_DEVICE = 1 } cublasPointerMode_t;
typedef enum { CUBLAS_ATOMICS_NOT_ALLOWED = 0, CUBLAS_ATOMICS_ALLOWED = 1 } cublasAtomicsMode_t;

/*For different GEMM algorithm */
typedef enum {
CUBLAS_GEMM_DFALT = -1,
CUBLAS_GEMM_DEFAULT = -1,
CUBLAS_GEMM_ALGO0 = 0,
CUBLAS_GEMM_ALGO1 = 1,
CUBLAS_GEMM_ALGO2 = 2,
CUBLAS_GEMM_ALGO3 = 3,
CUBLAS_GEMM_ALGO4 = 4,
CUBLAS_GEMM_ALGO5 = 5,
CUBLAS_GEMM_ALGO6 = 6,
CUBLAS_GEMM_ALGO7 = 7,
CUBLAS_GEMM_ALGO8 = 8,
CUBLAS_GEMM_ALGO9 = 9,
CUBLAS_GEMM_ALGO10 = 10,
CUBLAS_GEMM_ALGO11 = 11,
CUBLAS_GEMM_ALGO12 = 12,
CUBLAS_GEMM_ALGO13 = 13,
CUBLAS_GEMM_ALGO14 = 14,
CUBLAS_GEMM_ALGO15 = 15,
CUBLAS_GEMM_ALGO16 = 16,
CUBLAS_GEMM_ALGO17 = 17,
CUBLAS_GEMM_ALGO18 = 18,  // sliced 32x32
CUBLAS_GEMM_ALGO19 = 19,  // sliced 64x32
CUBLAS_GEMM_ALGO20 = 20,  // sliced 128x32
CUBLAS_GEMM_ALGO21 = 21,  // sliced 32x32  -splitK
CUBLAS_GEMM_ALGO22 = 22,  // sliced 64x32  -splitK
CUBLAS_GEMM_ALGO23 = 23,  // sliced 128x32 -splitK
CUBLAS_GEMM_DEFAULT_TENSOR_OP = 99,
CUBLAS_GEMM_DFALT_TENSOR_OP = 99,
CUBLAS_GEMM_ALGO0_TENSOR_OP = 100,
CUBLAS_GEMM_ALGO1_TENSOR_OP = 101,
CUBLAS_GEMM_ALGO2_TENSOR_OP = 102,
CUBLAS_GEMM_ALGO3_TENSOR_OP = 103,
CUBLAS_GEMM_ALGO4_TENSOR_OP = 104,
CUBLAS_GEMM_ALGO5_TENSOR_OP = 105,
CUBLAS_GEMM_ALGO6_TENSOR_OP = 106,
CUBLAS_GEMM_ALGO7_TENSOR_OP = 107,
CUBLAS_GEMM_ALGO8_TENSOR_OP = 108,
CUBLAS_GEMM_ALGO9_TENSOR_OP = 109,
CUBLAS_GEMM_ALGO10_TENSOR_OP = 110,
CUBLAS_GEMM_ALGO11_TENSOR_OP = 111,
CUBLAS_GEMM_ALGO12_TENSOR_OP = 112,
CUBLAS_GEMM_ALGO13_TENSOR_OP = 113,
CUBLAS_GEMM_ALGO14_TENSOR_OP = 114,
CUBLAS_GEMM_ALGO15_TENSOR_OP = 115
} cublasGemmAlgo_t;

/*Enum for default math mode/tensor operation*/
typedef enum {
CUBLAS_DEFAULT_MATH = 0,

/* deprecated, same effect as using CUBLAS_COMPUTE_32F_FAST_16F, will be removed in a future release */
CUBLAS_TENSOR_OP_MATH = 1,

/* same as using matching _PEDANTIC compute type when using cublas<T>routine calls or cublasEx() calls with
    cudaDataType as compute type */
CUBLAS_PEDANTIC_MATH = 2,

/* allow accelerating single precision routines using TF32 tensor cores */
CUBLAS_TF32_TENSOR_OP_MATH = 3,

/* flag to force any reductons to use the accumulator type and not output type in case of mixed precision routines
    with lower size output type */
CUBLAS_MATH_DISALLOW_REDUCED_PRECISION_REDUCTION = 16,
} cublasMath_t;

/* For backward compatibility purposes */
//typedef cudaDataType cublasDataType_t;

/* Enum for compute type
*
* - default types provide best available performance using all available hardware features
*   and guarantee internal storage precision with at least the same precision and range;
* - _PEDANTIC types ensure standard arithmetic and exact specified internal storage format;
* - _FAST types allow for some loss of precision to enable higher throughput arithmetic.
*/
typedef enum {
CUBLAS_COMPUTE_16F = 64,           /* half - default */
CUBLAS_COMPUTE_16F_PEDANTIC = 65,  /* half - pedantic */
CUBLAS_COMPUTE_32F = 68,           /* float - default */
CUBLAS_COMPUTE_32F_PEDANTIC = 69,  /* float - pedantic */
CUBLAS_COMPUTE_32F_FAST_16F = 74,  /* float - fast, allows down-converting inputs to half or TF32 */
CUBLAS_COMPUTE_32F_FAST_16BF = 75, /* float - fast, allows down-converting inputs to bfloat16 or TF32 */
CUBLAS_COMPUTE_32F_FAST_TF32 = 77, /* float - fast, allows down-converting inputs to TF32 */
CUBLAS_COMPUTE_64F = 70,           /* double - default */
CUBLAS_COMPUTE_64F_PEDANTIC = 71,  /* double - pedantic */
CUBLAS_COMPUTE_32I = 72,           /* signed 32-bit int - default */
CUBLAS_COMPUTE_32I_PEDANTIC = 73,  /* signed 32-bit int - pedantic */
} cublasComputeType_t;

/* Opaque structure holding CUBLAS library context */
struct cublasContext;
typedef struct cublasContext* cublasHandle_t;

cublasStatus_t  cublasCreate_v2(cublasHandle_t* handle);
cublasStatus_t  cublasDestroy_v2(cublasHandle_t handle);

cublasStatus_t  cublasGetVersion_v2(cublasHandle_t handle, int* version);
cublasStatus_t  cublasGetProperty(libraryPropertyType type, int* value);
size_t          cublasGetCudartVersion(void);

cublasStatus_t  cublasSetWorkspace_v2(cublasHandle_t handle,
                                                        void* workspace,
                                                        size_t workspaceSizeInBytes);

cublasStatus_t  cublasSetStream_v2(cublasHandle_t handle, cudaStream_t streamId);
cublasStatus_t  cublasGetStream_v2(cublasHandle_t handle, cudaStream_t* streamId);

cublasStatus_t  cublasGetPointerMode_v2(cublasHandle_t handle, cublasPointerMode_t* mode);
cublasStatus_t  cublasSetPointerMode_v2(cublasHandle_t handle, cublasPointerMode_t mode);

cublasStatus_t  cublasGetAtomicsMode(cublasHandle_t handle, cublasAtomicsMode_t* mode);
cublasStatus_t  cublasSetAtomicsMode(cublasHandle_t handle, cublasAtomicsMode_t mode);

cublasStatus_t  cublasGetMathMode(cublasHandle_t handle, cublasMath_t* mode);
cublasStatus_t  cublasSetMathMode(cublasHandle_t handle, cublasMath_t mode);

cublasStatus_t  cublasGetSmCountTarget(cublasHandle_t handle, int* smCountTarget);
cublasStatus_t  cublasSetSmCountTarget(cublasHandle_t handle, int smCountTarget);

const char*  cublasGetStatusName(cublasStatus_t status);
const char*  cublasGetStatusString(cublasStatus_t status);

/* Cublas logging */
typedef void (*cublasLogCallback)(const char* msg);

cublasStatus_t  cublasLoggerConfigure(int logIsOn,
                                                        int logToStdOut,
                                                        int logToStdErr,
                                                        const char* logFileName);
cublasStatus_t  cublasSetLoggerCallback(cublasLogCallback userCallback);
cublasStatus_t  cublasGetLoggerCallback(cublasLogCallback* userCallback);
cublasStatus_t  cublasSetVector(int n, int elemSize, const void* x, int incx, void* devicePtr, int incy);
cublasStatus_t  cublasGetVector(int n, int elemSize, const void* x, int incx, void* y, int incy);
cublasStatus_t  cublasSetMatrix(int rows, int cols, int elemSize, const void* A, int lda, void* B, int ldb);
cublasStatus_t  cublasGetMatrix(int rows, int cols, int elemSize, const void* A, int lda, void* B, int ldb);
cublasStatus_t  cublasSetVectorAsync(
int n, int elemSize, const void* hostPtr, int incx, void* devicePtr, int incy, cudaStream_t stream);
cublasStatus_t  cublasGetVectorAsync(
int n, int elemSize, const void* devicePtr, int incx, void* hostPtr, int incy, cudaStream_t stream);
cublasStatus_t 
cublasSetMatrixAsync(int rows, int cols, int elemSize, const void* A, int lda, void* B, int ldb, cudaStream_t stream);

/*
* cublasStatus_t
* cublasGetMatrixAsync (int rows, int cols, int elemSize, const void *A,
*                       int lda, void *B, int ldb, cudaStream_t stream)
*
* cublasGetMatrixAsync has the same functionnality as cublasGetMatrix
* but the transfer is done asynchronously within the CUDA stream passed
* in parameter.
*
* Return Values
* -------------
* CUBLAS_STATUS_NOT_INITIALIZED  if CUBLAS library has not been initialized
* CUBLAS_STATUS_INVALID_VALUE    if rows, cols, eleSize, lda, or ldb <= 0
* CUBLAS_STATUS_MAPPING_ERROR    if error occurred accessing GPU memory
* CUBLAS_STATUS_SUCCESS          if the operation completed successfully
*/
//cublasStatus_t 
//cublasGetMatrixAsync(int rows, int cols, int elemSize, const void* A, int lda, void* B, int ldb, cudaStream_t stream);

void  cublasXerbla(const char* srName, int info);
/* ---------------- CUBLAS BLAS1 functions ---------------- */
cublasStatus_t  cublasNrm2Ex(cublasHandle_t handle,
                                                int n,
                                                const void* x,
                                                cudaDataType xType,
                                                int incx,
                                                void* result,
                                                cudaDataType resultType,
                                                cudaDataType executionType); /* host or device pointer */

cublasStatus_t 
cublasSnrm2_v2(cublasHandle_t handle, int n, const float* x, int incx, float* result); /* host or device pointer */

cublasStatus_t 
cublasDnrm2_v2(cublasHandle_t handle, int n, const double* x, int incx, double* result); /* host or device pointer */

cublasStatus_t 
cublasScnrm2_v2(cublasHandle_t handle, int n, const cuComplex* x, int incx, float* result); /* host or device pointer */

cublasStatus_t  cublasDznrm2_v2(
cublasHandle_t handle, int n, const cuDoubleComplex* x, int incx, double* result); /* host or device pointer */

cublasStatus_t  cublasDotEx(cublasHandle_t handle,
                                                int n,
                                                const void* x,
                                                cudaDataType xType,
                                                int incx,
                                                const void* y,
                                                cudaDataType yType,
                                                int incy,
                                                void* result,
                                                cudaDataType resultType,
                                                cudaDataType executionType);

cublasStatus_t  cublasDotcEx(cublasHandle_t handle,
                                                int n,
                                                const void* x,
                                                cudaDataType xType,
                                                int incx,
                                                const void* y,
                                                cudaDataType yType,
                                                int incy,
                                                void* result,
                                                cudaDataType resultType,
                                                cudaDataType executionType);

cublasStatus_t  cublasSdot_v2(cublasHandle_t handle,
                                                int n,
                                                const float* x,
                                                int incx,
                                                const float* y,
                                                int incy,
                                                float* result); /* host or device pointer */

cublasStatus_t  cublasDdot_v2(cublasHandle_t handle,
                                                int n,
                                                const double* x,
                                                int incx,
                                                const double* y,
                                                int incy,
                                                double* result); /* host or device pointer */

                                                    
cublasStatus_t  cublasCdotu_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* result); /* host or device pointer */

cublasStatus_t  cublasCdotc_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* result); /* host or device pointer */
                                                    
cublasStatus_t  cublasZdotu_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* result); /* host or device pointer */

                                                    
cublasStatus_t  cublasZdotc_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* result); /* host or device pointer */

                                                    
cublasStatus_t  cublasScalEx(cublasHandle_t handle,
                                                int n,
                                                const void* alpha, /* host or device pointer */
                                                cudaDataType alphaType,
                                                void* x,
                                                cudaDataType xType,
                                                int incx,
                                                cudaDataType executionType);

cublasStatus_t  cublasSscal_v2(cublasHandle_t handle,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDscal_v2(cublasHandle_t handle,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCscal_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasCsscal_v2(cublasHandle_t handle,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZscal_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    cuDoubleComplex* x,
                                                    int incx);

cublasStatus_t  cublasZdscal_v2(cublasHandle_t handle,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    cuDoubleComplex* x,
                                                    int incx);

cublasStatus_t  cublasAxpyEx(cublasHandle_t handle,
                                                int n,
                                                const void* alpha, /* host or device pointer */
                                                cudaDataType alphaType,
                                                const void* x,
                                                cudaDataType xType,
                                                int incx,
                                                void* y,
                                                cudaDataType yType,
                                                int incy,
                                                cudaDataType executiontype);

cublasStatus_t  cublasSaxpy_v2(cublasHandle_t handle,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* x,
                                                    int incx,
                                                    float* y,
                                                    int incy);

cublasStatus_t  cublasDaxpy_v2(cublasHandle_t handle,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* x,
                                                    int incx,
                                                    double* y,
                                                    int incy);

cublasStatus_t  cublasCaxpy_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* x,
                                                    int incx,
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZaxpy_v2(cublasHandle_t handle,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    cuDoubleComplex* y,
                                                    int incy);

cublasStatus_t  cublasCopyEx(
cublasHandle_t handle, int n, const void* x, cudaDataType xType, int incx, void* y, cudaDataType yType, int incy);

cublasStatus_t 
cublasScopy_v2(cublasHandle_t handle, int n, const float* x, int incx, float* y, int incy);

cublasStatus_t 
cublasDcopy_v2(cublasHandle_t handle, int n, const double* x, int incx, double* y, int incy);

cublasStatus_t 
cublasCcopy_v2(cublasHandle_t handle, int n, const cuComplex* x, int incx, cuComplex* y, int incy);

cublasStatus_t 
cublasZcopy_v2(cublasHandle_t handle, int n, const cuDoubleComplex* x, int incx, cuDoubleComplex* y, int incy);

cublasStatus_t 
cublasSswap_v2(cublasHandle_t handle, int n, float* x, int incx, float* y, int incy);

cublasStatus_t 
cublasDswap_v2(cublasHandle_t handle, int n, double* x, int incx, double* y, int incy);

cublasStatus_t 
cublasCswap_v2(cublasHandle_t handle, int n, cuComplex* x, int incx, cuComplex* y, int incy);

cublasStatus_t 
cublasZswap_v2(cublasHandle_t handle, int n, cuDoubleComplex* x, int incx, cuDoubleComplex* y, int incy);

cublasStatus_t  cublasSwapEx(
cublasHandle_t handle, int n, void* x, cudaDataType xType, int incx, void* y, cudaDataType yType, int incy);

cublasStatus_t 
cublasIsamax_v2(cublasHandle_t handle, int n, const float* x, int incx, int* result); /* host or device pointer */

cublasStatus_t 
cublasIdamax_v2(cublasHandle_t handle, int n, const double* x, int incx, int* result); /* host or device pointer */

cublasStatus_t 
cublasIcamax_v2(cublasHandle_t handle, int n, const cuComplex* x, int incx, int* result); /* host or device pointer */

cublasStatus_t  cublasIzamax_v2(
cublasHandle_t handle, int n, const cuDoubleComplex* x, int incx, int* result); /* host or device pointer */

cublasStatus_t  cublasIamaxEx(
cublasHandle_t handle, int n, const void* x, cudaDataType xType, int incx, int* result /* host or device pointer */
);

cublasStatus_t 
cublasIsamin_v2(cublasHandle_t handle, int n, const float* x, int incx, int* result); /* host or device pointer */

cublasStatus_t 
cublasIdamin_v2(cublasHandle_t handle, int n, const double* x, int incx, int* result); /* host or device pointer */

cublasStatus_t 
cublasIcamin_v2(cublasHandle_t handle, int n, const cuComplex* x, int incx, int* result); /* host or device pointer */

cublasStatus_t  cublasIzamin_v2(
cublasHandle_t handle, int n, const cuDoubleComplex* x, int incx, int* result); /* host or device pointer */

cublasStatus_t  cublasIaminEx(
cublasHandle_t handle, int n, const void* x, cudaDataType xType, int incx, int* result /* host or device pointer */
);

cublasStatus_t  cublasAsumEx(cublasHandle_t handle,
                                                int n,
                                                const void* x,
                                                cudaDataType xType,
                                                int incx,
                                                void* result,
                                                cudaDataType resultType, /* host or device pointer */
                                                cudaDataType executiontype);

cublasStatus_t 
cublasSasum_v2(cublasHandle_t handle, int n, const float* x, int incx, float* result); /* host or device pointer */

cublasStatus_t 
cublasDasum_v2(cublasHandle_t handle, int n, const double* x, int incx, double* result); /* host or device pointer */

cublasStatus_t 
cublasScasum_v2(cublasHandle_t handle, int n, const cuComplex* x, int incx, float* result); /* host or device pointer */

cublasStatus_t  cublasDzasum_v2(
cublasHandle_t handle, int n, const cuDoubleComplex* x, int incx, double* result); /* host or device pointer */

cublasStatus_t  cublasSrot_v2(cublasHandle_t handle,
                                                int n,
                                                float* x,
                                                int incx,
                                                float* y,
                                                int incy,
                                                const float* c,  /* host or device pointer */
                                                const float* s); /* host or device pointer */

cublasStatus_t  cublasDrot_v2(cublasHandle_t handle,
                                                int n,
                                                double* x,
                                                int incx,
                                                double* y,
                                                int incy,
                                                const double* c,  /* host or device pointer */
                                                const double* s); /* host or device pointer */

cublasStatus_t  cublasCrot_v2(cublasHandle_t handle,
                                                int n,
                                                cuComplex* x,
                                                int incx,
                                                cuComplex* y,
                                                int incy,
                                                const float* c,      /* host or device pointer */
                                                const cuComplex* s); /* host or device pointer */

cublasStatus_t  cublasCsrot_v2(cublasHandle_t handle,
                                                    int n,
                                                    cuComplex* x,
                                                    int incx,
                                                    cuComplex* y,
                                                    int incy,
                                                    const float* c,  /* host or device pointer */
                                                    const float* s); /* host or device pointer */

cublasStatus_t  cublasZrot_v2(cublasHandle_t handle,
                                                int n,
                                                cuDoubleComplex* x,
                                                int incx,
                                                cuDoubleComplex* y,
                                                int incy,
                                                const double* c,           /* host or device pointer */
                                                const cuDoubleComplex* s); /* host or device pointer */

cublasStatus_t  cublasZdrot_v2(cublasHandle_t handle,
                                                    int n,
                                                    cuDoubleComplex* x,
                                                    int incx,
                                                    cuDoubleComplex* y,
                                                    int incy,
                                                    const double* c,  /* host or device pointer */
                                                    const double* s); /* host or device pointer */

cublasStatus_t  cublasRotEx(cublasHandle_t handle,
                                                int n,
                                                void* x,
                                                cudaDataType xType,
                                                int incx,
                                                void* y,
                                                cudaDataType yType,
                                                int incy,
                                                const void* c, /* host or device pointer */
                                                const void* s,
                                                cudaDataType csType,
                                                cudaDataType executiontype);

cublasStatus_t  cublasSrotg_v2(cublasHandle_t handle,
                                                    float* a,  /* host or device pointer */
                                                    float* b,  /* host or device pointer */
                                                    float* c,  /* host or device pointer */
                                                    float* s); /* host or device pointer */

cublasStatus_t  cublasDrotg_v2(cublasHandle_t handle,
                                                    double* a,  /* host or device pointer */
                                                    double* b,  /* host or device pointer */
                                                    double* c,  /* host or device pointer */
                                                    double* s); /* host or device pointer */

cublasStatus_t  cublasCrotg_v2(cublasHandle_t handle,
                                                    cuComplex* a,  /* host or device pointer */
                                                    cuComplex* b,  /* host or device pointer */
                                                    float* c,      /* host or device pointer */
                                                    cuComplex* s); /* host or device pointer */

cublasStatus_t  cublasZrotg_v2(cublasHandle_t handle,
                                                    cuDoubleComplex* a,  /* host or device pointer */
                                                    cuDoubleComplex* b,  /* host or device pointer */
                                                    double* c,           /* host or device pointer */
                                                    cuDoubleComplex* s); /* host or device pointer */

cublasStatus_t  cublasRotgEx(cublasHandle_t handle,
                                                void* a, /* host or device pointer */
                                                void* b, /* host or device pointer */
                                                cudaDataType abType,
                                                void* c, /* host or device pointer */
                                                void* s, /* host or device pointer */
                                                cudaDataType csType,
                                                cudaDataType executiontype);

cublasStatus_t  cublasSrotm_v2(cublasHandle_t handle,
                                                    int n,
                                                    float* x,
                                                    int incx,
                                                    float* y,
                                                    int incy,
                                                    const float* param); /* host or device pointer */

cublasStatus_t  cublasDrotm_v2(cublasHandle_t handle,
                                                    int n,
                                                    double* x,
                                                    int incx,
                                                    double* y,
                                                    int incy,
                                                    const double* param); /* host or device pointer */

cublasStatus_t  cublasRotmEx(cublasHandle_t handle,
                                                int n,
                                                void* x,
                                                cudaDataType xType,
                                                int incx,
                                                void* y,
                                                cudaDataType yType,
                                                int incy,
                                                const void* param, /* host or device pointer */
                                                cudaDataType paramType,
                                                cudaDataType executiontype);

cublasStatus_t  cublasSrotmg_v2(cublasHandle_t handle,
                                                    float* d1,       /* host or device pointer */
                                                    float* d2,       /* host or device pointer */
                                                    float* x1,       /* host or device pointer */
                                                    const float* y1, /* host or device pointer */
                                                    float* param);   /* host or device pointer */

cublasStatus_t  cublasDrotmg_v2(cublasHandle_t handle,
                                                    double* d1,       /* host or device pointer */
                                                    double* d2,       /* host or device pointer */
                                                    double* x1,       /* host or device pointer */
                                                    const double* y1, /* host or device pointer */
                                                    double* param);   /* host or device pointer */

cublasStatus_t  cublasRotmgEx(cublasHandle_t handle,
                                                void* d1, /* host or device pointer */
                                                cudaDataType d1Type,
                                                void* d2, /* host or device pointer */
                                                cudaDataType d2Type,
                                                void* x1, /* host or device pointer */
                                                cudaDataType x1Type,
                                                const void* y1, /* host or device pointer */
                                                cudaDataType y1Type,
                                                void* param, /* host or device pointer */
                                                cudaDataType paramType,
                                                cudaDataType executiontype);
/* --------------- CUBLAS BLAS2 functions  ---------------- */

/* GEMV */
cublasStatus_t  cublasSgemv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* x,
                                                    int incx,
                                                    const float* beta, /* host or device pointer */
                                                    float* y,
                                                    int incy);

cublasStatus_t  cublasDgemv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* x,
                                                    int incx,
                                                    const double* beta, /* host or device pointer */
                                                    double* y,
                                                    int incy);

cublasStatus_t  cublasCgemv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZgemv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* y,
                                                    int incy);
/* GBMV */
cublasStatus_t  cublasSgbmv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    int kl,
                                                    int ku,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* x,
                                                    int incx,
                                                    const float* beta, /* host or device pointer */
                                                    float* y,
                                                    int incy);

cublasStatus_t  cublasDgbmv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    int kl,
                                                    int ku,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* x,
                                                    int incx,
                                                    const double* beta, /* host or device pointer */
                                                    double* y,
                                                    int incy);

cublasStatus_t  cublasCgbmv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    int kl,
                                                    int ku,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZgbmv_v2(cublasHandle_t handle,
                                                    cublasOperation_t trans,
                                                    int m,
                                                    int n,
                                                    int kl,
                                                    int ku,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* y,
                                                    int incy);

/* TRMV */
cublasStatus_t  cublasStrmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const float* A,
                                                    int lda,
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDtrmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const double* A,
                                                    int lda,
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCtrmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuComplex* A,
                                                    int lda,
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZtrmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    cuDoubleComplex* x,
                                                    int incx);

/* TBMV */
cublasStatus_t  cublasStbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const float* A,
                                                    int lda,
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDtbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const double* A,
                                                    int lda,
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCtbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const cuComplex* A,
                                                    int lda,
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZtbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    cuDoubleComplex* x,
                                                    int incx);

/* TPMV */
cublasStatus_t  cublasStpmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const float* AP,
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDtpmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const double* AP,
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCtpmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuComplex* AP,
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZtpmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuDoubleComplex* AP,
                                                    cuDoubleComplex* x,
                                                    int incx);

/* TRSV */
cublasStatus_t  cublasStrsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const float* A,
                                                    int lda,
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDtrsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const double* A,
                                                    int lda,
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCtrsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuComplex* A,
                                                    int lda,
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZtrsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    cuDoubleComplex* x,
                                                    int incx);

/* TPSV */
cublasStatus_t  cublasStpsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const float* AP,
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDtpsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const double* AP,
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCtpsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuComplex* AP,
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZtpsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    const cuDoubleComplex* AP,
                                                    cuDoubleComplex* x,
                                                    int incx);
/* TBSV */
cublasStatus_t  cublasStbsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const float* A,
                                                    int lda,
                                                    float* x,
                                                    int incx);

cublasStatus_t  cublasDtbsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const double* A,
                                                    int lda,
                                                    double* x,
                                                    int incx);

cublasStatus_t  cublasCtbsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const cuComplex* A,
                                                    int lda,
                                                    cuComplex* x,
                                                    int incx);

cublasStatus_t  cublasZtbsv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    cuDoubleComplex* x,
                                                    int incx);

/* SYMV/HEMV */
cublasStatus_t  cublasSsymv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* x,
                                                    int incx,
                                                    const float* beta, /* host or device pointer */
                                                    float* y,
                                                    int incy);

cublasStatus_t  cublasDsymv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* x,
                                                    int incx,
                                                    const double* beta, /* host or device pointer */
                                                    double* y,
                                                    int incy);

cublasStatus_t  cublasCsymv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZsymv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* y,
                                                    int incy);

cublasStatus_t  cublasChemv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZhemv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* y,
                                                    int incy);

/* SBMV/HBMV */
cublasStatus_t  cublasSsbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    int k,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* x,
                                                    int incx,
                                                    const float* beta, /* host or device pointer */
                                                    float* y,
                                                    int incy);

cublasStatus_t  cublasDsbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    int k,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* x,
                                                    int incx,
                                                    const double* beta, /* host or device pointer */
                                                    double* y,
                                                    int incy);

cublasStatus_t  cublasChbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZhbmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* y,
                                                    int incy);

/* SPMV/HPMV */
cublasStatus_t  cublasSspmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* AP,
                                                    const float* x,
                                                    int incx,
                                                    const float* beta, /* host or device pointer */
                                                    float* y,
                                                    int incy);

cublasStatus_t  cublasDspmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* AP,
                                                    const double* x,
                                                    int incx,
                                                    const double* beta, /* host or device pointer */
                                                    double* y,
                                                    int incy);

cublasStatus_t  cublasChpmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* AP,
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* y,
                                                    int incy);

cublasStatus_t  cublasZhpmv_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* AP,
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* y,
                                                    int incy);

/* GER */
cublasStatus_t  cublasSger_v2(cublasHandle_t handle,
                                                int m,
                                                int n,
                                                const float* alpha, /* host or device pointer */
                                                const float* x,
                                                int incx,
                                                const float* y,
                                                int incy,
                                                float* A,
                                                int lda);

cublasStatus_t  cublasDger_v2(cublasHandle_t handle,
                                                int m,
                                                int n,
                                                const double* alpha, /* host or device pointer */
                                                const double* x,
                                                int incx,
                                                const double* y,
                                                int incy,
                                                double* A,
                                                int lda);

cublasStatus_t  cublasCgeru_v2(cublasHandle_t handle,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* A,
                                                    int lda);

cublasStatus_t  cublasCgerc_v2(cublasHandle_t handle,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* A,
                                                    int lda);

cublasStatus_t  cublasZgeru_v2(cublasHandle_t handle,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* A,
                                                    int lda);

cublasStatus_t  cublasZgerc_v2(cublasHandle_t handle,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* A,
                                                    int lda);

/* SYR/HER */
cublasStatus_t  cublasSsyr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const float* alpha, /* host or device pointer */
                                                const float* x,
                                                int incx,
                                                float* A,
                                                int lda);

cublasStatus_t  cublasDsyr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const double* alpha, /* host or device pointer */
                                                const double* x,
                                                int incx,
                                                double* A,
                                                int lda);

cublasStatus_t  cublasCsyr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const cuComplex* alpha, /* host or device pointer */
                                                const cuComplex* x,
                                                int incx,
                                                cuComplex* A,
                                                int lda);

cublasStatus_t  cublasZsyr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const cuDoubleComplex* alpha, /* host or device pointer */
                                                const cuDoubleComplex* x,
                                                int incx,
                                                cuDoubleComplex* A,
                                                int lda);

cublasStatus_t  cublasCher_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const float* alpha, /* host or device pointer */
                                                const cuComplex* x,
                                                int incx,
                                                cuComplex* A,
                                                int lda);

cublasStatus_t  cublasZher_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const double* alpha, /* host or device pointer */
                                                const cuDoubleComplex* x,
                                                int incx,
                                                cuDoubleComplex* A,
                                                int lda);

/* SPR/HPR */
cublasStatus_t  cublasSspr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const float* alpha, /* host or device pointer */
                                                const float* x,
                                                int incx,
                                                float* AP);

cublasStatus_t  cublasDspr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const double* alpha, /* host or device pointer */
                                                const double* x,
                                                int incx,
                                                double* AP);

cublasStatus_t  cublasChpr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const float* alpha, /* host or device pointer */
                                                const cuComplex* x,
                                                int incx,
                                                cuComplex* AP);

cublasStatus_t  cublasZhpr_v2(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                int n,
                                                const double* alpha, /* host or device pointer */
                                                const cuDoubleComplex* x,
                                                int incx,
                                                cuDoubleComplex* AP);

/* SYR2/HER2 */
cublasStatus_t  cublasSsyr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* x,
                                                    int incx,
                                                    const float* y,
                                                    int incy,
                                                    float* A,
                                                    int lda);

cublasStatus_t  cublasDsyr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* x,
                                                    int incx,
                                                    const double* y,
                                                    int incy,
                                                    double* A,
                                                    int lda);

cublasStatus_t  cublasCsyr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* A,
                                                    int lda);

cublasStatus_t  cublasZsyr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* A,
                                                    int lda);

cublasStatus_t  cublasCher2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* A,
                                                    int lda);

cublasStatus_t  cublasZher2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* A,
                                                    int lda);

/* SPR2/HPR2 */
cublasStatus_t  cublasSspr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* x,
                                                    int incx,
                                                    const float* y,
                                                    int incy,
                                                    float* AP);

cublasStatus_t  cublasDspr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* x,
                                                    int incx,
                                                    const double* y,
                                                    int incy,
                                                    double* AP);

cublasStatus_t  cublasChpr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* x,
                                                    int incx,
                                                    const cuComplex* y,
                                                    int incy,
                                                    cuComplex* AP);

cublasStatus_t  cublasZhpr2_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* x,
                                                    int incx,
                                                    const cuDoubleComplex* y,
                                                    int incy,
                                                    cuDoubleComplex* AP);

/* ---------------- CUBLAS BLAS3 functions ---------------- */

/* GEMM */
cublasStatus_t  cublasSgemm_v2(cublasHandle_t handle,
                                                    cublasOperation_t transa,
                                                    cublasOperation_t transb,
                                                    int m,
                                                    int n,
                                                    int k,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* B,
                                                    int ldb,
                                                    const float* beta, /* host or device pointer */
                                                    float* C,
                                                    int ldc);

cublasStatus_t  cublasDgemm_v2(cublasHandle_t handle,
                                                    cublasOperation_t transa,
                                                    cublasOperation_t transb,
                                                    int m,
                                                    int n,
                                                    int k,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* B,
                                                    int ldb,
                                                    const double* beta, /* host or device pointer */
                                                    double* C,
                                                    int ldc);

cublasStatus_t  cublasCgemm_v2(cublasHandle_t handle,
                                                    cublasOperation_t transa,
                                                    cublasOperation_t transb,
                                                    int m,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* B,
                                                    int ldb,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasCgemm3m(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                int k,
                                                const cuComplex* alpha, /* host or device pointer */
                                                const cuComplex* A,
                                                int lda,
                                                const cuComplex* B,
                                                int ldb,
                                                const cuComplex* beta, /* host or device pointer */
                                                cuComplex* C,
                                                int ldc);
cublasStatus_t  cublasCgemm3mEx(cublasHandle_t handle,
                                                    cublasOperation_t transa,
                                                    cublasOperation_t transb,
                                                    int m,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha,
                                                    const void* A,
                                                    cudaDataType Atype,
                                                    int lda,
                                                    const void* B,
                                                    cudaDataType Btype,
                                                    int ldb,
                                                    const cuComplex* beta,
                                                    void* C,
                                                    cudaDataType Ctype,
                                                    int ldc);

cublasStatus_t  cublasZgemm_v2(cublasHandle_t handle,
                                                    cublasOperation_t transa,
                                                    cublasOperation_t transb,
                                                    int m,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* B,
                                                    int ldb,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZgemm3m(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                int k,
                                                const cuDoubleComplex* alpha, /* host or device pointer */
                                                const cuDoubleComplex* A,
                                                int lda,
                                                const cuDoubleComplex* B,
                                                int ldb,
                                                const cuDoubleComplex* beta, /* host or device pointer */
                                                cuDoubleComplex* C,
                                                int ldc);

#if defined(__cplusplus)
cublasStatus_t  cublasHgemm(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                int k,
                                                const __half* alpha, /* host or device pointer */
                                                const __half* A,
                                                int lda,
                                                const __half* B,
                                                int ldb,
                                                const __half* beta, /* host or device pointer */
                                                __half* C,
                                                int ldc);
#endif
/* IO in FP16/FP32, computation in float */
cublasStatus_t  cublasSgemmEx(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                int k,
                                                const float* alpha, /* host or device pointer */
                                                const void* A,
                                                cudaDataType Atype,
                                                int lda,
                                                const void* B,
                                                cudaDataType Btype,
                                                int ldb,
                                                const float* beta, /* host or device pointer */
                                                void* C,
                                                cudaDataType Ctype,
                                                int ldc);

cublasStatus_t  cublasGemmEx(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                int k,
                                                const void* alpha, /* host or device pointer */
                                                const void* A,
                                                cudaDataType Atype,
                                                int lda,
                                                const void* B,
                                                cudaDataType Btype,
                                                int ldb,
                                                const void* beta, /* host or device pointer */
                                                void* C,
                                                cudaDataType Ctype,
                                                int ldc,
                                                cublasComputeType_t computeType,
                                                cublasGemmAlgo_t algo);

/* IO in Int8 complex/cuComplex, computation in cuComplex */
cublasStatus_t  cublasCgemmEx(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                int k,
                                                const cuComplex* alpha,
                                                const void* A,
                                                cudaDataType Atype,
                                                int lda,
                                                const void* B,
                                                cudaDataType Btype,
                                                int ldb,
                                                const cuComplex* beta,
                                                void* C,
                                                cudaDataType Ctype,
                                                int ldc);

cublasStatus_t  cublasUint8gemmBias(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        cublasOperation_t transc,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const unsigned char* A,
                                                        int A_bias,
                                                        int lda,
                                                        const unsigned char* B,
                                                        int B_bias,
                                                        int ldb,
                                                        unsigned char* C,
                                                        int C_bias,
                                                        int ldc,
                                                        int C_mult,
                                                        int C_shift);

/* SYRK */
cublasStatus_t  cublasSsyrk_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* beta, /* host or device pointer */
                                                    float* C,
                                                    int ldc);

cublasStatus_t  cublasDsyrk_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* beta, /* host or device pointer */
                                                    double* C,
                                                    int ldc);

cublasStatus_t  cublasCsyrk_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZsyrk_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);
/* IO in Int8 complex/cuComplex, computation in cuComplex */
cublasStatus_t  cublasCsyrkEx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const cuComplex* alpha, /* host or device pointer */
                                                const void* A,
                                                cudaDataType Atype,
                                                int lda,
                                                const cuComplex* beta, /* host or device pointer */
                                                void* C,
                                                cudaDataType Ctype,
                                                int ldc);

/* IO in Int8 complex/cuComplex, computation in cuComplex, Gaussian math */
cublasStatus_t  cublasCsyrk3mEx(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha,
                                                    const void* A,
                                                    cudaDataType Atype,
                                                    int lda,
                                                    const cuComplex* beta,
                                                    void* C,
                                                    cudaDataType Ctype,
                                                    int ldc);

/* HERK */
cublasStatus_t  cublasCherk_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const float* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const float* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZherk_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const double* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const double* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);

/* IO in Int8 complex/cuComplex, computation in cuComplex */
cublasStatus_t  cublasCherkEx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const float* alpha, /* host or device pointer */
                                                const void* A,
                                                cudaDataType Atype,
                                                int lda,
                                                const float* beta, /* host or device pointer */
                                                void* C,
                                                cudaDataType Ctype,
                                                int ldc);

/* IO in Int8 complex/cuComplex, computation in cuComplex, Gaussian math */
cublasStatus_t  cublasCherk3mEx(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const float* alpha,
                                                    const void* A,
                                                    cudaDataType Atype,
                                                    int lda,
                                                    const float* beta,
                                                    void* C,
                                                    cudaDataType Ctype,
                                                    int ldc);

/* SYR2K */
cublasStatus_t  cublasSsyr2k_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* B,
                                                    int ldb,
                                                    const float* beta, /* host or device pointer */
                                                    float* C,
                                                    int ldc);

cublasStatus_t  cublasDsyr2k_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* B,
                                                    int ldb,
                                                    const double* beta, /* host or device pointer */
                                                    double* C,
                                                    int ldc);

cublasStatus_t  cublasCsyr2k_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* B,
                                                    int ldb,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZsyr2k_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* B,
                                                    int ldb,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);
/* HER2K */
cublasStatus_t  cublasCher2k_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* B,
                                                    int ldb,
                                                    const float* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZher2k_v2(cublasHandle_t handle,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    int n,
                                                    int k,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* B,
                                                    int ldb,
                                                    const double* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);
/* SYRKX : eXtended SYRK*/
cublasStatus_t  cublasSsyrkx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const float* alpha, /* host or device pointer */
                                                const float* A,
                                                int lda,
                                                const float* B,
                                                int ldb,
                                                const float* beta, /* host or device pointer */
                                                float* C,
                                                int ldc);

cublasStatus_t  cublasDsyrkx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const double* alpha, /* host or device pointer */
                                                const double* A,
                                                int lda,
                                                const double* B,
                                                int ldb,
                                                const double* beta, /* host or device pointer */
                                                double* C,
                                                int ldc);

cublasStatus_t  cublasCsyrkx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const cuComplex* alpha, /* host or device pointer */
                                                const cuComplex* A,
                                                int lda,
                                                const cuComplex* B,
                                                int ldb,
                                                const cuComplex* beta, /* host or device pointer */
                                                cuComplex* C,
                                                int ldc);

cublasStatus_t  cublasZsyrkx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const cuDoubleComplex* alpha, /* host or device pointer */
                                                const cuDoubleComplex* A,
                                                int lda,
                                                const cuDoubleComplex* B,
                                                int ldb,
                                                const cuDoubleComplex* beta, /* host or device pointer */
                                                cuDoubleComplex* C,
                                                int ldc);
/* HERKX : eXtended HERK */
cublasStatus_t  cublasCherkx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const cuComplex* alpha, /* host or device pointer */
                                                const cuComplex* A,
                                                int lda,
                                                const cuComplex* B,
                                                int ldb,
                                                const float* beta, /* host or device pointer */
                                                cuComplex* C,
                                                int ldc);

cublasStatus_t  cublasZherkx(cublasHandle_t handle,
                                                cublasFillMode_t uplo,
                                                cublasOperation_t trans,
                                                int n,
                                                int k,
                                                const cuDoubleComplex* alpha, /* host or device pointer */
                                                const cuDoubleComplex* A,
                                                int lda,
                                                const cuDoubleComplex* B,
                                                int ldb,
                                                const double* beta, /* host or device pointer */
                                                cuDoubleComplex* C,
                                                int ldc);
/* SYMM */
cublasStatus_t  cublasSsymm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    int m,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* B,
                                                    int ldb,
                                                    const float* beta, /* host or device pointer */
                                                    float* C,
                                                    int ldc);

cublasStatus_t  cublasDsymm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    int m,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* B,
                                                    int ldb,
                                                    const double* beta, /* host or device pointer */
                                                    double* C,
                                                    int ldc);

cublasStatus_t  cublasCsymm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* B,
                                                    int ldb,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZsymm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* B,
                                                    int ldb,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);

/* HEMM */
cublasStatus_t  cublasChemm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* B,
                                                    int ldb,
                                                    const cuComplex* beta, /* host or device pointer */
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZhemm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* B,
                                                    int ldb,
                                                    const cuDoubleComplex* beta, /* host or device pointer */
                                                    cuDoubleComplex* C,
                                                    int ldc);

/* TRSM */
cublasStatus_t  cublasStrsm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    float* B,
                                                    int ldb);

cublasStatus_t  cublasDtrsm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    double* B,
                                                    int ldb);

cublasStatus_t  cublasCtrsm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    cuComplex* B,
                                                    int ldb);

cublasStatus_t  cublasZtrsm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    cuDoubleComplex* B,
                                                    int ldb);

/* TRMM */
cublasStatus_t  cublasStrmm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const float* alpha, /* host or device pointer */
                                                    const float* A,
                                                    int lda,
                                                    const float* B,
                                                    int ldb,
                                                    float* C,
                                                    int ldc);

cublasStatus_t  cublasDtrmm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const double* alpha, /* host or device pointer */
                                                    const double* A,
                                                    int lda,
                                                    const double* B,
                                                    int ldb,
                                                    double* C,
                                                    int ldc);

cublasStatus_t  cublasCtrmm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const cuComplex* alpha, /* host or device pointer */
                                                    const cuComplex* A,
                                                    int lda,
                                                    const cuComplex* B,
                                                    int ldb,
                                                    cuComplex* C,
                                                    int ldc);

cublasStatus_t  cublasZtrmm_v2(cublasHandle_t handle,
                                                    cublasSideMode_t side,
                                                    cublasFillMode_t uplo,
                                                    cublasOperation_t trans,
                                                    cublasDiagType_t diag,
                                                    int m,
                                                    int n,
                                                    const cuDoubleComplex* alpha, /* host or device pointer */
                                                    const cuDoubleComplex* A,
                                                    int lda,
                                                    const cuDoubleComplex* B,
                                                    int ldb,
                                                    cuDoubleComplex* C,
                                                    int ldc);
/* BATCH GEMM */

cublasStatus_t  cublasHgemmBatched(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const __half* alpha, /* host or device pointer */
                                                        const __half* const Aarray[],
                                                        int lda,
                                                        const __half* const Barray[],
                                                        int ldb,
                                                        const __half* beta, /* host or device pointer */
                                                        __half* const Carray[],
                                                        int ldc,
                                                        int batchCount);

cublasStatus_t  cublasSgemmBatched(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const float* alpha, /* host or device pointer */
                                                        const float* const Aarray[],
                                                        int lda,
                                                        const float* const Barray[],
                                                        int ldb,
                                                        const float* beta, /* host or device pointer */
                                                        float* const Carray[],
                                                        int ldc,
                                                        int batchCount);

cublasStatus_t  cublasDgemmBatched(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const double* alpha, /* host or device pointer */
                                                        const double* const Aarray[],
                                                        int lda,
                                                        const double* const Barray[],
                                                        int ldb,
                                                        const double* beta, /* host or device pointer */
                                                        double* const Carray[],
                                                        int ldc,
                                                        int batchCount);

cublasStatus_t  cublasCgemmBatched(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const cuComplex* alpha, /* host or device pointer */
                                                        const cuComplex* const Aarray[],
                                                        int lda,
                                                        const cuComplex* const Barray[],
                                                        int ldb,
                                                        const cuComplex* beta, /* host or device pointer */
                                                        cuComplex* const Carray[],
                                                        int ldc,
                                                        int batchCount);

cublasStatus_t  cublasCgemm3mBatched(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const cuComplex* alpha, /* host or device pointer */
                                                        const cuComplex* const Aarray[],
                                                        int lda,
                                                        const cuComplex* const Barray[],
                                                        int ldb,
                                                        const cuComplex* beta, /* host or device pointer */
                                                        cuComplex* const Carray[],
                                                        int ldc,
                                                        int batchCount);

cublasStatus_t  cublasZgemmBatched(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const cuDoubleComplex* alpha, /* host or device pointer */
                                                        const cuDoubleComplex* const Aarray[],
                                                        int lda,
                                                        const cuDoubleComplex* const Barray[],
                                                        int ldb,
                                                        const cuDoubleComplex* beta, /* host or device pointer */
                                                        cuDoubleComplex* const Carray[],
                                                        int ldc,
                                                        int batchCount);

cublasStatus_t  cublasGemmBatchedEx(cublasHandle_t handle,
                                                        cublasOperation_t transa,
                                                        cublasOperation_t transb,
                                                        int m,
                                                        int n,
                                                        int k,
                                                        const void* alpha, /* host or device pointer */
                                                        const void* const Aarray[],
                                                        cudaDataType Atype,
                                                        int lda,
                                                        const void* const Barray[],
                                                        cudaDataType Btype,
                                                        int ldb,
                                                        const void* beta, /* host or device pointer */
                                                        void* const Carray[],
                                                        cudaDataType Ctype,
                                                        int ldc,
                                                        int batchCount,
                                                        cublasComputeType_t computeType,
                                                        cublasGemmAlgo_t algo);

cublasStatus_t  cublasGemmStridedBatchedEx(cublasHandle_t handle,
                                                                cublasOperation_t transa,
                                                                cublasOperation_t transb,
                                                                int m,
                                                                int n,
                                                                int k,
                                                                const void* alpha, /* host or device pointer */
                                                                const void* A,
                                                                cudaDataType Atype,
                                                                int lda,
                                                                long long int strideA, /* purposely signed */
                                                                const void* B,
                                                                cudaDataType Btype,
                                                                int ldb,
                                                                long long int strideB,
                                                                const void* beta, /* host or device pointer */
                                                                void* C,
                                                                cudaDataType Ctype,
                                                                int ldc,
                                                                long long int strideC,
                                                                int batchCount,
                                                                cublasComputeType_t computeType,
                                                                cublasGemmAlgo_t algo);

cublasStatus_t  cublasSgemmStridedBatched(cublasHandle_t handle,
                                                            cublasOperation_t transa,
                                                            cublasOperation_t transb,
                                                            int m,
                                                            int n,
                                                            int k,
                                                            const float* alpha, /* host or device pointer */
                                                            const float* A,
                                                            int lda,
                                                            long long int strideA, /* purposely signed */
                                                            const float* B,
                                                            int ldb,
                                                            long long int strideB,
                                                            const float* beta, /* host or device pointer */
                                                            float* C,
                                                            int ldc,
                                                            long long int strideC,
                                                            int batchCount);

cublasStatus_t  cublasDgemmStridedBatched(cublasHandle_t handle,
                                                            cublasOperation_t transa,
                                                            cublasOperation_t transb,
                                                            int m,
                                                            int n,
                                                            int k,
                                                            const double* alpha, /* host or device pointer */
                                                            const double* A,
                                                            int lda,
                                                            long long int strideA, /* purposely signed */
                                                            const double* B,
                                                            int ldb,
                                                            long long int strideB,
                                                            const double* beta, /* host or device pointer */
                                                            double* C,
                                                            int ldc,
                                                            long long int strideC,
                                                            int batchCount);

cublasStatus_t  cublasCgemmStridedBatched(cublasHandle_t handle,
                                                            cublasOperation_t transa,
                                                            cublasOperation_t transb,
                                                            int m,
                                                            int n,
                                                            int k,
                                                            const cuComplex* alpha, /* host or device pointer */
                                                            const cuComplex* A,
                                                            int lda,
                                                            long long int strideA, /* purposely signed */
                                                            const cuComplex* B,
                                                            int ldb,
                                                            long long int strideB,
                                                            const cuComplex* beta, /* host or device pointer */
                                                            cuComplex* C,
                                                            int ldc,
                                                            long long int strideC,
                                                            int batchCount);

cublasStatus_t  cublasCgemm3mStridedBatched(cublasHandle_t handle,
                                                                cublasOperation_t transa,
                                                                cublasOperation_t transb,
                                                                int m,
                                                                int n,
                                                                int k,
                                                                const cuComplex* alpha, /* host or device pointer */
                                                                const cuComplex* A,
                                                                int lda,
                                                                long long int strideA, /* purposely signed */
                                                                const cuComplex* B,
                                                                int ldb,
                                                                long long int strideB,
                                                                const cuComplex* beta, /* host or device pointer */
                                                                cuComplex* C,
                                                                int ldc,
                                                                long long int strideC,
                                                                int batchCount);

cublasStatus_t 
cublasZgemmStridedBatched(cublasHandle_t handle,
                        cublasOperation_t transa,
                        cublasOperation_t transb,
                        int m,
                        int n,
                        int k,
                        const cuDoubleComplex* alpha, /* host or device pointer */
                        const cuDoubleComplex* A,
                        int lda,
                        long long int strideA, /* purposely signed */
                        const cuDoubleComplex* B,
                        int ldb,
                        long long int strideB,
                        const cuDoubleComplex* beta, /* host or device poi */
                        cuDoubleComplex* C,
                        int ldc,
                        long long int strideC,
                        int batchCount);


cublasStatus_t  cublasHgemmStridedBatched(cublasHandle_t handle,
                                                            cublasOperation_t transa,
                                                            cublasOperation_t transb,
                                                            int m,
                                                            int n,
                                                            int k,
                                                            const __half* alpha, /* host or device pointer */
                                                            const __half* A,
                                                            int lda,
                                                            long long int strideA, /* purposely signed */
                                                            const __half* B,
                                                            int ldb,
                                                            long long int strideB,
                                                            const __half* beta, /* host or device pointer */
                                                            __half* C,
                                                            int ldc,
                                                            long long int strideC,
                                                            int batchCount);

/* ---------------- CUBLAS BLAS-like extension ---------------- */
/* GEAM */
cublasStatus_t  cublasSgeam(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                const float* alpha, /* host or device pointer */
                                                const float* A,
                                                int lda,
                                                const float* beta, /* host or device pointer */
                                                const float* B,
                                                int ldb,
                                                float* C,
                                                int ldc);

cublasStatus_t  cublasDgeam(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                const double* alpha, /* host or device pointer */
                                                const double* A,
                                                int lda,
                                                const double* beta, /* host or device pointer */
                                                const double* B,
                                                int ldb,
                                                double* C,
                                                int ldc);

cublasStatus_t  cublasCgeam(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                const cuComplex* alpha, /* host or device pointer */
                                                const cuComplex* A,
                                                int lda,
                                                const cuComplex* beta, /* host or device pointer */
                                                const cuComplex* B,
                                                int ldb,
                                                cuComplex* C,
                                                int ldc);

cublasStatus_t  cublasZgeam(cublasHandle_t handle,
                                                cublasOperation_t transa,
                                                cublasOperation_t transb,
                                                int m,
                                                int n,
                                                const cuDoubleComplex* alpha, /* host or device pointer */
                                                const cuDoubleComplex* A,
                                                int lda,
                                                const cuDoubleComplex* beta, /* host or device pointer */
                                                const cuDoubleComplex* B,
                                                int ldb,
                                                cuDoubleComplex* C,
                                                int ldc);

/* Batched LU - GETRF*/
cublasStatus_t  cublasSgetrfBatched(cublasHandle_t handle,
                                                        int n,
                                                        float* const A[], /*Device pointer*/
                                                        int lda,
                                                        int* P,    /*Device Pointer*/
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

cublasStatus_t  cublasDgetrfBatched(cublasHandle_t handle,
                                                        int n,
                                                        double* const A[], /*Device pointer*/
                                                        int lda,
                                                        int* P,    /*Device Pointer*/
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

cublasStatus_t  cublasCgetrfBatched(cublasHandle_t handle,
                                                        int n,
                                                        cuComplex* const A[], /*Device pointer*/
                                                        int lda,
                                                        int* P,    /*Device Pointer*/
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

cublasStatus_t  cublasZgetrfBatched(cublasHandle_t handle,
                                                        int n,
                                                        cuDoubleComplex* const A[], /*Device pointer*/
                                                        int lda,
                                                        int* P,    /*Device Pointer*/
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

/* Batched inversion based on LU factorization from getrf */
cublasStatus_t  cublasSgetriBatched(cublasHandle_t handle,
                                                        int n,
                                                        const float* const A[], /*Device pointer*/
                                                        int lda,
                                                        const int* P,     /*Device pointer*/
                                                        float* const C[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasDgetriBatched(cublasHandle_t handle,
                                                        int n,
                                                        const double* const A[], /*Device pointer*/
                                                        int lda,
                                                        const int* P,      /*Device pointer*/
                                                        double* const C[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasCgetriBatched(cublasHandle_t handle,
                                                        int n,
                                                        const cuComplex* const A[], /*Device pointer*/
                                                        int lda,
                                                        const int* P,         /*Device pointer*/
                                                        cuComplex* const C[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasZgetriBatched(cublasHandle_t handle,
                                                        int n,
                                                        const cuDoubleComplex* const A[], /*Device pointer*/
                                                        int lda,
                                                        const int* P,               /*Device pointer*/
                                                        cuDoubleComplex* const C[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int batchSize);

/* Batched solver based on LU factorization from getrf */

cublasStatus_t  cublasSgetrsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int n,
                                                        int nrhs,
                                                        const float* const Aarray[],
                                                        int lda,
                                                        const int* devIpiv,
                                                        float* const Barray[],
                                                        int ldb,
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasDgetrsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int n,
                                                        int nrhs,
                                                        const double* const Aarray[],
                                                        int lda,
                                                        const int* devIpiv,
                                                        double* const Barray[],
                                                        int ldb,
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasCgetrsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int n,
                                                        int nrhs,
                                                        const cuComplex* const Aarray[],
                                                        int lda,
                                                        const int* devIpiv,
                                                        cuComplex* const Barray[],
                                                        int ldb,
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasZgetrsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int n,
                                                        int nrhs,
                                                        const cuDoubleComplex* const Aarray[],
                                                        int lda,
                                                        const int* devIpiv,
                                                        cuDoubleComplex* const Barray[],
                                                        int ldb,
                                                        int* info,
                                                        int batchSize);

/* TRSM - Batched Triangular Solver */
cublasStatus_t  cublasStrsmBatched(cublasHandle_t handle,
                                                        cublasSideMode_t side,
                                                        cublasFillMode_t uplo,
                                                        cublasOperation_t trans,
                                                        cublasDiagType_t diag,
                                                        int m,
                                                        int n,
                                                        const float* alpha, /*Host or Device Pointer*/
                                                        const float* const A[],
                                                        int lda,
                                                        float* const B[],
                                                        int ldb,
                                                        int batchCount);

cublasStatus_t  cublasDtrsmBatched(cublasHandle_t handle,
                                                        cublasSideMode_t side,
                                                        cublasFillMode_t uplo,
                                                        cublasOperation_t trans,
                                                        cublasDiagType_t diag,
                                                        int m,
                                                        int n,
                                                        const double* alpha, /*Host or Device Pointer*/
                                                        const double* const A[],
                                                        int lda,
                                                        double* const B[],
                                                        int ldb,
                                                        int batchCount);

cublasStatus_t  cublasCtrsmBatched(cublasHandle_t handle,
                                                        cublasSideMode_t side,
                                                        cublasFillMode_t uplo,
                                                        cublasOperation_t trans,
                                                        cublasDiagType_t diag,
                                                        int m,
                                                        int n,
                                                        const cuComplex* alpha, /*Host or Device Pointer*/
                                                        const cuComplex* const A[],
                                                        int lda,
                                                        cuComplex* const B[],
                                                        int ldb,
                                                        int batchCount);

cublasStatus_t  cublasZtrsmBatched(cublasHandle_t handle,
                                                        cublasSideMode_t side,
                                                        cublasFillMode_t uplo,
                                                        cublasOperation_t trans,
                                                        cublasDiagType_t diag,
                                                        int m,
                                                        int n,
                                                        const cuDoubleComplex* alpha, /*Host or Device Pointer*/
                                                        const cuDoubleComplex* const A[],
                                                        int lda,
                                                        cuDoubleComplex* const B[],
                                                        int ldb,
                                                        int batchCount);

/* Batched - MATINV*/
cublasStatus_t  cublasSmatinvBatched(cublasHandle_t handle,
                                                        int n,
                                                        const float* const A[], /*Device pointer*/
                                                        int lda,
                                                        float* const Ainv[], /*Device pointer*/
                                                        int lda_inv,
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

cublasStatus_t  cublasDmatinvBatched(cublasHandle_t handle,
                                                        int n,
                                                        const double* const A[], /*Device pointer*/
                                                        int lda,
                                                        double* const Ainv[], /*Device pointer*/
                                                        int lda_inv,
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

cublasStatus_t  cublasCmatinvBatched(cublasHandle_t handle,
                                                        int n,
                                                        const cuComplex* const A[], /*Device pointer*/
                                                        int lda,
                                                        cuComplex* const Ainv[], /*Device pointer*/
                                                        int lda_inv,
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

cublasStatus_t  cublasZmatinvBatched(cublasHandle_t handle,
                                                        int n,
                                                        const cuDoubleComplex* const A[], /*Device pointer*/
                                                        int lda,
                                                        cuDoubleComplex* const Ainv[], /*Device pointer*/
                                                        int lda_inv,
                                                        int* info, /*Device Pointer*/
                                                        int batchSize);

/* Batch QR Factorization */
cublasStatus_t  cublasSgeqrfBatched(cublasHandle_t handle,
                                                        int m,
                                                        int n,
                                                        float* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        float* const TauArray[], /*Device pointer*/
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasDgeqrfBatched(cublasHandle_t handle,
                                                        int m,
                                                        int n,
                                                        double* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        double* const TauArray[], /*Device pointer*/
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasCgeqrfBatched(cublasHandle_t handle,
                                                        int m,
                                                        int n,
                                                        cuComplex* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        cuComplex* const TauArray[], /*Device pointer*/
                                                        int* info,
                                                        int batchSize);

cublasStatus_t  cublasZgeqrfBatched(cublasHandle_t handle,
                                                        int m,
                                                        int n,
                                                        cuDoubleComplex* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        cuDoubleComplex* const TauArray[], /*Device pointer*/
                                                        int* info,
                                                        int batchSize);
/* Least Square Min only m >= n and Non-transpose supported */
cublasStatus_t  cublasSgelsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int m,
                                                        int n,
                                                        int nrhs,
                                                        float* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        float* const Carray[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int* devInfoArray, /*Device pointer*/
                                                        int batchSize);

cublasStatus_t  cublasDgelsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int m,
                                                        int n,
                                                        int nrhs,
                                                        double* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        double* const Carray[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int* devInfoArray, /*Device pointer*/
                                                        int batchSize);

cublasStatus_t  cublasCgelsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int m,
                                                        int n,
                                                        int nrhs,
                                                        cuComplex* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        cuComplex* const Carray[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int* devInfoArray,
                                                        int batchSize);

cublasStatus_t  cublasZgelsBatched(cublasHandle_t handle,
                                                        cublasOperation_t trans,
                                                        int m,
                                                        int n,
                                                        int nrhs,
                                                        cuDoubleComplex* const Aarray[], /*Device pointer*/
                                                        int lda,
                                                        cuDoubleComplex* const Carray[], /*Device pointer*/
                                                        int ldc,
                                                        int* info,
                                                        int* devInfoArray,
                                                        int batchSize);
/* DGMM */
cublasStatus_t  cublasSdgmm(cublasHandle_t handle,
                                                cublasSideMode_t mode,
                                                int m,
                                                int n,
                                                const float* A,
                                                int lda,
                                                const float* x,
                                                int incx,
                                                float* C,
                                                int ldc);

cublasStatus_t  cublasDdgmm(cublasHandle_t handle,
                                                cublasSideMode_t mode,
                                                int m,
                                                int n,
                                                const double* A,
                                                int lda,
                                                const double* x,
                                                int incx,
                                                double* C,
                                                int ldc);

cublasStatus_t  cublasCdgmm(cublasHandle_t handle,
                                                cublasSideMode_t mode,
                                                int m,
                                                int n,
                                                const cuComplex* A,
                                                int lda,
                                                const cuComplex* x,
                                                int incx,
                                                cuComplex* C,
                                                int ldc);

cublasStatus_t  cublasZdgmm(cublasHandle_t handle,
                                                cublasSideMode_t mode,
                                                int m,
                                                int n,
                                                const cuDoubleComplex* A,
                                                int lda,
                                                const cuDoubleComplex* x,
                                                int incx,
                                                cuDoubleComplex* C,
                                                int ldc);

/* TPTTR : Triangular Pack format to Triangular format */
cublasStatus_t 
cublasStpttr(cublasHandle_t handle, cublasFillMode_t uplo, int n, const float* AP, float* A, int lda);

cublasStatus_t 
cublasDtpttr(cublasHandle_t handle, cublasFillMode_t uplo, int n, const double* AP, double* A, int lda);

cublasStatus_t 
cublasCtpttr(cublasHandle_t handle, cublasFillMode_t uplo, int n, const cuComplex* AP, cuComplex* A, int lda);

cublasStatus_t  cublasZtpttr(
cublasHandle_t handle, cublasFillMode_t uplo, int n, const cuDoubleComplex* AP, cuDoubleComplex* A, int lda);
/* TRTTP : Triangular format to Triangular Pack format */
cublasStatus_t 
cublasStrttp(cublasHandle_t handle, cublasFillMode_t uplo, int n, const float* A, int lda, float* AP);

cublasStatus_t 
cublasDtrttp(cublasHandle_t handle, cublasFillMode_t uplo, int n, const double* A, int lda, double* AP);

cublasStatus_t 
cublasCtrttp(cublasHandle_t handle, cublasFillMode_t uplo, int n, const cuComplex* A, int lda, cuComplex* AP);

cublasStatus_t  cublasZtrttp(
cublasHandle_t handle, cublasFillMode_t uplo, int n, const cuDoubleComplex* A, int lda, cuDoubleComplex* AP);



%inline %{ 

Cublas  _cublas;
Cublas *cublas = &_cublas;

void CreateCublas() { 
    cublas = new Cublas();
    assert(cublas != NULL);
}
void DeleteCublas() { 
    if(cublas) delete cublas;
    cublas = NULL;
}
void synchronize() {
    cudaDeviceSynchronize();
}
unsigned seed = 0;

void set_seed()
{
    typedef std::chrono::high_resolution_clock myclock;
    myclock::time_point beginning = myclock::now();
    myclock::duration d = myclock::now() - beginning;    
    seed = d.count();
}


template<typename T>
float cast_float(T val) { return (float)val; }

template<typename T>
double cast_double(T val) { return (double)val; }

template<typename T>
int8_t cast_int8(T val) { return (int8_t)val; }

template<typename T>
uint8_t cast_uint8(T val) { return (uint8_t)val; }

template<typename T>
int16_t cast_int16(T val) { return (int16_t)val; }

template<typename T>
uint16_t cast_uint16(T val) { return (uint16_t)val; }

template<typename T>
int32_t cast_int32(T val) { return (int32_t)val; }

template<typename T>
uint32_t cast_uint32(T val) { return (uint32_t)val; }

template<typename T>
int64_t cast_int64(T val) { return (int64_t)val; }

template<typename T>
uint64_t cast_uint64(T val) { return (uint64_t)val; }

std::vector<float> vector_range(int start, int end, int inc=1) {
    std::vector<float> r;    
    for(int i = start; i <= end; i+=inc) {
        r.push_back((float)i);
    }
    return r;
}

%}

// lua only has double
%template(cast_double_float) cast_double<float>;