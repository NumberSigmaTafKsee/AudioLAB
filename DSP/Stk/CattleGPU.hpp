/*
#ifdef CATTL3_USE_CUDA
#include "parameter_initialization/gpu/ConstantGPUParameterInitialization.hpp"
#include "parameter_initialization/gpu/GaussianGPUParameterInitialization.hpp"
#include "parameter_initialization/gpu/GlorotGPUParameterInitialization.hpp"
#include "parameter_initialization/gpu/HeGPUParameterInitialization.hpp"
#include "parameter_initialization/gpu/LeCunGPUParameterInitialization.hpp"
#include "parameter_initialization/gpu/OneGPUParameterInitialization.hpp"
#include "parameter_initialization/gpu/ZeroGPUParameterInitialization.hpp"
#include "parameter_regularization/gpu/L2GPUParameterRegularization.hpp"
#endif
*/
#pragma once
#include <cublas_v2.h>
#include <cuda_runtime.h>
#include <cudnn.h>
#include <curand.h>

#define cublasErrorCheck(status) { _cublas_error_check(status, __FILE__, __LINE__); }
#define cublasAssert(status) { _cublas_assert(status, __FILE__, __LINE__); }

namespace cattle {
    namespace gpu
    {

        /**
        * A class representing a cuBLAS runtime error.
        */
        class CuBLASError : public std::runtime_error {
        public:
            /**
            * @param status The cuBLAS status code.
            * @param file The name of the file in which the error occurred.
            * @param line The number of the line at which the error occurred.
            */
            CuBLASError(cublasStatus_t status, const char* file, int line) :
                std::runtime_error("cuBLAS Error: " + cublas_status_to_string(status) + "; File: " +
                        std::string(file) + "; Line: " + std::to_string(line)) { }
        private:
            /**
            * It returns a string representation of the provided cuBLAS status code.
            *
            * @param status The cuBLAS status code.
            * @return A string describing the status code.
            */
            inline static std::string cublas_status_to_string(cublasStatus_t status) {
                switch (status) {
                    case CUBLAS_STATUS_SUCCESS: return "CUBLAS_STATUS_SUCCESS";
                    case CUBLAS_STATUS_NOT_INITIALIZED: return "CUBLAS_STATUS_NOT_INITIALIZED";
                    case CUBLAS_STATUS_ALLOC_FAILED: return "CUBLAS_STATUS_ALLOC_FAILED";
                    case CUBLAS_STATUS_INVALID_VALUE: return "CUBLAS_STATUS_INVALID_VALUE";
                    case CUBLAS_STATUS_ARCH_MISMATCH: return "CUBLAS_STATUS_ARCH_MISMATCH";
                    case CUBLAS_STATUS_MAPPING_ERROR: return "CUBLAS_STATUS_MAPPING_ERROR";
                    case CUBLAS_STATUS_EXECUTION_FAILED: return "CUBLAS_STATUS_EXECUTION_FAILED";
                    case CUBLAS_STATUS_INTERNAL_ERROR: return "CUBLAS_STATUS_INTERNAL_ERROR";
                    case CUBLAS_STATUS_NOT_SUPPORTED: return "CUBLAS_STATUS_NOT_SUPPORTED";
                    case CUBLAS_STATUS_LICENSE_ERROR: return "CUBLAS_STATUS_LICENSE_ERROR";
                    default: return "<unknown>";
                }
            }
        };

        namespace {

            __inline__ void _cublas_error_check(cublasStatus_t status, const char* file, int line) {
                if (status != CUBLAS_STATUS_SUCCESS)
                    throw CuBLASError(status, file, line);
            }

            __inline__ void _cublas_assert(cublasStatus_t status, const char* file, int line) {
                try {
                    _cublas_error_check(status, file, line);
                } catch (const CuBLASError& e) {
                    std::cout << e.what() << std::endl;
                    exit(-1);
                }
            }

        }

        /**
        * A singleton cuBLAS handle class.
        */
        class CuBLASHandle {
        public:
            inline CuBLASHandle(const CuBLASHandle&) = delete;
            inline ~CuBLASHandle() {
                cublasAssert(cublasDestroy(handle));
            }
            inline CuBLASHandle& operator=(const CuBLASHandle&) = delete;
            inline operator const cublasHandle_t&() const {
                return handle;
            }
            /**
            * @return A reference to the only instance of the class.
            */
            inline static const CuBLASHandle& get_instance() {
                static CuBLASHandle instance;
                return instance;
            }
        private:
            cublasHandle_t handle;
            inline CuBLASHandle() :
                    handle() {
                cublasAssert(cublasCreate(&handle));
            }
        };

        /**
        * A class representing a CUDA runtime error.
        */
        class CUDAError : public std::runtime_error {
        public:
            /**
            * @param code The CUDA error code.
            * @param file The name of the file in which the error occurred.
            * @param line The number of the line at which the error occurred.
            */
            CUDAError(cudaError_t code, const char* file, int line) :
                std::runtime_error("CUDA Error: " + std::string(cudaGetErrorString(code)) + "; File: " +
                        std::string(file) + "; Line: " + std::to_string(line)) { }
        };

        namespace {

            __inline__ void _cuda_error_check(const char* file, int line, cudaError_t code = cudaGetLastError()) {
                if (code != cudaSuccess)
                    throw CUDAError(code, file, line);
            }

            __inline__ void _cuda_assert(const char* file, int line, cudaError_t code = cudaGetLastError()) {
                try {
                    _cuda_error_check(file, line, code);
                } catch (const CUDAError& e) {
                    std::cout << e.what() << std::endl;
                    exit(-1);
                }
            }
        }

        /**
        * A template class for CUDA device arrays of different data types.
        */
        template<typename Scalar>
        class CUDAArray {
            typedef CUDAArray<Scalar> Self;
        public:
            /**
            * @param data The device array to wrap in a CUDAArray. The ownership of the
            * pointer is not transfered to the created CUDAArray instance.
            * @param size The number of elements the array is to have.
            */
            inline CUDAArray(Scalar* data, std::size_t size) :
                    _size(size),
                    _data(data),
                    _wrapper(true) {
                assert(data || size == 0);
            }
            /**
            * @param size The number of elements the array is to have.
            */
            inline CUDAArray(std::size_t size) :
                    _size(size),
                    _data(nullptr),
                    _wrapper(false) {
                if (size > 0)
                    cudaAssert(cudaMalloc(&_data, size * sizeof(Scalar)));
            }
            inline CUDAArray() :
                    CUDAArray(0u) { }
            inline CUDAArray(const Self& array) :
                    CUDAArray(array._size),
                    _wrapper(false) {
                if (_size > 0)
                    cudaAssert(cudaMemcpy(_data, array._data, _size * sizeof(Scalar), cudaMemcpyDeviceToDevice));
            }
            inline CUDAArray(Self&& array) :
                    CUDAArray() {
                swap(*this, array);
            }
            inline virtual ~CUDAArray() {
                if (_size > 0 && !_wrapper)
                    cudaAssert(cudaFree(_data));
            }
            inline Self& operator=(Self array) {
                swap(*this, array);
                return *this;
            }
            /**
            * @return The size of the array.
            */
            inline std::size_t size() const {
                return _size;
            }
            /**
            * @return A device memory address pointer pointing to the (constant) first element of
            * the array.
            */
            inline const Scalar* data() const {
                return _data;
            }
            /**
            * @return A device memory address pointer pointing to the first element of the array.
            */
            inline Scalar* data() {
                return _data;
            }
            /**
            * @return Whether the instance is just a wrapper over a device array that it does not
            * own.
            */
            inline bool wrapper() const {
                return _wrapper;
            }
            /**
            * @param value The integer value to which the values of the array are to be set.
            */
            inline void set_values(int value) {
                if (_size > 0)
                    cudaAssert(cudaMemset(_data, value, _size * sizeof(Scalar)));
            }
            /**
            * It populates the entire device array with data from the host memory.
            *
            * @param host_array A pointer to the first element of the host array.
            */
            inline void copy_from_host(const Scalar* host_array) {
                if (_size > 0)
                    cudaAssert(cudaMemcpy(_data, host_array, _size * sizeof(Scalar), cudaMemcpyHostToDevice));
            }
            /**
            * It copies the entire device array to the host memory.
            *
            * @param host_array A pointer pointing to the beginning of a contiguous host memory
            * block to which the device tensor is to be copied.
            */
            inline void copy_to_host(Scalar* host_array) const {
                if (_size > 0)
                    cudaAssert(cudaMemcpy(host_array, _data, _size * sizeof(Scalar), cudaMemcpyDeviceToHost));
            }
            inline friend void swap(Self& array1, Self& array2) {
                using std::swap;
                swap(array1._size, array2._size);
                swap(array1._data, array2._data);
                swap(array1._wrapper, array2._wrapper);
            }
        private:
            std::size_t _size;
            Scalar* _data;
            bool _wrapper;
        };

        namespace {

        template<typename Scalar>
        using AsumRoutine = cublasStatus_t (*)(cublasHandle_t, int, const Scalar*, int, Scalar*);
        template<typename Scalar> __inline__ AsumRoutine<Scalar> asum_routine() { return &cublasDasum; }
        template<> __inline__ AsumRoutine<float> asum_routine() { return &cublasSasum; }

        template<typename Scalar>
        using Nrm2Routine = cublasStatus_t (*)(cublasHandle_t, int, const Scalar*, int, Scalar*);
        template<typename Scalar> __inline__ Nrm2Routine<Scalar> nrm2_routine() { return &cublasDnrm2; }
        template<> __inline__ Nrm2Routine<float> nrm2_routine() { return &cublasSnrm2; }

        template<typename Scalar>
        using ScalRoutine = cublasStatus_t (*)(cublasHandle_t, int, const Scalar*, Scalar*, int);
        template<typename Scalar> __inline__ ScalRoutine<Scalar> scal_routine() { return &cublasDscal; }
        template<> __inline__ ScalRoutine<float> scal_routine() { return &cublasSscal; }

        template<typename Scalar>
        using AxpyRoutine = cublasStatus_t (*)(cublasHandle_t, int, const Scalar*, const Scalar*, int, Scalar*, int);
        template<typename Scalar> __inline__ AxpyRoutine<Scalar> axpy_routine() { return &cublasDaxpy; }
        template<> __inline__ AxpyRoutine<float> axpy_routine() { return &cublasSaxpy; }

        template<typename Scalar>
        using GemmRoutine = cublasStatus_t (*)(cublasHandle_t, cublasOperation_t, cublasOperation_t, int, int, int,
                const Scalar*, const Scalar*, int, const Scalar*, int, const Scalar*, Scalar*, int);
        template<typename Scalar> __inline__ GemmRoutine<Scalar> gemm_routine() { return &cublasDgemm; }
        template<> __inline__ GemmRoutine<float> gemm_routine() { return &cublasSgemm; }

        }

        /**
        * A template class for column-major cuBLAS device matrices of different data types.
        */
        template<typename Scalar>
        class CuBLASMatrix : public CUDAArray<Scalar> {
            static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
            typedef CUDAArray<Scalar> Base;
            typedef CuBLASMatrix<Scalar> Self;
        public:
            /**
            * @param data The device array to wrap in a CuBLASMatrix. The ownership of the
            * pointer is not transfered to the matrix.
            * @param rows The number of rows of the matrix.
            * @param cols The number of columns of the matrix.
            */
            inline CuBLASMatrix(Scalar* data, std::size_t rows, std::size_t cols) :
                    Base(data, rows * cols),
                    _rows(rows),
                    _cols(cols) { }
            /**
            * @param rows The number of rows of the matrix.
            * @param cols The number of columns of the matrix.
            */
            inline CuBLASMatrix(std::size_t rows, std::size_t cols) :
                    Base(rows * cols),
                    _rows(rows),
                    _cols(cols) { }
            inline CuBLASMatrix() :
                    CuBLASMatrix(0u, 0u) { }
            inline CuBLASMatrix(const Matrix<Scalar>& matrix) :
                    CuBLASMatrix(matrix.rows(), matrix.cols()) {
                if (Base::size() > 0)
                    Base::copy_from_host(matrix.data());
            }
            inline CuBLASMatrix(const Self& matrix) :
                    Base(matrix),
                    _rows(matrix._rows),
                    _cols(matrix._cols) { }
            inline CuBLASMatrix(Self&& matrix) :
                    CuBLASMatrix() {
                swap(*this, matrix);
            }
            virtual ~CuBLASMatrix() = default;
            inline Self& operator=(Self matrix) {
                swap(*this, matrix);
                return *this;
            }
            inline operator Matrix<Scalar>() const {
                if (Base::size() == 0)
                    return Matrix<Scalar>();
                Matrix<Scalar> out(_rows, _cols);
                Base::copy_to_host(out.data());
                return out;
            }
            /**
            * @return The number of rows of the matrix.
            */
            inline std::size_t rows() const {
                return _rows;
            }
            /**
            * @return The number of columns of the matrix.
            */
            inline std::size_t cols() const {
                return _cols;
            }
            /**
            * @return The L1 norm of the matrix.
            */
            inline Scalar l1_norm() const {
                Scalar res;
                asum(Base::size(), 1, *this, &res);
                return res;
            }
            /**
            * @return The L2 (Euclidian) norm of the matrix.
            */
            inline Scalar l2_norm() const {
                Scalar res;
                nrm2(Base::size(), 1, *this, &res);
                return res;
            }
            inline Self& operator+=(const Self& rhs) {
                axpy(Base::size(), 1, rhs, 1, 1, *this);
                return *this;
            }
            inline Self& operator-=(const Self& rhs) {
                axpy(Base::size(), 1, rhs, -1, 1, *this);
                return *this;
            }
            inline Self& operator*=(const Self& rhs) {
                gemm(*this, false, rhs, false, 1, 1, *this);
                return *this;
            }
            inline Self operator*=(Scalar rhs) {
                scal(Base::size(), rhs, 1, *this);
                return *this;
            }
            inline Self operator*=(Scalar rhs) {
                return *this * (1 / rhs);
            }
            inline friend Self operator+(Self lhs, const Self& rhs) {
                return lhs += rhs;
            }
            inline friend Self operator-(Self lhs, const Self& rhs) {
                return lhs -= rhs;
            }
            inline friend Self operator*(Self lhs, const Self& rhs) {
                return lhs *= rhs;
            }
            inline friend Self operator*(Self lhs, Scalar rhs) {
                return lhs *= rhs;
            }
            inline friend Self operator/(Self lhs, Scalar rhs) {
                return lhs /= rhs;
            }
            inline friend void swap(Self& matrix1, Self& matrix2) {
                using std::swap;
                swap(static_cast<Base&>(matrix1), static_cast<Base&>(matrix2));
                swap(matrix1._rows, matrix2._rows);
                swap(matrix1._cols, matrix2._cols);
            }
            /**
            * It computes the sum of the absolute values of the matrix's coefficients.
            *
            * \f$R = \sum\limits_{i = 1}^n \left|A_i\right|\f$
            *
            * @param n The number of elements whose absolute value is to be summed up.
            * @param inc_a The stride between elements of the matrix.
            * @param a The matrix.
            * @param result The result of the computation.
            */
            inline static void asum(int n, int inc_a, const CuBLASMatrix<Scalar>& a, /* out */ Scalar& result) {
                AsumRoutine<Scalar> asum = asum_routine<Scalar>();
                cublasAssert(asum(CuBLASHandle::get_instance(), n, a.data(), inc_a, &result));
            }
            /**
            * It computes the Euclidian norm of the matrix's coefficients.
            *
            * \f$R = \sqrt{\sum\limits_{i = 1}^n A_i^2}\f$
            *
            * @param n The number of elements whose second norm is to be calculated.
            * @param inc_a The stride between elements of the matrix.
            * @param a The matrix.
            * @param result The result of the computation.
            */
            inline static void nrm2(int n, int inc_a, const CuBLASMatrix<Scalar>& a, /* out */ Scalar& result) {
                Nrm2Routine<Scalar> nrm2 = nrm2_routine<Scalar>();
                cublasAssert(nrm2(CuBLASHandle::get_instance(), n, a.data(), inc_a, &result));
            }
            /**
            * It scales the matrix by the specified factor.
            *
            * \f$A_i = \alpha * A_i, \forall i \in \{0,...,n\}\f$
            *
            * @param n The number of elements on which the operation is to be performed.
            * @param alpha The scaling factor.
            * @param inc_a The stride between elements of the matrix.
            * @param a The matrix to be scaled.
            */
            inline static void scal(int n, Scalar alpha, int inc_a, /* in/out */ CuBLASMatrix<Scalar>& a) {
                ScalRoutine<Scalar> scal = scal_routine<Scalar>();
                cublasAssert(scal(CuBLASHandle::get_instance(), n, &alpha, a.data(), inc_a));
            }
            /**
            * It adds a scaled matrix to another matrix.
            *
            * \f$B_i = \alpha * A_i + B_i, \forall i \in \{0,...,n\}\f$
            *
            * @param n The number of elements on which the operation is to be performed.
            * @param inc_a The stride between elements of the scaled matrix.
            * @param a The matrix to scale and add to the other matrix.
            * @param alpha The scaling factor.
            * @param inc_b The stride between elements of the target matrix.
            * @param b The traget matrix.
            */
            inline static void axpy(int n, int inc_a, const CuBLASMatrix<Scalar>& a, Scalar alpha,
                    int inc_b, /* in/out */ CuBLASMatrix<Scalar>& b) {
                assert(a.rows() == b.rows() && a.cols() == b.cols());
                AxpyRoutine<Scalar> axpy = axpy_routine<Scalar>();
                cublasAssert(axpy(CuBLASHandle::get_instance(), n, &alpha, a.data(), inc_a, b.data(), inc_b));
            }
            /**
            * It computes the product of the matrix multiplication.
            *
            * \f$C = \alpha * (A \textrm{ x } B) + \beta * C\f$
            *
            * @param a The multiplicand matrix.
            * @param transpose_a Whether the multiplicand is to be transposed for the operation.
            * @param b The multiplier matrix.
            * @param transpose_b Whether the multiplier is to be transposed for the operation.
            * @param alpha The factor by which the result of the multiplication is to be scaled.
            * @param beta The factor by which c is to be scaled before the result of the multiplication
            * is added to it.
            * @param c The product of the matrix multiplication.
            */
            inline static void gemm(const CuBLASMatrix<Scalar>& a, bool transpose_a, const CuBLASMatrix<Scalar>& b,
                    bool transpose_b, Scalar alpha, Scalar beta, /* out */ CuBLASMatrix<Scalar>& c) {
                std::size_t a_rows, a_cols, b_rows, b_cols;
                if (transpose_a) {
                    a_rows = a.cols();
                    a_cols = a.rows();
                } else {
                    a_rows = a.rows();
                    a_cols = a.cols();
                }
                if (transpose_b) {
                    b_rows = b.cols();
                    b_cols = b.rows();
                } else {
                    b_rows = b.rows();
                    b_cols = b.cols();
                }
                assert(a_cols == b_rows);
                assert(c.rows() == a_rows() && c.cols() == b_cols);
                cublasOperation_t a_op = transpose_a ? CUBLAS_OP_T : CUBLAS_OP_N;
                cublasOperation_t b_op = transpose_b ? CUBLAS_OP_T : CUBLAS_OP_N;
                // Resolve the GEMM precision based on the scalar type.
                GemmRoutine<Scalar> gemm = gemm_routine<Scalar>();
                // Perform the matrix multiplication.
                cublasAssert(gemm(CuBLASHandle::get_instance(), a_op, b_op, a_rows, b_cols, a_cols, &alpha,
                        a.data(), a.rows(), b.data(), b.rows(), &beta, c.data(), a_rows));
            }
        private:
            std::size_t _rows, _cols;
        };

        /**
        * A class representing a cuDNN runtime error.
        */
        class CuDNNError : public std::runtime_error {
        public:
            /**
            * @param status The cuDNN status code.
            * @param file The name of the file in which the error occurred.
            * @param line The number of the line at which the error occurred.
            */
            CuDNNError(cudnnStatus_t status, const char* file, int line) :
                std::runtime_error("cuDNN Error: " + std::string(cudnnGetErrorString(status)) + "; File: " +
                        std::string(file) + "; Line: " + std::to_string(line)) { }
        };

        namespace {

        __inline__ void _cudnn_error_check(cudnnStatus_t status, const char* file, int line) {
            if (status != CUDNN_STATUS_SUCCESS)
                throw CuDNNError(status, file, line);
        }

        __inline__ void _cudnn_assert(cudnnStatus_t status, const char* file, int line) {
            try {
                _cudnn_error_check(status, file, line);
            } catch (const CuDNNError& e) {
                std::cout << e.what() << std::endl;
                exit(-1);
            }
        }

        }

        /**
        * A singleton utility class representing a handle to the cuDNN library.
        */
        class CuDNNHandle {
        public:
            CuDNNHandle(const CuDNNHandle&) = delete;
            ~CuDNNHandle() {
                cudnnAssert(cudnnDestroy(handle));
            }
            CuDNNHandle& operator=(const CuDNNHandle&) = delete;
            inline operator const cudnnHandle_t&() const {
                return handle;
            }
            /**
            * @return A reference to the only instance of the class.
            */
            inline static const CuDNNHandle& get_instance() {
                static CuDNNHandle instance;
                return instance;
            }
        private:
            inline CuDNNHandle() :
                    handle() {
                cudnnAssert(cudnnCreate(&handle));
            }
            cudnnHandle_t handle;
        };

        /**
        * A template class for representing row-major cuDNN device tensors of different data types.
        */
        template<typename Scalar>
        class CuDNNTensor : public CUDAArray<Scalar> {
            static_assert(std::is_floating_point<Scalar>::value, "non floating-point scalar type");
            typedef CUDAArray<Scalar> Base;
            typedef CuDNNTensor<Scalar> Self;
        public:
            static constexpr cudnnDataType_t DATA_TYPE = std::is_same<Scalar,float>::value ?
                    CUDNN_DATA_FLOAT : CUDNN_DATA_DOUBLE;
            static constexpr cudnnTensorFormat_t TENSOR_FORMAT = CUDNN_TENSOR_NCHW;
            static constexpr cudnnNanPropagation_t NAN_PROP = CUDNN_PROPAGATE_NAN;
            /**
            * @param data The device array to wrap in a CuDNNTensor. The ownership of the
            * pointer is not transfered to the tensor.
            * @param samples The batch size.
            * @param height The height.
            * @param width The width.
            * @param channels The number of channels.
            */
            inline CuDNNTensor(Scalar* data, std::size_t samples, std::size_t height, std::size_t width,
                    std::size_t channels) :
                        Base(data, samples * height * width * channels),
                        _samples(samples),
                        _height(height),
                        _width(width),
                        _channels(channels),
                        _desc(),
                        _filter_desc() {
                if (Base::size() > 0) {
                    create_tensor_descriptor(_desc, samples, height, width, channels);
                    create_filter_descriptor(_filter_desc, samples, height, width, channels);
                }
            }
            /**
            * @param samples The batch size.
            * @param height The height.
            * @param width The width.
            * @param channels The number of channels.
            */
            inline CuDNNTensor(std::size_t samples, std::size_t height, std::size_t width, std::size_t channels) :
                    Base(samples * height * width * channels),
                    _samples(samples),
                    _height(height),
                    _width(width),
                    _channels(channels),
                    _desc(),
                    _filter_desc() {
                if (Base::size() > 0) {
                    create_tensor_descriptor(_desc, samples, height, width, channels);
                    create_filter_descriptor(_filter_desc, samples, height, width, channels);
                }
            }
            inline CuDNNTensor() :
                    CuDNNTensor(0u, 0u, 0u, 0u) { }
            inline CuDNNTensor(const Tensor<Scalar,4>& tensor) :
                    CuDNNTensor(tensor.dimension(0), tensor.dimension(1), tensor.dimension(2), tensor.dimension(3)) {
                if (Base::size() > 0) {
                    static std::array<std::size_t,4> eigen_to_cudnn_layout({ 2u, 1u, 2u, 0u });
                    Tensor<Scalar,4> shuffled_tensor = tensor.shuffle(eigen_to_cudnn_layout);
                    Base::copy_from_host(shuffled_tensor.data());
                }
            }
            inline CuDNNTensor(const Self& tensor) :
                    Base(tensor),
                    _samples(tensor._samples),
                    _height(tensor._height),
                    _width(tensor._width),
                    _channels(tensor._channels),
                    _desc(tensor._desc),
                    _filter_desc(tensor._filter_desc) { }
            inline CuDNNTensor(Self&& tensor) :
                    CuDNNTensor() {
                swap(*this, tensor);
            }
            inline ~CuDNNTensor() {
                if (Base::size() > 0) {
                    destroy_tensor_descriptor(_desc);
                    destroy_filter_descriptor(_filter_desc);
                }
            }
            inline Self& operator=(Self tensor) {
                swap(*this, tensor);
                return *this;
            }
            inline operator Tensor<Scalar,4>() const {
                if (Base::size() == 0)
                    return Tensor<Scalar,4>();
                Tensor<Scalar,4> out(_width, _height, _channels, _samples);
                Base::copy_to_host(out.data());
                static std::array<std::size_t,4> cudnn_to_eigen_layout({ 3u, 1u, 0u, 2u });
                return out.shuffle(cudnn_to_eigen_layout);
            }
            /**
            * @return The batch size of the tensor.
            */
            inline std::size_t samples() const {
                return _samples;
            }
            /**
            * @return The height of the tensor.
            */
            inline std::size_t height() const {
                return _height;
            }
            /**
            * @return The width of the tensor.
            */
            inline std::size_t width() const {
                return _width;
            }
            /**
            * @return The number of channels of the tensor.
            */
            inline std::size_t channels() const {
                return _channels;
            }
            /**
            * @return A constant reference to the tensor descriptor.
            */
            inline const cudnnTensorDescriptor_t& desc() const {
                return _desc;
            }
            /**
            * @return A constant reference to the filter descriptor.
            */
            inline const cudnnFilterDescriptor_t& filter_desc() const {
                return _filter_desc;
            }
            /**
            * @param value The value to which all elemetnts of the tensor are to be set.
            */
            inline void set_values(Scalar value) {
                cudnnAssert(cudnnSetTensor(CuDNNHandle::get_instance(), _desc, _data, value));
            }
            /**
            * Performs a reduction along all ranks of the tensor.
            *
            * @param op_type The reduction operation type.
            * @return The result of the reduction.
            */
            inline Scalar reduce(cudnnReduceTensorOp_t op_type) const {
                Self reduced_tensor(1u, 1u, 1u, 1u);
                reduce_op(1, *this, op_type, 0, reduced_tensor);
                Scalar res;
                reduced_tensor.copy_to_host(&res);
                return res;
            }
            /**
            * Performs a reduction along the specified ranks of the tensor.
            *
            * @param op_type The reduction operation type.
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The result of the reduction.
            */
            inline Self reduce(cudnnReduceTensorOp_t op_type, const std::array<bool,4>& ranks) const {
                Self reduced_tensor(ranks[0] ? 1u : _samples, ranks[1] ? 1u : _height,
                        ranks[2] ? 1u : _width, ranks[3] ? 1u : _channels);
                reduce_op(1, *this, op_type, 0, reduced_tensor);
                return reduced_tensor;
            }
            /**
            * @return The sum of all elements of the tensor.
            */
            inline Scalar sum() const {
                return reduce(CUDNN_REDUCE_TENSOR_ADD);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self sum(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_ADD, ranks);
            }
            /**
            * @return The mean of all elements of the tensor.
            */
            inline Scalar avg() const {
                return reduce(CUDNN_REDUCE_TENSOR_AVG);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self avg(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_AVG, ranks);
            }
            /**
            * @return The minimum of all elements of the tensor.
            */
            inline Scalar min() const {
                return reduce(CUDNN_REDUCE_TENSOR_MIN);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self min(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_MIN, ranks);
            }
            /**
            * @return The maximum of all elements of the tensor.
            */
            inline Scalar max() const {
                return reduce(CUDNN_REDUCE_TENSOR_MAX);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self max(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_MAX, ranks);
            }
            /**
            * @return The absolute maximum of all elements of the tensor.
            */
            inline Scalar abs_max() const {
                return reduce(CUDNN_REDUCE_TENSOR_AMAX);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self abs_max(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_AMAX, ranks);
            }
            /**
            * @return The L1 norm of the tensor.
            */
            inline Scalar l1_norm() const {
                return reduce(CUDNN_REDUCE_TENSOR_NORM1);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self l1_norm(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_NORM1, ranks);
            }
            /**
            * @return The L2 norm of the tensor.
            */
            inline Scalar l2_norm() const {
                return reduce(CUDNN_REDUCE_TENSOR_NORM2);
            }
            /**
            * @param ranks A size 4 boolean array denoting which ranks are to be reduced.
            * @return The reduced tensor.
            */
            inline Self l2_norm(const std::array<bool,4>& ranks) const {
                return reduce(CUDNN_REDUCE_TENSOR_NORM2, ranks);
            }
            /**
            * @return A string representation of the tensor.
            */
            std::string to_string() const {
                std::stringstream strm;
                strm << "data type: " << DATA_TYPE << "; format: " << TENSOR_FORMAT << "; " <<
                        "[N:" << _samples << ", C:" << _channels << ", H:" << _height << ", W:" << _width << "]";
                return strm.str();
            }
            inline Self& operator+=(const Self& rhs) {
                add(1, rhs, 1, *this);
                return *this;
            }
            inline Self& operator-=(const Self& rhs) {
                add(-1, rhs, 1, *this);
                return *this;
            }
            inline Self& operator*=(const Self& rhs) {
                op(*this, 1, rhs, 1, CUDNN_OP_TENSOR_MUL, 1, *this);
                return *this;
            }
            inline Self& operator+=(Scalar rhs) {
                Self rhs_tensor(1u, 1u, 1u, 1u);
                rhs_tensor.copy_from_host(&rhs);
                op(*this, 1, rhs_tensor, 1, CUDNN_OP_TENSOR_ADD, 1, *this);
                return *this;
            }
            inline Self& operator-=(Scalar rhs) {
                return *this += -rhs;
            }
            inline Self& operator*=(Scalar rhs) {
                scale(rhs, *this);
                return *this;
            }
            inline Self& operator/=(Scalar rhs) {
                return *this *= (1 / rhs);
            }
            inline friend Self operator+(Self lhs, const Self& rhs) {
                return lhs += rhs;
            }
            inline friend Self operator-(Self lhs, const Self& rhs) {
                return lhs -= rhs;
            }
            inline friend Self operator*(Self lhs, const Self& rhs) {
                return lhs *= rhs;
            }
            inline friend Self operator+(Self lhs, Scalar rhs) {
                return lhs += rhs;
            }
            inline friend Self operator-(Self lhs, Scalar rhs) {
                return lhs -= rhs;
            }
            inline friend Self operator*(Self lhs, Scalar rhs) {
                return lhs *= rhs;
            }
            inline friend Self operator/(Self lhs, Scalar rhs) {
                return lhs /= rhs;
            }
            inline friend std::ostream& operator<<(std::ostream& os, const Self& tensor) {
                return os << tensor.to_string() << std::endl;
            }
            inline friend void swap(Self& tensor1, Self& tensor2) {
                using std::swap;
                swap(static_cast<Base&>(tensor1), static_cast<Base&>(tensor2));
                swap(tensor1._samples, tensor2._samples);
                swap(tensor1._height, tensor2._height);
                swap(tensor1._width, tensor2._width);
                swap(tensor1._channels, tensor2._channels);
                swap(tensor1._filter, tensor2._filter);
                swap(tensor1._desc, tensor2._desc);
                swap(tensor1._filter_desc, tensor2._filter_desc);
            }
            /**
            * @param desc A reference to the tensor descri ptor object.
            * @param samples The batch size.
            * @param height The height.
            * @param width The width.
            * @param channels The number of channels.
            */
            inline static void create_tensor_descriptor(cudnnTensorDescriptor_t& desc, std::size_t samples,
                    std::size_t height, std::size_t width, std::size_t channels) {
                cudnnAssert(cudnnCreateTensorDescriptor(&desc));
                cudnnAssert(cudnnSetTensor4dDescriptor(desc, TENSOR_FORMAT, DATA_TYPE, samples, channels,
                        height, width));
            }
            /**
            * @param desc A constant reference to the tensor descriptor object.
            */
            inline static void destroy_tensor_descriptor(const cudnnTensorDescriptor_t& desc) {
                cudnnAssert(cudnnDestroyTensorDescriptor(desc));
            }
            /**
            * @param filter_desc A reference to the filter descriptor object.
            * @param samples The batch size.
            * @param height The height.
            * @param width The width.
            * @param channels The number of channels.
            */
            inline static void create_filter_descriptor(cudnnFilterDescriptor_t& filter_desc, std::size_t samples,
                    std::size_t height, std::size_t width, std::size_t channels) {
                cudnnAssert(cudnnCreateFilterDescriptor(&filter_desc));
                cudnnAssert(cudnnSetFilter4dDescriptor(filter_desc, Base::DATA_TYPE, Base::TENSOR_FORMAT, samples,
                        channels, height, width));
            }
            /**
            * @param filter_desc A constant reference to the filter descriptor object.
            */
            inline static void destroy_filter_descriptor(const cudnnFilterDescriptor_t& filter_desc) {
                cudnnAssert(cudnnDestroyFilterDescriptor(filter_desc));
            }
            /**
            * It scales the specified tensor by a certain factor.
            *
            * \f$A = \alpha * B\f$
            *
            * @param alpha The factor by which the tensor is to be scaled.
            * @param a The tensor to scale.
            */
            inline static void scale(Scalar alpha, /* in/out */ Self& a) {
                cudnnAssert(cudnnScaleTensor(CuDNNHandle::get_instance(), a.desc(), a.data(), &alpha));
            }
            /**
            * It adds tensor a to tensor b.
            *
            * \f$B = \alpha * A + \beta * B\f$
            *
            * @param alpha The scaling factor of the tensor to add to the other one.
            * @param a The tensor to add to the other tensor.
            * @param beta The scaling factor of the target tensor.
            * @param b The target tensor.
            */
            inline static void add(Scalar alpha, const Self& a, Scalar beta, /* in/out */ Self& b) {
                cudnnAssert(cudnnAddTensor(CuDNNHandle::get_instance(), &alpha, a.desc(), a.data(), &beta,
                        b.desc(), b.data()));
            }
            /**
            * It performs the specified operation on tensors a and b and saves the result in c.
            *
            * \f$C = op(\alpha * A, \beta * B) + \gamma * C\f$
            *
            * @param a The first operand.
            * @param alpha The scaling factor of the first operand.
            * @param b The second operand.
            * @param beta The scaling factor of the second operand.
            * @param op_type The operation type.
            * @param gamma The scaling factor of the result tensor.
            * @param c The result tensor.
            */
            inline static void op(const Self& a, Scalar alpha, const Self& b, Scalar beta, cudnnOpTensorOp_t op_type,
                    Scalar gamma, /* in/out */ Self& c) {
                cudnnOpTensorDescriptor_t desc;
                cudnnAssert(cudnnCreateOpTensorDescriptor(&desc));
                cudnnAssert(cudnnSetOpTensorDescriptor(desc, op_type, DATA_TYPE, NAN_PROP));
                cudnnAssert(cudnnOpTensor(desc, alpha, a.desc(), a.data(), beta, b.desc(), b.data(),
                        gamma, c.desc(), c.data()));
                cudnnAssert(cudnnDestroyOpTensorDescriptor(desc));
            }
            /**
            * It performs the specified reduction operation on tensor a and adds it to tensor b.
            *
            * \f$B = \alpha * reduce_op(A) + \beta * B\f$
            *
            * @param alpha The factor by which the result of the reduction is to be scaled.
            * @param a The tensor to reduce.
            * @param op_type The reduction operation type.
            * @param beta The scaling factor of the target tensor.
            * @param b The target tensor.
            */
            inline static void reduce_op(Scalar alpha, const Self& a, cudnnReduceTensorOp_t op_type,
                    Scalar beta, /* in/out */ Self& b) {
                // Create the reduction operation descriptor.
                cudnnReduceTensorDescriptor_t desc;
                cudnnAssert(createReduceTensorDescriptor_t(&desc));
                cudnnAssert(cudnnSetReduceTensorDescriptor(desc, op_type, DATA_TYPE, NAN_PROP,
                        CUDNN_REDUCE_TENSOR_NO_INDICES, CUDNN_32BIT_INDICES));
                // Calculate the array size needed for the indices (should be 0).
                std::size_t indices_size;
                cudnnAssert(cudnnGetReductionIndicesSize(CuDNNHandle::get_instance(), desc, a.desc(),
                        b.desc(), &indices_size));
                Base indices(static_cast<std::size_t>(ceil(static_cast<Scalar>(indices_size) / sizeof(Scalar))));
                // Calculate the workspace size.
                std::size_t workspace_size;
                cudnnAssert(cudnnGetReductionWorkspaceSize(CuDNNHandle::get_instance(), desc, a.desc(),
                        b.desc(), &workspace_size));
                Base workspace(static_cast<std::size_t>(ceil(static_cast<Scalar>(workspace_size) / sizeof(Scalar))));
                // Perform the reduction.
                cudnnAssert(cudnnReduceTensor(CuDNNHandle::get_instance(), indices.data(), indices_size,
                        workspace.data(), workspace_size, &alpha, a.desc(), a.data(), &beta, b.desc(), b.data()));
                // Free resources.
                cudnnAssert(cudnnDestroyReduceTensorDescriptor(desc));
            }
        private:
            std::size_t _samples, _height, _width, _channels;
            cudnnTensorDescriptor_t _desc;
            cudnnFilterDescriptor_t _filter_desc;
        };

        /**
        * A class representing a cuRAND runtime error.
        */
        class CuRANDError : public std::runtime_error {
        public:
            /**
            * @param status The cuRAND status code.
            * @param file The name of the file in which the error occurred.
            * @param line The number of the line at which the error occurred.
            */
            CuRANDError(curandStatus_t status, const char* file, int line) :
                std::runtime_error("cuRAND Error: " + curand_status_to_string(status) + "; File: " +
                        std::string(file) + "; Line: " + std::to_string(line)) { }
        private:
            /**
            * It returns a string representation of the provided cuRAND status code.
            *
            * @param status The cuRAND status code.
            * @return A string describing the status code.
            */
            inline static std::string curand_status_to_string(curandStatus_t status) {
                switch (status) {
                    case CURAND_STATUS_SUCCESS: return "CURAND_STATUS_SUCCESS";
                    case CURAND_STATUS_VERSION_MISMATCH: return "CURAND_STATUS_VERSION_MISMATCH";
                    case CURAND_STATUS_NOT_INITIALIZED: return "CURAND_STATUS_NOT_INITIALIZED";
                    case CURAND_STATUS_ALLOCATION_FAILED: return "CURAND_STATUS_ALLOCATION_FAILED";
                    case CURAND_STATUS_TYPE_ERROR: return "CURAND_STATUS_TYPE_ERROR";
                    case CURAND_STATUS_OUT_OF_RANGE: return "CURAND_STATUS_OUT_OF_RANGE";
                    case CURAND_STATUS_LENGTH_NOT_MULTIPLE: return "CURAND_STATUS_LENGTH_NOT_MULTIPLE";
                    case CURAND_STATUS_DOUBLE_PRECISION_REQUIRED: return "CURAND_STATUS_DOUBLE_PRECISION_REQUIRED";
                    case CURAND_STATUS_LAUNCH_FAILURE: return "CURAND_STATUS_LAUNCH_FAILURE";
                    case CURAND_STATUS_PREEXISTING_FAILURE: return "CURAND_STATUS_PREEXISTING_FAILURE";
                    case CURAND_STATUS_INITIALIZATION_FAILED: return "CURAND_STATUS_INITIALIZATION_FAILED";
                    case CURAND_STATUS_ARCH_MISMATCH: return "CURAND_STATUS_ARCH_MISMATCH";
                    case CURAND_STATUS_INTERNAL_ERROR: return "CURAND_STATUS_INTERNAL_ERROR";
                    default: return "<unknown>";
                }
            }
        };

        namespace {

        __inline__ void _curand_error_check(curandStatus_t status, const char* file, int line) {
            if (status != CURAND_STATUS_SUCCESS)
                throw CuRANDError(status, file, line);
        }

        __inline__ void _curand_assert(curandStatus_t status, const char* file, int line) {
            try {
                _curand_error_check(status, file, line);
            } catch (const CuRANDError& e) {
                std::cout << e.what() << std::endl;
                exit(-1);
            }
        }

        }

        namespace {

        template<typename Scalar>
        using NormalGenerationRoutine = curandStatus_t (*)(curandGenerator_t, Scalar*, std::size_t, Scalar, Scalar);
        template<typename Scalar>
        __inline__ NormalGenerationRoutine<Scalar> normal_gen_routine() { return &curandGenerateNormalDouble; }
        template<> __inline__ NormalGenerationRoutine<float> normal_gen_routine() { return &curandGenerateNormal; }

        template<typename Scalar>
        using UniformGenerationRoutine = curandStatus_t (*)(curandGenerator_t, Scalar*, std::size_t);
        template<typename Scalar>
        __inline__ UniformGenerationRoutine<Scalar> uniform_gen_routine() { return &curandGenerateUniformDouble; }
        template<> __inline__ UniformGenerationRoutine<float> uniform_gen_routine() { return &curandGenerateUniform; }

        }

        /**
        * A template class for generating normally or uniformly distributed random numbers
        * using the cuRAND generator.
        */
        template<typename Scalar>
        class CuRANDGenerator {
        public:
            static constexpr curandRngType_t RNG_TYPE = CURAND_RNG_PSEUDO_DEFAULT;
            inline CuRANDGenerator() :
                    gen() {
                curandAssert(curandCreateGenerator(&gen, RNG_TYPE));
            }
            inline ~CuRANDGenerator() {
                curandAssert(curandDestroyGenerator(gen));
            }
            /**
            * @param mean The mean of the distribution.
            * @param sd The standard deviation of the distribution.
            * @param array The device array to fill with the randomly generated numbers.
            */
            inline void generate_normal(Scalar mean, Scalar sd, CUDAArray<Scalar>& array) const {
                NormalGenerationRoutine<Scalar> norm_gen = normal_gen_routine<Scalar>();
                curandAssert(norm_gen(gen, array.data(), array.size() * sizeof(Scalar), mean, sd));
            }
            /**
            * @param array The device array to fill with the randomly generated numbers.
            */
            inline void generate_uniform(CUDAArray<Scalar>& array) const {
                UniformGenerationRoutine<Scalar> uni_gen = uniform_gen_routine<Scalar>();
                curandAssert(uni_gen(gen, array.data(), array.size() * sizeof(Scalar)));
            }
        private:
            curandGenerator_t gen;
        };


        template<typename Scalar>
        class GPUParameters : public virtual Parameters<Scalar> {
        public:
            virtual ~GPUParameters() = default;
            virtual GPUParameters<Scalar>* gpu_clone() const = 0;
            virtual std::size_t samples() const = 0;
            virtual std::size_t height() const = 0;
            virtual std::size_t width() const = 0;
            virtual std::size_t channels() const = 0;
            virtual const CuDNNTensor<Scalar>& get_gpu_values() const = 0;
            virtual void set_gpu_values(CuDNNTensor<Scalar> values) = 0;
            virtual const CuDNNTensor<Scalar>& get_gpu_grad() const = 0;
            virtual void accumulate_gpu_grad(const CuDNNTensor<Scalar>& grad) = 0;
            inline Parameters<Scalar>* clone() const {
                return gpu_clone();
            }
        };

        template<typename Scalar, std::size_t Rank>
        class GPULayer : public virtual Layer<Scalar,Rank> {
            typedef Layer<Scalar,Rank> Base;
        protected:
            typedef Dimensions<std::size_t,3> GPUDims;
        public:
            virtual ~GPULayer() = default;
            virtual GPULayer<Scalar,Rank>* gpu_clone() const = 0;
            virtual GPULayer<Scalar,Rank>* gpu_clone_with_shared_params() = 0;
            virtual const GPUDims& get_gpu_input_dims() const = 0;
            virtual const GPUDims& get_gpu_output_dims() const = 0;
            virtual std::vector<const GPUParameters<Scalar>*> get_gpu_params() const = 0;
            virtual std::vector<GPUParameters<Scalar>*> get_gpu_params() = 0;
            virtual CuDNNTensor<Scalar> pass_forward(CuDNNTensor<Scalar> in, bool training) = 0;
            virtual CuDNNTensor<Scalar> pass_back(CuDNNTensor<Scalar> out_grad) = 0;
            inline Layer<Scalar,Rank>* clone() const {
                return gpu_clone();
            }
            inline Layer<Scalar,Rank>* clone_with_shared_params() {
                return gpu_clone_with_shared_params();
            }
            inline typename Base::Data pass_forward(typename Base::Data in, bool training) {
                auto rows = in.dimension(0);
                auto in_gpu_dims = get_gpu_input_dims().template extend<>();
                in_gpu_dims(0) = rows;
                Tensor<Scalar,4> out_extended = pass_forward(CuDNNTensor<Scalar>(TensorMap<Scalar,4>(in.data(),
                        in_gpu_dims)), training);
                auto out_dims = get_output_dims().template extend<>();
                out_dims(0) = rows;
                return TensorMap<Scalar,Base::DATA_RANK>(out_extended.data(), out_dims);
            }
            inline typename Base::Data pass_back(typename Base::Data out_grad) {
                auto rows = out_grad.dimension(0);
                auto out_gpu_dims = get_gpu_output_dims().template extend<>();
                out_gpu_dims(0) = rows;
                Tensor<Scalar,4> prev_out_grad_extended = pass_back(TensorMap<Scalar,4>(out_grad.data(),
                        out_gpu_dims));
                auto in_dims = get_input_dims().template extend<>();
                in_dims(0) = rows;
                return TensorMap<Scalar,Base::DATA_RANK>(prev_out_grad_extended.data(), in_dims);
            }
        };

        template<typename Scalar, std::size_t Rank, bool Sequential>
        class GPUNeuralNetwork : public virtual NeuralNetwork<Scalar,Rank,Sequential> {
            typedef NeuralNetwork<Scalar,Rank,Sequential> Base;
        protected:
            typedef Dimensions<std::size_t,3> GPUDims;
        public:
            virtual const GPUDims& get_gpu_input_dims() const = 0;
            virtual const GPUDims& get_gpu_output_dims() const = 0;
            virtual std::vector<const GPULayer<Scalar,Rank>*> get_gpu_layers() const = 0;
            virtual std::vector<GPULayer<Scalar,Rank>*> get_gpu_layers() = 0;
            virtual CuDNNTensor<Scalar> propagate(CuDNNTensor<Scalar> input, bool training) = 0;
            virtual CuDNNTensor<Scalar> backpropagate(CuDNNTensor<Scalar> out_grad) = 0;
            inline typename Base::Data propagate(typename Base::Data input, bool training) {
                auto rows = input.dimension(0);
                auto in_gpu_dims = get_gpu_input_dims().template extend<>();
                in_gpu_dims(0) = rows;
                Tensor<Scalar,4> out_extended = propagate(CuDNNTensor<Scalar>(TensorMap<Scalar,4>(input.data(),
                        in_gpu_dims)), training);
                auto out_dims = get_output_dims().template extend<>();
                out_dims(0) = rows;
                return TensorMap<Scalar,Base::DATA_RANK>(out_extended.data(), out_dims);
            }
            inline typename Base::Data backpropagate(typename Base::Data out_grad) {
                auto rows = out_grad.dimension(0);
                auto out_gpu_dims = get_gpu_output_dims().template extend<>();
                out_gpu_dims(0) = rows;
                Tensor<Scalar,4> prev_out_grad_extended = backpropagate(CuDNNTensor<Scalar>(
                        TensorMap<Scalar,4>(out_grad.data(), out_gpu_dims)));
                auto in_dims = get_input_dims().template extend<>();
                in_dims(0) = rows;
                return TensorMap<Scalar,Base::DATA_RANK>(prev_out_grad_extended.data(), in_dims);
            }
        };        

        template<typename Scalar>
        class GPUParameterInitialization : public virtual ParameterInitialization<Scalar> {
        public:
            virtual void apply(CuBLASMatrix<Scalar>& params) const = 0;
            inline void apply(Matrix<Scalar>& params) const {
                CuBLASMatrix<Scalar> gpu_params = params;
                apply(gpu_params);
                gpu_params.copy_to_host(params.data());
            }
        };

        template<typename Scalar>
        class GPUParameterRegularization : public virtual ParameterRegularization<Scalar> {
        public:
            virtual Scalar function(const CuBLASMatrix<Scalar>& params) const;
            virtual CuBLASMatrix<Scalar> d_function(const CuBLASMatrix<Scalar>& params) const;
            inline Scalar function(const Matrix<Scalar>& params) const {
                return function(CuBLASMatrix<Scalar>(params));
            }
            inline Matrix<Scalar> d_function(const Matrix<Scalar>& params) const {
                return d_function(CuBLASMatrix<Scalar>(params));
            }
        };

        template<typename Scalar>
        class ConstantGPUParameterInitialization : public GPUParameterInitialization<Scalar> {
        public:
            /**
            * @param constant The value to which all elements of the parameter matrix are to be
            * initialized.
            */
            ConstantGPUParameterInitialization(Scalar constant) :
                    constant(constant) { }
            inline void apply(CuBLASMatrix<Scalar>& params) const {
                CuDNNTensor<Scalar> cudnn_params(params.data(), params.cols(), params.rows(), 1u, 1u);
                cudnn_params.set_values(constant);
            }
        private:
            Scalar constant;
        };

        template<typename Scalar>
        class GaussianGPUParameterInitialization : public GPUParameterInitialization<Scalar> {
        public:
            /**
            * @param mean The mean of the distribution.
            * @param sd_scaling_factor The standard deviation scaling factor.
            */
            GaussianGPUParameterInitialization(Scalar mean = 0, Scalar sd_scaling_factor = 1) :
                    mean(mean),
                    sd_scaling_factor(sd_scaling_factor) {
                assert(sd_scaling_factor > 0);
            }
            inline void apply(CuBLASMatrix<Scalar>& params) const {
                CuRANDGenerator<Scalar> gen;
                gen.generate_normal(0, sd_scaling_factor * _sd(params.rows(), params.cols()), params);
            }
        protected:
            /**
            * It computes the standard deviation of the distribution to sample from.
            *
            * @param fan_ins The input size of the kernel (the number of rows in the weight matrix
            * excluding the bias row).
            * @param fan_outs The output size of the kernel (the number of columns in the weight
            * matrix).
            * @return The standard deviation of the normal distribution from which the values
            * of the initialized weight matrix are to be sampled.
            */
            inline virtual Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
                return 1;
            }
        private:
            const Scalar mean, sd_scaling_factor;
        };        

        template<typename Scalar>
        class GlorotGPUParameterInitialization : public GaussianGPUParameterInitialization<Scalar> {
        public:
            /**
            * @param sd_scaling_factor The standard deviation scaling factor.
            */
            GlorotGPUParameterInitialization(Scalar sd_scaling_factor = 1) :
                    GaussianGPUParameterInitialization<Scalar>(0, sd_scaling_factor) { }
        protected:
            inline Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
                return sqrt(2 / (Scalar) (fan_ins + fan_outs));
            }
        };      


        template<typename Scalar>
        class HeGPUParameterInitialization : public GaussianGPUParameterInitialization<Scalar> {
        public:
            /**
            * @param sd_scaling_factor The standard deviation scaling factor.
            */
            HeGPUParameterInitialization(Scalar sd_scaling_factor = 1) :
                    GaussianGPUParameterInitialization<Scalar>(0, sd_scaling_factor) { }
        protected:
            inline Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
                return sqrt(2 / (Scalar) fan_ins);
            }
        };  

        template<typename Scalar>
        class LeCunGPUParameterInitialization : public GaussianGPUParameterInitialization<Scalar> {
        public:
            /**
            * @param sd_scaling_factor The standard deviation scaling factor.
            */
            LeCunGPUParameterInitialization(Scalar sd_scaling_factor = 1) :
                    GaussianGPUParameterInitialization<Scalar>(0, sd_scaling_factor) { }
        protected:
            inline Scalar _sd(unsigned fan_ins, unsigned fan_outs) const {
                return sqrt(1 / (Scalar) fan_ins);
            }
        };

        template<typename Scalar>
        class OneGPUParameterInitialization : public GPUParameterInitialization<Scalar> {
        public:
            inline void apply(CuBLASMatrix<Scalar>& params) const {
                params.set_values(1);
            }
        };

        template<typename Scalar>
        class ZeroGPUParameterInitialization : public GPUParameterInitialization<Scalar> {
        public:
            inline void apply(CuBLASMatrix<Scalar>& params) const {
                params.set_values(0);
            }
        };

    } /* namespace gpu */
} /* namespace cattle */