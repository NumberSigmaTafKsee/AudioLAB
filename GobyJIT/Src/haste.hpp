#pragma once

// GENERAL NOTES:
// No pointers may be null unless otherwise specified.
// All pointers are expected to point to device memory.
// The square brackets below describe tensor shapes, where
//     T = number of RNN time steps
//     N = batch size
//     C = input size
//     H = hidden size
// and the rightmost dimension changes the fastest.


#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <cuda_fp16.h>

template<typename T>
struct blas {
  struct set_pointer_mode {
    set_pointer_mode(cublasHandle_t handle) : handle_(handle) {
      cublasGetPointerMode(handle_, &old_mode_);
      cublasSetPointerMode(handle_, CUBLAS_POINTER_MODE_HOST);
    }
    ~set_pointer_mode() {
      cublasSetPointerMode(handle_, old_mode_);
    }
    private:
      cublasHandle_t handle_;
      cublasPointerMode_t old_mode_;
  };
  struct enable_tensor_cores {
    enable_tensor_cores(cublasHandle_t handle) : handle_(handle) {
      cublasGetMathMode(handle_, &old_mode_);
      cublasSetMathMode(handle_, CUBLAS_TENSOR_OP_MATH);
    }
    ~enable_tensor_cores() {
      cublasSetMathMode(handle_, old_mode_);
    }
    private:
      cublasHandle_t handle_;
      cublasMath_t old_mode_;
  };
};

template<>
struct blas<__half> {
  static constexpr decltype(cublasHgemm)* gemm = &cublasHgemm;
};

template<>
struct blas<float> {
  static constexpr decltype(cublasSgemm)* gemm = &cublasSgemm;
};

template<>
struct blas<double> {
  static constexpr decltype(cublasDgemm)* gemm = &cublasDgemm;
};


namespace haste {
namespace v0 {
namespace gru {

template<typename T>
class ForwardPass {
  public:
    // training: `true` if the caller intends to perform a backward pass to compute gradients.
    // batch_size: the number of training/inference inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    ForwardPass(
        const bool training,
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~ForwardPass();

    // Performs one forward iteration of the GRU cell.
    //
    // W: [C,H*3] the input weight matrix.
    // R: [H,H*3] the recurrent weight matrix.
    // bx: [H*3] the bias for the input weight matrix.
    // br: [H*3] the bias for the recurrent weight matrix.
    // x: [N,C] the GRU input for this iteration (N vectors, each with dimension C).
    // h: [N,H] the t-1 iteration's `h_out` or the initial hidden state if this is the
    //     t=0 iteration (typically zeros).
    // h_out: [N,H] the GRU's output, and the input to the next iteration's `h`. This
    //     pointer may be the same as `h`. Each iteration may reuse the same memory region.
    // v: [N,H*4] if `training` is `false`, this can be a null pointer. If `training` is
    //     `true`, this vector will contain intermediate activations for this iteration which
    //     must be provided as-is to the corresponding backward iteration. The caller must
    //     provide a new memory region for each iteration.
    // tmp_Wx: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector, and must provide a new memory region for
    //     each iteration.
    // tmp_Rh: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. The same memory region may be provided
    //     for each iteration.
    // zoneout_prob: 0.0 <= zoneout_prob <= 1.0; specifies the probability of a hidden
    //     activation being randomly zoned out. If zoneout was used during training, this
    //     parameter must also be specified during inference with the same value.
    // zoneout_mask: [N,H] may be null to disable zoneout. This is a random binary mask
    //     following a Bernoulli(1-zoneout_prob) distribution. A different mask is typically
    //     used for each iteration.
    void Iterate(
        const T* W,
        const T* R,
        const T* bx,
        const T* br,
        const T* x,
        const T* h,
        T* h_out,
        T* v,
        T* tmp_Wx,
        T* tmp_Rh,
        const float zoneout_prob,
        const T* zoneout_mask);

    void Run(
        const int steps,
        const T* W,
        const T* R,
        const T* bx,
        const T* br,
        const T* x,
        T* h,
        T* v,
        T* tmp_Wx,
        T* tmp_Rh,
        const float zoneout_prob,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R,
        const T* bx,
        const T* br,
        const T* h,
        T* h_out,
        T* v,
        T* tmp_Wx,
        T* tmp_Rh,
        const float zoneout_prob,
        const T* zoneout_mask);

    struct private_data;
    private_data* data_;
};

template<typename T>
class BackwardPass {
  public:
    // batch_size: the number of training inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    BackwardPass(
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~BackwardPass();

    // Performs one backward iteration of the GRU cell.
    //
    // Note that BackwardPass must be iterated in the reverse order as ForwardPass.
    // If ForwardPass iterates from 0 to T-1, BackwardPass needs to iterate from
    // T-1 down to 0. When iteration numbers are described, they will be based on the
    // iteration index (i.e., the T-1'th iteration of the forward pass is the last call
    // to ForwardPass::Iterate, whereas it is the first call to BackwardPass::Iterate).
    //
    // W_t: [H*3,C] the transpose of the input weight matrix.
    // R_t: [H*3,H] the transpose of the recurrent weight matrix.
    // bx: [H*3] the bias vector for the input weight matrix.
    // br: [H*3] the bias vector for the recurrent weight matrix.
    // x_t: [C,N] the transpose of the GRU input for this iteration.
    // h: [N,H] the t-1 iteration's `h_out` or the initial hidden state if this is the t=0
    //     iteration (typically zeros).
    // v: [N,H*4] the same vector as returned by ForwardPass::Iterate on its corresponding
    //     iteration.
    // dh_new: [N,H] the gradient of `h_out` with respect to the loss at this iteration.
    // dx: [N,C] the gradient of the input at this time step with respect to the loss.
    // dW: [C,H*3] the gradient of the input weight matrix with respect to the loss.
    // dR: [H,H*3] the gradient of the recurrent weight matrix with respect to the loss.
    // dbx: [H*3] the gradient of the bias vector for the input weight matrix with respect to
    //     the loss.
    // dbr: [H*3] the gradient of the bias vector for the recurrent weight matrix with respect
    //     to the loss.
    // dh: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros
    //     for the T-1'th iteration and the same pointer should be passed in for each
    //     iteration. After a complete backward pass, this vector will contain the gradient
    //     of the initial hidden state with respect to the loss.
    // dp: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. A new memory region must be provided
    //     for each iteration.
    // dq: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. A new memory region must be provided
    //     for each iteration.
    // zoneout_mask: [N,H] may be null if zoneout was disabled in the forward pass. This vector
    //     must be the same as the one provided during the corresponding forward iteration.
    void Iterate(
        const T* W_t,
        const T* R_t,
        const T* bx,
        const T* br,
        const T* x_t,
        const T* h,
        const T* v,
        const T* dh_new,
        T* dx,
        T* dW,
        T* dR,
        T* dbx,
        T* dbr,
        T* dh,
        T* dp,
        T* dq,
        const T* zoneout_mask);

    void Run(
        const int steps,
        const T* W_t,
        const T* R_t,
        const T* bx,
        const T* br,
        const T* x_t,
        const T* h,
        const T* v,
        const T* dh_new,
        T* dx,
        T* dW,
        T* dR,
        T* dbx,
        T* dbr,
        T* dh,
        T* dp,
        T* dq,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R_t,
        const T* h,
        const T* v,
        const T* dh_new,
        T* dbx,
        T* dbr,
        T* dh,
        T* dp,
        T* dq,
        const T* zoneout_mask);

    struct private_data;
    private_data* data_;
};

}  // namespace gru
}  // namespace v0
}  // namespace haste


namespace haste {
namespace v0 {
namespace indrnn {

template<typename T>
class ForwardPass {
  public:
    ForwardPass(
        const bool training,
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    ~ForwardPass();

    void Run(
        const int steps,
        const T* W,
        const T* u,
        const T* b,
        const T* x,
        T* h,
        T* workspace,
        const float zoneout_prob,
        const T* zoneout_mask);

  private:
    struct private_data;
    private_data* data_;
};

template<typename T>
class BackwardPass {
  public:
    BackwardPass(
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    ~BackwardPass();

    void Run(
        const int steps,
        const T* W_t,
        const T* u,
        const T* b,
        const T* x_t,
        const T* h,
        const T* dh_new,
        T* dx,
        T* dW,
        T* du,
        T* db,
        T* dh,
        T* workspace,
        const T* zoneout_mask);

  private:
    struct private_data;
    private_data* data_;
};

}  // namespace indrnn
}  // namespace v0
}  // namespace haste

namespace haste {
namespace v0 {
namespace layer_norm {

template<typename T>
class ForwardPass {
  public:
    // gamma: [H]
    // beta: [H]
    // cache: [N,2]
    ForwardPass(
        const int batch_size,
        const int hidden_size,
        const T* gamma,
        const T* beta,
        T* cache);

    // Computes the layer norm of an input tensor `x` over its innermost (fastest changing)
    // dimension. The layer norm is defined as: \(\frac{x-\mu}{\sigma} \gamma + \beta\)
    // where `\gamma` and `\beta` are trainable parameters.
    //
    // x: [N,H]
    // y: [N,H]
    void Run(const cudaStream_t& stream, const T* x, T* y);

    void RunPartial(
        const cudaStream_t& stream,
        const int minibatch,
        const T* x,
        T* y);

  private:
    const int batch_size_;
    const int hidden_size_;
    const T* gamma_;
    const T* beta_;
    T* cache_;
    int partial_;
};

template<typename T>
class BackwardPass {
  public:
    BackwardPass(
        const int batch_size,
        const int hidden_size,
        const T* gamma,
        const T* beta,
        const T* x,
        T* dgamma,
        T* dbeta,
        T* cache);

    void Run(const cudaStream_t& stream, const T* dy, T* dx);

    void RunPartial(
        const cudaStream_t& stream,
        const int minibatch,
        const T* dy,
        T* dx);

  private:
    const int batch_size_;
    const int hidden_size_;
    const T* gamma_;
    const T* beta_;
    const T* x_;
    T* dgamma_;
    T* dbeta_;
    T* cache_;
    int partial_;
};

}  // namespace layer_norm
}  // namespace v0
}  // namespace haste

namespace haste {
namespace v0 {
namespace layer_norm_gru {

template<typename T>
class ForwardPass {
  public:
    // training: `true` if the caller intends to perform a backward pass to compute gradients.
    // batch_size: the number of training/inference inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    ForwardPass(
        const bool training,
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~ForwardPass();

    // Performs one forward iteration of the GRU cell.
    //
    // W: [C,H*3] the input weight matrix.
    // R: [H,H*3] the recurrent weight matrix.
    // bx: [H*3] the bias for the input weight matrix.
    // br: [H*3] the bias for the recurrent weight matrix.
    // x: [N,C] the GRU input for this iteration (N vectors, each with dimension C).
    // h: [N,H] the t-1 iteration's `h_out` or the initial hidden state if this is the
    //     t=0 iteration (typically zeros).
    // h_out: [N,H] the GRU's output, and the input to the next iteration's `h`. This
    //     pointer may be the same as `h`. Each iteration may reuse the same memory region.
    // v: [N,H*4] if `training` is `false`, this can be a null pointer. If `training` is
    //     `true`, this vector will contain intermediate activations for this iteration which
    //     must be provided as-is to the corresponding backward iteration. The caller must
    //     provide a new memory region for each iteration.
    // tmp_Wx: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector, and must provide a new memory region for
    //     each iteration.
    // tmp_Rh: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. The same memory region may be provided
    //     for each iteration.
    // zoneout_prob: 0.0 <= zoneout_prob <= 1.0; specifies the probability of a hidden
    //     activation being randomly zoned out. If zoneout was used during training, this
    //     parameter must also be specified during inference with the same value.
    // zoneout_mask: [N,H] may be null to disable zoneout. This is a random binary mask
    //     following a Bernoulli(1-zoneout_prob) distribution. A different mask is typically
    //     used for each iteration.
    void Run(
        const int steps,
        const T* W,
        const T* R,
        const T* bx,
        const T* br,
        const T* x,
        T* h,
        T* v,
        T* act_Wx,
        layer_norm::ForwardPass<T>& layer_norm1,
        T* tmp_Wx_norm,
        T* act_Rh,
        layer_norm::ForwardPass<T>& layer_norm2,
        T* tmp_Rh_norm,
        const float zoneout_prob,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R,
        const T* bx,
        const T* br,
        const T* h,
        T* h_out,
        T* v,
        T* tmp_Wx_norm,
        T* act_Rh,
        layer_norm::ForwardPass<T>& layer_norm2,
        T* tmp_Rh_norm,
        const float zoneout_prob,
        const T* zoneout_mask);

    struct private_data;
    private_data* data_;
};

template<typename T>
class BackwardPass {
  public:
    // batch_size: the number of training inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    BackwardPass(
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~BackwardPass();

    // Performs one backward iteration of the GRU cell.
    //
    // Note that BackwardPass must be iterated in the reverse order as ForwardPass.
    // If ForwardPass iterates from 0 to T-1, BackwardPass needs to iterate from
    // T-1 down to 0. When iteration numbers are described, they will be based on the
    // iteration index (i.e., the T-1'th iteration of the forward pass is the last call
    // to ForwardPass::Iterate, whereas it is the first call to BackwardPass::Iterate).
    //
    // W_t: [H*3,C] the transpose of the input weight matrix.
    // R_t: [H*3,H] the transpose of the recurrent weight matrix.
    // bx: [H*3] the bias vector for the input weight matrix.
    // br: [H*3] the bias vector for the recurrent weight matrix.
    // x_t: [C,N] the transpose of the GRU input for this iteration.
    // h: [N,H] the t-1 iteration's `h_out` or the initial hidden state if this is the t=0
    //     iteration (typically zeros).
    // v: [N,H*4] the same vector as returned by ForwardPass::Iterate on its corresponding
    //     iteration.
    // dh_new: [N,H] the gradient of `h_out` with respect to the loss at this iteration.
    // dx: [N,C] the gradient of the input at this time step with respect to the loss.
    // dW: [C,H*3] the gradient of the input weight matrix with respect to the loss.
    // dR: [H,H*3] the gradient of the recurrent weight matrix with respect to the loss.
    // dbx: [H*3] the gradient of the bias vector for the input weight matrix with respect to
    //     the loss.
    // dbr: [H*3] the gradient of the bias vector for the recurrent weight matrix with respect
    //     to the loss.
    // dh: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros
    //     for the T-1'th iteration and the same pointer should be passed in for each
    //     iteration. After a complete backward pass, this vector will contain the gradient
    //     of the initial hidden state with respect to the loss.
    // dp: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. A new memory region must be provided
    //     for each iteration.
    // dq: [N,H*3] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. A new memory region must be provided
    //     for each iteration.
    // zoneout_mask: [N,H] may be null if zoneout was disabled in the forward pass. This vector
    //     must be the same as the one provided during the corresponding forward iteration.
    void Run(
        const int steps,
        const T* W_t,
        const T* R_t,
        const T* bx,
        const T* br,
        const T* x_t,
        const T* h,
        const T* v,
        const T* dh_new,
        T* dx,
        T* dW,
        T* dR,
        T* dbx,
        T* dbr,
        T* dh,
        T* dp,
        T* dq,
        layer_norm::BackwardPass<T>& layer_norm1,
        layer_norm::BackwardPass<T>& layer_norm2,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R_t,
        const T* h,
        const T* v,
        const T* dh_new,
        T* dbx,
        T* dbr,
        T* dh,
        T* dp,
        T* dq,
        layer_norm::BackwardPass<T>& layer_norm2,
        const T* zoneout_mask);

    struct private_data;
    private_data* data_;
};

}  // namespace layer_norm_gru
}  // namespace v0
}  // namespace haste

namespace haste {
namespace v0 {
namespace layer_norm_indrnn {

template<typename T>
class ForwardPass {
  public:
    ForwardPass(
        const bool training,
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    ~ForwardPass();

    void Run(
        const int steps,
        const T* W,
        const T* u,
        const T* b,
        const T* x,
        T* h,
        T* workspace,
        T* act_Wx,
        layer_norm::ForwardPass<T>& layer_norm1,
        const float zoneout_prob,
        const T* zoneout_mask);

  private:
    struct private_data;
    private_data* data_;
};

template<typename T>
class BackwardPass {
  public:
    BackwardPass(
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    ~BackwardPass();

    void Run(
        const int steps,
        const T* W_t,
        const T* u,
        const T* b,
        const T* x_t,
        const T* h,
        const T* dh_new,
        T* dx,
        T* dW,
        T* du,
        T* db,
        T* dh,
        T* workspace,
        layer_norm::BackwardPass<T>& layer_norm1,
        const T* zoneout_mask);

  private:
    struct private_data;
    private_data* data_;
};

}  // namespace layer_norm_indrnn
}  // namespace v0
}  // namespace haste

namespace haste {
namespace v0 {
namespace layer_norm_lstm {

template<typename T>
class ForwardPass {
  public:
    // training: `true` if the caller intends to perform a backward pass to compute gradients.
    // batch_size: the number of training/inference inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    ForwardPass(
        const bool training,
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~ForwardPass();

    // Runs the LSTM over all time steps. This method is faster than using a per-step
    // `Iterate` but requires that the entire input sequence be available upfront. In some
    // situations, this constraint may not be satisfiable (e.g. autoregressive models).
    // Users should prefer calling `Run` over `Iterate` whenever possible.
    //
    // steps: the number of iterations to run (i.e. T).
    // W: [C,H*4] the input weight matrix.
    // R: [H,H*4] the recurrent weight matrix.
    // b: [H*4] the bias vector.
    // x: [T,N,C] the LSTM input for this iteration (N vectors, each with dimension C).
    // h: [T+1,N,H] the hidden state vectors across all time steps. The t=0'th vector should
    //      be set to the desired initial hidden state (typically zeros). The rest of the
    //      vectors will be set by this function. `h[1:,:,:]` forms the output of this LSTM
    //      layer.
    // c: [T+1,N,H] the cell state vectors across all time steps. The t=0'th vector should be
    //      set to the desired initial cell state (typically zeros). The rest of the vectors
    //      will be set by this function.
    // v: [T,N,H*4] if `training` is `false`, this is scratch space and should not be used by
    //     the caller. If `training` is `true`, this parameter will contain intermediate
    //     activations which must be provided as-is to `BackwardPass::Run` or manually urolled
    //     for `BackwardPass::Iterate`.
    // tmp_Rh: [N,H*4] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. The same memory region may be provided
    //     for each iteration.
    // zoneout_prob: 0.0 <= zoneout_prob <= 1.0; specifies the probability of a hidden
    //     activation being randomly zoned out. If zoneout was used during training, this
    //     parameter must also be specified during inference with the same value.
    // zoneout_mask: [T,N,H] may be null to disable zoneout. This is a random binary mask
    //     following a Bernoulli(1-zoneout_prob) distribution. A different mask is typically
    //     used for each iteration.
    void Run(
        const int steps,
        const T* W,
        const T* R,
        const T* b,
        const T* x,
        T* h,
        T* c,
        T* act_Wx,
        T* tmp_Rh,
        layer_norm::ForwardPass<T>& layer_norm1,
        T* act_Wx_norm,
        T* act_Rh,
        layer_norm::ForwardPass<T>& layer_norm2,
        layer_norm::ForwardPass<T>& layer_norm3,
        T* act_c_norm,
        const float zoneout_prob,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R,
        const T* b,
        const T* h,
        const T* c,
        T* h_out,
        T* c_out,
        T* v,
        T* tmp_Rh,
        T* act_Rh,
        layer_norm::ForwardPass<T>& layer_norm2,
        layer_norm::ForwardPass<T>& layer_norm3,
        T* act_c_norm,
        const float zoneout_prob,
        const T* zoneout_mask);

    struct private_data;
    private_data* data_;
};

template<typename T>
class BackwardPass {
  public:
    // batch_size: the number of training inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    BackwardPass(
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~BackwardPass();

    // Runs the LSTM backward pass over all time steps. This method is faster than using a
    // per-step `Iterate` but requires that the entire input sequence be available upfront.
    // In some situations, this constraint may not be satisfiable (e.g. autoregressive models).
    // Users should prefer calling `Run` over `Iterate` whenever possible.
    //
    // steps: the number of iterations to run (i.e. T).
    // W_t: [H*4,C] the transpose of the input weight matrix.
    // R_t: [H*4,H] the transpose of the recurrent weight matrix.
    // b: [H*4] the bias vector.
    // x_t: [C,T,N] the transpose of the LSTM input for this iteration.
    // h: [T+1,N,H] the hidden state vectors after running `ForwardPass::Run`.
    // c: [T+1,N,H] the cell state vectors after running `ForwardPass::Run`.
    // dh_new: [T+1,N,H] the gradient of the loss with respect to `h`.
    // dc_new: [T+1,N,H] the gradient of the loss with respect to `c`.
    // dx: [T,N,C] the gradient of the loss with respect to the input.
    // dW: [C,H*4] the gradient of the loss with respect to the input weight matrix.
    // dR: [H,H*4] the gradient of the loss with respect to the recurrent weight matrix.
    // db: [H*4] the gradient of the loss with respect to the bias vector.
    // dh: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros.
    //     When this function returns, `dh` will contain the gradient of the loss with respect
    //     to the initial hidden state.
    // dc: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros.
    //     When this function returns, `dc` will contain the gradient of the loss with respect
    //     to the initial cell state.
    // v: [T,N,H*4] the same tensor that was passed to `ForwardPass::Run`.
    // zoneout_mask: [T,N,H] may be null if zoneout was disabled in the forward pass. This
    //     vector must be the same as the one provided during the forward pass.
    void Run(
        const int steps,
        const T* W_t,
        const T* R_t,
        const T* b,
        const T* x_t,
        const T* h,
        const T* c,
        const T* dh_new,
        const T* dc_new,
        T* dx,
        T* dW,
        T* dR,
        T* db,
        T* dh,
        T* dc,
        T* act_Wx,
        layer_norm::BackwardPass<T>& layer_norm1,
        T* act_Wx_norm,
        T* act_Rh,
        layer_norm::BackwardPass<T>& layer_norm2,
        layer_norm::BackwardPass<T>& layer_norm3,
        T* act_c_norm,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R_t,
        const T* c,
        const T* c_new,
        const T* dh_new,
        const T* dc_new,
        T* db,
        T* dh,
        T* dc,
        T* v,
        T* act_Rh,
        layer_norm::BackwardPass<T>& layer_norm2,
        layer_norm::BackwardPass<T>& layer_norm3,
        T* act_c_norm,
        const T* zoneout_mask);
    struct private_data;
    private_data* data_;
};

}  // namespace layer_norm_lstm
}  // namespace v0
}  // namespace haste

namespace haste {
namespace v0 {
namespace lstm {

template<typename T>
class ForwardPass {
  public:
    // training: `true` if the caller intends to perform a backward pass to compute gradients.
    // batch_size: the number of training/inference inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    ForwardPass(
        const bool training,
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~ForwardPass();

    // Performs one forward iteration of the LSTM cell.
    //
    // W: [C,H*4] the input weight matrix.
    // R: [H,H*4] the recurrent weight matrix.
    // b: [H*4] the bias vector.
    // x: [N,C] the LSTM input for this iteration (N vectors, each with dimension C).
    // h: [N,H] the t-1 iteration's `h_out` or the initial hidden state if this is the
    //     t=0 iteration (typically zeros).
    // c: [N,H] the t-1 iteration's `c_out` or the initial cell state if this is the
    //     t=0 iteration (typically zeros).
    // h_out: [N,H] the LSTM's output, and the input to the next iteration's `h`. This
    //     pointer may be the same as `h`. Each iteration may reuse the same memory region.
    // c_out: [N,H] the LSTM's internal cell state after this iteration is complete. This
    //     will become the input to the next iteration's `c`.
    // v: [N,H*4] if `training` is `false`, this is scratch space and should not be used by
    //     the caller. If `training` is `true`, this vector will contain intermediate
    //     activations for this iteration which must be provided as-is to the corresponding
    //     backward iteration. In either case, a new memory region must be provided for each
    //     iteration.
    // tmp_Rh: [N,H*4] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. The same memory region may be provided
    //     for each iteration.
    // zoneout_prob: 0.0 <= zoneout_prob <= 1.0; specifies the probability of a hidden
    //     activation being randomly zoned out. If zoneout was used during training, this
    //     parameter must also be specified during inference with the same value.
    // zoneout_mask: [N,H] may be null to disable zoneout. This is a random binary mask
    //     following a Bernoulli(1-zoneout_prob) distribution. A different mask is typically
    //     used for each iteration.
    void Iterate(
        const cudaStream_t& stream,
        const T* W,
        const T* R,
        const T* b,
        const T* x,
        const T* h,
        const T* c,
        T* h_out,
        T* c_out,
        T* v,
        T* tmp_Rh,
        const float zoneout_prob,
        const T* zoneout_mask);

    // Runs the LSTM over all time steps. This method is faster than using a per-step
    // `Iterate` but requires that the entire input sequence be available upfront. In some
    // situations, this constraint may not be satisfiable (e.g. autoregressive models).
    // Users should prefer calling `Run` over `Iterate` whenever possible.
    //
    // steps: the number of iterations to run (i.e. T).
    // W: [C,H*4] the input weight matrix.
    // R: [H,H*4] the recurrent weight matrix.
    // b: [H*4] the bias vector.
    // x: [T,N,C] the LSTM input for this iteration (N vectors, each with dimension C).
    // h: [T+1,N,H] the hidden state vectors across all time steps. The t=0'th vector should
    //      be set to the desired initial hidden state (typically zeros). The rest of the
    //      vectors will be set by this function. `h[1:,:,:]` forms the output of this LSTM
    //      layer.
    // c: [T+1,N,H] the cell state vectors across all time steps. The t=0'th vector should be
    //      set to the desired initial cell state (typically zeros). The rest of the vectors
    //      will be set by this function.
    // v: [T,N,H*4] if `training` is `false`, this is scratch space and should not be used by
    //     the caller. If `training` is `true`, this parameter will contain intermediate
    //     activations which must be provided as-is to `BackwardPass::Run` or manually urolled
    //     for `BackwardPass::Iterate`.
    // tmp_Rh: [N,H*4] additional temporary work space required for this iteration. The caller
    //     should not use the contents of this vector. The same memory region may be provided
    //     for each iteration.
    // zoneout_prob: 0.0 <= zoneout_prob <= 1.0; specifies the probability of a hidden
    //     activation being randomly zoned out. If zoneout was used during training, this
    //     parameter must also be specified during inference with the same value.
    // zoneout_mask: [T,N,H] may be null to disable zoneout. This is a random binary mask
    //     following a Bernoulli(1-zoneout_prob) distribution. A different mask is typically
    //     used for each iteration.
    void Run(
        const int steps,
        const T* W,
        const T* R,
        const T* b,
        const T* x,
        T* h,
        T* c,
        T* v,
        T* tmp_Rh,
        const float zoneout_prob,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R,
        const T* b,
        const T* h,
        const T* c,
        T* h_out,
        T* c_out,
        T* v,
        T* tmp_Rh,
        const float zoneout_prob,
        const T* zoneout_mask);

    struct private_data;
    private_data* data_;
};

template<typename T>
class BackwardPass {
  public:
    // batch_size: the number of training inputs provided in each tensor.
    // input_size: the dimension of each input vector.
    // hidden_size: the expected dimension of each output vector.
    // blas_handle: an initialized cuBLAS handle (see `cublasCreate`).
    BackwardPass(
        const int batch_size,
        const int input_size,
        const int hidden_size,
        const cublasHandle_t& blas_handle,
        const cudaStream_t& stream = 0);

    // Releases internal resources.
    // Blocks until all iterations have completed executing on the GPU.
    ~BackwardPass();

    // Performs one backward iteration of the LSTM cell.
    //
    // Note that BackwardPass must be iterated in the reverse order as ForwardPass.
    // If ForwardPass iterates from 0 to T-1, BackwardPass needs to iterate from
    // T-1 down to 0. When iteration numbers are described, they will be based on the
    // iteration index (i.e., the T-1'th iteration of the forward pass is the last call
    // to ForwardPass::Iterate, whereas it is the first call to BackwardPass::Iterate).
    //
    // W_t: [H*4,C] the transpose of the input weight matrix.
    // R_t: [H*4,H] the transpose of the recurrent weight matrix.
    // b: [H*4] the bias vector.
    // x_t: [C,N] the transpose of the LSTM input for this iteration.
    // h: [N,H] the hidden state of the t'th iteration or the initial hidden state if this is
    //     the t=0 iteration (typically zeros).
    // c: [N,H] the t-1'th forward iteration's `c_out` or the initial cell state if this is
    //     the t=0 iteration (typically zeros).
    // c_new: [N,H] the t'th forward iteration's `c_out` vector.
    // dh_new: [N,H] the gradient of the loss with respect to `h_out` at this iteration.
    // dc_new: [N,H] the gradient of the loss with respect to `c_out` at this iteration.
    // dx: [N,C] the gradient of the loss with respect to the input at this time step.
    // dW: [C,H*4] the gradient of the loss with respect to the input weight matrix.
    // dR: [H,H*4] the gradient of the loss with respect to the recurrent weight matrix.
    // db: [H*4] the gradient of the loss with respect to the bias vector.
    // dh: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros
    //     for the T-1'th iteration and the same pointer should be passed in for each
    //     iteration. After a complete backward pass, this vector will contain the gradient
    //     of the loss with respect to the initial hidden state.
    // dc: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros
    //     for the T-1'th iteration and the same pointer should be passed in for each
    //     iteration. After a complete backward pass, this vector will contain the gradient
    //     of the loss with respect to the initial cell state.
    // v: [N,H*4] the same tensor that was passed to `ForwardPass::Iterate` on its corresponding
    //     iteration.
    // zoneout_mask: [N,H] may be null if zoneout was disabled in the forward pass. This vector
    //     must be the same as the one provided during the corresponding forward iteration.
    void Iterate(
        const cudaStream_t& stream,
        const T* W_t,
        const T* R_t,
        const T* b,
        const T* x_t,
        const T* h,
        const T* c,
        const T* c_new,
        const T* dh_new,
        const T* dc_new,
        T* dx,
        T* dW,
        T* dR,
        T* db,
        T* dh,
        T* dc,
        T* v,
        const T* zoneout_mask);

    // Runs the LSTM backward pass over all time steps. This method is faster than using a
    // per-step `Iterate` but requires that the entire input sequence be available upfront.
    // In some situations, this constraint may not be satisfiable (e.g. autoregressive models).
    // Users should prefer calling `Run` over `Iterate` whenever possible.
    //
    // steps: the number of iterations to run (i.e. T).
    // W_t: [H*4,C] the transpose of the input weight matrix.
    // R_t: [H*4,H] the transpose of the recurrent weight matrix.
    // b: [H*4] the bias vector.
    // x_t: [C,T,N] the transpose of the LSTM input for this iteration.
    // h: [T+1,N,H] the hidden state vectors after running `ForwardPass::Run`.
    // c: [T+1,N,H] the cell state vectors after running `ForwardPass::Run`.
    // dh_new: [T+1,N,H] the gradient of the loss with respect to `h`.
    // dc_new: [T+1,N,H] the gradient of the loss with respect to `c`.
    // dx: [T,N,C] the gradient of the loss with respect to the input.
    // dW: [C,H*4] the gradient of the loss with respect to the input weight matrix.
    // dR: [H,H*4] the gradient of the loss with respect to the recurrent weight matrix.
    // db: [H*4] the gradient of the loss with respect to the bias vector.
    // dh: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros.
    //     When this function returns, `dh` will contain the gradient of the loss with respect
    //     to the initial hidden state.
    // dc: [N,H] NOTE: this is an input and output parameter. Should be initialized to zeros.
    //     When this function returns, `dc` will contain the gradient of the loss with respect
    //     to the initial cell state.
    // v: [T,N,H*4] the same tensor that was passed to `ForwardPass::Run`.
    // zoneout_mask: [T,N,H] may be null if zoneout was disabled in the forward pass. This
    //     vector must be the same as the one provided during the forward pass.
    void Run(
        const int steps,
        const T* W_t,
        const T* R_t,
        const T* b,
        const T* x_t,
        const T* h,
        const T* c,
        const T* dh_new,
        const T* dc_new,
        T* dx,
        T* dW,
        T* dR,
        T* db,
        T* dh,
        T* dc,
        T* v,
        const T* zoneout_mask);

  private:
    void IterateInternal(
        const T* R_t,
        const T* c,
        const T* c_new,
        const T* dh_new,
        const T* dc_new,
        T* db,
        T* dh,
        T* dc,
        T* v,
        const T* zoneout_mask);
    struct private_data;
    private_data* data_;
};

}  // namespace lstm
}  // namespace v0
}  // namespace haste


extern "C"
__host__ __device__
void __assertfail(
    const char * __assertion,
    const char *__file,
    unsigned int __line,
    const char *__function,
    size_t charsize);

#define device_assert_fail(msg) \
      __assertfail((msg), __FILE__, __LINE__, __PRETTY_FUNCTION__, sizeof(char))



template<typename T>
__device__ __forceinline__
T sigmoid(const T x) {
  return static_cast<T>(1.0) / (static_cast<T>(1.0) + exp(-x));
}

template<typename T>
__device__ __forceinline__
T tanh(const T x) {
  return std::tanh(x);
}

template<typename T>
__device__ __forceinline__
T d_sigmoid(const T sigmoid_output) {
  return sigmoid_output * (static_cast<T>(1.0) - sigmoid_output);
}

template<typename T>
__device__ __forceinline__
T d_tanh(const T tanh_output) {
  return (static_cast<T>(1.0) - tanh_output * tanh_output);
}

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)

__device__ __forceinline__
double atomicAdd(double* address, double val) {
  unsigned long long int* address_as_ull = (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

#endif

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ >= 600)

template<>
__device__ __forceinline__
half sigmoid(const half x) {
  return static_cast<half>(1.0) / (static_cast<half>(1.0) + hexp(-x));
}

template<>
__device__ __forceinline__
half tanh(const half x) {
  return std::tanh(float(x));
}
#endif


namespace haste {
namespace v0 {
namespace gru {

namespace {

template<typename T, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* h,
                         const T* v,
                         const T* dh_new,
                         T* dbx_out,
                         T* dbr_out,
                         T* dh_inout,
                         T* dp_out,
                         T* dq_out,
                         const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  const int base_idx = col * hidden_dim + row;

  T dh_total = dh_new[base_idx] + dh_inout[base_idx];

  const int stride4_base_idx = col * (hidden_dim * 4) + row;
  const int z_idx = stride4_base_idx + 0 * hidden_dim;
  const int r_idx = stride4_base_idx + 1 * hidden_dim;
  const int g_idx = stride4_base_idx + 2 * hidden_dim;
  const int q_g_idx = stride4_base_idx + 3 * hidden_dim;

  const T z = v[z_idx];
  const T r = v[r_idx];
  const T g = v[g_idx];
  const T q_g = v[q_g_idx];

  if (ApplyZoneout) {
    const T mask = zoneout_mask[base_idx];
    dh_inout[base_idx] = (static_cast<T>(1.0) - mask) * dh_total;
    dh_total = mask * dh_total;
    dh_inout[base_idx] += z * dh_total;
  } else {
    dh_inout[base_idx] = z * dh_total;
  }

  const T dg = (static_cast<T>(1.0) - z) * dh_total;
  const T dz = (h[base_idx] - g) * dh_total;
  const T dp_g = d_tanh(g) * dg;
  const T dq_g = dp_g * r;
  const T dr = dp_g * q_g;
  const T dp_r = d_sigmoid(r) * dr;
  const T dq_r = dp_r;
  const T dp_z = d_sigmoid(z) * dz;
  const T dq_z = dp_z;

  const int idx = col * (hidden_dim * 3) + row;

  dp_out[idx + 0 * hidden_dim] = dp_z;
  dp_out[idx + 1 * hidden_dim] = dp_r;
  dp_out[idx + 2 * hidden_dim] = dp_g;

  dq_out[idx + 0 * hidden_dim] = dq_z;
  dq_out[idx + 1 * hidden_dim] = dq_r;
  dq_out[idx + 2 * hidden_dim] = dq_g;

  atomicAdd(&dbx_out[row + 0 * hidden_dim], dp_z);
  atomicAdd(&dbx_out[row + 1 * hidden_dim], dp_r);
  atomicAdd(&dbx_out[row + 2 * hidden_dim], dp_g);

  atomicAdd(&dbr_out[row + 0 * hidden_dim], dq_z);
  atomicAdd(&dbr_out[row + 1 * hidden_dim], dq_r);
  atomicAdd(&dbr_out[row + 2 * hidden_dim], dq_g);
}

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 700)
template<typename T, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const half* h,
                         const half* v,
                         const half* dh_new,
                         half* dbx_out,
                         half* dbr_out,
                         half* dh_inout,
                         half* dp_out,
                         half* dq_out,
                         const half* zoneout_mask) {
  device_assert_fail("FP16 is not supported on compute capability < 7.0.");
}
#endif

}  // anonymous namespace

template<typename T>
struct BackwardPass<T>::private_data {
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[2];
  cudaEvent_t event;
  cudaStream_t sync_stream;
};

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
BackwardPass<T>::~BackwardPass() {
  if (data_->sync_stream) {
    cudaEventRecord(data_->event, data_->stream[1]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[0]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
  } else {
    cudaStreamSynchronize(data_->stream[1]);
    cudaStreamSynchronize(data_->stream[0]);
  }
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void BackwardPass<T>::Iterate(
    const T* W_t,     // [H*3,C]
    const T* R_t,     // [H*3,H]
    const T* bx,      // [H*3]
    const T* br,      // [H*3]
    const T* x_t,     // [C,N]
    const T* h,       // [N,H]
    const T* v,       // [N,H*4]
    const T* dh_new,  // [N,H]
    T* dx,            // [N,C]
    T* dW,            // [C,H*3]
    T* dR,            // [H,H*3]
    T* dbx,           // [H*3]
    T* dbr,           // [H*3]
    T* dh,            // [N,H]
    T* dp,            // [N,H*3]
    T* dq,            // [N,H*3]
    const T* zoneout_mask) {  // [N,H]
  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);
  const T beta_assign = static_cast<T>(0.0);

  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const int input_size = data_->input_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  IterateInternal(
      R_t,
      h,
      v,
      dh_new,
      dbx,
      dbr,
      dh,
      dp,
      dq,
      zoneout_mask);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, input_size, batch_size,
      &alpha,
      dp, hidden_size * 3,
      x_t, batch_size,
      &beta_sum,
      dW, hidden_size * 3);

  // Wait for pointwise operations to complete since there's a
  // data dependency between its output (`dp`, `dq`) and the following matmuls.
  cudaStreamWaitEvent(stream2, event, 0);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, batch_size, hidden_size * 3,
      &alpha,
      W_t, input_size,
      dp, hidden_size * 3,
      &beta_assign,
      dx, input_size);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_T,
      hidden_size * 3, hidden_size, batch_size,
      &alpha,
      dq, hidden_size * 3,
      h, hidden_size,
      &beta_sum,
      dR, hidden_size * 3);

  cublasSetStream(blas_handle, save_stream);
}

template<typename T>
void BackwardPass<T>::IterateInternal(
    const T* R_t,     // [H*3,H]
    const T* h,       // [N,H]
    const T* v,       // [N,H*4]
    const T* dh_new,  // [N,H]
    T* dbx,           // [H*3]
    T* dbr,           // [H*3]
    T* dh,            // [N,H]
    T* dp,            // [N,H*3]
    T* dq,            // [N,H*3]
    const T* zoneout_mask) {  // [N,H]
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);

  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(32, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  if (zoneout_mask) {
    PointwiseOperations<T, true><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        h,
        v,
        dh_new,
        dbx,
        dbr,
        dh,
        dp,
        dq,
        zoneout_mask
    );
  } else {
    PointwiseOperations<T, false><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        h,
        v,
        dh_new,
        dbx,
        dbr,
        dh,
        dp,
        dq,
        nullptr
    );
  }
  cudaEventRecord(event, stream1);

  cublasSetStream(blas_handle,  stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, batch_size, hidden_size * 3,
      &alpha,
      R_t, hidden_size,
      dq, hidden_size * 3,
      &beta_sum,
      dh, hidden_size);
}

template<typename T>
void BackwardPass<T>::Run(
    const int steps,
    const T* W_t,
    const T* R_t,
    const T* bx,
    const T* br,
    const T* x_t,
    const T* h,
    const T* v,
    const T* dh_new,
    T* dx,
    T* dW,
    T* dR,
    T* dbx,
    T* dbr,
    T* dh,
    T* dp,
    T* dq,
    const T* zoneout_mask) {
  const blas<void>::enable_tensor_cores scoped0(data_->blas_handle);
  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);
  const T beta_assign = static_cast<T>(0.0);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  const int NH = batch_size * hidden_size;
  for (int i = steps - 1; i >= 0; --i) {
    IterateInternal(
        R_t,
        h + i * NH,
        v + i * NH * 4,
        dh_new + (i + 1) * NH,
        dbx,
        dbr,
        dh,
        dp + i * NH * 3,
        dq + i * NH * 3,
        zoneout_mask ? zoneout_mask + i * NH : nullptr );
  }

  // Wait for pointwise operations to complete since there's a
  // data dependency between its output (`dp`, `dq`) and the following matmuls.
  cudaStreamWaitEvent(stream2, event, 0);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, batch_size * steps, hidden_size * 3,
      &alpha,
      W_t, input_size,
      dp, hidden_size * 3,
      &beta_assign,
      dx, input_size);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_T,
      hidden_size * 3, hidden_size, batch_size * steps,
      &alpha,
      dq, hidden_size * 3,
      h, hidden_size,
      &beta_sum,
      dR, hidden_size * 3);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, input_size, batch_size * steps,
      &alpha,
      dp, hidden_size * 3,
      x_t, batch_size * steps,
      &beta_sum,
      dW, hidden_size * 3);

  cublasSetStream(blas_handle, save_stream);
}

template struct BackwardPass<half>;
template struct BackwardPass<float>;
template struct BackwardPass<double>;

}  // namespace gru
}  // namespace v0
}  // namespace haste




namespace haste {
namespace v0 {
namespace gru {

namespace {

template<typename T, bool Training, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* Wx,
                         const T* Rh,
                         const T* bx,
                         const T* br,
                         const T* h,
                         T* h_out,
                         T* v,
                         const T zoneout_prob,
                         const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  const int weight_idx = col * (hidden_dim * 3) + row;

  // Index into the `h` and `h_out` vectors (they have a stride of `hidden_dim`).
  const int output_idx = col * hidden_dim + row;

  // Indicies into the Wx and Rh matrices (for each of the u, r, and e components).
  const int z_idx = weight_idx + 0 * hidden_dim;
  const int r_idx = weight_idx + 1 * hidden_dim;
  const int g_idx = weight_idx + 2 * hidden_dim;

  // Indices into the bias vectors (for each of the u, r, and e components).
  const int bz_idx = row + 0 * hidden_dim;
  const int br_idx = row + 1 * hidden_dim;
  const int bg_idx = row + 2 * hidden_dim;

  const T z = sigmoid(Wx[z_idx] + Rh[z_idx] + bx[bz_idx] + br[bz_idx]);
  const T r = sigmoid(Wx[r_idx] + Rh[r_idx] + bx[br_idx] + br[br_idx]);
  const T g = tanh   (Wx[g_idx] + r * (Rh[g_idx] + br[bg_idx]) + bx[bg_idx]);

  // Store internal activations if we're eventually going to backprop.
  if (Training) {
    const int base_v_idx = col * (hidden_dim * 4) + row;
    v[base_v_idx + 0 * hidden_dim] = z;
    v[base_v_idx + 1 * hidden_dim] = r;
    v[base_v_idx + 2 * hidden_dim] = g;
    v[base_v_idx + 3 * hidden_dim] = Rh[g_idx] + br[bg_idx];
  }

  T cur_h_value = z * h[output_idx] + (static_cast<T>(1.0) - z) * g;

  if (ApplyZoneout) {
    if (Training) {
      cur_h_value = (cur_h_value - h[output_idx]) * zoneout_mask[output_idx] + h[output_idx];
    } else {
      cur_h_value = (zoneout_prob * h[output_idx]) + ((static_cast<T>(1.0) - zoneout_prob) * cur_h_value);
    }
  }

  h_out[output_idx] = cur_h_value;
}

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 700)
template<typename T, bool Training, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const half* Wx,
                         const half* Rh,
                         const half* bx,
                         const half* br,
                         const half* h,
                         half* h_out,
                         half* v,
                         const half zoneout_prob,
                         const half* zoneout_mask) {
  device_assert_fail("FP16 is not supported on compute capability < 7.0.");
}
#endif
}  // anonymous namespace

template<typename T>
struct ForwardPass<T>::private_data {
  bool training;
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[2];
  cudaEvent_t event;
  cudaStream_t sync_stream;
};

template<typename T>
ForwardPass<T>::ForwardPass(
    const bool training,
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->training = training;
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
ForwardPass<T>::~ForwardPass() {
  if (data_->sync_stream) {
    cudaEventRecord(data_->event, data_->stream[1]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[0]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
  } else {
    cudaStreamSynchronize(data_->stream[1]);
    cudaStreamSynchronize(data_->stream[0]);
  }
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void ForwardPass<T>::Iterate(
    const T* W,  // [C,H*3]
    const T* R,  // [H,H*3]
    const T* bx, // [H*3]
    const T* br, // [H*3]
    const T* x,  // [N,C]
    const T* h,  // [N,H]
    T* h_out,    // [N,H]
    T* v,        // [N,H*4]
    T* tmp_Wx,   // [N,H*3]
    T* tmp_Rh,   // [N,H*3]
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, batch_size, input_size,
      &alpha,
      W, hidden_size * 3,
      x, input_size,
      &beta,
      tmp_Wx, hidden_size * 3);
  cudaEventRecord(event, stream2);

  IterateInternal(
      R,
      bx,
      br,
      h,
      h_out,
      v,
      tmp_Wx,
      tmp_Rh,
      zoneout_prob,
      zoneout_mask);

  cublasSetStream(blas_handle, save_stream);
}

template<typename T>
void ForwardPass<T>::IterateInternal(
    const T* R,  // [H,H*3]
    const T* bx, // [H*3]
    const T* br, // [H*3]
    const T* h,  // [N,H]
    T* h_out,    // [N,H]
    T* v,        // [N,H*4]
    T* tmp_Wx,   // [N,H*3]
    T* tmp_Rh,   // [N,H*3]
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  // Constants for GEMM
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const bool training = data_->training;
  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, batch_size, hidden_size,
      &alpha,
      R, hidden_size * 3,
      h, hidden_size,
      &beta,
      tmp_Rh, hidden_size * 3);

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(32, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  cudaStreamWaitEvent(stream1, event, 0);

  if (training) {
    if (zoneout_prob && zoneout_mask) {
      PointwiseOperations<T, true, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx,
          tmp_Rh,
          bx,
          br,
          h,
          h_out,
          v,
          zoneout_prob,
          zoneout_mask);
    } else {
      PointwiseOperations<T, true, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx,
          tmp_Rh,
          bx,
          br,
          h,
          h_out,
          v,
          0.0f,
          nullptr);
    }
  } else {
    if (zoneout_prob && zoneout_mask) {
      PointwiseOperations<T, false, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx,
          tmp_Rh,
          bx,
          br,
          h,
          h_out,
          nullptr,
          zoneout_prob,
          zoneout_mask);
    } else {
      PointwiseOperations<T, false, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx,
          tmp_Rh,
          bx,
          br,
          h,
          h_out,
          nullptr,
          0.0f,
          nullptr);
    }
  }
}

template<typename T>
void ForwardPass<T>::Run(
    const int steps,
    const T* W,  // [C,H*3]
    const T* R,  // [H,H*3]
    const T* bx, // [H*3]
    const T* br, // [H*3]
    const T* x,  // [N,C]
    T* h,        // [N,H]
    T* v,        // [N,H*4]
    T* tmp_Wx,   // [N,H*3]
    T* tmp_Rh,   // [N,H*3]
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::enable_tensor_cores scoped0(data_->blas_handle);
  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, steps * batch_size, input_size,
      &alpha,
      W, hidden_size * 3,
      x, input_size,
      &beta,
      tmp_Wx, hidden_size * 3);
  cudaEventRecord(event, stream2);

  const int NH = batch_size * hidden_size;
  for (int i = 0; i < steps; ++i) {
    IterateInternal(
        R,
        bx,
        br,
        h + i * NH,
        h + (i + 1) * NH,
        v + i * NH * 4,
        tmp_Wx + i * NH * 3,
        tmp_Rh,
        zoneout_prob,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }

  cublasSetStream(blas_handle, save_stream);
}

template struct ForwardPass<half>;
template struct ForwardPass<float>;
template struct ForwardPass<double>;

}  // namespace gru
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool ApplyZoneout>
__global__
void IndrnnBwdOps(
    const int steps,
    const int batch_size,
    const int hidden_size,
    const T* u,
    const T* h_prev,
    const T* h,
    const T* dh_new,
    T* du_out,
    T* db_out,
    T* dh_inout,
    T* dk_out,
    const T* zoneout_mask) {
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  const int NH = batch_size * hidden_size;
  const int idx = col * hidden_size + row;

  const T u_row = u[row];
  T dh_inout_idx = dh_inout[idx];
  T du_sum = static_cast<T>(0.0);
  T db_sum = static_cast<T>(0.0);

  for (int i = (steps - 1) * NH; i >= 0; i -= NH) {
    T dh_total = dh_new[idx + i] + dh_inout_idx;
    T dh = static_cast<T>(0.0);
    if (ApplyZoneout) {
      const T mask = zoneout_mask[idx + i];
      dh = (static_cast<T>(1.0) - mask) * dh_total;
      dh_total = mask * dh_total;
    }

    const T dk = d_tanh(h[idx + i]) * dh_total;

    dk_out[idx + i] = dk;
    dh_inout_idx = dh + u_row * dk;
    du_sum += h_prev[idx + i] * dk;
    db_sum += dk;
  }

  dh_inout[idx] = dh_inout_idx;
  atomicAdd(&du_out[row], du_sum);
  atomicAdd(&db_out[row], db_sum);
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace indrnn {

template<typename T>
struct BackwardPass<T>::private_data {
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream;
  cudaStream_t sync_stream;
};

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream);
}

template<typename T>
BackwardPass<T>::~BackwardPass() {
  if (data_->sync_stream) {
    cudaEvent_t event;
    cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
    cudaEventRecord(event, data_->stream);
    cudaStreamWaitEvent(data_->sync_stream, event, 0);
    cudaEventDestroy(event);
  } else {
    cudaStreamSynchronize(data_->stream);
  }
  cudaStreamDestroy(data_->stream);
  delete data_;
}

template<typename T>
void BackwardPass<T>::Run(
    const int steps,
    const T* W_t,
    const T* u,
    const T* b,
    const T* x_t,
    const T* h,
    const T* dh_new,
    T* dx,
    T* dW,
    T* du,
    T* db,
    T* dh,
    T* workspace,
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream = data_->stream;

  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);
  const int NH = batch_size * hidden_size;
  if (zoneout_mask) {
    IndrnnBwdOps<T, true><<<gridDim, blockDim, 0, stream>>>(
        steps,
        batch_size,
        hidden_size,
        u,
        h,
        h + NH,
        dh_new + NH,
        du,
        db,
        dh,
        workspace,
        zoneout_mask);
  } else {
    IndrnnBwdOps<T, false><<<gridDim, blockDim, 0, stream>>>(
        steps,
        batch_size,
        hidden_size,
        u,
        h,
        h + NH,
        dh_new + NH,
        du,
        db,
        dh,
        workspace,
        nullptr);
  }

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, input_size, batch_size * steps,
      &alpha,
      workspace, hidden_size,
      x_t, batch_size * steps,
      &beta,
      dW, hidden_size);

  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, steps * batch_size, hidden_size,
      &alpha,
      W_t, input_size,
      workspace, hidden_size,
      &beta,
      dx, input_size);

  cublasSetStream(blas_handle, save_stream);
}

template class BackwardPass<float>;
template class BackwardPass<double>;

}  // namespace indrnn
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool Training, bool ApplyZoneout>
__global__
void IndrnnFwdOps(
    const int steps,
    const int batch_size,
    const int hidden_size,
    const T* Wx,
    const T* u,
    const T* b,
    const T* h,
    T* h_out,
    const float zoneout_prob,
    const T* zoneout_mask) {
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  const int idx = col * hidden_size + row;
  const int NH = batch_size * hidden_size;
  const T u_row = u[row];
  const T b_row = b[row];

  for (int i = 0; i < steps * NH; i += NH) {
    const T a = Wx[idx + i] + u_row * h[idx + i] + b_row;
    T cur_h_value = tanh(a);

    if (ApplyZoneout) {
      if (Training) {
        cur_h_value = (cur_h_value - h[idx + i]) * zoneout_mask[idx + i] + h[idx + i];
      } else {
        cur_h_value = (zoneout_prob * h[idx + i]) + ((1.0f - zoneout_prob) * cur_h_value);
      }
    }

    h_out[idx + i] = cur_h_value;
  }
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace indrnn {

template<typename T>
struct ForwardPass<T>::private_data {
  bool training;
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream;
  cudaStream_t sync_stream;
};

template<typename T>
ForwardPass<T>::ForwardPass(
    const bool training,
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->training = training;
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream);
}

template<typename T>
ForwardPass<T>::~ForwardPass() {
  if (data_->sync_stream) {
    cudaEvent_t event;
    cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
    cudaEventRecord(event, data_->stream);
    cudaStreamWaitEvent(data_->sync_stream, event, 0);
    cudaEventDestroy(event);
  } else {
    cudaStreamSynchronize(data_->stream);
  }
  cudaStreamDestroy(data_->stream);
  delete data_;
}

template<typename T>
void ForwardPass<T>::Run(
    const int steps,
    const T* W,
    const T* u,
    const T* b,
    const T* x,
    T* h,
    T* workspace,
    const float zoneout_prob,
    const T* zoneout_mask) {
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const bool training = data_->training;
  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream = data_->stream;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, steps * batch_size, input_size,
      &alpha,
      W, hidden_size,
      x, input_size,
      &beta,
      workspace, hidden_size);

  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);
  const int NH = batch_size * hidden_size;
  if (training) {
    if (zoneout_prob && zoneout_mask) {
      IndrnnFwdOps<T, true, true><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          zoneout_prob,
          zoneout_mask);
    } else {
      IndrnnFwdOps<T, true, false><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          0.0f,
          nullptr);
    }
  } else {
    if (zoneout_prob && zoneout_mask) {
      IndrnnFwdOps<T, false, true><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          zoneout_prob,
          zoneout_mask);
    } else {
      IndrnnFwdOps<T, false, false><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          0.0f,
          nullptr);
    }
  }

  cublasSetStream(blas_handle, save_stream);
}

template class ForwardPass<float>;
template class ForwardPass<double>;

}  // namespace indrnn
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool ApplyBeta>
__global__
void LayerNormGrad(
    const int batch_size,
    const int hidden_size,
    const T* gamma,
    const T* x,
    const T* dy,
    T* dgamma,
    T* dbeta,
    T* dx,
    T* cache) {
  const int batch = blockDim.x * blockIdx.x + threadIdx.x;
  if (batch >= batch_size)
    return;

  extern __shared__ int shared_var[];
  T* shared = reinterpret_cast<T*>(shared_var);
  const int index = threadIdx.y;
  const int stride = blockDim.y;
  const int batch_idx = batch * hidden_size;
  const int batch_block_idx = threadIdx.x * stride * 3;

  const T mean   = cache[batch * 2 + 0];
  const T invstd = cache[batch * 2 + 1];

  T dsigma_tmp = static_cast<T>(0.0);
  T dmu1_tmp = static_cast<T>(0.0);
  T dmu2_tmp = static_cast<T>(0.0);
  for (int i = index; i < hidden_size; i += stride) {
    const T cur_dy = dy[batch_idx + i];
    const T centered_x = x[batch_idx + i] - mean;
    const T z = centered_x * invstd;

    atomicAdd(&dgamma[i], z * cur_dy);
    if (ApplyBeta)
      atomicAdd(&dbeta[i], cur_dy);

    const T db = gamma[i] * cur_dy;
    dsigma_tmp += centered_x * db;
    dmu1_tmp += centered_x;
    dmu2_tmp += db;
  }
  shared[batch_block_idx + index * 3 + 0] = dsigma_tmp;
  shared[batch_block_idx + index * 3 + 1] = dmu1_tmp;
  shared[batch_block_idx + index * 3 + 2] = dmu2_tmp;
  __syncthreads();

  for (int s = stride / 2; s > 0; s >>= 1) {
    if (index < s) {
      shared[batch_block_idx + index * 3 + 0] += shared[batch_block_idx + (index + s) * 3 + 0];
      shared[batch_block_idx + index * 3 + 1] += shared[batch_block_idx + (index + s) * 3 + 1];
      shared[batch_block_idx + index * 3 + 2] += shared[batch_block_idx + (index + s) * 3 + 2];
    }
    __syncthreads();
  }

  const T dsigma = static_cast<T>(-0.5) * shared[batch_block_idx + 0] * invstd * invstd * invstd;
  const T dmu = (static_cast<T>(-2.0) * shared[batch_block_idx + 1] * dsigma / hidden_size) -
                (shared[batch_block_idx + 2] * invstd);

  for (int i = index; i < hidden_size; i += stride) {
    const T cur_dy = dy[batch_idx + i];
    const T centered_x = x[batch_idx + i] - mean;

    const T db = gamma[i] * cur_dy;
    dx[batch_idx + i] = (static_cast<T>(2.0) * centered_x * dsigma / hidden_size) +
                        (invstd * db) +
                        (dmu / hidden_size);
  }
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace layer_norm {

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int hidden_size,
    const T* gamma,
    const T* beta,
    const T* x,
    T* dgamma,
    T* dbeta,
    T* cache)
        : batch_size_(batch_size),
          hidden_size_(hidden_size),
          gamma_(gamma),
          beta_(beta),
          x_(x),
          dgamma_(dgamma),
          dbeta_(dbeta),
          cache_(cache),
          partial_(batch_size) {
}

template<typename T>
void BackwardPass<T>::Run(const cudaStream_t& stream, const T* dy, T* dx) {
  RunPartial(stream, batch_size_, dy, dx);
}

template<typename T>
void BackwardPass<T>::RunPartial(
    const cudaStream_t& stream,
    const int minibatch,
    const T* dy,
    T* dx) {
  assert(partial_ - minibatch >= 0);

  dim3 blockDim(4, 256);
  dim3 gridDim;
  gridDim.x = (minibatch + blockDim.x - 1) / blockDim.x;
  const int shared_mem_size = sizeof(T) * blockDim.x * blockDim.y * 3;

  if (beta_ && dbeta_) {
    LayerNormGrad<T, true><<<gridDim, blockDim, shared_mem_size, stream>>>(
        minibatch,
        hidden_size_,
        gamma_,
        x_ + (partial_ - minibatch) * hidden_size_,
        dy,
        dgamma_,
        dbeta_,
        dx,
        cache_ + (partial_ - minibatch) * 2);
  } else {
    LayerNormGrad<T, false><<<gridDim, blockDim, shared_mem_size, stream>>>(
        minibatch,
        hidden_size_,
        gamma_,
        x_ + (partial_ - minibatch) * hidden_size_,
        dy,
        dgamma_,
        nullptr,
        dx,
        cache_ + (partial_ - minibatch) * 2);
  }

  partial_ -= minibatch;
}

template class BackwardPass<float>;
template class BackwardPass<double>;

}  // namespace layer_norm
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool ApplyBeta>
__global__
void LayerNorm(
    const int batch_size,
    const int hidden_size,
    const T* gamma,
    const T* beta,
    const T* x,
    T* y,
    T* cache) {
  const int batch = blockDim.x * blockIdx.x + threadIdx.x;
  if (batch >= batch_size)
    return;

  extern __shared__ int shared_var[];
  T* shared = reinterpret_cast<T*>(shared_var);
  const int index = threadIdx.y;
  const int stride = blockDim.y;
  const int batch_idx = batch * hidden_size;
  const int batch_block_idx = threadIdx.x * stride;

  // TODO: use parallel single-pass algorithm to compute moments.

  // Reduce sum
  T sum = static_cast<T>(0.0);
  for (int i = index; i < hidden_size; i += stride)
    sum += x[batch_idx + i];
  shared[batch_block_idx + index] = sum;
  __syncthreads();

  for (int s = stride / 2; s > 0; s >>= 1) {
    if (index < s)
      shared[batch_block_idx + index] += shared[batch_block_idx + index + s];
    __syncthreads();
  }

  // Make sure this read completes before we start writing to `shared` again below.
  const T mean = shared[batch_block_idx] / hidden_size;
  __syncthreads();

  // Reduce squared difference
  T sumsq = static_cast<T>(0.0);
  for (int i = index; i < hidden_size; i += stride) {
    const T diff = x[batch_idx + i] - mean;
    sumsq += diff * diff;
  }
  shared[batch_block_idx + index] = sumsq;
  __syncthreads();

  for (int s = stride / 2; s > 0; s >>= 1) {
    if (index < s)
      shared[batch_block_idx + index] += shared[batch_block_idx + index + s];
    __syncthreads();
  }

  const T invstd = rsqrt(shared[batch_block_idx] / hidden_size + static_cast<T>(1e-5));

  for (int i = index; i < hidden_size; i += stride) {
    if (ApplyBeta)
      y[batch_idx + i] = (x[batch_idx + i] - mean) * invstd * gamma[i] + beta[i];
    else
      y[batch_idx + i] = (x[batch_idx + i] - mean) * invstd * gamma[i];
  }

  cache[batch * 2 + 0] = mean;
  cache[batch * 2 + 1] = invstd;
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace layer_norm {

template<typename T>
ForwardPass<T>::ForwardPass(
    const int batch_size,
    const int hidden_size,
    const T* gamma,
    const T* beta,
    T* cache)
        : batch_size_(batch_size),
          hidden_size_(hidden_size),
          gamma_(gamma),
          beta_(beta),
          cache_(cache),
          partial_(0) {
}

template<typename T>
void ForwardPass<T>::Run(const cudaStream_t& stream, const T* x, T* y) {
  RunPartial(stream, batch_size_, x, y);
}

template<typename T>
void ForwardPass<T>::RunPartial(
    const cudaStream_t& stream,
    const int minibatch,
    const T* x,
    T* y) {
  assert(partial_ + minibatch <= batch_size_);

  dim3 blockDim(4, 256);
  dim3 gridDim;
  gridDim.x = (minibatch + blockDim.x - 1) / blockDim.x;
  const int shared_mem_size = sizeof(T) * blockDim.x * blockDim.y;

  if (beta_) {
    LayerNorm<T, true><<<gridDim, blockDim, shared_mem_size, stream>>>(
        minibatch,
        hidden_size_,
        gamma_,
        beta_,
        x,
        y,
        cache_ + partial_ * 2);
  } else {
    LayerNorm<T, false><<<gridDim, blockDim, shared_mem_size, stream>>>(
        minibatch,
        hidden_size_,
        gamma_,
        nullptr,
        x,
        y,
        cache_ + partial_ * 2);
  }

  partial_ += minibatch;
}

template class ForwardPass<float>;
template class ForwardPass<double>;

}  // namespace layer_norm
}  // namespace v0
}  // namespace haste


namespace haste {
namespace v0 {
namespace layer_norm_gru {

namespace {

template<typename T, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* h,
                         const T* v,
                         const T* dh_new,
                         T* dbx_out,
                         T* dbr_out,
                         T* dh_inout,
                         T* dp_out,
                         T* dq_out,
                         const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  const int base_idx = col * hidden_dim + row;

  T dh_total = dh_new[base_idx] + dh_inout[base_idx];

  const int stride4_base_idx = col * (hidden_dim * 4) + row;
  const int z_idx = stride4_base_idx + 0 * hidden_dim;
  const int r_idx = stride4_base_idx + 1 * hidden_dim;
  const int g_idx = stride4_base_idx + 2 * hidden_dim;
  const int q_g_idx = stride4_base_idx + 3 * hidden_dim;

  const T z = v[z_idx];
  const T r = v[r_idx];
  const T g = v[g_idx];
  const T q_g = v[q_g_idx];

  if (ApplyZoneout) {
    const T mask = zoneout_mask[base_idx];
    dh_inout[base_idx] = (static_cast<T>(1.0) - mask) * dh_total;
    dh_total = mask * dh_total;
    dh_inout[base_idx] += z * dh_total;
  } else {
    dh_inout[base_idx] = z * dh_total;
  }

  const T dg = (static_cast<T>(1.0) - z) * dh_total;
  const T dz = (h[base_idx] - g) * dh_total;
  const T dp_g = d_tanh(g) * dg;
  const T dq_g = dp_g * r;
  const T dr = dp_g * q_g;
  const T dp_r = d_sigmoid(r) * dr;
  const T dq_r = dp_r;
  const T dp_z = d_sigmoid(z) * dz;
  const T dq_z = dp_z;

  const int idx = col * (hidden_dim * 3) + row;

  dp_out[idx + 0 * hidden_dim] = dp_z;
  dp_out[idx + 1 * hidden_dim] = dp_r;
  dp_out[idx + 2 * hidden_dim] = dp_g;

  dq_out[idx + 0 * hidden_dim] = dq_z;
  dq_out[idx + 1 * hidden_dim] = dq_r;
  dq_out[idx + 2 * hidden_dim] = dq_g;

  atomicAdd(&dbx_out[row + 0 * hidden_dim], dp_z);
  atomicAdd(&dbx_out[row + 1 * hidden_dim], dp_r);
  atomicAdd(&dbx_out[row + 2 * hidden_dim], dp_g);

  atomicAdd(&dbr_out[row + 0 * hidden_dim], dq_z);
  atomicAdd(&dbr_out[row + 1 * hidden_dim], dq_r);
  atomicAdd(&dbr_out[row + 2 * hidden_dim], dq_g);
}

}  // anonymous namespace

template<typename T>
struct BackwardPass<T>::private_data {
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[2];
  cudaEvent_t event;
  cudaStream_t sync_stream;
};

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
BackwardPass<T>::~BackwardPass() {
  if (data_->sync_stream) {
    cudaEventRecord(data_->event, data_->stream[1]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[0]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
  } else {
    cudaStreamSynchronize(data_->stream[1]);
    cudaStreamSynchronize(data_->stream[0]);
  }
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void BackwardPass<T>::IterateInternal(
    const T* R_t,     // [H*3,H]
    const T* h,       // [N,H]
    const T* v,       // [N,H*4]
    const T* dh_new,  // [N,H]
    T* dbx,           // [H*3]
    T* dbr,           // [H*3]
    T* dh,            // [N,H]
    T* dp,            // [N,H*3]
    T* dq,            // [N,H*3]
    layer_norm::BackwardPass<T>& layer_norm2,
    const T* zoneout_mask) {  // [N,H]
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);

  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(32, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  if (zoneout_mask) {
    PointwiseOperations<T, true><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        h,
        v,
        dh_new,
        dbx,
        dbr,
        dh,
        dp,
        dq,
        zoneout_mask
    );
  } else {
    PointwiseOperations<T, false><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        h,
        v,
        dh_new,
        dbx,
        dbr,
        dh,
        dp,
        dq,
        nullptr
    );
  }
  cudaEventRecord(event, stream1);

  cublasSetStream(blas_handle,  stream1);
  layer_norm2.RunPartial(stream1, batch_size, dq, dq);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, batch_size, hidden_size * 3,
      &alpha,
      R_t, hidden_size,
      dq, hidden_size * 3,
      &beta_sum,
      dh, hidden_size);
}

template<typename T>
void BackwardPass<T>::Run(
    const int steps,
    const T* W_t,
    const T* R_t,
    const T* bx,
    const T* br,
    const T* x_t,
    const T* h,
    const T* v,
    const T* dh_new,
    T* dx,
    T* dW,
    T* dR,
    T* dbx,
    T* dbr,
    T* dh,
    T* dp,
    T* dq,
    layer_norm::BackwardPass<T>& layer_norm1,
    layer_norm::BackwardPass<T>& layer_norm2,
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);
  const T beta_assign = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  const int NH = batch_size * hidden_size;
  for (int i = steps - 1; i >= 0; --i) {
    IterateInternal(
        R_t,
        h + i * NH,
        v + i * NH * 4,
        dh_new + (i + 1) * NH,
        dbx,
        dbr,
        dh,
        dp + i * NH * 3,
        dq + i * NH * 3,
        layer_norm2,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }

  // Wait for pointwise operations to complete since there's a
  // data dependency between its output (`dp`, `dq`) and the following matmuls.
  cudaStreamWaitEvent(stream2, event, 0);

  cublasSetStream(blas_handle, stream2);
  layer_norm1.Run(stream2, dp, dp);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, batch_size * steps, hidden_size * 3,
      &alpha,
      W_t, input_size,
      dp, hidden_size * 3,
      &beta_assign,
      dx, input_size);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_T,
      hidden_size * 3, hidden_size, batch_size * steps,
      &alpha,
      dq, hidden_size * 3,
      h, hidden_size,
      &beta_sum,
      dR, hidden_size * 3);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, input_size, batch_size * steps,
      &alpha,
      dp, hidden_size * 3,
      x_t, batch_size * steps,
      &beta_sum,
      dW, hidden_size * 3);

  cublasSetStream(blas_handle, save_stream);
}

template struct BackwardPass<float>;
template struct BackwardPass<double>;

}  // namespace layer_norm_gru
}  // namespace v0
}  // namespace haste


namespace haste {
namespace v0 {
namespace layer_norm_gru {

namespace {

template<typename T, bool Training, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* Wx,
                         const T* Rh,
                         const T* bx,
                         const T* br,
                         const T* h,
                         T* h_out,
                         T* v,
                         const float zoneout_prob,
                         const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  const int weight_idx = col * (hidden_dim * 3) + row;

  // Index into the `h` and `h_out` vectors (they have a stride of `hidden_dim`).
  const int output_idx = col * hidden_dim + row;

  // Indicies into the Wx and Rh matrices (for each of the u, r, and e components).
  const int z_idx = weight_idx + 0 * hidden_dim;
  const int r_idx = weight_idx + 1 * hidden_dim;
  const int g_idx = weight_idx + 2 * hidden_dim;

  // Indices into the bias vectors (for each of the u, r, and e components).
  const int bz_idx = row + 0 * hidden_dim;
  const int br_idx = row + 1 * hidden_dim;
  const int bg_idx = row + 2 * hidden_dim;

  const T z = sigmoid(Wx[z_idx] + Rh[z_idx] + bx[bz_idx] + br[bz_idx]);
  const T r = sigmoid(Wx[r_idx] + Rh[r_idx] + bx[br_idx] + br[br_idx]);
  const T g = tanh   (Wx[g_idx] + r * (Rh[g_idx] + br[bg_idx]) + bx[bg_idx]);

  // Store internal activations if we're eventually going to backprop.
  if (Training) {
    const int base_v_idx = col * (hidden_dim * 4) + row;
    v[base_v_idx + 0 * hidden_dim] = z;
    v[base_v_idx + 1 * hidden_dim] = r;
    v[base_v_idx + 2 * hidden_dim] = g;
    v[base_v_idx + 3 * hidden_dim] = Rh[g_idx] + br[bg_idx];
  }

  T cur_h_value = z * h[output_idx] + (static_cast<T>(1.0) - z) * g;

  if (ApplyZoneout) {
    if (Training) {
      cur_h_value = (cur_h_value - h[output_idx]) * zoneout_mask[output_idx] + h[output_idx];
    } else {
      cur_h_value = (zoneout_prob * h[output_idx]) + ((1.0f - zoneout_prob) * cur_h_value);
    }
  }

  h_out[output_idx] = cur_h_value;
}

}  // anonymous namespace

template<typename T>
struct ForwardPass<T>::private_data {
  bool training;
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[2];
  cudaEvent_t event;
  cudaStream_t sync_stream;
};

template<typename T>
ForwardPass<T>::ForwardPass(
    const bool training,
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->training = training;
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
ForwardPass<T>::~ForwardPass() {
  if (data_->sync_stream) {
    cudaEventRecord(data_->event, data_->stream[1]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[0]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
  } else {
    cudaStreamSynchronize(data_->stream[1]);
    cudaStreamSynchronize(data_->stream[0]);
  }
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void ForwardPass<T>::IterateInternal(
    const T* R,     // [H,H*3]
    const T* bx,    // [H*3]
    const T* br,    // [H*3]
    const T* h,     // [N,H]
    T* h_out,       // [N,H]
    T* v,           // [N,H*4]
    T* tmp_Wx_norm, // [N,H*3]
    T* act_Rh,      // [N,H*3]
    layer_norm::ForwardPass<T>& layer_norm2,
    T* tmp_Rh_norm,
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  // Constants for GEMM
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const bool training = data_->training;
  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, batch_size, hidden_size,
      &alpha,
      R, hidden_size * 3,
      h, hidden_size,
      &beta,
      act_Rh, hidden_size * 3);
  layer_norm2.RunPartial(stream1, batch_size, act_Rh, tmp_Rh_norm);

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(32, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  cudaStreamWaitEvent(stream1, event, 0);

  if (training) {
    if (zoneout_prob && zoneout_mask) {
      PointwiseOperations<T, true, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx_norm,
          tmp_Rh_norm,
          bx,
          br,
          h,
          h_out,
          v,
          zoneout_prob,
          zoneout_mask);
    } else {
      PointwiseOperations<T, true, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx_norm,
          tmp_Rh_norm,
          bx,
          br,
          h,
          h_out,
          v,
          0.0f,
          nullptr);
    }
  } else {
    if (zoneout_prob && zoneout_mask) {
      PointwiseOperations<T, false, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx_norm,
          tmp_Rh_norm,
          bx,
          br,
          h,
          h_out,
          nullptr,
          zoneout_prob,
          zoneout_mask);
    } else {
      PointwiseOperations<T, false, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          tmp_Wx_norm,
          tmp_Rh_norm,
          bx,
          br,
          h,
          h_out,
          nullptr,
          0.0f,
          nullptr);
    }
  }
}

template<typename T>
void ForwardPass<T>::Run(
    const int steps,
    const T* W,  // [C,H*3]
    const T* R,  // [H,H*3]
    const T* bx, // [H*3]
    const T* br, // [H*3]
    const T* x,  // [N,C]
    T* h,        // [N,H]
    T* v,        // [N,H*4]
    T* act_Wx,   // [N,H*3]
    layer_norm::ForwardPass<T>& layer_norm1,
    T* tmp_Wx_norm,
    T* act_Rh,   // [N,H*3]
    layer_norm::ForwardPass<T>& layer_norm2,
    T* tmp_Rh_norm,
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 3, steps * batch_size, input_size,
      &alpha,
      W, hidden_size * 3,
      x, input_size,
      &beta,
      act_Wx, hidden_size * 3);
  layer_norm1.Run(stream2, act_Wx, tmp_Wx_norm);
  cudaEventRecord(event, stream2);

  const int NH = batch_size * hidden_size;
  for (int i = 0; i < steps; ++i) {
    IterateInternal(
        R,
        bx,
        br,
        h + i * NH,
        h + (i + 1) * NH,
        v + i * NH * 4,
        tmp_Wx_norm + i * NH * 3,
        act_Rh + i * NH * 3,
        layer_norm2,
        tmp_Rh_norm,
        zoneout_prob,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }

  cublasSetStream(blas_handle, save_stream);
}

template struct ForwardPass<float>;
template struct ForwardPass<double>;

}  // namespace layer_norm_gru
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool ApplyZoneout>
__global__
void LayerNormIndrnnBwdOps(
    const int steps,
    const int batch_size,
    const int hidden_size,
    const T* u,
    const T* h_prev,
    const T* h,
    const T* dh_new,
    T* du_out,
    T* db_out,
    T* dh_inout,
    T* dk_out,
    const T* zoneout_mask) {
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  const int NH = batch_size * hidden_size;
  const int idx = col * hidden_size + row;

  const T u_row = u[row];
  T dh_inout_idx = dh_inout[idx];
  T du_sum = static_cast<T>(0.0);
  T db_sum = static_cast<T>(0.0);

  for (int i = (steps - 1) * NH; i >= 0; i -= NH) {
    T dh_total = dh_new[idx + i] + dh_inout_idx;
    T dh = static_cast<T>(0.0);
    if (ApplyZoneout) {
      const T mask = zoneout_mask[idx + i];
      dh = (static_cast<T>(1.0) - mask) * dh_total;
      dh_total = mask * dh_total;
    }

    const T dk = d_tanh(h[idx + i]) * dh_total;

    dk_out[idx + i] = dk;
    dh_inout_idx = dh + u_row * dk;
    du_sum += h_prev[idx + i] * dk;
    db_sum += dk;
  }

  dh_inout[idx] = dh_inout_idx;
  atomicAdd(&du_out[row], du_sum);
  atomicAdd(&db_out[row], db_sum);
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace layer_norm_indrnn {

template<typename T>
struct BackwardPass<T>::private_data {
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream;
  cudaStream_t sync_stream;
};

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream);
}

template<typename T>
BackwardPass<T>::~BackwardPass() {
  if (data_->sync_stream) {
    cudaEvent_t event;
    cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
    cudaEventRecord(event, data_->stream);
    cudaStreamWaitEvent(data_->sync_stream, event, 0);
    cudaEventDestroy(event);
  } else {
    cudaStreamSynchronize(data_->stream);
  }
  cudaStreamDestroy(data_->stream);
  delete data_;
}

template<typename T>
void BackwardPass<T>::Run(
    const int steps,
    const T* W_t,
    const T* u,
    const T* b,
    const T* x_t,
    const T* h,
    const T* dh_new,
    T* dx,
    T* dW,
    T* du,
    T* db,
    T* dh,
    T* workspace,
    layer_norm::BackwardPass<T>& layer_norm1,
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream = data_->stream;

  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);
  const int NH = batch_size * hidden_size;
  if (zoneout_mask) {
    LayerNormIndrnnBwdOps<T, true><<<gridDim, blockDim, 0, stream>>>(
        steps,
        batch_size,
        hidden_size,
        u,
        h,
        h + NH,
        dh_new + NH,
        du,
        db,
        dh,
        workspace,
        zoneout_mask);
  } else {
    LayerNormIndrnnBwdOps<T, false><<<gridDim, blockDim, 0, stream>>>(
        steps,
        batch_size,
        hidden_size,
        u,
        h,
        h + NH,
        dh_new + NH,
        du,
        db,
        dh,
        workspace,
        nullptr);
  }

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream);
  layer_norm1.Run(stream, workspace, workspace);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, input_size, batch_size * steps,
      &alpha,
      workspace, hidden_size,
      x_t, batch_size * steps,
      &beta,
      dW, hidden_size);

  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, steps * batch_size, hidden_size,
      &alpha,
      W_t, input_size,
      workspace, hidden_size,
      &beta,
      dx, input_size);

  cublasSetStream(blas_handle, save_stream);
}

template class BackwardPass<float>;
template class BackwardPass<double>;

}  // namespace layer_norm_indrnn
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool Training, bool ApplyZoneout>
__global__
void LayerNormIndrnnFwdOps(
    const int steps,
    const int batch_size,
    const int hidden_size,
    const T* Wx,
    const T* u,
    const T* b,
    const T* h,
    T* h_out,
    const float zoneout_prob,
    const T* zoneout_mask) {
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  const int idx = col * hidden_size + row;
  const int NH = batch_size * hidden_size;
  const T u_row = u[row];
  const T b_row = b[row];

  for (int i = 0; i < steps * NH; i += NH) {
    const T a = Wx[idx + i] + u_row * h[idx + i] + b_row;
    T cur_h_value = tanh(a);

    if (ApplyZoneout) {
      if (Training) {
        cur_h_value = (cur_h_value - h[idx + i]) * zoneout_mask[idx + i] + h[idx + i];
      } else {
        cur_h_value = (zoneout_prob * h[idx + i]) + ((1.0f - zoneout_prob) * cur_h_value);
      }
    }

    h_out[idx + i] = cur_h_value;
  }
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace layer_norm_indrnn {

template<typename T>
struct ForwardPass<T>::private_data {
  bool training;
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream;
  cudaStream_t sync_stream;
};

template<typename T>
ForwardPass<T>::ForwardPass(
    const bool training,
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->training = training;
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream);
}

template<typename T>
ForwardPass<T>::~ForwardPass() {
  if (data_->sync_stream) {
    cudaEvent_t event;
    cudaEventCreateWithFlags(&event, cudaEventDisableTiming);
    cudaEventRecord(event, data_->stream);
    cudaStreamWaitEvent(data_->sync_stream, event, 0);
    cudaEventDestroy(event);
  } else {
    cudaStreamSynchronize(data_->stream);
  }
  cudaStreamDestroy(data_->stream);
  delete data_;
}

template<typename T>
void ForwardPass<T>::Run(
    const int steps,
    const T* W,
    const T* u,
    const T* b,
    const T* x,
    T* h,
    T* workspace,
    T* act_Wx,
    layer_norm::ForwardPass<T>& layer_norm1,
    const float zoneout_prob,
    const T* zoneout_mask) {
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const bool training = data_->training;
  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream = data_->stream;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, steps * batch_size, input_size,
      &alpha,
      W, hidden_size,
      x, input_size,
      &beta,
      act_Wx, hidden_size);
  layer_norm1.Run(stream, act_Wx, workspace);

  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);
  const int NH = batch_size * hidden_size;
  if (training) {
    if (zoneout_prob && zoneout_mask) {
      LayerNormIndrnnFwdOps<T, true, true><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          zoneout_prob,
          zoneout_mask);
    } else {
      LayerNormIndrnnFwdOps<T, true, false><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          0.0f,
          nullptr);
    }
  } else {
    if (zoneout_prob && zoneout_mask) {
      LayerNormIndrnnFwdOps<T, false, true><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          zoneout_prob,
          zoneout_mask);
    } else {
      LayerNormIndrnnFwdOps<T, false, false><<<gridDim, blockDim, 0, stream>>>(
          steps,
          batch_size,
          hidden_size,
          workspace,
          u,
          b,
          h,
          h + NH,
          0.0f,
          nullptr);
    }
  }

  cublasSetStream(blas_handle, save_stream);
}

template class ForwardPass<float>;
template class ForwardPass<double>;

}  // namespace layer_norm_indrnn
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool ApplyZoneout>
__global__
void ComputeOutputGrad(
    const int batch_size,
    const int hidden_size,
    const T* act_c_norm,
    const T* dh_new,
    T* dh_inout,
    T* dlayer_norm,
    T* v,
    const T* zoneout_mask) {
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  const int base_idx = col * hidden_size + row;
  const int stride4_base_idx = col * (hidden_size * 4) + row;
  const int o_idx = stride4_base_idx + 3 * hidden_size;

  T dh_total = dh_new[base_idx] + dh_inout[base_idx];
  if (ApplyZoneout) {
    const T mask = zoneout_mask[base_idx];
    dh_inout[base_idx] = (static_cast<T>(1.0) - mask) * dh_total;
    dh_total = mask * dh_total;
  } else {
    dh_inout[base_idx] = static_cast<T>(0.0);
  }

  const T c_tanh = tanh(act_c_norm[base_idx]);
  const T o = v[o_idx];

  const T do_ = c_tanh * dh_total;
  const T dc_tanh = o * dh_total;

  dlayer_norm[base_idx] = d_tanh(c_tanh) * dc_tanh;
  v[o_idx] = d_sigmoid(o) * do_;
}

template<typename T>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* c,
                         const T* v,
                         const T* dc_new,
                         const T* dlayer_norm,
                         T* db_out,
                         T* dc_inout,
                         T* dv_out) {
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  const int base_idx = col * hidden_dim + row;
  const int stride4_base_idx = col * (hidden_dim * 4) + row;
  const int i_idx = stride4_base_idx + 0 * hidden_dim;
  const int g_idx = stride4_base_idx + 1 * hidden_dim;
  const int f_idx = stride4_base_idx + 2 * hidden_dim;
  const int o_idx = stride4_base_idx + 3 * hidden_dim;

  const T i = v[i_idx];
  const T g = v[g_idx];
  const T f = v[f_idx];
  const T o = v[o_idx];

  const T dc_total = dc_new[base_idx] + dc_inout[base_idx] + dlayer_norm[base_idx];
  const T df = c[base_idx] * dc_total;
  const T dc = f * dc_total;
  const T di = g * dc_total;
  const T dg = i * dc_total;
  const T dv_g = d_tanh(g) * dg;
  const T dv_o = o;
  const T dv_i = d_sigmoid(i) * di;
  const T dv_f = d_sigmoid(f) * df;

  // TODO: performance optimization opportunity on this reduce operation.
  atomicAdd(&db_out[row + 0 * hidden_dim], dv_i);
  atomicAdd(&db_out[row + 1 * hidden_dim], dv_g);
  atomicAdd(&db_out[row + 2 * hidden_dim], dv_f);
  atomicAdd(&db_out[row + 3 * hidden_dim], dv_o);

  dc_inout[base_idx] = dc;

  dv_out[i_idx] = dv_i;
  dv_out[g_idx] = dv_g;
  dv_out[f_idx] = dv_f;
  dv_out[o_idx] = dv_o;
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace layer_norm_lstm {

template<typename T>
struct BackwardPass<T>::private_data {
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[3];
  cudaEvent_t event;
  cudaStream_t sync_stream;
};

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaStreamCreate(&data_->stream[2]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
BackwardPass<T>::~BackwardPass() {
  if (data_->sync_stream) {
    cudaEventRecord(data_->event, data_->stream[2]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[1]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[0]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
  } else {
    cudaStreamSynchronize(data_->stream[2]);
    cudaStreamSynchronize(data_->stream[1]);
    cudaStreamSynchronize(data_->stream[0]);
  }
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[2]);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void BackwardPass<T>::IterateInternal(
    const T* R_t,     // [H*4,H]
    const T* c,       // [N,H]
    const T* c_new,   // [N,H]
    const T* dh_new,  // [N,H]
    const T* dc_new,  // [N,H]
    T* db,            // [H*4]
    T* dh,            // [N,H]
    T* dc,            // [N,H]
    T* v,             // [N,H*4]
    T* act_Rh,
    layer_norm::BackwardPass<T>& layer_norm2,
    layer_norm::BackwardPass<T>& layer_norm3,
    T* act_c_norm,
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);  // Accumulate into output matrix!

  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  if (zoneout_mask) {
    ComputeOutputGrad<T, true><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        act_c_norm,
        dh_new,
        dh,
        act_c_norm,
        v,
        zoneout_mask);
  } else {
    ComputeOutputGrad<T, false><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        act_c_norm,
        dh_new,
        dh,
        act_c_norm,
        v,
        nullptr);
  }
  layer_norm3.RunPartial(stream1, batch_size, act_c_norm, act_c_norm);
  PointwiseOperations<T><<<gridDim, blockDim, 0, stream1>>>(
      batch_size,
      hidden_size,
      c,
      v,
      dc_new,
      act_c_norm,
      db,
      dc,
      v);

  // Signal completion of pointwise operations for data-dependent streams.
  cudaEventRecord(event, stream1);

  cublasSetStream(blas_handle, stream1);
  layer_norm2.RunPartial(stream1, batch_size, v, act_Rh);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, batch_size, hidden_size * 4,
      &alpha,
      R_t, hidden_size,
      act_Rh, hidden_size * 4,
      &beta_sum,
      dh, hidden_size);
}

template<typename T>
void BackwardPass<T>::Run(
    const int steps,
    const T* W_t,     // [H*4,C]
    const T* R_t,     // [H*4,H]
    const T* b,       // [H*4]
    const T* x_t,     // [C,T,N]
    const T* h,       // [T+1,N,H]
    const T* c,       // [T+1,N,H]
    const T* dh_new,  // [T+1,N,H]
    const T* dc_new,  // [T+1,N,H]
    T* dx,            // [T,N,C]
    T* dW,            // [C,H*4]
    T* dR,            // [H,H*4]
    T* db,            // [H*4]
    T* dh,            // [N,H]
    T* dc,            // [N,H]
    T* act_Wx,        // [T,N,H*4]
    layer_norm::BackwardPass<T>& layer_norm1,
    T* act_Wx_norm,   // [T,N,H*4]
    T* act_Rh,
    layer_norm::BackwardPass<T>& layer_norm2,
    layer_norm::BackwardPass<T>& layer_norm3,
    T* act_c_norm,
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);  // Accumulate into output matrix!
  const T beta_assign = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaStream_t stream3 = data_->stream[2];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  const int NH = batch_size * hidden_size;
  for (int i = steps - 1; i >= 0; --i) {
    IterateInternal(
        R_t,
        c + i * NH,
        c + (i + 1) * NH,
        dh_new + (i + 1) * NH,
        dc_new + (i + 1) * NH,
        db,
        dh,
        dc,
        act_Wx_norm + i * NH * 4,
        act_Rh + i * NH * 4,
        layer_norm2,
        layer_norm3,
        act_c_norm + i * NH,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }
  cudaEventRecord(event, stream1);

  cudaStreamWaitEvent(stream2, event, 0);
  layer_norm1.Run(stream2, act_Wx_norm, act_Wx);
  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, input_size, batch_size * steps,
      &alpha,
      act_Wx, hidden_size * 4,
      x_t, batch_size * steps,
      &beta_sum,
      dW, hidden_size * 4);

  cudaStreamWaitEvent(stream3, event, 0);
  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_T,
      hidden_size * 4, hidden_size, batch_size * steps,
      &alpha,
      act_Rh, hidden_size * 4,
      h, hidden_size,
      &beta_sum,
      dR, hidden_size * 4);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, steps * batch_size, hidden_size * 4,
      &alpha,
      W_t, input_size,
      act_Wx, hidden_size * 4,
      &beta_assign,
      dx, input_size);

  cublasSetStream(blas_handle, save_stream);
}

template struct BackwardPass<float>;
template struct BackwardPass<double>;

}  // namespace layer_norm_lstm
}  // namespace v0
}  // namespace haste

namespace {

// `c` and `c_out` may be aliased.
template<typename T, bool Training>
__global__
void ComputeCellState(
    const int batch_size,
    const int hidden_size,
    const T* Wx,  // Precomputed (Wx) vector
    const T* Rh,  // Precomputed (Rh) vector
    const T* b,   // Bias for gates
    const T* c,   // Input cell state
    T* c_out,     // Output cell state
    T* v_out) {   // Output vector v (Wx + Rh + b)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  // Base index into the Wx and Rh matrices.
  const int weight_idx = col * (hidden_size * 4) + row;

  // Base index into the output matrix. This is different from `weight_idx` because
  // the number of rows are different between the two sets of matrices.
  const int output_idx = col * hidden_size + row;

  const int i_idx = weight_idx + 0 * hidden_size;
  const int g_idx = weight_idx + 1 * hidden_size;
  const int f_idx = weight_idx + 2 * hidden_size;
  const int o_idx = weight_idx + 3 * hidden_size;

  const T i = sigmoid(Wx[i_idx] + Rh[i_idx] + b[row + 0 * hidden_size]);
  const T g = tanh   (Wx[g_idx] + Rh[g_idx] + b[row + 1 * hidden_size]);
  const T f = sigmoid(Wx[f_idx] + Rh[f_idx] + b[row + 2 * hidden_size]);
  const T o = sigmoid(Wx[o_idx] + Rh[o_idx] + b[row + 3 * hidden_size]);

  // Compile-time constant branch should be eliminated by compiler so we have
  // straight-through code.
  if (Training) {
    v_out[i_idx] = i;
    v_out[g_idx] = g;
    v_out[f_idx] = f;
    v_out[o_idx] = o;
  } else {
    v_out[o_idx] = o;
  }

  c_out[output_idx] = (f * c[output_idx]) + (i * g);
}

// `h` and `h_out` may be aliased.
template<typename T, bool Training, bool ApplyZoneout>
__global__
void ComputeCellOutput(
    const int batch_size,
    const int hidden_size,
    const T* h,   // Input recurrent state
    const T* c,   // Input cell state
    const T* v,
    T* h_out,     // Output recurrent state
    const float zoneout_prob,
    const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_size || col >= batch_size)
    return;

  const int weight_idx = col * (hidden_size * 4) + row;
  const int output_idx = col * hidden_size + row;

  const T o = v[weight_idx + 3 * hidden_size];
  const T cur_c_value = c[output_idx];

  T cur_h_value = o * tanh(cur_c_value);

  if (ApplyZoneout) {
    if (Training) {
      cur_h_value = (cur_h_value - h[output_idx]) * zoneout_mask[output_idx] + h[output_idx];
    } else {
      cur_h_value = (zoneout_prob * h[output_idx]) + ((1.0f - zoneout_prob) * cur_h_value);
    }
  }

  h_out[output_idx] = cur_h_value;
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace layer_norm_lstm {

template<typename T>
struct ForwardPass<T>::private_data {
  bool training;
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[2];
  cudaEvent_t event;
  cudaStream_t sync_stream;
};

template<typename T>
ForwardPass<T>::ForwardPass(
    const bool training,
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->training = training;
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  data_->sync_stream = stream;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
ForwardPass<T>::~ForwardPass() {
  if (data_->sync_stream) {
    cudaEventRecord(data_->event, data_->stream[1]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
    cudaEventRecord(data_->event, data_->stream[0]);
    cudaStreamWaitEvent(data_->sync_stream, data_->event, 0);
  } else {
    cudaStreamSynchronize(data_->stream[1]);
    cudaStreamSynchronize(data_->stream[0]);
  }
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void ForwardPass<T>::IterateInternal(
    const T* R,  // Weight matrix for recurrent state (Rh) [H,H*4]
    const T* b,  // Bias for gates (Wx + Rh + b) [H*4]
    const T* h,  // Recurrent state [N,H]
    const T* c,  // Cell state [N,H]
    T* h_out,    // Output recurrent state [N,H]
    T* c_out,    // Output cell state [N,H]
    T* v,        // Output vector (Wx + Rh + b) [N,H*4]
    T* tmp_Rh,   // Temporary storage for Rh vector [N,H*4]
    T* act_Rh,
    layer_norm::ForwardPass<T>& layer_norm2,
    layer_norm::ForwardPass<T>& layer_norm3,
    T* act_c_norm,
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const bool training = data_->training;
  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, batch_size, hidden_size,
      &alpha,
      R, hidden_size * 4,
      h, hidden_size,
      &beta,
      act_Rh, hidden_size * 4);
  layer_norm2.RunPartial(stream1, batch_size, act_Rh, tmp_Rh);
  cudaStreamWaitEvent(stream1, event, 0);

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  if (training) {
    ComputeCellState<T, true><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        v,
        tmp_Rh,
        b,
        c,
        c_out,
        v);
    layer_norm3.RunPartial(stream1, batch_size, c_out, act_c_norm);
    if (zoneout_prob && zoneout_mask) {
      ComputeCellOutput<T, true, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          h,
          act_c_norm,
          v,
          h_out,
          zoneout_prob,
          zoneout_mask);
    } else {
      ComputeCellOutput<T, true, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          h,
          act_c_norm,
          v,
          h_out,
          0.0f,
          nullptr);
    }
  } else {
    ComputeCellState<T, false><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        v,
        tmp_Rh,
        b,
        c,
        c_out,
        v);
    layer_norm3.RunPartial(stream1, batch_size, c_out, act_c_norm);
    if (zoneout_prob && zoneout_mask) {
      ComputeCellOutput<T, false, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          h,
          act_c_norm,
          v,
          h_out,
          zoneout_prob,
          zoneout_mask);
    } else {
      ComputeCellOutput<T, false, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          h,
          act_c_norm,
          v,
          h_out,
          0.0f,
          nullptr);
    }
  }
}

template<typename T>
void ForwardPass<T>::Run(
    const int steps,
    const T* W,  // Weight matrix for input (Wx) [C,H*4]
    const T* R,  // Weight matrix for recurrent state (Rh) [H,H*4]
    const T* b,  // Bias for gates (Wx + Rh + b) [H*4]
    const T* x,  // Input vector [T,N,C]
    T* h,        // Recurrent state [T+1,N,H]
    T* c,        // Cell state [T+1,N,H]
    T* act_Wx,   // Output vector (Wx + Rh + b) [T,N,H*4]
    T* tmp_Rh,   // Temporary storage for Rh vector [N,H*4]
    layer_norm::ForwardPass<T>& layer_norm1,
    T* act_Wx_norm,
    T* act_Rh,
    layer_norm::ForwardPass<T>& layer_norm2,
    layer_norm::ForwardPass<T>& layer_norm3,
    T* act_c_norm,
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [T,N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, steps * batch_size, input_size,
      &alpha,
      W, hidden_size * 4,
      x, input_size,
      &beta,
      act_Wx, hidden_size * 4);
  layer_norm1.Run(stream1, act_Wx, act_Wx_norm);

  for (int i = 0; i < steps; ++i) {
    const int NH = batch_size * hidden_size;
    IterateInternal(
        R,
        b,
        h + i * NH,
        c + i * NH,
        h + (i + 1) * NH,
        c + (i + 1) * NH,
        act_Wx_norm + i * NH * 4,
        tmp_Rh,
        act_Rh + i * NH * 4,
        layer_norm2,
        layer_norm3,
        act_c_norm + i * NH,
        zoneout_prob,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }

  cublasSetStream(blas_handle, save_stream);
}

template struct ForwardPass<float>;
template struct ForwardPass<double>;

}  // namespace layer_norm_lstm
}  // namespace v0
}  // namespace haste

namespace {

template<typename T, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* c,
                         const T* v,
                         const T* c_new,
                         const T* dh_new,
                         const T* dc_new,
                         T* db_out,
                         T* dh_inout,
                         T* dc_inout,
                         T* dv_out,
                         const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  const int base_idx = col * hidden_dim + row;

        T dc_total = dc_new[base_idx] + dc_inout[base_idx];
        T dh_total = dh_new[base_idx] + dh_inout[base_idx];
  const T c_tanh = tanh(c_new[base_idx]);

  const int stride4_base_idx = col * (hidden_dim * 4) + row;
  const int i_idx = stride4_base_idx + 0 * hidden_dim;
  const int g_idx = stride4_base_idx + 1 * hidden_dim;
  const int f_idx = stride4_base_idx + 2 * hidden_dim;
  const int o_idx = stride4_base_idx + 3 * hidden_dim;

  const T i = v[i_idx];
  const T g = v[g_idx];
  const T f = v[f_idx];
  const T o = v[o_idx];

  if (ApplyZoneout) {
    const T mask = zoneout_mask[base_idx];
    dh_inout[base_idx] = (static_cast<T>(1.0) - mask) * dh_total;
    dh_total = mask * dh_total;
  } else {
    dh_inout[base_idx] = static_cast<T>(0.0);
  }

  const T do_ = c_tanh * dh_total;
  const T dc_tanh = o * dh_total;
          dc_total += d_tanh(c_tanh) * dc_tanh;
  const T df = c[base_idx] * dc_total;
  const T dc = f * dc_total;
  const T di = g * dc_total;
  const T dg = i * dc_total;
  const T dv_g = d_tanh(g) * dg;
  const T dv_o = d_sigmoid(o) * do_;
  const T dv_i = d_sigmoid(i) * di;
  const T dv_f = d_sigmoid(f) * df;

  // TODO: performance optimization opportunity on this reduce operation.
  atomicAdd(&db_out[row + 0 * hidden_dim], dv_i);
  atomicAdd(&db_out[row + 1 * hidden_dim], dv_g);
  atomicAdd(&db_out[row + 2 * hidden_dim], dv_f);
  atomicAdd(&db_out[row + 3 * hidden_dim], dv_o);

  dc_inout[base_idx] = dc;

  dv_out[i_idx] = dv_i;
  dv_out[g_idx] = dv_g;
  dv_out[f_idx] = dv_f;
  dv_out[o_idx] = dv_o;
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace lstm {

template<typename T>
struct BackwardPass<T>::private_data {
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[3];
  cudaEvent_t event;
};

template<typename T>
BackwardPass<T>::BackwardPass(
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaStreamCreate(&data_->stream[2]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
}

template<typename T>
BackwardPass<T>::~BackwardPass() {
  cudaStreamSynchronize(data_->stream[2]);
  cudaStreamSynchronize(data_->stream[1]);
  cudaStreamSynchronize(data_->stream[0]);
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[2]);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void BackwardPass<T>::Iterate(
    const cudaStream_t& stream,
    const T* W_t,     // [H*4,C]
    const T* R_t,     // [H*4,H]
    const T* b,       // [H*4]
    const T* x_t,     // [C,N]
    const T* h,       // [N,H]
    const T* c,       // [N,H]
    const T* c_new,   // [N,H]
    const T* dh_new,  // [N,H]
    const T* dc_new,  // [N,H]
    T* dx,            // [N,C]
    T* dW,            // [C,H*4]
    T* dR,            // [H,H*4]
    T* db,            // [H*4]
    T* dh,            // [N,H]
    T* dc,            // [N,H]
    T* v,             // [N,H*4]
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);  // Accumulate into output matrix!
  const T beta_assign = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaStream_t stream3 = data_->stream[2];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  // Make sure inputs are ready before using them.
  if (stream) {
    cudaEventRecord(event, stream);
    cudaStreamWaitEvent(stream1, event, 0);
  }

  IterateInternal(
      R_t,
      c,
      c_new,
      dh_new,
      dc_new,
      db,
      dh,
      dc,
      v,
      zoneout_mask);

  // Wait for pointwise operations to complete since there's a
  // data dependency between its output (`v`) and the following matmuls.
  cudaStreamWaitEvent(stream2, event, 0);
  cudaStreamWaitEvent(stream3, event, 0);

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, batch_size, hidden_size * 4,
      &alpha,
      W_t, input_size,
      v, hidden_size * 4,
      &beta_assign,
      dx, input_size);

  // We can get away with only waiting for the `dx` and `dh` outputs and
  // let the `dR` and `dW` matrices complete whenever they complete. It's
  // a little unsafe, but we make the assumption that callers won't have
  // upstream data-dependencies on those matrices.
  if (stream) {
    cudaEventRecord(event, stream2);
    cudaStreamWaitEvent(stream, event, 0);
  }

  cublasSetStream(blas_handle, stream3);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_T,
      hidden_size * 4, hidden_size, batch_size,
      &alpha,
      v, hidden_size * 4,
      h, hidden_size,
      &beta_sum,
      dR, hidden_size * 4);

  cublasSetStream(blas_handle, stream3);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, input_size, batch_size,
      &alpha,
      v, hidden_size * 4,
      x_t, batch_size,
      &beta_sum,
      dW, hidden_size * 4);

  cublasSetStream(blas_handle, save_stream);
}

template<typename T>
void BackwardPass<T>::IterateInternal(
    const T* R_t,     // [H*4,H]
    const T* c,       // [N,H]
    const T* c_new,   // [N,H]
    const T* dh_new,  // [N,H]
    const T* dc_new,  // [N,H]
    T* db,            // [H*4]
    T* dh,            // [N,H]
    T* dc,            // [N,H]
    T* v,             // [N,H*4]
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);  // Accumulate into output matrix!

  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  if (zoneout_mask) {
    PointwiseOperations<T, true><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        c,
        v,
        c_new,
        dh_new,
        dc_new,
        db,
        dh,
        dc,
        v,
        zoneout_mask
    );
  } else {
    PointwiseOperations<T, false><<<gridDim, blockDim, 0, stream1>>>(
        batch_size,
        hidden_size,
        c,
        v,
        c_new,
        dh_new,
        dc_new,
        db,
        dh,
        dc,
        v,
        nullptr
    );
  }

  // Signal completion of pointwise operations for data-dependent streams.
  cudaEventRecord(event, stream1);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size, batch_size, hidden_size * 4,
      &alpha,
      R_t, hidden_size,
      v, hidden_size * 4,
      &beta_sum,
      dh, hidden_size);
}

template<typename T>
void BackwardPass<T>::Run(
    const int steps,
    const T* W_t,     // [H*4,C]
    const T* R_t,     // [H*4,H]
    const T* b,       // [H*4]
    const T* x_t,     // [C,T,N]
    const T* h,       // [T+1,N,H]
    const T* c,       // [T+1,N,H]
    const T* dh_new,  // [T+1,N,H]
    const T* dc_new,  // [T+1,N,H]
    T* dx,            // [T,N,C]
    T* dW,            // [C,H*4]
    T* dR,            // [H,H*4]
    T* db,            // [H*4]
    T* dh,            // [N,H]
    T* dc,            // [N,H]
    T* v,            // [T,N,H*4]
    const T* zoneout_mask) {
  const T alpha = static_cast<T>(1.0);
  const T beta_sum = static_cast<T>(1.0);  // Accumulate into output matrix!
  const T beta_assign = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaStream_t stream3 = data_->stream[2];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  const int NH = batch_size * hidden_size;
  for (int i = steps - 1; i >= 0; --i) {
    IterateInternal(
        R_t,
        c + i * NH,
        c + (i + 1) * NH,
        dh_new + (i + 1) * NH,
        dc_new + (i + 1) * NH,
        db,
        dh,
        dc,
        v + i * NH * 4,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }
  cudaEventRecord(event, stream1);

  cudaStreamWaitEvent(stream2, event, 0);
  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, input_size, batch_size * steps,
      &alpha,
      v, hidden_size * 4,
      x_t, batch_size * steps,
      &beta_sum,
      dW, hidden_size * 4);

  cudaStreamWaitEvent(stream3, event, 0);
  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_T,
      hidden_size * 4, hidden_size, batch_size * steps,
      &alpha,
      v, hidden_size * 4,
      h, hidden_size,
      &beta_sum,
      dR, hidden_size * 4);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      input_size, steps * batch_size, hidden_size * 4,
      &alpha,
      W_t, input_size,
      v, hidden_size * 4,
      &beta_assign,
      dx, input_size);

  cublasSetStream(blas_handle, save_stream);
}

template struct BackwardPass<float>;
template struct BackwardPass<double>;

}  // namespace lstm
}  // namespace v0
}  // namespace haste

namespace {

// `h` and `h_out` may be aliased.
// `c` and `c_out` may be aliased.
template<typename T, bool Training, bool ApplyZoneout>
__global__
void PointwiseOperations(const int batch_dim,
                         const int hidden_dim,
                         const T* Wx,  // Precomputed (Wx) vector
                         const T* Rh,  // Precomputed (Rh) vector
                         const T* b,   // Bias for gates
                         const T* h,   // Input recurrent state
                         const T* c,   // Input cell state
                         T* h_out,     // Output recurrent state
                         T* c_out,     // Output cell state
                         T* v_out,     // Output vector v (Wx + Rh + b) (only used if Training==true)
                         const float zoneout_prob,
                         const T* zoneout_mask) {  // Zoneout mask (only used if ApplyZoneout==true)
  // We're in column-major order here, so increase x => increase row.
  const int row = blockDim.x * blockIdx.x + threadIdx.x;
  const int col = blockDim.y * blockIdx.y + threadIdx.y;

  if (row >= hidden_dim || col >= batch_dim)
    return;

  // Base index into the Wx and Rh matrices.
  const int weight_idx = col * (hidden_dim * 4) + row;

  // Base index into the output matrix. This is different from `weight_idx` because
  // the number of rows are different between the two sets of matrices.
  const int output_idx = col * hidden_dim + row;

  const int i_idx = weight_idx + 0 * hidden_dim;
  const int g_idx = weight_idx + 1 * hidden_dim;
  const int f_idx = weight_idx + 2 * hidden_dim;
  const int o_idx = weight_idx + 3 * hidden_dim;

  const T i = sigmoid(Wx[i_idx] + Rh[i_idx] + b[row + 0 * hidden_dim]);
  const T g = tanh   (Wx[g_idx] + Rh[g_idx] + b[row + 1 * hidden_dim]);
  const T f = sigmoid(Wx[f_idx] + Rh[f_idx] + b[row + 2 * hidden_dim]);
  const T o = sigmoid(Wx[o_idx] + Rh[o_idx] + b[row + 3 * hidden_dim]);

  // Compile-time constant branch should be eliminated by compiler so we have
  // straight-through code.
  if (Training) {
    v_out[i_idx] = i;
    v_out[g_idx] = g;
    v_out[f_idx] = f;
    v_out[o_idx] = o;
  }

  T cur_c_value = (f * c[output_idx]) + (i * g);
  T cur_h_value = o * tanh(cur_c_value);

  if (ApplyZoneout) {
    if (Training) {
      cur_h_value = (cur_h_value - h[output_idx]) * zoneout_mask[output_idx] + h[output_idx];
    } else {
      cur_h_value = (zoneout_prob * h[output_idx]) + ((1.0f - zoneout_prob) * cur_h_value);
    }
  }

  c_out[output_idx] = cur_c_value;
  h_out[output_idx] = cur_h_value;
}

}  // anonymous namespace

namespace haste {
namespace v0 {
namespace lstm {

template<typename T>
struct ForwardPass<T>::private_data {
  bool training;
  int batch_size;
  int input_size;
  int hidden_size;
  cublasHandle_t blas_handle;
  cudaStream_t stream[2];
  cudaEvent_t event;
  cudaEvent_t ready_event;
  cudaEvent_t finished_event;
};

template<typename T>
ForwardPass<T>::ForwardPass(
    const bool training,
    const int batch_size,
    const int input_size,
    const int hidden_size,
    const cublasHandle_t& blas_handle,
    const cudaStream_t& stream) : data_(new private_data) {
  data_->training = training;
  data_->batch_size = batch_size;
  data_->input_size = input_size;
  data_->hidden_size = hidden_size;
  data_->blas_handle = blas_handle;
  cudaStreamCreate(&data_->stream[0]);
  cudaStreamCreate(&data_->stream[1]);
  cudaEventCreateWithFlags(&data_->event, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&data_->ready_event, cudaEventDisableTiming);
  cudaEventCreateWithFlags(&data_->finished_event, cudaEventDisableTiming);
}

template<typename T>
ForwardPass<T>::~ForwardPass() {
  cudaStreamSynchronize(data_->stream[1]);
  cudaStreamSynchronize(data_->stream[0]);
  cudaEventDestroy(data_->finished_event);
  cudaEventDestroy(data_->ready_event);
  cudaEventDestroy(data_->event);
  cudaStreamDestroy(data_->stream[1]);
  cudaStreamDestroy(data_->stream[0]);
  delete data_;
}

template<typename T>
void ForwardPass<T>::Iterate(
    const cudaStream_t& stream,
    const T* W,  // Weight matrix for input (Wx) [C,H*4]
    const T* R,  // Weight matrix for recurrent state (Rh) [H,H*4]
    const T* b,  // Bias for gates (Wx + Rh + b) [H*4]
    const T* x,  // Input vector [N,C]
    const T* h,  // Recurrent state [N,H]
    const T* c,  // Cell state [N,H]
    T* h_out,    // Output recurrent state [N,H]
    T* c_out,    // Output cell state [N,H]
    T* v,        // Output vector (Wx + Rh + b) [N,H*4]
    T* tmp_Rh,   // Temporary storage for Rh vector [N,H*4]
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  // Constants for GEMM
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaStream_t stream2 = data_->stream[1];
  const cudaEvent_t event = data_->event;

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  // Make sure inputs are ready before we use them.
  if (stream) {
    cudaEventRecord(event, stream);
    cudaStreamWaitEvent(stream2, event, 0);
  }

  cublasSetStream(blas_handle, stream2);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, batch_size, input_size,
      &alpha,
      W, hidden_size * 4,
      x, input_size,
      &beta,
      v, hidden_size * 4);
  cudaEventRecord(event, stream2);

  IterateInternal(
      R,
      b,
      h,
      c,
      h_out,
      c_out,
      v,
      tmp_Rh,
      zoneout_prob,
      zoneout_mask);

  // Make sure outputs have settled.
  if (stream) {
    cudaEventRecord(event, stream1);
    cudaStreamWaitEvent(stream, event, 0);
  }

  cublasSetStream(blas_handle, save_stream);
}

template<typename T>
void ForwardPass<T>::IterateInternal(
    const T* R,  // Weight matrix for recurrent state (Rh) [H,H*4]
    const T* b,  // Bias for gates (Wx + Rh + b) [H*4]
    const T* h,  // Recurrent state [N,H]
    const T* c,  // Cell state [N,H]
    T* h_out,    // Output recurrent state [N,H]
    T* c_out,    // Output cell state [N,H]
    T* v,        // Output vector (Wx + Rh + b) [N,H*4]
    T* tmp_Rh,   // Temporary storage for Rh vector [N,H*4]
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const bool training = data_->training;
  const int batch_size = data_->batch_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];
  const cudaEvent_t event = data_->event;

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, batch_size, hidden_size,
      &alpha,
      R, hidden_size * 4,
      h, hidden_size,
      &beta,
      tmp_Rh, hidden_size * 4);

  cudaStreamWaitEvent(stream1, event, 0);

  // Compute launch configuration for pointwise operations kernel.
  const dim3 blockDim(64, 16);
  const dim3 gridDim(
      (hidden_size + blockDim.x - 1) / blockDim.x,
      (batch_size + blockDim.y - 1) / blockDim.y);

  if (training) {
    if (zoneout_prob && zoneout_mask) {
      PointwiseOperations<T, true, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          v,
          tmp_Rh,
          b,
          h,
          c,
          h_out,
          c_out,
          v,
          zoneout_prob,
          zoneout_mask);
    } else {
      PointwiseOperations<T, true, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          v,
          tmp_Rh,
          b,
          h,
          c,
          h_out,
          c_out,
          v,
          0.0f,
          nullptr);
    }
  } else {
    if (zoneout_prob && zoneout_mask) {
      PointwiseOperations<T, false, true><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          v,
          tmp_Rh,
          b,
          h,
          c,
          h_out,
          c_out,
          nullptr,
          zoneout_prob,
          zoneout_mask);
    } else {
      PointwiseOperations<T, false, false><<<gridDim, blockDim, 0, stream1>>>(
          batch_size,
          hidden_size,
          v,
          tmp_Rh,
          b,
          h,
          c,
          h_out,
          c_out,
          nullptr,
          0.0f,
          nullptr);
    }
  }
}

template<typename T>
void ForwardPass<T>::Run(
    const int steps,
    const T* W,  // Weight matrix for input (Wx) [C,H*4]
    const T* R,  // Weight matrix for recurrent state (Rh) [H,H*4]
    const T* b,  // Bias for gates (Wx + Rh + b) [H*4]
    const T* x,  // Input vector [T,N,C]
    T* h,        // Recurrent state [T+1,N,H]
    T* c,        // Cell state [T+1,N,H]
    T* v,        // Output vector (Wx + Rh + b) [T,N,H*4]
    T* tmp_Rh,   // Temporary storage for Rh vector [N,H*4]
    const float zoneout_prob,
    const T* zoneout_mask) { // Zoneout mask [T,N,H]
  static const T alpha = static_cast<T>(1.0);
  static const T beta = static_cast<T>(0.0);

  const blas<void>::set_pointer_mode scoped1(data_->blas_handle);

  const int batch_size = data_->batch_size;
  const int input_size = data_->input_size;
  const int hidden_size = data_->hidden_size;
  const cublasHandle_t blas_handle = data_->blas_handle;
  const cudaStream_t stream1 = data_->stream[0];

  cudaStream_t save_stream;
  cublasGetStream(blas_handle, &save_stream);

  cublasSetStream(blas_handle, stream1);
  blas<T>::gemm(blas_handle,
      CUBLAS_OP_N, CUBLAS_OP_N,
      hidden_size * 4, steps * batch_size, input_size,
      &alpha,
      W, hidden_size * 4,
      x, input_size,
      &beta,
      v, hidden_size * 4);

  for (int i = 0; i < steps; ++i) {
    const int NH = batch_size * hidden_size;
    IterateInternal(
        R,
        b,
        h + i * NH,
        c + i * NH,
        h + (i + 1) * NH,
        c + (i + 1) * NH,
        v + i * NH * 4,
        tmp_Rh,
        zoneout_prob,
        zoneout_mask ? zoneout_mask + i * NH : nullptr);
  }

  cublasSetStream(blas_handle, save_stream);
}

template struct ForwardPass<float>;
template struct ForwardPass<double>;

}  // namespace lstm
}  // namespace v0
}  // namespace haste