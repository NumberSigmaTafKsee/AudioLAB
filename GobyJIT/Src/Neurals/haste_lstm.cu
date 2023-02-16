#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <string>
#include <unsupported/Eigen/CXX11/Tensor>
#include <vector>


#include "haste.hpp"

#include <cassert>
#include <cuda.h>
#include <cuda_runtime_api.h>

template<typename T>
struct device_ptr {
  static constexpr size_t ElemSize = sizeof(typename T::Scalar);

  static device_ptr<T> NewByteSized(size_t bytes) {
    return device_ptr<T>((bytes + ElemSize - 1) / ElemSize);
  }

  explicit device_ptr(size_t size_)
      : data(nullptr), size(size_) {
    void* tmp;
    cudaMalloc(&tmp, size * ElemSize);
    data = static_cast<typename T::Scalar*>(tmp);
  }

  explicit device_ptr(const T& elem)
      : data(nullptr), size(elem.size()) {
    void* tmp;
    cudaMalloc(&tmp, size * ElemSize);
    data = static_cast<typename T::Scalar*>(tmp);
    ToDevice(elem);
  }

  device_ptr(device_ptr<T>&& other) : data(other.data), size(other.size) {
    other.data = nullptr;
    other.size = 0;
  }

  device_ptr& operator=(const device_ptr<T>&& other) {
    if (&other != this) {
      data = other.data;
      size = other.size;
      other.data = nullptr;
      other.size = 0;
    }
    return *this;
  }

  device_ptr(const device_ptr<T>& other) = delete;
  device_ptr& operator=(const device_ptr<T>& other) = delete;

  void ToDevice(const T& src) {
    assert(size == src.size());
    cudaMemcpy(data, src.data(), src.size() * ElemSize, cudaMemcpyHostToDevice);
  }

  void ToHost(T& target) const {
    assert(size == target.size());
    cudaMemcpy(target.data(), data, target.size() * ElemSize, cudaMemcpyDeviceToHost);
  }

  size_t Size() const {
    return size;
  }

  void zero() {
    cudaMemset(data, 0, size * ElemSize);
  }

  ~device_ptr() {
    cudaFree(data);
  }

  typename T::Scalar* data;
  size_t size;
};

using haste::v0::lstm::BackwardPass;
using haste::v0::lstm::ForwardPass;
using std::string;

using Tensor1 = Eigen::Tensor<float, 1>;
using Tensor2 = Eigen::Tensor<float, 2>;
using Tensor3 = Eigen::Tensor<float, 3>;

constexpr int BATCH_SIZE = 64;
constexpr int SEQUENCE_LEN = 1000;
constexpr int HIDDEN_DIMS = 512;
constexpr int INPUT_DIMS = 512;

static cublasHandle_t g_blas_handle;

class ScopeTimer {
  public:
    ScopeTimer(const string& msg) : msg_(msg) {
      cudaEventCreate(&start_);
      cudaEventCreate(&stop_);
      cudaDeviceSynchronize();
      cudaEventRecord(start_);
    }

    ~ScopeTimer() {
      float elapsed_ms;
      cudaEventRecord(stop_);
      cudaEventSynchronize(stop_);
      cudaEventElapsedTime(&elapsed_ms, start_, stop_);
      printf("%s %.1fms\n", msg_.c_str(), elapsed_ms);
      cudaEventDestroy(start_);
      cudaEventDestroy(stop_);
    }

  private:
    string msg_;
    cudaEvent_t start_, stop_;
};

void LstmInference(const Tensor2& W, const Tensor2& R, const Tensor1& b, const Tensor3& x) {
  const int time_steps = x.dimension(2);
  const int batch_size = x.dimension(1);
  const int input_size = x.dimension(0);
  const int hidden_size = R.dimension(1);

  // Copy weights over to GPU.
  device_ptr<Tensor2> W_dev(W);
  device_ptr<Tensor2> R_dev(R);
  device_ptr<Tensor1> b_dev(b);
  device_ptr<Tensor3> x_dev(x);

  device_ptr<Tensor2> h_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor2> c_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> v_dev(time_steps * batch_size * hidden_size * 4);
  device_ptr<Tensor2> tmp_Rh_dev(batch_size * hidden_size * 4);

  h_dev.zero();
  c_dev.zero();

  ScopeTimer t("Inference:");

  ForwardPass<float> forward(
      false,  // training
      batch_size,
      input_size,
      hidden_size,
      g_blas_handle);

  forward.Run(
      time_steps,
      W_dev.data,
      R_dev.data,
      b_dev.data,
      x_dev.data,
      h_dev.data,
      c_dev.data,
      v_dev.data,
      tmp_Rh_dev.data,
      0.0f,      // zoneout prob
      nullptr);  // zoneout mask
}

void LstmTrain(const Tensor2& W, const Tensor2& R, const Tensor1& b, const Tensor3& x,
               const Tensor3& dh, const Tensor3& dc) {
  const int time_steps = x.dimension(2);
  const int batch_size = x.dimension(1);
  const int input_size = x.dimension(0);
  const int hidden_size = R.dimension(1);

  // Copy weights over to GPU.
  device_ptr<Tensor2> W_dev(W);
  device_ptr<Tensor2> R_dev(R);
  device_ptr<Tensor1> b_dev(b);
  device_ptr<Tensor3> x_dev(x);

  // This is nearly the same as the inference code except we have an extra dimension
  // for h and c. We'll store those outputs of the cell for all time steps and use
  // them during the backward pass below.
  device_ptr<Tensor3> h_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> c_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> v_dev(batch_size * hidden_size * 4 * time_steps);
  device_ptr<Tensor2> tmp_Rh_dev(batch_size * hidden_size * 4);

  h_dev.zero();
  c_dev.zero();

  {
    ScopeTimer t("Train forward:");
    ForwardPass<float> forward(
        true,  // training
        batch_size,
        input_size,
        hidden_size,
        g_blas_handle);

    forward.Run(
        time_steps,
        W_dev.data,
        R_dev.data,
        b_dev.data,
        x_dev.data,
        h_dev.data,
        c_dev.data,
        v_dev.data,
        tmp_Rh_dev.data,
        0.0f,      // zoneout prob
        nullptr);  // zoneout mask
  }

  Eigen::array<int, 3> transpose_x({ 1, 2, 0 });
  Tensor3 x_t = x.shuffle(transpose_x);

  Eigen::array<int, 2> transpose({ 1, 0 });
  Tensor2 W_t = W.shuffle(transpose);
  Tensor2 R_t = R.shuffle(transpose);

  device_ptr<Tensor3> x_t_dev(x_t);
  device_ptr<Tensor2> W_t_dev(W_t);
  device_ptr<Tensor2> R_t_dev(R_t);

  // These gradients should actually come "from above" but we're just allocating
  // a bunch of uninitialized memory and passing it in.
  device_ptr<Tensor3> dh_new_dev(dh);
  device_ptr<Tensor3> dc_new_dev(dc);

  device_ptr<Tensor3> dx_dev(time_steps * batch_size * input_size);
  device_ptr<Tensor2> dW_dev(input_size * hidden_size * 4);
  device_ptr<Tensor2> dR_dev(hidden_size * hidden_size * 4);
  device_ptr<Tensor2> db_dev(hidden_size * 4);
  device_ptr<Tensor2> dh_dev(batch_size * hidden_size);
  device_ptr<Tensor2> dc_dev(batch_size * hidden_size);

  dW_dev.zero();
  dR_dev.zero();
  db_dev.zero();
  dh_dev.zero();
  dc_dev.zero();

  {
    ScopeTimer t("Train backward:");
    BackwardPass<float> backward(
        batch_size,
        input_size,
        hidden_size,
        g_blas_handle);

    backward.Run(
        time_steps,
        W_t_dev.data,
        R_t_dev.data,
        b_dev.data,
        x_t_dev.data,
        h_dev.data,
        c_dev.data,
        dh_new_dev.data,
        dc_new_dev.data,
        dx_dev.data,
        dW_dev.data,
        dR_dev.data,
        db_dev.data,
        dh_dev.data,
        dc_dev.data,
        v_dev.data,
        nullptr);
  }
}

void LstmTrainIterative(const Tensor2& W, const Tensor2& R, const Tensor1& b, const Tensor3& x,
                        const Tensor3& dh, const Tensor3& dc) {
  const int time_steps = x.dimension(2);
  const int batch_size = x.dimension(1);
  const int input_size = x.dimension(0);
  const int hidden_size = R.dimension(1);

  // Copy weights over to GPU.
  device_ptr<Tensor2> W_dev(W);
  device_ptr<Tensor2> R_dev(R);
  device_ptr<Tensor1> b_dev(b);
  device_ptr<Tensor3> x_dev(x);

  device_ptr<Tensor3> h_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> c_dev((time_steps + 1) * batch_size * hidden_size);
  device_ptr<Tensor3> v_dev(time_steps * batch_size * hidden_size * 4);
  device_ptr<Tensor2> tmp_Rh_dev(batch_size * hidden_size * 4);

  h_dev.zero();
  c_dev.zero();

  {
    ScopeTimer t("Train forward (iterative):");
    ForwardPass<float> forward(
        true,  // training
        batch_size,
        input_size,
        hidden_size,
        g_blas_handle);

    const int NC = batch_size * input_size;
    const int NH = batch_size * hidden_size;
    for (int t = 0; t < time_steps; ++t) {
      forward.Iterate(
          0,
          W_dev.data,
          R_dev.data,
          b_dev.data,
          x_dev.data + t * NC,
          h_dev.data + t * NH,
          c_dev.data + t * NH,
          h_dev.data + (t + 1) * NH,
          c_dev.data + (t + 1) * NH,
          v_dev.data + t * NH * 4,
          tmp_Rh_dev.data,
          0.0f,      // zoneout prob
          nullptr);  // zoneout mask
    }
  }

  Eigen::array<int, 3> transpose_x({ 1, 2, 0 });
  Tensor3 x_t = x.shuffle(transpose_x);

  Eigen::array<int, 2> transpose({ 1, 0 });
  Tensor2 W_t = W.shuffle(transpose);
  Tensor2 R_t = R.shuffle(transpose);

  device_ptr<Tensor3> x_t_dev(x_t);
  device_ptr<Tensor2> W_t_dev(W_t);
  device_ptr<Tensor2> R_t_dev(R_t);

  // These gradients should actually come "from above" but we're just allocating
  // a bunch of uninitialized memory and passing it in.
  device_ptr<Tensor3> dh_new_dev(dh);
  device_ptr<Tensor3> dc_new_dev(dc);

  device_ptr<Tensor3> dx_dev(time_steps * batch_size * input_size);
  device_ptr<Tensor2> dW_dev(input_size * hidden_size * 4);
  device_ptr<Tensor2> dR_dev(hidden_size * hidden_size * 4);
  device_ptr<Tensor2> db_dev(hidden_size * 4);
  device_ptr<Tensor2> dh_dev(batch_size * hidden_size);
  device_ptr<Tensor2> dc_dev(batch_size * hidden_size);

  dW_dev.zero();
  dR_dev.zero();
  db_dev.zero();
  dh_dev.zero();
  dc_dev.zero();

  {
    ScopeTimer t("Train backward (iterative):");
    BackwardPass<float> backward(
        batch_size,
        input_size,
        hidden_size,
        g_blas_handle);

    const int NC = batch_size * input_size;
    const int NH = batch_size * hidden_size;
    for (int t = time_steps - 1; t >= 0; --t) {
      backward.Iterate(
          0,
          W_t_dev.data,
          R_t_dev.data,
          b_dev.data,
          x_t_dev.data + t * NC,
          h_dev.data + t * NH,
          c_dev.data + t * NH,
          c_dev.data + (t + 1) * NH,
          dh_new_dev.data + t * NH,
          dc_new_dev.data + t * NH,
          dx_dev.data + t * NC,
          dW_dev.data,
          dR_dev.data,
          db_dev.data,
          dh_dev.data,
          dc_dev.data,
          v_dev.data + t * NH * 4,
          nullptr);
    }
  }
}

int main() {
  srand(time(0));

  cublasCreate(&g_blas_handle);

  // Weights.
  // W: input weight matrix
  // R: recurrent weight matrix
  // b: bias
  Tensor2 W(HIDDEN_DIMS * 4, INPUT_DIMS);
  Tensor2 R(HIDDEN_DIMS * 4, HIDDEN_DIMS);
  Tensor1 b(HIDDEN_DIMS * 4);

  // Input.
  Tensor3 x(INPUT_DIMS, BATCH_SIZE, SEQUENCE_LEN);

  // Gradients from upstream layers.
  Tensor3 dh(HIDDEN_DIMS, BATCH_SIZE, SEQUENCE_LEN + 1);
  Tensor3 dc(HIDDEN_DIMS, BATCH_SIZE, SEQUENCE_LEN + 1);

  W.setRandom();
  R.setRandom();
  b.setRandom();
  x.setRandom();
  dh.setRandom();
  dc.setRandom();

  LstmInference(W, R, b, x);
  LstmTrain(W, R, b, x, dh, dc);
  LstmTrainIterative(W, R, b, x, dh, dc);

  cublasDestroy(g_blas_handle);

  return 0;
}

