
/* C++ wrapper for real-to-real transforms of FFTW 3.
 * For example, type 1 discrete cosine transform (DCT):
 * Initialize an object
 *   r2r::r2r dct1(FFTW_REDFT00);
 * Transform with
 *   Y = dct1(X);
 * or
 *   dct1.inplace(X,Y);
 * GTC – 1/10/21
*/
#pragma once
#include <fftw3.h>
#include <unordered_map>
#include <cstring>

// define inverse kinds
const fftw_r2r_kind INV_FFTW_R2HC     =  FFTW_HC2R;
const fftw_r2r_kind INV_FFTW_HC2R     =  FFTW_R2HC;
const fftw_r2r_kind INV_FFTW_DHT      =  FFTW_DHT;
const fftw_r2r_kind INV_FFTW_REDFT00  =  FFTW_REDFT00;
const fftw_r2r_kind INV_FFTW_REDFT10  =  FFTW_REDFT01;
const fftw_r2r_kind INV_FFTW_REDFT01  =  FFTW_REDFT10;
const fftw_r2r_kind INV_FFTW_REDFT11  =  FFTW_REDFT11;
const fftw_r2r_kind INV_FFTW_RODFT00  =  FFTW_RODFT00;
const fftw_r2r_kind INV_FFTW_RODFT10  =  FFTW_RODFT01;
const fftw_r2r_kind INV_FFTW_RODFT01  =  FFTW_RODFT10;
const fftw_r2r_kind INV_FFTW_RODFT11  =  FFTW_RODFT11;

namespace r2r {

  class r2r {
    private:
      std::unordered_map<int, fftw_plan> plans;
      void initialize_plan(int);
      fftw_r2r_kind kind;
    public:
      r2r(fftw_r2r_kind X=FFTW_REDFT00){kind=X;};
      ~r2r();
      template <class vec_t>
      void inplace(vec_t const &, vec_t &);
      template <class vec_t>
      vec_t operator()(vec_t const &);
  };

  r2r::~r2r() {
    for (auto plan : plans) fftw_destroy_plan(plan.second);
    plans.clear();
  }
  
  void r2r::initialize_plan(int N) {
    if (!(plans.count(N))) {
      #pragma omp critical
      if (!(plans.count(N))) {
      double* input;
      input = (double*) fftw_malloc(sizeof(double) * N);
      double* output;
      output = (double*) fftw_malloc(sizeof(double) * N);
      plans[N] = fftw_plan_r2r_1d(N, input, output, kind, FFTW_MEASURE);
      fftw_free(input);
      fftw_free(output);
      //std::cerr << "FFTW planning for size: " << N << std::endl;
      }
    }
  }
  
  template <class vec_t>
  void r2r::inplace(vec_t const &X, vec_t &ans) {
    int N = X.size();
    initialize_plan(N);
    double* input;
    input = (double*) fftw_malloc(sizeof(double) * N);
    double* output;
    output = (double*) fftw_malloc(sizeof(double) * N);
    memcpy(input, &X[0], sizeof(double) * N);
    fftw_execute_r2r(plans[N], input, output);
    memcpy(&ans[0], output, sizeof(double) * N);
    fftw_free(input);
    fftw_free(output);
  }

  template <class vec_t>
  vec_t r2r::operator()(vec_t const &X) {
    vec_t ans = X;
    inplace(X,ans);
    return ans;
  }

}
