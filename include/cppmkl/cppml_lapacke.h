#pragma once
#include <assert.h>
#include <mkl_lapacke.h>
#include "cppmkl/cppmkl_type_utils.h"

namespace cppmkl
{

  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const float *X, const lapack_int lda, lapack_int *ipiv)
  {
      return LACPACK_sgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const MKL_Complex8 *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LACPACK_cgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const std::complex<float> *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LACPACK_cgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const double *X, const lapack_int lda, lapack_int *ipiv)
  { 
    return LACPACK_dgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const MKL_Complex16 *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LACPACK_zgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const std::complex<double> *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LACPACK_zgetrf(layout, M, N, X, lda, ipiv);
  }

  //always return a double, auto conversion to float in calling code
  template <typename VECTOR_T>
  double cblas_asum(const VECTOR_T& v, const MKL_INT incX=1)
  {
    return cblas_asum(v.size()/incX, ptr_to_first(v), incX);
  }
