#pragma once
#include <assert.h>
#include <mkl_lapacke.h>
#include "cppmkl/cppmkl_type_utils.h"

namespace cppmkl
{

  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const float *X, const lapack_int lda, lapack_int *ipiv)
  {
      return LAPACKE_sgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const MKL_Complex8 *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LAPACKE_cgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const std::complex<float> *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LAPACKE_cgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const double *X, const lapack_int lda, lapack_int *ipiv)
  { 
    return LAPACKE_dgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const MKL_Complex16 *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LAPACKE_zgetrf(layout, M, N, X, lda, ipiv);
  }
  inline lapack_int lapacke_getrf(int layout, lapack_int M, lapack_int N, const std::complex<double> *X, const lapack_int lda, lapack_int *ipiv)
  {
    return LAPACKE_zgetrf(layout, M, N, X, lda, ipiv);
  }

  //always return a double, auto conversion to float in calling code
  template <typename MATRIX_T>
  lapack_int LAPACKE_getrf(int layout, lapack_int M, lapack_int N, const MATRIX_T& X, const lapack_int lda, lapack_int *ipiv)
  {
    return lapacke_getrf(layout,M,N,X.data(),lda,ipiv);
  }
