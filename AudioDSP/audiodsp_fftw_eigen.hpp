#pragma once

#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <fstream>
#include <tuple>
#include <iomanip>
#include <vector>
#include <cmath>
#include <fftw3.h>

#define EIGEN_DEFAULT_TO_ROW_MAJOR
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigen>


namespace EigenFFT
{
	/*
	 * FFTW3 Function reference:
	 * 
	 * Complex DFT
	 * 
	 *  fftw_plan fftw_plan_dft_1d (int           n,
	 *                              fftw_complex* in,   fftw_complex* out,
	 *                              int           sign, unsigned      flags);
	 * 
	 *  fftw_plan fftw_plan_dft_2d (int           n0,   int           n1,
	 *                              fftw_complex* in,   fftw_complex* out,
	 *                              int           sign, unsigned      flags);
	 * 
	 *  fftw_plan fftw_plan_dft_3d (int           n0,   int           n1,        int           n2,
	 *                              fftw_complex* in,   fftw_complex* out,
	 *                              int           sign, unsigned      flags);
	 * 
	 *  fftw_plan fftw_plan_dft    (int           rank, const int*    n,
	 *                              fftw_complex* in,   fftw_complex* out,
	 *                              int           sign, unsigned      flags);
	 * 
	 * 
	 * Real DFT
	 * 
	 *  fftw_plan fftw_plan_dft_r2c_1d(int n,
	 *                                 double*  in,   fftw_complex* out,
	 *                                 unsigned flags                  );
	 * 
	 *  fftw_plan fftw_plan_dft_r2c_2d(int      n0,   int           n1,
	 *                                 double*  in,   fftw_complex* out,
	 *                                 unsigned flags                  );
	 * 
	 *  fftw_plan fftw_plan_dft_r2c_3d(int      n0,   int           n1, int         n2,
	 *                                 double*  in,   fftw_complex* out,
	 *                                 unsigned flags                  );
	 * 
	 *  fftw_plan fftw_plan_dft_r2c   (int      rank, const int*    n,
	 *                                 double*  in,   fftw_complex* out,
	 *                                 unsigned flags                  );
	 * 
	 * 
	 * Real-to-real DFT
	 * 
	 *  fftw_plan fftw_plan_r2r_1d(int           n, 
	 *                             double*       in,    double*       out,
	 *                             fftw_r2r_kind kind,               
	 *                             unsigned      flags                  );
	 *                                                               
	 *  fftw_plan fftw_plan_r2r_2d(int           n0,    int           n1, 
	 *                             double*       in,    double*       out,
	 *                             fftw_r2r_kind kind0, fftw_r2r_kind kind1, 
	 *                             unsigned      flags);
	 * 
	 *  fftw_plan fftw_plan_r2r_3d(int           n0,     int       n1, int       n2,
	 *                             double*       in,     double*   out,
	 *                             fftw_r2r_kind kind0,
	 *                             fftw_r2r_kind kind1,
	 *                             fftw_r2r_kind kind2,
	 *                             unsigned      flags);
	 * 
	 *  fftw_plan fftw_plan_r2r   (int           rank, const int*      n, 
	 *                             double*       in,   double*         out, 
	 *                             const fftw_r2r_kind*  kind, 
	 *                             unsigned      flags);
	 * 
	 * 
	 * For more detailed info: https://www.cnblogs.com/aiguona/p/9407425.html
	 */




	namespace ionizing {
	  using namespace Eigen;

	  static_assert(sizeof(size_t) == sizeof(uint64_t),
	      "Uncompatible compiler: sizeof(size_t) != 8 !");
	  static_assert(sizeof(double) == 8,
	      "Uncompatible compiler: sizeof(double) != 8 !");
	  static_assert(sizeof(float)  == 4,
	      "Uncompatible compiler: sizeof( float) != 4 !");
	  static_assert(sizeof(int)    == 4,
	      "Uncompatible compiler: sizeof(   int) != 4 !");
	  static_assert(sizeof(std::complex<double>) == 16, 
	      "Uncompatible compiler: sizeof(complex<double>) != 16 !");
	  static_assert(sizeof(std::complex< float>) ==  8, 
	      "Uncompatible compiler: sizeof(complex< float>) !=  8 !");

	// global type aliasing
	  using Vecd   = VectorXd;
	  using Veccd  = VectorXcd;
	  using Matd   = Matrix<              double, Dynamic, Dynamic, RowMajor>;
	  using Matcd  = Matrix<std::complex<double>, Dynamic, Dynamic, RowMajor>;
	  using Cubd   = Tensor<              double,                3, RowMajor>;
	  using Cubcd  = Tensor<std::complex<double>,                3, RowMajor>;
	  using Mat33d = Matrix<double,                     3,       3, RowMajor>;
	  using MatX3d = Matrix<double,               Dynamic,       3, RowMajor>;

	  template <typename T>
	    using MatT    = Matrix<T, Dynamic, Dynamic, RowMajor>;

	  template <typename T>
	    using ColVecT = Matrix<T, Dynamic,       1, RowMajor>;

	  template <typename T>
	    using RowVecT = Matrix<T,       1, Dynamic, RowMajor>;
	  
	namespace FFTW {
	  /*
	   * using Vecd   = VectorXd;
	   * using Veccd  = VectorXcd;
	   * using Matd   = Matrix<              double, Dynamic, Dynamic>;
	   * using Matcd  = Matrix<std::complex<double>, Dynamic, Dynamic>;
	   * using Cubd   = Tensor<              double,                3, RowMajor>;
	   * using Cubcd  = Tensor<std::complex<double>,                3, RowMajor>;
	   * using Mat33d = Matrix<double,                     3,       3>;
	   * using MatX3d = Matrix<double,               Dynamic,       3>;
	   */

	  Veccd  fft_1d(const Veccd& vec);
	  Veccd  fft_1d(const  Vecd& vec);
	  Veccd ifft_1d(const Veccd& vec);
	  Veccd ifft_1d(const  Vecd& vec);

	  Matcd  fft_2d(const Matcd& mat);
	  Matcd  fft_2d(const  Matd& mat);
	  Matcd ifft_2d(const Matcd& mat);
	  Matcd ifft_2d(const  Matd& mat);

	  Cubcd  fft_3d(const Cubcd& cub);
	  Cubcd  fft_3d(const  Cubd& cub);
	  Cubcd ifft_3d(const Cubcd& cub);
	  Cubcd ifft_3d(const  Cubd& cub);

	  Veccd  rfft_1d(const  Vecd& vec);
	  Vecd  irfft_1d(       Veccd vec);
		   
	  Matcd  rfft_2d(const  Matd& mat);
	  Matd  irfft_2d(       Matcd mat);
		   
	  Cubcd  rfft_3d(const  Cubd& cub);
	  Cubd  irfft_3d(       Cubcd cub);
	}

	#include <fft.hpp>

	static IOFormat CommaInitFmt{StreamPrecision, DontAlignCols, ", ", ", ", "", "", " << ", ";"};
	static IOFormat CleanFmt{4, 0, ", ", "\n", "[", "]"};
	static IOFormat OctaveFmt{StreamPrecision, 0, ", ", ";\n", "", "", "[", "]"};
	static IOFormat HeavyFmt{FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]"};

	namespace FFTW {
	  /*
	   * using Vecd   = VectorXd;
	   * using Veccd  = VectorXcd;
	   * using Matd   = Matrix<              double, Dynamic, Dynamic>;
	   * using Matcd  = Matrix<std::complex<double>, Dynamic, Dynamic>;
	   * using Cubd   = Tensor<              double,                3, RowMajor>;
	   * using Cubcd  = Tensor<std::complex<double>,                3, RowMajor>;
	   * using Mat33d = Matrix<double,                     3,       3>;
	   * using MatX3d = Matrix<double,               Dynamic,       3>;
	   */

	  Veccd fft_1d(const Veccd& vec){
	    if (vec.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input array is 0" << std::endl;
	      std::abort();
	    }
	    Veccd out(vec.size());
	    auto plan = fftw_plan_dft_1d(vec.size(),
		(fftw_complex *)vec.data(), (fftw_complex *)out.data(),
		FFTW_FORWARD, FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    return out;
	  }

	  Veccd fft_1d(const Vecd& vec){
	    Veccd tmp = vec.cast<std::complex<double>>();
	    // std::cout << __FUNCTION__ << ": tmp = \n" << tmp << std::endl;
	    return fft_1d(tmp);
	  }

	  Veccd ifft_1d(const Veccd& vec) {
	    if (vec.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input array is 0" << std::endl;
	      std::abort();
	    }
	    Veccd out(vec.size());
	    auto plan = fftw_plan_dft_1d(vec.size(),
		(fftw_complex *)vec.data(), (fftw_complex *)out.data(),
		FFTW_BACKWARD, FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    out /= out.size();
	    return out;
	  }

	  Veccd ifft_1d(const Vecd& vec){
	    Veccd tmp = vec.cast<std::complex<double>>();
	    // std::cout << __FUNCTION__ << ": tmp = \n" << tmp << std::endl;
	    return ifft_1d(tmp);
	  }


	  Matcd  fft_2d(const Matcd& mat) {
	    if (mat.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input matrix is 0" << std::endl;
	      std::abort();
	    }
	    Matcd out(mat.rows(), mat.cols());
	    auto plan = fftw_plan_dft_2d(mat.rows(), mat.cols(),
		(fftw_complex *)mat.data(), (fftw_complex *)out.data(),
		FFTW_FORWARD, FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    return out;
	  }



	  Matcd  fft_2d(const  Matd& mat) {
	    Matcd tmp = mat.cast<std::complex<double>>();
	    // std::cout << __FUNCTION__ << ": tmp = \n" << tmp << std::endl;
	    return fft_2d(tmp);
	  }
	  
	  Matcd ifft_2d(const Matcd& mat) {
	    if (mat.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input matrix is 0" << std::endl;
	      std::abort();
	    }
	    Matcd out(mat.rows(), mat.cols());
	    auto plan = fftw_plan_dft_2d(mat.rows(), mat.cols(),
		(fftw_complex *)mat.data(), (fftw_complex *)out.data(),
		FFTW_BACKWARD, FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    out /= out.size();
	    return out;
	  }

	  Matcd ifft_2d(const  Matd& mat) {
	    Matcd tmp = mat.cast<std::complex<double>>();
	    // std::cout << __FUNCTION__ << ": tmp = \n" << tmp << std::endl;
	    return ifft_2d(tmp);
	  }

	  Cubcd  fft_3d(const Cubcd& cub) {
	    if (cub.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input cub is 0" << std::endl;
	      std::abort();
	    }

	    Cubcd out(cub.dimension(0), cub.dimension(1), cub.dimension(2));
	    auto plan = fftw_plan_dft_3d(
		cub.dimension(0), cub.dimension(1), cub.dimension(2), 
		(fftw_complex *)cub.data(), (fftw_complex *)out.data(),
		FFTW_FORWARD, FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    return out;
	  }

	  Cubcd  fft_3d(const  Cubd& cub) {
	    Cubcd tmp = cub.cast<std::complex<double>>();
	    // std::cout << __FUNCTION__ << ": tmp = \n" << tmp << std::endl;
	    return fft_3d(tmp);
	  }

	  Cubcd ifft_3d(const Cubcd& cub) {
	    if (cub.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input cub is 0" << std::endl;
	      std::abort();
	    }
	    Cubcd out(cub.dimension(0), cub.dimension(1), cub.dimension(2));
	    auto plan = fftw_plan_dft_3d(
		cub.dimension(0), cub.dimension(1), cub.dimension(2), 
		(fftw_complex *)cub.data(), (fftw_complex *)out.data(),
		FFTW_BACKWARD, FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    // std::cout << __FUNCTION__ << "out.size() == " << out.size() << std::endl;
	    out /= out.constant(out.size());
	    return out;
	  }

	  Cubcd ifft_3d(const  Cubd& cub) {
	    Cubcd tmp = cub.cast<std::complex<double>>();
	    // std::cout << __FUNCTION__ << ": tmp = \n" << tmp << std::endl;
	    return ifft_3d(tmp);
	  }


	}


	/************************************************/
	/************* For real numbers *****************/
	/************************************************/

	  /*
	   * using Vecd   = VectorXd;
	   * using Veccd  = VectorXcd;
	   * using Matd   = Matrix<              double, Dynamic, Dynamic>;
	   * using Matcd  = Matrix<std::complex<double>, Dynamic, Dynamic>;
	   * using Cubd   = Tensor<              double,                3, RowMajor>;
	   * using Cubcd  = Tensor<std::complex<double>,                3, RowMajor>;
	   * using Mat33d = Matrix<double,                     3,       3>;
	   * using MatX3d = Matrix<double,               Dynamic,       3>;
	   */

	namespace ionizing {

	/*
	 * For real to complex FFT, out array only need to have a
	 * size of [in.size() / 2 + 1].
	 *  in array shape: m
	 * out array shape: m/2 + 1
	 */
	  Veccd rfft_1d(const Vecd& vec){
	    if (vec.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __FUNCTION__ 
		<< ": " << "size of input vector is 0" << std::endl;
	      std::abort();
	    }

	    Veccd out(vec.size() / 2 + 1);
	    auto plan = fftw_plan_dft_r2c_1d(vec.size(),
		(double *)vec.data(), (fftw_complex *)out.data(),
		FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    return out;
	  }

	/*
	 * fftw_c2r destroyes the original vector, so in this function
	 * you cannot pass `const T& in`
	 *  in array shape:  m
	 * out array shape:  (m-1) * 2
	 */
	  Vecd irfft_1d(Veccd vec) {
	    if (vec.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __FUNCTION__ 
		<< ": " << "size of input vector is 0" << std::endl;
	      std::abort();
	    }
	    // if original vector length is odd, size = in.size() * 2 - 1,
	    //                             else, size = in.size() * 2 - 2;
	    int size = fabs(vec.tail<1>()(0).imag()) > 1e-5 ?
	      vec.size() * 2 - 1 : vec.size() * 2 - 2 ;
	    Vecd out(size);
	    auto plan = fftw_plan_dft_c2r_1d(out.size(),
		(fftw_complex *)vec.data(), (double *)out.data(),
		FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    out /= out.size();
	    return out;
	  }



	/*
	 * For real to comlex 2D-FFT:
	 *  in array shape:  m x n,
	 * out array shape:  m x (n/2 + 1);
	 */
	  Matcd  rfft_2d(const Matd& mat) {
	    if (mat.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input matrix is 0" << std::endl;
	      std::abort();
	    }

	    Matcd out(mat.rows(), mat.cols() / 2 + 1);
	    auto plan = fftw_plan_dft_r2c_2d(mat.rows(), mat.cols(),
		(double *)mat.data(), (fftw_complex *)out.data(),
		FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    return out;
	  }
	  
	/*
	 * fftw_c2r destroyes the original matrix, you cannot pass
	 * `const T& in`
	 *  in array shape: m x n
	 * out array shape: m x ((n-1) * 2)
	 */
	  Matd irfft_2d(Matcd mat) {
	    if (mat.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __LINE__ << ": "
		<< "size of input matrix is 0" << std::endl;
	      std::abort();
	    }

	    int cols = fabs(mat.rightCols<1>().head<1>()(0).imag()) > 1e-5 ?
	      mat.cols() * 2 - 1 : mat.cols() * 2 - 2;
	    Matd out(mat.rows(), cols);
	    auto plan = fftw_plan_dft_c2r_2d(out.rows(), out.cols(),
		(fftw_complex *)mat.data(), (double *)out.data(),
		FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    out /= out.size();
	    return out;
	  }

	/*
	 * For real to comlex 3D-FFT:
	 *  in array shape:  m x n x q,
	 * out array shape:  m x n x (q/2 + 1);
	 */
	  Cubcd  rfft_3d(const Cubd& cub) {
	    if (cub.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __FUNCTION__ << ": "
		<< "size of input cub is 0" << std::endl;
	      std::abort();
	    }

	    Cubcd out(cub.dimension(0), cub.dimension(1), cub.dimension(2) / 2 + 1);
	    auto plan = fftw_plan_dft_r2c_3d(
		cub.dimension(0), cub.dimension(1), cub.dimension(2), 
		(double *)cub.data(), (fftw_complex *)out.data(),
		FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    return out;
	  }

	/*
	 *  in array shape: m x n x q
	 * out array shape: m x n x ((q-1) * 2)
	 */
	  Cubd irfft_3d(Cubcd cub) {
	    if (cub.size() == 0) {
	      std::cerr << std::endl << __FILE__ << ":" << __FUNCTION__ << ": "
		<< "size of input cub is 0" << std::endl;
	      std::abort();
	    }

	/*
	 * if in cube is [m x n x q]
	 *  tmp_last_elem = cube[0, 0, -1].imag()
	 *  this variable can be used to determine the
	 *  if the last dimension of original cube is
	 *  odd or even.
	 */
	    double tmp_last_elem = cub(0, 0, cub.dimension(2) - 1).imag();
	    int dim3 = fabs(tmp_last_elem) > 1e-5 ?
		cub.dimension(2) * 2 - 1 :
		cub.dimension(2) * 2 - 2 ;
	    Cubd out(cub.dimension(0), cub.dimension(1), dim3);
	    auto plan = fftw_plan_dft_c2r_3d(
		out.dimension(0), out.dimension(1), out.dimension(2), 
		(fftw_complex *)cub.data(), (double *)out.data(),
		FFTW_ESTIMATE);
	    fftw_execute(plan);
	    fftw_destroy_plan(plan);
	    out /= out.constant(out.size());
	    return out;
	  }
}

