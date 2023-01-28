#pragma once

//#define EIGEN2_SUPPORT

#define EIGEN_ARRAYBASE_PLUGIN "loudia_cwise_addons.h"
#define EIGEN_FUNCTORS_PLUGIN  "louda_functors_addons.h"

//#define EIGEN_DEFAULT_TO_ROW_MAJOR

#include "Debug.h"

#include <cmath>
#include <sstream>
#include <cassert>
#include <exception>
#include <cstring>
#include <limits>
#include <iostream>

#include <Eigen/Eigen>
#include <Eigen/LU>

#include <fftw3.h>
#include <samplerate.h>

using namespace std;
using namespace Eigen;

#define MAX_AUDIO_FRAME_SIZE 192000

#ifndef __STDC_CONSTANT_MACROS
#define __STDC_CONSTANT_MACROS
#endif

extern "C" {
    #ifndef __STDC_CONSTANT_MACROS
    #define __STDC_CONSTANT_MACROS
    #endif

    #ifdef LOUDIA_OLD_FFMPEG
    #include <ffmpeg/avcodec.h>
    #include <ffmpeg/avformat.h>
    #else
    #include <libavcodec/avcodec.h>
    #include <libavformat/avformat.h>
    #include <libavresample/avresample.h>
    #include <libavutil/opt.h>
    #include <libavutil/mathematics.h>
    #endif
}


namespace Loudia
{
    // Types for scalar values
    typedef int Integer;
    typedef float Real;
    typedef std::complex< Real > Complex;

    typedef std::complex< float > complex_float;
    typedef std::complex< double > complex_double;

    // Types for vector values
    typedef Eigen::Matrix< Integer, 1, Eigen::Dynamic > RowXI;
    typedef Eigen::Matrix< Real, 1, Eigen::Dynamic > RowXR;
    typedef Eigen::Matrix< Complex, 1, Eigen::Dynamic > RowXC;

    typedef Eigen::Matrix< Integer, Eigen::Dynamic, 1 > ColXI;
    typedef Eigen::Matrix< Real, Eigen::Dynamic, 1 > ColXR;
    typedef Eigen::Matrix< Complex, Eigen::Dynamic, 1 > ColXC;

    typedef ColXI VectorXI;
    typedef ColXR VectorXR;
    typedef ColXC VectorXC;

    // Types for matrix values
    typedef Eigen::Matrix< Integer, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXI;
    typedef Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXR;
    typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXC;

    // Types for mapping Scipy matrices (these are RowMajor)
    typedef Eigen::Matrix< Integer, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXIscipy;
    typedef Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXRscipy;
    typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXCscipy;

    typedef Eigen::Matrix< float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXfscipy;
    typedef Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXdscipy;

    typedef Eigen::Matrix< std::complex< float >, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXcfscipy;
    typedef Eigen::Matrix< std::complex< double >, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXcdscipy;

    template <typename T>
    inline T square(const T& val) {
    return val*val;
    }

    /**
    * Exception class for Loudia. It has a whole slew of different constructors
    * to make it as easy as possible to throw an exception with a descriptive
    * message.
    */
    class LoudiaException : public std::exception {

    public:
    LoudiaException(const char* msg) : exception(), _msg(msg) {}
    LoudiaException(const std::string& msg) : exception(), _msg(msg) {}
    LoudiaException(const std::ostringstream& msg) : exception(), _msg(msg.str()) {}

    template <typename T, typename U>
    LoudiaException(const T& a, const U& b) : exception() {
        std::ostringstream oss; oss << a << b; _msg = oss.str();
    }

    template <typename T, typename U, typename V>
    LoudiaException(const T& a, const U& b, const V& c) : exception() {
        std::ostringstream oss; oss << a << b << c; _msg = oss.str();
    }

    virtual ~LoudiaException() throw() {}
    virtual const char* what() const throw() { return _msg.c_str(); }

    protected:
    std::string _msg;

    };

    #define LOUDIA_ERROR(msg) ostringstream loudiaErrorMessage__;  loudiaErrorMessage__ << msg; throw LoudiaException(loudiaErrorMessage__);
    #define LOUDIA_WARNING(msg) std::cerr << msg << std::endl;


    // If none of the NO_DEBUG are defined we enable debugging
    #if (!defined(LOUDIA_NO_DEBUG) && !defined(NDEBUG))

    #include <iostream>
    #define LOUDIA_DEBUG(msg) std::cerr << msg << std::endl;

    // else we do nothing
    #else

    #define LOUDIA_DEBUG(msg)

    #endif

    /**
    * Given a matrix of polynomes (one per column)
    * returns a matrix of roots (a vector of roots per column)
    */
    void roots(const MatrixXR& poly, MatrixXC* result);

    /**
    * Given a matrix of roots (a vector of roots per column)
    * returns a matrix of polynomes (a polynome per vector of roots)
    */
    void poly(const MatrixXC& roots, MatrixXC* result);

    /**
    * Given two row matrices
    * returns the convolution of both
    */
    void convolve(const MatrixXC& a, const MatrixXC& b, MatrixXC* c);
    void convolve(const MatrixXR& a, const MatrixXR& b, MatrixXR* c);


    /**
    * Given two row matrices
    * returns the correlation of both
    */
    void correlate(const MatrixXC& a, const MatrixXC& b, MatrixXC* c,
                int _minlag = -std::numeric_limits<int>::infinity(),
                int _maxlag = std::numeric_limits<int>::infinity());

    void correlate(const MatrixXR& a, const MatrixXR& b, MatrixXR* c,
                int _minlag = -std::numeric_limits<int>::infinity(),
                int _maxlag = std::numeric_limits<int>::infinity());

    /**
    * Given a row matrix
    * returns the autocorrelation
    */
    void autocorrelate(const MatrixXR& a, MatrixXR* c,
                    int _minlag = 0,
                    int _maxlag = std::numeric_limits<int>::infinity());

    void autocorrelate(const MatrixXC& a, MatrixXC* c,
                    int _minlag = 0,
                    int _maxlag = std::numeric_limits<int>::infinity());


    /**
    * Reverse in place the order of the columns
    */
    void reverseCols(MatrixXC* in);
    void reverseCols(MatrixXR* in);


    /**
    * Calculate inplace the cumulative sum
    */
    void rowCumsum(MatrixXR* in);
    void colCumsum(MatrixXR* in);

    /**
    * Calculate inplace shift of a matrix
    */
    void rowShift(MatrixXR* in, int num);
    void colShift(MatrixXR* in, int num);

    /**
    * Calculate inplace range matrix
    */
    void range(Real start, Real end, int steps, MatrixXC* in);
    void range(Real start, Real end, int steps, int rows, MatrixXC* in);
    void range(Real start, Real end, int steps, MatrixXR* in);
    void range(Real start, Real end, int steps, int rows, MatrixXR* in);
    void range(Real start, Real end, int steps, MatrixXI* in);
    void range(Real start, Real end, int steps, int rows, MatrixXI* in);

    /**
    * Create a matrix of complex numbers given the polar coordinates
    */
    void polar(const MatrixXR& mag, const MatrixXR& phase, MatrixXC* complex);

    /**
    * Calculate the combinations of N elements in groups of k
    *
    */
    int combination(int N, int k);

    /**
    * Calculate the aliased cardinal sine defined as:
    *
    *   asinc(M, T, x) = sin(M * pi * x * T) / sin(pi * x * T)
    */
    Real asinc(int M, Real omega);

    /**
    * Calculate the Fourier transform of a hamming window
    */
    void raisedCosTransform(Real position, Real magnitude,
                            int windowSize, int fftSize,
                            Real alpha, Real beta,
                            MatrixXR* spectrum, int* begin, int* end, int bandwidth);

    void raisedCosTransform(Real position, Real magnitude,
                            int windowSize, int fftSize,
                            Real alpha, Real beta,
                            MatrixXR* spectrum, int bandwidth);

    void hannTransform(Real position, Real magnitude,
                    int windowSize, int fftSize,
                    MatrixXR* spectrum, int bandwidth = 4);

    void hannTransform(Real position, Real magnitude,
                    int windowSize, int fftSize,
                    MatrixXR* spectrum, int* begin, int* end, int bandwidth = 4);


    void hammingTransform(Real position, Real magnitude,
                        int windowSize, int fftSize,
                        MatrixXR* spectrum, int bandwidth = 4);

    void hammingTransform(Real position, Real magnitude,
                        int windowSize, int fftSize,
                        MatrixXR* spectrum, int* begin, int* end, int bandwidth = 4);

    void dbToMag(const MatrixXR& db, MatrixXR* mag);

    void magToDb(const MatrixXR& mag, MatrixXR* db, Real minMag = 0.0001 );

    void unwrap(const MatrixXR& phases, MatrixXR* unwrapped);

    void freqz(const MatrixXR& b, const MatrixXR& a, const MatrixXR& w, MatrixXC* resp);
    void freqz(const MatrixXR& b, const MatrixXR& w, MatrixXC* resp);

    void derivate(const MatrixXR& a, MatrixXR* b);

    int nextPowerOf2(Real a, int factor = 0);

    Real gaussian(Real x, Real mu, Real fi);

    void gaussian(Real x, MatrixXR mu, MatrixXR fi, MatrixXR* result);

    void pseudoInverse(const MatrixXR& a, MatrixXR* result, Real epsilon = 1e-6);

    /**
    * Create the zeros, poles and gain of an analog prototype of a Chebyshev Type I filter.
    */
    void chebyshev1(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

    /**
    * Create the  zeros, poles and gain of an analog prototype of a Chebyshev Type II filter.
    */
    void chebyshev2(int order, Real rippleDB, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

    /**
    * Create the zeros, poles and gain of an analog prototype of a Butterworth filter.
    */
    void butterworth(int order, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

    /**
    * Create the zeros, poles and gain of an analog prototype of a Bessel filter.
    */
    void bessel(int order, int channels, MatrixXC* zeros, MatrixXC* poles, Real* gain);

    /**
    * Convert from the b and a coefficients of an IIR filter to the
    * zeros, poles and gain of the filter
    */
    void coeffsToZpk(const MatrixXR& b, const MatrixXR& a, MatrixXC* zeros, MatrixXC* poles, Real* gain);

    /**
    * Convert from zeros, poles and gain of an IIR filter to the
    * the b and a coefficients of the filter
    */
    void zpkToCoeffs(const MatrixXC& zeros, const MatrixXC& poles, Real gain, MatrixXC*  b, MatrixXC*  a);

    /**
    * Convert from the b and a coefficients from low pass to low pass of an IIR filter with critical frequency 1.0
    * to the coefficients with the critical frequency passed as argument
    */
    void lowPassToLowPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout);

    /**
    * Convert from the b and a coefficients from low pass to high pass of an IIR filter with critical frequency 1.0
    * to the coefficients with the critical frequency passed as argument
    */
    void lowPassToHighPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout);

    /**
    * Convert from the b and a coefficients from low pass to band pass of an IIR filter with critical frequency 1.0
    * to the coefficients with the critical frequency passed as argument
    */
    void lowPassToBandPass(const MatrixXC& b, const MatrixXC& a, Real freq, Real freqStop, MatrixXC*  bout, MatrixXC*  aout);

    /**
    * Convert from the b and a coefficients from low pass to band stop of an IIR filter with critical frequency 1.0
    * to the coefficients with the critical frequency passed as argument
    */
    void lowPassToBandStop(const MatrixXC& b, const MatrixXC& a, Real freq, Real freqStop, MatrixXC*  bout, MatrixXC*  aout);

    /**
    * Normalize to a first coefficient
    * 
    */
    void normalize(MatrixXC& b, MatrixXC& a);

    /**
    * Apply the biliniear transformations to a set of coefficients
    * 
    */
    void bilinear(const MatrixXC& b, const MatrixXC& a, Real fs, MatrixXR*  bout, MatrixXR*  aout);
}

#include "loudia_aok.hpp"
#include "loudia_audioloader.hpp"
#include "loudia_window.hpp"
#include "loudia_spectrum.hpp"
#include "loudia_correlation.hpp"
#include "loudia_filters.hpp"    
#include "loudia_bands.hpp"    
#include "loudia_framecutter.hpp"    
#include "loudia_inmf.hpp"
#include "loudia_lpc.hpp"    
#include "loudia_mel.hpp"    
#include "loudia_nmf.hpp"
#include "loudia_peak_detection.hpp"
#include "loudia_spectralodf.hpp"
#include "loudia_onset.hpp"    
#include "loudia_peak_interpolation.hpp"
#include "loudia_peak.hpp"    
#include "loudia_pitch.hpp"    
#include "loudia_resampler.hpp"    
#include "loudia_spectral.hpp"    
#include "loudia_voice.hpp"
    



