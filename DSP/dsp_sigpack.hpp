#pragma once
#include <armadillo>

#define SP_VERSION_MAJOR 1
#define SP_VERSION_MINOR 2
#define SP_VERSION_PATCH 7

#include "base/base.h"
#include "window/window.h"
#include "filter/filter.h"
#include "resampling/resampling.h"
#include "spectrum/spectrum.h"
#include "timing/timing.h"
#include "gplot/gplot.h"
#include "parser/parser.h"
#ifdef HAVE_FFTW
  #include "fftw/fftw.h"
#endif
#include "image/image.h"
#include "kalman/kalman.h"


namespace sp
{
    ///
    /// @defgroup math Math
    /// \brief Math functions.
    /// @{

    const double PI   = 3.14159265358979323846;     ///< ... _or use arma::datum::pi_
    const double PI_2 = 6.28318530717958647692;

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A sinc, sin(x)/x, function.
    /// @param x The angle in radians
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline double sinc( double x )
    {
        if(x==0.0)
            return 1.0;
        else
            return std::sin(PI*x)/(PI*x);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief A sinc, sin(x)/x, function.
    /// @param x The angle in radians
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec sinc(const arma::vec& x)
    {
        arma::vec out;
        out.copy_size(x);
        for (unsigned int n = 0; n < out.size(); n++)
        {
            out(n) = sinc(x(n));
        }
        return out;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Modified first kind bessel function order zero.
    ///
    /// See bessel functions on [Wikipedia](https://en.wikipedia.org/wiki/Bessel_function)
    /// @param x
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline double besseli0( double x )
    {
        double y=1.0,s=1.0,x2=x*x,n=1.0;
        while (s > y*1.0e-9)
        {
            s *= x2/4.0/(n*n);
            y += s;
            n += 1;
        }
        return y;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Calculates angle in radians for complex input.
    /// @param x Complex input value
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    double angle( const std::complex<T>& x )
    {
        return std::arg(x);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Calculates angle in radians for complex input.
    /// @param x Complex input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec angle( const arma::cx_vec& x )
    {
        arma::vec P;
        P.copy_size(x);
        for(unsigned int r=0; r<x.n_rows; r++)
            P(r) = std::arg(x(r));
        return P;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Calculates angle in radians for complex input.
    /// @param x Complex input matrix
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::mat angle( const arma::cx_mat& x )
    {
        arma::mat P;
        P.copy_size(x);
        for(unsigned int r=0; r<x.n_rows; r++)
            for(unsigned int c=0; c<x.n_cols; c++)
                P(r,c) = std::arg(x(r,c));
        return P;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Unwraps the angle vector x, accumulates phase.
    /// @param x Complex input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec unwrap( const arma::vec& x )
    {
        arma::vec P;
        double pacc = 0, pdiff = 0;
        const double thr=PI*170/180;
        P.copy_size(x);
        P(0)=x(0);
        for(unsigned int r=1; r<x.n_rows; r++)
        {
            pdiff = x(r)-x(r-1);
            if( pdiff >= thr ) pacc += -PI_2;
            if( pdiff <= -thr) pacc +=  PI_2;
            P(r)=pacc+x(r);
        }
        return P;
    }
    /// @} // END math


    ///
    /// @defgroup data Data
    /// \brief Data generation/manipulation ...
    /// @{

    /// \brief Generates a linear time vector with specified sample rate. Delta time=1/Fs.
    /// @param N  Number of data points
    /// @param Fs Sample rate
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec timevec( const int N, const double Fs )
    {
        return arma::regspace(0,N-1.0)/Fs;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 1D FFT shift.
    /// @returns Circular shifted FFT
    /// @param Pxx Complex FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Col<T> fftshift(const arma::Col<T>& Pxx)
    {
        arma::Col<T> x(Pxx.n_elem);
        x = shift(Pxx, floor(Pxx.n_elem / 2));
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 1D FFT inverse/reverse shift.
    /// @returns Circular shifted FFT
    /// @param Pxx Complex FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Col<T> ifftshift(const arma::Col<T>& Pxx)
    {
        arma::Col<T> x(Pxx.n_elem);
        x = shift(Pxx, -ceil(Pxx.n_elem / 2));
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 2D FFT shift.
    /// @returns Circular shifted FFT
    /// @param Pxx FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Mat<T> fftshift(const arma::Mat<T>& Pxx)
    {
        arma::uword R = Pxx.n_rows;
        arma::uword C = Pxx.n_cols;
        arma::Mat<T> x(R, C);
        x = arma::shift(Pxx, static_cast<arma::sword>(floor(R / 2)), 0);
        x = arma::shift(x, static_cast<arma::sword>(floor(C / 2)), 1);
        return x;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief 2D FFT inverse/reverse shift.
    /// @returns Circular shifted FFT
    /// @param Pxx FFT
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <typename T>
    arma::Mat<T> ifftshift(const arma::Mat<T>& Pxx)
    {
        arma::uword R = Pxx.n_rows;
        arma::uword C = Pxx.n_cols;
        arma::Mat<T> x(R, C);
        x = shift(Pxx, -ceil(R / 2), 0);
        x = shift(x, -ceil(C / 2), 1);
        return x;
    }

    /// @} // END data

    ///
    /// @defgroup misc Misc
    /// \brief Misc functions, error handling etc.
    /// @{

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief SigPack version string
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline std::string sp_version(void)
    {
        return std::to_string(SP_VERSION_MAJOR)+"."+std::to_string(SP_VERSION_MINOR)+"."+std::to_string(SP_VERSION_PATCH);
    }

    ///////////////////////////////////
    // err_handler("Error string")
    //      Prints an error message, waits for input and
    //      then exits with error
#define err_handler(msg) \
    { \
        std::cout << "SigPack Error [" << __FILE__  << "@" << __LINE__ << "]: " << msg << std::endl; \
        std::cin.get(); \
        exit(EXIT_FAILURE);\
    }

    ///////////////////////////////////
    // wrn_handler("Warning string")
    //      Prints an warning message
#define wrn_handler(msg)  \
    { \
        std::cout << "SigPack warning [" << __FILE__ << "@" << __LINE__ << "]: " << msg << std::endl;\
    }
    /// @} // END misc

    ///
    /// @defgroup window Window
    /// \brief Window functions.
    ///
    /// See window functions at [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    /// @{

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Generic fifth order symmetric cos window.
    ///
    /// \f$ w_i = a_0-a_1\ cos(2\pi i /(N-1))+a_2\ cos(4\pi i /(N-1))-a_3\ cos(6\pi i /(N-1))+a_4\ cos(8\pi i /(N-1))\f$
    /// @returns The cosinus window based on the <b>a</b> vector
    /// @param N Number of window taps
    /// @param a A vector of cosinus coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec cos_win( const arma::uword N, const arma::vec& a )
    {
        arma::vec h(N);
        for(arma::uword i=0; i<N; i++)
        {
            h[i] = a[0] - a[1]*std::cos(1.0*PI_2*i/(N-1)) + a[2]*std::cos(2.0*PI_2*i/(N-1)) \
                   - a[3]*std::cos(3.0*PI_2*i/(N-1)) + a[4]*std::cos(4.0*PI_2*i/(N-1));
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Hamming window.
    ///
    /// \f$ w_i = 0.54-0.46\ cos(2\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec hamming( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.54;
        a[1] = 0.46;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Hann window.
    ///
    /// \f$ w_i = 0.5-0.5\ cos(2\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec hann( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.5;
        a[1] = 0.5;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Blackman window.
    ///
    /// \f$ w_i = 0.42-0.5\ cos(2\pi i /(N-1))+0.08\ cos(4\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec blackman( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.42; // 7938/18608.0
        a[1] = 0.5;  // 9240/18608.0
        a[2] = 0.08; // 1430/18608.0
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Blackman-Harris window.
    /// Symmetric BH4 window
    ///
    /// \f$ w_i = 0.359-0.488\ cos(2\pi i /(N-1))+0.141\ cos(4\pi i /(N-1))-0.011\ cos(6\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec blackmanharris( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.35875;
        a[1] = 0.48829;
        a[2] = 0.14128;
        a[3] = 0.01168;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Flattop window.
    ///
    /// \f$ w_i = 0.216-0.417\ cos(2\pi i /(N-1))+0.277\ cos(4\pi i /(N-1))-0.084\ cos(6\pi i /(N-1))+0.007\ cos(8\pi i /(N-1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec flattopwin( const arma::uword N )
    {
        arma::vec a=arma::zeros<arma::vec>(5);
        a[0] = 0.21557895;
        a[1] = 0.41663158;
        a[2] = 0.277263158;
        a[3] = 0.083578947;
        a[4] = 0.006947368;
        return cos_win(N,a);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Hanning window.
    ///
    /// \f$ w_i = 0.5-0.5\ cos(2\pi (i+1) /(N+1))\f$
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec hanning( const arma::uword N )
    {
        arma::vec h(N);
        for(arma::uword i=0; i<N; i++)
        {
            h[i] = 0.5-0.5*std::cos(PI_2*(i+1)/(N+1));
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Kaiser window.
    ///
    /// See Kaiser window at [Wikipedia](https://en.wikipedia.org/wiki/Window_function#Kaiser_window)
    /// @param N Nr of taps
    /// @param beta Beta factor
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec kaiser( const arma::uword N, double beta )
    {
        arma::vec h(N);
        double bb = besseli0(beta);
        for( arma::uword i=0; i<N; i++)
        {
            h[i] = besseli0(beta*sqrt(4.0*i*(N-1-i))/(N-1))/bb;
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Triangle window.
    ///
    /// See Triangle window at [Wikipedia](https://en.wikipedia.org/wiki/Window_function#Triangular_window)
    /// @param N Nr of taps
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec triang( const arma::uword N )
    {
        arma::vec h(N);
        if(N%2)    // Odd
        {
            for(arma::uword i=0; i<(N-1)/2; i++)
            {
                h[i]     = 2.0*(i+1)/(N+1);
                h[N-i-1] = h[i];
            }
            h[(N-1)/2] = 1.0;
        }
        else      // Even
        {
            for(arma::uword i=0; i<N/2; i++)
            {
                h[i]     = (2.0*i+1)/N;
                h[N-i-1] = h[i];
            }
        }
        return h;
    }
    /// @}

   ///
    /// @defgroup filter Filter
    /// \brief FIR/MA and IIR/ARMA filter functions.
    /// @{

    ///
    /// \brief FIR/MA filter class.
    ///
    /// Implements FIR/MA filter functions as \f[  y(n) = \sum_{k=0}^{M-1}{b_kx(n-k)}=b_0x(n)+b_1x(n-1)+...+b_{M-1}x(n-(M-1))\f]
    /// where M is the number of taps in the FIR filter. The filter order is M-1.
    /// Adaptive update of filter is possible with LMS or NLMS algorithms
    template <class T1, class T2, class T3>
    class FIR_filt
    {
    private:
        // Ordinary FIR filter
        arma::uword M;           ///< Nr of filter taps
        arma::uword cur_p;       ///< Pointer to current sample in buffer
        arma::Mat<T1> buf;       ///< Signal buffer
        arma::Mat<T2> b;         ///< Filter coefficients
        // Adaptive LMS FIR filter
        double mu;               ///< Adaptive filter step size
        arma::uword L;           ///< Adaptive filter block length
        arma::uword blk_ctr;     ///< Adaptive filter block length counter
        T2 c;                    ///< Adaptive filter NLMS regulation const.
        arma::Mat<T1> P;         ///< Adaptive filter Inverse corr matrix (estimated accuracy)
        arma::Mat<T1> Q;         ///< Adaptive filter Process noise
        arma::Mat<T1> R;         ///< Adaptive filter Measurement noise
        arma::Mat<T1> K;         ///< Adaptive filter gain vector
        double lmd;              ///< Adaptive filter RLS forgetting factor
        arma::Mat<T1> X_tpz;     ///< Adaptive filter Toeplitz for Corr matrix calc.
        arma::uword do_adapt;    ///< Adaptive filter enable flag
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        FIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        ~FIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Clears the internal states and pointer.
        ////////////////////////////////////////////////////////////////////////////////////////////
        void clear(void)
        {
            buf.zeros();
            cur_p = 0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in FIR filter.
        /// The internal state and pointers are cleared
        /// @param _b Filter coefficients \f$ [b_0 ..b_{M-1}]^T \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(const arma::Mat<T2> &_b)
        {
            M = _b.n_elem;
            buf.set_size(M,1);
            this->clear();
            b.set_size(M,1);
            b = _b;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in FIR filter (col format)
        /// The internal state and pointers are cleared
        /// @param _b_col Filter coefficients \f$ [b_0 ..b_{M-1}]^T \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(const arma::Col<T2> &_b_col)
        {
          arma::Mat<T2> b_mat = arma::conv_to<arma::Mat<T2> >::from(_b_col);
          set_coeffs(b_mat);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Get coefficients from FIR filter.
        /// @return b Filter coefficients \f$ [b_0 ..b_{M-1}]^T \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Col<T2> get_coeffs()
        {
           return arma::conv_to<arma::Col<T2> >::from(b);
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Updates coefficients in FIR filter without clearing the internal states.
        /// @param _b Filter coefficients \f$ [b_0 ..b_{M-1}] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void update_coeffs(const arma::Mat<T2> &_b)
        {
            b = _b;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter operator.
        /// @return Filtered output
        /// @param in Input sample
        ////////////////////////////////////////////////////////////////////////////////////////////
        T3 operator()(const T1 & in)
        {
            T3 out=0;
            arma::uword p = 0;
            buf[cur_p] = in;                    // Insert new sample
            for( arma::uword m = cur_p; m < M; m++)
                out += b[p++]*buf[m];           // Calc upper part
            for( arma::uword m = 0; m < cur_p; m++)
                out += b[p++]*buf[m];           // ... and lower

            // Move insertion point
            if(cur_p == 0)
                cur_p = M-1;
            else
                cur_p--;

            return out;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter function.
        /// @return Filtered output
        /// @param in Input vector
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Mat<T3> filter(const arma::Mat<T1> & in)
        {
            arma::uword sz = in.n_elem;
            arma::Mat<T3> out(sz,1);
            for( arma::uword n=0;n<sz;n++)
                out[n] = this->operator()(in[n]);
            return out;
        }
        arma::Col<T3> filter(const arma::Col<T1> & in)
        {
           arma::Mat<T1> in_col = arma::conv_to<arma::Mat<T1> >::from(in);
           return arma::conv_to<arma::Col<T3> >::from(filter(in_col));
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief LMS Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _mu Step size
        /// @param _L Block length
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_lms(const arma::uword _N, const double _mu, const  arma::uword _L=1)
        {
            M  = _N;
            mu = _mu;
            L  = _L;
            buf.set_size(M,1);buf.zeros();
            b.set_size(M,1);b.zeros();
            K.set_size(M,1);K.zeros();
            cur_p = 0;
            blk_ctr = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief LMS Filter update function
        ///
        ///  The LMS filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+2\mu\mathbf{x(n)}err(n) \f$ <br>
        ///  where <br> \f$ err(n) = d(n)-\mathbf{b(n-1)^Tx(n)} \f$
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void lms_adapt(const T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(M,1);
                for(arma::uword m=0; m<M; m++)
                {
                    buf_tmp(m) = buf((cur_p+m+1)%M);
                }

                // Accumulate
                K += _err*buf_tmp;

                // Block update
                if(blk_ctr++%L==0)
                {
                      b+=2*mu*K/double(L);
                      K.zeros();
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief NLMS Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _mu Step size
        /// @param _c Regularization factor
        /// @param _L Block length
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_nlms(const  arma::uword _N, const double _mu, const T2 _c, const  arma::uword _L=1)
        {
            M  = _N;
            mu = _mu;
            L  = _L;
            c  = _c;
            buf.set_size(M,1);buf.zeros();
            b.set_size(M,1);b.zeros();
            K.set_size(M,1);K.zeros();
            cur_p = 0;
            blk_ctr = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief NLMS Filter update function
        ///
        ///  The NLMS filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+2\mu\frac{\mathbf{x(n)}err(n)}{c+\mathbf{x(n)^Tx(n)}} \f$ <br>
        ///  where <br> \f$ err(n) = d(n)-\mathbf{b(n-1)^Tx(n)} \f$
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void nlms_adapt(const T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(M,1);
                for(arma::uword m=0; m<M; m++)
                {
                    buf_tmp(m) = buf((cur_p+m+1)%M);
                }

                // Accumulate
                T1 S = c + arma::as_scalar(buf_tmp.t()*buf_tmp);
                K += _err*buf_tmp/S;

                // Block update
                if(blk_ctr++%L==0)
                {
                      b+=2*mu*K/double(L);
                      K.zeros();
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief LMS-Newton Filter function setup. (Affine Projection Algorithm)
        /// @param _N  Number of filter taps
        /// @param _mu Step size
        /// @param _c Regularization factor
        /// @param _L Block length
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_newt(const  arma::uword _N, const double _mu, const T2 _c, const  arma::uword _L=1)
        {
            M  = _N;
            mu = _mu;
            L  = _L;
            c  = _c;
            buf.set_size(M,1);buf.zeros();
            b.set_size(M,1);b.zeros();
            K.set_size(M,1);K.zeros();
            X_tpz.set_size(M,L);X_tpz.zeros();
            cur_p = 0;
            blk_ctr = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief LMS-Newton Filter update function
        ///
        ///  The LMS-Newton filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+2\mu\frac{\mathbf{x(n)}err(n)}{c+\mathbf{R_{xx}}} \f$ <br>
        ///  where <br> \f$ err(n) = d(n)-\mathbf{b(n-1)^Tx(n)} \f$ <br>
        ///  and \f$ \mathbf{R_{xx}} \f$ is the correlation matrix
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void newt_adapt(const T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(M,1);
                for(arma::uword m=0; m<M; m++)
                {
                    buf_tmp(m) = buf((cur_p+m+1)%M);
                }

                // Accumulate in buf
                K += _err*buf_tmp;
                X_tpz.col(blk_ctr%L) = buf_tmp;

                // Block update
                if(blk_ctr++%L==0)
                {
                      // Correlation matrix estimate
                      arma::Mat<T1> Rxx = X_tpz*X_tpz.t()/L;
                      arma::Mat<T1> I; I.eye(M,M);
                      b+=mu*pinv(Rxx+c*I)*K/L;
                      K.zeros();
                }
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief RLS Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _lmd Lambda
        /// @param _P0 Inverse corr matrix initializer
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_rls(const arma::uword _N, const double _lmd,const double _P0)
        {
            M  = _N;
            lmd  = _lmd;
            L = 1;
            P.eye(M,M);
            P =_P0*P;
            K.set_size(M,1);K.zeros();
            buf.set_size(M,1);buf.zeros();
            b.set_size(M,1);b.zeros();
            cur_p = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief RLS Filter update function
        ///
        ///  The RLS filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+\mathbf{Kx(n)}\f$ <br>
        ///  where <br>\f$ \mathbf{K} =\frac{\mathbf{Px}}{\lambda+\mathbf{x^TPx}} \f$ <br>
        ///  and <br>\f$ \mathbf{P^+} =\frac{\mathbf{P^-+xP^-x^T }}{\lambda} \f$ <br>
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void rls_adapt(const T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(M,1);
                for(arma::uword m=0; m<M; m++)
                {
                    buf_tmp(m) = buf((cur_p+m+1)%M);
                }

                // Update P
                T1 S = lmd + arma::as_scalar(buf_tmp.t()*P*buf_tmp);
                K = P*buf_tmp/S;
                P = (P-K*buf_tmp.t()*P)/lmd;

                // Update coeffs
                b = b + K*_err;
            }
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Kalman Filter function setup.
        /// @param _N  Number of filter taps
        /// @param _P0 Inverse corr matrix initializer
        /// @param _Q0 Process noise matrix initializer
        /// @param _R0 Measurement noise matrix initializer
        ////////////////////////////////////////////////////////////////////////////////////////////
        void setup_kalman(const arma::uword _N, const double _P0, const double _Q0, const double _R0)
        {
            M  = _N;
            L = 1;
            P.eye(M,M);
            P =_P0*P;
            Q.eye(M,M);
            Q =_Q0*Q;
            R.ones(1,1);
            R =_R0*R;
            K.set_size(M,1);K.zeros();
            buf.set_size(M,1);buf.zeros();
            b.set_size(M,1);b.zeros();
            cur_p = 0;
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Kalman Filter update function
        ///
        ///  The Kalman filter is updated as <br>
        ///  \f$ \mathbf{b(n)} = \mathbf{b(n-1)}+\mathbf{Kx(n)}\f$ <br>
        ///  where <br>\f$ \mathbf{K} =\frac{\mathbf{Px}}{R+\mathbf{x^TPx}} \f$ <br>
        ///  and <br>\f$ \mathbf{P^+} =\mathbf{P^-+xP^-x^T }+Q \f$ <br>
        /// @param _err  Feedback error
        ////////////////////////////////////////////////////////////////////////////////////////////
        void kalman_adapt(const T3 _err)
        {
            if(do_adapt)
            {
                // Reshape buf
                arma::Mat<T1> buf_tmp(M,1);
                for(arma::uword m=0; m<M; m++)
                {
                    buf_tmp(m) = buf((cur_p+m+1)%M);
                }

                // Innovation/error covariance
                T1 S = arma::as_scalar(R+buf_tmp.t()*P*buf_tmp);

                // Kalman gain
                K = P*buf_tmp/S;

                // Update coeffs/state
                b = b + K*_err;

                // Update estimate covariance
                P = P-K*buf_tmp.t()*P+Q;
            }
        }


        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Get step size
        /// @return Step size mu
        ////////////////////////////////////////////////////////////////////////////////////////////
        double get_step_size(void)
        {
            return mu;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Get P
        /// @return P
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Mat<T1> get_P(void)
        {
            return P;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Get K
        /// @return K
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Mat<T1> get_K(void)
        {
            return K;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Set step size
        /// @param _mu Step size mu
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_step_size(const double _mu)
        {
            mu = _mu;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Start adapt
        ////////////////////////////////////////////////////////////////////////////////////////////
        void adapt_enable(void)
        {
            do_adapt = 1;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Stop adapt
        ////////////////////////////////////////////////////////////////////////////////////////////
        void adapt_disble(void)
        {
            do_adapt = 0;
        }

    };


    ///
    /// \brief IIR/ARMA filter class.
    ///
    /// Implements IIR/ARMA filter functions as \f[  a_0y(n) = b_0x(n)+b_1x(n-1)+...+b_{M-1}x(n-(M-1))-a_1y(n-1)-...-a_{M-1}y(n-(M-1))\f]
    /// where M is the number of taps in the FIR filter part and M is the number of taps in the IIR filter. The filter order is (M-1,M-1)
    ///
    template <class T1, class T2, class T3>
    class IIR_filt
    {
    private:
        arma::uword M;                ///< Nr of MA filter taps
        arma::uword N;                ///< Nr of AR filter taps
        arma::uword b_cur_p;          ///< Pointer to current sample in MA buffer
        arma::uword a_cur_p;          ///< Pointer to current sample in AR buffer
        arma::Col<T2> b;      ///< MA Filter coefficients
        arma::Col<T2> a;      ///< AR Filter coefficients
        arma::Col<T1> b_buf;  ///< MA Signal buffer
        arma::Col<T1> a_buf;  ///< AR Signal buffer
    public:
        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Constructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        IIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Destructor.
        ////////////////////////////////////////////////////////////////////////////////////////////
        ~IIR_filt(){}

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Clears the internal states and pointers.
        ////////////////////////////////////////////////////////////////////////////////////////////
        void clear(void)
        {
            b_buf.zeros();
            a_buf.zeros();
            b_cur_p = 0;
            a_cur_p = 0;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Sets coefficients in IIR filter.
        /// The internal state and pointers are cleared
        /// @param _b Filter coefficients \f$ [b_0 ..b_M] \f$
        /// @param _a Filter coefficients \f$ [a_0 ..a_N] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void set_coeffs(const arma::Col<T2> &_b,const arma::Col<T2> &_a)
        {
            M = _b.size();
            N = _a.size();
            b_buf.set_size(M);
            a_buf.set_size(N);
            this->clear();
            b = _b/_a[0];
            a = _a/_a[0];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Updates coefficients in filter without clearing the internal states.
        /// @param _b Filter coefficients \f$ [b_0 ..b_M] \f$
        /// @param _a Filter coefficients \f$ [a_0 ..a_N] \f$
        ////////////////////////////////////////////////////////////////////////////////////////////
        void update_coeffs(const arma::Col<T2> &_b,const arma::Col<T2> &_a)
        {
            b = _b/_a[0];
            a = _a/_a[0];
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter operator.
        /// @return Filtered output
        /// @param in Input sample
        ////////////////////////////////////////////////////////////////////////////////////////////
        T3 operator()(const T1 & in)
        {
            T3 out=0;
            arma::uword p = 0;

            // MA part
            b_buf[b_cur_p] = in;                // Insert new sample
            for(arma::uword m = b_cur_p; m < M; m++)
                out += b[p++]*b_buf[m];         // Calc upper part
            for(arma::uword m = 0; m < b_cur_p; m++)
                out += b[p++]*b_buf[m];         // ... and lower

            // Move insertion point
            if(b_cur_p == 0)
                b_cur_p = M-1;
            else
                b_cur_p--;

            // AR part
            p=1;
            for(arma::uword n = a_cur_p+1; n < N; n++)
                out -= a[p++]*a_buf[n];         // Calc upper part
            for(arma::uword n = 0; n < a_cur_p; n++)
                out -= a[p++]*a_buf[n];         // ... and lower

            a_buf[a_cur_p] = out;		        // Insert output

            // Move insertion point
            if(a_cur_p == 0)
                a_cur_p = N-1;
            else
                a_cur_p--;

            return out;
        }

        ////////////////////////////////////////////////////////////////////////////////////////////
        /// \brief Filter function.
        /// @return Filtered output
        /// @param in Input vector
        ////////////////////////////////////////////////////////////////////////////////////////////
        arma::Col<T3> filter(const arma::Col<T1> & in)
        {
            arma::uword sz = in.size();
            arma::Col<T3> out(sz);
            for( arma::uword n=0;n<sz;n++)
                out[n] = this->operator()(in[n]);
            return out;
        }
    };


    ///
    /// Filter design functions
    ///

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief FIR lowpass design function.
    /// FIR lowpassdesign using windows method (hamming window).
    /// NB! Returns size M+1
    /// @return b Filter coefficients \f$ [b_0 ..b_N] \f$
    /// @param M Filter order
    /// @param f0 Filter cutoff frequency in interval [0..1]
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec fir1(const arma::uword M, const double f0)
    {
        arma::vec b(M+1), h(M+1);
        h = hamming(M+1);
        for (arma::uword m=0;m<M+1;m++)
        {
            b[m] = f0*h[m]*sinc(f0*(m-M/2.0));
        }
        return b/arma::sum(b);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief FIR highpass design function.
    /// FIR design using windows method (hamming window).
    /// NB! Returns size M+1
    /// @return b Filter coefficients \f$ [b_0 ..b_N] \f$
    /// @param M Filter order (must be even)
    /// @param f0 Filter cutoff frequency in interval [0..1]
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec fir1_hp(const arma::uword M, const double f0)
    {
        if(M%2 != 0)
            err_handler("Filter order must be even");
        arma::vec b(M+1), h(M+1);
        h = hamming(M+1);
        for (arma::uword m=0;m<M+1;m++)
        {
            b[m] = h[m]*(sinc(m-M/2.0)-f0*sinc(f0*(m-M/2.0)));
        }

        // Scale
        std::complex<double> i(0,1);
        double nrm;
        arma::vec fv=arma::regspace(0,double(M));
        nrm = abs(arma::sum(exp(-i*fv*PI)%b));

        return b/nrm;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief FIR bandpass design function.
    /// FIR design using windows method (hamming window).
    /// NB! Returns size M+1
    /// @return b Filter coefficients \f$ [b_0 ..b_N] \f$
    /// @param M Filter order
    /// @param f0 Filter low cutoff frequency in interval [0..1]
    /// @param f1 Filter high cutoff frequency in interval [0..1]
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec fir1_bp(const arma::uword M, const double f0, const double f1)
    {
        if(f1<=f0)
            err_handler("Frequencies must be [0 < f0 < f1 < 1]");

        arma::vec b(M+1), h(M+1);
        h = hamming(M+1);
        for (arma::uword m=0;m<M+1;m++)
        {
            b[m] = h[m]*(f1*sinc(f1*(m-M/2.0))-f0*sinc(f0*(m-M/2.0)));
        }

        // Scale
        double fc = (f0+f1)/2;
        std::complex<double> i(0,1);
        double nrm;
        arma::vec fv=arma::regspace(0,double(M));
        nrm = abs(arma::sum(exp(-i*fv*PI*fc)%b));

        return b/nrm;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief FIR bandstop design function.
    /// FIR design using windows method (hamming window).
    /// NB! Returns size M+1
    /// @return b Filter coefficients \f$ [b_0 ..b_N] \f$
    /// @param M Filter order (must be even)
    /// @param f0 Filter low cutoff frequency in interval [0..1]
    /// @param f1 Filter high cutoff frequency in interval [0..1]
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec fir1_bs(const arma::uword M, const double f0, const double f1)
    {
        if(M%2 != 0)
            err_handler("Filter order must be even");
        if(f1<=f0)
            err_handler("Frequencies must be [0 < f0 < f1 < 1]");

        arma::vec b(M+1), h(M+1);
        h = hamming(M+1);
        for (arma::uword m=0;m<M+1;m++)
        {
            b[m] = h[m]*(sinc(m-M/2.0)-f1*sinc(f1*(m-M/2.0))+f0*sinc(f0*(m-M/2.0)));
        }

        return b/arma::sum(b);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Fractional delay function.
    /// Fractional delay filter design using windowed sinc method.
    /// Actual delay is M/2+fd samples for even nr of taps and
    /// (M-1)/2+fd for odd nr of taps
    /// Best performance if -1 < fd < 1
    /// @param M Filter length
    /// @param fd Fractional delay
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec fd_filter( const arma::uword M, double fd )
    {
        arma::vec h(M);
        arma::vec w = blackmanharris(M);
        if( M % 2 == 1 ) fd = fd-0.5; // Offset for odd nr of taps
        for(arma::uword m=0;m<M;m++)
        {
            h(m) = w(m)*sinc(m-M/2.0-fd);
        }
        h = h/arma::sum(h);  // Normalize gain

        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency response function.
    /// Calculates the frequency response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param K Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::cx_vec freq( const arma::vec b, const arma::vec a, const arma::uword K=512)
    {
        arma::cx_vec h(K);
        arma::uword M = b.size();
        arma::uword N = a.size();
        std::complex<double> b_tmp,a_tmp,i(0,1);
        for(arma::uword k=0;k<K;k++)
        {
            b_tmp=std::complex<double>(b(0),0);
            for(arma::uword m=1;m<M;m++)
                b_tmp+= b(m)*(cos(m*PI*k/K)-i*sin(m*PI*k/K));
            a_tmp=std::complex<double>(a(0),0);
            for(arma::uword n=1;n<N;n++)
                a_tmp+= a(n)*(cos(n*PI*k/K)-i*sin(n*PI*k/K));
            h(k) = b_tmp/a_tmp;
        }
        return h;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency magnitude response function.
    /// Calculates the frequency magnitude response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param K Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec freqz( const arma::vec b, const arma::vec a, const arma::uword K=512)
    {
        arma::cx_vec f = freq(b,a,K);
        return abs(f);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Frequency phase response function.
    /// Calculates the frequency phase response
    /// @param b FIR/MA filter coefficients
    /// @param a IIR/AR filter coefficients
    /// @param K Number of evaluation points, Default 512
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::vec phasez( const arma::vec b, const arma::vec a, const arma::uword K=512)
    {
        arma::cx_vec f = freq(b,a,K);
        return angle(f);
    }
    /// @}

///
/// @defgroup resampling Resampling
/// \brief Resampling functions.
/// @{

////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Upsampling without anti alias filtering.
/// @returns A vector with p-1 zeros inserted in the input vector [x0,0,0,..,x1,0,0,..-..,xN,0,0,..]
/// @param x Input vector
/// @param p Upsampling factor
////////////////////////////////////////////////////////////////////////////////////////////
template <class T1>
arma::Col<T1> upsample(const arma::Col<T1>& x, const int p )
{
    long int N = x.size();
    arma::Col<T1> y;
    y.set_size(p*N);
    y.zeros();
    for(long int n=0; n<N; n++)
        y[p*n] = x[n];
    return y;
}

////////////////////////////////////////////////////////////////////////////////////////////
/// \brief Downsampling without anti alias filtering.
/// @returns A vector with every q:th value from the input vector
/// @param x Input vector
/// @param q Downsampling factor
////////////////////////////////////////////////////////////////////////////////////////////
template <class T1>
arma::Col<T1> downsample(const arma::Col<T1>& x, const int q )
{
    arma::Col<T1> y;
    int N = int(floor(1.0*x.size()/q));
    y.set_size(N);
    for(long int n=0; n<N; n++)
        y[n] = x[n*q];
    return y;
}

///
/// \brief A resampling class.
///
/// Implements up/downsampling functions
///
template <class T1>
class resampling
{
private:
    FIR_filt<T1,double,T1> aa_filt;
    arma::vec H;         ///< Filter coefficients
    arma::vec K;         ///< Number of filter coefficients
    arma::uword P;       ///< Upsampling rate
    arma::uword Q;       ///< Downsampling rate

public:
    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Constructor.
    /// @param _P Upsampling rate
    /// @param _Q Downsampling rate
    /// @param _H FIR filter coefficients
    ////////////////////////////////////////////////////////////////////////////////////////////
    resampling(const arma::uword _P,const arma::uword _Q,const arma::vec _H)
    {
        P = _P;
        Q = _Q;
        H = _H;
        K = H.n_elem;
        aa_filt.set_coeffs(H);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Constructor using a #fir1 filter with 8*M+1 taps and cutoff 1/M where M=max(P,Q)
    /// @param _P Upsampling rate
    /// @param _Q Downsampling rate
    ////////////////////////////////////////////////////////////////////////////////////////////
    resampling(const arma::uword _P,const arma::uword _Q)
    {
        P = _P;
        Q = _Q;
        arma::uword M=(P>Q)?P:Q;
        H = fir1(8*M,1.0f/M);
        K = H.n_elem;
        aa_filt.set_coeffs(H);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Destructor.
    ////////////////////////////////////////////////////////////////////////////////////////////
    ~resampling() {}


    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Downsampling with anti alias filter.
    ///
    /// @param in  Input vector
    /// @param out Output vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    void downfir(const arma::Col<T1>& in, arma::Col<T1>& out)
    {
        arma::uword sz = in.n_elem;
        for( arma::uword n=0; n<sz; n++)
        {
            T1 tmp = aa_filt(in[n]);
            if(n%Q==0)
                out[n/Q] = tmp;
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Upsampling with anti alias filter.
    ///
    /// @param in  Input vector
    /// @param out Output vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    void upfir(const arma::Col<T1>& in, arma::Col<T1>& out)
    {
        arma::uword sz = P*in.n_elem;
        for( arma::uword n=0; n<sz; n++)
        {
            if(n%P==0)
                out[n] = P*aa_filt(in[n/P]);
            else
                out[n] = P*aa_filt(0.0);
        }
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Resampling by a rational P/Q with anti alias filtering.
    ///
    /// The caller needs to allocate the input and output vector so that length(out)==length(in)*P/Q
    /// @param in  Input vector
    /// @param out Output vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    void upfirdown(const arma::Col<T1>& in, arma::Col<T1>& out)
    {
        arma::uword sz = P*in.n_elem;
        T1 tmp;
        for( arma::uword n=0; n<sz; n++)
        {
            if(n%P==0)
                tmp = aa_filt(in[n/P]);
            else
                tmp = aa_filt(0.0);
            if(n%Q==0)
                out[n/Q] = P*tmp;
        }
    }

      ///
    /// @defgroup spectrum Spectrum
    /// \brief Spectrum functions.
    /// @{

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Windowed spectrum calculation.
    ///
    /// The spectrum is calculated using the fast fourier transform of the windowed input data vector
    /// @returns A complex spectrum vector
    /// @param x Input vector
    /// @param W Window function vector. NB! Must be same size as input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::cx_vec spectrum(const arma::Col<T1>& x, const arma::vec& W)
    {
        arma::cx_vec Pxx(x.size());
        double wc = sum(W);     // Window correction factor
        Pxx = fft(x % W)/wc;    // FFT calc
        return Pxx;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrum density calculation using windowed data.
    /// @returns A real valued PSD vector
    /// @param x Input vector
    /// @param W Window function vector. NB! Must be same size as input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec psd(const arma::Col<T1>& x, const arma::vec& W)
    {
        arma::cx_vec X(x.size());
        arma::vec Pxx(x.size());
        X = spectrum(x,W);          // FFT calc
        Pxx = real(X % conj(X));    // Calc power spectra
        return Pxx;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrum density calculation using Hamming windowed data.
    /// @returns A real valued PSD vector
    /// @param x Input vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec psd(const arma::Col<T1>& x)
    {
        arma::vec W;
        W = hamming(x.size());
        return psd(x,W);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Spectrogram calculation using Hamming windowed data.
    ///
    /// See spectrogram at [Wikipedia](https://en.wikipedia.org/wiki/Spectrogram)
    /// @returns A complex spectrogram matrix
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::cx_mat specgram_cx(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::cx_mat Pw;

        //Def params
        arma::uword N = x.size();
        arma::uword D = Nfft-Noverl;
        arma::uword m = 0;
        if(N > Nfft)
        {
            arma::Col<T1> xk(Nfft);
            arma::vec W(Nfft);

            W = hamming(Nfft);
            arma::uword U = static_cast<arma::uword>(floor((N-Noverl)/double(D)));
            Pw.set_size(Nfft,U);
            Pw.zeros();

            // Avg loop
            for(arma::uword k=0; k<N-Nfft; k+=D)
            {
                xk = x.rows(k,k+Nfft-1);       // Pick out chunk
                Pw.col(m++) = spectrum(xk,W);  // Calculate spectrum
            }
        }
        else
        {
            arma::vec W(N);
            W = hamming(N);
            Pw.set_size(N,1);
            Pw = spectrum(x,W);
        }
        return Pw;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrogram calculation.
    ///
    /// See spectrogram at [Wikipedia](https://en.wikipedia.org/wiki/Spectrogram)
    /// @returns A power spectrogram matrix
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::mat specgram(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::cx_mat Pw;
        arma::mat Sg;
        Pw = specgram_cx(x,Nfft,Noverl);
        Sg = real(Pw % conj(Pw));              // Calculate power spectrum
        return Sg;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Phase spectrogram calculation.
    ///
    /// See spectrogram at [Wikipedia](https://en.wikipedia.org/wiki/Spectrogram)
    /// @returns A phase spectrogram matrix
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::mat specgram_ph(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::cx_mat Pw;
        arma::mat Sg;
        Pw = specgram_cx(x,Nfft,Noverl);
        Sg = angle(Pw);                        // Calculate phase spectrum
        return Sg;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Phase spectrum calculation using Welch's method.
    ///
    /// See Welch's method at [Wikipedia](https://en.wikipedia.org/wiki/Welch%27s_method)
    /// @returns A phase spectrum vector
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec pwelch_ph(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::mat Ph;
        Ph  = specgram_ph(x,Nfft,Noverl);
        return arma::mean(Ph,1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Power spectrum calculation using Welch's method.
    ///
    /// _abs(pwelch(x,Nfft,Noverl))_ is equivalent to Matlab's: _pwelch(x,Nfft,Noverl,'twosided','power')_ <br>
    /// See Welch's method at [Wikipedia](https://en.wikipedia.org/wiki/Welch%27s_method)
    /// @returns A power spectrum vector
    /// @param x Input vector
    /// @param Nfft  FFT size
    /// @param Noverl FFT overlap size
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::vec pwelch(const arma::Col<T1>& x, const arma::uword Nfft=512, const arma::uword Noverl=256)
    {
        arma::mat Pxx;
        Pxx = specgram(x,Nfft,Noverl);
        return arma::mean(Pxx,1);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief DFT calculation of a single frequency using Goertzel's method.
    ///
    /// For more details see [Sysel and Rajmic](https://asp-eurasipjournals.springeropen.com/track/pdf/10.1186/1687-6180-2012-56?site=asp.eurasipjournals.springeropen.com)
    /// @returns The DFT of frequency f
    /// @param x Input vector
    /// @param f Frequency index
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    std::complex<double> goertzel(const arma::Col<T1>& x, const double f)
    {
        // Constants
        arma::uword N = x.size();
        double      Q = f/N;
        double      A = PI_2*Q;
        double      B = 2*cos(A);
        std::complex<double> C(cos(A),-sin(A));
        // States
        T1 s0 = 0;
        T1 s1 = 0;
        T1 s2 = 0;

        // Accumulate data
        for (arma::uword n=0;n<N;n++)
        {
            // Update filter
            s0 = x(n)+B*s1-s2;

            // Shift buffer
            s2 = s1;
            s1 = s0;
        }
        // Update output state
        s0 = B*s1-s2;

        // Return the complex DFT output
        return s0-s1*C;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief DFT calculation of a vector of frequencies using Goertzel's method.
    ///
    /// For more details see [Sysel and Rajmic](https://asp-eurasipjournals.springeropen.com/track/pdf/10.1186/1687-6180-2012-56?site=asp.eurasipjournals.springeropen.com)
    /// @returns The DFT of frequency f
    /// @param x Input vector
    /// @param f Frequency index vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    template <class T1>
    arma::cx_vec goertzel(const arma::Col<T1>& x, const arma::vec f)
    {
        arma::uword N = f.size();
        arma::cx_vec P(N);
        for (arma::uword n=0;n<N;n++)
        {
            P(n) = goertzel(x,f(n));
        }

        // Return the complex DFT output vector
        return P;
    }

    /// @}
        ///
    /// @defgroup timing Timing
    /// \brief Timing functions.
    /// @{

    ///
    /// \brief A delay class.
    ///
    /// Implements different timing related functions such as delay
    ///
    template <class T1>
    class Delay
    {
        private:
            arma::uword D;        ///< The delay value
            arma::uword cur_p;    ///< Pointer to current sample in buffer
            arma::Col<T1> buf;    ///< Signal buffer
        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            Delay()
            {
                cur_p = 0;
                D = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor with delay input.
            /// @param _D delay
            ////////////////////////////////////////////////////////////////////////////////////////////
            Delay(const arma::uword _D)
            {
                set_delay(_D);
                clear();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~Delay() {}

            ////////////////////////////////////////////////////////////////////////////////////////////
            ///  \brief Clears internal state.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void clear(void)
            {
                buf.zeros();
                cur_p = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Sets delay.
            /// @param _D delay
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_delay(const arma::uword _D)
            {
                D = _D+1;
                buf.set_size(D);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief A delay operator.
            /// @param in sample input
            ////////////////////////////////////////////////////////////////////////////////////////////
            T1 operator()(const T1& in)
            {
                buf[cur_p] = in;                    // Insert new sample
                // Move insertion point
                if (cur_p == 0)
                    cur_p = D-1;
                else
                    cur_p--;

                return buf[cur_p];
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief A delay operator (vector version).
            /// @param in vector input
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::Col<T1> delay(const arma::Col<T1>& in)
            {
                arma::uword sz = in.size();
                arma::Col<T1> out(sz);
                for(arma::uword n=0; n<sz; n++)
                    out[n] = this->operator()(in[n]);
                return out;
            }
    };
    /// @}
        ///
    /// @defgroup gplot GPlot
    /// \brief Collection of Gnuplot functions
    /// @{

    ///
    /// \brief Gnuplot class.
    ///
    /// Implements a class for streaming data to Gnuplot using a pipe.
    /// Inspiration from https://code.google.com/p/gnuplot-cpp/
    ///
    /// Requires Gnuplot version > 5.0
    /// \note In Windows only one class is allowed. Using multiple figures are controlled by a figure number. In Linux we may use one instance per figure.
    ///
    class gplot
    {
        private:
            FILE*           gnucmd;          ///< File handle to pipe
            std::string     term;
            int             fig_ix;
            int             plot_ix;

            struct plot_data_s
            {
                std::string label;
                std::string linespec;
            };

            std::vector<plot_data_s> plotlist;

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot y vs. x.
            /// @param x x vector
            /// @param y y vector
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T1, typename T2>
            void plot_str2( const T1& x, const T2& y)
            {
                std::ostringstream tmp_s;
                std::string s;
                tmp_s << "$Dxy" << plot_ix << " << EOD \n";
                arma::uword Nelem = x.n_elem;
                for(arma::uword n=0; n<Nelem; n++)
                {
                    tmp_s << x(n) << " " << y(n);
                    s = tmp_s.str();
                    send2gp(s.c_str());
                    tmp_s.str(""); // Clear buffer
                }
                send2gp("EOD");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get type
            /// @param x     x input type
            /// @returns     Type name
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            std::string get_type(T x)
            {
                if (typeid(x) == typeid(int8_t))     return "int8";
                if (typeid(x) == typeid(uint8_t))    return "uint8";
                if (typeid(x) == typeid(int16_t))    return "int16";
                if (typeid(x) == typeid(uint16_t))   return "uint16";
                if (typeid(x) == typeid(int32_t))    return "int32";
                if (typeid(x) == typeid(uint32_t))   return "uint32";
                if (typeid(x) == typeid(int64_t))    return "int64";
                if (typeid(x) == typeid(arma::sword))return "int64";
                if (typeid(x) == typeid(uint64_t))   return "uint64";
                if (typeid(x) == typeid(arma::uword))return "uint64";
                if (typeid(x) == typeid(float))      return "float32";
                if (typeid(x) == typeid(double))     return "float64";
                err_handler("Unknown type");
            }

        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ///
            /// Opens a pipe to gnuplot program. Make sure it is found/accessable by the system.
            ////////////////////////////////////////////////////////////////////////////////////////////
            gplot()
            {
#if defined(_MSC_VER)
                gnucmd = _popen("gnuplot -persist 2> NUL","wb");
                term = "win";
//#elif defined(_APPLE_)
//            gnucmd = popen("gnuplot -persist &> /dev/null","w");
//#define term "aqua"
#else
                gnucmd = popen("gnuplot -persist","w");
                term = "x11";
#endif
                if(!gnucmd)
                {
                    err_handler("Could not start gnuplot");
                }

                // Set global params
                plot_ix   = 0;
                fig_ix    = 0;
                plotlist.clear();

            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~gplot()
            {
#if defined(_MSC_VER)
                _pclose(gnucmd);
#else
                pclose(gnucmd);
#endif
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Send command to Gnuplot pipe.
            /// @param cmdstr  Command string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void send2gp(const char* cmdstr)
            {
                std::string s_in(cmdstr);
                std::string tmp=s_in+"\n";
                std::fputs(tmp.c_str(), gnucmd );
//                std::cout << tmp.c_str() << std::endl;
            }

			////////////////////////////////////////////////////////////////////////////////////////////
			/// \brief Flush command buffer to Gnuplot pipe.
			////////////////////////////////////////////////////////////////////////////////////////////
			void flush_cmd_buf(void)
			{
				std::fflush(gnucmd);
			}

			////////////////////////////////////////////////////////////////////////////////////////////
			/// \brief Updates gnuplot instantly. (Flushes the command buffer)
			////////////////////////////////////////////////////////////////////////////////////////////
			void draw_now(void)
			{
				std::fflush(gnucmd);
			}

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Sets the active figure.
            /// @param fig  Figure number
            ////////////////////////////////////////////////////////////////////////////////////////////
            void figure(const int fig)
            {
                fig_ix = fig;
                std::ostringstream tmp_s;
                tmp_s << "set term " << term << " " << fig;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                send2gp("reset");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Configure the figure used Windows environment.
            /// @param fig     Figure number
            /// @param name    Window name
            /// @param x       x position of upper left corner
            /// @param y       y position of upper left corner
            /// @param width   width of window
            /// @param height  height of window
            ////////////////////////////////////////////////////////////////////////////////////////////
            void window(const int fig, const char* name,const int x,const int y,const int width,const int height)
            {
                fig_ix = fig;
                std::ostringstream tmp_s;
                tmp_s << "set term " << term << " " << fig << " title \"" << name << "\" position " << x << "," << y << " size " << width << "," << height;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                send2gp("reset");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Configure the figure/window - used in Linux environment where no figure numbers are needed.
            /// @param name    Window name
            /// @param x       x position of upper left corner
            /// @param y       y position of upper left corner
            /// @param width   width of window
            /// @param height  height of window
            ////////////////////////////////////////////////////////////////////////////////////////////
            void window(const char* name,const int x,const int y,const int width,const int height)
            {
                window(0,name,x,y,width,height);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Close window
            ////////////////////////////////////////////////////////////////////////////////////////////
            void close_window(void)
            {
                std::ostringstream tmp_s;
                tmp_s << "set term " << term << " close";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                send2gp("reset");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set grid
            ////////////////////////////////////////////////////////////////////////////////////////////
            void grid_on(void)
            {
                send2gp("set grid");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set grid
            ////////////////////////////////////////////////////////////////////////////////////////////
            void grid_off(void)
            {
                send2gp("unset grid");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set label for X-axis.
            /// @param label label string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void xlabel(const char* label)
            {
                std::ostringstream tmp_s;
                tmp_s << "set xlabel \"" << label << "\" ";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set label for X-axis.
            /// @param label label string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ylabel(const char* label)
            {
                std::ostringstream tmp_s;
                tmp_s << "set ylabel \"" << label << "\" ";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set label at position x,y.
            /// @param x x value
            /// @param y y value
            /// @param label label string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void label(const double& x, const double& y, const char* label)
            {
                std::ostringstream tmp_s;
                tmp_s << "set label \"" << label << "\" at " << x << "," << y;
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set windowtitle.
            /// @param name title string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void title(const char* name)
            {
                std::ostringstream tmp_s;
                tmp_s << "set title \"" << name << " \" ";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set X-axis range.
            /// @param xmin xmin
            /// @param xmax xmax
            ////////////////////////////////////////////////////////////////////////////////////////////
            void xlim(const double xmin, const double xmax)
            {
                std::ostringstream tmp_s;
                tmp_s << "set xrange [" << xmin << ":" << xmax << "]";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set Y-axis range.
            /// @param ymin ymin
            /// @param ymax ymax
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ylim(const double ymin, const double ymax)
            {
                std::ostringstream tmp_s;
                tmp_s << "set yrange [" << ymin << ":" << ymax << "]";
                std::string s = tmp_s.str();
                send2gp(s.c_str());
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Push plot y vs. x with label and linespec
            /// @param x      x vector
            /// @param y      y vector
            /// @param lb     label
            /// @param ls     line spec
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T1, typename T2>
            void plot_add( const T1& x, const T2& y, const std::string lb, const std::string ls="lines")
            {
                plot_data_s pd;

                pd.linespec = ls;
                pd.label    = lb;

                plotlist.push_back(pd);
                plot_str2(x,y);
                plot_ix++;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Push plot y vs. x with label and linespec
            /// @param y      y vector
            /// @param lb     label
            /// @param ls     line spec
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T1>
            void plot_add( const T1& y, const std::string lb, const std::string ls="lines")
            {
                arma::vec x=arma::regspace(0,double(y.n_elem-1));
                plot_data_s pd;

                pd.linespec = ls;
                pd.label    = lb;

                plotlist.push_back(pd);
                plot_str2(x,y);
                plot_ix++;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Push multiple plot, each row gives a plot without label
            /// @param y      y matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot_add_mat( const arma::mat& y)
            {
                arma::vec x=arma::regspace(0,double(y.n_cols-1));
                plot_data_s pd;

                pd.linespec = "lines";
                pd.label    = "";
                for(arma::uword r=0;r<y.n_rows;r++)
                {
                    plotlist.push_back(pd);
                    plot_str2(x,y.row(r));
                    plot_ix++;
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Push multiple plot, each row gives a plot with prefix label
            /// @param y      y matrix
            /// @param p_lb   Label prefix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot_add_mat( const arma::mat& y, const std::string p_lb)
            {
                arma::vec x=arma::regspace(0,double(y.n_cols-1));
                plot_data_s pd;
                pd.linespec = "lines";

                for(arma::uword r=0;r<y.n_rows;r++)
                {
                    std::ostringstream tmp_s;
                    tmp_s << p_lb << r;
                    std::string s = tmp_s.str();
                    pd.label = s;
                    plotlist.push_back(pd);
                    plot_str2(x,y.row(r));
                    plot_ix++;
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Show plots
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot_show(void)
            {
                std::ostringstream tmp_s;

                tmp_s << "plot $Dxy0 title \"" << plotlist[0].label << "\" with " << plotlist[0].linespec;
                for(int r=1; r<plot_ix; r++)
                {
                    tmp_s << " ,$Dxy" << r <<" title \"" << plotlist[r].label << "\" with " << plotlist[r].linespec;
                }
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                plotlist.clear();
                plot_ix = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Clear plots
            ////////////////////////////////////////////////////////////////////////////////////////////
            void plot_clear(void)
            {
                plotlist.clear();
                plot_ix = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot vector, fast version
            /// @param x     x vector
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            void fast_plot( const arma::Col<T>& x, const std::string fmt_args = "with lines")
            {
                std::string fmt=get_type(x.at(0));
                std::string s;
                send2gp("unset key");

                s= "plot '-' binary format='%"+
                   fmt + "' array=("+std::to_string(x.n_elem)+") "+fmt_args;
                send2gp(s.c_str());
                std::fwrite(x.memptr(), sizeof(x.at(0)),x.n_elem,gnucmd );
                flush_cmd_buf();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot vector, fast version
			///
			/// x and y needs to have the same type
            /// @param x     x vector
			/// @param y     y vector
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            void fast_plot( const arma::Col<T>& x, const arma::Col<T>& y, const std::string fmt_args = "with lines")
            {
                std::string fmt1=get_type(x.at(0));
                std::string s;
                send2gp("unset key");
                const arma::uword N=x.n_elem;
                arma::Col<T> v(2*N);
                for(arma::uword n=0;n<N;n++)
                {
                    v.at(2*n) = x.at(n);
                    v.at(2*n+1) = y.at(n);
                }

                s= "plot '-' binary format='%"+fmt1+"' record=("+std::to_string(x.n_elem)+ ") "+fmt_args;
                send2gp(s.c_str());
                std::fwrite(v.memptr(), sizeof(x.at(0)),v.n_elem,gnucmd );
                flush_cmd_buf();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot mat as image
            /// @param x     x matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            void image( const arma::Mat<T>& x)
            {
                xlim(-0.5,x.n_cols-0.5);
                ylim(x.n_rows-0.5,-0.5);
                send2gp("unset key");
                std::string fmt=get_type(x.at(0));
                std::string s;

                s= "plot '-' binary array=(" +
                   std::to_string(x.n_cols) + "," +
                   std::to_string(x.n_rows) + ") scan=yx format='%" +
                   fmt + "' w image";
                send2gp(s.c_str());
                std::fwrite(x.memptr(), sizeof(x.at(0)),x.n_elem,gnucmd );
                flush_cmd_buf();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot cube as image
            /// @param x     x matrix (R,G,B)
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            void image( const arma::Cube<T>& x)
            {
                xlim(-0.5,x.n_cols-0.5);
                ylim(x.n_rows-0.5,-0.5);
                send2gp("unset key");
                std::string fmt=get_type(x.at(0));
                std::string s;

                // Conv cube to gnuplot rgb array
                arma::Cube<T> gp_im(arma::size(x));
                T* ptr=gp_im.memptr();
                for(arma::uword c=0;c<x.n_cols; c++ )
                {
                    for(arma::uword r=0;r<x.n_rows;  r++ )
                    {
                        *ptr++ = x.at(r,c,0);    // R
                        *ptr++ = x.at(r,c,1);    // G
                        *ptr++ = x.at(r,c,2);    // B
                    }
                }

                s= "plot '-' binary array=(" +
                   std::to_string(x.n_cols) + "," +
                   std::to_string(x.n_rows) + ") scan=yx format='%" +
                   fmt+ "' w rgbimage";
                send2gp(s.c_str());
                std::fwrite(gp_im.memptr(), sizeof(x.at(0)),gp_im.n_elem,gnucmd );
                flush_cmd_buf();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot mat as mesh
            /// @param x     x matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            void mesh( const arma::Mat<T>& x)
            {
                send2gp("unset key");
                send2gp("set hidden3d");
                std::string fmt=get_type(x.at(0));
                std::string s= "splot '-' binary array=(" +
                   std::to_string(x.n_cols) + "," +
                   std::to_string(x.n_rows) + ") scan=yx format='%" +
                   fmt+ "' w lines";
                send2gp(s.c_str());
                std::fwrite(x.memptr(),sizeof(x.at(0)),x.n_elem,gnucmd );
                flush_cmd_buf();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Plot mat as surf
            /// @param x     x matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            void surf( const arma::Mat<T>& x)
            {
                send2gp("set pm3d");
                mesh(x);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set linetype to Matlab 'parula' NB! doesn't work with X11 -terminal
            /// Data from https://github.com/Gnuplotting/gnuplot-palettes
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_parula_line(void)
            {
                send2gp("set linetype 1 lc rgb '#0072bd' "); // blue
                send2gp("set linetype 2 lc rgb '#d95319' "); // orange
                send2gp("set linetype 3 lc rgb '#edb120' "); // yellow
                send2gp("set linetype 4 lc rgb '#7e2f8e' "); // purple
                send2gp("set linetype 5 lc rgb '#77ac30' "); // green
                send2gp("set linetype 6 lc rgb '#4dbeee' "); // light-blue
                send2gp("set linetype 7 lc rgb '#a2142f' "); // red
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set linetype to Matlab 'jet' NB! doesn't work with X11 -terminal
            /// Data from https://github.com/Gnuplotting/gnuplot-palettes
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_jet_line(void)
            {
                send2gp("set linetype 1 lc rgb '#0000ff' "); // blue
                send2gp("set linetype 2 lc rgb '#007f00' "); // green
                send2gp("set linetype 3 lc rgb '#ff0000' "); // red
                send2gp("set linetype 4 lc rgb '#00bfbf' "); // cyan
                send2gp("set linetype 5 lc rgb '#bf00bf' "); // pink
                send2gp("set linetype 6 lc rgb '#bfbf00' "); // yellow
                send2gp("set linetype 7 lc rgb '#3f3f3f' "); // black
            }
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set linetype to Matlab 'parula' NB! doesn't work with X11 -terminal
            /// Data from https://github.com/Gnuplotting/gnuplot-palettes
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_set1_line(void)
            {
                send2gp("set linetype 1 lc rgb '#E41A1C' ");// red
                send2gp("set linetype 2 lc rgb '#377EB8' ");// blue
                send2gp("set linetype 3 lc rgb '#4DAF4A' ");// green
                send2gp("set linetype 4 lc rgb '#984EA3' ");// purple
                send2gp("set linetype 5 lc rgb '#FF7F00' ");// orange
                send2gp("set linetype 6 lc rgb '#FFFF33' ");// yellow
                send2gp("set linetype 7 lc rgb '#A65628' ");// brown
                send2gp("set linetype 8 lc rgb '#F781BF' ");// pink

                send2gp("set palette maxcolors 8");
                char str[] ="set palette defined ( \
                      0 '#E41A1C',\
                      1 '#377EB8',\
                      2 '#4DAF4A',\
                      3 '#984EA3',\
                      4 '#FF7F00',\
                      5 '#FFFF33',\
                      6 '#A65628',\
                      7 '#F781BF')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to Matlab 'jet'
            /// Data from https://github.com/Gnuplotting/gnuplot-palettes
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_jet_palette(void)
            {
                char str[] ="set palette defined ( \
                      0 '#000090',\
                      1 '#000fff',\
                      2 '#0090ff',\
                      3 '#0fffee',\
                      4 '#90ff70',\
                      5 '#ffee00',\
                      6 '#ff7000',\
                      7 '#ee0000',\
                      8 '#7f0000')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to Matlab 'parula'
            /// Data from https://github.com/Gnuplotting/gnuplot-palettes
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_parula_palette(void)
            {
                char str[] ="set palette defined (\
                      0 '#352a87',\
                      1 '#0363e1',\
                      2 '#1485d4',\
                      3 '#06a7c6',\
                      4 '#38b99e',\
                      5 '#92bf73',\
                      6 '#d9ba56',\
                      7 '#fcce2e',\
                      8 '#f9fb0e')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to 'cool-warm'
            // See http://www.kennethmoreland.com/color-advice/
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_coolwarm_palette(void)
            {
                char str[] = "set palette defined (\
                      0 '#5548C1', \
                      1 '#7D87EF', \
                      2 '#A6B9FF', \
                      3 '#CDD7F0', \
                      4 '#EBD1C2', \
                      5 '#F3A889', \
                      6 '#DE6A53', \
                      7 '#B10127')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set palette to 'black body'
            // See http://www.kennethmoreland.com/color-advice/
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_blackbody_palette(void)
            {
                char str[] = "set palette defined (\
                      0 '#000000', \
                      1 '#2B0F6B', \
                      2 '#5D00CB', \
                      3 '#C60074', \
                      4 '#EB533C', \
                      5 '#F59730', \
                      6 '#E9D839', \
                      7 '#FFFFFF')";
                send2gp(str);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Save plot to file.
            /// @param name filename
            ///
            /// Extensions that are supported:
            /// - png
            /// - ps
            /// - eps
            /// - tex
            /// - pdf
            /// - svg
            /// - emf
            /// - gif
            ///
            /// \note When 'latex' output is used the '\' must be escaped by '\\\\' e.g set_xlabel("Frequency $\\\\omega = 2 \\\\pi f$")
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_output(const char* name)
            {
                std::string name_s(name);
                size_t found = name_s.find_last_of(".");
                std::string ext;
                ext = name_s.substr(found + 1);
                std::ostringstream tmp_s;

                if (ext.compare("png")==0)
                {
                    tmp_s << "set terminal pngcairo enhanced font 'Verdana,10'";
                }
                else if (ext.compare("ps") == 0)
                {
                    tmp_s << "set terminal postscript enhanced color";
                }
                else if (ext.compare("eps") == 0)
                {
                    tmp_s << "set terminal postscript eps enhanced color";
                }
                else if (ext.compare("tex") == 0)
                {
                    tmp_s << "set terminal cairolatex eps color enhanced";
                }
                else if (ext.compare("pdf") == 0)
                {
                    tmp_s << "set terminal pdfcairo color enhanced";
                }
                else if (ext.compare("svg") == 0)
                {
                    tmp_s << "set terminal svg enhanced";
                }
                else if (ext.compare("emf") == 0)
                {
                    tmp_s << "set terminal emf color enhanced";
                }
                else if (ext.compare("gif") == 0)
                {
                    tmp_s << "set terminal gif enhanced";
                }
                //else if (ext.compare("jpg") == 0)
                //{
                //	tmp_s << "set terminal jpeg ";
                //}
                else
                {
                    tmp_s << "set terminal " << term;
                }
                std::string s = tmp_s.str();
                send2gp(s.c_str());
                tmp_s.str("");  // Clear buffer
                tmp_s << "set output '" << name_s << "'";
                s = tmp_s.str();
                send2gp(s.c_str());
            }


            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Reset output terminal.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void reset_term(void)
            {
                send2gp("reset session");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Set output terminal.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_term(const char* ttype)
            {
                std::ostringstream tmp_s;
                tmp_s << "set terminal " << ttype;
                std::string s = tmp_s.str();
                term = s;
                send2gp(s.c_str());
            }

    }; // End Gnuplot Class


    /// @}
        ///
    /// @defgroup parser Parser
    /// \brief Parameter file parser functions.
    /// @{

    ///
    /// \brief A parser class.
    ///
    /// Implements parsing from text file for different types
    ///
    class parser
    {
        private:
            std::map<std::string, std::string> par_map; ///< Map structure to store parameter and value

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Check if key is valid.
            /// @returns TRUE if key is valid FALSE otherwise
            /// @param key Key string
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool valid_key(const std::string& key)
            {
                if(par_map.find(key) == par_map.end() )
                {
                    std::cout << "SigPack: Parameter "+ key + " not found!" <<std::endl;
                    return false;
                }
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Parse complex value.
            ///
            /// String may be in format e.g. "1+1i","1+1j" or entirely real or imaginary "1", "1i" or "1j"
            /// @returns Complex value
            /// @param str Complex notation string
            ////////////////////////////////////////////////////////////////////////////////////////////
            std::complex<double> parse_cx(std::string str)
            {
                double re,im;
                char i_ch;
                std::stringstream iss(str);

                // Parse full ...
                if(iss >> re >> im >> i_ch && (i_ch=='i'|| i_ch=='j')) return std::complex<double>(re,im);

                // ... or only imag
                iss.clear();
                iss.seekg(0,iss.beg);
                if(iss >> im >> i_ch && (i_ch=='i'|| i_ch=='j')) return std::complex<double>(0.0,im);

                // .. or only real
                iss.clear();
                iss.seekg(0,iss.beg);
                if(iss >> re) return std::complex<double>(re,0.0);

                // ... otherwise
                err_handler("Could not parse complex number!");
            }

        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ///
            /// Opens parameter file and puts the key and value in the map structure
            /// @param fname Parameter file name
            ////////////////////////////////////////////////////////////////////////////////////////////
            parser(const std::string& fname)
            {
                std::ifstream fh;
                std::string line;
                size_t mark = 0;

                // Clear parameter map
                par_map.clear();

                // Open file
                fh.open(fname.c_str());
                if (!fh)
                {
                    wrn_handler("Could not find " + fname);
                }
                else
                {
                    // Parse
                    while (std::getline(fh, line))
                    {
                        std::string keyS="";
                        std::string dataS="";

                        // Skip empty lines
                        if (line.empty())
                            continue;

                        // Skip lines with only whitespace
                        if(line.find_first_not_of("\t ")==std::string::npos)
                            continue;

                        // Remove comment
                        mark = line.find("%");
                        if(mark!=std::string::npos)
                            line.erase(mark,line.length());

                        // Do we have a '='
                        mark = line.find("=");
                        if(mark!=std::string::npos)
                        {
                            // Find key
                            keyS = line.substr(line.find_first_not_of("\t "),mark-line.find_first_not_of("\t "));
                            keyS = keyS.substr(0,keyS.find_last_not_of("\t ")+1);

                            // Find data
                            dataS = line.substr(mark+1,line.length());
                            dataS = dataS.substr(0,dataS.find_last_not_of("\t ")+1);

                            // Do we have a string
                            mark = dataS.find("\"");
                            if(mark!=std::string::npos)
                            {
                                dataS = dataS.substr(mark+1,dataS.length());
                                dataS = dataS.substr(0,dataS.find_last_of("\""));
                            }
                            // Do we have a vector/matrix
                            mark = dataS.find("[");
                            if(mark!=std::string::npos)
                            {
                                dataS = dataS.substr(mark+1,dataS.length());
                                dataS = dataS.substr(0,dataS.find_last_of("]"));
                            }

                            // Insert to map
                            par_map.insert(std::pair<std::string, std::string>(keyS, dataS));
                        }
                    }

                    // Close file
                    fh.close();
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~parser() {}

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Generic type get function.
            /// @returns Value of type T, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            T  getParam(const std::string key,const T def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default "<< def_val  <<std::endl;
                    return def_val;
                }
                std::istringstream iss(par_map.find(key)->second);
                T out;
                iss >> out;
                return out;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief String type get function.
            /// @returns String value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            std::string getString(const std::string key,const std::string def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default " << def_val <<std::endl;
                    return def_val;
                }
                return par_map.find(key)->second;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Col type get function.
            /// @returns Col value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            arma::Col<T> getCol(const std::string key,const arma::Col<T> def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string row,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),';')+1);
                arma::Col<T> x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, row, ';');
                    std::stringstream iss(row);
                    iss >> x(k);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief cx_vec type get function.
            /// @returns cx_vec value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_vec getCxCol(const std::string key,const arma::cx_vec def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string row,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),';')+1);
                arma::cx_vec x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, row, ';');
                    x(k) = parse_cx(row);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Row type get function.
            /// @returns Row value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            arma::Row<T> getRow(const std::string key,const arma::Row<T> def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string col,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
                arma::Row<T> x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, col, ',');
                    std::stringstream iss(col);
                    iss >> x(k);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief cx_rowvec type get function.
            /// @returns cx_rowvec value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_rowvec getCxRow(const std::string key,const arma::cx_rowvec def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }

                std::string col,str=par_map.find(key)->second;
                std::istringstream full_str(str);
                int K = static_cast<int>(std::count(str.begin(),str.end(),',')+1);
                arma::cx_rowvec x(K);
                for(int k=0; k<K; k++)
                {
                    std::getline(full_str, col, ',');
                    x(k) = parse_cx(col);
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Mat type get function.
            /// @returns Mat value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            template <typename T>
            arma::Mat<T> getMat(const std::string key,const arma::Mat<T> def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }
                std::string full_str,row,col;
                std::istringstream iss_full;

                full_str = par_map.find(key)->second;
                int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

                iss_full.str(full_str);
                std::getline(iss_full, row, ';');
                int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

                arma::Mat<T> x(R,C);

                iss_full.seekg(0,iss_full.beg);
                for(int r=0; r<R; r++)
                {
                    std::getline(iss_full, row, ';');
                    std::istringstream iss_row(row);
                    for(int c=0; c<C; c++)
                    {
                        std::getline(iss_row, col, ',');
                        std::istringstream iss_col(col);
                        iss_col >> x(r,c);
                    }
                }
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief cx_mat type get function.
            /// @returns cx_mat value, if it not exists in file the default value is returned
            /// @param key Parameter name
            /// @param def_val Default value if key was not found in file
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_mat getCxMat(const std::string key,const arma::cx_mat def_val)
            {
                if(!valid_key(key))
                {
                    std::cout << "Setting default \n"<< def_val  <<std::endl;
                    return def_val;
                }
                std::string full_str,row,col;
                std::istringstream iss_full;

                full_str = par_map.find(key)->second;
                int R = static_cast<int>(std::count(full_str.begin(),full_str.end(),';')+1);

                iss_full.str(full_str);
                std::getline(iss_full, row, ';');
                int C = static_cast<int>(std::count(row.begin(),row.end(),',')+1);

                arma::cx_mat x(R,C);

                iss_full.seekg(0,iss_full.beg);
                for(int r=0; r<R; r++)
                {
                    std::getline(iss_full, row, ';');
                    std::istringstream iss_row(row);
                    for(int c=0; c<C; c++)
                    {
                        std::getline(iss_row, col, ',');
                        x(r,c)=parse_cx(col);
                    }
                }
                return x;
            }

    }; // end Class
    /// @}
       /// @}
        ///
    /// @defgroup image Image
    /// \brief Image functions.
    /// @{

    ///
    /// \brief Portable anymap format class.
    ///
    /// Implements portable anymap image functions
    /// Supports .pbm, .pgm and .ppm plain and raw
    ///

    class PNM
    {
        private:
            std::ifstream ifs;          ///< Input stream handle
            std::ofstream ofs;          ///< Output stream handle
            arma::uword cols;                   ///< Nr of columns in image
            arma::uword rows;                   ///< Nr of rows in image
            int maxval;                 ///< Maximum pixel value in image
        public:
            enum imtype { NOTUSED, PBM_A, PGM_A, PPM_A, PBM_B, PGM_B, PPM_B } type; ///< Image format

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            PNM()
            {
                type   = NOTUSED;
                cols   = 0;
                rows   = 0;
                maxval = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~PNM() {}

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Clears the internal variables.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void clear(void)
            {
                type   = NOTUSED;
                cols   = 0;
                rows   = 0;
                maxval = 0;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Reads the .pnm header.
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////
            void read_header()
            {
                std::string str;
                while (ifs >> str ) // Read until we have maxval
                {
                    // Remove comments
                    size_t mark = str.find_first_of("#");
                    if ( mark!= std::string::npos)
                    {
                        ifs.ignore(256, '\n');
                        str.erase(mark, 1);

                        if (str.empty()) continue;
                    }

                    if (type == NOTUSED)
                    {
                        type = static_cast<imtype>(str.at(1)-48);  // Conv char to int
                        if(str.at(0)!='P' || type>PPM_B) err_handler("Wrong type!");
                    }
                    else if (cols   == 0)
                    {
                        cols = atoi(str.c_str());
                    }
                    else if (rows   == 0)
                    {
                        rows = atoi(str.c_str());
                        if(type==PBM_A || type==PBM_B)
                        {
                            maxval = 1;  // Set maxval=1 for .PBM types
                            break;
                        }
                    }
                    else if (maxval == 0)
                    {
                        maxval = atoi(str.c_str());
                        break;
                    }
                }
                ifs.ignore(1); // Skip one char before data
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Write the .pnm header.
            ///
            /// @param _type     Image type
            /// @param _rows     Nr of rows
            /// @param _cols     Nr of cols
            /// @param _maxval   Maxval
            /// @param comments Comments
            ////////////////////////////////////////////////////////////////////////////////////////////
            void write_header(const imtype _type, const arma::uword _rows, const arma::uword _cols, const int _maxval, const std::string comments)
            {
                type   = _type;
                rows   = _rows;
                cols   = _cols;
                maxval = _maxval;
                ofs << "P" << type << std::endl;
                ofs << "# " << comments << std::endl;
                ofs << cols << " " << rows << std::endl;
                if(!((type==PBM_A) || (type==PBM_B))) ofs << maxval << std::endl;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Write the .pnm file.
            ///
            /// @returns true if success
            /// @param fname  File name
            /// @param _type  File name
            /// @param img    Image data
            /// @param info   File comments
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool write(std::string fname, const imtype _type, const arma::cube& img, const std::string info="")
            {
                // Open file
                ofs.open(fname.c_str(), std::ofstream::binary);
                if (!ofs.good())
                {
                    std::cout << "Could not open " << fname << std::endl;
                    return false;
                }
                else
                {
                    write_header(_type, img.n_rows, img.n_cols, (int)img.max(),info);
                    //            get_info();

                    // Write data
                    if(type==PPM_A ) // Plain (ASCII )type
                    {
                        for(arma::uword r=0; r<rows; r++)
                        {
                            arma::uword i = 0;
                            for(arma::uword c=0; c<cols; c++)
                            {
                                ofs << img(r,c,0) << " " << img(r,c,1) << " " << img(r,c,2) << " "; // R G B
                                if(++i%5==0) ofs << std::endl; // Max len is 70 chars/line
                            }
                        }
                    }
                    else if(type==PPM_B)
                    {
                        for(arma::uword r=0; r<rows; r++)
                            for(arma::uword c=0; c<cols; c++)
                            {
                                unsigned char bb;
                                bb= static_cast<unsigned char>(img(r,c,0));   // R
                                ofs.write(reinterpret_cast<char*>(&bb),1);
                                bb= static_cast<unsigned char>(img(r,c,1));   // G
                                ofs.write(reinterpret_cast<char*>(&bb),1);
                                bb= static_cast<unsigned char>(img(r,c,2));   // B
                                ofs.write(reinterpret_cast<char*>(&bb),1);
                            }
                    }
                }
                ofs.close();
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Write the .pnm file.
            ///
            /// @returns true if success
            /// @param fname  File name
            /// @param _type  Image type
            /// @param img    Image data
            /// @param info   File comments
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool write(std::string fname, const imtype _type, arma::mat& img, const std::string info="")
            {
                // Open file
                ofs.open(fname.c_str(), std::ofstream::binary);
                if (!ofs.good())
                {
                    std::cout << "Could not open " << fname << std::endl;
                    return false;
                }
                else
                {
                    write_header(_type, img.n_rows, img.n_cols, (int)img.max(),info);
                    //            get_info();

                    // Write data
                    if(type==PBM_A || type ==PGM_A ) // Plain (ASCII )type
                    {
                        arma::uword i=0;
                        for(arma::mat::iterator ii=img.begin(); ii!=img.end(); ++ii)
                        {
                            ofs << *ii << " ";
                            if(++i%11==0) ofs << std::endl; // Max len is 70 chars/line
                        }
                    }
                    else if(type == PBM_B) // Raw .pbm
                    {
                        std::bitset<8> b;
                        for(arma::uword r=0; r<rows; r++)
                            for(arma::uword c=0; c<cols; c++)
                            {
                                arma::uword ix = 7-(c%8);
                                b[ix] =  (img(r,c)>0);
                                if(ix==0 || c==cols-1)
                                {
                                    ofs.write(reinterpret_cast<char*>(&b),1);
                                    b.reset();
                                }
                            }
                    }
                    else if(type == PGM_B) // Raw .pgm
                    {
                        if(maxval<=255)
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned char bb= static_cast<unsigned char>(img(r,c));
                                    ofs.write(reinterpret_cast<char*>(&bb),1);
                                }
                        }
                        else
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned int bb;
                                    bb = ((static_cast<unsigned int>(img(r,c)))>> 8) & 0x00ff;   // Write MSB first
                                    ofs.write(reinterpret_cast<char*>(&bb),1);
                                    bb = static_cast<unsigned int>(img(r,c)) & 0x00ff;
                                    ofs.write(reinterpret_cast<char*>(&bb),1);
                                }
                        }
                    }

                }

                ofs.close();
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Prints header info
            ///
            ////////////////////////////////////////////////////////////////////////////////////////////
            void get_info()
            {
                std::cout << "Type:  P"  << type   << std::endl;
                std::cout << "cols:   "  << cols   << std::endl;
                std::cout << "rows:   "  << rows   << std::endl;
                if(type==PGM_A || type==PGM_B) std::cout << "Maxval: "  << maxval << std::endl;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get nr of rows
            ///
            /// @returns number of rows
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::uword get_rows()
            {
                return rows;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get nr of cols
            ///
            /// @returns number of columns
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::uword get_cols()
            {
                return cols;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Get maxval
            ///
            /// @returns Maximum value in image
            ////////////////////////////////////////////////////////////////////////////////////////////
            int get_maxval()
            {
                return maxval;
            }


            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Read image
            ///
            /// @returns true if success
            /// @param fname File name
            /// @param img   Image data
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool read(std::string fname, arma::cube& img)
            {
                // Open file
                ifs.open(fname.c_str(), std::ifstream::binary);
                if (!ifs.good())
                {
                    std::cout << "Could not open " << fname << std::endl;
                    return false;
                }
                else
                {
                    read_header();
                    //            get_info();

                    img.set_size(rows,cols,3);
                    arma::uword r = 0, c = 0;
                    // Get the data
                    if (type==PPM_A )  // Plain .PPM
                    {
                        std::string str;
                        arma::uword i=0;
                        while (ifs >> str && r<rows ) // Read until eof
                        {
                            // Remove comments
                            size_t mark = str.find_first_of("#");
                            if ( mark!= std::string::npos)
                            {
                                ifs.ignore(256, '\n');
                                str.erase(mark, 1);

                                if (str.empty()) continue;
                            }
                            int pix= atoi(str.c_str());  // Convert to int

                            img(r, c, i%3) = pix;
                            i++;
                            if(i%3==0)
                                if (++c == cols)
                                {
                                    c = 0;
                                    r++;
                                }
                        }
                    }
                    else if (type==PPM_B )  // Raw .PPM
                    {
                        for(arma::uword r=0; r<rows; r++)
                            for(arma::uword c=0; c<cols; c++)
                            {
                                unsigned char bb;
                                ifs.read(reinterpret_cast<char*>(&bb),1);    // R
                                img(r,c,0) = bb;
                                ifs.read(reinterpret_cast<char*>(&bb),1);    // G
                                img(r,c,1) = bb;
                                ifs.read(reinterpret_cast<char*>(&bb),1);    // B
                                img(r,c,2) = bb;
                            }
                    }
                }
                ifs.close();
                return true;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Read image
            ///
            /// @returns true if success
            /// @param fname File name
            /// @param img Image data
            ////////////////////////////////////////////////////////////////////////////////////////////
            bool read(std::string fname, arma::mat& img)
            {
                // Open file
                ifs.open(fname.c_str(), std::ifstream::binary);
                if (!ifs.good())
                {
                    std::cout << "Could not open " << fname << std::endl;
                    return false;
                }
                else
                {
                    read_header();
                    get_info();

                    img.set_size(rows,cols);
                    arma::uword r = 0, c = 0;
                    // Get the data
                    if (type==PBM_A || type == PGM_A)  // Plain .PGM or .PBM
                    {
                        std::string str;
                        while (ifs >> str && r<rows ) // Read until eof
                        {
                            // Remove comments
                            size_t mark = str.find_first_of("#");
                            if ( mark!= std::string::npos)
                            {
                                ifs.ignore(256, '\n');
                                str.erase(mark, 1);

                                if (str.empty()) continue;
                            }
                            int pix= atoi(str.c_str());  // Convert to int
                            img(r, c) = pix;

                            if (++c == cols)
                            {
                                c = 0;
                                r++;
                            }
                        }
                    }
                    else if(type== PBM_B) // Raw PBM
                    {
                        unsigned char ch;
                        while (ifs.read(reinterpret_cast<char*>(&ch),1) && r<rows)  // Read until eof
                        {
                            std::bitset<8> pix(ch);
                            for(int b=7; b>=0; b--)
                            {
                                img(r,c) = pix[b];
                                if (++c >= cols)
                                {
                                    c = 0;
                                    r++;
                                    break;
                                }
                            }
                        }
                    }
                    else if(type==PGM_B) // Raw PGM
                    {
                        if(maxval<=255)
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned char bb;
                                    ifs.read(reinterpret_cast<char*>(&bb),1);
                                    img(r,c) = bb;
                                }
                        }
                        else
                        {
                            for(arma::uword r=0; r<rows; r++)
                                for(arma::uword c=0; c<cols; c++)
                                {
                                    unsigned char bb[2];
                                    ifs.read(reinterpret_cast<char*>(bb),2);
                                    img(r,c) = (bb[0]<<8)+bb[1];
                                }
                        }
                    }

                }
                ifs.close();
                return true;
            }
    };
    /// @}
        #define FCN_XUW [=](arma::mat x,arma::mat u,arma::mat w)     // Lambda function f(x,u,w) ([capture] by copy)
    using fcn_t = std::function< double(arma::mat,arma::mat,arma::mat) >;
    using fcn_v = std::vector<fcn_t>;
    using fcn_m = std::vector<fcn_v>;

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Evaluate function f(x,u,w)
    /// @param  f Function pointer vector
    /// @param  x Input vector
    /// @param  u Input vector
    /// @param  w Input vector
    /// @return y Output column vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::mat eval_fcn( const fcn_v f, const arma::mat& x, const arma::mat& u, const arma::mat& w)
    {
        arma::mat y((arma::uword)(f.size()),1,arma::fill::zeros);
        for( arma::uword n=0; n<y.n_rows;n++)
            y(n,0)    = f[n](x,u,w);
        return y;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Evaluate function f(x,u,w=0)
    /// @param  f Function pointer vector
    /// @param  x Input vector
    /// @param  u Input vector
    /// @return y Output column vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::mat eval_fcn( const fcn_v f, const arma::mat& x, const arma::mat& u)
    {
        arma::mat w0(0,0,arma::fill::zeros);
        return eval_fcn(f,x,u,w0);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Evaluate function f(x,u=0,w=0)
    /// @param  f Function pointer vector
    /// @param  x Input vector
    /// @return y Output column vector
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline arma::mat eval_fcn( const fcn_v f, const arma::mat& x)
    {
        arma::mat w0(0,0,arma::fill::zeros);
        arma::mat u0(w0);
        return eval_fcn(f,x,u0,w0);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////
    /// \brief Discretize to a state transition and state noise cov matrix from an LTI system.
    ///
    ///  Discretize a continious LTI system using van Loan method. The model is in the form
    ///
    ///     dx/dt = F x + Ww,  w ~ N(0,Qc)
    ///
    ///  Result of discretization is the model
    ///
    ///     x[k] = A x[k-1] + q, q ~ N(0,Q)
    ///
    ///  See http://becs.aalto.fi/en/research/bayes/ekfukf/
    /// @param F  LTI system model matrix
    /// @param W  LTI system noise model matrix
    /// @param Qc LTI power spectra density matrix
    /// @param dT Discretization delta time
    /// @param A  Output - discrete system model
    /// @param Q  Output - discrete state noise cov matrix
    ////////////////////////////////////////////////////////////////////////////////////////////
    arma_inline void lti2discr(const arma::mat& F,const arma::mat& W, const arma::mat& Qc, const double dT, arma::mat& A,arma::mat& Q)
    {
        arma::uword M = F.n_rows;

        // Solve A
        A = arma::expmat(F*dT);

        // Solve Q by using matrix fraction decomposition
        arma::mat AB = arma::zeros(2*M,M);
        arma::mat CD = arma::zeros(2*M,2*M);
        arma::mat EF = arma::zeros(2*M,M);
        EF.submat(M,0, 2*M-1,M-1)  = arma::eye(M,M);
        CD.submat(0,0, M-1,M-1)    = F;
        CD.submat(M,M,2*M-1,2*M-1) = -F.t();
        CD.submat(0,M,M-1,2*M-1)   = W*Qc*W.t();

        AB = arma::expmat(CD*dT)*EF;

        Q = AB.rows(0,M-1)*arma::inv(AB.rows(M,2*M-1));
    }

    ///
    /// @defgroup kalman Kalman
    /// \brief Kalman predictor/filter functions.
    /// @{

    ////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// \brief Kalman filter class.
    ///
    /// Implements Kalman functions for the discrete system
    /// \f[ x_k = Ax_{k-1}+Bu_{k-1} + w_{k-1}  \f]
    /// with measurements
    /// \f[ z_k = Hx_k + v_k  \f]
    /// The predicting stage is
    /// \f[ \hat{x}^-_k = A\hat{x}_{k-1}+Bu_{k-1} \f]
    /// \f[ P^-_k = AP_{k-1}A^T+Q \f]
    /// and the updates stage
    /// \f[ K_k = P^-_kH^T(HP^-_kH^T+R)^{-1} \f]
    /// \f[ \hat{x}_k = \hat{x}^-_k + K_k(z_k-H\hat{x}^-_k) \f]
    /// \f[ P_k = (I-K_kH)P^-_k \f]
    ///
    /// For detailed info see: http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
    ////////////////////////////////////////////////////////////////////////////////////////////
    class KF
    {
        protected:
            arma::uword N;        ///< Number of states
            arma::uword M;        ///< Number of inputs
            arma::uword L;        ///< Number of measurements/observations
            bool lin_proc;        ///< Linearity flag for process
            bool lin_meas;        ///< Linearity flag for measurement
            arma::mat x;          ///< State vector
            arma::mat z_err;      ///< Prediction error
            arma::mat A;          ///< State transition matrix
            arma::mat B;          ///< Input matrix
            arma::mat H;          ///< Measurement matrix
            arma::mat P;          ///< Error covariance matrix (estimated accuracy)
            arma::mat Q;          ///< Process noise
            arma::mat R;          ///< Measurement noise
            arma::mat K;          ///< Kalman gain vector
            fcn_v f;              ///< Vector of Kalman state transition functions
            fcn_v h;              ///< Vector of Kalman measurement functions
        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            KF(arma::uword _N,arma::uword _M,arma::uword _L)
            {
                N = _N;   // Nr of states
                M = _M;   // Nr of measurements/observations
                L = _L;   // Nr of inputs
                lin_proc = true;
                lin_meas = true;
                x.set_size(N,1); x.zeros();
                z_err.set_size(M,1); z_err.zeros();
                A.set_size(N,N); A.eye();
                B.set_size(N,L); B.zeros();
                H.set_size(M,N); H.zeros();
                P.set_size(N,N); P.eye();
                Q.set_size(N,N); Q.eye();
                R.set_size(M,M); R.eye();
                K.set_size(N,M); K.zeros();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~KF() {}

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Clear the internal states and pointer.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void clear(void)
            {
                K.zeros();
                P.eye();
                x.zeros();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            //  Set/get functions
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_state_vec(const arma::mat& _x)       { x = _x;   }    // Set state vector.[Nx1]
            void set_trans_mat(const arma::mat& _A)       { A = _A;   }    // Set state transition matrix.[NxN]
            void set_control_mat(const arma::mat& _B)     { B = _B;   }    // Set input matrix.[NxL]
            void set_meas_mat(const arma::mat& _H)        { H = _H;   }    // Set measurement matrix.[MxN]
            void set_err_cov(const arma::mat& _P)         { P = _P;   }    // Set error covariance matrix.[NxN]
            void set_proc_noise(const arma::mat& _Q)      { Q = _Q;   }    // Set process noise cov. matrix.
            void set_meas_noise(const arma::mat& _R)      { R = _R;   }    // Set meas noise cov. matrix.
            void set_kalman_gain(const arma::mat& _K)     { K = _K;   }    // Set Kalman gain matrix.[NxM]
            void set_trans_fcn(fcn_v _f)    // Set state transition functions
            {
                f = _f;
                lin_proc = false;
            }
            void set_meas_fcn(fcn_v _h)     // Set measurement functions
            {
                h = _h;
                lin_meas = false;
            }

            arma::mat get_state_vec(void)        { return x;       }   // Get states [Nx1]
            arma::mat get_err(void)              { return z_err;   }   // Get pred error [Mx1]
            arma::mat get_kalman_gain(void)      { return K;       }   // Get Kalman gain [NxM]
            arma::mat get_err_cov(void)          { return P;       }   // Get error cov matrix [NxN]

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Predict the internal states using a control input.
            /// @param u Input/control signal
            ////////////////////////////////////////////////////////////////////////////////////////////
            void predict(const arma::mat u )
            {
                x = A*x+B*u;      // New state
                P = A*P*A.t()+Q;  // New error covariance
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Predict the internal states, no control.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void predict(void)
            {
                 arma::mat u0(L,1,arma::fill::zeros);
                 predict(u0);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Correct and update the internal states.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void update(const arma::mat z )
            {
                // Compute the Kalman gain
                K = P*H.t()*inv(H*P*H.t()+R);

                // Update estimate with measurement z
                z_err = z-H*x;
                x = x+K*z_err;

                // Joseph’s form covariance update
                arma::mat Jf = arma::eye<arma::mat>(N,N)-K*H;
                P = Jf*P*Jf.t() + K*R*K.t();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Rauch-Tung-Striebel smoother.
            /// See http://www.lce.hut.fi/~ssarkka/course_k2011/pdf/course_booklet_2011.pdf
            ////////////////////////////////////////////////////////////////////////////////////////////
            void rts_smooth(const arma::mat& Xf, const arma::cube& Pf, arma::mat& Xs, arma::cube& Ps )
            {
                arma::uword Nf = Xf.n_cols;
                arma::mat X_pred(N,1);
                arma::mat P_pred(N,N);
                arma::mat C(N,N);

                // Copy forward data
                Xs = Xf;
                Ps = Pf;

                // Backward filter
                for(arma::uword n=Nf-2; n>0; n--)
                {
                    // Project state and error covariance
                    X_pred = A*Xf.col(n);
                    P_pred = A*Pf.slice(n)*A.t()+Q;

                    // Update
                    C = Pf.slice(n)*A.t()*inv(P_pred);
                    Xs.col(n)   += C*(Xs.col(n+1)-X_pred);
                    Ps.slice(n) += C*(Ps.slice(n+1)-P_pred)*C.t();
                }
            }
    }; // End class KF

    ////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// \brief First order Extended Kalman filter class
    ///
    /// Implements Kalman functions for the discrete system with additive noise
    /// \f[ x_k = f(x_{k-1})+Bu_{k-1} + w_{k-1}  \f]
    /// and with measurements
    /// \f[ z_k = h(x_k) + v_k  \f]
    /// where f(x) and h(x) may be nonlinear functions.
    /// The predicting stage is
    /// \f[ \hat{x}^-_k = A\hat{x}_{k-1}+Bu_{k-1} \f]
    /// \f[ P^-_k = AP_{k-1}A^T+Q \f]
    /// and the updates stage
    /// \f[ K_k = P^-_kH^T(HP^-_kH^T+R)^{-1} \f]
    /// \f[ \hat{x}_k = \hat{x}^-_k + K_k(z_k-H\hat{x}^-_k) \f]
    /// \f[ P_k = (I-K_kH)P^-_k \f]
    ///
    /// Where A and H is the Jacobians of f(x) and h(x) functions.
    ///
    /// For detailed info see: http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf
    ////////////////////////////////////////////////////////////////////////////////////////////
    class EKF: public KF
    {
        protected:
            fcn_m f_jac;            ///< Matrix of Extended Kalman state transition jacobian
            fcn_m h_jac;            ///< Matrix of Extended Kalman measurement transition jacobian
            double dx;              ///< Finite difference approximation step size
        public:

            EKF(arma::uword _N,arma::uword _M,arma::uword _L): KF(_N,_M,_L)
            {
                dx = 1e-7;  // Default diff step size
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            //  Set/get functions
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_diff_step(double _dx)  { dx  = _dx;    }   // Set diff step size
            void set_state_jac(fcn_m _f)    { f_jac = _f;   }   // Set EKF state transition jacobian
            void set_meas_jac(fcn_m _h)     { h_jac = _h;   }   // Set EKF measurement transition jacobian

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Calculate and evaluate Jacobian matrix using finite difference approximation
            /// @param  _F Jacobian matrix d/dx evaluated at x
            /// @param  _f Function vector [ f0(x,u,w) ...  fN(x,u,w)]
            /// @param  _x State vector
            ///
            /// Alternative: Complex Step Diff: http://blogs.mathworks.com/cleve/2013/10/14/complex-step-differentiation/
            ////////////////////////////////////////////////////////////////////////////////////////////
            void jacobian_diff(arma::mat& _F, fcn_v _f, const arma::mat& _x)
            {
                arma::uword nC = _F.n_cols;
                arma::uword nR = static_cast<arma::uword>(_f.size());
                arma::mat z0(nC,1,arma::fill::zeros); // Zero matrix, assume dim u and w <= states

                if(nR==0 || nR!=_F.n_rows) err_handler("Function list is empty or wrong size");

                for(arma::uword c=0; c<nC; c++)
                {
                    arma::mat xp(_x);
                    arma::mat xm(_x);
                    xp(c,0) += dx;
                    xm(c,0) -= dx;

                    // Finite diff approx, evaluate at x,u=0,w=0
                    for(arma::uword r=0; r<nR; r++)
                        _F(r,c) = (_f[r](xp,z0,z0)-_f[r](xm,z0,z0))/(2*dx);
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Calculate and evaluate Jacobian matrix using finite difference approximation
            /// @param  _F Jacobian matrix d/dx evaluated at x
            /// @param  _f Function vector [ f0(x,u,w) ...  fN(x,u,w)]
            ////////////////////////////////////////////////////////////////////////////////////////////
            void jacobian_diff(arma::mat& _F, fcn_v _f)
            {
               jacobian_diff(_F,_f,x); // Use current state from object
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Evaluate Jacobian matrix using analytical jacobian
            /// @param  _F Jacobian matrix d/dx evaluated at x
            /// @param  _f_m Jacobian function matrix
            /// @param  _x State vector
            ////////////////////////////////////////////////////////////////////////////////////////////
            void jacobian_analytical(arma::mat& _F, fcn_m _f_m, const arma::mat& _x)
            {
                arma::uword nC = _F.n_cols;
                arma::uword nR = static_cast<arma::uword>(_f_m.size());

                if(nR==0 || nR!=_F.n_rows) err_handler("Function list is empty or wrong size");

                arma::mat z0(nC,1,arma::fill::zeros); // Zero matrix, assume dim u and w <= states
                for(arma::uword c=0; c<nC; c++)
                {
                    for(arma::uword r=0; r<nR; r++)
                        _F(r,c) = _f_m[r][c](_x,z0,z0);
                }
            }
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Evaluate Jacobian matrix using analytical jacobian
            /// @param  _F Jacobian matrix d/dx evaluated at x
            /// @param  _f_m Jacobian function matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void jacobian_analytical(arma::mat& _F, fcn_m _f_m)
            {
               jacobian_analytical(_F,_f_m,x); // Use current state from object
            }


            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Predict the internal states using a control input.
            /// @param u Input/control signal
            ////////////////////////////////////////////////////////////////////////////////////////////
            void predict(const arma::mat u )
            {
                if( !lin_proc )
                {
                    // Update A with jacobian approx or analytical if set
                    if(f_jac.size()>0)
                        jacobian_analytical(A,f_jac);
                    else
                        jacobian_diff(A,f);

                    // Predict state   x+ = f(x,u,0)
                    x = eval_fcn(f,x,u);
                }
                else  // Linear process
                    x = A*x+B*u;

                // Project error covariance
                P = A*P*A.t()+Q;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Predict the internal states, no control.
            ////////////////////////////////////////////////////////////////////////////////////////////
            void predict(void)
            {
                 arma::mat u0(L,1,arma::fill::zeros);
                 predict(u0);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Correct and update the internal states. EKF
            ////////////////////////////////////////////////////////////////////////////////////////////
            void update(const arma::mat z )
            {
                arma::mat z_hat(M,1);

                if(!lin_meas) // Nonlinear measurement
                {
                    // Update H with jacobian approx or analytical if set
                    if( h_jac.size()>0)
                        jacobian_analytical(H,h_jac);
                    else
                        jacobian_diff(H,h);

                    // Update measurement
                    z_hat = eval_fcn(h,x);
                }
                else  // Linear meas
                    z_hat = H*x;

                // Calc residual
                z_err = z-z_hat;

                // Compute the Kalman gain
                K = P*H.t()*inv(H*P*H.t()+R);

                // Update estimate with measurement residual
                x = x+K*z_err;

                // Joseph’s form covariance update
                arma::mat Jf = arma::eye<arma::mat>(N,N)-K*H;
                P = Jf*P*Jf.t()+K*R*K.t();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Rauch-Tung-Striebel smoother.
            /// See http://www.lce.hut.fi/~ssarkka/course_k2011/pdf/course_booklet_2011.pdf
            ////////////////////////////////////////////////////////////////////////////////////////////
            void rts_smooth(const arma::mat& Xf, const arma::cube& Pf, arma::mat& Xs, arma::cube& Ps )
            {
                arma::uword Nf = Xf.n_cols;
                arma::mat X_pred(N,1);
                arma::mat P_pred(N,N);
                arma::mat C(N,N);

                // Copy forward data
                Xs = Xf;
                Ps = Pf;

                // Backward filter
                for(arma::uword n=Nf-2; n>0; n--)
                {
                    if( !lin_proc )
                    {
                        // Update A with jacobian approx or analytical if set
                        if(f_jac.size()>0)
                            jacobian_analytical(A,f_jac,Xf.col(n));
                        else
                            jacobian_diff(A,f,Xf.col(n));

                        // Project state
                        X_pred = eval_fcn(f,Xf.col(n));
                    }
                    else  // Linear process
                    {
                        // Project state
                        X_pred = A*Xf.col(n);
                    }

                    // Project error covariance
                    P_pred = A*Pf.slice(n)*A.t()+Q;

                    // Update
                    C = Pf.slice(n)*A.t()*inv(P_pred);
                    Xs.col(n)   += C*(Xs.col(n+1)-X_pred);
                    Ps.slice(n) += C*(Ps.slice(n+1)-P_pred)*C.t();
                }

            }

    }; // End class EKF


    ////////////////////////////////////////////////////////////////////////////////////////////
    ///
    /// \brief Unscented Kalman filter class
    ///
    /// Implements Kalman functions for the discrete system with additive noise
    /// \f[ x_k = f(x_{k-1})+Bu_{k-1} + w_{k-1}  \f]
    /// and with measurements
    /// \f[ z_k = h(x_k) + v_k  \f]
    /// where f(x) and h(x) may be nonlinear functions.
    ///
    /// The predict and update stage is using the unscented transform of the sigma points
    /// of the input and/or the measurements. <br>
    /// For detailed info see: http://www.lce.hut.fi/~ssarkka/course_k2011/pdf/course_booklet_2011.pdf
    ////////////////////////////////////////////////////////////////////////////////////////////
    class UKF: public KF
    {
        protected:
            double    alpha;     ///< Spread factor of sigma points
            double    beta;      ///< x distr. prior knowledge factor
            double    kappa;     ///< Scaling par.
            double    lambda;

            arma::mat X;         ///< Sigma points
            arma::mat S;         ///< Output covariance
            arma::mat C;         ///< Cross covariance input-output
            arma::vec Wx;        ///< Weights states
            arma::vec Wp;        ///< Weights covariance
        public:

            UKF(arma::uword _N,arma::uword _M,arma::uword _L): KF(_N,_M,_L)
            {
                alpha  = 1e-3;
                beta   = 2.0;
                kappa  = 0;
                lambda = alpha*alpha*(_N+kappa)-_N;
                X.set_size(_N,2*_N+1);X.zeros();
                Wx.set_size(2*_N+1);Wx.zeros();
                Wp.set_size(2*_N+1);Wp.zeros();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            //  Set/get functions
            ////////////////////////////////////////////////////////////////////////////////////////////
            void set_alpha(double _a)  { alpha   = _a;    }   // Set alpha
            void set_beta(double _b)   { beta    = _b;    }   // Set beta
            void set_kappa(double _k)  { kappa   = _k;    }   // Set kappa
            void set_lambda(double _l) { lambda  = _l;    }   // Set lambda

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Calculate sigma point weights
            ////////////////////////////////////////////////////////////////////////////////////////////
            void update_weights( void )
            {
                // Update lambda
                lambda = alpha*alpha*(N+kappa)-N;

                // Update weights
                Wx(0) = lambda/(N+lambda);
                Wp(0) = lambda/(N+lambda)+(1-alpha*alpha+beta);

                for(arma::uword n=1;n<=2*N;n++)
                {
                    Wx(n) = 1/(2*(N+lambda));
                    Wp(n) = Wx(n);
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Calculate sigma points around a reference point
            /// @param _x State matrix
            /// @param _P Covariance matrix
            ////////////////////////////////////////////////////////////////////////////////////////////
            void update_sigma(const arma::mat& _x, const arma::mat& _P )
            {
                // Update sigma points using Cholesky decomposition
                arma::mat _A = sqrt(N + lambda)*arma::chol(_P,"lower");

                X = arma::repmat(_x,1,2*N+1);
                X.cols(1  ,  N) += _A;
                X.cols(N+1,2*N) -= _A;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Calculate unscented transform
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::mat ut( const arma::mat& _x, const arma::mat& _P, const fcn_v _f )
            {
                arma::uword Ny = static_cast<arma::uword>(_f.size());
                arma::mat y(Ny,1);
                S.set_size(Ny,Ny);
                C.set_size(N,Ny);

                update_weights();
                update_sigma(_x,_P);

                // Propagate sigma points through nonlinear function
                arma::mat Xy(Ny,2*N+1);
                for(arma::uword n=0;n<2*N+1;n++)
                    Xy.col(n) = eval_fcn(_f,X.col(n));

                // New mean
                y = Xy*Wx;

                // New cov
                S = (Xy-arma::repmat(y,1,2*N+1))*arma::diagmat(Wp)*(Xy-arma::repmat(y,1,2*N+1)).t();

                // New crosscov
                C = (X-arma::repmat(_x,1,2*N+1))*arma::diagmat(Wp)*(Xy-arma::repmat(y,1,2*N+1)).t();

                return y;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Predict the internal states using a control input.
            /// @param u Input/control signal
            ////////////////////////////////////////////////////////////////////////////////////////////
            void predict(const arma::mat u )
            {
                if(!lin_proc) // Nonlinear process
                {
                    // Do the Unscented Transform
                    x = ut(x,P,f)+B*u;

                    // Add process noise cov
                    P = S + Q;
                }
                else  // Linear process
                {
                    // Project state
                    x = A*x+B*u;

                    // Project error covariance
                    P = A*P*A.t()+Q;
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Predict the internal states, no control. Convenient function
            ////////////////////////////////////////////////////////////////////////////////////////////
            void predict(void)
            {
                 arma::mat u0(L,1,arma::fill::zeros);
                 predict(u0);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Correct and update the internal states. UKF
            ////////////////////////////////////////////////////////////////////////////////////////////
            void update(const arma::mat z )
            {
                arma::mat z_hat(M,1);
                if(!lin_meas) // Nonlinear measurement
                {
                    // Do the Unscented Transform
                    z_hat = ut(x,P,h);

                    // Add measurement noise cov
                    S = S + R;

                    // Compute the Kalman gain
                    K = C *arma::inv(S);

                    // Update estimate with measurement residual
                    z_err = z-z_hat;
                    x = x+K*z_err;

                    // Update covariance, TODO: improve numerical perf. Josephs form?
                    P = P-K*S*K.t();
                }
                else  // Linear measurement
                {
                    // Calc residual
                    z_err = z-H*x;

                    // Compute the Kalman gain
                    K = P*H.t()*inv(H*P*H.t()+R);

                    // Update estimate with measurement residual
                    x = x+K*z_err;

                    // Joseph’s form covariance update
                    arma::mat Jf = arma::eye<arma::mat>(N,N)-K*H;
                    P = Jf*P*Jf.t()+K*R*K.t();
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Rauch-Tung-Striebel smoother.
            /// See http://www.lce.hut.fi/~ssarkka/course_k2011/pdf/course_booklet_2011.pdf
            ////////////////////////////////////////////////////////////////////////////////////////////
            void rts_smooth(const arma::mat& Xf, const arma::cube& Pf, arma::mat& Xs, arma::cube& Ps )
            {
                arma::uword Nf = Xf.n_cols;
                arma::mat X_pred(N,1,arma::fill::zeros);
                arma::mat P_pred(N,N,arma::fill::zeros);
                arma::mat D_pred(N,N,arma::fill::zeros);

                // Copy forward data
                Xs = Xf;
                Ps = Pf;

                // Backward filter
                for(arma::uword k=Nf-2; k>0; k--)
                {
                    if( !lin_proc )
                    {
                        // Do the unscented transform
                        X_pred = ut(Xf.col(k),Pf.slice(k),f);
                        P_pred = S+Q;

                        // Update
                        D_pred = C*inv(P_pred);
                        Xs.col(k)   += D_pred*(Xs.col(k+1)-X_pred);
                        Ps.slice(k) += D_pred*(Ps.slice(k+1)-P_pred)*D_pred.t();
                    }
                    else  // Linear process
                    {
                        // Project state
                        X_pred = A*Xf.col(k);

                        // Project error covariance
                        P_pred = A*Pf.slice(k)*A.t()+Q;

                        // Update
                        D_pred = Pf.slice(k)*A.t()*inv(P_pred);
                        Xs.col(k)   += D_pred*(Xs.col(k+1)-X_pred);
                        Ps.slice(k) += D_pred*(Ps.slice(k+1)-P_pred)*D_pred.t();
                    }
                }
            }
    }; // End class UKF
}


namespace fftw
{
    ///
    /// @defgroup fftw FFTW
    /// \brief One dimensional FFT functions using FFTW3 library.
    ///
    /// \note If a single FFT is to be used the Armadillo version is faster.
    /// FFTW takes longer time at the first calculation but is faster in the following loops
    /// @{


    ///
    /// \brief FFTW class.
    ///
    /// Implements FFT functions for Armadillo types. For more info see [fftw.org](http://fftw.org/)
    ///
    class FFTW
    {
        private:
            fftw_plan pl_fft;     ///< Real FFTW plan
            fftw_plan pl_ifft;    ///< Real IFFTW plan
            fftw_plan pl_fft_cx;  ///< Complex FFTW plan
            fftw_plan pl_ifft_cx; ///< Complex IFFTW plan
            unsigned int N;       ///< FFT length
            unsigned int R,C;     ///< FFT 2D dims
            int alg;              ///< One of FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE, FFTW_WISDOM_ONLY see [FFTW plans](http://fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags)
            int export_alg;       ///< Alg used for exporting wisdom
        public:
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            /// @param _N FFT length
            /// @param _alg FFTW algorithm selection
            ////////////////////////////////////////////////////////////////////////////////////////////
            FFTW(unsigned int _N, int _alg = FFTW_ESTIMATE)
            {
                N = _N;
                R = 0;
                C = 0;
                alg = _alg;
                export_alg = FFTW_PATIENT;
                pl_fft = NULL;
                pl_ifft = NULL;
                pl_fft_cx = NULL;
                pl_ifft_cx = NULL;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Constructor.
            /// @param _R FFT Nr of rows
            /// @param _C FFT Nr of cols
            /// @param _alg FFTW algorithm selection
            ////////////////////////////////////////////////////////////////////////////////////////////
            FFTW(unsigned int _R, unsigned int _C, int _alg )
            {
                R = _R;
                C = _C;
                N = 0;
                alg = _alg;
                export_alg = FFTW_PATIENT;
                pl_fft = NULL;
                pl_ifft = NULL;
                pl_fft_cx = NULL;
                pl_ifft_cx = NULL;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Destructor.
            ////////////////////////////////////////////////////////////////////////////////////////////
            ~FFTW()
            {
                if (pl_fft != NULL) fftw_destroy_plan(pl_fft);
                if (pl_ifft != NULL) fftw_destroy_plan(pl_ifft);
                if (pl_fft_cx != NULL) fftw_destroy_plan(pl_fft_cx);
                if (pl_ifft_cx != NULL) fftw_destroy_plan(pl_ifft_cx);
                fftw_cleanup();
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief FFT of complex input.
            /// @param x Complex input data
            /// @param[out] Pxx Vector to hold complex FFT of length N
            ////////////////////////////////////////////////////////////////////////////////////////////
            void fft_cx(arma::cx_vec& x, arma::cx_vec& Pxx)
            {
                fftw_complex*  in = reinterpret_cast<fftw_complex*>(x.memptr());
                fftw_complex* out = reinterpret_cast<fftw_complex*>(Pxx.memptr());
                if (pl_fft_cx == NULL)
                {
                    pl_fft_cx = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, alg);
                    if (pl_fft_cx == NULL)
                    {
                        err_handler("Unable to create complex data FFTW plan");
                    }
                }
                fftw_execute_dft(pl_fft_cx, in, out);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief FFT of complex input.
            /// @returns Complex FFT of length N
            /// @param x Complex input data
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_vec fft_cx(arma::cx_vec& x)
            {
                arma::cx_vec Pxx(N);
                fft_cx(x, Pxx);
                return Pxx;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Inverse FFT.
            /// @param Pxx Complex FFT
            /// @param[out] x Vector to hold complex data of length N
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ifft_cx( arma::cx_vec& Pxx, arma::cx_vec& x)
            {
                fftw_complex*  in = reinterpret_cast<fftw_complex*>(Pxx.memptr());
                fftw_complex* out = reinterpret_cast<fftw_complex*>(x.memptr());
                if (pl_ifft_cx == NULL)
                {
                    pl_ifft_cx = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, alg);
                    if (pl_ifft_cx == NULL)
                    {
                        err_handler("Unable to create complex data IFFTW plan");
                    }
                }
                fftw_execute_dft(pl_ifft_cx, in, out);
                x /= N;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Inverse FFT.
            /// @returns Complex data vector of length N
            /// @param Pxx Complex FFT
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::cx_vec ifft_cx( arma::cx_vec& Pxx)
            {
                arma::cx_vec x(N);
                ifft_cx(Pxx, x);
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief FFT of real input.
            /// @param x Input data
            /// @param[out] Pxx Vector to hold complex FFT of length N
            ////////////////////////////////////////////////////////////////////////////////////////////
            void fft( arma::vec& x, arma::cx_vec& Pxx)
            {
                double*        in = x.memptr();
                fftw_complex* out = reinterpret_cast<fftw_complex*>(Pxx.memptr());
                if (pl_fft == NULL)
                {
                    pl_fft = fftw_plan_dft_r2c_1d(N, in, out, alg);
                    if (pl_fft == NULL)
                    {
                        err_handler("Unable to create real data FFTW plan");
                    }
                }

                fftw_execute_dft_r2c(pl_fft, in, out);
                int offset = static_cast<int>(ceil(N / 2.0));
                int n_elem = N - offset;
                for (int i = 0; i < n_elem; ++i)
                {
                    Pxx(offset + i) = std::conj(Pxx(n_elem - i));
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief FFT of real input.
            /// @returns Complex FFT of length N
            /// @param x Real input data
            arma::cx_vec fft( arma::vec& x)
            {
                arma::cx_vec Pxx(N);
                fft(x, Pxx);
                return Pxx;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Inverse FFT.
            /// @param Pxx Complex FFT
            /// @param[out] x Vector to hold real data of length N
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ifft( arma::cx_vec& Pxx, arma::vec& x)
            {
                fftw_complex* in = reinterpret_cast<fftw_complex*>(Pxx.memptr());
                double*      out = x.memptr();
                if (pl_ifft == NULL)
                {
                    pl_ifft = fftw_plan_dft_c2r_1d(N, in, out, alg);
                    if (pl_ifft == NULL)
                    {
                        err_handler("Unable to create real data IFFTW plan");
                    }
                }
                fftw_execute_dft_c2r(pl_ifft, in, out);
                x /= N;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Inverse FFT.
            /// @returns Real data vector of length N
            /// @param Pxx Complex FFT
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::vec ifft( arma::cx_vec& Pxx)
            {
                arma::vec x(N);
                ifft(Pxx, x);
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief FFT of real 2D input.
            /// @param x Input data matrix
            /// @param[out] Pxx Matrix to hold complex FFT of length [RxC]
            ////////////////////////////////////////////////////////////////////////////////////////////
            void fft2( arma::mat& x,  arma::cx_mat& Pxx)
            {
                arma::cx_mat Ptmp(R / 2 + 1, C,arma::fill::ones);
                double*        in = x.memptr();
                fftw_complex* out = reinterpret_cast<fftw_complex*>(Ptmp.memptr());
                if (pl_fft == NULL)
                {
                    pl_fft = fftw_plan_dft_r2c_2d( C,R, in, out, alg); // Column to row-major order trick: switch C and R
                    if (pl_fft == NULL)
                    {
                        err_handler("Unable to create real data FFTW plan");
                    }
                }

                fftw_execute_dft_r2c(pl_fft, in, out);

                // Reshape Pxx - upper half
                std::complex<double>* ptr= reinterpret_cast<std::complex<double>*>(out);
                const unsigned int Roff = R / 2 + 1;
                for (unsigned int r = 0; r < Roff; r++)
                {
                    for (unsigned int c = 0; c < C; c++)
                    {
                        Pxx(r, c) = ptr[r + c*Roff];
                    }
                }
                // Reshape Pxx - conj symmetry
                for (unsigned int r = Roff; r < R; r++)
                {
                    Pxx(r, 0) = conj(ptr[R - r]);
                    for (unsigned int c = 1; c < C; c++)
                    {
                        Pxx(r, c) = conj(ptr[R - r + (C - c)*Roff]);
                    }
                }
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief FFT of real 2D input.
            /// @returns Complex FFT of size [RxC]
            /// @param x Real input matrix
            arma::cx_mat fft2( arma::mat& x)
            {
                arma::cx_mat Pxx(R,C,arma::fill::ones);
                fft2(x, Pxx);
                return Pxx;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Inverse 2D FFT.
            /// @param Pxx Complex FFT
            /// @param[out] x Matrix to hold real data of size[RxC]
            ////////////////////////////////////////////////////////////////////////////////////////////
            void ifft2( arma::cx_mat& Pxx, arma::mat& x)
            {
                // Reshape to row-major format
                unsigned int Roff = R / 2 + 1;
                arma::cx_mat Ptmp(Roff, C);
                std::complex<double>* ptr = reinterpret_cast<std::complex<double>*>(Ptmp.memptr());
                for (unsigned int r = 0; r < Roff; r++)
                {
                    for (unsigned int c = 0; c < C; c++)
                    {
                        ptr[r + c*Roff] = Pxx(r,c);
                    }
                }

                fftw_complex* in = reinterpret_cast<fftw_complex*>(Ptmp.memptr());
                double*      out = x.memptr();
                if (pl_ifft == NULL)
                {
                    pl_ifft = fftw_plan_dft_c2r_2d(C,R, in, out, alg);
                    if (pl_ifft == NULL)
                    {
                        err_handler("Unable to create real data IFFTW plan");
                    }
                }
                fftw_execute_dft_c2r(pl_ifft, in, out);
                x /= (R*C);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Inverse FFT.
            /// @returns Real data vector of length N
            /// @param Pxx Complex FFT
            ////////////////////////////////////////////////////////////////////////////////////////////
            arma::mat ifft2( arma::cx_mat& Pxx)
            {
                arma::mat x(R,C);
                ifft2(Pxx, x);
                return x;
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Import wisdom from string.
            /// @param wisd Wisdom string
            ////////////////////////////////////////////////////////////////////////////////////////////
            void import_wisdom_string(const std::string wisd)
            {
                int res = fftw_import_wisdom_from_string(wisd.c_str());
                if (res == 0)
                    err_handler("Unable to import wisdom from string!");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Import wisdom from file.
            /// @param fname File name
            ////////////////////////////////////////////////////////////////////////////////////////////
            void import_wisdom_file(const std::string fname)
            {
                int res = fftw_import_wisdom_from_filename(fname.c_str());
                if (res == 0)
                    err_handler("Unable to import wisdom from file!");
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Export real FFT wisdom to file.
            /// @param fname File name
            ////////////////////////////////////////////////////////////////////////////////////////////
            void export_wisdom_fft(const std::string fname)
            {
                fftw_plan pl_w = NULL;
                double* x_r;
                fftw_complex* x_cx1;

                if(R==0 || C==0)   // 1D
                {
                    x_r   = fftw_alloc_real(N);
                    x_cx1 = fftw_alloc_complex(N);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_r2c_1d(N, x_r, x_cx1, export_alg);

                }
                else             // 2D
                {
                    x_r   = fftw_alloc_real(R*C);
                    x_cx1 = fftw_alloc_complex(R*C);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_r2c_2d(C,R, x_r, x_cx1, export_alg);
                }
                if (pl_w == NULL)
                {
                    err_handler("Unable to create real data FFTW plan");
                }

                // Export
                if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
                {
                    err_handler("Could not export wisdom to file!");
                }

                fftw_destroy_plan(pl_w);
                fftw_free(x_r);
                fftw_free(x_cx1);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Export real IFFT wisdom to file.
            /// @param fname File name
            ////////////////////////////////////////////////////////////////////////////////////////////
            void export_wisdom_ifft(const std::string fname)
            {
                fftw_plan pl_w = NULL;
                double* x_r;
                fftw_complex* x_cx1;

                if(R==0 || C==0)   // 1D
                {
                    x_r   = fftw_alloc_real(N);
                    x_cx1 = fftw_alloc_complex(N);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_c2r_1d(N, x_cx1, x_r, export_alg);
                }
                else             // 2D
                {
                    x_r   = fftw_alloc_real(R*C);
                    x_cx1 = fftw_alloc_complex(R*C);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_c2r_2d(C,R, x_cx1, x_r, export_alg);
                }

                if (pl_w == NULL)
                {
                    err_handler("Unable to create real data FFTW plan");
                }

                // Export
                if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
                {
                    err_handler("Could not export wisdom to file!");
                }

                fftw_destroy_plan(pl_w);
                fftw_free(x_r);
                fftw_free(x_cx1);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Export complex FFT wisdom to file.
            /// @param fname File name
            ////////////////////////////////////////////////////////////////////////////////////////////
            void export_wisdom_fft_cx(const std::string fname)
            {
                fftw_plan pl_w = NULL;
                fftw_complex* x_cx1, *x_cx2;

                if(R==0 || C==0)      // 1D
                {
                    x_cx1 = fftw_alloc_complex(N);
                    x_cx2 = fftw_alloc_complex(N);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_1d(N, x_cx1, x_cx2, FFTW_FORWARD, export_alg);
                }
                else
                {
                    x_cx1 = fftw_alloc_complex(R*C);
                    x_cx2 = fftw_alloc_complex(R*C);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_2d(C, R, x_cx1, x_cx2, FFTW_FORWARD, export_alg);
                }

                if (pl_w == NULL)
                {
                    err_handler("Unable to create complex data FFTW plan");
                }

                // Export
                if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
                {
                    err_handler("Could not export wisdom to file!");
                }

                fftw_destroy_plan(pl_w);
                fftw_free(x_cx1);
                fftw_free(x_cx2);
            }
            ////////////////////////////////////////////////////////////////////////////////////////////
            /// \brief Export complex IFFT wisdom to file.
            /// @param fname File name
            ////////////////////////////////////////////////////////////////////////////////////////////
            void export_wisdom_ifft_cx(const std::string fname)
            {
                fftw_plan pl_w = NULL;
                fftw_complex* x_cx1, *x_cx2;

                if(R==0 || C==0)      // 1D
                {
                    x_cx1 = fftw_alloc_complex(N);
                    x_cx2 = fftw_alloc_complex(N);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_1d(N, x_cx2, x_cx1, FFTW_BACKWARD, export_alg);
                }
                else
                {
                    x_cx1 = fftw_alloc_complex(R*C);
                    x_cx2 = fftw_alloc_complex(R*C);

                    // Replan using wisdom
                    pl_w = fftw_plan_dft_2d(C, R, x_cx2, x_cx1, FFTW_BACKWARD, export_alg);
                }
                if (pl_w == NULL)
                {
                    err_handler("Unable to create complex data IFFTW plan");
                }

                // Export
                if (fftw_export_wisdom_to_filename(fname.c_str()) == 0)
                {
                    err_handler("Could not export wisdom to file!");
                }

                fftw_destroy_plan(pl_w);
                fftw_free(x_cx1);
                fftw_free(x_cx2);
            }
    };
} // end namespace