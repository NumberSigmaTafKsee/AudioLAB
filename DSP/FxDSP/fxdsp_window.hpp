#pragma once

#include <cstdlib>
#include <cstddef>

namespace FXDSP
{

    /* Modified Bessel function of the first kind */
    template<typename T>
    T modZeroBessel(T x)
    {
        T x_2 = x/2;
        T num = 1;
        T fact = 1;
        T result = 1;
        
        unsigned i;
        #pragma omp simd
        for (i=1 ; i<20 ; i++) 
        {
            num *= x_2 * x_2;
            fact *= i;
            result += num / (fact * fact);
        }
        return result;
    }


    
    /* Chebyshev Polynomial */

    template<typename T>
    T chebyshev_poly(int n, T x)
    {
        T y;
        if (std::fabs(x) <= 1)
        {
            y = std::cos(n * std::acos(x));
        }
        else
        {
            y = std::cosh(n * std::acosh(x));
        }
        
        return y;
    }

    

    /*******************************************************************************
    Hamming */
    template<typename T>
    void hamming(unsigned n, T* dest)
    {
        // do it manually. it seems like vDSP_hamm_window is wrong.
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            *dest++ = 0.54 - 0.46 * std::cos((2 * M_PI * buf_idx) / (n - 1));
        }        
    }

    template<typename T>
    void blackman(unsigned n, T a, T* dest)
    {        
        T a0 = (1 - a) / 2;
        T a1 = 0.5;
        T a2 = a / 2;
        
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            *dest++ = a0 - a1 * std::cos((2 * M_PI * buf_idx) / (n - 1)) + a2 * cosf((4 * M_PI * buf_idx) / (n - 1));
        }        
    }


    template<typename T>
    void tukey(unsigned n, T a, T* dest)
    {
        T term = a * (n - 1) / 2;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            if (buf_idx <= term)
            {
                dest[buf_idx] = 0.5 * (1 + std::cos(M_PI * ((2 * buf_idx) /
                                                        (a * (n - 1)) - 1)));
            }
            else if (term <= buf_idx && buf_idx <= (n - 1) * (1 - a / 2))
            {
                dest[buf_idx] = 1.0;
            }
            else
            {
                dest[buf_idx] = 0.5 * (1 + std::cos(M_PI *((2 * buf_idx) /
                                                    (a * (n - 1)) - (2 / a) + 1)));
            }                    
        }        
    }

    template<typename T>
    void cosine(unsigned n, T* dest)
    {
        T N = n - 1.0;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            *dest++ = std::sin((M_PI * buf_idx) / N);
        }      
    }

    template<typename T>
    void lanczos(unsigned n, T* dest)
    {
        T N = n - 1.0;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            T term = M_PI * ((2 * n)/ N) - 1.0;
            *dest++ = std::sin(term) / term;
        }
        return NOERR;
    }

    template<typename T>
    void bartlett(unsigned n, T* dest)
    {
        unsigned buf_idx;        
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            if (buf_idx <= (n - 1) / 2)
            {
                *dest++ = (T)(2 * buf_idx)/(n - 1);
            }
            else
            {
                *dest ++ = 2.0 - (T)(2 * buf_idx) / (n - 1);
            }
        }    
    }

    template<typename T>
    void gaussian(unsigned n, T sigma, T* dest)
    {
        T N = n - 1;
        T L = N / 2.0;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {   
            *dest++ = std::exp(-0.5 * std::pow((buf_idx - L)/(sigma * L),2));
        }        
    }

    template<typename T>
    void bartlett_hann(unsigned n, T* dest)
    {
        T N = n - 1.0;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            T term = ((buf_idx / N) - 0.5);
            *dest++ = 0.62 - 0.48 * std::fab(term) + 0.38 * std::cos(2 * M_PI * term);
        }      
    }

    template<typename T>
    void kaiser(unsigned n, T a, T* dest)
    {
        // Pre-calc
        T beta = M_PI * a;
        T m_2 = (T)(n - 1.0) / 2.0;
        T denom = modZeroBessel(beta);
        
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            T val = ((buf_idx) - m_2) / m_2;
            val = 1 - (val * val);
            *dest++ = modZeroBessel(beta * std::sqrt(val)) / denom;
        }      
    }

    template<typename T>
    void nuttall(unsigned n, T* dest)
    {
        T term;
        T N = n - 1.0;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            term = 2 * M_PI * (buf_idx / N);
            *dest++ = 0.3635819  - 0.4891775 * std::cos(term) + 0.1365995 *
            std::cos(2 * term) - 0.0106411 * std::cos(3 * term);
        }        
    }

    template<typename T>
    void blackman_harris(unsigned n, T* dest)
    {
        T term;
        unsigned buf_idx
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            term = (2 * M_PI * buf_idx) / (n - 1);
            *dest++ = 0.35875 - 0.48829 * std::cos(term)+ 0.14128 * std::cos(2 * term) -
            0.01168 * std::cos(3 * term);
        }        
    }

    template<typename T>
    blackman_nuttall(unsigned n, T* dest)
    {
        T term;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            term = (2 * M_PI * buf_idx) / (n - 1);
            *dest++ = 0.3635819 - 0.4891775 * cosf(term)+ 0.1365995 * std::cos(2 * term) - 0.0106411 * std::cos(3 * term);
        }
    }
    template<typename T>
    void flat_top(unsigned n, T* dest)
    {
        T N = n - 1.0;
        T term;
        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            term = (2 * M_PI * buf_idx) / N;
            *dest++ = 0.21557895 - 0.41663158 * std::cos(term)+ 0.277263158 *
            std::cos(2 * term) - 0.083578947 * std::cos(3 * term) + 0.006947368 *
            std::cos(4 * term);
        }        
    }
    template<typename T>
    void poisson(unsigned n, T D, T* dest)
    {
        T term = (n - 1) / 2;
        T tau_inv = 1. / ((n / 2) * (8.69 / D));

        unsigned buf_idx;
        #pragma omp simd
        for (buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            *dest++ = std::exp(-std::fab(buf_idx - term) * tau_inv);
        }
        return NOERR;    
    }
    template<typename T>
    void chebyshev(unsigned n, T A, T *dest)
    {
        T max = 0;
        T N = n - 1.0;
        T M = N / 2;
        T tg = std::pow(10, A / 20.0);
        T x0 = std::cosh((1.0 / N) * std::acosh(tg));
        
        #pragma omp simd
        for(unsigned buf_idx=0; buf_idx<(n/2+1); ++buf_idx)
        {
            T y = buf_idx - M;
            T sum = 0;
            for(unsigned i=1; i<=M; i++){
                sum += chebyshev_poly(N, x0 * std::cos(M_PI * i / n)) *
                std::cos( 2.0 * y * M_PI * i / n);
            }
            dest[buf_idx] = tg + 2 * sum;
            dest[(unsigned)N - buf_idx] = dest[buf_idx];
            if(dest[buf_idx] > max)
            {
                max = dest[buf_idx];
            }
        }
        #pragma omp simd
        for(unsigned buf_idx = 0; buf_idx < n; ++buf_idx)
        {
            dest[buf_idx] /= max;
        }
    }

    template<typename T>
    void boxcar(unsigned n, T* dest)
    {
        FillBuffer(dest, n, 1.0);        
    }

    /** Window function type */
    typedef enum _Window_t
    {
        /** Rectangular window */
        BOXCAR,

        /** Hann window */
        HANN,

        /** Hamming window */
        HAMMING,

        /** Blackman window */
        BLACKMAN,

        /** Tukey window */
        TUKEY,

        /** Cosine window */
        COSINE,

        /** Lanczos window */
        LANCZOS,

        /** Bartlett window */
        BARTLETT,

        /** Gauss window */
        GAUSSIAN,

        /** Bartlett-Hann window */
        BARTLETT_HANN,

        /** Kaiser window */
        KAISER,

        /** Nuttall window */
        NUTTALL,

        /** Blackaman-Harris window */
        BLACKMAN_HARRIS,

        /** Blackman-Nuttall window */
        BLACKMAN_NUTTALL,

        /** Flat top window */
        FLATTOP,

        /** Poisson window */
        POISSON,

        /** The number of window types */
        N_WINDOWTYPES
    } Window_t;

    template<typename T>
    struct WindowFunction
    {
        T*      window;
        unsigned    length;
        Window_t    type;

        WindowFunction(size_t n, Window_t type)
        {
            length = n;
            window = new T[n]; 
            this->type = type;
            
            switch (type)
            {
                case BOXCAR:
                    boxcar(length, window);
                    break;
                case HANN:
                    hann(length, window);
                    break;
                case HAMMING:
                    hamming(length, window);
                    break;
                case BLACKMAN:
                    blackman(length, 0.16, window);
                    break;
                case TUKEY:
                    tukey(length, 0.5, window);
                    break;
                case COSINE:
                    cosine(length, window);
                    break;
                case LANCZOS:
                    lanczos(length, window);
                    break;
                case BARTLETT:
                    bartlett(length, window);
                    break;
                case GAUSSIAN:
                    gaussian(length, 0.4, window);
                    break;
                case BARTLETT_HANN:
                    bartlett_hann(length, window);
                    break;
                case KAISER:
                    kaiser(length, 0.5, window);
                    break;
                case NUTTALL:
                    nuttall(length, window);
                    break;
                case BLACKMAN_HARRIS:
                    blackman_harris(length, window);
                    break;
                case BLACKMAN_NUTTALL:
                    blackman_nuttall(length, window);
                    break;
                case FLATTOP:
                    flat_top(length, window);
                    break;
                case POISSON:
                    poisson(length, 8.69, window);
                    
                default:
                    boxcar(length, window);
                    break;
            }
        }
        ~WindowFunction() {
            if(window) delete [] window;
        }
        void ProcessBlock(size_t n_samples, T * inBuffer, T * outBuffer)
        {
            VectorVectorMultiply(outBuffer, inBuffer, window, n_samples);
        }
    };
}