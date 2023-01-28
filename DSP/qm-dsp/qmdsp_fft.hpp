#pragma once

#include "qmdsp_math.hpp"
#include "kiss_fft.h"
#include "kiss_fftr.h"
#include <vector>
#include <cmath>

#include <stdexcept>


namespace qmdsp
{
    class FFT  
    {
    public:
        /**
        * Construct an FFT object to carry out complex-to-complex
        * transforms of size nsamples. nsamples does not have to be a
        * power of two.
        */
        FFT(int nsamples);
        ~FFT();

        /**
        * Carry out a forward or inverse transform (depending on the
        * value of inverse) of size nsamples, where nsamples is the value
        * provided to the constructor above.
        *
        * realIn and (where present) imagIn should contain nsamples each,
        * and realOut and imagOut should point to enough space to receive
        * nsamples each.
        *
        * imagIn may be NULL if the signal is real, but the other
        * pointers must be valid.
        *
        * The inverse transform is scaled by 1/nsamples.
        */
        void process(bool inverse,
                    const double *realIn, const double *imagIn,
                    double *realOut, double *imagOut);
        
    private:
        class D;
        D *m_d;
    };

    class FFTReal
    {
    public:
        /**
        * Construct an FFT object to carry out real-to-complex transforms
        * of size nsamples. nsamples does not have to be a power of two,
        * but it does have to be even. (Use the complex-complex FFT above
        * if you need an odd FFT size. This constructor will throw
        * std::invalid_argument if nsamples is odd.)
        */
        FFTReal(int nsamples);
        ~FFTReal();

        /**
        * Carry out a forward real-to-complex transform of size nsamples,
        * where nsamples is the value provided to the constructor above.
        *
        * realIn, realOut, and imagOut must point to (enough space for)
        * nsamples values. For consistency with the FFT class above, and
        * compatibility with existing code, the conjugate half of the
        * output is returned even though it is redundant.
        */
        void forward(const double *realIn,
                    double *realOut, double *imagOut);

        /**
        * Carry out a forward real-to-complex transform of size nsamples,
        * where nsamples is the value provided to the constructor
        * above. Return only the magnitudes of the complex output values.
        *
        * realIn and magOut must point to (enough space for) nsamples
        * values. For consistency with the FFT class above, and
        * compatibility with existing code, the conjugate half of the
        * output is returned even though it is redundant.
        */
        void forwardMagnitude(const double *realIn, double *magOut);

        /**
        * Carry out an inverse real transform (i.e. complex-to-real) of
        * size nsamples, where nsamples is the value provided to the
        * constructor above.
        *
        * realIn and imagIn should point to at least nsamples/2+1 values;
        * if more are provided, only the first nsamples/2+1 values of
        * each will be used (the conjugate half will always be deduced
        * from the first nsamples/2+1 rather than being read from the
        * input data).  realOut should point to enough space to receive
        * nsamples values.
        *
        * The inverse transform is scaled by 1/nsamples.
        */
        void inverse(const double *realIn, const double *imagIn,
                    double *realOut);

    private:
        class D;
        D *m_d;
    };    

    class DCT
    {
    public:
        /**
        * Construct a DCT object to calculate the Discrete Cosine
        * Transform given input of size n samples. The transform is
        * implemented using an FFT of size 4n, for simplicity.
        */
        DCT(int n);

        ~DCT();

        /**
        * Carry out a type-II DCT of size n, where n is the value
        * provided to the constructor above.
        *
        * The in and out pointers must point to (enough space for) n
        * values.
        */
        void forward(const double *in, double *out);

        /**
        * Carry out a type-II unitary DCT of size n, where n is the value
        * provided to the constructor above. This is a scaled version of
        * the type-II DCT corresponding to a transform with an orthogonal
        * matrix. This is the transform implemented by the dct() function
        * in MATLAB.
        *
        * The in and out pointers must point to (enough space for) n
        * values.
        */
        void forwardUnitary(const double *in, double *out);

        /**
        * Carry out a type-III (inverse) DCT of size n, where n is the
        * value provided to the constructor above.
        *
        * The in and out pointers must point to (enough space for) n
        * values.
        */
        void inverse(const double *in, double *out);

        /**
        * Carry out a type-III (inverse) unitary DCT of size n, where n
        * is the value provided to the constructor above. This is the
        * inverse of forwardUnitary().
        *
        * The in and out pointers must point to (enough space for) n
        * values.
        */
        void inverseUnitary(const double *in, double *out);

    private:
        int m_n;
        double m_scale;
        std::vector<double> m_scaled;
        std::vector<double> m_time2;
        std::vector<double> m_freq2r;
        std::vector<double> m_freq2i;
        FFTReal m_fft;
    };

    DCT::DCT(int n) :
        m_n(n),
        m_scaled(n, 0.0),
        m_time2(n*4, 0.0),
        m_freq2r(n*4, 0.0),
        m_freq2i(n*4, 0.0),
        m_fft(n*4)
    {
        m_scale = m_n * sqrt(2.0 / m_n);
    }

    DCT::~DCT()
    {
    }

    void
    DCT::forward(const double *in, double *out)
    {
        for (int i = 0; i < m_n; ++i) {
            m_time2[i*2 + 1] = in[i];
            m_time2[m_n*4 - i*2 - 1] = in[i];
        }

        m_fft.forward(m_time2.data(), m_freq2r.data(), m_freq2i.data());

        for (int i = 0; i < m_n; ++i) {
            out[i] = m_freq2r[i];
        }
    }

    void
    DCT::forwardUnitary(const double *in, double *out)
    {
        forward(in, out);
        for (int i = 0; i < m_n; ++i) {
            out[i] /= m_scale;
        }
        out[0] /= sqrt(2.0);
    }

    void
    DCT::inverse(const double *in, double *out)
    {
        for (int i = 0; i < m_n; ++i) {
            m_freq2r[i] = in[i];
        }
        for (int i = 0; i < m_n; ++i) {
            m_freq2r[m_n*2 - i] = -in[i];
        }
        m_freq2r[m_n] = 0.0;

        for (int i = 0; i <= m_n*2; ++i) {
            m_freq2i[i] = 0.0;
        }
        
        m_fft.inverse(m_freq2r.data(), m_freq2i.data(), m_time2.data());

        for (int i = 0; i < m_n; ++i) {
            out[i] = m_time2[i*2 + 1];
        }
    }

    void
    DCT::inverseUnitary(const double *in, double *out)
    {
        for (int i = 0; i < m_n; ++i) {
            m_scaled[i] = in[i] * m_scale;
        }
        m_scaled[0] *= sqrt(2.0);
        inverse(m_scaled.data(), out);
    }

    class FFT::D
    {
    public:
        D(int n) : m_n(n) {
            m_planf = kiss_fft_alloc(m_n, 0, NULL, NULL);
            m_plani = kiss_fft_alloc(m_n, 1, NULL, NULL);
            m_kin = new kiss_fft_cpx[m_n];
            m_kout = new kiss_fft_cpx[m_n];
        }

        ~D() {
            kiss_fft_free(m_planf);
            kiss_fft_free(m_plani);
            delete[] m_kin;
            delete[] m_kout;
        }

        void process(bool inverse,
                    const double *ri,
                    const double *ii,
                    double *ro,
                    double *io) {

            for (int i = 0; i < m_n; ++i) {
                m_kin[i].r = ri[i];
                m_kin[i].i = (ii ? ii[i] : 0.0);
            }

            if (!inverse) {

                kiss_fft(m_planf, m_kin, m_kout);

                for (int i = 0; i < m_n; ++i) {
                    ro[i] = m_kout[i].r;
                    io[i] = m_kout[i].i;
                }

            } else {

                kiss_fft(m_plani, m_kin, m_kout);

                double scale = 1.0 / m_n;

                for (int i = 0; i < m_n; ++i) {
                    ro[i] = m_kout[i].r * scale;
                    io[i] = m_kout[i].i * scale;
                }
            }
        }
        
    private:
        int m_n;
        kiss_fft_cfg m_planf;
        kiss_fft_cfg m_plani;
        kiss_fft_cpx *m_kin;
        kiss_fft_cpx *m_kout;
    };        

    FFT::FFT(int n) :
        m_d(new D(n))
    {
    }

    FFT::~FFT()
    {
        delete m_d;
    }

    void
    FFT::process(bool inverse,
                const double *p_lpRealIn, const double *p_lpImagIn,
                double *p_lpRealOut, double *p_lpImagOut)
    {
        m_d->process(inverse,
                    p_lpRealIn, p_lpImagIn,
                    p_lpRealOut, p_lpImagOut);
    }
        
    class FFTReal::D
    {
    public:
        D(int n) : m_n(n) {
            if (n % 2) {
                throw std::invalid_argument
                    ("nsamples must be even in FFTReal constructor");
            }
            m_planf = kiss_fftr_alloc(m_n, 0, NULL, NULL);
            m_plani = kiss_fftr_alloc(m_n, 1, NULL, NULL);
            m_c = new kiss_fft_cpx[m_n];
        }

        ~D() {
            kiss_fftr_free(m_planf);
            kiss_fftr_free(m_plani);
            delete[] m_c;
        }

        void forward(const double *ri, double *ro, double *io) {

            kiss_fftr(m_planf, ri, m_c);

            for (int i = 0; i <= m_n/2; ++i) {
                ro[i] = m_c[i].r;
                io[i] = m_c[i].i;
            }

            for (int i = 0; i + 1 < m_n/2; ++i) {
                ro[m_n - i - 1] =  ro[i + 1];
                io[m_n - i - 1] = -io[i + 1];
            }
        }

        void forwardMagnitude(const double *ri, double *mo) {

            double *io = new double[m_n];

            forward(ri, mo, io);

            for (int i = 0; i < m_n; ++i) {
                mo[i] = sqrt(mo[i] * mo[i] + io[i] * io[i]);
            }

            delete[] io;
        }

        void inverse(const double *ri, const double *ii, double *ro) {

            // kiss_fftr.h says
            // "input freqdata has nfft/2+1 complex points"

            for (int i = 0; i < m_n/2 + 1; ++i) {
                m_c[i].r = ri[i];
                m_c[i].i = ii[i];
            }
            
            kiss_fftri(m_plani, m_c, ro);

            double scale = 1.0 / m_n;

            for (int i = 0; i < m_n; ++i) {
                ro[i] *= scale;
            }
        }

    private:
        int m_n;
        kiss_fftr_cfg m_planf;
        kiss_fftr_cfg m_plani;
        kiss_fft_cpx *m_c;
    };

    FFTReal::FFTReal(int n) :
        m_d(new D(n)) 
    {
    }

    FFTReal::~FFTReal()
    {
        delete m_d;
    }

    void
    FFTReal::forward(const double *ri, double *ro, double *io)
    {
        m_d->forward(ri, ro, io);
    }

    void
    FFTReal::forwardMagnitude(const double *ri, double *mo)
    {
        m_d->forwardMagnitude(ri, mo);
    }

    void
    FFTReal::inverse(const double *ri, const double *ii, double *ro)
    {
        m_d->inverse(ri, ii, ro);
    }

}