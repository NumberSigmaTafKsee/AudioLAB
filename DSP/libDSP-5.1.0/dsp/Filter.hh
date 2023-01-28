/*

    FFT filter
    Copyright (C) 1999-2005 Jussi Laako

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/


#include <dsp/DSPConfig.hh>
#include <dsp/DSPOp.hh>
#include <dsp/ReBuffer.hh>


#ifndef FILTER_HH
    #define FILTER_HH


    /**
        Filter smoothing window types.
    */
    enum
    {
        FILTER_SMOOTH_NONE = 0,
        FILTER_SMOOTH_KAISER = 1,
        FILTER_SMOOTH_KAISER_BESSEL = 2,
        FILTER_SMOOTH_DOLPH_CHEBYSHEV = 3
    };

    /**
        Class implementing FFT-based FIR-filter with design functions.
        Filtering is done in complex-plane.
    */
    class clFilter : private clDSPOp
    {
            bool bInitialized;
            long lFFTSize;
            long lHalfSize;
            long lOldSize;
            long lNewSize;
            long lSpectPoints;
            clDSPAlloc CoeffWin;
            clDSPAlloc Prev;
            clDSPAlloc Proc;
            clDSPAlloc CCoeffs;
            clDSPAlloc CProc;
            clReBuffer InBuf;
            clReBuffer OutBuf;
            /**
                Initialize complex coefficients.

                \f[\left\vert\begin{array}{l}\Re_{C}(i)=1\\
                \Im_{C}(i)=0\end{array}\right.\f]
            */
            void InitCoeffsS ();
            /// \overload
            void InitCoeffsD ();
            /**
                Create complex impulse response from specified frequency
                response.

                This is done by doing \f$FFT^{-1}(X)\f$, then multiplying
                it with desired window function and doing \f$FFT(x)\f$.
            */
            void ReadyFilterS ();
            /// \overload
            void ReadyFilterD ();
            /**
                Calculate beta value for Kaiser window.
                
                \f[\beta=\left{\begin{array}{ll}0.1102(\alpha_{s}-8.7)&, \alpha_{s}>50\\
                0.5842(\alpha_{s}-21)^{0.4}+0.07886(\alpha_{s}-21&, 21\leq\alpha_{s}\leq{50}\\
                0&, \alpha_{s}<21\end{array}\right.\f]
                
                \param Alpha Stopband rippleratio in dB
                \return Beta for Kaiser window
            */
            float GetKaiserBeta (float fAlpha);
            /// \overload
            double GetKaiserBeta (double dAlpha);
        public:
            clFilter ();
            virtual ~clFilter ();
            /**
                Initialize filter.

                \param WindowSize Number of FFT points used in filter
                \param FiltCoeffs Optional coefficient vector, can be NULL
                \param Overlap Normalized overlap factor
                \param Beta Optional beta value for window
                \param SmoothWindow Type of smoothing window to use
                \return Success
            */
            bool Initialize (long, const float *,
                float = DSP_FILT_DEF_OVERLAPF,
                float = DSP_FILT_DEF_BETAF,
                int = FILTER_SMOOTH_DOLPH_CHEBYSHEV);
            /// \overload
            bool Initialize (long, const double *,
                double = DSP_FILT_DEF_OVERLAP,
                double = DSP_FILT_DEF_BETA,
                int = FILTER_SMOOTH_DOLPH_CHEBYSHEV);
            /// \overload
            bool InitFloat (long);
            /// \overload
            bool InitDouble (long);
            /**
                Initialize low-pass or high-pass filter.
                
                \f[N_{o}=\frac{\alpha_{s}-8}{2.285\delta\omega}\f]
                
                \param PassBand Normalized passband
                \param StopBand Normalized stopband
                \param RippleRatio Stopband ripple ratio in dB
                \param Overlap Normalized overlap factor
                \return Success
            */
            bool InitializeLP (float, float, float,
                float = DSP_FILT_DEF_OVERLAPF);
            /// \overload
            bool InitializeLP (double, double, double,
                double = DSP_FILT_DEF_OVERLAP);
            /// \overload
            bool InitializeHP (float, float, float,
                float = DSP_FILT_DEF_OVERLAPF);
            /// \overload
            bool InitializeHP (double, double, double,
                double = DSP_FILT_DEF_OVERLAP);
            /**
                Uninitialize filter.
            */
            void Uninitialize ();
            /**
                Set coefficients.

                You can feed the filter with desired frequency response.
                Length of this dataset is specified window size / 2 + 1.
                Filter is then created from this coefficient vector.

                \param FiltCoeffs Coefficients for designing filter
                \param Smooth Smooth filter response using smoothing window
            */
            void SetCoeffs (const float *);
            /// \overload
            void SetCoeffs (const double *);
            /// \overload
            void SetCoeffs (const stpSCplx, bool = false);
            /// \overload
            void SetCoeffs (const stpDCplx, bool = false);
            /**
                Get coefficients.

                Get filter's frequency response.

                \param FiltCoeffs Filter coefficients
            */
            void GetCoeffs (float *);
            /// \overload
            void GetCoeffs (double *);
            /// \overload
            void GetCoeffs (stpSCplx);
            /// \overload
            void GetCoeffs (stpDCplx);
            /**
                Get filter delay in samples.
                
                Delay caused by the filter in number of samples.
                
                \return Delay (in samples)
            */
            long GetDelay ();
            /**
                Feed data into filter.

                Filtering is done by first doing \f$FFT(x)\f$, then doing
                complex multiply with complex impulse response, and finally
                doing \f$FFT^{-1}(X)\f$. Windows overlap with specified
                overlap factor and 1 - overlap factor of window length is used.

                \param SrcData Source data
                \param SrcCount Source data count
            */
            void Put (const float *, long);
            /// \overload
            void Put (const double *, long);
            /**
                Feed data into filter using complex coefficients.

                \param SrcData Source data
                \param SrcCount Source data count
                \param Coeffs Complex coefficients
            */
            void Put (const float *, long, const stpSCplx);
            /// \overload
            void Put (const double *, long, const stpDCplx);
            /**
                Get data from filter.

                \param DestData Destination buffer
                \param DestCount Number of samples
                \return Data was available?
            */
            bool Get (float *, long);
            /// \overload
            bool Get (double *, long);
            /**
                Design low-pass filter.

                \param Corner Corner frequency
                \param SampleRate Sampling frequency
                \param DCBlock Include DC-blocking
            */
            void DesignLP (float *, bool = false);
            /// \overload
            void DesignLP (double *, bool = false);
            /// \overload
            void DesignLP (float *, float, bool = false);
            /// \overload
            void DesignLP (double *, double, bool = false);
            /**
                Design high-pass filter.
            */
            void DesignHP (float *);
            /// \overload
            void DesignHP (double *);
            /// \overload
            void DesignHP (float *, float);
            /// \overload
            void DesignHP (double *, double);
            /**
                Design band-pass filter.

                \param LowCorner Lower corner frequency
                \param HighCorner Higher corner frequency
                \param SampleRate Sampling frequency
            */
            void DesignBP (float *, float *);
            /// \overload
            void DesignBP (double *, double *);
            /// \overload
            void DesignBP (float *, float *, float);
            /// \overload
            void DesignBP (double *, double *, double);
            /**
                Design band-reject filter.
            */
            void DesignBR (float *, float *);
            /// \overload
            void DesignBR (double *, double *);
            /// \overload
            void DesignBR (float *, float *, float);
            /// \overload
            void DesignBR (double *, double *, double);
    };

#endif
