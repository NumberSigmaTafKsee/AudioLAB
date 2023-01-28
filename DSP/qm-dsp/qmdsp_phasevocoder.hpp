#pragma once

#include "qmdsp_math.hpp"
#include "qmdsp_fft.hpp"

namespace qmdsp
{
    class PhaseVocoder  
    {
    public:
        PhaseVocoder(int size, int hop);
        virtual ~PhaseVocoder();

        /**
        * Given one frame of time-domain samples, FFT and return the
        * magnitudes, instantaneous phases, and unwrapped phases.
        *
        * src must have size values (where size is the frame size value
        * as passed to the PhaseVocoder constructor), and should have
        * been windowed as necessary by the caller (but not fft-shifted).
        *
        * mag, phase, and unwrapped must each be non-NULL and point to
        * enough space for size/2 + 1 values. The redundant conjugate
        * half of the output is not returned.
        */
        void processTimeDomain(const double *src,
                            double *mag, double *phase, double *unwrapped);

        /**
        * Given one frame of frequency-domain samples, return the
        * magnitudes, instantaneous phases, and unwrapped phases.
        *
        * reals and imags must each contain size/2+1 values (where size
        * is the frame size value as passed to the PhaseVocoder
        * constructor).
        *
        * mag, phase, and unwrapped must each be non-NULL and point to
        * enough space for size/2+1 values.
        */
        void processFrequencyDomain(const double *reals, const double *imags,
                                    double *mag, double *phase, double *unwrapped);

        /**
        * Reset the stored phases to zero. Note that this may be
        * necessary occasionally (depending on the application) to avoid
        * loss of floating-point precision in the accumulated unwrapped
        * phase values as they grow.
        */
        void reset();

    protected:
        void FFTShift(double *src);
        void getMagnitudes(double *mag);
        void getPhases(double *theta);
        void unwrapPhases(double *theta, double *unwrapped);

        int m_n;
        int m_hop;
        FFTReal *m_fft;
        double *m_time;
        double *m_imag;
        double *m_real;
        double *m_phase;
        double *m_unwrapped;
    };

    PhaseVocoder::PhaseVocoder(int n, int hop) :
        m_n(n),
        m_hop(hop)
    {
        m_fft = new FFTReal(m_n);
        m_time = new double[m_n];
        m_real = new double[m_n];
        m_imag = new double[m_n];
        m_phase = new double[m_n/2 + 1];
        m_unwrapped = new double[m_n/2 + 1];

        for (int i = 0; i < m_n/2 + 1; ++i) {
            m_phase[i] = 0.0;
            m_unwrapped[i] = 0.0;
        }

        reset();
    }

    PhaseVocoder::~PhaseVocoder()
    {
        delete[] m_unwrapped;
        delete[] m_phase;
        delete[] m_real;
        delete[] m_imag;
        delete[] m_time;
        delete m_fft;
    }

    void PhaseVocoder::FFTShift(double *src)
    {
        const int hs = m_n/2;
        for (int i = 0; i < hs; ++i) {
            double tmp = src[i];
            src[i] = src[i + hs];
            src[i + hs] = tmp;
        }
    }

    void PhaseVocoder::processTimeDomain(const double *src,
                                        double *mag, double *theta,
                                        double *unwrapped)
    {
        for (int i = 0; i < m_n; ++i) {
            m_time[i] = src[i];
        }
        FFTShift(m_time);
        m_fft->forward(m_time, m_real, m_imag);
        getMagnitudes(mag);
        getPhases(theta);
        unwrapPhases(theta, unwrapped);
    }

    void PhaseVocoder::processFrequencyDomain(const double *reals, 
                                            const double *imags,
                                            double *mag, double *theta,
                                            double *unwrapped)
    {
        for (int i = 0; i < m_n/2 + 1; ++i) {
            m_real[i] = reals[i];
            m_imag[i] = imags[i];
        }
        getMagnitudes(mag);
        getPhases(theta);
        unwrapPhases(theta, unwrapped);
    }

    void PhaseVocoder::reset()
    {
        for (int i = 0; i < m_n/2 + 1; ++i) {
            // m_phase stores the "previous" phase, so set to one step
            // behind so that a signal with initial phase at zero matches
            // the expected values. This is completely unnecessary for any
            // analytical purpose, it's just tidier.
            double omega = (2 * M_PI * m_hop * i) / m_n;
            m_phase[i] = -omega;
            m_unwrapped[i] = -omega;
        }
    }

    void PhaseVocoder::getMagnitudes(double *mag)
    {       
        for (int i = 0; i < m_n/2 + 1; i++) {
            mag[i] = sqrt(m_real[i] * m_real[i] + m_imag[i] * m_imag[i]);
        }
    }

    void PhaseVocoder::getPhases(double *theta)
    {
        for (int i = 0; i < m_n/2 + 1; i++) {
            theta[i] = atan2(m_imag[i], m_real[i]);
        }   
    }

    void PhaseVocoder::unwrapPhases(double *theta, double *unwrapped)
    {
        for (int i = 0; i < m_n/2 + 1; ++i) {

            double omega = (2 * M_PI * m_hop * i) / m_n;
            double expected = m_phase[i] + omega;
            double error = MathUtilities::princarg(theta[i] - expected);

            unwrapped[i] = m_unwrapped[i] + omega + error;

            m_phase[i] = theta[i];
            m_unwrapped[i] = unwrapped[i];
        }
    }

}