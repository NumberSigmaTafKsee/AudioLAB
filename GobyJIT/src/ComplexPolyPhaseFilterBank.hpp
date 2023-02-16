//https://github.com/alexranaldi/PolyphaseFilterBank/blob/main/PolyPhaseFilterBank.hpp
#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>


// Radix-2 FFT forward transform
class FFT_R2 {

public:

    // Creates an FFT object
    FFT_R2(size_t N_log2) : N(1 << N_log2), oddPow(N_log2 & 1) {
        permute = new size_t[N];
        for (size_t n = 0, r = 0; n < N; ++n, r = 0) {
            for (size_t m = 0; m < N_log2; ++m) r = (r << 1) | ((n >> m) & 1);
            permute[n] = r; // reversed bits
        }
        phasors = new std::complex<double>[N/2];
        for (size_t n = 0; n < N/2; ++n) {
            double x = 2.0*PI*(double)n/(double)N;
            phasors[n] = std::complex<double>(cos(x),-sin(x)); // conjugate
        }
        workspace = new std::complex<double>[N];
    }

    // Destroy the FFT object
    ~FFT_R2(void) {
        delete[] permute;
        delete[] phasors;
        delete[] workspace;
    }

    // Perform the forward transform (out of place) on complex input
    //  The two pointers must hold at least N samples each and not overlap!
    void forwardOOP(std::complex<double> *in, std::complex<double> *out) {
        if (2 < N) {
            // We'll be using the workspace vector and the output vector to
            // alternate computations. We need to choose the right one to use
            // for the first stage so that the final stage is the output.
            std::complex<double> *x = out, *y = out, *swap;
            if (oddPow) {
                x = workspace;
            } else {
                y = workspace;
            }
            // First stage
            for (size_t n = 0; n < N; n += 2) {
                y[n]   = in[permute[n]] + in[permute[n+1]];
                y[n+1] = in[permute[n]] - in[permute[n+1]];
            }
            // Remaining stages
            std::complex<double> temp;
            for (size_t M = 4; M <= N; M *= 2) {
                swap = x; x = y; y = swap; // swap ptrs
                for (size_t n = 0; n < N; n += M) {
                    y[n]     = x[n] + x[n+M/2];
                    y[n+M/2] = x[n] - x[n+M/2];
                    for (size_t m = n+1, k = N/M; m < n+M/2; ++m, k += N/M) {
                        temp = x[m+M/2]*phasors[k];
                        y[m]     = x[m] + temp;
                        y[m+M/2] = x[m] - temp;
                    }
                }
            }
        } else if (2 == N) {
            out[0] = in[0] + in[1];
            out[1] = in[0] - in[1];
        } else if (1 == N) {
            out[0] = in[0];
        }
    }

private:

    //
    size_t N;
    bool oddPow;
    size_t *permute;
    std::complex<double> *phasors;
    std::complex<double> *workspace;

}; // class FFT_R2


// An efficient implementation of a poly-phase filter bank
//   - Uses a power-of-two based sinc filters (many coeffs are 0).
//   - Uses a sample buffer to aid in processing consecutive samples (can be
//   cleared if consecutive calls contain data that is not contiguous).
//   - The output samples are delayed by FiltOrd/(NumChan*2)-1 samples.
//   - Note that the output sample rate is 1/NumChan of the input rate.
//   - Filter bandwidth scaling can be used for filter overlap
class PolyPhaseFilterBank {

public:

    // Constructor
    PolyPhaseFilterBank(
            size_t NumChan_log2,
            size_t FiltOrd_log2,
            double FiltScale = 1.0);

    // Destructor
    ~PolyPhaseFilterBank(void);

    // Consumes the next NumChan samples from the provided pointer and
    // generates the NumChan outputs (delayed and unnormalized).
    // In and out pointers can be identical (in-place operation).
    void Iterate(
            std::complex<double> *in,
            std::complex<double> *out);

    void ProcessBlock(size_t n, float * in, float * out)
    {
        std::vector<std::complex<double>> I(n);
        std::vector<std::complex<double>> O(n);
        for(size_t i = 0; i < n; i++) {
            O[i] = 0;
            I[i].real(in[i]);
            I[i].imag(0);
        }
        Iterate(I.data(),O.data());
        for(size_t i = 0; i < n; i++) {
            out[i] = std::abs(O[i]);

        }        
    }
    
    // Resets the internal buffers to zero, losing history
    void ResetBuffer(void);

    // Returns a scalar which should be multiplied by the output as:
    //  output/sqrt(GetNorm())
    // It is not multiplied internally for computational savings
    inline double GetNorm(void) const {return norm;}

    // Returns the fractional bandwidth (1/gain) of the filters
    inline double GetBandwidth(void) const {return Scale/NumChan;}

    // Returns the output delay of this filter (at the output rate)
    inline size_t GetDelay(void) const {return FiltOrd/NumChan/2-1;}

    // Returns the number of channels
    inline size_t GetNumChan(void) const {return NumChan;}

private:

    // The number of channels
    const size_t NumChan;
    // The low-pass filter order
    const size_t FiltOrd;
    // Scalar for filter bandwidth
    const double Scale;
    // The normalization scalar
    double norm;
    // The weighted low-pass filter
    double *filter;
    // Buffer for samples
    std::complex<double> *samples;
    // Index to place next sample
    size_t next;
    // The FFT object
    FFT_R2 fft;
    // FFT workspace
    std::complex<double> *workspace;

}; // class PolyPhaseFilterBank


PolyPhaseFilterBank::PolyPhaseFilterBank(
        size_t NumChan_log2, size_t FiltOrd_log2, double FiltScale) :
    NumChan(1 << NumChan_log2), FiltOrd(1 << FiltOrd_log2),
    Scale(FiltScale), next(0), fft(NumChan_log2)
{

    // Check sizes
    if (1 + NumChan_log2 > FiltOrd_log2)
        throw std::invalid_argument("FiltOrd/NumChan must be >= 2");
    if (16 < FiltOrd_log2)
        throw std::invalid_argument("log2(FiltOrd) must be <= 16");
    if (1.0 > FiltScale || 2.0 < FiltScale)
        throw std::invalid_argument("FiltScale must be in [1,2]");

    // Create the filter taps and compute the normalization factor.
    //  Note that the number of taps in the filter is order+1 but since
    //  the last coeff would be 0 so we don't need the +1.
    //  The first tap is 0, the middle tap is 1, and the taps are symmetrical.
    //  Additionally, many taps are 0 (every NumChan/Scale)
    filter = new double[FiltOrd];
    filter[0] = 0.0;
    filter[FiltOrd/2] = 1.0;
    double sumY = 1.0, sumZ2 = 1.0, sumY2 = 1.0; // 1.0 for middle tap
    double xd = Scale/NumChan, x0 = -xd*FiltOrd/2;
    for (size_t n = 1; n < FiltOrd/2; ++n) { // filter[0] = 0 so skip
        double xPi = (x0+xd*n)*PI;
        double y = sin(xPi)/xPi; // un-weighted
        double z = y * (0.54-0.46*cos(2.0*PI*(double)n/(double)FiltOrd));
        filter[n] = z;
        filter[FiltOrd-n] = z;
        sumY  += 2.0*y;
        sumY2 += 2.0*y*y; // x2 for left and right sides
        sumZ2 += 2.0*z*z; // x2 for left and right sides
    }

    // Compute the normalization scalar
    norm = sqrt(sumY/NumChan/Scale) *   // sinc extent
            sumZ2/sumY2 *               // Hamming Window
            NumChan;                    // FFT
    if (2 == FiltOrd/NumChan) norm *= 0.83176;

    // Create and clear the sample buffer
    samples = new std::complex<double>[FiltOrd];
    ResetBuffer();

    // Create the FFT workspace
    workspace = new std::complex<double>[NumChan];

}


PolyPhaseFilterBank::~PolyPhaseFilterBank(void) {
    delete[] filter;
    delete[] samples;
    delete[] workspace;
}


void PolyPhaseFilterBank::ResetBuffer(void) {
    std::fill(samples, samples+FiltOrd, std::complex<double>(0.0, 0.0));
}


void PolyPhaseFilterBank::Iterate(
        std::complex<double> *in, std::complex<double> *out)
{

    // Rather than performing modulo operations, we can use a bit-mask since
    // our base is a power of two (significant CPU savings)
    const size_t MASK = FiltOrd-1;

    // Put the new samples in the buffer
    std::copy(in, in+NumChan, samples+next);
    next = (next + NumChan) & MASK; // "next" now points to 1st sample

    // Filter the samples and add the blocks together
    for (size_t n = 0; n < NumChan; ++n) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t m = n; m < FiltOrd; m += NumChan)
            sum += filter[m]*samples[(m+next) & MASK];
        workspace[n] = sum;
    }

    // Perform forward FFT
    fft.forwardOOP(workspace, out);

}