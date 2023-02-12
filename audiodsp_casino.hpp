#pragma once

#include "DSP/Carlo/carlo_casino.hpp"

namespace AudioDSP
{
    using ArrayXf    = Casino::IPP::IPPArray<float>;
    using ArrayXd    = Casino::IPP::IPPArray<double>;
    using ArrayXcf    = Casino::IPP::IPPArray<std::complex<float>>;
    using ArrayXcd    = Casino::IPP::IPPArray<std::complex<double>>;
    using MatrixXcf  = Casino::MKL::ComplexMatrix<float>;
    using MatrixXcd  = Casino::MKL::ComplexMatrix<double>;
    using MatrixXf   = Casino::MKL::Matrix<float>;
    using MatrixXd   = Casino::MKL::Matrix<double>;
    using VectorXcf  = Casino::MKL::ComplexVector<float>;
    using VectorXcd  = Casino::MKL::ComplexVector<double>;
    using VectorXf   = Casino::MKL::Vector<float>;
    using VectorXd   = Casino::MKL::Vector<double>;

    using autocorr_f32 = Casino::IPP::AutoCorr<float>;
    using autocorr_f64 = Casino::IPP::AutoCorr<double>;
    using xcorr_f32 = Casino::IPP::CrossCorr<float>;
    using xcorr_f64 = Casino::IPP::CrossCorr<double>;
    using convolve_f32 = Casino::IPP::Convolver<float>;
    using convolve_f64 = Casino::IPP::Convolver<double>;
    using convfilter_f32 = Casino::IPP::ConvolutionFilter<float>;
    using convfilter_f64 = Casino::IPP::ConvolutionFilter<double>;
    using dct_f32 = Casino::IPP::DCT<float>;
    using dct_f64 = Casino::IPP::DCT<double>;
    using dft_f32 = Casino::IPP::DFT<float>;
    using dft_f64 = Casino::IPP::DFT<double>;
    using fft_f32 = Casino::IPP::CFFT<float>;
    using fft_f64 = Casino::IPP::CFFT<double>;
    using firlms_f32 = Casino::IPP::FIRLMS32;
    using firmr_f32  = Casino::IPP::FIRMR<float>;
    using firmr_f64  = Casino::IPP::FIRMR<double>;
    using firsr_f32  = Casino::IPP::FIRSR<float>;
    using firsr_f64  = Casino::IPP::FIRSR<double>;
    using hilbert_f32 = Casino::IPP::HilbertTransform32;
    using hilbert_f64 = Casino::IPP::HilbertTransform64;
    using iir_f32 = Casino::IPP::IIR<float>;
    using iir_f64 = Casino::IPP::IIR<double>;
    using biquad_f32 = Casino::IPP::IIRBiquad<float>;
    using biquad_f64 = Casino::IPP::IIRBiquad<double>;
    using random_f32 = Casino::IPP::RandomUniform<float>;
    using random_f64 = Casino::IPP::RandomUniform<double>;
    using resample_f32 = Casino::IPP::Resample32;
}