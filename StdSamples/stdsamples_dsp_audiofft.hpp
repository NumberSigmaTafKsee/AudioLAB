#pragma once

#include "AudioFFT/AudioFFT.h"
#include "AudioFFT/FFTConvolver.h"


namespace AudioDSP
{
    ////////////////////////////////////////////////////////////////
    // AudioFFT
    ////////////////////////////////////////////////////////////////
    complex_vector<float> audiofft_forward(sample_vector<float> & in)
    {
        audiofft::AudioFFT fft;
        size_t s = in.size();
        if(s % 2 != 0)
        {
            s = (size_t)std::pow(2,std::log((double)s)+1.0);
        }
        sample_vector<float> temp(s);
        sample_vector<float> real(s);
        sample_vector<float> imag(s);
        zeros(temp);
        zeros(real);
        zeros(imag);
        memcpy(temp.data(),in.data(),in.size()*sizeof(float));
        fft.init(in.size());
        fft.fft(temp.data(),real.data(),imag.data());
        complex_vector<float> out(s);
        for(size_t i = 0; i < s; i++)
        {
            out[i].real(real[i]);
            out[i].imag(imag[i]);
        }
        return out;
    }

    sample_vector<float> audiofft_inverse(complex_vector<float> & in)
    {
        audiofft::AudioFFT fft;
        size_t s = in.size();
        if(s % 2 != 0)
        {
            s = (size_t)std::pow(2,std::log((double)s)+1.0);
        }
        sample_vector<float> temp(s);
        sample_vector<float> real(s);
        sample_vector<float> imag(s);
        zeros(temp);
        zeros(real);
        zeros(imag);
        for(size_t i = 0; i < in.size(); i++)
        {
            real[i] = in[i].real();
            imag[i] = in[i].imag();
        }
        fft.init(in.size());
        fft.ifft(temp.data(),real.data(),imag.data());        
        return temp;
    }

    ////////////////////////////////////////////////////////////////
    // FFTConvolver
    ////////////////////////////////////////////////////////////////
    sample_vector<float> fftconvolver(sample_vector<float>& in, size_t block_size, sample_vector<float> & ir)
    {
        sample_vector<float> out(in.size());
        fftconvolver::FFTConvolver conv;
        conv.init(block_size,ir.data(),ir.size());
        conv.process(in.data(),out.data(),in.size());
        return out;
    }
 }