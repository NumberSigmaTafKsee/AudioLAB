#pragma once

#include "AudioFFT/AudioFFT.hpp"
#include "AudioFFT/FFTConvolver.hpp"

namespace KfrDSP1
{
    template<typename T>
    void zeros(kfr::univector<T> & v) {
        memset(v.data(),0,v.size()*sizeof(T));
    }
    
    ///////////////////////////////////////////////////////////////
    // AudioFFT
    ///////////////////////////////////////////////////////////////

    kfr::univector<kfr::complex<float>> audio_fft_forward(kfr::univector<float> in)
    {
        audiofft::AudioFFT fft;
        size_t s = in.size();
        if(s % 2 != 0)
        {
            s = (size_t)std::pow(2,std::log((double)s)+1.0);
        }
        kfr::univector<float> temp(s);
        kfr::univector<float> real(s);
        kfr::univector<float> imag(s);
        zeros(temp);
        zeros(real);
        zeros(imag);
        memcpy(temp.data(),in.data(),in.size()*sizeof(float));
        fft.init(in.size());
        fft.fft(temp.data(),real.data(),imag.data());
        kfr::univector<kfr::complex<float>> out(s);
        for(size_t i = 0; i < s; i++)
        {
            out[i].real(real[i]);
            out[i].imag(imag[i]);
        }
        return out;
    }

    kfr::univector<float> audio_fft_inverse(kfr::univector<kfr::complex<float>> in)
    {
        audiofft::AudioFFT fft;
        size_t s = in.size();
        if(s % 2 != 0)
        {
            s = (size_t)std::pow(2,std::log((double)s)+1.0);
        }
        kfr::univector<float> temp(s);
        kfr::univector<float> real(s);
        kfr::univector<float> imag(s);
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


    ///////////////////////////////////////////////////////////////
    // FFT Convolver
    ///////////////////////////////////////////////////////////////

    kfr::univector<float> audio_convolve(kfr::univector<float> in, size_t block_size, kfr::univector<float> ir)
    {
        kfr::univector<float> out(in.size());
        fftconvolver::FFTConvolver conv;
        conv.init(block_size,ir.data(),ir.size());
        conv.process(in.data(),out.data(),in.size());
        return out;
    }
}
