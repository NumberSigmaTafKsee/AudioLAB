#pragma once

namespace Gamma::FFT
{
    struct CFFTPadProcessor
    {
        gam::CFFT<DspFloatType> fft;
        std::vector<std::complex<DspFloatType>> temp;
        size_t block_size;

        CFFTPadProcessor(size_t blocks, size_t zeros=1) : fft(blocks*zeros)
        {
            temp.resize(blocks*2*zeros);
            memset((DspFloatType*)temp.data(),0,2*blocks*zeros*sizeof(DspFloatType)*2);
            block_size = blocks;
        }
        void forward(std::vector<std::complex<DspFloatType>> & input,std::vector<std::complex<DspFloatType>> & output, bool normalize=true, DspFloatType nrmGain=1.0f )
        {
            temp = input;
            for(size_t i = input.size(); i < temp.size(); i++) temp[i] = std::complex<DspFloatType>(0,0);
            fft.forward(static_cast<DspFloatType*>(temp.data()),normalize,nrmGain);
            for(size_t i = 0; i < block_size; i++) output[i] = temp[i];
        }
        void inverse(std::vector<std::complex<DspFloatType>> & input,std::vector<std::complex<DspFloatType>> & output)
        {
            temp = input;
            for(size_t i = input.size(); i < temp.size(); i++) temp[i] = std::complex<DspFloatType>(0,0);
            fft.inverse(static_cast<DspFloatType*>(temp.data()));
            for(size_t i = 0; i < block_size; i++) output[i] = temp[i];
        }
    };

    struct RFFTPadProcessor
    {
        gam::RFFT<DspFloatType> fft;
        std::vector<DspFloatType> temp;
        std::vector<std::complex<DspFloatType>> tempc;
        size_t block_size;

        RFFTPadProcessor(size_t blocks, size_t zeros=1) : fft(blocks*zeros)
        {
            temp.resize(blocks*2*zeros);
            memset((DspFloatType*)temp.data(),0,2*blocks*zeros*sizeof(DspFloatType));
            block_size = blocks;
        }
        void forward(std::vector<DspFloatType> & input,std::vector<std::complex<DspFloatType>> & output, bool normalize=true, DspFloatType nrmGain=1.0f )
        {
            temp = input;
            for(size_t i = input.size(); i < temp.size(); i++) temp[i] = std::complex<DspFloatType>(0,0);
            fft.forward(static_cast<DspFloatType*>(temp.data()),normalize,nrmGain);
            for(size_t i = 0; i < block_size; i++) {
                output[i].real(temp[2*i]);
                output[i].imag(temp[2*i+1]);
            }
        }
        void inverse(std::vector<std::complex<DspFloatType>> & input,std::vector<DspFloatType> & output)
        {
            tempc = input;
            for(size_t i = input.size(); i < tempc.size(); i++) tempc[i] = std::complex<DspFloatType>(0,0);
            fft.inverse(static_cast<DspFloatType*>(tempc.data()));
            memcpy(output,(DspFloatType*)tempc.data(),block_size*sizeof(DspFloatType));            
        }
    };
}