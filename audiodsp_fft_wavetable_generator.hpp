#pragma once
#include <vector>
#include <iostream>
#include <fftw3.h>

namespace AudioDSP
{
    struct C2RF
    {
        fftwf_complex * in;    
        float * out;
        size_t size;
        fftwf_plan p;

        C2RF() {
            in = NULL;
            out = NULL;
            size = 0;
        }
        C2RF(size_t n) {
            in = NULL;
            out = NULL;
            size = 0;
            init(n);
        }
        ~C2RF() {
            fftwf_destroy_plan(p);
            fftwf_free(in);
            fftwf_free(out);    
        }
        void init(size_t n) {
            if(in != NULL) fftwf_destroy_plan(p);
            if(in != NULL) fftwf_free(in);
            if(out != NULL) fftwf_free(out);
            in = NULL;
            out = NULL;
            size = n;
            in = fftwf_alloc_complex(n);
            out= fftwf_alloc_real(n);                    
            p = fftwf_plan_dft_c2r_1d(n, in, out, FFTW_ESTIMATE);
        }
        void set_input(std::vector<std::complex<float>> & input) {
            for(size_t i = 0; i < size; i++) {
                in[i][0] = input[i].real();
                in[i][1] = input[i].imag();
            }
        }        
        std::vector<float> get_output() {
            std::vector<float> r(size);
            memcpy(r.data(),out, size*sizeof(float));
            return r;
        }        
        void get_output( std::vector<float> & output)
        {
            if(output.size() != size) output.resize(size);
            for(size_t i = 0; i < size; i++ )
            {
                output[i] = out[i];
            }
        }
        void normalize()
        {
            for(size_t i = 0; i < size; i++) 
                out[i] /= (float)size;                
        }
        void Execute() {
            fftwf_execute(p);            
        }
    };
    struct FFTWaveTableGenerator
    {
        // number of harmonics for note 
        // 0 = DC
        // 44100/4096 = 10.766
        // 1 = f0
        // 2 = f1
        // 4 = f2
        // 8 = f3

        static std::vector<DspFloatType> sawtooth(DspFloatType f, DspFloatType sr)
        {
            std::vector<std::complex<DspFloatType>> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/2)/f;
            //std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm; i++)
            {
                DspFloatType n = 1/(DspFloatType)i;
                buffer[i] = std::complex<DspFloatType>(0,n);
            }                        
            std::vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 1/(M_PI);            
            return out;
        }
        static std::vector<DspFloatType> reverse_sawtooth(DspFloatType f, DspFloatType sr)
        {
            std::vector<std::complex<DspFloatType>> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            size_t harm = (sr/2)/f;            
            //std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm; i++)
            {
                DspFloatType n = 1/(DspFloatType)i;
                buffer[i] = std::complex<DspFloatType>(0,-n);
            }                        
            std::vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 1/(M_PI);            
            return out;
        }
        static std::vector<DspFloatType> square(DspFloatType f, DspFloatType sr)
        {
            std::vector<std::complex<DspFloatType>> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/2)/f;
            //std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm; i+=2)
            {
                DspFloatType n = 1/(DspFloatType)i;
                buffer[i] = std::complex<DspFloatType>(0,-n);
            }                        
            std::vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 2.0/(M_PI);
            return out;
        }
        static std::vector<DspFloatType> triangle(DspFloatType f, DspFloatType sr)
        {
            std::vector<std::complex<DspFloatType>> buffer;
            std::complex<DspFloatType> temp(0,-M_PI);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/2)/f;
            //std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm; i+=2)
            {
                DspFloatType n = 1.0/(DspFloatType)(i*i);                       
                buffer[i] = std::complex<DspFloatType>(n,0)*exp(temp);
            }                        
            std::vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 4.0/(M_PI*M_PI);            
            return out;
        }
        static std::vector<DspFloatType> sine(DspFloatType f, DspFloatType sr)
        {
            std::vector<std::complex<DspFloatType>> buffer;
            size_t size = 4096;
            buffer.resize(size);                        
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            buffer[1] = std::complex<DspFloatType>(0,-1);
            std::vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            return out;
        }
        static std::vector<DspFloatType> cyclize(std::vector<std::complex<DspFloatType>> & c)
        {            
            C2RF inverse(c.size());
            inverse.set_input(c);
            inverse.Execute();
            std::vector<DspFloatType> out;
            inverse.get_output(out);
            return out;
        }
    };
}