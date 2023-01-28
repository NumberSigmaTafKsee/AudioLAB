#pragma once

namespace FFT
{
    struct FFTWaveTableGenerator
    {
        // number of harmonics for note 
        // 0 = DC
        // 44100/4096 = 10.766
        // 1 = f0
        // 2 = f1
        // 4 = f2
        // 8 = f3

        static sample_vector<DspFloatType> sawtooth(DspFloatType f, DspFloatType sr)
        {
            complex_vector<DspFloatType> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/f)/(sr/size);
            std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm-1; i++)
            {
                DspFloatType n = 1/(DspFloatType)i;
                buffer[i] = std::complex<DspFloatType>(0,n);
            }                        
            sample_vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 1/(M_PI);            
            return out;
        }
        static sample_vector<DspFloatType> reverse_sawtooth(DspFloatType f, DspFloatType sr)
        {
            complex_vector<DspFloatType> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/f)/(sr/size);
            std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm-1; i++)
            {
                DspFloatType n = 1/(DspFloatType)i;
                buffer[i] = std::complex<DspFloatType>(0,-n);
            }                        
            sample_vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 1/(M_PI);            
            return out;
        }
        static sample_vector<DspFloatType> square(DspFloatType f, DspFloatType sr)
        {
            complex_vector<DspFloatType> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/f)/(sr/size);
            std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm-1; i+=2)
            {
                DspFloatType n = 1/(DspFloatType)i;
                buffer[i] = std::complex<DspFloatType>(0,-n);
            }                        
            sample_vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 2.0/(M_PI);
            return out;
        }
        static sample_vector<DspFloatType> triangle(DspFloatType f, DspFloatType sr)
        {
            complex_vector<DspFloatType> buffer;
            std::complex<DspFloatType> temp(0,-M_PI);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/f)/(sr/size);
            std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            for(size_t i=1; i < harm-1; i+=2)
            {
                DspFloatType n = 1.0/(DspFloatType)(i*i);                       
                buffer[i] = std::complex<DspFloatType>(n,0)*exp(temp);
            }                        
            sample_vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            for(size_t i = 0; i < size; i++) out[i] *= 4.0/(M_PI*M_PI);            
            return out;
        }
        static sample_vector<DspFloatType> sine(DspFloatType f, DspFloatType sr)
        {
            complex_vector<DspFloatType> buffer;
            std::complex<DspFloatType> temp = (0,-1);
            size_t size = 4096;
            buffer.resize(size);
            size_t harm = (sr/f)/(sr/size);
            std::cout << harm << std::endl;
            memset(buffer.data(),0,size*sizeof(std::complex<DspFloatType>));            
            buffer[1] = std::complex<DspFloatType>(0,-1);
            sample_vector<DspFloatType> out(4096);
            C2RF inverse(4096);
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(out);
            return out;
        }
        static sample_vector<DspFloatType> cyclize(complex_vector<DspFloatType> & c)
        {            
            C2RF inverse(c.size());
            inverse.set_input(c);
            inverse.Execute();
            sample_vector<DspFloatType> out;
            inverse.get_output(out);
            return out;
        }
    };
}