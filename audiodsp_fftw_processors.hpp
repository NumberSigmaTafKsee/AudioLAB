#pragma once

namespace Spectrum::FFT
{
    /////////////////////////////////////////////////////////////////////////////////////
    // Pad FFT with zeros
    /////////////////////////////////////////////////////////////////////////////////////
    struct FFTPadProcessor
    {
        sample_vector<float> buffer;    
        int32_t bufptr;
        uint32_t blocksize;
        Spectrum::R2CF forward;
        Spectrum::C2RF inverse;
        FFTPadProcessor(uint32_t blocks, uint32_t zero=2) {            
            bufptr = 0;        
            blocksize = blocks;        
            buffer.resize(blocksize*2*zero);
            forward.init(blocksize*2*zero);
            inverse.init(blocksize*2*zero);        
            zeros(buffer);
        }
        void ProcessBlock(size_t hopsize, float * input, float * output) {
            Spectrum::complex_vector<float> x;
            sample_vector<float> out;                
            memcpy(buffer.data(),input,hopsize*sizeof(float));
            forward.set_input(buffer);
            forward.Execute();
            forward.normalize();
            x = forward.get_output();
            inverse.set_input(x);
            inverse.Execute();
            out = inverse.get_output();     

            for(size_t i = 0; i < hopsize; i++)
            {            
                output[i]  = out[i];
            }                                                         
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

    
    // experiment
    struct FFTOverlapAdd
    {
        Spectrum::complex_vector<float> buffer;
        sample_vector<float> temp,ola;
        Spectrum::Rectangle<float> window;
        Spectrum::R2CF forward;
        Spectrum::C2RF inverse;
        size_t hop;

        FFTOverlapAdd(size_t size=1024, size_t hop_size=128) : window(size),hop(hop_size) {                
            ola.resize(size);
            buffer.resize(size);
            temp.resize(size);
            forward.init(size);
            inverse.init(size);
            zeros(buffer);
            zeros(temp);            
        }
        void ProcessBlock(size_t n, float * input, float * output)
        {
            memcpy(temp.data(),input,n*sizeof(float));
            
            for(size_t i = 0; i < n; i++)
                ola.push_back(input[i]);
            for(size_t i = 0; i < hop; i++)
            {                
                temp[i] += ola[i];                
            }
            ola.erase(ola.begin(),ola.begin()+hop);
            for(size_t i = n; i < ola.size(); i++) ola[i] = 0;
            for(size_t i = 0; i < temp.size(); i++)
            {            
                temp[i] *= window[i];            
            }        
            forward.set_input(temp);        
            forward.Execute();
            forward.normalize();
            forward.get_output(buffer);                                    
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(temp.data());
            for(size_t i = 0; i < n; i++)       
            {     
                output[i] = temp[i];
            }            
        }
    };

    // experiment
    struct FFTShifter
    {
        Spectrum::complex_vector<float> buffer;
        sample_vector<float> temp,ola;
        Spectrum::Rectangle<float> window;
        Spectrum::R2CF forward;
        Spectrum::C2RF inverse;

        FFTShifter(size_t size=1024) : window(size) {                
            ola.resize(size);
            buffer.resize(size);
            temp.resize(size);
            forward.init(size);
            inverse.init(size);
            zeros(buffer);
            zeros(temp);            
        }
        void ProcessBlock(size_t n, float * input, float * output)
        {
            for(size_t i = 0; i < n; i++)
                ola[ola.size()-n+i] = input[i];
            
            memcpy(temp.data(),ola.data(),temp.size()*sizeof(float));
            
            for(size_t i = 0; i < temp.size(); i++)
            {            
                temp[i] *= window[i];            
            }        
            forward.set_input(temp);        
            forward.Execute();
            forward.normalize();
            forward.get_output(buffer);                                    
            inverse.set_input(buffer);
            inverse.Execute();
            inverse.get_output(temp.data());
            for(size_t i = 0; i < n; i++)                   
                output[i] = temp[i];            
            for(size_t i = n; i < ola.size(); i++)
                ola[i-n] = ola[i];
        }
    };
}