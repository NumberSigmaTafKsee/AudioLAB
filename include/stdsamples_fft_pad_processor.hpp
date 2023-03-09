#pragma once

namespace FFT
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
}