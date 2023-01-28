#pragma once
#include "SoundPipe.hpp"

namespace SoundPipe
{
    struct GBuzz : public SoundPipe
    {    
        sp_gbuzz * gbuzz;
        // have to fill this table with the impulse response
        
        GBuzz(sp_data * data) : SoundPipe(data)        
        {                    
            sp_gbuzz_create(&gbuzz);
            sp_gbuzz_init(data,gbuzz);            
        }
        ~GBuzz() {
            if(conv) sp_gbuzz_destroy(&gbuzz);            
        }
        void setFreq(float f) {
            gbuzz->freq = f;
        }
        void setAmp(float a) {
            gbuzz->amp = a;
        }
        void setNharm(float v) {
            gbuzz->nharm =v;
        }
        void setLharm(float v) {
            gbuzz->lharm =v;
        }
        void setMul(float v) {
            gbuzz->mul =v;
        }
        void setIphs(float v) {
            gbuzz->iphs =v;
        }
        float Tick(float I,float A=1, float X = 1, float Y = 1) {
            float in = 0;
            float out = 0;    
            float c = filt->freq;
            filt->freq = filt->freq*X*Y;                    
            sp_gbuzz_compute(sp,filt,&in,&out);
            filt->freq = c;
            return A * out;
        }
    };
}