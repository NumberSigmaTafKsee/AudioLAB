#pragma once
#include "SoundPipe.hpp"

namespace SoundPipe
{
    struct Filt : public SoundPipe
    {    
        sp_filt * filt;
        // have to fill this table with the impulse response
        
        Filt(sp_data * data) : SoundPipe(data)        
        {                    
            sp_filt_create(&filt);
            sp_filt_init(data,filt);
            
        }
        ~Filt() {
            if(conv) sp_filt_destroy(&conv);            
        }
        void setFreq(float f) {
            filt->freq = f;
        }
        void setAtk(float a) {
            filt->atk = a;
        }
        void setDec(float d) {
            filt->dec = d;
        }

        float Tick(float I,float A=1, float X = 1, float Y = 1) {
            float in = 0;
            float out = 0;    
            float c = filt->freq;
            filt->freq = filt->freq*X*Y;                    
            sp_filt_compute(sp,filt,&in,&out);
            filt->freq = c;
            return A * out;
        }
    };
}