#pragma once
#include "SoundPipe.hpp"

namespace SoundPipe
{
    struct FOsc : public SoundPipe
    {    
        sp_fosc * fosc;
        // have to fill this table with the impulse response
        
        FOsc(sp_data * data, sp_ftbl *ft) : SoundPipe(data)        
        {                    
            sp_fosc_create(&fosc);
            sp_fosc_init(data,osc, ft);
            
        }
        ~FOsc() {
            if(fosc) sp_filt_destroy(&fosc);            
        }
        void setFreq(float f) {
            fosc->freq = f;
        }
        void setAmp(float a) {
            fosc->amp = a;
        }
        void setCar(float v) {
            fosc->car =v;
        }
        void setMod(float v) {
            fosc->mod =v;
        }
        void setIndx(float v) {
            fosc->indx =v;
        }
        void setIphs(float v) {
            fosc->iphs =v;
        }

        float Tick(float I,float A=1, float X = 1, float Y = 1) {
            float in = 0;
            float out = 0;                
            sp_fosc_compute(sp,fosc,&in,&out);            
            return A * out;
        }
    };
}