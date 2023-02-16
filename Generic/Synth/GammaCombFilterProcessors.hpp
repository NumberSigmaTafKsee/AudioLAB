#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::Delay
{
    // todo: other interpolations
    struct Comb : public FunctionProcessor, public gam::Comb<DspFloatType>
    {
        Comb() : FunctionProcessor(),gam::Comb<DspFloatType>()
        {

        }
        Comb(float delay, float ffd=0.0, float fbk=0.0) 
        : FunctionProcessor(),
          gam::Comb<DspFloatType>(delay,ffd,fbk)
        {

        }
        Comb(float maxDelay, float delay, float ffd=0.0, float fbk=0.0) 
        : FunctionProcessor(),
          gam::Comb<DspFloatType>(maxDelay,delay,ffd,fbk)
        {

        }
        enum {
            PORT_ALLPASS,
            PORT_FBK,
            PORT_FFD,  
            PORT_DECAY,          
            PORT_NEXTFBK,
            PORT_NORM,
            PORT_NORMFBK,
            PORT_NORMFFD,    
            //PORT_FBW,        
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_ALLPASS: this->allPass(v); break;
                case PORT_FBK: this->fbk(v); break;
                case PORT_FFD: this->ffd(v); break;
                case PORT_DECAY: this->decay(v); break;
            }
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_FBK: return this->fbk(); 
                //case PORT_FBW: return this->fbw();     
                case PORT_NORM: return this->norm(); 
                case PORT_NORMFBK: return this->normFbk(); 
                case PORT_NORMFFD: return this->normFfd(); 
            }
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {
            return (*this)(I);
        }
    };
}