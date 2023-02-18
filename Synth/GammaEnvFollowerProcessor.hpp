#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"




namespace Gamma::Analysis
{

    // ports
    // value
    // done


// ports
    // value
    // done


    // Y >= 1 or done(X)  callback()

    struct EnvFollower : public FunctionProcessor, public gam::EnvFollow<DspFloatType,DspFloatType>
    {        
        EnvFollower(DspFloatType f=10.0) :
         FunctionProcessor(),
         gam::EnvFollow<DspFloatType,DspFloatType>(f)
        {

        }
        enum {
            PORT_VALUE,
            PORT_LAG,
            PORT_DONE,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_LAG: this->lag(v); break;                
            }
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_VALUE: return this->value();
                case PORT_DONE: return this->done();
            }
            return 0.0;
        }

        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=0.001, DspFloatType Y=0)
        {        
            return A*(*this)(I);
        }
    };

}