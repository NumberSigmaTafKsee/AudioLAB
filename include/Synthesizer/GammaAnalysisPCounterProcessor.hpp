#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"

namespace Gamma::Analysis
{
    // ports
    // period n
    // period
    // cycled
    // count
    // reset

    // X > 0  period(X)
    // Y >=1  reset()

    struct PCounter : public GeneratorProcessor, public gam::PCounter
    {
        std::function<void (PCounter &)> callback;

        PCounter(unsigned period = 256)
        : GeneratorProcessor(),gam::PCounter(period)
        {
            callback =  [](PCounter &) {};
        }
        enum {
            PORT_PERIOD,
            PORT_CYCLED,
            PORT_COUNT,
            PORT_RESET,
        };
        void setPort(int port, DspFloatType v) {
            if(port == PORT_RESET && v == 1.0) this->reset();
        }
        DspFloatType getPort(int port) {
            switch(port) {
                case PORT_PERIOD: return this->period();
                case PORT_CYCLED: return this->cycled();
                case PORT_COUNT: return this->count();                
            }
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0)
        {           
            return (*this)();
        }
    };
}