#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct Noise : public GeneratorProcessorPlugin<stk::Noise>
    {
        Noise() : GeneratorProcessorPlugin<stk::Noise>()
        {

        }
        enum {
            PORT_SEED,
            PORT_LASTOUT,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_SEED: this->setSeed(v); break;
                default: printf("No port %d\n",port);
            }
        }
        DspFloatType getPort(int port) {
            if(port == PORT_LASTOUT) return this->lastOut();
            printf("No port %d\n",port);
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };
}