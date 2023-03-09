#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct FormSwep : public GeneratorProcessorPlugin<stk::FormSwep>
    {
        FormSwep() : GeneratorProcessorPlugin<stk::FormSwep>()
        {

        }
        enum {
            PORT_RESONANCE,
            PORT_STATES,
            PORT_TARGETS,
            PORT_SWEEPRATE,
            PORT_SWEEPTIME,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {                
                case PORT_SWEEPRATE: this->setSweepRate(v); break;
                case PORT_SWEEPTIME: this->setSweepTime(v); break;
            }
        }
        void setPort2(int port, DspFloatType a, DspFloatType b)
        {
            if(port != PORT_RESONANCE) return;
            this->setResonance(a,b);
        }
        void setPortV(int port, const std::vector<DspFloatType> & v) {
            switch(port) {
                case PORT_STATES: this->setStates(v[0],v[1],v[2]); break;
                case PORT_TARGETS: this->setTargets(v[0],v[1],v[2]); break;
            }
        }
        DspFloatType Tick(DspFloatType I, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return A*this->tick(I);
        }
    };
}