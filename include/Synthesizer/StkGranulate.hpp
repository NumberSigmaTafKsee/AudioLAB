#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    // todo: ports
    struct Granulate : public GeneratorProcessorPlugin<stk::Granulate>
    {
        Granulate(const char * filename, bool raw=false) : GeneratorProcessorPlugin<stk::Granulate>()
        {
            this->openFile(filename, raw);
        }
        enum {
            PORT_RESET,
            PORT_VOICES,
            PORT_STRETCH,
            PORT_GRAINPARAMETERS,
            PORT_RANDOMFACTOR,
            PORT_LASTOUT,
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_RESET: this->reset(); break;
                case PORT_VOICES: this->setVoices(v); break;
                case PORT_STRETCH: this->setStretch(v); break;
                case PORT_RANDOMFACTOR: this->setRandomFactor(v); break;
            }
        }
        void setPortV(int port, const std::vector<DspFloatType> & v) {
            switch(port) {
                case PORT_GRAINPARAMETERS: this->setGrainParameters(v[0],v[1],v[2],v[3]); break;                
            }
        }
        DspFloatType getPort(int port) {
            if(port == PORT_LASTOUT) return this->lastOut();
            return 0.0;
        }
        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };
}