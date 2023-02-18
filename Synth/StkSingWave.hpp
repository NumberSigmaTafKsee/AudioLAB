#pragma once
#include "StkHeaders.hpp"
namespace Stk::Generators
{
    struct SingWave : public GeneratorProcessor, public stk::SingWave
    {
        SingWave(const std::string& filename, bool raw=false) 
        : GeneratorProcessor(),stk::SingWave(filename,raw)
        {

        }

        DspFloatType Tick(DspFloatType I=0, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            return this->tick();
        }
    };

}