#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::Spatial
{

    struct LoopGain : public FunctionProcessorPlugin<gam::LoopGain<DspFloatType>>
    {

    };

    struct Loop1P : public FunctionProcessorPlugin<gam::Loop1P<DspFloatType>>
    {

    };

    struct Loop1P1Z : public FunctionProcessorPlugin<gam::Loop1P1Z<DspFloatType>>
    {

    };

    struct Echo : public FunctionProcessorPlugin<gam::Echo<DspFloatType>>
    {

    };
    struct EchoCSine : public FunctionProcessorPlugin<gam::EchoCSine<DspFloatType>>
    {

    };    
    struct ReverbMS : public FunctionProcessorPlugin<gam::ReverbMS<DspFloatType>>
    {

    };
    struct Dist : public FunctionProcessorPlugin<gam::Dist<2,DspFloatType>>
    {

    };
    
}
