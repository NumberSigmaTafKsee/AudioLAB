#pragma once

#define GAMMA_H_INC_ALL
#include "Gamma/Gamma.h"
#include "Gamma/Voices.h"


namespace Gamma::Oscillators
{

    struct Sweep : public GeneratorProcessorPlugin<gam::Sweep>
    {

    };

    struct WaveTable : public GeneratorProcessorPlugin<gam::Osc<DspFloatType>>
    {

    };

    struct CSine : public GeneratorProcessorPlugin<gam::CSine<DspFloatType>>
    {

    };

    struct Sine : public GeneratorProcessorPlugin<gam::Sine<DspFloatType>>
    {

    };

    struct SineR : public GeneratorProcessorPlugin<gam::SineR<DspFloatType>>
    {

    };

    struct SineRs : public GeneratorProcessorPlugin<gam::SineRs<DspFloatType>>
    {

    };

    struct SineD : public GeneratorProcessorPlugin<gam::SineD<DspFloatType>>
    {

    };

    struct SineDs : public GeneratorProcessorPlugin<gam::SineDs<DspFloatType>>
    {

    };

    struct Chirplet : public GeneratorProcessorPlugin<gam::Chriplet<DspFloatType>>
    {

    };

    struct LFO : public GeneratorProcessorPlugin<gam::LFO>
    {

    };

    struct DWO : public GeneratorProcessorPlugin<gam::DWO>
    {

    };

    struct Buzz : public GeneratorProcessorPlugin<gam::Buzz<DspFloatType>>
    {

    };

    struct Impulse : public GeneratorProcessorPlugin<gam::Impulse<DspFloatType>>
    {

    };

    struct Saw : public GeneratorProcessorPlugin<gam::Saw<DspFloatType>>
    {

    };

    struct Square : public GeneratorProcessorPlugin<gam::Square<DspFloatType>>
    {

    };

    struct DSF : public GeneratorProcessorPlugin<gam::DSF<DspFloatType>>
    {

    };

    struct Upsample : public GeneratorProcessorPlugin<gam::Upsample>
    {

    };

    struct Ramped : public GeneratorProcessorPlugin<gam::Ramped<DspFloatType>>
    {

    };
}
