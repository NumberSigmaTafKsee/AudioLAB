
//------------------------------------------------------------------------------
// This file was generated using the Faust compiler (https://faust.grame.fr),
// and the Faust post-processor (https://github.com/jpcima/faustpp).
//
// Source: mooghalfladder.dsp
// Name: MoogHalfLadder
// Author: Eric Tarr
// Copyright: 
// License: MIT-style STK-4.3 license
// Version: 
//------------------------------------------------------------------------------

#pragma once

#include <memory>

namespace Filters::Moog
{
    struct ParameterRange {
        DspFloatType init;
        DspFloatType min;
        DspFloatType max;
    };

    class MoogHalfLadder {
    public:
        MoogHalfLadder();
        ~MoogHalfLadder();

        void init(DspFloatType sample_rate);
        void clear() noexcept;

        void process(
            const DspFloatType *in0,
            DspFloatType *out0,
            unsigned count) noexcept;

        enum { NumInputs = 1 };
        enum { NumOutputs = 1 };
        enum { NumActives = 2 };
        enum { NumPassives = 0 };
        enum { NumParameters = 2 };

        enum Parameter {
            p_cutoff,
            p_q,
            
        };


        static const char *parameter_label(unsigned index) noexcept;
        static const char *parameter_short_label(unsigned index) noexcept;
        static const char *parameter_symbol(unsigned index) noexcept;
        static const char *parameter_unit(unsigned index) noexcept;
        static const ParameterRange *parameter_range(unsigned index) noexcept;
        static bool parameter_is_trigger(unsigned index) noexcept;
        static bool parameter_is_boolean(unsigned index) noexcept;
        static bool parameter_is_integer(unsigned index) noexcept;
        static bool parameter_is_logarithmic(unsigned index) noexcept;

        DspFloatType get_parameter(unsigned index) const noexcept;
        void set_parameter(unsigned index, DspFloatType value) noexcept;

        
        DspFloatType get_cutoff() const noexcept;
        
        DspFloatType get_q() const noexcept;
        
        
        void set_cutoff(DspFloatType value) noexcept;
        
        void set_q(DspFloatType value) noexcept;
        
        enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
        };
        void setPort(int port, double v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                set_cutoff(v);
                break;
            case PORT_RESONANCE:
                set_q(v);
                break;		
            }
        }
        void Process(size_t n, DspFloatType * input, DspFloatType * output) {
            process(input,output,n);
        }
        void Process(DspFloatType * input, size_t n ) {
            process(input,input,n);
        }
        DspFloatType Tick(DspFloatType input) {
            DspFloatType r = 0.0;
            Process(1,&input,&r);
            return r;	
        }
    public:
        class BasicDsp;

    private:
        std::unique_ptr<BasicDsp> fDsp;



    };
}

#include "VAMoogHalfLadderFilter.cpp"