
//------------------------------------------------------------------------------
// This file was generated using the Faust compiler (https://faust.grame.fr),
// and the Faust post-processor (https://github.com/jpcima/faustpp).
//
// Source: diodeladder.dsp
// Name: DiodeLadder
// Author: Christopher Arndt
// Copyright: 
// License: MIT-style STK-4.3 license
// Version: 
//------------------------------------------------------------------------------
#pragma once

#include <memory>
#include <vector>
#include "DSP/Decimators.hpp"

namespace Analog::Filters::DiodeLadder
{
    struct ParameterRange {
        DspFloatType init;
        DspFloatType min;
        DspFloatType max;
    };

    class DiodeLadder 
    {
    public:
        DiodeLadder();
        ~DiodeLadder();

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
        
        enum {
            PORT_CUTOFF,
            PORT_Q,
        };
        void setPort(int port, double v) {
            switch(port) {
                case PORT_CUTOFF: set_cutoff(v); break;
                case PORT_Q: set_q(v); break;
                default: printf("No port %d\n", port);
            }
        }
        DspFloatType Tick(DspFloatType in) {
            DspFloatType r = 0.0;
            process(&in,&r,1);
            return r;
        }
        void Process(size_t n, DspFloatType * input, DspFloatType * output) {
            process(input,output,n);
        }


    public:
        class BasicDsp;

    private:
        std::unique_ptr<BasicDsp> fDsp;

    };

}
#include "VADiodeLadderFilter.cpp"