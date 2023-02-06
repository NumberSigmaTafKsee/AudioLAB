
//------------------------------------------------------------------------------
// This file was generated using the Faust compiler (https://faust.grame.fr),
// and the Faust post-processor (https://github.com/jpcima/faustpp).
//
// Source: korg35lpf.dsp
// Name: Korg35LPF
// Author: Christopher Arndt
// Copyright: 
// License: MIT-style STK-4.3 license
// Version: 
//------------------------------------------------------------------------------

#pragma once

#include <algorithm>
#include <memory>

namespace Analog::Filters::Korg35
{
    struct ParameterRange {
        DspFloatType init;
        DspFloatType min;
        DspFloatType max;
    };

    class Korg35LPF {
    public:
        Korg35LPF();
        ~Korg35LPF();

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
            PORT_RESONANCE,
        };
        void setPort(int port, double v) {
            switch(port) {
                case PORT_CUTOFF: set_cutoff(v); break;
                case PORT_RESONANCE: set_q(v); break;
                default: printf("No port %d\n",port);
            }
        }
        double Tick(double input, double A=1, double X=1, double Y=1) {
            DspFloatType r = input;
            process(&r,&r,1);
            return r;
        }
        void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
            process(in,out,n);
        }
    public:
        class BasicDsp;

    private:
        std::unique_ptr<BasicDsp> fDsp;
    };    
}
#include "VAKorg35LPFFilter.cpp"



