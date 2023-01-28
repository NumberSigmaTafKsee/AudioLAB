
//------------------------------------------------------------------------------
// This file was generated using the Faust compiler (https://faust.grame.fr),
// and the Faust post-processor (https://github.com/jpcima/faustpp).
//
// Source: oberheim.dsp
// Name: Oberheim
// Author: Christopher Arndt
// Copyright: 
// License: MIT-style STK-4.3 license
// Version: 
//------------------------------------------------------------------------------

#pragma once

#include <memory>

namespace Analog::Filters::Oberheim
{
    struct ParameterRange {
            DspFloatType init;
            DspFloatType min;
            DspFloatType max;
        };    
    class Oberheim {
    public:
        Oberheim();
        ~Oberheim();

        void init(DspFloatType sample_rate);
        void clear() noexcept;

        void process(
            const DspFloatType *in0,
            DspFloatType *out0,DspFloatType *out1,DspFloatType *out2,DspFloatType *out3,
            unsigned count) noexcept;

        enum { NumInputs = 1 };
        enum { NumOutputs = 4 };
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
        
        
        enum filterOut {
            OUT_1,
            OUT_2,
            OUT_3,
            OUT_4,
        } filterout = OUT_1;
        enum {
            PORT_CUTOFF,
            PORT_Q,
            PORT_OUTPUT,
        };
        void setPort(int port, double v) {
            switch(port) {
                case PORT_CUTOFF: set_cutoff(v); break;
                case PORT_Q: set_q(v); break;
                case PORT_OUTPUT: filterout = (filterOut)v; break;
            }
        }
        DspFloatType Tick(DspFloatType input) {
            DspFloatType out1,out2,out3,out4;
            process(&input,&out1,&out2,&out3,&out4,1);
            switch(filterout) {
                case OUT_1: return out1;
                case OUT_2: return out2;
                case OUT_3: return out3;
                case OUT_4: return out4;
            }
            return 0;
        }
        void Process(size_t n, DspFloatType * input, DspFloatType * output) {
            for(size_t i = 0; i < n; i++) output[i] = Tick(input[i]);
        }
        void Process(DspFloatType * samples, size_t n ) {
            for(size_t i = 0; i < n; i++) samples[i] = Tick(samples[i]);
        }

    public:
        class BasicDsp;

    private:
        std::unique_ptr<BasicDsp> fDsp;



    };
}

#include "VAOberheimFilter.cpp"