#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    //////////////////////////////////////////////
    // Old Blit it works
    //////////////////////////////////////////////
    template<typename DSP>
    struct BlitSaw : public GSSoundProcessor<DSP>
    {
        //! Class constructor.
        BlitSaw( DSP frequency = 220.0, DSP sampleRate=44100 ) : GSSoundProcessor<DSP>()
        {
            nHarmonics_ = 0;
            offset = 0;
            this->sampleRate = sampleRate;
            reset();
            setFrequency( frequency );
            block.setFc(10.0f/sampleRate);
            gain = 1;
        }

        //! Class destructor.
        ~BlitSaw() = default;

        //! Resets the oscillator state and phase to 0.
        void reset()
        {
            phase_ = 0.0f;
            state_ = 0.0;
            y = 0.0;
        }

        //! Set the sawtooth oscillator rate in terms of a frequency in Hz.
        void setFrequency( DSP frequency )
        {
            f = frequency;
            p_ = sampleRate / frequency;
            C2_ = 1 / p_;
            rate_ = M_PI * C2_;
            updateHarmonics();
        }

        void setHarmonics( unsigned int nHarmonics = 0 )
        {
            nHarmonics_ = nHarmonics;
            this->updateHarmonics();        
            state_ = -0.5 * a_;
        }
        void setGain(DSP g) {
            gain = g;
        }
        void setPhaseOffset(DSP o) {
            phase_ = o;    
        }
        DSP getPhase() { return phase_; }
        void updateHarmonics( void )
        {            
            nHarmonics_ = (unsigned int) floor(  0.5 * p_ );
            m_ = 2 * nHarmonics_ + 1;        
            a_ = m_ / p_;
        }

        enum {
            PORT_FREQ,
            PORT_HARMONICS,
            PORT_GAIN,
            PORT_PHASE
        };
        void setPort(int port, DSP v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_HARMONICS: setHarmonics(v); break;
                case PORT_GAIN: setGain(v); break;
                case PORT_PHASE: setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        //! Return the last computed output value.
        DSP lastOut( void ) const { return y; };

        
        //! Compute and return one output sample.
        
        // blit = sin(m * phase) / (p * sin(phase));

        DSP Tick( DSP I=1, DSP A=1, DSP X=0, DSP Y=0 )
        {     
            // I = index
            // X = FM
            // Y = PM
            DSP tmp, denominator = sin( phase_ );
            if ( fabs(denominator) <= std::numeric_limits<DSP>::epsilon() )
                tmp = a_;
            else {
                tmp =  sin( m_ * phase_ );
                tmp /= p_ * denominator;
            }

            tmp += (state_ - C2_);
            state_ = tmp * 0.995;            
            phase_ += rate_;
            if ( phase_ >= M_PI ) phase_ -= M_PI;
            y = tmp;
            return 2*y;
        }

        void ProcessSIMD(size_t n, DSP * in, DSP * out) {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                // I = index
                // X = FM
                // Y = PM
                DSP tmp, denominator = sin( phase_ );
                if ( fabs(denominator) <= std::numeric_limits<DSP>::epsilon() )
                    tmp = a_;
                else {
                    tmp =  sin( m_ * phase_ );
                    tmp /= p_ * denominator;
                }

                tmp += (state_ - C2_);
                state_ = tmp * 0.995;            
                phase_ += rate_;
                if ( phase_ >= M_PI ) phase_ -= M_PI;
                y = tmp;
                out[i] = 2*y;     
                if(in) out[i] *= in[i];       
            }
        }
        void ProcessBlock(size_t n, DSP * in, DSP * out) {
            ProcessSIMD(n,in,out);
        }
        void ProcessInplace(size_t n, DSP * out) {
            ProcessSIMD(n,nullptr,out);
        }

        FX::Filters::OnePole     block;    
        unsigned int nHarmonics_;
        unsigned int m_;
        DSP f;
        DSP rate_;
        DSP phase_;
        DSP offset;
        DSP p_;
        DSP C2_;
        DSP a_;
        DSP state_;
        DSP y;
        DSP gain;
        DSP sampleRate=44100;
    };
}
