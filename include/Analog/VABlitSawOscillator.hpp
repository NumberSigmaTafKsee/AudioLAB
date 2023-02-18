#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    //////////////////////////////////////////////
    // Old Blit it works
    //////////////////////////////////////////////
    struct BlitSaw : public OscillatorProcessor
    {
        //! Class constructor.
        BlitSaw( DspFloatType frequency = 220.0, DspFloatType sampleRate=44100 ) : OscillatorProcessor()
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
        void setFrequency( DspFloatType frequency )
        {
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
        void setGain(DspFloatType g) {
            gain = g;
        }
        void setPhaseOffset(DspFloatType o) {
            phase_ = o;    
        }
        DspFloatType getPhase() { return phase_; }
        void updateHarmonics( void )
        {
            if ( nHarmonics_ <= 0 ) {
                unsigned int maxHarmonics = (unsigned int) floor( 0.5 * p_ );
                m_ = 2 * maxHarmonics + 1;
            }
            else
                m_ = 2 * nHarmonics_ + 1;

            a_ = m_ / p_;
        }

        enum {
            PORT_FREQ,
            PORT_HARMONICS,
            PORT_GAIN,
            PORT_PHASE
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_HARMONICS: setHarmonics(v); break;
                case PORT_GAIN: setGain(v); break;
                case PORT_PHASE: setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        //! Return the last computed output value.
        DspFloatType lastOut( void ) const { return y; };

        
        //! Compute and return one output sample.
        
        // blit = sin(m * phase) / (p * sin(phase));

        DspFloatType Tick( DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0 )
        {     
            // I = index
            // X = FM
            // Y = PM
            DspFloatType tmp, denominator = sin( phase_ );
            if ( fabs(denominator) <= std::numeric_limits<DspFloatType>::epsilon() )
                tmp = a_;
            else {
                tmp =  sin( m_ * phase_ );
                tmp /= p_ * denominator;
            }

            tmp += state_ - C2_;
            state_ = tmp * 0.995;            
            phase_ += rate_;
            if ( phase_ >= M_PI ) phase_ -= M_PI;
            y = tmp - block.process(tmp);
            return 1.8*y;
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType tmp, denominator = sin( phase_ );
                if ( fabs(denominator) <= std::numeric_limits<DspFloatType>::epsilon() )
                    tmp = a_;
                else {
                    tmp =  sin( m_ * phase_ );
                    tmp /= p_ * denominator;
                }

                tmp += state_ - C2_;
                state_ = tmp * 0.995;            
                phase_ += rate_;
                if ( phase_ >= M_PI ) phase_ -= M_PI;
                y = tmp - block.process(tmp);
                out[i] = 1.8*y;
                if(in) out[i] *= in[i];
            }            
        }
        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
            
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,nullptr,input);
        }

        FX::Filters::OnePole     block;    
        unsigned int nHarmonics_;
        unsigned int m_;
        DspFloatType rate_;
        DspFloatType phase_;
        DspFloatType offset;
        DspFloatType p_;
        DspFloatType C2_;
        DspFloatType a_;
        DspFloatType state_;
        DspFloatType y;
        DspFloatType gain;
        DspFloatType sampleRate=44100;
    };
}
