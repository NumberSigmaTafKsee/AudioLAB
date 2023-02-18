#pragma once

#include "SoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"


namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // The new blit
    ///////////////////////////////////////////////////////////////////////////////////////
    DspFloatType BlitDSF(DspFloatType phase,DspFloatType m,DspFloatType p, DspFloatType a) 
    {
        DspFloatType tmp, denominator = sin( phase );
        if ( fabs(denominator) <= std::numeric_limits<DspFloatType>::epsilon() )
            tmp = a;
        else {
            tmp =  sin( m * phase );
            tmp /= p * denominator;
        }
        return tmp;
    }

    struct blitSaw : public OscillatorProcessor
    {
        FX::Filters::OnePole block;
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
        DspFloatType sampleRate=44100;
        
        blitSaw(DspFloatType sampleRate=44100.0f, DspFloatType frequency=440.0f)
        : OscillatorProcessor()
        {
            this->sampleRate = sampleRate;
            nHarmonics_ = 0;
            offset = 0;
            reset();
            setFrequency( frequency );
            block.setFc(10.0f/sampleRate);        
        }

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
            p_      = (sampleRate) / frequency;
            C2_     = 1 / p_;
            rate_   = M_PI * C2_;
            updateHarmonics();
        }

        void setHarmonics( unsigned int nHarmonics = 0 )
        {
            nHarmonics_ = nHarmonics;
            this->updateHarmonics();        
            state_ = -0.5 * a_;
        }

        DspFloatType getPhase() { 
            return phase_; 
        }

        void setPhaseOffset(DspFloatType o) {
            phase_ = o;    
        }

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
            //PORT_GAIN,
            PORT_PHASE
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_HARMONICS: setHarmonics(v); break;
                //case PORT_GAIN: setGain(v); break;
                case PORT_PHASE: setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        //! Return the last computed output value.
        DspFloatType lastOut( void ) const  { 
            return y; 
        };
        
        
        // blit = sin(m * phase) / (p * sin(phase));
        DspFloatType Tick( DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0 )
        {    
            DspFloatType tmp = BlitDSF(phase_,m_,p_,a_);        
            tmp += state_ - C2_;
            state_ = tmp * 0.995;        
            phase_ += rate_;
            if ( phase_ >= M_PI ) phase_ -= M_PI;
            y = clamp(tmp,-1,1);        
            y -= block.process(y);
            return 1.9*y;
        }                

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType tmp = BlitDSF(phase_,m_,p_,a_);        
                tmp += state_ - C2_;
                state_ = tmp * 0.995;        
                phase_ += rate_;
                if ( phase_ >= M_PI ) phase_ -= M_PI;
                y = clamp(tmp,-1,1);        
                y -= block.process(y);
                out[i] = 1.9*y;
                if(in) out[i] *= in[i];
            }
        }

        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
        
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,input,input);
        }

    };

    struct blitSquare : public OscillatorProcessor
    {
        FX::Filters::OnePole block;
        unsigned int nHarmonics_;
        unsigned int m_;
        DspFloatType f;
        DspFloatType rate_;
        DspFloatType phase_;
        DspFloatType offset;
        DspFloatType p_;
        DspFloatType C2_;
        DspFloatType a_;
        DspFloatType state_;
        DspFloatType y;
        DspFloatType D;
        DspFloatType sampleRate=44100;
        
        blitSquare(DspFloatType sampleRate=44100.0f, DspFloatType frequency=440.0f)
        : OscillatorProcessor()
        {
            this->sampleRate = sampleRate;
            nHarmonics_ = 0;
            offset = 0;
            reset();
            setFrequency( frequency );
            block.setFc(10.0f/sampleRate);        
            D = 0.5;
        }

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
            f       = frequency;
            p_      = (sampleRate) / frequency;
            C2_     = 1 / p_;
            rate_   = M_PI * C2_;
            updateHarmonics();
        }

        void setHarmonics( unsigned int nHarmonics = 0 )
        {
            nHarmonics_ = nHarmonics;
            this->updateHarmonics();        
            state_ = -0.5 * a_;
        }
        void setDuty( DspFloatType d) {
            D = d;
        }
        DspFloatType getPhase() { 
            return phase_; 
        }

        void setPhaseOffset(DspFloatType o) {
            phase_ = o;    
        }

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
            //PORT_GAIN,
            PORT_PHASE
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_HARMONICS: setHarmonics(v); break;
                //case PORT_GAIN: setGain(v); break;
                case PORT_PHASE: setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        //! Return the last computed output value.
        DspFloatType lastOut( void ) const  { 
            return y; 
        };
        
        
        // blit = sin(m * phase) / (p * sin(phase));
        DspFloatType Tick( DspFloatType I=1, DspFloatType A=1, DspFloatType X=0, DspFloatType Y=0 )
        {    
            DspFloatType tmp = BlitDSF(phase_,m_,p_,a_);        
            DspFloatType tmp2= BlitDSF(phase_+D*M_PI,m_,p_,a_);
            tmp      = tmp - tmp2;
            //tmp     += state_ - C2_;        
            state_ += tmp * 0.995;
            phase_ += rate_;
            if ( phase_ >= 2*M_PI ) phase_ -= 2*M_PI;
            y = state_;
            y -= block.process(y);
            return y;
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType tmp = BlitDSF(phase_,m_,p_,a_);        
                DspFloatType tmp2= BlitDSF(phase_+D*M_PI,m_,p_,a_);
                tmp      = tmp - tmp2;
                //tmp     += state_ - C2_;        
                state_ += tmp * 0.995;
                phase_ += rate_;
                if ( phase_ >= 2*M_PI ) phase_ -= 2*M_PI;
                y = state_;
                y -= block.process(y);
                out[i] = y;
                if(in) out[i] *= in[i];
            }
        }

        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
        
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,nullptr,input);
        }

    };

    struct blitTriangle : public OscillatorProcessor
    {
        blitSquare sqr;
        FX::Filters::OnePole b1;
        DspFloatType sampleRate=44100;

        DspFloatType triangle;
        blitTriangle(DspFloatType sampleRate=44100.0f, DspFloatType frequency=440.0f) : 
        OscillatorProcessor(),
        sqr(sampleRate,frequency)
        {
            this->sampleRate = sampleRate;
            b1.setFc(10.0f/sampleRate);
            triangle = 0;
        }
        void reset() {
            triangle = 0;
        }
        void setDuty(DspFloatType d) {
            sqr.setDuty(d);
        }
        void setFrequency(DspFloatType f) {
            sqr.setFrequency(f);
        }
        enum {
            PORT_FREQ,
            PORT_DUTY,
            PORT_RESET,
            PORT_HARMONICS,
            //PORT_GAIN,
            PORT_PHASE
        };
        void setPort(int port, DspFloatType v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_DUTY: setDuty(v); break;
                case PORT_RESET: reset(); break;
                case PORT_HARMONICS: sqr.setHarmonics(v); break;
                //case PORT_GAIN: sqr.setGain(v); break;
                case PORT_PHASE: sqr.setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        DspFloatType Tick(DspFloatType I=1, DspFloatType A=1, DspFloatType X=1, DspFloatType Y=1)
        {
            DspFloatType x = sqr.Tick();
            DspFloatType a = 1.0 - 0.1*std::fmin(1,sqr.f/1000.0);
            triangle = a*triangle + x/sqr.p_;
            DspFloatType kaka = b1.process(triangle);
            triangle -= kaka;
            return 4*triangle;
        }

        void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out)
        {
            #pragma omp simd aligned(in,out)
            for(size_t i = 0; i < n; i++)
            {
                DspFloatType x = sqr.Tick();
                DspFloatType a = 1.0 - 0.1*std::fmin(1,sqr.f/1000.0);
                triangle = a*triangle + x/sqr.p_;
                DspFloatType kaka = b1.process(triangle);
                triangle -= kaka;                
                out[i] = 4*triangle;
                if(in) out[i] *= in[i];
            }
        }

        void ProcessBlock(size_t n, DspFloatType * input, DspFloatType * output) {
            ProcessSIMD(n,input,output);
        }
        
        void ProcessInplace(size_t n, DspFloatType * input) {
            ProcessBlock(n,input,input);
        }

    };
}
