#pragma once

#include "GenericSoundObject.hpp"
#include "FX/OnePole.hpp"
#include "FX/Filters.h"

namespace Analog::Oscillators
{
    ///////////////////////////////////////////////////////////////////////////////////////
    // The new blit
    ///////////////////////////////////////////////////////////////////////////////////////

    double BlitDSF(double phase, double m, double p, double a) 
    {
        double tmp, denominator = sin( phase );
        if ( std::fabs(denominator) <= std::numeric_limits<double>::epsilon() )
            tmp = a;
        else {
            tmp =  std::sin( m * phase );
            tmp /= p * denominator;
        }
        return tmp;
    }


///////////////////////////////////////////////////////////////////////////////////////
// blitSaw
///////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    struct Blit2SawOscillator : public GSSoundProcessor<T>
    {
        FX::Filters::OnePole block;
        unsigned int nHarmonics_;
        unsigned int m_;
        
        T rate_;
        T phase_;
        T offset;
        T p_;
        T C2_;
        T a_;
        T state_;
        T y;
        T sampleRate=44100;
        
        Blit2SawOscillator(T sampleRate=44100.0f, T frequency=440.0f)
        : GSSoundProcessor<T>()
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
        void setFrequency( T frequency )
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

        T getPhase() { 
            return phase_; 
        }

        void setPhaseOffset(T o) {
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
        void setPort(int port, T v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_HARMONICS: setHarmonics(v); break;
                //case PORT_GAIN: setGain(v); break;
                case PORT_PHASE: setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        //! Return the last computed output value.
        T lastOut( void ) const  { 
            return y; 
        };
        
        
        // blit = sin(m * phase) / (p * sin(phase));
        T Tick( T I=1, T A=1, T X=0, T Y=0 )
        {    
            T tmp = BlitDSF(phase_,m_,p_,a_);        
            tmp += state_ - C2_;
            state_ = tmp * 0.995;        
            phase_ += rate_;
            if ( phase_ >= M_PI ) phase_ -= M_PI;
            y = clamp(tmp,-1,1);        
            //y -= block.process(y);
            return 2*y;
        }

        
        void ProcessSIMD(size_t n, T * out) {
            
            #ifdef USE_GPU
            #pragma omp target teams distribute parallel for
            #else
            #pragma omp simd
            #endif
            for(size_t i = 0; i < n; i++)
            {
                T tmp = BlitDSF(phase_,m_,p_,a_);        
                tmp += state_ - C2_;
                state_ = tmp * 0.995;        
                phase_ += rate_;
                if ( phase_ >= M_PI ) phase_ -= M_PI;
                y = clamp(tmp,-1,1);                    
                out[i] = 2*y;
            }
        }
        T operator()() {
            return Tick();
        }                
    };

///////////////////////////////////////////////////////////////////////////////////////
// blitSquare
///////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    struct Blit2SquareOscillator : public GSSoundProcessor<T>
    {
        FX::Filters::OnePole block;
        unsigned int nHarmonics_;
        unsigned int m_;
        T f;
        T rate_;
        T phase_;
        T offset;
        T p_;
        T C2_;
        T a_;
        T state_;
        T y;
        T D;
        T sampleRate=44100;
        
        Blit2SquareOscillator(T sampleRate=44100.0f, T frequency=440.0f)
        : GSSoundProcessor<T>()
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
        void setFrequency( T frequency )
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
        void setDuty( T d) {
            D = d;
        }
        T getPhase() { 
            return phase_; 
        }

        void setPhaseOffset(T o) {
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
        void setPort(int port, T v) {
            switch(port) {
                case PORT_FREQ: setFrequency(v); break;
                case PORT_HARMONICS: setHarmonics(v); break;
                //case PORT_GAIN: setGain(v); break;
                case PORT_PHASE: setPhaseOffset(v); break;
                default: printf("No port %d\n", port);
            }
        }

        //! Return the last computed output value.
        T lastOut( void ) const  { 
            return y; 
        };
        
        
        // blit = sin(m * phase) / (p * sin(phase));
        T Tick( T I=1, T A=1, T X=0, T Y=0 )
        {    
            T tmp = BlitDSF(phase_,m_,p_,a_);        
            T tmp2= BlitDSF(phase_+D*M_PI,m_,p_,a_);
            tmp      = tmp - tmp2;
            //tmp     += state_ - C2_;        
            state_ += tmp * 0.995;
            phase_ += rate_;
            if ( phase_ >= 2*M_PI ) phase_ -= 2*M_PI;
            y = state_;            
            y = clamp(y,-1,1);
            return y;
        }
        void ProcessSIMD(size_t n, T * out)
        {
            #ifdef USE_GPU
            #pragma omp target teams distribute parallel for
            #else
            #pragma omp simd
            #endif
            for(size_t i = 0; i < n; i++)
            {
                T tmp = BlitDSF(phase_,m_,p_,a_);        
                T tmp2= BlitDSF(phase_+D*M_PI,m_,p_,a_);
                tmp      = tmp - tmp2;
                //tmp     += state_ - C2_;        
                state_ += tmp * 0.995;
                phase_ += rate_;
                if ( phase_ >= 2*M_PI ) phase_ -= 2*M_PI;
                y = state_;            
                y = clamp(y,-1,1);
                out[i] = y;
            }
        }
        T operator()() {
            return Tick();
        }
    };

///////////////////////////////////////////////////////////////////////////////////////
// blitTriangle
///////////////////////////////////////////////////////////////////////////////////////

    template<typename T>
    struct Blit2TriangleOscillator : public GSSoundProcessor<T>
    {
        Blit2SquareOscillator<T> sqr;
        FX::Filters::OnePole b1;
        T sampleRate=44100;

        T triangle;
        Blit2TriangleOscillator(T sampleRate=44100.0f, T frequency=440.0f) : 
        GSSoundProcessor<T>(),
        sqr(sampleRate,frequency)
        {
            this->sampleRate = sampleRate;
            b1.setFc(10.0f/sampleRate);
            triangle = 0;
        }
        void reset() {
            triangle = 0;
        }
        void setDuty(T d) {
            sqr.setDuty(d);
        }
        void setFrequency(T f) {
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
        void setPort(int port, T v) {
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

        T Tick(T I=1, T A=1, T X=1, T Y=1)
        {
            T x = sqr.Tick();
            T a = 1.0 - 0.01*std::fmin(1,sqr.f/1000.0);
            triangle = a*triangle + x/sqr.p_;            
            triangle -= b1.process(triangle);
            return 4*triangle;
        }
        void ProcessSIMD(size_t n, T * out)
        {
            #ifdef USE_GPU
            #pragma omp target teams distribue parallel for
            #else
            #pragma omp simd
            #endif
            for(size_t i = 0; i < n; i++)
            {
                T x = sqr.Tick();
                T a = 1.0 - 0.01*std::fmin(1,sqr.f/1000.0);
                triangle = a*triangle + x/sqr.p_;            
                triangle -= b1.process(triangle);
                out[i] = 4*triangle;
            }
        }    
        T operator()() {
            return Tick();
        }
    };
}
