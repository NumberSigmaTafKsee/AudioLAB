
//---------------------------------------------------------------------------
//
//                                3 Band EQ :)
//
// EQ.H - Header file for 3 band EQ
//
// (c) Neil C / Etanza Systems / 2K6
//
// Shouts / Loves / Moans = etanza at lycos dot co dot uk
//
// This work is hereby placed in the public domain for all purposes, including
// use in commercial applications.
//
// The author assumes NO RESPONSIBILITY for any problems caused by the use of
// this software.
//
//----------------------------------------------------------------------------

#pragma once 
#include <cmath>
#include "Undenormal.hpp"
#include "GenericSoundObject.hpp"

// ------------
//| Structures |
// ------------

template<typename T>
struct ToneStack : public GSSoundProcessor<T>
{
    // Filter #1 (Low band)

    T  lf;       // Frequency
    T  f1p0;     // Poles ...
    T  f1p1;
    T  f1p2;
    T  f1p3;

    // Filter #2 (High band)

    T  hf;       // Frequency
    T  f2p0;     // Poles ...
    T  f2p1;
    T  f2p2;
    T  f2p3;

    // Sample history buffer

    T  sdm1;     // Sample data minus 1
    T  sdm2;     //                   2
    T  sdm3;     //                   3

    // Gain Controls

    T  lg;       // low  gain
    T  mg;       // mid  gain
    T  hg;       // high gain

    T lowfreq,highfreq,mixfreq,input_gain,output_gain;
    
    ToneStack(T lowfreq=880.0f,T highfreq=5000.0f, T mixfreq=44100.0f) 
    : GSSoundProcessor<T>()
    {
        lg = 1.0;
        mg = 1.0;
        hg = 1.0;
        this->lowfreq = lowfreq;
        this->highfreq = highfreq;
        this->mixfreq = mixfreq;
        input_gain = output_gain = 1.0f;
        // Calculate filter cutoff frequencies

        lf = 2 * std::sin(M_PI * ((T)lowfreq / (T)mixfreq));
        hf = 2 * std::sin(M_PI * ((T)highfreq / (T)mixfreq));
    }
    T Tick(T sample, T A = 1, T X = 0, T Y = 0)
    {
            T  l,m,h;      // Low / Mid / High - Sample Values
            Undenormal denormal;

            
            T tl = lowfreq;
            T th = highfreq;

            lf = 2 * std::sin(M_PI * ((T)(lowfreq + X * lowfreq)/ (T)mixfreq));
            hf = 2 * std::sin(M_PI * ((T)(highfreq + Y * highfreq)/ (T)mixfreq));

            // Filter #1 (lowpass)
            sample *= input_gain;

            f1p0  += (lf * (sample - f1p0));
            f1p1  += (lf * (f1p0 - f1p1));
            f1p2  += (lf * (f1p1 - f1p2));
            f1p3  += (lf * (f1p2 - f1p3));

            l          = f1p3;

            // Filter #2 (highpass)
            f2p0  += (hf * (sample - f2p0));
            f2p1  += (hf * (f2p0 - f2p1));
            f2p2  += (hf * (f2p1 - f2p2));
            f2p3  += (hf * (f2p2 - f2p3));

            h          = sdm3 - f2p3;

            // Calculate midrange (signal - (low + high))
            m          = sdm3 - (h + l);

            // Scale, Combine and store

            l         *= lg;
            m         *= mg;
            h         *= hg;

            // Shuffle history buffer

            sdm3   = sdm2;
            sdm2   = sdm1;
            sdm1   = output_gain*sample;

            lf = tl;
            hf = th;

            // Return result
            return A*(l + m + h);
    }
	void ProcessSIMD(size_t n, T * in, T * out) {
		#pragma omp simd
		for(size_t i = 0; i < n; i++) out[i] = Tick(in[i]);
	}
	void ProcessBlock(size_t n, T * in, T * out) {
		ProcessSIMD(n,in,out);
	}
	void ProcessInplace(size_t n, T * buffer) {
		ProcessSIMD(n,buffer,buffer);
	}
};




