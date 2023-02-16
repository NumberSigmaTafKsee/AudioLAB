#pragma once 

#include "chorus_1.h"

namespace SoundWave
{
    enum t
    {
        // Parameters Tags
        kRate = 0,
        kWidth,
        kDelay,
        kWetLevel,

        kNumParams
    };

    class Chorus 
    {
    public:
        Chorus();
        ~Chorus();
        
        void resetBuffer();
        void resetCoeffs();
        
        void setParameter(int index, float value);
        void setSampleRate (float sampleRate);
        
        void setRate ();
        void setWidth ();
        void setDelay ();
        void setWetLevel ();
        void setSweep();
        
        void process(float &in, float & outL, float & outR);
        void process(float &inL, float &inR, float & outL, float & outR);
        void process(double &in, double & outL, double & outR);
        void process(double &inL, double &inR, double & outL, double & outR);

        void Process(size_t n, float * in, float * out);
        void Process(size_t n, float * inL, float * inR, float * outL, float * outR);
        void Process(size_t n, float * in_stereo, float * out_stereo);

        void Process(size_t n, double * in, double * out);
        void Process(size_t n, double * inL, double * inR, double * outL, double * outR);
        void Process(size_t n, double * in_stereo, double * out_stereo);

    private :
        float *parameter_;
        float paramSweepRate;
        float paramWidth;
        float paramDelay;
        float paramWetLevel;	
        
        double _minSweepSamples;	// lower bound of calculated sweep range, calc'd by setSweep from rate, width, delay
        double _maxSweepSamples;	// upper bound, ditto
        
        double _sweepRate;			// actual calculated sweep rate
        int		_sweepSamples;			// sweep width in # of samples
        int	   _delaySamples;		// number of samples to run behind filling pointer

        double *_bufferL;		// stored sound
        double *_bufferR;
        int	   writeIndex;					// fill/write pointer
        double _step;				// amount to step the sweep each sample
        double _sweep;
    };

}

