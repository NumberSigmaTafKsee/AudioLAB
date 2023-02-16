#pragma once

#include <string>

#define	MAXBUFFERSIZE 8192		// must be about 1/5 of a second at given sample rate
#define ROUND(n)		((int)((double)(n)+0.5))

#define	MODF(n,i,f) ((i) = (int)(n), (f) = (n) - (double)(i))
//	modf - vaguely related to the library routine modf(), this macro breaks a double into
//	integer and fractional components i and f respectively.
//	n - input number, a double
//	i - integer portion, an integer (the input number integer portion should fit)
//	f - fractional portion, a double


#define	NUM_MIX_MODES	7
#define	NUM_DELAYS	11

namespace Fx::CHORUS
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



	class Chorus : public StereoFXProcessor
	{
	public:
		Chorus();
		~Chorus();
		
		void resetBuffer();
		void resetCoeffs();
		
		void setParameter(int index, double value);
		void setSampleRate (double sampleRate);
		
		void setRate ();
		void setWidth ();
		void setDelay ();
		void setWetLevel ();
		void setSweep();
		
		void process(double &in);
		void process(double &inL, double &inR);
		void process(double &in);
		void process(double &inL, double &inR);

		double Tick(double iL, double iR, double & L, double &R, double A=1, double X=0, double Y=0)
		{
			double ps = paramSweepRate;
			double pw = paramWidth;
			double pd = paramDelay;
			double pwet = paramWetLevel;

			paramSweepRate = paramSweepRate + (X+Y);
			paramWidth = paramWidth + X;
			paramDelay = paramDelay + Y;
			paramWetLevel = paramWetLevel*A;
			setParam
			L = iL;
			R = iR;
			process(L,R);
			paramSweepRate = ps;
			paramWidth = pw;
			paramDelay = pd;
			paramWetLevel = pwet;
			return 0.5*(L+R);
		}
		void ProcessBlock(size_t n, double ** in, double ** out) {
			for(size_t i = 0; i < n; i++) {
				out[0][i] = in[0][i];
				out[1][i] = in[1][i];
				process(out[0][i],out[1][i]);
			}
		}

	private :
		double *parameter_;
		double paramSweepRate;
		double paramWidth;
		double paramDelay;
		double paramWetLevel;	
		
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
		double sampleRate_;
		std::string name_;
	};
}