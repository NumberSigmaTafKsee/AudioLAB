#include <math.h>
#include "Chorus.h"

using namespace CHORUS;
	
Chorus::Chorus()
{
	name_ = "Chorus";

	parameter_ = new float[kNumParams];

	for(int i = 0; i < kNumParams; ++i)
	{
		parameter_[i] = 0;
	}

	reset();
	//constructor
	paramSweepRate = 0.2f;
	paramWidth = 0.3f;
	paramDelay = 0.2f;
	paramWetLevel = 0.5f;
	
	_sweepRate = 0.2;
	_sweepSamples = 0;
	_delaySamples = 10;
	writeIndex = 0;
	_sweep = 0.0;

	// allocate the buffer
	_bufferL = new double[MAXBUFFERSIZE];
	_bufferR = new double[MAXBUFFERSIZE];
}

Chorus::~Chorus ()
{
	//destructor
	if(parameter_)
		delete[] parameter_;
	parameter_ = 0;
}

void Chorus::resetBuffer()
{

}

void Chorus::resetCoeffs()
{

}

void Chorus::setSampleRate(float sampleRate)
{
	sampleRate_ = sampleRate;
}

void Chorus::setParameter(int index, float value)
{
	switch(index)
	{
	case kRate :
		parameter_[kRate] = value;
		paramSweepRate = parameter_[kRate];
		break;
	
	case kWidth :
		parameter_[kWidth] = value;
		paramWidth = parameter_[kWidth];
		break;

	case kDelay :
		parameter_[kDelay] = value;
		paramDelay = parameter_[kDelay];
		break;

	case kWetLevel :
		parameter_[kWetLevel] = value;
		paramWetLevel = parameter_[kWetLevel];
		break;
	}
}

void Chorus::setRate ()
{
	// map into param onto desired sweep range with log curve
	_sweepRate = std::pow(10.0,(double)paramSweepRate);
	_sweepRate  -= 0.1f;
	_sweepRate  *= 1.1f;
	_sweepRate  += 0.1f;

	// finish setup
	setSweep();
}

void Chorus::setWidth ()
{
	// map so that we can spec between 0ms and 50ms
	_sweepSamples = std::round(paramWidth * 0.05 * sampleRate_);

	// finish setup
	setSweep();
}

void Chorus::setDelay ()
{	
	double delay = pow(10.0, (double)paramDelay * 2.0)/1000.0;		// map logarithmically and convert to seconds
	_delaySamples = std::round(delay * sampleRate_);

	// finish setup
	setSweep();
}

void Chorus::setWetLevel()
{
	
}

void Chorus::setSweep()
{
	// calc # of samples per second we'll need to move to achieve spec'd sweep rate
	_step = (double)(_sweepSamples * 2.0 * _sweepRate) / sampleRate_;
	
	if( _step <= 1.0 ) // NEED CALC
		_step = 1.0;	// NEED CALC
		
	// calc min and max sweep now
	_minSweepSamples = _delaySamples;
	_maxSweepSamples = _delaySamples + _sweepSamples;

	// set intial sweep pointer to midrange
	_sweep = (_minSweepSamples + _maxSweepSamples) / 2;
}

void Chorus::process(float &inL, float &inR, float & outL, float & outR)
{
		// assemble input value and store it in circle queue
		_bufferL[writeIndex] = inL;
		_bufferR[writeIndex] = inR;
		
		writeIndex = (writeIndex + 1) & (MAXBUFFERSIZE-1);

		// build the two emptying pointers and do linear interpolation
		int ep1, ep2;
		double w1, w2;
		double ep = writeIndex - _sweep;
		
		w2 = std::fmod(ep, ep1);
		
		ep1 &= (MAXBUFFERSIZE-1);
		ep2 = ep1 + 1;
		ep2 &= (MAXBUFFERSIZE-1);
		
		w1 = 1.0 - w2;
		
		outL = _bufferL[ep1] * w1 + _bufferL[ep2] * w2;
		outR = _bufferL[ep1] *w1 + _bufferR[ep2] * w2;

		// develop output wet
		inL = (float)(paramWetLevel * outL + (1 - paramWetLevel) * inL);
		inR = (float)(paramWetLevel * outR + (1 - paramWetLevel) * inR);

		// increment the sweep
		_sweep += _step;
		if( _sweep >= _maxSweepSamples || _sweep <= _minSweepSamples )
		{
			_step = -_step;
		}

}

void Chorus::process(float &in, float & outL, float & outR)
{
		// assemble input value and store it in circle queue
		_bufferL[writeIndex] = in;
		
		writeIndex = (writeIndex + 1) & (MAXBUFFERSIZE-1);

		// build the two emptying pointers and do linear interpolation
		int ep1, ep2;
		double w1, w2;
		double ep = writeIndex - _sweep;
		
		w2 = std::fmod(ep, ep1);
		
		ep1 &= (MAXBUFFERSIZE-1);
		ep2 = ep1 + 1;
		ep2 &= (MAXBUFFERSIZE-1);
		
		w1 = 1.0 - w2;
		
		outL = _bufferL[ep1] * w1 + _bufferL[ep2] * w2;
		outR = _bufferL[ep1] *w1 + _bufferR[ep2] * w2;

		// develop output wet
		in = (float)(paramWetLevel * outL + (1 - paramWetLevel) * in);

		// increment the sweep
		_sweep += _step;
		if( _sweep >= _maxSweepSamples || _sweep <= _minSweepSamples )
		{
			_step = -_step;
		}
}

void Chorus::process(double &in, double & outL, double & outR)
{
		// assemble input value and store it in circle queue
		_bufferL[writeIndex] = in;
		
		writeIndex = (writeIndex + 1) & (MAXBUFFERSIZE-1);

		// build the two emptying pointers and do linear interpolation
		int ep1, ep2;
		double w1, w2;
		double ep = writeIndex - _sweep;
		
		w2 = std::fmod(ep, ep1);
		
		ep1 &= (MAXBUFFERSIZE-1);
		ep2 = ep1 + 1;
		ep2 &= (MAXBUFFERSIZE-1);
		
		w1 = 1.0 - w2;
		
		outL = _bufferL[ep1] * w1 + _bufferL[ep2] * w2;
		outR = _bufferL[ep1] *w1 + _bufferR[ep2] * w2;

		// develop output wet
		in = paramWetLevel * outL + (1 - paramWetLevel) * in;

		// increment the sweep
		_sweep += _step;
		if( _sweep >= _maxSweepSamples || _sweep <= _minSweepSamples )
		{
			_step = -_step;
		}
}
void Chorus::process(double &inL, double &inR, double & outL, double & outR)
{
		// assemble input value and store it in circle queue
		_bufferL[writeIndex] = inL;
		_bufferR[writeIndex] = inR;
		
		writeIndex = (writeIndex + 1) & (MAXBUFFERSIZE-1);

		// build the two emptying pointers and do linear interpolation
		int ep1, ep2;
		double w1, w2;
		double ep = writeIndex - _sweep;
		
		w2 = std::fmod(ep, ep1);
		
		ep1 &= (MAXBUFFERSIZE-1);
		ep2 = ep1 + 1;
		ep2 &= (MAXBUFFERSIZE-1);
		
		w1 = 1.0 - w2;
		
		outL = _bufferL[ep1] * w1 + _bufferL[ep2] * w2;
		outR = _bufferL[ep1] *w1 + _bufferR[ep2] * w2;

		// develop output wet
		inL = paramWetLevel * outL + (1 - paramWetLevel) * inL;
		inR = paramWetLevel * outR + (1 - paramWetLevel) * inR;

		// increment the sweep
		_sweep += _step;
		if( _sweep >= _maxSweepSamples || _sweep <= _minSweepSamples )
		{
			_step = -_step;
		}

}