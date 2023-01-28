#include "math.h"
//#include "FxBase.h"
//#include "ModuleBase.h"

namespace WAVESHAPER {
	enum INVERT {OFF, ON};
	enum t
	{
		kAtanPos = 0,
		kAtanNeg,
		kStages,
		kInvertStages,
		kNumParams
	};
};

class WaveShaper
{
public:
	WaveShaper();
	~WaveShaper();

	void setSampleRate(float sampleRate);
	void setParameter(int index, float value);

	void resetBuffer();
	void resetCoeffs();

	void setAtanPos(float value);
	void setAtanNeg(float value);
	void setStages(float value);
	void setInvertStages(float value);

	void process(float &in);
    void process(float &inL, float &inR);
	void process(double &in);
    void process(double &inL, double &inR);

private:
	float sampleRate_;

	float arcTanPos;
	float arcTanNeg;

	unsigned int invertStages;
	float stages_;

	float *parameter_;
};

WaveShaper::WaveShaper()
{
	//name_ = "WaveShaper";

	parameter_ = new float[kNumParams];

	for(int i = 0; i < kNumParams; ++i)
	{
		parameter_[i] = 0;
	}

	resetBuffer();
	resetCoeffs();
}

WaveShaper::~WaveShaper()
{
	if(parameter_)
		delete[] parameter_;
	parameter_ = 0;
}

void WaveShaper::resetBuffer()
{
	
}

void WaveShaper::resetCoeffs()
{

}

void WaveShaper::setSampleRate(float sampleRate)
{
	sampleRate_ = sampleRate;
}

void WaveShaper::setAtanPos(float value)
{
	arcTanPos = 0.1f + value *20.f;
}

void WaveShaper::setAtanNeg(float value)
{
	arcTanNeg = 0.1f + value*20.f;
}

void WaveShaper::setStages(float value)
{
	stages_ = value;
}

void WaveShaper::setInvertStages(float value)
{
	if (value <= 0.5f)
			invertStages = ON;
		else 
			invertStages = OFF;
}

void WaveShaper::setParameter(int index, float value)
{
	switch(index)
	{
	case kAtanPos :
		parameter_[kAtanPos] = value;
		arcTanPos = 0.1f + parameter_[kAtanPos]*20.f;
		break;
	
	case kAtanNeg :
		parameter_[kAtanNeg] = value;
		arcTanNeg = 0.1f + parameter_[kAtanNeg]*20.f;
		break;

	case kStages :
		parameter_[kStages] = value;
		stages_ = parameter_[kStages];
		break;

	case kInvertStages :
		parameter_[kInvertStages] = value;
		if (parameter_[kInvertStages] <= 0.5f)
			invertStages = ON;
		else 
			invertStages = OFF;
		break;
	
	default : return;
	}
}


void WaveShaper::process(float &in)
{
	float f_xn = in;
	for(int i=0; i < stages_; i++)
	{
		if(f_xn >= 0.f)
			f_xn = (1.f/atan(arcTanPos))*atan(arcTanPos*f_xn);
		else
			f_xn = (1.f/atan(arcTanNeg))*atan(arcTanNeg*f_xn);

		if(invertStages == ON && i%2 == 0)
			f_xn *= -1.f;
	}

	in = f_xn;
}

void WaveShaper::process(float &inL, float &inR)
{
	float f_xnL = inL;
	for(int i=0; i < stages_; i++)
	{
		if(f_xnL >= 0.f)
			f_xnL = (1.f/atan(arcTanPos))*atan(arcTanPos*f_xnL);
		else
			f_xnL = (1.f/atan(arcTanNeg))*atan(arcTanNeg*f_xnL);

		if(invertStages == ON && i%2 == 0)
			f_xnL *= -1.f;
	}

	inL = f_xnL;
	inR = inL;
}
void WaveShaper::process(double &in)
{
	float f_xn = (float)in;
	for(int i=0; i < stages_; i++)
	{
		if(f_xn >= 0.f)
			f_xn = (1.f/atan(arcTanPos))*atan(arcTanPos*f_xn);
		else
			f_xn = (1.f/atan(arcTanNeg))*atan(arcTanNeg*f_xn);

		if(invertStages == ON && i%2 == 0)
			f_xn *= -1.f;
	}

	in = (double)f_xn;
}

void WaveShaper::process(double &inL, double &inR)
{
	float f_xnL = (float)inL;
	for(int i=0; i < stages_; i++)
	{
		if(f_xnL >= 0.f)
			f_xnL = (1.f/atan(arcTanPos))*atan(arcTanPos*f_xnL);
		else
			f_xnL = (1.f/atan(arcTanNeg))*atan(arcTanNeg*f_xnL);

		if(invertStages == ON && i%2 == 0)
			f_xnL *= -1.f;
	}

	inL =(double) f_xnL;
	inR =(double) inL;
}