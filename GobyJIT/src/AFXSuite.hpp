#pragma once

#include <iostream>
#include <cmath>
#include <cstdint>

//==============================================================================
/**
	A Base class for delay based digital effects. Provides the basic methods
	that are shared amongst Flanger, Delay, Chorus and Phaser
    
    @version 0.1
 
 @see DelayEffectBase
 */
class DelayEffectBase
{
public: // Methods
	//==============================================================================
	/** Constructor. */
    DelayEffectBase();
    DelayEffectBase(int bufferSizeSamples);

	/** Destructor. */
    ~DelayEffectBase();
	
	//==============================================================================
	/** Main process block for applying audio effect
		@param inputSample The input audio sample for the effect to be applied to
		
		@returns an audio sample as a double with effect applied
	 */
	virtual double process(double inputSample) = 0;

    /**
     <#Description#>

     @param bufferSizeSamples <#bufferSizeSamples description#>
     */
    void setupDelayEffectBase(const int bufferSizeSamples);
private:	//Methods
	//==============================================================================
	/** Sets the internal lagrange interpolation table. Ideally it should be shared
	 amongst all
	 
	 @returns    false and an error description if there's a problem,
	 or sets the interpolation lookup table and returns true
	 */
	static double** setInterpolationTable();
	
	/***/
	void printInterpTable();
protected:
	//==============================================================================
	/**this should be informed by what kind off effect will be implemented*/
	//	virtual bool setMaxBufferSize() = 0;
	
	/** Allocates memory for delay buffer and initialises all elements to 0
	 
	 @param bufferSizeSamples		the size of delay buffer to be used
	 
	 @returns   false and error description if there's a problem,
	 or sets the internal delay buffer and returns true
	 */
	bool setDelayBuffer(int bufferSizeSamples);
	
	/** store input sample into the delay buffer
	 
	 @param inputSample sample to be stored for delay (double)
	 */
	void storeSample(double inputSample);
	
	/** Increments the currentDelayWriteBufferIndex by 1
	 */
	void incDelayBuffWriteIndex();
	
	/** Increments the currentDelayBufferReadIndex by indexInc
		@param indexInc The amount to increment the delay buffer index
	 */
	void incDelayBuffReadIndex(double indexInc);
	
	/** sets the currentDelayBufferReadIndex by indexInc (Currently no wrapping)
		@param index the read index index required 
	 */
	void setDelayBuffReadIndex(double index);
	
	/** store input sample into the delay buffer and increment currentDelayWriteIndex
		for tracking when to loop back to start of buffer
	 
	 @param inputSample sample to be stored for delay (double)
	 */
	void delaySample(double inputSample);
	/** get the value of the requested buffer index by interpolating other points
		@param bufferIndex	The required buffer index
		
		@returns interpolated output
	 */
	double getInterpolatedOut(double bufferIndex);
	
protected:// member variables
	//==============================================================================
	/** Table of interpolation values as a 2D array indexed by
		interpolationTable[pointIndex][alphaIndex]
	 */
	static double** interpolationTable;

	/** buffer to stored audio buffer for delay effects*/
	double* delayBuffer = 0;
	
	/** Maximum number of samples that can be stored in delayBuffer*/
	int maxDelayBufferSize = 441000;
	
	/** the delay time of signal in samples*/
	int delayTimeSamples = 44100;
	
	/***/
	int currentDelayWriteIndex = 0;
	
	/***/
	double currentDelayReadIndex = 0;
	
	/***/
	static const int interpOrder = 4;
	
	/***/
	static const int interpResolution = 1000;
	
	/** internal class error boolean*/
	bool error;
};


//==============================================================================
/**
	A Base class for filter based effects including methods for simple
	high, low and band pass filtering
 
 @see FilterEffectBase
 */
class FilterEffectBase
{
public:		// methods
	//==========================================================================
	/** Constructor. */
    FilterEffectBase();
	/** Destructor. */
    ~FilterEffectBase();
	//==========================================================================
	/** with the current filter coefficients this method filters a
	 sample then stores it the sample Buffer and increments the index
	 
	 @param sampVal is the sample to be processed
	 
	 @returns filtered audio sample
	 */
	double applyFilter(double sampVal);
	/** detect the envelop of an incoming signal
	 
	 @param sample		the incoming signal sample value
	 
	 @returns returns envelope dection signal sample value
	 */
	 double envelope(double sample);
	//==========================================================================
	/***/
	void printBuffers();
	/***/
	void printCoefs();
	//==========================================================================
	/** changes the current Chebyshev type 1 coefficients without altering the filter
		order. This allows for use in an audio process thread as it avoids dynamic allocation
		of memory. Filter sample and coefficient buffers are unaltered
		
		for initialising a chebyshev type 1 filter @see setChebyICoefficients
		
		@params cutFreq		normalised cutoff frequency (0 < x < .5)
		@params shelfType	bool filter shelf type, false = low pass, true = high pass
		@params ripple		percentage ripple (<.2929)
		
		@returns boolean false on error and true on success
	 */
	bool setChebyICoefficients(double cutFreq, bool shelfType, double ripple);
	//==========================================================================
protected:	// methods
	//==========================================================================
	/** set firCoefficients and iirCoefficients
		for required chebyshev type I filter
		sampleBuffer memory is also set
		
		@params cutFreq		normalised cutoff frequency (0 < x < .5)
		@params shelfType	bool filter shelf type, false = low pass, true = high pass
		@params ripple		percentage ripple (<.2929)
		@params poles		number of poles
		
		@returns boolean false on error and true on success
	 */
	bool changeChebyICoefficients(double cutFreq, bool shelfType, double ripple,int poles);
	/** a simple normalised fir low pass filter
		@params order	number of delay coefficients
		*/
	//==========================================================================
	bool setSimpleLpf(int order);
private:	// methods
	//==========================================================================
	/** increment the buffer index and wrap it to the filter order*/
	void incBufferIndex();
	
	//==========================================================================
	/** checks internal memory storage of filter coeffcients and deletes if
		required
		*/
	void clearMemory();
	
	/** will allocate memory to a buffer given the current filter order and set
		all values == 0.00
		*/
	void allocateBufferMemory();
	/** root mean square of signal over a specific sample window */
	double rms(double sample);
	
public:		// variables
	//==========================================================================
protected:  // variables
	//==========================================================================
	/** Numerator coefficients in delay filter
		firCoefficients[0] z^0 coeffcieint
		firCoefficients[1] z^-1 coefficient
	 */
	double* firCoefficients = 0;
	/** Denomiator coefficients in delay filter
		@see firCoefficients
	 */
	double* iirCoefficients = 0;
	/** hold temporary values for fir coeffcient buffer*/
	double* firTemp = 0;
	/** hold temporary values for iir coeffcient buffer*/
	double* iirTemp = 0;
	/** buffer to hold forward delay sample data
	 */	
	double* firBuffer = 0;
	/** buffer to hold backward delay sample data
	 */
	double* iirBuffer = 0;
	/**current buffer index*/
	int bufferIndex = 0;
	/** order of delay filter including the zero delay coefficients*/
	int filterOrder = 0;
	/***/
	int samplingRate = 0;
	/** window size in samples of rms window*/
	const int rmsWindowSize = 128;
	/** current write index of rmsBuffer */
	int rmsBufferIndex = 0;
	/** RMS window buffer */
	double* rmsBuffer = new double[rmsWindowSize];
private:	// variables
	//==========================================================================
	/***/
	constexpr static const double pi = 3.141592653589793;
};

/**
 Class provides a wave table that can be populated with a number of preallocated
 waveforms. These can be used to generate audio in themselves or to modulate
 The parameters of another effect.
 
 Class initialised with sample rate.
 */
class ModulationBaseClass {
public:
    //==============================================================================
    /** Constructor */
    ModulationBaseClass();
    ModulationBaseClass(double extSampRate);
    /** Destructor */
    ~ModulationBaseClass();
    //==============================================================================
    /**
     setup the class with a given sample rate. Basically reperforming the
     constructor

     @param extSampRate External sample rate
     */
    void setupModulationBaseClass(double extSampRate);
    
    /**
     sets wavetable to one period of a triangle wave
     */
    void setTriangle();
    /**
     sets wavetable to one period of a square wave
     */
    void setSquare();
    /**
     sets wavetable to one period of a sawtooth wave
     */
    void setSawtooth();
    /**
     sets wavetable to one period of a sine wave oscillating between -1 and 1
     */
    void setSine();
    /**
     sets wavetable to one period of a sine wave oscillating between 0 and 1
     */
    void setOffSine();
    /**
     sets wavetable to DC one
     */
    void setDC();
    /** set wave table to be a ramp from 0 to 1 */
    void setRamp();
    /**
     reads out white noise
     
     @return random number between -1 and 1
     */
    double readNoise();
    /**
     clip wave table values with a tanh function. Effect change with a variable
     amp to control intensity.
     
     @param amp amount to multiply signal before being fed through a tanh function
     */
    void clipWave(double amp);
    //==============================================================================
    /**
     Read through values in waveTable as a given frequency
     
     @param freq read speed in Hz: essentially the number of samples
     jumped between reads
     @return value from table as double
     */
    double readTable(double freq);
    //==============================================================================
    void printInterpTable();
    //==============================================================================
    /**
     populates the internal interpolation table
     
     @return return tru on success, else false
     */
    bool setInterpTable();
private:
    //==============================================================================
    /**
     allocate memory to internal wave table based on sample rate
     
     @return returns true on success or false on failure
     */
    bool allocateMemory();
    //==============================================================================
    /**
     get the interpolated output of the waveTable from the
     given buffer index
     
     @param bufferIndex buffer index as double
     @return returns interpolated value from surrounding wavtable indices
     */
    double getInterpOut(double bufferIndex);
    
    /**
     get a cubic spline interpolated out from the wave table
     
     Derived from code by Alec Wright at repository:
     https://github.com/Alec-Wright/Chorus
     
     @authors Matthew Hamilton, Alec Wright
     @param bufferIndex the required buffer index
     @param freq (speed) that the table is being read through
     @return returns interpolated value as double
     */
    double getSplineOut(double bufferIndex, int freq);
public:
    //==============================================================================
    /** current table read index */
    double tableIndex = 0;
    /** Internal Sample Rate */
    int sampleRate;
    /** time between samples: 1/sampRate */
    double timeStep;
    /** store modulation signal as*/
    double* waveTable;
    
private:
    //==============================================================================
    static const int order = 4;
    static const int res = 100;
    double interpTable[order][res] = {1};
};



class SimpleDelay;
//==============================================================================
/**
	Simple Delay effect consiting of a single tap delay with Effect Gain and
	feed back controls
	
	Constructor requires internal delay in samples
	
	@see process
 */
class SimpleDelay : public DelayEffectBase
{
public: // Methods
	//==============================================================================
    /** Constructor: DigitalEffect Base Must Be initialised
     
     @param delayInSamples Set the amount of delay in samples
     
     @see DelayEffectBase constructor
     */
    SimpleDelay(int maxDelayInSamples, int samplingRate);
	
	/** Destructor. */
	~SimpleDelay();
	
	/** setDelayGain: sets the delay gain to a value between 1 and -1
	 
		@param gain		required delay gain. Values beyond 1 and -1
	 are capped to the maximum to avoid idiocy.
	 Negative velus invoke a phase inversion.
	 
		@see setFeedbackGain
	 */
	void setDelayGain(double gain);
	
	/** setDelayGain: sets the feedback gain to a value between 1 and -1
	 
		@param gain		required delay gain. Values beyond 1 and -1
	 are capped to the maximum to avoid idiocy.
	 Negative velus invoke a phase inversion.
	 
		@see setDelayGain
	 */
	void setFeedbackGain(double gain);
	
	//==============================================================================
	/**
	 Apply delay and return input sample along with delay buffer signal

	 @param inputSample <#inputSample description#>
	 @return <#return value description#>
	 */
	double process(double inputSample) override;
	
    /**
     <#Description#>

     @param delayInSamples <#delayInSamples description#>
     */
    void setupSimpleDelay(int delayInSamples);
    /**
     <#Description#>

     @param delayInSample <#delayInSample description#>
     */
    void setDelayTime(double delayInSamples);
    /**
     <#Description#>

     @param seconds <#seconds description#>
     */
    void setDelayTransitionTime(double seconds);
private: //Methods
	/** capGain: caps gain to a range of 1 and -1;
	*	@param gain address of gain value
	*/
	void capGain(double& gain);
    //==============================================================================
    /**
     get a cubic spline interpolated out from the wave table
     
     Derived from code by Alec Wright at repository:
     https://github.com/Alec-Wright/Chorus
     
     @authors Matthew Hamilton, Alec Wright
     @param bufferIndex the required buffer index
     @return returns interpolated value as double
     */
    double getSplineOut(double bufferIndex);
    //==============================================================================
private: //member vairables
	/***/
    double delayGain = .707;
    /***/
    double feedbackGain = 0.;
    /***/
    double readHeadIndex;
    /***/
    unsigned int writeHeadIndex;
    /***/
    double currentDelaySamples;
    /***/
    double targetDelaySamples;
    /** increment when transition from current to target delay per sample set by delayTransitionTime*/
    double delayIncrement;
    /** inverse of delay increment: for working out whole number cyles to reach target delay*/
    double invDelayIncrement;
    /** time in seconds to transition from one delay to another*/
    double delayTransitionTime;
    double delayTransitionTimeInSamples;
    /***/
    const int sampleRate;
    int count = 0;
    bool delayTimeChanged = false;
};


/** Delay effect that filters the repeat delay */
class FilteredDelay: public DelayEffectBase,
					 public FilterEffectBase
{
public:
	/** Constructor */
	FilteredDelay(int delayInSamples) : DelayEffectBase(44100)
	{
	delayTimeSamples = delayInSamples;
	changeChebyICoefficients(.05, true, .1, 4);
	};
	/** Destructor */
	~FilteredDelay(){};
public:
  	/** setDelayGain: sets the delay gain to a value between 1 and -1
	 
		@param gain		required delay gain. Values beyond 1 and -1
	 are capped to the maximum to avoid idiocy.
	 Negative velus invoke a phase inversion.
	 
		@see setFeedbackGain
	 */
	void setDelayGain(double gain);
	
	/** setDelayGain: sets the feedback gain to a value between 1 and -1
	 
		@param gain		required delay gain. Values beyond 1 and -1
	 are capped to the maximum to avoid idiocy.
	 Negative velus invoke a phase inversion.
	 
		@see setDelayGain
	 */
	void setFeedbackGain(double gain);
	
	//==============================================================================
	/**apply the DSP effect*/
	double process(double inputSample) override;
	
private: //Methods
	/** capGain: caps gain to a range of 1 and -1;
	*	@param gain address of gain value
	*/
	void capGain(double& gain);

private: //member vairables
	/***/
	double delayGain=.707, feedbackGain=0.0;

};


/** Class for applying a FIR lowpass filter
	
	Intialised with the order of filter and the normalised cutoff frequency
	
	object declared with order of low pass filter
	*/
class SimpleLPF : public FilterEffectBase
{
public:
/** Constructor: Intialised with the order of FIR filter
	
	minimum cutoff frequency is 1/sampleRate
	
	@param cutoff	normalised cutoff frequency [0, 1);
    @param order	filter order
	*/
	SimpleLPF(double cutoff, int order)
	{
		changeChebyICoefficients(cutoff, false, .1, order);
	};
public:
	/** destructor*/
	~SimpleLPF(){};
};


class EnvelopeFilter : public FilterEffectBase
{
public:
	/** Constructor */
	EnvelopeFilter(): envelopeFollower(.00006, 4)
	{
	// NOTE: Initialising chebyshev coeffcients allocates memory, perhaps alter so that memory is already pre allocated
	changeChebyICoefficients(.01, false, .1, 4);
	envelopeFollower.setChebyICoefficients(.00006, false, 0);
	};
	/** Desructor */
	~EnvelopeFilter(){};
	/** main process method: applies an envelope filter
		to the incoming signal sample
		
		@params sample		incoming signal sample value
		
		@returns processed sample value
		*/
		double process(double sample);
private:
	/** this follows the signal envelope and alters the internal
		low pass filter cutoff
	*/
	SimpleLPF envelopeFollower;
};


/** Simple Chorus effect with a single delay voice and mono output
 
	Chorus is effective between 15 and 20 miliseconds delay of original
	audio. Requires the sample rate when initialising.*/
class SimpleChorus : public DelayEffectBase,
                     public ModulationBaseClass,
                     public SimpleLPF
{
public:
    /** Constructor: initialises the effect parameters
     Since the class inherits DelayEffectBase it must set a
     delay buffer size when initialising.
     
     */
    SimpleChorus() : SimpleLPF(0.0001,4)
    {
    }
    SimpleChorus(int extSampleRate) : ModulationBaseClass(extSampleRate),
                                      DelayEffectBase(double(extSampleRate)*.031),
                                      SimpleLPF(0.0001,4)
    {
        swing = 0.005*sampleRate;
        base = 0.015*sampleRate;
        if (sampleRate != 0)
            setRandLfo();
    }
	//==============================================================================
    /**
     apply chorus effect to audio sample

     @param inputSample input audio sample
     @return processed audio sample
     */
    double process(double inputSample);
    
    /**
     set parameters and internal sample rate

     @param extSampleRate external sample rate
     */
    void setupChorus(double extSampleRate);
    
    /**
     set the 'swing' of the chorus: The amount intensity of vibrato in the delay
     signal

     @param swingAmount <#swingAmount description#>
     */
    void setSwing(double swingAmount);
    
    
    /**
     sets the 'base' of the chorus: the minimum delay in the signal.

     @param baseAmount <#baseAmount description#>
     */
    void setBase(double baseAmount);
private:
    //==============================================================================
	/** swing range of delayed audio index in samples
	 typical maximum swing would be 15 milliseconds
		(i.e. 0.015*sampleRate)*/
	double swing;
	/** minimum delay in samples. Typically 10 milliseconds */
	double base;
	/** the minimum value of the lowpassed random modulation signal*/
	double modMin = .5;
	/** the maximum value of the lowpassed random modulation signal*/
	double modMax = .5;
	/** the normalising factor of the random delay signal*/
	double modNorm = 1/(modMax - modMin);
	
	/** modulation signal scaling equation: ((n - modMin)*modNorm*swing) + base
	 
		modulates a random white noise signal by lowpass filtering then scaling
		to a range between 15 to 25 miliseconds of delay.*/
	double getModSignal();
    
    const double readSpeed = ((readNoise()+1)*.5)*.0005;
    
    void setRandLfo();
};


class SimpleFlanger;
//==============================================================================
/**
	Simple Flanger Effect Consistig of a single voice flanger
	The flanger has an effective range between 0 and 15 miliseconds
	in this case dleay buffer should be set to sampleRate*3/200
	
	Constructor requires internal delay in samples
	
	@see process
 */
class SimpleFlanger : public DelayEffectBase
{
public: // Methods
	//==============================================================================
	/** Constructor: DigitalEffect Base Must Be initialised
		
		@see DelayEffectBase constructor
	 */
    SimpleFlanger();
    SimpleFlanger(double extSampleRate);
	
	/** Destructor. */
    ~SimpleFlanger();
	
	/** setEffectGain: sets the effect gain to a value between 1 and -1
	 
		@param gain		required delay gain. Values beyond 1 and -1
	 are capped to the maximum to avoid idiocy.
	 Negative velus invoke a phase inversion.
	 */
	void setEffectGain(double gain);
    
    /**
     <#Description#>

     @param depth <#depth description#>
     */
    void setDepth(const double depth);
    
    /**
     <#Description#>

     @param rate <#rate description#>
     */
    void setRate(const double rate);
	
	/** setEffectGain: sets the parameters for effect
		@param	gain	effect gain
		@param	depth	depth of modulation in samples
		@param	rate	rate of modulation in Hz
	 */
	void setEffectParams(double gain, double depth, double rate);
	//==============================================================================
	/**Apply the DSP effect*/
	double process(double inputSample) override;
	
    void setupSimpleFlanger(double extSampleRate);
private: //Methods
	/** capGain: caps gain to a range of 1 and -1;
	 *	@param gain address of gain value
	 */
	double capGain(double gain);
	
	/** setAngleDelta: sets the angleDelta for delay modulation*/
    void setAngleDelta();
	
	/** updateModulation: updates the modulationIndex by the correct
		increment*/
	void updateModulation();
	
private: //member vairables
	/** internal class declaration of pi
		it would likely make sense to have this moved
		to a higher class
		*/
	constexpr static const double internal_Pi = 3.141592653589793;
	
	/***/
	double modulationDepth = 1000, modulationRate = 0, effectGain = .01;
	
	//	M1 = M01*(1+sin(2*pi*f1.*n/FS))
	/***/
	double modulationIndex = 0;
	
	/** 1/sampleRate: The time in seconds between samples*/
	double timeStep = 1./44100.;
	
	/** 2 * pi * modulationRate / sampleRate */
	double modulationConstant, modulationAngle = 0;
	
	//const double cyclesPerSample = modulationRate * timeStep;
	/**increment value for modulation signal*/
	double angleDelta = 2*internal_Pi*timeStep;
};


DelayEffectBase::DelayEffectBase()
{
}

DelayEffectBase::DelayEffectBase(int bufferSizeSamples)
{
    error = setDelayBuffer(bufferSizeSamples);
    delayTimeSamples = bufferSizeSamples;
}

DelayEffectBase::~DelayEffectBase()
{
    delete[] delayBuffer;
}

void DelayEffectBase::setupDelayEffectBase(const int bufferSizeSamples)
{
//    delete [] delayBuffer;
    error = setDelayBuffer(bufferSizeSamples);
    delayTimeSamples = bufferSizeSamples;
}
//==============================================================================
double** DelayEffectBase::setInterpolationTable()
{
	const int order = interpOrder;
	const int res = interpResolution;
	double** interpTable = new double*[order];
	if(!interpTable){return NULL;}
	
	for(int i=0;i<order;i++)
	{
		interpTable[i] = new double[res+1];
		if(!interpTable[i]){return NULL;}
		std::fill(interpTable[i], interpTable[i]+res, 1);
	}
	
	double *polynomial_normaliser = new double [order];
	if(!polynomial_normaliser){return NULL;}
	std::fill(polynomial_normaliser, polynomial_normaliser+order, 1);
	double *alphas = new double [res];
	if(!alphas){return NULL;}
	
	for(int i = 0; i<res;i++)
	{
		alphas[i] = (i/float(res)) - 0.5;
	}
	
	double *anchors = new double [order];
	
	if ((order % 2)==0)
	{
		for(int i = 0; i<order;i++)
		{
			anchors[i] = -(double(order) - 1)*0.5 + double(i);
			
		}
	}
	else
	{
		for(int i = 0; i<order;i++)
		{
			anchors[i] = (-(double(order))*0.5) + double(i);
		}
	}
	
	for (int q = 0; q<res;q++) // loop for every value of alpha
	{
		for (int j = 0; j<order;j++) // loop for sub polynomial
		{
			for (int m = 0; m < order; m++) //loop for each point in subpoly
			{
				if (m != j)
				{
					if (q == 0)
					{
						polynomial_normaliser[j] = polynomial_normaliser[j]*(anchors[j]-anchors[m]);
					}
					interpTable[j][q] *= (alphas[q]-anchors[m]);
					
				}
			}
			interpTable[j][q] /= polynomial_normaliser[j];
		}
	}
	delete[] polynomial_normaliser;
	delete[] alphas;
	delete[] anchors;
	return interpTable;
}
//==============================================================================

double DelayEffectBase::getInterpolatedOut(double bufferIndex)
{
	const int order = interpOrder;
	const int orderHalf = order*.5;
	const int res = interpResolution;
	double interpOut = 0;
	int intBufferIndex = floor(bufferIndex);
	int alphaIndex = int(floor((bufferIndex-intBufferIndex)*res));
	
	for(int i = 0; i<order;i++)
	{
		
		int interpIndex = (i+1 - orderHalf) + intBufferIndex;
		
		if(interpIndex < 0 || (interpIndex >= maxDelayBufferSize))
		{
			if(interpIndex < 0){interpIndex = maxDelayBufferSize+interpIndex;}
			else{interpIndex = interpIndex-maxDelayBufferSize;}
		}
		
		interpOut += (interpolationTable[i][alphaIndex]) * (delayBuffer[interpIndex]);
	}
	return interpOut;
}

//==============================================================================
bool DelayEffectBase::setDelayBuffer(int bufferSizeSamples)
{
	maxDelayBufferSize = bufferSizeSamples;
	delayBuffer = new double [maxDelayBufferSize];
	if(!delayBuffer){return false;}
	std::fill(delayBuffer, delayBuffer+maxDelayBufferSize, 0);
	return true;
}
//==============================================================================
void DelayEffectBase::storeSample(double inputSample)
{
	delayBuffer[currentDelayWriteIndex] = inputSample;
}
void DelayEffectBase::incDelayBuffWriteIndex()
{
	currentDelayWriteIndex++;
    currentDelayWriteIndex %= delayTimeSamples;
//	if(currentDelayWriteIndex>=delayTimeSamples){currentDelayWriteIndex=0;}
}

void DelayEffectBase::incDelayBuffReadIndex(double indexInc)
{
	currentDelayReadIndex += indexInc;
	if(currentDelayReadIndex>=double(delayTimeSamples)){currentDelayReadIndex=0;}
	if(currentDelayReadIndex<0){currentDelayReadIndex=0;}
}

void DelayEffectBase::setDelayBuffReadIndex(double index)
{
	currentDelayReadIndex = index;
	if(currentDelayReadIndex>=double(delayTimeSamples)){currentDelayReadIndex=0;}
	if(currentDelayReadIndex<0){currentDelayReadIndex=0;}
}

void DelayEffectBase::delaySample(double inputSample)
{
	storeSample(inputSample);
	incDelayBuffWriteIndex();
}

//==============================================================================
void DelayEffectBase::printInterpTable()
{
	for (int j = 0; j<interpResolution; j++)
	{
		for (int i = 0; i<interpOrder;i++)
		{
			printf("index %d: %.2f \t",i,interpolationTable[i][j]);
		}
		printf("\n");
	}
}

//==============================================================================
// set static variables
//==============================================================================
double** DelayEffectBase::interpolationTable = DelayEffectBase::setInterpolationTable();


FilterEffectBase::FilterEffectBase()
{
    std::fill(rmsBuffer, rmsBuffer+rmsWindowSize, 0);
}

FilterEffectBase::~FilterEffectBase()
{
    
}

//==============================================================================
double FilterEffectBase::applyFilter(double sampVal)
{
	double outSample = 0;
	firBuffer[bufferIndex] = sampVal;
	
	for(int j = 0; j<filterOrder; j++)
	{
		int i = ((bufferIndex-j)+filterOrder)%filterOrder;
		outSample += firCoefficients[j]*firBuffer[i] + iirCoefficients[j]*iirBuffer[i];
	}
	
	iirBuffer[bufferIndex] = outSample;
	incBufferIndex();
	
	return outSample;
}

//==============================================================================
double FilterEffectBase::rms(double sample)
{
	rmsBuffer[rmsBufferIndex] = sample;
	double rmsValue = 0;
	for(int j = 0; j < rmsBufferIndex; j++)
	{
	int i = ((rmsBufferIndex-j)+rmsWindowSize)%rmsWindowSize;
	rmsValue += rmsBuffer[i]*rmsBuffer[i];
	}
	
//	printf("samp: %e\tsum: %e\n", sample, rmsValue);
	
	rmsValue /= rmsWindowSize;
	rmsValue = sqrt(rmsValue);
	
	rmsBufferIndex++;
	rmsBufferIndex%= rmsWindowSize;
	
	return rmsValue;
}
//==============================================================================
double FilterEffectBase::envelope(double sample)
{
	return applyFilter(rms(sample));
}

//==============================================================================

void FilterEffectBase::incBufferIndex()
{
	bufferIndex++;
	bufferIndex%=filterOrder;
}

//==============================================================================

bool FilterEffectBase::setChebyICoefficients(double cutFreq, bool shelfType, double ripple)
{
	//NOTE: coefficient buffers must be cleared as are additive in the following
	//		code
	std::fill(firCoefficients, firCoefficients+22, 0);
	std::fill(iirCoefficients, iirCoefficients+22, 0);
	
	double poles = (double)filterOrder-1;
	int order = (int)poles;

	firCoefficients[2] = 1;
	iirCoefficients[2] = 1;
	
	double Es,Vx,Kx;
	if(ripple!=0)
	{
		Es = sqrt( pow(1/(1-ripple),2) - 1);
		Vx = (1/poles)*log(1/Es + sqrt(1/(pow(Es,2))+1));
		Kx = (1/poles)*log(1/Es + sqrt(1/(pow(Es,2))-1));
		Kx = cosh(Kx);
	}
	else
	{
		Vx = 1;
		Kx = 1;
	}
	
	const double T = 2*tan(.5);
	const double W = 2*pi*cutFreq;
	
	double K;
	
	if(shelfType==0) //// if low pass
	{
		K = sin(.5 - W/2)/sin(.5 + W/2);
	}
	else //// if high pass
	{
		
		K = -cos(.5 + W/2)/cos(W/2 - .5);
	}
	
	////// main algorithm
	for (int i = 0; i<(order/2); i++)
	{
		////// Sub routine
		const double alpha = pi/(2*poles) + (i-1)*(pi/poles);
		
		double Rp, Ip;
		if(ripple!=0)
		{
			Rp = -cos(alpha)*sinh(Vx)/Kx;
			Ip = sin(alpha)*cosh(Vx)/Kx;
		}
		else
		{
			Rp = -cos(alpha);
			Ip = sin(alpha);
		}
		
		const double M = pow(Rp,2) + pow(Ip,2);
		const double D = 4 - 4*Rp*T + M*T;
		
		const double X0 = (pow(T,2))/D;
		const double X1 = (2*pow(T,2))/D;
		const double X2 = X0;
		
		const double Y1 = (8-(2*M*pow(T,2)))/D;
		const double Y2 = (-4 - 4*Rp*T - M*T)/D;
		
		const double D1 = 1/(1 + Y1*K - Y2*pow(K,2)); // renamed and inverted from original algorithm
		
		const double A0 =  (X0 - X1*K + X2*pow(K,2))*D1;
		double A1 =  (-2*X0*K + X1 + X1*pow(K,2) - 2*X2*K)*D1;
		const double A2 =  (X0*pow(K,2) - X1*K + X2)*D1;
		
		double B1 = (2*K + Y1 +Y1*pow(K,2) - 2*Y2*K)*D1;
		const double B2 = (-(pow(K,2)) - Y1*K + Y2)*D1;
		
		if(shelfType==1){A1 = -A1; B1 = -B1;}
		
		for (int j = 0; j<22; j++)
		{
			firTemp[j] = firCoefficients[j];
			iirTemp[j] = iirCoefficients[j];
		}
		for (int j = 2; j<22; j++)
		{
			firCoefficients[j] = A0*firTemp[j] + A1*firTemp[j-1] + A2*firTemp[j-2];
			iirCoefficients[j] =    iirTemp[j] - B1*iirTemp[j-1] - B2*iirTemp[j-2];
		}
	} //// end for i
	
	iirCoefficients[2] = 0;
	for (int j = 0; j<filterOrder; j++)
	{
		firCoefficients[j] = firCoefficients[j+2];
		iirCoefficients[j] = -iirCoefficients[j+2];
	}
	//////// Normalising
	double SA = 0; double SB = 0;
	if (shelfType==0)
	{
		for (int j = 0; j<filterOrder; j++)
		{
			SA += firCoefficients[j];
			SB += iirCoefficients[j];
		}
	}
	else
	{
		for (int j = 0; j<order; j++)
		{
			SA += firCoefficients[j]*pow(-1,j);
			SB += iirCoefficients[j]*pow(-1,j);
		}
	}
	
	const double gain = SA/(1-SB);
	for (int j = 0; j<filterOrder; j++)
	{
		firCoefficients[j] /= gain;
	}
	
	return true;
}

bool FilterEffectBase::changeChebyICoefficients(double cutFreq, bool shelfType, double ripple,int poles)
{
	filterOrder = poles+1;
	clearMemory();
	allocateBufferMemory();
	setChebyICoefficients(cutFreq, shelfType, ripple);

	return true;
}

//==============================================================================
bool FilterEffectBase::setSimpleLpf(int order)
{
	filterOrder = order;
	clearMemory();
	allocateBufferMemory();
	firCoefficients = new double[filterOrder];
	iirCoefficients = new double[filterOrder];
	std::fill(iirCoefficients, iirCoefficients+filterOrder, 0);
	int coef = 1;
	double gain = 0;
	for(int j = 0; j<filterOrder; j++)
	{
		if(j==0)
		{coef = 1;}
		else
		{coef = coef*(filterOrder-j)/j;}
		
		firCoefficients[j] = (double)coef;
		gain += firCoefficients[j];
	}
	
	for(int j = 0; j<=filterOrder; j++)
	{
		firCoefficients[j] /= gain;
	}
	
	return true;
}

//==============================================================================

void FilterEffectBase::clearMemory()
{
	
	if(firCoefficients)
	{
		delete[] firCoefficients;
	}
	
	if(iirCoefficients)
	{
		delete[] iirCoefficients;
	}
}

//==============================================================================

void FilterEffectBase::allocateBufferMemory()
{
	if(firBuffer)
	{
		delete[] firBuffer;
	}
	
	if(iirBuffer)
	{
		delete[] iirBuffer;
	}
	firBuffer = new double[filterOrder];
	iirBuffer = new double[filterOrder];
	std::fill(firBuffer, firBuffer+filterOrder, 0);
	std::fill(iirBuffer, iirBuffer+filterOrder, 0);
	
		if(firCoefficients)
	{
		delete[] firCoefficients;
	}
	
	if(iirCoefficients)
	{
		delete[] iirCoefficients;
	}
	
		if(firTemp)
	{
		delete[] firTemp;
	}
	
	if(iirTemp)
	{
		delete[] iirTemp;
	}
	
	firCoefficients = new double[22];
	iirCoefficients = new double[22];
	firTemp = new double[22];
	iirTemp = new double[22];
	std::fill(firCoefficients, firCoefficients+22, 0);
	std::fill(iirCoefficients, iirCoefficients+22, 0);
	std::fill(firTemp, firTemp+22, 0);
	std::fill(iirTemp, iirTemp+22, 0);

}
//==============================================================================

void FilterEffectBase::printBuffers()
{
	printf("FIRb\t\tIIRb\n");
 for (int i = 0; i<filterOrder;i++)
 {
	 printf("%.4e\t%.4e\n",firBuffer[i],iirBuffer[i]);
 }
 printf("\n");
}

void FilterEffectBase::printCoefs()
{
	printf("FIR\t\tIIR\n");
 for (int i = 0; i<filterOrder;i++)
 {
	 printf("%.4e\t%.4e\n",firCoefficients[i],iirCoefficients[i]);
 }
 printf("\n");
}


ModulationBaseClass::ModulationBaseClass()
{
    srand (static_cast <unsigned> (time(0)));
}

ModulationBaseClass::ModulationBaseClass(double extSampRate)
{
    sampleRate = extSampRate;
    timeStep = 1./extSampRate;
    allocateMemory();
    //                setInterpTable();
    srand (static_cast <unsigned> (time(0)));
}

ModulationBaseClass::~ModulationBaseClass()
{
}

void ModulationBaseClass::setupModulationBaseClass(double extSampRate)
{
    if (waveTable != nullptr)
    {
//        delete [] waveTable;
    }
    sampleRate = extSampRate;
    timeStep = 1./extSampRate;
    allocateMemory();
}

//==============================================================================
double ModulationBaseClass::readTable(double freq)
{
    if (freq > 0)
    {
        //    const double out = getInterpOut(tableIndex);
        const double out = getSplineOut(tableIndex, int(freq));
        tableIndex += freq;
        if (tableIndex-sampleRate > 0)
            tableIndex -= sampleRate;
        
        return out;
    }
    else
    {
        return 0.;
    }
}
//==============================================================================
bool ModulationBaseClass::allocateMemory()
{
    waveTable = new double[sampleRate];
    if(!waveTable){
        return false;
    }
    std::fill(waveTable, waveTable+sampleRate, 0);
    return true;
}

//==============================================================================

bool ModulationBaseClass::setInterpTable()
{
    double *polynomial_normaliser = new double [order];
    if(!polynomial_normaliser){return false;}
    std::fill(polynomial_normaliser, polynomial_normaliser+order, 1);
    double *alphas = new double [res];
    if(!alphas){return false;}
    
    for(int i = 0; i<res;i++)
    {
        alphas[i] = (i/float(res)) - 0.5;
    }
    
    double *anchors = new double [order];
    
    if ((order % 2)==0)
    {
        for(int i = 0; i<order;i++)
        {
            anchors[i] = -(double(order) - 1)*0.5 + double(i);
            std::fill(interpTable[i], interpTable[i]+res, 1);
        }
    }
    
    else
    {
        for(int i = 0; i<order;i++)
        {
            anchors[i] = (-(double(order))*0.5) + double(i);
        }
    }
    
    for (int q = 0; q<res;q++) // loop for every value of alpha
    {
        for (int j = 0; j<order;j++) // loop for sub polynomial
        {
            for (int m = 0; m < order; m++) //loop for each point in subpoly
            {
                if (m != j)
                {
                    if (q == 0)
                    {
                        polynomial_normaliser[j] = polynomial_normaliser[j]*(anchors[j]-anchors[m]);
                    }
                    interpTable[j][q] *= (alphas[q]-anchors[m]);
                }
            }
            interpTable[j][q] /= polynomial_normaliser[j];
        }
    }
    
    
    delete[] polynomial_normaliser;
    delete[] alphas;
    delete[] anchors;
    return true;
}

double ModulationBaseClass::getInterpOut(double bufferIndex)
{
    const int order = 4;
    const int orderHalf = order*.5;
    const int res = 100;
    double interpOut = 0;
    int intBufferIndex = floor(bufferIndex);
    int alphaIndex = int(floor((bufferIndex-intBufferIndex)*res));
    
    for(int i = 0; i<order;i++)
    {
        int interpIndex = (((i+1 - orderHalf) + intBufferIndex)+sampleRate)%sampleRate;
        interpOut += (interpTable[i][alphaIndex]) * (waveTable[interpIndex]);
    }
    return interpOut;
}

void ModulationBaseClass::printInterpTable()
{
    for (int j = 0; j < res; j++){
        for (int i = 0; i < order; i++) {
            std::cout << interpTable[i][j] << '\t';
        }
        std::cout << '\n';
    }
}
//==============================================================================
double ModulationBaseClass::getSplineOut(double bufferIndex, int freq)
{
    
    if (freq<1) {
        freq = 1;
    }
    const int n0 = floor(bufferIndex);
    const int n1 = (n0+freq)%sampleRate;
    const int n2 = (n0+(2*freq))%sampleRate;
    const double alpha = bufferIndex - n0;
    const double a = waveTable[n1];
    const double c = ((3*(waveTable[n2] - waveTable[n1])) -  (3*(waveTable[n1] - waveTable[n0]))) *.25;
    const double b = (waveTable[n2] - waveTable[n1]) - ((2*c))/3;
    const double d = (-c)/3;
    return a + (b*alpha) + (c*alpha*alpha) + (d*alpha*alpha*alpha);
}

//==============================================================================

void ModulationBaseClass::setSine()
{
    const double radPerSec = 2*3.1415926536*timeStep;
    for(int i = 0; i < sampleRate; i++)
        waveTable[i] = sin(i*radPerSec);
}
//==============================================================================
void ModulationBaseClass::setOffSine()
{
    const double radPerSec = 2*3.1415926536*timeStep;
    for(int i = 0; i < sampleRate; i++)
        waveTable[i] = (sin(i*radPerSec) + 1)*.5;
}

//==============================================================================

void ModulationBaseClass::setSawtooth()
{
    std::fill(waveTable, waveTable+sampleRate, 0);
    const double radPerSec = 2*3.1415926536*timeStep;
    for(int i = 0; i < sampleRate; i++)
    {
        for (int j = 1; j <11; j+=1)
            waveTable[i] += pow(-1,j)*sin(j*radPerSec*i)/double(j);
    }
}

//==============================================================================
void ModulationBaseClass::setSquare()
{
    std::fill(waveTable, waveTable+sampleRate, 0);
    const double radPerSec = 2*3.1415926536*timeStep;
    for(int i = 0; i < sampleRate; i++)
    {
        for (int j = 0; j <35; j+=1)
            waveTable[i] += (sin((2*j + 1)*i*radPerSec))/(2*j +1);
    }
}

//==============================================================================
void ModulationBaseClass::setTriangle()
{
    std::fill(waveTable, waveTable+sampleRate, 0);
    const double radPerSec = 2*3.1415926536*timeStep;
    for(int i = 0; i < sampleRate; i++)
    {
        for (int j = 0; j <35; j+=1)
            waveTable[i] += pow(-1.,j)*(sin((2.*double(j) + 1)*i*radPerSec))/(2.*double(j) +1);
    }
}

//==============================================================================

void ModulationBaseClass::setDC()
{
    for (int i = 0; i < sampleRate; i++)
        waveTable[i] = 1.0;
}

//==============================================================================

void ModulationBaseClass::setRamp()
{
    for (int i = 0; i < sampleRate; i++)
        waveTable[i] = i / double(sampleRate);
}

//==============================================================================

void ModulationBaseClass::clipWave(double amp)
{
    if (amp < .01)
    {
        amp = .01;
    }
    
    for(int i = 0; i < sampleRate; i++)
        waveTable[i] = tanh(amp*waveTable[i])/tanh(amp);
}

double ModulationBaseClass::readNoise()
{
    const double lo = -1.;
    const double hi =  1.;
    return lo + static_cast <double> (rand()) /( static_cast <double> (RAND_MAX/(hi-lo)));
}



//==============================================================================

double FilteredDelay::process(double inputSample)
{
	delaySample(applyFilter((inputSample*delayGain)+ feedbackGain*getInterpolatedOut(currentDelayWriteIndex)));
	const double out = getInterpolatedOut(currentDelayWriteIndex) + inputSample;
	return out;
}

//==============================================================================

void FilteredDelay::capGain(double& gain)
{
	if (gain > 1.)
	{
		gain = 1.;
	}
	else if (gain < -1.)
	{
		gain = -1.;
	}
	return;
}

void FilteredDelay::setDelayGain(double gain)
{
	capGain(gain);
	delayGain = gain;
}

void FilteredDelay::setFeedbackGain(double gain)
{
	capGain(gain);
	feedbackGain = gain;
}
//==============================================================================


void SimpleChorus::setSwing(double swingAmount)
{
    swing = swingAmount*sampleRate;
}

void SimpleChorus::setBase(double baseAmount)
{
    base = baseAmount*sampleRate;
}

void SimpleChorus::setupChorus(double extSampleRate)
{
    setupModulationBaseClass(extSampleRate);
    setupDelayEffectBase(double(extSampleRate)*.1);
    //    SimpleLPF(0.0004,4)
    setChebyICoefficients(0.00005, false, 0);
    
    swing = readSpeed*extSampleRate*5;
    base = readSpeed*extSampleRate*20;
//    delete [] waveTable;
    setRandLfo();
}

//==============================================================================
double SimpleChorus::process(double inputSample)
{
    delaySample(inputSample);
    const double waveDelay = getModSignal();
    const double delayAmount = ((int(currentDelayWriteIndex - waveDelay) + delayTimeSamples) % delayTimeSamples)
                                + ((currentDelayWriteIndex - waveDelay) - trunc(currentDelayWriteIndex - waveDelay));
    
    const double out = .0*inputSample + 1.*getInterpolatedOut(delayAmount);
    
    return out;
}
//==============================================================================
void SimpleChorus::setRandLfo()
{
    std::fill(iirBuffer, iirBuffer+filterOrder, .5);
    for (int i = 0; i < sampleRate; i++)
    {
        waveTable[i] = (readNoise()+1)*.5;
//        waveTable[i] = applyFilter((readNoise()+1)*.5);
        if (waveTable[i] < modMin)
            modMin = waveTable[i];
        if (waveTable[i] > modMax)
        {
            modMax = waveTable[i];
        }
    }
    
    modNorm = 1/(modMax-modMin);
    
    // normalises the delay signal
    for (int i = 0; i < sampleRate; i++)
    {
        waveTable[i] -= modMin;
        waveTable[i] *= modNorm;
    }
    
//    setOffSine();
    
    // this code fades out at the end and fades in at the start
    // to avoid any discontinuities int the signal.
//    const int fadeSize = 10000;
//    const double fadeSpeed = 2*M_PI/fadeSize;
//    for (int i = 0; i < fadeSize; i++)
//    {
//        const int fadeIndex = ((sampleRate-fadeSize/2)+i)%sampleRate;
//        waveTable[fadeIndex] *= (1+cos(fadeSpeed*i))*.5;
//    }
}
//==============================================================================
double SimpleChorus::getModSignal()
{
    return (readTable(readSpeed)*swing) + base;
}

//==============================================================================


//==============================================================================
SimpleDelay::SimpleDelay(int maxDelayInSamples, int samplingRate) : DelayEffectBase(maxDelayInSamples), sampleRate(samplingRate)
{
    writeHeadIndex = 0;
    readHeadIndex = 1;
    currentDelaySamples = maxDelayInSamples;
    targetDelaySamples = maxDelayInSamples;
    setDelayTransitionTime(0.5);
}

SimpleDelay::~SimpleDelay()
{
}
//==============================================================================

void SimpleDelay::setupSimpleDelay(int delayInSamples)
{
    setupDelayEffectBase(delayInSamples);
}

//==============================================================================

double SimpleDelay::process(double inputSample)
{
    // write sample
    delayBuffer[writeHeadIndex] = inputSample;
    writeHeadIndex++;
    writeHeadIndex %= maxDelayBufferSize;
    
    // read sample
    double outSample = getSplineOut(readHeadIndex) + (inputSample * 1);
    if (delayTimeChanged)
    {
        count++;
        const double difference = (currentDelaySamples - targetDelaySamples);
        const double increment = delayIncrement * (difference / fabs(difference));
        currentDelaySamples -= increment;
        readHeadIndex += 1 + increment;
        readHeadIndex = std::fmod(readHeadIndex, maxDelayBufferSize);
        if (count > floor(delayTransitionTimeInSamples))
        {
            currentDelaySamples = targetDelaySamples;
            readHeadIndex = floor(readHeadIndex);
            delayTimeChanged = false;
        }
    }
    else
    {
        readHeadIndex++;
        readHeadIndex = std::fmod(readHeadIndex, maxDelayBufferSize);
    }
    return outSample;
}

//==============================================================================

void SimpleDelay::capGain(double& gain)
{
    if (gain > 1.)
    {
        gain = 1.;
    }
    else if (gain < -1.)
    {
        gain = -1.;
    }
    return;
}

void SimpleDelay::setDelayGain(double gain)
{
    capGain(gain);
    delayGain = gain;
    
}

void SimpleDelay::setFeedbackGain(double gain)
{
    capGain(gain);
    feedbackGain = gain;
}

//==============================================================================
double SimpleDelay::getSplineOut(double bufferIndex)
{
    const int n0 = floor(bufferIndex);
    const int n1 = (n0 + 1) % maxDelayBufferSize;
    const int n2 = (n0 + 2) % maxDelayBufferSize;
    const double alpha = bufferIndex - n0;
    
    const double a = delayBuffer[n1];
    const double c = ((3 * (delayBuffer[n2] - delayBuffer[n1])) -  (3 * (delayBuffer[n1] - delayBuffer[n0]))) * 0.25;
    const double b = (delayBuffer[n2] - delayBuffer[n1]) - (2 * c * 0.33333);
    const double d = (-c) * 0.33333;
    return a + (b * alpha) + (c * alpha * alpha) + (d * alpha * alpha * alpha);
}
//==============================================================================
void SimpleDelay::setDelayTransitionTime(double seconds)
{
    delayTransitionTime = seconds;
    delayTransitionTimeInSamples = seconds * sampleRate;
//    std::cout << delayTransitionTime << '\n' << delayTransitionTimeInSamples << '\n' << '\n';
}

void SimpleDelay::setDelayTime(double delayInSamples)
{
    delayTimeChanged = true;
    targetDelaySamples = delayInSamples;
    const double delayTimeDifference = currentDelaySamples - targetDelaySamples;
    delayIncrement = delayTimeDifference / delayTransitionTimeInSamples;
    count = 0;
}

SimpleFlanger::SimpleFlanger()
{
}

SimpleFlanger::SimpleFlanger(double extSampleRate) : DelayEffectBase(static_cast<int>(extSampleRate*0.02))
{
}

/** Destructor. */
SimpleFlanger::~SimpleFlanger()
{
}

void SimpleFlanger::setupSimpleFlanger(double extSampleRate)
{
    setupDelayEffectBase(extSampleRate*.02);
    timeStep = 1./extSampleRate;
    setEffectParams(.707, extSampleRate*.02, .1);
}
//==============================================================================

double SimpleFlanger::process(double inputSample)
{
	delaySample(inputSample);
	const double out = ((1-fabs(effectGain*.2))*(inputSample) + (effectGain * getInterpolatedOut(modulationIndex)));
	updateModulation();
	return out;
}
//==============================================================================
void SimpleFlanger::updateModulation() //TODO: swap for usage of ModulationBaseClassInheritance
{
	modulationAngle += angleDelta;
	modulationIndex = (currentDelayWriteIndex-(modulationDepth*(1+(sin(modulationAngle))))) - 12;
    modulationIndex = ( (int(modulationIndex) + delayTimeSamples) % delayTimeSamples)
                      + (modulationIndex - floor(modulationIndex) );
}
//==============================================================================

double SimpleFlanger::capGain(double gain)
{
	if (gain > 1.)
	{
		gain = 1.;
	}
	else if (gain < -1.)
	{
		gain = -1.;
	}
	return gain;
}

//==============================================================================
void SimpleFlanger::setEffectGain(const double gain)
{
    effectGain = capGain(gain);
}


void SimpleFlanger::setDepth(const double depth)
{
    if (depth > double(delayTimeSamples))
        modulationDepth = double(delayTimeSamples)-1;
    else
        modulationDepth = depth;
}

void SimpleFlanger::setRate(const double rate)
{
    modulationRate = rate;
    setAngleDelta();
}
//==============================================================================
void SimpleFlanger::setAngleDelta()
{
	const double cyclesPerSample = modulationRate * timeStep;
	angleDelta = cyclesPerSample * 2.0 * internal_Pi;
}

void SimpleFlanger::setEffectParams(const double gain, const double depth, const double rate)
{
    setEffectGain(gain);
    setDepth(depth);
    setRate(rate);
}


double EnvelopeFilter::process(double sample)
{
	setChebyICoefficients(0.001+envelopeFollower.envelope(2*sample),false,.1); //Offset avoids zero cutoff value
	return applyFilter(sample);
}
