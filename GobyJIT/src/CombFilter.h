/**
\struct CombFilterParameters
\ingroup FX-Objects
\brief
Custom parameter structure for the CombFilter object. Used for reverb algorithms in book.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
struct CombFilterParameters
{
	CombFilterParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	CombFilterParameters& operator=(const CombFilterParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		delayTime_mSec = params.delayTime_mSec;
		RT60Time_mSec = params.RT60Time_mSec;
		enableLPF = params.enableLPF;
		lpf_g = params.lpf_g;
		interpolate = params.interpolate;
		return *this;
	}

	// --- individual parameters
	double delayTime_mSec = 0.0;	///< delay time in mSec
	double RT60Time_mSec = 0.0;		///< RT 60 time ini mSec
	bool enableLPF = false;			///< enable LPF flag
	double lpf_g = 0.0;				///< gain value for LPF (if enabled)
	bool interpolate = false;		///< interpolation flag (diagnostics)
};


/**
\class CombFilter
\ingroup FX-Objects
\brief
The CombFilter object implements a comb filter with optional LPF in feedback loop. Used for reverb algorithms in book.
Audio I/O:
- Processes mono input to mono output.
Control I/F:
- Use CombFilterParameters structure to get/set object params.
\author Will Pirkle http://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed. by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class CombFilter : public IAudioSignalProcessor
{
public:
	CombFilter(void) {}		/* C-TOR */
	~CombFilter(void) {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- flush
		lpf_state = 0.0;

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffer(sampleRate, bufferLength_mSec);

		return true;
	}

	/** process CombFilter */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		double yn = delay.readDelay();
		double input = 0.0;

		// --- form input & write
		if (combFilterParameters.enableLPF)
		{
			// --- apply simple 1st order pole LPF
			double g2 = lpf_g*(1.0 - comb_g); // see book for equation 11.27 (old book)
			double filteredSignal = yn + g2*lpf_state;
			input = xn + comb_g*(filteredSignal);
			lpf_state = filteredSignal;
		}
		else
		{
			input = xn + comb_g*yn;
		}

		delay.writeDelay(input);

		// --- done
		return yn;
	}

	/** return false: this object only processes samples */
	virtual bool canProcessAudioFrame() { return false; }

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return CombFilterParameters custom data structure
	*/
	CombFilterParameters getParameters()
	{
		return combFilterParameters;
	}

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param CombFilterParameters custom data structure
	*/
	void setParameters(const CombFilterParameters& params)
	{
		combFilterParameters = params;

		// --- update the delay line parameters first
		SimpleDelayParameters delayParams = delay.getParameters();
		delayParams.delayTime_mSec = combFilterParameters.delayTime_mSec;
		delayParams.interpolate = combFilterParameters.interpolate;
		delay.setParameters(delayParams); // this will set the delay time in samples

										  // --- calculate g with RT60 time (requires updated delay above^^)
		double exponent = -3.0*delayParams.delay_Samples*(1.0 / sampleRate);
		double rt60_mSec = combFilterParameters.RT60Time_mSec / 1000.0; // RT is in mSec!
		comb_g = pow(10.0, exponent / rt60_mSec);

		// --- set LPF g
		lpf_g = combFilterParameters.lpf_g;
	}

	/** create new buffers */
	void createDelayBuffer(double _sampleRate, double delay_mSec)
	{
		sampleRate = _sampleRate;
		bufferLength_mSec = delay_mSec;

		// --- create new buffer, will store sample rate and length(mSec)
		delay.createDelayBuffer(_sampleRate, delay_mSec);
	}

private:
	CombFilterParameters combFilterParameters; ///< object parameters
	double sampleRate = 0.0;	///< sample rate
	double comb_g = 0.0;		///< g value for comb filter (feedback, not %)
	double bufferLength_mSec = 0.0; ///< total buffer length

	// --- LPF support
	double lpf_g = 0.0;		///< LPF g value
	double lpf_state = 0.0; ///< state register for LPF (z^-1)

	// --- delay buffer of doubles
	SimpleDelay delay;		///< delay for comb filter
};
