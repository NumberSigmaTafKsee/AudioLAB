/*
MoogLadder Filter class based on an implementation by Aaron Krajeski
This code is Unlicensed (i.e. public domain)
Aaron Krajeski has stated in an email: "That work is under no copyright. You may use it however you might like."
Source: http://song-swap.com/MUMT618/aaron/Presentation/demo.html
*/

// moogladder.h: header file for defining the moog ladder filter

// Thermal voltage (26 milliwats at room temperature)

#include <cmath>
#include "GenericSoundProcessor.hpp"

namespace Analog::Flters::Moog::MoogLadder
{
	#define VT 0.026
	template<typename DSP>
	class MoogLadder : public GSSoundProcessor<DSP> {

	public:
		// Constructor
		MoogLadder();
		
		// Constructor specifying a sample rate
		MoogLadder(DSP sampleRate);
		
		// Set the sample rate
		void setSampleRate(DSP rate);
		
		// Set the cutoff frequency 
		void setFrequency(DSP cutoff);
		
		// Set the filter Q
		void setQ(DSP resonance);
		
		void setDrive(DSP drive);
		
		// Reset previous history of filter
		void reset();
		
		enum
        {
            PORT_CUTOFF,
            PORT_RESONANCE,			
			PORT_DRIVE,
			PORT_RESET,
        };
        void setPort(int port, DSP v)
        {
            switch (port)
            {
            case PORT_CUTOFF:
                setFrequency(v);
                break;
            case PORT_RESONANCE:
                setQ(v);
                break;
			case PORT_DRIVE:
				setDrive(v);
				break;
			case PORT_RESET:
				reset();
				break;
            }
        }
		// Calculate the next sample of output, changing the envelope
		// state as needed
		DSP Tick(DSP input, DSP A=1, DSP X=1, DSP Y=1); 
		
		void ProcessSIMD(size_t n, DSP * in, DSP * out);
		
		void ProcessBlock(size_t n, DSP * in, DSP * out) {
			ProcessSIMD(n,in,out);
		}
		void ProcessInplace(size_t n, DSP * p) {
			ProcessSIMD(n,p,p);
		}
		
		// Destructor
		~MoogLadder();

	private:

		// State variables, not accessible to the outside world
		DSP sampleRate_;	// Filter sample rate
		DSP cutoff_;		// Filter cutoff frequency
		DSP resonance_;	// Filter resonance
		DSP drive_;	// Filter drive.
		DSP gComp_;	// Compensation factor, used to decide how much of the input signal is added into the feedback loop.
		
		DSP state[5];// Array for holding the output of each stage of the ladder
		DSP delay[5];// Array for holding the delayed components
		DSP wc;		// The angular frequency of the cutoff.
		DSP g;		// A derived parameter for the cutoff frequency
		DSP gRes;	// A similar derived parameter for resonance.

		
	};


	// Constructor
	MoogLadder::MoogLadder() : MoogLadder(44100.0) {}

	// Constructor specifying a sample rate
	MoogLadder::MoogLadder(DSP sampleRate)
	{
		// Set the samplerate and inital values for drive, compensation, cutoff and q/resonance
		setSampleRate(sampleRate);
		for(int i = 0; i < 5; i++){
			state[i] = 0;
			delay[i] = 0;
		}
		
		drive_ = 1.0;
		gComp_ = 0.5;
		
		setFrequency(1000.0f);
		setQ(0.10f);
	}
		
	// Set the sample rate, used for all calculations
	void MoogLadder::setSampleRate(DSP rate)
	{
		sampleRate_ = rate;	
	}

	// Set the cutoff frequency 
	void MoogLadder::setFrequency(DSP c)
	{
		//Set the filter cutoff
		cutoff_ = c;
		
		// Calculate angular frequency
		wc = 2 * M_PI * cutoff_ / sampleRate_;
		
		// Calculate Huovilainen derived g parameter controlled by cutoff
		g = 0.9892 * wc - 0.4342 * std::pow(wc, 2) + 0.1381 * std::pow(wc, 3) - 0.0202 * std::pow(wc, 4);
	}
		
	// Set the Q and recalculate the coefficients
	void MoogLadder::setQ(DSP resonance)
	{
		// Set the filter resonance
		resonance_ = resonance;
		// Calculate derived resonance value
		gRes = resonance_ * (1.0029 + 0.0526 * wc - 0.926 * std::pow(wc, 2) + 0.0218 * std::pow(wc, 3));
	}

	// Set the sample rate, used for all calculations
	void MoogLadder::setDrive(DSP drive)
	{
		drive_ = drive;	
	}
		
	// Calculate the next sample of output
	DSP MoogLadder::Tick(DSP input, DSP A=1, DSP X=1, DSP Y=1); 
	{
		Undenormal denormals;
		//Next input to the filter is a combination of the last output and current sample, 
		// scaled by the resonance, filter drive and compensation value
		state[0] = tanh(drive_ * (input - 4 * gRes * (state[4] - gComp_ * input)));
		
		// Loop through each pole of the ladder filter		
		for(int i = 0; i < 4; i++)
		{
			// Equation implementing the Huovilainen (2006) one pole circuit diagram
			state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
			
			// Add output to feedback loop
			delay[i] = state[i];
		}
		
		// Filtered sample is the output of the final pole
		DSP output = state[4];
		
		return output;
	}
	void MoogLadder::ProcessSIMD(size_t n, DSP * in, DSP * out) {
		Undenormal denormals;
		#pragma omp simd aligned(in,out)
		for(size_t i = 0; i < n; i++)
		{
			//Next input to the filter is a combination of the last output and current sample, 
			// scaled by the resonance, filter drive and compensation value
			state[0] = std::tanh(drive_ * (in[i] - 4 * gRes * (state[4] - gComp_ * input)));
			
			// Loop through each pole of the ladder filter		
			for(int i = 0; i < 4; i++)
			{
				// Equation implementing the Huovilainen (2006) one pole circuit diagram
				state[i+1] = g * (0.3 / 1.3 * state[i] + 1 / 1.3 * delay[i] - state[i + 1]) + state[i + 1];
				
				// Add output to feedback loop
				delay[i] = state[i];
			}						
			out[i] =state[4];
		}		
	}
	// Reset the filter poles and feedback state
	void MoogLadder::reset()
	{
		for(int i = 0; i < 5; i++){
			state[i] = 0;
			delay[i] = 0;
		}	
	}
		
	// Destructor
	MoogLadder::~MoogLadder()
	{
		// Nothing to do here
	}
}
#undef VT
