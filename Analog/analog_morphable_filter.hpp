/*
  ==============================================================================

    MorphableFilter.h
    Created: 8 Oct 2021 4:20:46pm
    Author:  Eoin Roe

  ==============================================================================
*/

#pragma once


#include <cmath>
#include <array>
#include <vector>
#include <functional>
#include "GenericSoundObject.hpp"

namespace Analog::Filters
{
	template<typename DSP>
	class FirstOrderFilter {
	public:
		// Using member intializer list to initialize the reference member variables
		FirstOrderFilter(DSP& a1, DSP& b0, DSP& b1, DSP& g) : a1(a1), b0(b0), b1(b1), g(g) {}
		DSP process(DSP input);
	private:
		// Using reference_wrapper since it is more flexible
		// and allows the object to be copied if necessary.
		std::reference_wrapper<DSP> a1;
		std::reference_wrapper<DSP> b0;
		std::reference_wrapper<DSP> b1;
		std::reference_wrapper<DSP> g;
		
		// Make sure previous input and output are initialized.

		// x[n-1]
		DSP previousInput {};
		
		// y[n-1]
		DSP previousOutput {};
	};

	DSP FirstOrderFilter::process(DSP in)
	{
		DSP out = ((b0*in + b1*previousInput - a1*previousOutput) * g) + a1*previousOutput;

		previousInput = in;
		previousOutput = out;

		return out;
	}
	
	template<typename DSP>
	class MorphableFilter : public GSSoundProcessor<DSP> {
	public:
		MorphableFilter();
		DSP process(DSP input);
		void calculate_coefficients(DSP sampleRate, DSP frequencyHz, DSP resonance);
		void set_filter_type(int type);
	private:
		// Coefficients
		DSP a1 {};
		DSP b0 {};
		DSP b1 {};
		DSP g  {};
		
		std::vector<FirstOrderFilter<DSP>> filters;
		
		// Fixed by design
		const DSP gComp = 0.5f;
		DSP gRes = 0.75f;
		
		DSP delayed_filter_output {};

		enum filter_type
		{
			two_pole_low_pass,
			two_pole_band_pass,
			two_pole_high_pass,
			four_pole_low_pass,
			four_pole_band_pass,
			four_pole_high_pass
		};
						
		filter_type filter_t = two_pole_low_pass;
		
		std::array<std::array<DSP, 5>, 6> filter_weights;
	};



	template<typename DSP>
	MorphableFilter<DSP>::MorphableFilter() : GSSoundProcessor<DSP>()
	{
		// Create four first-order filters and add them to the vector
		for (int i = 0; i < 4; ++i)
			filters.push_back({a1, b0, b1, g});
		
		filter_weights[two_pole_low_pass]   = { {0,  0,  1,  0, 0} };
		filter_weights[two_pole_band_pass]  = { {0,  2, -2,  0, 0} };
		filter_weights[two_pole_high_pass]  = { {1, -2,  1,  0, 0} };
		
		filter_weights[four_pole_low_pass]  = { {0,  0,  0,  0, 1} };
		filter_weights[four_pole_band_pass] = { {0,  0,  4, -8, 4} };
		filter_weights[four_pole_high_pass] = { {1, -4,  6, -4, 1} };
	}

	template<typename DSP>
	void MorphableFilter<DSP>::set_filter_type(int type)
	{
		// Convert integer to the enumeration type
		filter_t = static_cast<filter_type>(type);
	}

	template<typename DSP>
	void MorphableFilter<DSP>::calculate_coefficients(DSP sampleRate, DSP frequencyHz, DSP resonance)
	{
		a1 = 1.0f;
		b0 = 1.0f / 1.3f;
		b1 = 0.3f / 1.3f;
		
		// Angular frequency
		DSP w =  2.0f * M_PI * frequencyHz / sampleRate;
		
		// Polynomial model which improves the fit
		g = 0.9892f*w - 0.4342f*powf(w, 2.0f) + 0.1381f*powf(w, 3.0f) - 0.0202f*powf(w, 4.0f);
		
		// Calculate a cubic polynomial to correct for the error
		gRes = resonance * (1.0029f + 0.0526f*w - 0.0926f*powf(w, 2.0f) + 0.0218f*powf(w, 3.0f));
	}

	template<typename DSP>
	DSP MorphableFilter<DSP>::process(DSP input)
	{
		DSP feedback_loop = input - ( (delayed_filter_output - (input * gComp)) * gRes * 4.0f );
		
		// Use the hyperbolic tangent function to approximate the nonlinearity of an analog circuit
		DSP nonlinearity = tanhf (feedback_loop);

		DSP filter_out {};
		DSP summing_junction {};
		
		// Use filter_type to index into the two dimensional array
		for (unsigned int i = 0; i < filters.size(); ++i) {
			if (i == 0) {
				summing_junction += nonlinearity * filter_weights[filter_t][i];
				filter_out = filters[i].process(nonlinearity);
			}
			else {
				summing_junction += filter_out * filter_weights[filter_t][i];
				filter_out = filters[i].process(filter_out);
			}
		}
		
		// This happens just before the multiplication by the coefficient E
		delayed_filter_output = filter_out;
		
		// Add the final weighting of the output outside of the for loop
		summing_junction += filter_out * filter_weights[filter_t][4];
		
		return summing_junction;
	}
}
