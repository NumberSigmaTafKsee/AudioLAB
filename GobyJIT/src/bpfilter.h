/* bpfilter.h
 * assignment 3 : Formant Filter
 * ECS7012 Music and Audio Programming
 * Student name : Selim Sheta
 * Student ID : 190295033
 * This code runs on the Bela embedded audio platform (bela.io).
 * This is a header file containing the declaration for the
 * lpfilter, hpfilter and bpfilter classes, for use outside of the
 * class.
 */

#ifndef _BPFILTER_H
#define _BPFILTER_H

// Low-pass filter
class lpfilter
{
	private:
		// Filter coefficients
		double gCoefY[3] = {0,0,0};		// Filter coefficients for outputs
		double gCoefX[3] = {0,0,0};		// Filter coefficients for inputs
		
		// Filter states for y and x. First row is for filter 1, 
		// second row is for filter 2.
		double gStateY[2][3] = {{0,0,0},{0,0,0}};	// Filter states for outputs
		double gStateX[2][3] = {{0,0,0},{0,0,0}};	// Filter states for inputs
		
		double update_filter(int n);
	
	public:
		lpfilter();
		void calculate_coefficients(float sampleRate, float frequency, float q);
		float process(float sample);
};

// High-pass filter
class hpfilter
{
	private:
		// Filter coefficients
		double gCoefY[3] = {0,0,0};		// Filter coefficients for outputs
		double gCoefX[3] = {0,0,0};		// Filter coefficients for inputs
		
		// Filter states for y and x. First row is for filter 1, 
		// second row is for filter 2.
		double gStateY[2][3] = {{0,0,0},{0,0,0}};	// Filter states for outputs
		double gStateX[2][3] = {{0,0,0},{0,0,0}};	// Filter states for inputs
		
		double update_filter(int n);
		
	public:
		hpfilter();
		void calculate_coefficients(float sampleRate, float frequency, float q);
		float process(float sample);
};

// Band-pass filter
class bpfilter
{
	private: 
		lpfilter lowPassFilter; 	// low-pass filter
		hpfilter highPassFilter;	// high-pass filter
		float q;					// Q-factor
		float centerFrequency;		// Centre frequency of the pass band
		float gaindB;				// Gain in dB

	public:
		bpfilter();
		void setParameters(float sampleRate, float frequency, float filterQ, float filterGaindB);
		float process(float sample);
};


#endif