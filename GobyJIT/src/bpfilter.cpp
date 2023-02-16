/* biquadfilter.cpp
 * assignment 3 : Formant Filter
 * ECS7012 Music and Audio Programming
 * Student name : Selim Sheta
 * Student ID : 190295033
 * This code runs on the Bela embedded audio platform (bela.io).
 * It contains the functioh definitions for the lpfilter, hpfilter
 * and bpfilter classes declared in the header file bpfilter.h.
 */
 
#include <cmath>
#include <bpfilter.h>


/********************** DEFINTIONS FOR 4th ORDER LOW-PASS FILTER CLASS ***************************
 This Class is a 4th order low pass filter and is used in the the band-pass filter class below
 The filter is made of two cascaded 2nd order butterworth filters
*/
		
lpfilter::lpfilter() {}

// Calculate 2nd-order lowpass filter coefficients given specifications
void lpfilter::calculate_coefficients(float sampleRate, float frequency, float q)
{
	double t = 1/sampleRate;		// Convert sampling rate to sampling period
	double w = 2*M_PI*frequency;	// Convert cutoff frequency to rad
	w = (2/t) * tanf(w*t/2);		// pre-warp
	
	// Coefficients for output
	gCoefY[0] = 4*q + 2*t*w + q*pow(t*w,2); 			// a0
	gCoefY[1] = (-8*q + 2*q*pow(t*w,2))/gCoefY[0];		// a1
	gCoefY[2] = (4*q - 2*t*w + q*pow(t*w,2))/gCoefY[0];	// a2
	
	// Coefficients for input
	gCoefX[0] = (q*pow(t*w,2))/gCoefY[0];	// b0
	gCoefX[1] = (2*q*pow(t*w,2))/gCoefY[0];	// b1
	gCoefX[2] = gCoefX[0];	// b2
	
	// The coefficients are normalised by dividing by a0, then a0 is set to 1
	gCoefY[0] = 1; 
}

// This function takes an integer n which specifies the filter to update, and returns
// the output of that filter using the difference equation of the filter.
// n = 0 -> Update the first filter 
// n = 1 -> Update the second filter 
double lpfilter::update_filter(int n)
{
	// Using the difference equation to calculate the Normalised output
	gStateY[n][0] = (gCoefX[0]*gStateX[n][0] + gCoefX[1]*gStateX[n][1] + gCoefX[2]*gStateX[n][2] 
		- gCoefY[1]*gStateY[n][1] - gCoefY[2]*gStateY[n][2]);
		
	// updating the x registers
	gStateX[n][2] = gStateX[n][1];
	gStateX[n][1] = gStateX[n][0];
	// updating the y registers
	gStateY[n][2] = gStateY[n][1];
	gStateY[n][1] = gStateY[n][0];
	
	// Return the output of the filter
	return gStateY[n][0];
}

float lpfilter::process(float sample){
	// Pass signal to input of first filter
	gStateX[0][0] = sample;
	
	// Pass output of first filter to input of second filter
	gStateX[1][0] = update_filter(0);
	
	// Return the output of the second filter
	float out = update_filter(1);
	return out;
}

/********************** DEFINTIONS FOR 4th ORDER HIGH-PASS FILTER CLASS ***************************
 This Class is a 4th order high pass filter and is used in the the band-pass filter class below
 The filter is made of two cascaded 2nd order butterworth filters
*/
		
hpfilter::hpfilter() {}

// Calculate 2nd-order highpass filter coefficients given specifications
void hpfilter::calculate_coefficients(float sampleRate, float frequency, float q)
{
	double t = 1/sampleRate;		// Convert sampling rate to sampling period
	double w = 2*M_PI*frequency;	// Convert cutoff frequency to rad and pre-warp
	w = (2/t) * tanf(w*t/2);		// pre-warp
	
	// Coefficients for output
	gCoefY[0] = 4*q + 2*t*w + q*pow(t*w,2); 			// a0
	gCoefY[1] = (-8*q + 2*q*pow(t*w,2))/gCoefY[0];		// a1
	gCoefY[2] = (4*q - 2*t*w + q*pow(t*w,2))/gCoefY[0];	// a2
	
	// Coefficients for input
	gCoefX[0] = (4*q)/gCoefY[0];	// b0
	gCoefX[1] = (-8*q)/gCoefY[0];	// b1
	gCoefX[2] = gCoefX[0];	// b2
	
	// The coefficients are normalised by dividing by a0 then a0 is set to 1
	gCoefY[0] = 1; 
}

// This function takes an integer n which specifies the filter to update, and returns
// the output of that filter using the difference equation of the filter.
// n = 0 -> Update the first filter 
// n = 1 -> Update the second filter 
double hpfilter::update_filter(int n)
{
	// Using the difference equation to calculate the normalised output
	gStateY[n][0] = (gCoefX[0]*gStateX[n][0] + gCoefX[1]*gStateX[n][1] + gCoefX[2]*gStateX[n][2] 
		- gCoefY[1]*gStateY[n][1] - gCoefY[2]*gStateY[n][2]);
		
	// updating the x registers
	gStateX[n][2] = gStateX[n][1];
	gStateX[n][1] = gStateX[n][0];
	// updating the y registers
	gStateY[n][2] = gStateY[n][1];
	gStateY[n][1] = gStateY[n][0];
	
	// Return the output of the filter
	return gStateY[n][0];
}

float hpfilter::process(float sample){
	// Pass signal to input of first filter
	gStateX[0][0] = sample;
	
	// Pass output of first filter to input of second filter
	gStateX[1][0] = update_filter(0);
	
	// Return the output of the second filter
	float out = update_filter(1);
	return out;
}

/********************** DEFINTIONS FOR 4th ORDER BAND-PASS FILTER CLASS ***************************
 This Class is a 4th order band pass filter and uses the low-pass and high-pass filter classes above.
*/

bpfilter::bpfilter(){
	lowPassFilter = lpfilter(); 	// Low pass filter
	highPassFilter = hpfilter();	// High pass filter
	centerFrequency = 1000;			// Center frequency of the pass band
	q = 2;							// Q factor
	gaindB = 0;						// gain in dB
}

// This function updates the variables and the coefficients of the filters to match the specifications
void bpfilter::setParameters(float sampleRate, float frequency, float filterQ, float filterGaindB){
	q = filterQ;
	gaindB = filterGaindB;
	centerFrequency = frequency;
	lowPassFilter.calculate_coefficients(sampleRate,centerFrequency*(1+1/q),0.707);
	highPassFilter.calculate_coefficients(sampleRate,centerFrequency*(1-1/q),0.707);
}

// This function takes as input an audio sample and returned the processed version of that sample. 
// It also updates the registers of the low and high-pass filters.
float bpfilter::process(float sample){
	// Process input signal with high pass filter
	float xhp = highPassFilter.process(sample);
	// Process signal with low pass filter
	float out = lowPassFilter.process(xhp);
	// gain
	float g = pow(10,gaindB/20);
	// Return filter output
	return g*out;
}