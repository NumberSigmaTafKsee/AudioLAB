#pragma once
//helper functions
//-------------------------------------------------------------------------
// numSamplesFromMSf :
// Determine the number of samples from the number of milliseconds delay
// passed to function
//-------------------------------------------------------------------------
inline float numSamplesFromMSf(const int sr, const float d_ms){
    return sr * d_ms * .001;
}

//----------------------------------------------------------------------
//  Linear Interpolation Function
//  x1  :  weighting 1
//  x2  :  weighting 2
//  y1  :  output y(n) at read index in buffer
//  y2  :  output y(n-1) at read index minus 1 in buffer
//  x   :  fractional value between samples
//----------------------------------------------------------------------
inline float linInterp(float x1, float x2, float y1, float y2, float x){
   	float denom = x2 - x1;
	if(denom == 0)
		return y1; // should not ever happen

	// calculate decimal position of x
	float dx = (x - x1)/(x2 - x1);

	// use weighted sum method of interpolating
	float result = dx*y2 + (1-dx)*y1;

	return result; 
}

//helper functions
//-------------------------------------------------------------------------
// numSamplesFromMSf :
// Determine the number of samples from the number of milliseconds delay
// passed to function
//-------------------------------------------------------------------------
inline float tapNumSamplesFromMSf(const int sr, const float d_ms){
    return sr * d_ms;
}

//helper functions
inline float calcCombGain(const float d_ms, const float rt60){
    return pow(10.0, ((-3.0 * d_ms) / (rt60 * 1000.0)));
    //(((return pow(10.0, (-3.0 * (d_ms*44100*.001) * (1.0/44100)) / (rt60 * 1000.0));

}

inline float calcAPGain(const float d_ms, const float rt60){
    return pow(10.0, ((-3.0 * d_ms) / (rt60 * 1000.0)));
}
