#pragma once
#include <cassert>

#define M_LN2	   0.69314718055994530942

float clamp(float x, float a, float b) {
  return x < a? a: x > b? b: x;
  /*
  if(x < a) return a;
  if(x > b) return b;
  return x;
  */
}
struct ResonantLowpassFilter
{

  double fs,fc,q,g;
  double hp,bp,lp;
  double buf0,buf1;
  
  ResonantLowpassFilter(const double SampleRate,
                        const double cutoff,
                        const double res,
                        const double gain)
  { 
    fc = cutoff;
    fs = SampleRate;
    q  = clamp(res,0.0001,0.9995);
    g  = gain;                
  }

  ~ResonantLowpassFilter() {
    
  }
 
  void setCutoff(float c) {
    fc = c;
  }
  void setResonance(float r) {
    q  = r;
  }
  
  float Tick(float in, float A = 1, float Fc=0, float R=0) {            	
      //set feedback amount given f and q between 0 and 1
      float f  = 2.0*sin(2*M_PI*fc/fs);
      float fb = q + q/(1.0 - f);

      hp = in - buf0;
      bp = buf0 - buf1;
          
      buf0 = buf0 + f * (hp + fb*bp);
      buf1 = buf1 + f * (buf0 - buf1);
      lp   = buf1;  

      return lp;
  }

  void Process(size_t n, float * input, float * output) {    
    for(size_t i = 0; i < n; i++) output[i] = Tick(input[i]);
  }
};  
 