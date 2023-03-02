#pragma once
#include <cstring>
#include <cmath>

//==============================================================================
/**
   This class implements a LUT-based LFO with various waveforms and linear 
   interpolation.
   
   It uses 32-bit fixed-point phase and increment, where the 8 MSB
   represent the integer part of the number and the 24 LSB the fractionnal part
   
   @author		Remy Muller
   @date		20030822
   modified by macho charlie 1993
*/
//==============================================================================


class LFO
{
public:

  /** type enumeration used to select the different waveforms*/
  //typedef 
  typedef enum {triangle, sinus, sawtooth, square, exponent, kNumWave} waveform_t;
  
  /** phase type */
  typedef unsigned int phase_t;


  /**  @param samplerate the samplerate in Hz */
  LFO(DspFloatType samplerate);
  virtual ~LFO() {}

  /** increments the phase and outputs the new LFO value.
      @return the new LFO value between [-1;+1] */ 
  DspFloatType tick();
  void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out );
  void ProcessBlock(size_t n, DspFloatType * in, DspFloatType * out) {
    ProcessSIMD(n,in,out);
  }
  void ProcessInplace(size_t n, DspFloatType * in) {
    ProcessSIMD(n,in);
  }
  /** change the current rate
      @param rate new rate in Hz */
  void setRate(const DspFloatType rate);

  /** change the current samplerate
      @param samplerate new samplerate in Hz */
  void setSampleRate(const DspFloatType samplerate_) {samplerate = (samplerate_>0.0) ? samplerate : 44100.0f;}

  /** select the desired waveform for the LFO
      @param index tag of the waveform
   */
  void setWaveform(waveform_t index);

  /** @return the waveform's name as a std::string*/
  const std::string getWaveformName(long index){return waveNames[index];}

  /** @return the waveform's name as a C-string (char *) */
  const char * get_C_WaveformName(long index){return waveNames[index].c_str();}

private:
  /** names of the waveforms for display purpose*/
  static const std::string waveNames[kNumWave];

  DspFloatType samplerate;
  
  /** phase and phase increment
      the 8 Most Significant Bits represent the integer part,
      the 24 Least Significant Bits represent the fractionnal part.
      that way we can automatically adress the table with the MSB
      between 0-255 (which overflow automatically) and use the LSB 
      to determine the fractionnal part for linear interpolation with a precision of 
      \f[ 2^-24 \f] */
  phase_t phase,inc;
  
  /** table length is 256+1, with table[0] = table[256] 
      that way we can perform linear interpolation:
      \f[ val = (1-frac)*u[n] + frac*u[n+1] \f]
      even with n = 255.
      For n higher than 255, n is automatically  wrapped to 0-255*/
  DspFloatType table[257]; 
};


const std::string LFO::waveNames[] = {"triangle", "sinus", "sawtooth", "square", "exponent"};

inline LFO::LFO(DspFloatType samplerate)
  : samplerate(samplerate),
    phase(0),
    inc(0)
{
  setWaveform(LFO::triangle);   
  setRate(1.0f); //1Hz
}

const DspFloatType k1Div24lowerBits = 1.0f/(DspFloatType)(1<<24);

inline DspFloatType LFO::tick()
{
  // the 8 MSB are the index in the table in the range 0-255 
  int i = phase >> 24; 
  
  // and the 24 LSB are the fractionnal part
  DspFloatType frac = (phase & 0x00FFFFFF) * k1Div24lowerBits;

  // increment the phase for the next tick
  phase += inc; // the phase overflow itself

  return table[i]*(1.0f-frac) + table[i+1]*frac; // linear interpolation
}
void ProcessSIMD(size_t n, DspFloatType * in, DspFloatType * out ) {
  #pragma omp simd aligned(in,out)
  for(size_t i = 0; i < n; i++) {
    // the 8 MSB are the index in the table in the range 0-255 
    int i = phase >> 24;     
    // and the 24 LSB are the fractionnal part
    DspFloatType frac = (phase & 0x00FFFFFF) * k1Div24lowerBits;
    // increment the phase for the next tick
    phase += inc; // the phase overflow itself
    out[i] =  table[i]*(1.0f-frac) + table[i+1]*frac; // linear interpolation
    if(in) out[i] *= in[i];
  }
}
inline void LFO::setRate(DspFloatType rate)
{
  /** the rate in Hz is converted to a phase increment with the following formula
      \f[ inc = (256*rate/samplerate) * 2^24 \f] */
  inc =  (unsigned int)((256.0f * rate / samplerate) * (DspFloatType)(1<<24));
}

inline void LFO::setWaveform(waveform_t index)
{
  switch(index)
    {
    case sinus:
      {
    	DspFloatType pi = 4.0 * atan(1.0);

	int i;
  #pragma omp simd
	for(i=0;i<=256;i++)
	  table[i] = sin(2.0f*pi*(i/256.0f));

	break;
      }
    case triangle:
      {
	int i;
  #pragma omp simd
	for(i=0;i<64;i++)
	  {
	    table[i]     =        i / 64.0f;
	    table[i+64]  =   (64-i) / 64.0f;
	    table[i+128] =      - i / 64.0f;
	    table[i+192] = - (64-i) / 64.0f;
	  }
	table[256] = 0.0f;
	break;
      }
    case sawtooth:
      {
	int i;
  #pragma omp simd
	for(i=0;i<256;i++)
	  {
	    table[i] = 2.0f*(i/255.0f) - 1.0f;
	  }
	table[256] = -1.0f;
	break;
      }
    case square:
      {
	int i;
  #pragma omp simd
	for(i=0;i<128;i++)
	  {
	    table[i]     =  1.0f;
	    table[i+128] = -1.0f;
	  }
	table[256] = 1.0f;
	break;
      }
    case exponent:
      {
	/* symetric exponent similar to triangle */
	int i;
	DspFloatType e = (DspFloatType)exp(1.0f);
  #pragma omp simd
	for(i=0;i<128;i++)
	  {
	    table[i] = 2.0f * ((exp(i/128.0f) - 1.0f) / (e - 1.0f)) - 1.0f  ;
	    table[i+128] = 2.0f * ((exp((128-i)/128.0f) - 1.0f) / (e - 1.0f)) - 1.0f  ;
	  }
	table[256] = -1.0f;
	break;
      }
    default:
      {
	break;
      }
    } 
}
