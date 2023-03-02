#pragma once

// https://github.com/janne808/kocmoc-rack-modules/blob/master/src/svfilter.cpp

#include "KocMocIIR.hpp"


// filter modes
enum SVFFilterMode {
   SVF_LOWPASS_MODE,
   SVF_BANDPASS_MODE,
   SVF_HIGHPASS_MODE
};

// integration methods
enum SVFIntegrationMethod {
   SVF_SEMI_IMPLICIT_EULER,
   SVF_PREDICTOR_CORRECTOR,
   SVF_TRAPEZOIDAL,
   SVF_INV_TRAPEZOIDAL
};

class SVFilter{
public:
  // constructor/destructor
  SVFilter(double newCutoff, double newResonance, int newOversamplingFactor,
	   SVFFilterMode newFilterMode, double newSampleRate,
	   SVFIntegrationMethod newIntegrationMethod, int newDecimatorOrder);
  SVFilter();
  ~SVFilter();

  // set filter parameters
  void SetFilterCutoff(double newCutoff);
  void SetFilterResonance(double newResonance);
  void SetFilterMode(SVFFilterMode newFilterMode);
  void SetFilterSampleRate(double newSampleRate);
  void SetFilterIntegrationMethod(SVFIntegrationMethod method);
  void SetFilterOversamplingFactor(int newOversamplingFactor);
  void SetFilterDecimatorOrder(int decimatorOrder);
    
  // get filter parameters
  double GetFilterCutoff();
  double GetFilterResonance();
  SVFFilterMode GetFilterMode();  
  double GetFilterSampleRate();
  SVFIntegrationMethod GetFilterIntegrationMethod();
  int GetFilterOversamplingFactor();  
  int GetFilterDecimatorOrder();
  
  // tick filter state
  void filter(double input);

  // get filter responses
  double GetFilterLowpass();
  double GetFilterBandpass();
  double GetFilterHighpass();

  // get filter output
  double GetFilterOutput();

  // reset state
  void ResetFilterState();
  
private:
  // set integration rate
  void SetFilterIntegrationRate();

  // pade approximant functions for hyperbolic functions
  // filter parameters
  double cutoffFrequency;
  double Resonance;
  SVFFilterMode filterMode;
  SVFIntegrationMethod integrationMethod;
  double dt;
  double sampleRate;
  int oversamplingFactor;
  int decimatorOrder;
  
  // filter state
  double lp;
  double bp;
  double hp;
  double u_t1;
  
  // filter output
  double out;

  // IIR downsampling filter
  IIRLowpass *iir;
};




// steepness of downsample filter response
#define IIR_DOWNSAMPLE_ORDER 16

// downsampling passthrough bandwidth
#define IIR_DOWNSAMPLING_BANDWIDTH 0.9

// maximum newton-raphson iteration steps
#define SVF_MAX_NEWTON_STEPS 8

// check for newton-raphson breaking limit
#define SVF_NEWTON_BREAKING_LIMIT 1

// constructor
SVFilter::SVFilter(double newCutoff, double newResonance, int newOversamplingFactor,
		   SVFFilterMode newFilterMode, double newSampleRate,
		   SVFIntegrationMethod newIntegrationMethod, int newDecimatorOrder){
  // initialize filter parameters
  cutoffFrequency = newCutoff;
  Resonance = newResonance;
  filterMode = newFilterMode;
  sampleRate = newSampleRate;
  oversamplingFactor = newOversamplingFactor;
  decimatorOrder = newDecimatorOrder;

  SetFilterIntegrationRate();

  // initialize filter state
  hp = bp = lp = out = u_t1 = 0.0;
  
  integrationMethod = newIntegrationMethod;
  
  // instantiate downsampling filter
  iir = new IIRLowpass(sampleRate * oversamplingFactor,
		       IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0,
		       decimatorOrder);
}

// default constructor
SVFilter::SVFilter(){
  // initialize filter parameters
  cutoffFrequency = 0.25;
  Resonance = 0.5;
  filterMode = SVF_LOWPASS_MODE;
  sampleRate = 44100.0;
  oversamplingFactor = 2;
  decimatorOrder = IIR_DOWNSAMPLE_ORDER;
  
  SetFilterIntegrationRate();
  
  // initialize filter state
  hp = bp = lp = out = u_t1 = 0.0;
  
  integrationMethod = SVF_TRAPEZOIDAL;
  
  // instantiate downsampling filter
  iir = new IIRLowpass(sampleRate * oversamplingFactor, IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0, decimatorOrder);
}

// default destructor
SVFilter::~SVFilter(){
  delete iir;
}

void SVFilter::ResetFilterState(){
  // initialize filter parameters
  cutoffFrequency = 0.25;
  Resonance = 0.5;

  SetFilterIntegrationRate();
  
  // initialize filter state
  hp = bp = lp = out = u_t1 = 0.0;
  
  // set oversampling
  iir->SetFilterSamplerate(sampleRate * oversamplingFactor);
  iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
  iir->SetFilterOrder(decimatorOrder);
}

void SVFilter::SetFilterCutoff(double newCutoff){
  cutoffFrequency = newCutoff;

  SetFilterIntegrationRate();
}

void SVFilter::SetFilterResonance(double newResonance){
  Resonance = newResonance;
}

void SVFilter::SetFilterMode(SVFFilterMode newFilterMode){
  filterMode = newFilterMode;
}

void SVFilter::SetFilterSampleRate(double newSampleRate){
  sampleRate = newSampleRate;
  iir->SetFilterSamplerate(sampleRate * (double)(oversamplingFactor));
  iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
  iir->SetFilterOrder(decimatorOrder);

  SetFilterIntegrationRate();
}

void SVFilter::SetFilterIntegrationMethod(SVFIntegrationMethod method){
  integrationMethod = method;
  ResetFilterState();
}

void SVFilter::SetFilterOversamplingFactor(int newOversamplingFactor){
  oversamplingFactor = newOversamplingFactor;
  iir->SetFilterSamplerate(sampleRate * oversamplingFactor);
  iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
  iir->SetFilterOrder(decimatorOrder);

  SetFilterIntegrationRate();
}

void SVFilter::SetFilterDecimatorOrder(int newDecimatorOrder){
  decimatorOrder = newDecimatorOrder;
  iir->SetFilterOrder(decimatorOrder);
}

void SVFilter::SetFilterIntegrationRate(){
  // normalize cutoff freq to samplerate
  dt = 44100.0 / (sampleRate * (double)(oversamplingFactor)) * cutoffFrequency;

  // clamp integration rate
  if(dt < 0.0){
    dt=0.0;
  }
}

double SVFilter::GetFilterCutoff(){
  return cutoffFrequency;
}

double SVFilter::GetFilterResonance(){
  return Resonance;
}

double SVFilter::GetFilterOutput(){
  return out;
}

SVFFilterMode SVFilter::GetFilterMode(){
  return filterMode;
}

double SVFilter::GetFilterSampleRate(){
  return sampleRate;
}

int SVFilter::GetFilterOversamplingFactor(){
  return oversamplingFactor;
}

int SVFilter::GetFilterDecimatorOrder(){
  return decimatorOrder;
}

SVFIntegrationMethod SVFilter::GetFilterIntegrationMethod(){
  return integrationMethod;
}

void SVFilter::filter(double input){
  // noise term
  double noise;

  // feedback amount variables
  double fb = 1.0 - (3.5*Resonance);

  // integration rate
  double dt2 = dt;
  
  // update noise terms
  noise = static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
  noise = 1.0e-6 * 2.0 * (noise - 0.5);

  input += noise;

  // clamp integration rate
  switch(integrationMethod){
  case SVF_TRAPEZOIDAL:
    if(dt2 > 0.65){
      dt2 = 0.65;
    }
    break;
  case SVF_INV_TRAPEZOIDAL:
    if(dt2 > 1.0){
      dt2 = 1.0;
    }
    break;
  default:
    if(dt2 > 0.25){
      dt2 = 0.25;
    }
    break;
  }
  
  // integrate filter state
  // with oversampling
  for(int nn = 0; nn < oversamplingFactor; nn++){
    // switch integration method
    switch(integrationMethod){
    case SVF_SEMI_IMPLICIT_EULER:
      {
	// loss factor
	double beta = 1.0 - (0.0075/oversamplingFactor);

       	hp = input - lp - fb*bp - SinhPade54(bp);
	bp += dt2*hp;
	bp *= beta;
	lp += dt2*bp;
      }
      break;
    case SVF_TRAPEZOIDAL:
      // trapezoidal integration
      {
	double alpha = dt2/2.0;
	double beta = 1.0 - (0.0075/oversamplingFactor);
	double alpha2 = dt2*dt2/4.0 + fb*alpha;
	double D_t = (1.0 - dt2*dt2/4.0)*bp +
	              alpha*(u_t1 + input - 2.0*lp - fb*bp - SinhPade54(bp));
	double x_k, x_k2;

	// starting point is last output
	x_k = bp;
	
	// newton-raphson
	for(int ii=0; ii < SVF_MAX_NEWTON_STEPS; ii++) {
	  x_k2 = x_k - (x_k + alpha*SinhPade54(x_k) + alpha2*x_k - D_t)/
	                  (1.0 + alpha*CoshPade54(x_k) + alpha2);

#ifdef SVF_NEWTON_BREAKING_LIMIT
	  // breaking limit
	  if(abs(x_k2 - x_k) < 1.0e-9) {
	    x_k = x_k2;
	    break;
	  }
#endif
	  x_k = x_k2;
	}

	lp += alpha*bp;
	bp = beta*x_k;
	lp += alpha*bp;
      	hp = input - lp - fb*bp;
      }
      break;
    case SVF_INV_TRAPEZOIDAL:
      // inverse trapezoidal integration
      {
	double alpha = dt2/2.0;
	double beta = 1.0 - (0.0075/oversamplingFactor);
	double alpha2 = dt2*dt2/4.0 + fb*alpha;
	double D_t = (1.0 - dt2*dt2/4.0)*bp +
	              alpha*(u_t1 + input - 2.0*lp - fb*bp - sinh(bp));
	double y_k, y_k2;

	// starting point is last output
	y_k = sinh(bp);
	
	// newton-raphson
	for(int ii=0; ii < SVF_MAX_NEWTON_STEPS; ii++) {
	  y_k2 = y_k - (alpha*y_k + ASinhPade54(y_k)*(1.0 + alpha2) - D_t)/
	                  (alpha + (1.0 + alpha2)*dASinhPade54(y_k));

#ifdef SVF_NEWTON_BREAKING_LIMIT
	  // breaking limit
	  if(abs(y_k2 - y_k) < 1.0e-9) {
	    y_k = y_k2;
	    break;
	  }
#endif
	  
	  y_k = y_k2;
	}

     	lp += alpha*bp;
	bp = beta*asinh(y_k);
	lp += alpha*bp;
      	hp = input - lp - fb*bp;
      }
      break;
    default:
      break;
    }
    
    switch(filterMode){
    case SVF_LOWPASS_MODE:
      out = lp;
      break;
    case SVF_BANDPASS_MODE:
      out = bp;
      break;
    case SVF_HIGHPASS_MODE:
      out = hp;
      break;
    default:
      out = 0.0;
    }
    
    // downsampling filter
    if(oversamplingFactor > 1){
      out = iir->IIRfilter(out);
    }
  }
  
  // set input at t-1
  u_t1 = input;    
}

double SVFilter::GetFilterLowpass(){
  return lp;
}

double SVFilter::GetFilterBandpass(){
  return bp;
}

double SVFilter::GetFilterHighpass(){
  return hp;
}


