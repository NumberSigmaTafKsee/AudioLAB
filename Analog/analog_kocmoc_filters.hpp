#pragma once
#include <cmath>
#include <cstdlib>
#include <cmath>
#include "StdSamples/Core/core_fastmath.hpp"


#define IIR_MAX_ORDER 32

// steepness of downsample filter response
#define IIR_DOWNSAMPLE_ORDER 16

// downsampling passthrough bandwidth
#define IIR_DOWNSAMPLING_BANDWIDTH 0.9

// maximum newton-raphson iteration steps
#define SKF_MAX_NEWTON_STEPS 8

// check for newton-raphson breaking limit
#define SKF_NEWTON_BREAKING_LIMIT 1

namespace Analog::KocMoc
{
  class IIRLowpass
  {
  public:
    // constructor/destructor
    IIRLowpass(DspFloatType newSamplerate, DspFloatType newCutoff, int newOrder);
    IIRLowpass();
    ~IIRLowpass();

    // set filter parameters
    void SetFilterOrder(int newOrder);
    void SetFilterSamplerate(DspFloatType newSamplerate);
    void SetFilterCutoff(DspFloatType newCutoff);

    // initialize biquad cascade delayline
    void InitializeBiquadCascade();
    
    // IIR filter signal 
    DspFloatType IIRfilter(DspFloatType input);

    // get filter coefficients
    DspFloatType* GetFilterCoeffA1();
    DspFloatType* GetFilterCoeffA2();
    DspFloatType* GetFilterCoeffK();
    
  private:
    // compute biquad cascade coefficients
    void ComputeCoefficients();

    // filter design variables
    DspFloatType samplerate;
    DspFloatType cutoff;
    int order;
    
    // dsp variables
    DspFloatType *a1;
    DspFloatType *a2;
    DspFloatType *K;
    DspFloatType *pa_real;
    DspFloatType *pa_imag;
    DspFloatType *p_real;
    DspFloatType *p_imag;
    
    // cascaded biquad buffers
    DspFloatType *z;
  };


  // constructor
  IIRLowpass::IIRLowpass(DspFloatType newSamplerate, DspFloatType newCutoff, int newOrder)
  {
    // initialize filter design parameters
    samplerate = newSamplerate;
    cutoff = newCutoff;
    order = newOrder;

    // allocate dsp vectors
    a1 = new DspFloatType[IIR_MAX_ORDER/2];
    a2 = new DspFloatType[IIR_MAX_ORDER/2];
    K = new DspFloatType[IIR_MAX_ORDER/2];
    pa_real = new DspFloatType[IIR_MAX_ORDER/2];
    pa_imag = new DspFloatType[IIR_MAX_ORDER/2];
    p_real = new DspFloatType[IIR_MAX_ORDER/2];
    p_imag = new DspFloatType[IIR_MAX_ORDER/2];

    // allocate cascaded biquad buffer
    z = new DspFloatType[IIR_MAX_ORDER];
    
    // initialize cascade delayline
    InitializeBiquadCascade();
    
    // compute impulse response
    ComputeCoefficients();
  }

  // default constructor
  IIRLowpass::IIRLowpass()
  {
    // set default design parameters
    samplerate=(DspFloatType)(44100.0);
    cutoff=(DspFloatType)(440.0);
    order=IIR_MAX_ORDER;
    
    // allocate dsp vectors
    a1 = new DspFloatType[IIR_MAX_ORDER/2];
    a2 = new DspFloatType[IIR_MAX_ORDER/2];
    K = new DspFloatType[IIR_MAX_ORDER/2];
    pa_real = new DspFloatType[IIR_MAX_ORDER/2];
    pa_imag = new DspFloatType[IIR_MAX_ORDER/2];
    p_real = new DspFloatType[IIR_MAX_ORDER/2];
    p_imag = new DspFloatType[IIR_MAX_ORDER/2];

    // allocate cascaded biquad buffer
    z = new DspFloatType[IIR_MAX_ORDER];
    
    // initialize cascade delayline
    InitializeBiquadCascade();
    
    // compute impulse response
    ComputeCoefficients();
  }

  // destructor
  IIRLowpass::~IIRLowpass(){
    // free dsp vectors
    delete[] a1;
    delete[] a2;
    delete[] K;
    delete[] pa_real;
    delete[] pa_imag;
    delete[] p_real;
    delete[] p_imag;
    
    // free cascaded biquad buffer
    delete[] z;
  }

  void IIRLowpass::SetFilterOrder(int newOrder){
    if(newOrder > IIR_MAX_ORDER){
      order = IIR_MAX_ORDER;
    }
    else{
      order = newOrder;
    }
    
    // initialize cascade delayline
    InitializeBiquadCascade();
    
    // compute new impulse response
    ComputeCoefficients();
  }

  void IIRLowpass::SetFilterSamplerate(DspFloatType newSamplerate){
    samplerate = newSamplerate;

    // initialize cascade delayline
    InitializeBiquadCascade();
    
    // compute new cascade coefficients
    ComputeCoefficients();
  }

  void IIRLowpass::SetFilterCutoff(DspFloatType newCutoff){
    cutoff = newCutoff;

    // initialize cascade delayline
    InitializeBiquadCascade();
    
    // compute new cascade coefficients
    ComputeCoefficients();
  }

  void IIRLowpass::InitializeBiquadCascade(){
    for(int ii=0; ii<order/2; ii++){
      z[ii*2+1] = 0.0;
      z[ii*2] = 0.0;
    }
  }

  DspFloatType IIRLowpass::IIRfilter(DspFloatType input){
    DspFloatType out=input;
    DspFloatType in;

    // process biquad cascade
    for(int ii=0; ii<order/2; ii++) {
      // compute biquad input
      in = K[ii]*out - a1[ii]*z[ii*2] - a2[ii]*z[ii*2+1];
        
      // compute biquad output
      out = in + 2.0*z[ii*2] + z[ii*2+1];
      
      // update delays
      z[ii*2+1] = z[ii*2];
      z[ii*2] = in;
    }
    
    return out;
  }

  DspFloatType* IIRLowpass::GetFilterCoeffA1(){
    return a1;
  }

  DspFloatType* IIRLowpass::GetFilterCoeffA2(){
    return a2;
  }

  DspFloatType* IIRLowpass::GetFilterCoeffK(){
    return K;
  }

  void IIRLowpass::ComputeCoefficients(){
    // place butterworth style analog filter poles
    DspFloatType theta;

    for(int ii = 0; ii<order/2; ii++) {
      int k = order/2 - ii;
      theta = (2.0*(DspFloatType)(k) - 1.0)*M_PI/(2.0*(DspFloatType)(order));
      
      pa_real[ii] = -1.0*sin(theta);
      pa_imag[ii] = cos(theta);
    }

    // prewarp and scale poles
    DspFloatType Fc = samplerate/M_PI*tan(M_PI*cutoff/samplerate);  
    for(int ii = 0; ii<order/2; ii++) {
      pa_real[ii] *= 2.0*M_PI*Fc; 
      pa_imag[ii] *= 2.0*M_PI*Fc; 
    }

    // bilinear transform to z-plane
    for(int ii = 0; ii<order/2; ii++) {
      // complex division
      DspFloatType u = (2.0*samplerate+pa_real[ii])/(2.0*samplerate); 
      DspFloatType v = pa_imag[ii]/(2.0*samplerate); 
      DspFloatType x = (2.0*samplerate-pa_real[ii])/(2.0*samplerate); 
      DspFloatType y = -1.0*pa_imag[ii]/(2.0*samplerate);
      
      DspFloatType c = 1.0/(x*x + y*y);
      
      p_real[ii] = c*(u*x + v*y);
      p_imag[ii] = c*(v*x - u*y);
    }
    
    // compute cascade coefficients
    for(int ii = 0; ii<order/2; ii++) {
      a1[ii] = -2.0*p_real[ii];
      a2[ii] = p_real[ii]*p_real[ii] + p_imag[ii]*p_imag[ii];
      K[ii] = (1.0 + a1[ii] + a2[ii])/4.0;
    }
  }

  // https://github.com/janne808/kocmoc-rack-modules/blob/master/src/fir.cpp
  class FIRLowpass{
  public:
    // constructor/destructor
    FIRLowpass(DspFloatType newSamplerate, DspFloatType newCutoff, int newOrder);
    FIRLowpass();
    ~FIRLowpass();

    // set filter parameters
    void SetFilterOrder(int newOrder);
    void SetFilterSamplerate(DspFloatType newSamplerate);
    void SetFilterCutoff(DspFloatType newCutoff);
    
    // FIR filter signal 
    DspFloatType FIRfilter(DspFloatType input);

    // get impulse response
    DspFloatType* GetImpulseResponse();

    // clean filter ring buffer
    void InitializeRingbuffer();
    
  private:
    // compute windowed ideal impulse response
    void ComputeImpulseResponse();

    // filter design variables
    DspFloatType samplerate;
    DspFloatType cutoff;
    int order;
    
    // dsp variables
    DspFloatType omega_c;
    DspFloatType *h_d;
    DspFloatType *h;
    DspFloatType *w;

    // ring buffer delay line
    DspFloatType *ringBuffer;
    int ringBufferIndex;
  };

  // constructor
  FIRLowpass::FIRLowpass(DspFloatType newSamplerate, DspFloatType newCutoff, int newOrder)
  {
    // initialize filter design parameters
    samplerate = newSamplerate;
    cutoff = newCutoff;
    order = newOrder;
    
    // allocate dsp vectors
    h_d = new DspFloatType[order];
    h = new DspFloatType[order];
    w = new DspFloatType[order];

    // initialize ring buffer delay line
    InitializeRingbuffer();

    // compute impulse response
    ComputeImpulseResponse();
  }

  // default constructor
  FIRLowpass::FIRLowpass()
  {
    // set default design parameters
    samplerate=(DspFloatType)(44100.0);
    cutoff=(DspFloatType)(440.0);
    order=128;
    
    // allocate dsp vectors
    h_d = new DspFloatType[order];
    h = new DspFloatType[order];
    w = new DspFloatType[order];

    // initialize ring buffer delay line
    InitializeRingbuffer();

    // compute impulse response
    ComputeImpulseResponse();
  }

  // destructor
  FIRLowpass::~FIRLowpass(){
    // free dsp vectors
    delete[] h_d;
    delete[] h;
    delete[] w;
  }

  void FIRLowpass::SetFilterOrder(int newOrder){
    order = newOrder;

    // free dsp vectors
    delete[] h_d;
    delete[] h;
    delete[] w;
    
    // allocate dsp vectors
    h_d = new DspFloatType[order];
    h = new DspFloatType[order];
    w = new DspFloatType[order];

    // compute new impulse response
    ComputeImpulseResponse();

    // initialize ring buffer delay line
    InitializeRingbuffer();
  }

  void FIRLowpass::SetFilterSamplerate(DspFloatType newSamplerate){
    samplerate = newSamplerate;

    // compute new impulse response
    ComputeImpulseResponse();

    // initialize ring buffer delay line
    InitializeRingbuffer();
  }

  void FIRLowpass::SetFilterCutoff(DspFloatType newCutoff){
    cutoff = newCutoff;

    // compute new impulse response
    ComputeImpulseResponse();

  }

  DspFloatType* FIRLowpass::GetImpulseResponse(){
    return h;
  }

  void FIRLowpass::InitializeRingbuffer(){
    // initialize ring buffer delay line
    ringBufferIndex = 0;
    ringBuffer = new DspFloatType[order];

    for(int n=0; n<order; n++){
      ringBuffer[n] = 0.0;
    }
  }

  DspFloatType FIRLowpass::FIRfilter(DspFloatType input){
    // update delay line
    ringBuffer[ringBufferIndex++] = input;
    if(ringBufferIndex > (order - 1)){
      ringBufferIndex -= order;
    }

    // compute convolution
    DspFloatType output = 0.0;
    int ii;
    for(int n = 0; n < order; n++){
      // wrap around index
      ii = (ringBufferIndex - (n + 1));
      if(ii < 0)
        ii += order;
      
      // multiply and accumulate
      output += h[n] * ringBuffer[ii]; 
    }

    return output;
  }

  void FIRLowpass::ComputeImpulseResponse(){
    // index as -M..M
    DspFloatType ii;
    
    // set cutoff frequency in radians
    omega_c = cutoff/samplerate;

    // compute truncated ideal impulse response
    for(int n=0; n<order; n++){
      // compute index as -M..M and avoid NaN at impulse peak
      ii = (DspFloatType)(n) - 1.0 - (DspFloatType)(floor((DspFloatType)(order)/2.0)) + 1.0e-9;

      // sample sinc function
      h_d[n] = std::sin((DspFloatType)(2.0 * M_PI * omega_c * ii))/(DspFloatType)(2.0 * M_PI * omega_c * ii);
    }

    // compute windowing function
    for(int n=0; n<order; n++){
      // compute index as -M..M and avoid NaN at impulse peak
      ii = (DspFloatType)(n) - 1.0 - (DspFloatType)(floor((DspFloatType)(order) / 2.0)) + 1.0e-9;

      // hanning window function
      w[n] = std::cos(M_PI * ii/(DspFloatType)(order));
      w[n] *= w[n];
    }

    // compute windowed ideal impulse function
    for(int n=0; n<order; n++){
      // window truncated ideal impulse response
      h[n] = w[n] * h_d[n];
    }
  }

  /////////////////////////////////////////////////////
  // Moog Ladder
  /////////////////////////////////////////////////////

  // filter modes
  enum LadderFilterMode {
    LADDER_LOWPASS_MODE,
    LADDER_BANDPASS_MODE,
    LADDER_HIGHPASS_MODE
  };

  // integration methods
  enum LadderIntegrationMethod {
    LADDER_EULER_FULL_TANH,
    LADDER_PREDICTOR_CORRECTOR_FULL_TANH,
    LADDER_PREDICTOR_CORRECTOR_FEEDBACK_TANH,
    LADDER_TRAPEZOIDAL_FEEDBACK_TANH
  };

  class Ladder{
  public:
    // constructor/destructor
    Ladder(DspFloatType newCutoff, DspFloatType newResonance, int newOversamplingFactor,
    LadderFilterMode newFilterMode, DspFloatType newSampleRate,
    LadderIntegrationMethod newIntegrationMethod, int newDecimatorOrder);
    Ladder();
    ~Ladder();

    // set filter parameters
    void SetFilterCutoff(DspFloatType newCutoff);
    void SetFilterResonance(DspFloatType newResonance);
    void SetFilterMode(LadderFilterMode newFilterMode);
    void SetFilterSampleRate(DspFloatType newSampleRate);
    void SetFilterIntegrationMethod(LadderIntegrationMethod method);
    void SetFilterOversamplingFactor(int newOversamplingFactor);
    void SetFilterDecimatorOrder(int decimatorOrder);
    
    // get filter parameters
    DspFloatType GetFilterCutoff();
    DspFloatType GetFilterResonance();
    LadderFilterMode GetFilterMode();  
    DspFloatType GetFilterSampleRate();
    LadderIntegrationMethod GetFilterIntegrationMethod();
    int GetFilterOversamplingFactor();  
    int GetFilterDecimatorOrder();
    
    // tick filter state
    void LadderFilter(DspFloatType input);

    // get filter responses
    DspFloatType GetFilterLowpass();
    DspFloatType GetFilterBandpass();
    DspFloatType GetFilterHighpass();

    // get filter output
    DspFloatType GetFilterOutput();

    // reset state
    void ResetFilterState();

  private:
    // set integration rate
    void SetFilterIntegrationRate();

    // filter parameters
    DspFloatType cutoffFrequency;
    DspFloatType Resonance;
    LadderFilterMode filterMode;
    DspFloatType sampleRate;
    DspFloatType dt;
    LadderIntegrationMethod integrationMethod;
    int oversamplingFactor;
    int decimatorOrder;
    
    // filter state
    DspFloatType p0, p1, p2, p3;
    DspFloatType ut_1;
    
    // filter output
    DspFloatType out;

    // IIR downsampling filter
    IIRLowpass *iir;
  };


  // constructor
  Ladder::Ladder(DspFloatType newCutoff, DspFloatType newResonance, int newOversamplingFactor,
          LadderFilterMode newFilterMode, DspFloatType newSampleRate,
          LadderIntegrationMethod newIntegrationMethod, int newDecimatorOrder){
    // initialize filter parameters
    cutoffFrequency = newCutoff;
    Resonance = newResonance;
    filterMode = newFilterMode;
    sampleRate = newSampleRate;
    oversamplingFactor = newOversamplingFactor;
    decimatorOrder = newDecimatorOrder;
    
    SetFilterIntegrationRate();

    // initialize filter state
    p0 = p1 = p2 = p3 = out = ut_1 = 0.0;
    
    integrationMethod = newIntegrationMethod;
    
    // instantiate downsampling filter
    iir = new IIRLowpass(sampleRate * oversamplingFactor, IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0, decimatorOrder);
  }

  // default constructor
  Ladder::Ladder(){
    // initialize filter parameters
    cutoffFrequency = 0.25;
    Resonance = 0.5;
    filterMode = LADDER_LOWPASS_MODE;
    sampleRate = 44100.0;
    oversamplingFactor = 2;
    decimatorOrder = IIR_DOWNSAMPLE_ORDER;
    
    SetFilterIntegrationRate();
    
    // initialize filter state
    p0 = p1 = p2 = p3 = out = ut_1 = 0.0;
    
    integrationMethod = LADDER_PREDICTOR_CORRECTOR_FULL_TANH;
    
    // instantiate downsampling filter
    iir = new IIRLowpass(sampleRate * oversamplingFactor, IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0, decimatorOrder);
  }

  // default destructor
  Ladder::~Ladder(){
    delete iir;
  }

  void Ladder::ResetFilterState(){
    // initialize filter parameters
    cutoffFrequency = 0.25;
    Resonance = 0.0;

    SetFilterIntegrationRate();
    
    // initialize filter state
    p0 = p1 = p2 = p3 = out = ut_1 = 0.0;
    
    // set oversampling
    iir->SetFilterSamplerate(sampleRate * oversamplingFactor);
    iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
    iir->SetFilterOrder(decimatorOrder);  
  }

  void Ladder::SetFilterCutoff(DspFloatType newCutoff){
    cutoffFrequency = newCutoff;

    SetFilterIntegrationRate();
  }

  void Ladder::SetFilterResonance(DspFloatType newResonance){
    Resonance = newResonance;
  }

  void Ladder::SetFilterMode(LadderFilterMode newFilterMode){
    filterMode = newFilterMode;
  }

  void Ladder::SetFilterSampleRate(DspFloatType newSampleRate){
    sampleRate = newSampleRate;
    iir->SetFilterSamplerate(sampleRate * (DspFloatType)(oversamplingFactor));
    iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
    iir->SetFilterOrder(decimatorOrder);

    SetFilterIntegrationRate();
  }

  void Ladder::SetFilterIntegrationMethod(LadderIntegrationMethod method){
    integrationMethod = method;
  }

  void Ladder::SetFilterOversamplingFactor(int newOversamplingFactor){
    oversamplingFactor = newOversamplingFactor;
    iir->SetFilterSamplerate(sampleRate * oversamplingFactor);
    iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
    iir->SetFilterOrder(decimatorOrder);

    SetFilterIntegrationRate();
  }

  void Ladder::SetFilterDecimatorOrder(int newDecimatorOrder){
    decimatorOrder = newDecimatorOrder;
    iir->SetFilterOrder(decimatorOrder);
  }

  void Ladder::SetFilterIntegrationRate(){
    // normalize cutoff freq to samplerate
    dt = 44100.0 / (sampleRate * oversamplingFactor) * cutoffFrequency;

    // clamp integration rate
    if(dt < 0.0){
      dt = 0.0;
    }
    else if(dt > 0.6){
      dt = 0.6;
    }
  }

  DspFloatType Ladder::GetFilterCutoff(){
    return cutoffFrequency;
  }

  DspFloatType Ladder::GetFilterResonance(){
    return Resonance;
  }

  int Ladder::GetFilterOversamplingFactor(){
    return oversamplingFactor;
  }

  int Ladder::GetFilterDecimatorOrder(){
    return decimatorOrder;
  }

  DspFloatType Ladder::GetFilterOutput(){
    return out;
  }

  LadderFilterMode Ladder::GetFilterMode(){
    return filterMode;
  }

  DspFloatType Ladder::GetFilterSampleRate(){
    return sampleRate;
  }

  LadderIntegrationMethod Ladder::GetFilterIntegrationMethod(){
    return integrationMethod;
  }

  void Ladder::LadderFilter(DspFloatType input){
    // noise term
    DspFloatType noise;

    // feedback amount
    DspFloatType fb = 8.0*Resonance;

    // update noise terms
    noise = static_cast <DspFloatType> (rand()) / static_cast <DspFloatType> (RAND_MAX);
    noise = 1.0e-6 * 2.0 * (noise - 0.5);

    input += noise;
    
    // integrate filter state
    // with oversampling
    for(int nn = 0; nn < oversamplingFactor; nn++){
      // switch integration method
      switch(integrationMethod){
      case LADDER_EULER_FULL_TANH:
        // semi-implicit euler integration
        // with full tanh stages
        {
    p0 = p0 + dt*(TanhPade32(input - fb*p3) - TanhPade32(p0));
    p1 = p1 + dt*(TanhPade32(p0) - TanhPade32(p1));
    p2 = p2 + dt*(TanhPade32(p1) - TanhPade32(p2));
    p3 = p3 + dt*(TanhPade32(p2) - TanhPade32(p3));
        }
        break;
        
      case LADDER_PREDICTOR_CORRECTOR_FULL_TANH:
        // predictor-corrector integration
        // with full tanh stages
        {
    DspFloatType p0_prime, p1_prime, p2_prime, p3_prime, p3t_1;

    // predictor
    p0_prime = p0 + dt*(TanhPade32(ut_1 - fb*p3) - TanhPade32(p0));
    p1_prime = p1 + dt*(TanhPade32(p0) - TanhPade32(p1));
    p2_prime = p2 + dt*(TanhPade32(p1) - TanhPade32(p2));
    p3_prime = p3 + dt*(TanhPade32(p2) - TanhPade32(p3));

    // corrector
    p3t_1 = p3;
    p3 = p3 + 0.5*dt*((TanhPade32(p2) - TanhPade32(p3)) + (TanhPade32(p2_prime) - TanhPade32(p3_prime)));
    p2 = p2 + 0.5*dt*((TanhPade32(p1) - TanhPade32(p2)) + (TanhPade32(p1_prime) - TanhPade32(p2_prime)));
    p1 = p1 + 0.5*dt*((TanhPade32(p0) - TanhPade32(p1)) + (TanhPade32(p0_prime) - TanhPade32(p1_prime)));
    p0 = p0 + 0.5*dt*((TanhPade32(ut_1 - fb*p3t_1) - TanhPade32(p0)) + (TanhPade32(input - fb*p3) - TanhPade32(p0_prime)));
        }
        break;
        
      case LADDER_PREDICTOR_CORRECTOR_FEEDBACK_TANH:
        // predictor-corrector integration
        // with feedback tanh stage only
        {
    DspFloatType p0_prime, p1_prime, p2_prime, p3_prime, p3t_1;

    // predictor
    p0_prime = p0 + dt*(TanhPade32(ut_1 - fb*p3) - p0);
    p1_prime = p1 + dt*(p0 - p1);
    p2_prime = p2 + dt*(p1 - p2);
    p3_prime = p3 + dt*(p2 - p3);

    // corrector
    p3t_1 = p3;
    p3 = p3 + 0.5*dt*((p2 - p3) + (p2_prime - p3_prime));
    p2 = p2 + 0.5*dt*((p1 - p2) + (p1_prime - p2_prime));
    p1 = p1 + 0.5*dt*((p0 - p1) + (p0_prime - p1_prime));
    p0 = p0 + 0.5*dt*((TanhPade32(ut_1 - fb*p3t_1) - p0) +
          (TanhPade32(input - fb*p3) - p0_prime));
        }
        break;
        
      case LADDER_TRAPEZOIDAL_FEEDBACK_TANH:
        // implicit trapezoidal integration
        // with feedback tanh stage only
        {
    DspFloatType x_k, x_k2, g, b, c, C_t, D_t, ut, ut_2;
    DspFloatType p0_prime, p1_prime, p2_prime, p3_prime;

    ut = TanhPade32(ut_1 - fb*p3);
        b = (0.5*dt)/(1.0 + 0.5*dt);
    c = (1.0 - 0.5*dt)/(1.0 + 0.5*dt);
    g = -fb*b*b*b*b;
    x_k = ut;
    D_t = c*p3 + (b + c*b)*p2 + (b*b+b*b*c)*p1 +
                  (b*b*b+b*b*b*c)*p0 + b*b*b*b*ut;
    C_t = TanhPade32(input - fb*D_t);

    // newton-raphson 
    for(int ii=0; ii < LADDER_MAX_NEWTON_STEPS; ii++) {
      DspFloatType tanh_g_xk, tanh_g_xk2;
      
      tanh_g_xk = TanhPade32(g*x_k);
      tanh_g_xk2 = g*(1.0 - TanhPade32(g*x_k)*TanhPade32(g*x_k));
      
      x_k2 = x_k - (x_k + x_k*tanh_g_xk*C_t - tanh_g_xk - C_t) /
                    (1.0 + C_t*(tanh_g_xk + x_k*tanh_g_xk2) - tanh_g_xk2);
      
  #ifdef LADDER_NEWTON_BREAKING_LIMIT
      // breaking limit
      if(abs(x_k2 - x_k) < 1.0e-9) {
        x_k = x_k2;
        break;
      }
  #endif	  
      x_k = x_k2;
    }
    
    ut_2 = x_k;

    p0_prime = p0;
    p1_prime = p1;
    p2_prime = p2;
    p3_prime = p3;

    p0 = c*p0_prime + b*(ut + ut_2);
    p1 = c*p1_prime + b*(p0_prime + p0);
    p2 = c*p2_prime + b*(p1_prime + p1);
    p3 = c*p3_prime + b*(p2_prime + p2);
        }
        break;
        
      default:
        break;
      }

      // input at t-1
      ut_1 = input;

      //switch filter mode
      switch(filterMode){
      case LADDER_LOWPASS_MODE:
        out = p3;
        break;
      case LADDER_BANDPASS_MODE:
        out = p1 - p3;
        break;
      case LADDER_HIGHPASS_MODE:
        out = TanhPade32(input - p0 - fb*p3);
        break;
      default:
        out = 0.0;
      }

      // downsampling filter
      if(oversamplingFactor > 1){
        out = iir->IIRfilter(out);
      }
    }
  }

  DspFloatType Ladder::GetFilterLowpass(){
    return p3;
  }

  DspFloatType Ladder::GetFilterBandpass(){
    return 0.0;
  }

  DspFloatType Ladder::GetFilterHighpass(){
    return 0.0;
  }

  /////////////////////////////////////////////////////
  // Sallen-Key (Korg35)
  /////////////////////////////////////////////////////

  // filter modes
  enum SKFilterMode {
    SK_LOWPASS_MODE,
    SK_BANDPASS_MODE,
    SK_HIGHPASS_MODE
  };

  // integration methods
  enum SKIntegrationMethod {
    SK_SEMI_IMPLICIT_EULER,
    SK_PREDICTOR_CORRECTOR,
    SK_TRAPEZOIDAL
  };

  class SKFilter{
  public:
    // constructor/destructor
    SKFilter(DspFloatType newCutoff, DspFloatType newResonance, int newOversamplingFactor,
      SKFilterMode newFilterMode, DspFloatType newSampleRate,
      SKIntegrationMethod newIntegrationMethod, int newDecimatorOrder);
    SKFilter();
    ~SKFilter();

    // set filter parameters
    void SetFilterCutoff(DspFloatType newCutoff);
    void SetFilterResonance(DspFloatType newResonance);
    void SetFilterMode(SKFilterMode newFilterMode);
    void SetFilterSampleRate(DspFloatType newSampleRate);
    void SetFilterIntegrationMethod(SKIntegrationMethod method);
    void SetFilterOversamplingFactor(int newOversamplingFactor);
    void SetFilterDecimatorOrder(int decimatorOrder);
    
    // get filter parameters
    DspFloatType GetFilterCutoff();
    DspFloatType GetFilterResonance();
    SKFilterMode GetFilterMode();  
    DspFloatType GetFilterSampleRate();
    SKIntegrationMethod GetFilterIntegrationMethod();
    int GetFilterOversamplingFactor();  
    int GetFilterDecimatorOrder();
    
    // tick filter state
    void filter(DspFloatType input);

    // set filter inputs
    void SetFilterLowpassInput(DspFloatType input);
    void SetFilterBandpassInput(DspFloatType input);
    void SetFilterHighpassInput(DspFloatType input);

    // get filter output
    DspFloatType GetFilterOutput();

    // reset state
    void ResetFilterState();
    
  private:
    // set integration rate
    void SetFilterIntegrationRate();

    // filter parameters
    DspFloatType cutoffFrequency;
    DspFloatType Resonance;
    SKFilterMode filterMode;
    DspFloatType sampleRate;
    DspFloatType dt;
    SKIntegrationMethod integrationMethod;
    int oversamplingFactor;
    int decimatorOrder;
    
    // filter state
    DspFloatType p0;
    DspFloatType p1;

    // filter input
    DspFloatType input_lp;
    DspFloatType input_bp;
    DspFloatType input_hp;
    DspFloatType input_lp_t1;
    DspFloatType input_bp_t1;
    DspFloatType input_hp_t1;
    
    // filter output
    DspFloatType out;

    // IIR downsampling filter
    IIRLowpass *iir;
  };


  /*
  *  (C) 2021 Janne Heikkarainen <janne808@radiofreerobotron.net>
  *
  *  All rights reserved.
  *
  *  This file is part of Kocmoc VCV Rack plugin.
  *
  *  Kocmoc VCV Rack plugin is free software: you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation, either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  Kocmoc VCV Rack plugin is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with Kocmoc VCV Rack plugin.  If not, see <http://www.gnu.org/licenses/>.
  */


  // constructor
  SKFilter::SKFilter(DspFloatType newCutoff, DspFloatType newResonance, int newOversamplingFactor,
        SKFilterMode newFilterMode, DspFloatType newSampleRate,
        SKIntegrationMethod newIntegrationMethod, int newDecimatorOrder){
    // initialize filter parameters
    cutoffFrequency = newCutoff;
    Resonance = newResonance;
    filterMode = newFilterMode;
    sampleRate = newSampleRate;
    oversamplingFactor = newOversamplingFactor;
    decimatorOrder = newDecimatorOrder;

    SetFilterIntegrationRate();

    // initialize filter state
    p0 = p1 = out = 0.0;

    // initialize filter inputs
    input_lp = input_bp = input_hp = 0.0;
    input_lp_t1 = input_bp_t1 = input_hp_t1 = 0.0;
    
    integrationMethod = newIntegrationMethod;
    
    // instantiate downsampling filter
    iir = new IIRLowpass(sampleRate * oversamplingFactor,
            IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0,
            decimatorOrder);
  }

  // default constructor
  SKFilter::SKFilter(){
    // initialize filter parameters
    cutoffFrequency = 0.25;
    Resonance = 0.5;
    filterMode = SK_LOWPASS_MODE;
    sampleRate = 44100.0;
    oversamplingFactor = 2;
    decimatorOrder = IIR_DOWNSAMPLE_ORDER;

    SetFilterIntegrationRate();
    
    // initialize filter state
    p0 = p1 = out = 0.0;

    // initialize filter inputs
    input_lp = input_bp = input_hp = 0.0;
    input_lp_t1 = input_bp_t1 = input_hp_t1 = 0.0;
    
    integrationMethod = SK_TRAPEZOIDAL;
    
    // instantiate downsampling filter
    iir = new IIRLowpass(sampleRate * oversamplingFactor,
            IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0,
            decimatorOrder);
  }

  // default destructor
  SKFilter::~SKFilter(){
    delete iir;
  }

  void SKFilter::ResetFilterState(){
    // initialize filter parameters
    cutoffFrequency = 0.25;
    Resonance = 0.5;

    SetFilterIntegrationRate();
    
    // initialize filter state
    p0 = p1 = out = 0.0;

    // initialize filter inputs
    input_lp = input_bp = input_hp = 0.0;
    input_lp_t1 = input_bp_t1 = input_hp_t1 = 0.0;
    
    // set oversampling
    iir->SetFilterSamplerate(sampleRate * oversamplingFactor);
    iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
    iir->SetFilterOrder(decimatorOrder);
  }

  void SKFilter::SetFilterCutoff(DspFloatType newCutoff){
    cutoffFrequency = newCutoff;

    SetFilterIntegrationRate();
  }

  void SKFilter::SetFilterResonance(DspFloatType newResonance){
    Resonance = newResonance;
  }

  void SKFilter::SetFilterMode(SKFilterMode newFilterMode){
    filterMode = newFilterMode;
  }

  void SKFilter::SetFilterSampleRate(DspFloatType newSampleRate){
    sampleRate = newSampleRate;
    iir->SetFilterSamplerate(sampleRate * (DspFloatType)(oversamplingFactor));
    iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);

    SetFilterIntegrationRate();
  }

  void SKFilter::SetFilterIntegrationMethod(SKIntegrationMethod method){
    integrationMethod = method;
  }

  void SKFilter::SetFilterOversamplingFactor(int newOversamplingFactor){
    oversamplingFactor = newOversamplingFactor;
    iir->SetFilterSamplerate(sampleRate * oversamplingFactor);
    iir->SetFilterCutoff(IIR_DOWNSAMPLING_BANDWIDTH*sampleRate/2.0);
    iir->SetFilterOrder(decimatorOrder);

    SetFilterIntegrationRate();
  }

  void SKFilter::SetFilterDecimatorOrder(int newDecimatorOrder){
    decimatorOrder = newDecimatorOrder;
    iir->SetFilterOrder(decimatorOrder);
  }

  void SKFilter::SetFilterIntegrationRate(){
    // normalize cutoff freq to samplerate
    dt = 44100.0 / (sampleRate * oversamplingFactor) * cutoffFrequency;

    // clamp integration rate
    if(dt < 0.0){
      dt = 0.0;
    }
    else if(dt > 0.55){
      dt = 0.55;
    }
  }

  DspFloatType SKFilter::GetFilterCutoff(){
    return cutoffFrequency;
  }

  DspFloatType SKFilter::GetFilterResonance(){
    return Resonance;
  }

  int SKFilter::GetFilterOversamplingFactor(){
    return oversamplingFactor;
  }

  int SKFilter::GetFilterDecimatorOrder(){
    return decimatorOrder;
  }

  DspFloatType SKFilter::GetFilterOutput(){
    return out;
  }

  SKFilterMode SKFilter::GetFilterMode(){
    return filterMode;
  }

  DspFloatType SKFilter::GetFilterSampleRate(){
    return sampleRate;
  }

  SKIntegrationMethod SKFilter::GetFilterIntegrationMethod(){
    return integrationMethod;
  }

  void SKFilter::filter(DspFloatType input){
    // noise term
    DspFloatType noise;

    // feedback amount variables
    DspFloatType res=4.0*Resonance;
    DspFloatType fb=0.0;

    // update noise terms
    noise = static_cast <DspFloatType> (rand()) / static_cast <DspFloatType> (RAND_MAX);
    noise = 1.0e-6 * 2.0 * (noise - 0.5);

    input += noise;

    // set filter mode
    switch(filterMode){
    case SK_LOWPASS_MODE:
      input_lp = input;
      input_bp = 0.0;
      input_hp = 0.0;
      break;
    case SK_BANDPASS_MODE:
      input_lp = 0.0;
      input_bp = input;
      input_hp = 0.0;
      break;
    case SK_HIGHPASS_MODE:
      input_lp = 0.0;
      input_bp = 0.0;
      input_hp = input;
      break;
    default:
      input_lp = 0.0;
      input_bp = 0.0;
      input_hp = 0.0;
    }
      
    // integrate filter state
    // with oversampling
    for(int nn = 0; nn < oversamplingFactor; nn++){
      // switch integration method
      switch(integrationMethod){
      case SK_SEMI_IMPLICIT_EULER:
        // semi-implicit euler integration
        {
    fb = input_bp + res*p1;
    p0 += dt*(input_lp - p0 - fb);
          p1 += dt*(p0 + fb - p1 - 1.0/4.0*SinhPade34(p0*4.0));
          out = p1;
        }
        break;
      case SK_PREDICTOR_CORRECTOR:
        // predictor-corrector integration
        {
    DspFloatType p0_prime, p1_prime, fb_prime;
      
    fb = input_bp_t1 + res*p1;
    p0_prime = p0 + dt*(input_lp_t1 - p0 - fb);
          p1_prime = p1 + dt*(p0 + fb - p1 - 1.0/4.0*SinhPade34(p1*4.0));	
    fb_prime = input_bp + res*p1_prime;
    
          p1 += 0.5*dt*((p0 + fb - p1 - 1.0/4.0*SinhPade34(p1*4.0)) +
            (p0_prime + fb_prime - p1_prime - 1.0/4.0*SinhPade34(p1*4.0)));
    p0 += 0.5*dt*((input_lp_t1 - p0 - fb) +
            (input_lp - p0_prime - fb_prime));

    out = p1;
        }
        break;
      case SK_TRAPEZOIDAL:
        // trapezoidal integration
        {
    DspFloatType x_k, x_k2;
    DspFloatType fb_t = input_bp_t1 + res*p1;
    DspFloatType alpha = dt/2.0;
    DspFloatType A = p0 + fb_t - p1 - 1.0/4.0*SinhPade54(4.0*p1) +
              p0/(1.0 + alpha) + alpha/(1 + alpha)*(input_lp_t1 - p0 - fb_t + input_lp);
    DspFloatType c = 1.0 - (alpha - alpha*alpha/(1.0 + alpha))*res + alpha;
    DspFloatType D_n = p1 + alpha*A + (alpha - alpha*alpha/(1.0 + alpha))*input_bp;

    x_k = p1;
    
    // newton-raphson
    for(int ii=0; ii < SKF_MAX_NEWTON_STEPS; ii++) {
      x_k2 = x_k - (c*x_k + alpha*1.0/4.0*SinhPade54(4.0*x_k) - D_n)/(c + alpha*CoshPade54(4.0*x_k));
      
  #ifdef SKF_NEWTON_BREAKING_LIMIT
      // breaking limit
      if(abs(x_k2 - x_k) < 1.0e-9) {
        x_k = x_k2;
        break;
      }
  #endif	  
      x_k = x_k2;
    }
    
    p1 = x_k;
    fb = input_bp + res*p1;
    p0 = p0/(1.0 + alpha) + alpha/(1.0 + alpha)*(input_lp_t1 - p0 - fb_t + input_lp - fb);
    out = p1;
        }
        break;
      default:
        break;
      }

      // downsampling filter
      if(oversamplingFactor > 1){
        out = iir->IIRfilter(out);
      }
    }
    
    // set input at t-1
    input_lp_t1 = input_lp;    
    input_bp_t1 = input_bp;    
    input_hp_t1 = input_hp;    
  }

  void SKFilter::SetFilterLowpassInput(DspFloatType input){
    input_lp = input;
  }

  void SKFilter::SetFilterBandpassInput(DspFloatType input){
    input_bp = input;
  }

  void SKFilter::SetFilterHighpassInput(DspFloatType input){
    input_hp = input;
  }


  /////////////////////////////////////////////////////
  // State Variable Filter
  /////////////////////////////////////////////////////
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
    SVFilter(DspFloatType newCutoff, DspFloatType newResonance, int newOversamplingFactor,
      SVFFilterMode newFilterMode, DspFloatType newSampleRate,
      SVFIntegrationMethod newIntegrationMethod, int newDecimatorOrder);
    SVFilter();
    ~SVFilter();

    // set filter parameters
    void SetFilterCutoff(DspFloatType newCutoff);
    void SetFilterResonance(DspFloatType newResonance);
    void SetFilterMode(SVFFilterMode newFilterMode);
    void SetFilterSampleRate(DspFloatType newSampleRate);
    void SetFilterIntegrationMethod(SVFIntegrationMethod method);
    void SetFilterOversamplingFactor(int newOversamplingFactor);
    void SetFilterDecimatorOrder(int decimatorOrder);
      
    // get filter parameters
    DspFloatType GetFilterCutoff();
    DspFloatType GetFilterResonance();
    SVFFilterMode GetFilterMode();  
    DspFloatType GetFilterSampleRate();
    SVFIntegrationMethod GetFilterIntegrationMethod();
    int GetFilterOversamplingFactor();  
    int GetFilterDecimatorOrder();
    
    // tick filter state
    void filter(DspFloatType input);

    // get filter responses
    DspFloatType GetFilterLowpass();
    DspFloatType GetFilterBandpass();
    DspFloatType GetFilterHighpass();

    // get filter output
    DspFloatType GetFilterOutput();

    // reset state
    void ResetFilterState();
    
  private:
    // set integration rate
    void SetFilterIntegrationRate();

    // pade approximant functions for hyperbolic functions
    // filter parameters
    DspFloatType cutoffFrequency;
    DspFloatType Resonance;
    SVFFilterMode filterMode;
    SVFIntegrationMethod integrationMethod;
    DspFloatType dt;
    DspFloatType sampleRate;
    int oversamplingFactor;
    int decimatorOrder;
    
    // filter state
    DspFloatType lp;
    DspFloatType bp;
    DspFloatType hp;
    DspFloatType u_t1;
    
    // filter output
    DspFloatType out;

    // IIR downsampling filter
    IIRLowpass *iir;
  };

  // constructor
  SVFilter::SVFilter(DspFloatType newCutoff, DspFloatType newResonance, int newOversamplingFactor,
        SVFFilterMode newFilterMode, DspFloatType newSampleRate,
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

  void SVFilter::SetFilterCutoff(DspFloatType newCutoff){
    cutoffFrequency = newCutoff;

    SetFilterIntegrationRate();
  }

  void SVFilter::SetFilterResonance(DspFloatType newResonance){
    Resonance = newResonance;
  }

  void SVFilter::SetFilterMode(SVFFilterMode newFilterMode){
    filterMode = newFilterMode;
  }

  void SVFilter::SetFilterSampleRate(DspFloatType newSampleRate){
    sampleRate = newSampleRate;
    iir->SetFilterSamplerate(sampleRate * (DspFloatType)(oversamplingFactor));
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
    dt = 44100.0 / (sampleRate * (DspFloatType)(oversamplingFactor)) * cutoffFrequency;

    // clamp integration rate
    if(dt < 0.0){
      dt=0.0;
    }
  }

  DspFloatType SVFilter::GetFilterCutoff(){
    return cutoffFrequency;
  }

  DspFloatType SVFilter::GetFilterResonance(){
    return Resonance;
  }

  DspFloatType SVFilter::GetFilterOutput(){
    return out;
  }

  SVFFilterMode SVFilter::GetFilterMode(){
    return filterMode;
  }

  DspFloatType SVFilter::GetFilterSampleRate(){
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

  void SVFilter::filter(DspFloatType input){
    // noise term
    DspFloatType noise;

    // feedback amount variables
    DspFloatType fb = 1.0 - (3.5*Resonance);

    // integration rate
    DspFloatType dt2 = dt;
    
    // update noise terms
    noise = static_cast <DspFloatType> (rand()) / static_cast <DspFloatType> (RAND_MAX);
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
    DspFloatType beta = 1.0 - (0.0075/oversamplingFactor);

          hp = input - lp - fb*bp - SinhPade54(bp);
    bp += dt2*hp;
    bp *= beta;
    lp += dt2*bp;
        }
        break;
      case SVF_TRAPEZOIDAL:
        // trapezoidal integration
        {
    DspFloatType alpha = dt2/2.0;
    DspFloatType beta = 1.0 - (0.0075/oversamplingFactor);
    DspFloatType alpha2 = dt2*dt2/4.0 + fb*alpha;
    DspFloatType D_t = (1.0 - dt2*dt2/4.0)*bp +
                  alpha*(u_t1 + input - 2.0*lp - fb*bp - SinhPade54(bp));
    DspFloatType x_k, x_k2;

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
    DspFloatType alpha = dt2/2.0;
    DspFloatType beta = 1.0 - (0.0075/oversamplingFactor);
    DspFloatType alpha2 = dt2*dt2/4.0 + fb*alpha;
    DspFloatType D_t = (1.0 - dt2*dt2/4.0)*bp +
                  alpha*(u_t1 + input - 2.0*lp - fb*bp - sinh(bp));
    DspFloatType y_k, y_k2;

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

  DspFloatType SVFilter::GetFilterLowpass(){
    return lp;
  }

  DspFloatType SVFilter::GetFilterBandpass(){
    return bp;
  }

  DspFloatType SVFilter::GetFilterHighpass(){
    return hp;
  }
}
