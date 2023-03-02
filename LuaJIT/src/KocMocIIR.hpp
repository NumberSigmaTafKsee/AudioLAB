class IIRLowpass{
public:
  // constructor/destructor
  IIRLowpass(double newSamplerate, double newCutoff, int newOrder);
  IIRLowpass();
  ~IIRLowpass();

  // set filter parameters
  void SetFilterOrder(int newOrder);
  void SetFilterSamplerate(double newSamplerate);
  void SetFilterCutoff(double newCutoff);

  // initialize biquad cascade delayline
  void InitializeBiquadCascade();
  
  // IIR filter signal 
  double IIRfilter(double input);

  // get filter coefficients
  double* GetFilterCoeffA1();
  double* GetFilterCoeffA2();
  double* GetFilterCoeffK();
  
private:
  // compute biquad cascade coefficients
  void ComputeCoefficients();

  // filter design variables
  double samplerate;
  double cutoff;
  int order;
  
  // dsp variables
  double *a1;
  double *a2;
  double *K;
  double *pa_real;
  double *pa_imag;
  double *p_real;
  double *p_imag;
  
  // cascaded biquad buffers
  double *z;
};

#define IIR_MAX_ORDER 32

// constructor
IIRLowpass::IIRLowpass(double newSamplerate, double newCutoff, int newOrder)
{
  // initialize filter design parameters
  samplerate = newSamplerate;
  cutoff = newCutoff;
  order = newOrder;

  // allocate dsp vectors
  a1 = new double[IIR_MAX_ORDER/2];
  a2 = new double[IIR_MAX_ORDER/2];
  K = new double[IIR_MAX_ORDER/2];
  pa_real = new double[IIR_MAX_ORDER/2];
  pa_imag = new double[IIR_MAX_ORDER/2];
  p_real = new double[IIR_MAX_ORDER/2];
  p_imag = new double[IIR_MAX_ORDER/2];

  // allocate cascaded biquad buffer
  z = new double[IIR_MAX_ORDER];
  
  // initialize cascade delayline
  InitializeBiquadCascade();
  
  // compute impulse response
  ComputeCoefficients();
}

// default constructor
IIRLowpass::IIRLowpass()
{
  // set default design parameters
  samplerate=(double)(44100.0);
  cutoff=(double)(440.0);
  order=IIR_MAX_ORDER;
  
  // allocate dsp vectors
  a1 = new double[IIR_MAX_ORDER/2];
  a2 = new double[IIR_MAX_ORDER/2];
  K = new double[IIR_MAX_ORDER/2];
  pa_real = new double[IIR_MAX_ORDER/2];
  pa_imag = new double[IIR_MAX_ORDER/2];
  p_real = new double[IIR_MAX_ORDER/2];
  p_imag = new double[IIR_MAX_ORDER/2];

  // allocate cascaded biquad buffer
  z = new double[IIR_MAX_ORDER];
  
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

void IIRLowpass::SetFilterSamplerate(double newSamplerate){
  samplerate = newSamplerate;

  // initialize cascade delayline
  InitializeBiquadCascade();
  
  // compute new cascade coefficients
  ComputeCoefficients();
}

void IIRLowpass::SetFilterCutoff(double newCutoff){
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

double IIRLowpass::IIRfilter(double input){
  double out=input;
  double in;

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

double* IIRLowpass::GetFilterCoeffA1(){
  return a1;
}

double* IIRLowpass::GetFilterCoeffA2(){
  return a2;
}

double* IIRLowpass::GetFilterCoeffK(){
  return K;
}

void IIRLowpass::ComputeCoefficients(){
  // place butterworth style analog filter poles
  double theta;

  for(int ii = 0; ii<order/2; ii++) {
    int k = order/2 - ii;
    theta = (2.0*(double)(k) - 1.0)*M_PI/(2.0*(double)(order));
    
    pa_real[ii] = -1.0*sin(theta);
    pa_imag[ii] = cos(theta);
  }

  // prewarp and scale poles
  double Fc = samplerate/M_PI*tan(M_PI*cutoff/samplerate);  
  for(int ii = 0; ii<order/2; ii++) {
    pa_real[ii] *= 2.0*M_PI*Fc; 
    pa_imag[ii] *= 2.0*M_PI*Fc; 
  }

  // bilinear transform to z-plane
  for(int ii = 0; ii<order/2; ii++) {
    // complex division
    double u = (2.0*samplerate+pa_real[ii])/(2.0*samplerate); 
    double v = pa_imag[ii]/(2.0*samplerate); 
    double x = (2.0*samplerate-pa_real[ii])/(2.0*samplerate); 
    double y = -1.0*pa_imag[ii]/(2.0*samplerate);
    
    double c = 1.0/(x*x + y*y);
    
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
