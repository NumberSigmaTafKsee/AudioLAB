#pragma once

namespace Loudia
{
    class OnsetComplex {
    protected:
        // Internal parameters
        int _frameLength;
        int _fftLength;
        bool _zeroPhase;
        Window::WindowType _windowType;
        
        // Internal variables
        Window _window;
        FFT _fft;
        SpectralODFComplex _odf;

        MatrixXR _windowed;
        MatrixXC _ffted;
        
    public:
        OnsetComplex(int frameLength, int fftLength, Window::WindowType windowType = Window::RECTANGULAR, bool zeroPhase = true);

        ~OnsetComplex();

        void setup();

        void process(const MatrixXR& samples, MatrixXR* odfValue);
        void reset();
    };

OnsetComplex::OnsetComplex(int frameLength, int fftLength, Window::WindowType windowType, bool zeroPhase) : 
  _frameLength( frameLength ),
  _fftLength( fftLength ),
  _zeroPhase( zeroPhase ),
  _windowType( windowType ),
  _window( frameLength, windowType ), 
  _fft( fftLength, zeroPhase ),
  _odf( fftLength )
{
  
  LOUDIA_DEBUG("OnsetComplex: Constructor frameLength: " << frameLength << 
        ", fftLength: " << fftLength);
   
  setup();
}

OnsetComplex::~OnsetComplex() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void OnsetComplex::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("OnsetComplex: Setting up...");

  _window.setup();
  _fft.setup();
  _odf.setup();

  reset();

  LOUDIA_DEBUG("OnsetComplex: Finished set up...");
}


void OnsetComplex::process(const MatrixXR& samples, MatrixXR* odfValue) {
  LOUDIA_DEBUG("OnsetComplex: Processing windowed");
 
  _window.process(samples, &_windowed);

  _fft.process(_windowed, &_ffted);
  
  _odf.process(_ffted, odfValue);
    
  LOUDIA_DEBUG("OnsetComplex: Finished Processing");
}

void OnsetComplex::reset() {
  // Initial values
  _window.reset();
  _fft.reset();
  _odf.reset();
}
}