#pragma once

namespace Loudia
{
    /**
    * @class SpectralNoiseSuppression
    *
    * @brief Algorithm to remove the non-harmonic part of the spectrums magnitudes represented as vectors of Real values.
    *
    * This class represents an object to perform spectral noise suppresion on vectors
    * of Real values.  Which is a useful technique to keep only the peaks of a spectrum magnitude
    * in harmonic sounds.
    *
    * This implementation consists in estimating the spectral noise by performing a
    * moving average on the power warped spectrum magnitude using varying bandwidths.
    * The spectral noise is then removed from the original spectrum, clipping the result to zero to avoid
    * negative values.
    *
    * The sampleRate and FFT size of the input spectrum are specified using setSampleRate() and
    * setFftSize().
    *
    * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
    * setHighFrequency().
    *
    * @author Ricard Marxer
    *
    * @sa MelBands, Bands, PeakDetection
    */
    class SpectralNoiseSuppression {
    protected:
        int _fftSize;
        Real _sampleRate;

        Real _lowFrequency;
        Real _highFrequency;

        int _k0;
        int _k1;

        MatrixXR _g;

        Bands _bands;

    public:
        /**
            Constructs a spectral noise suppression object with the specified @a lowFrequency, @a highFrequency,
            @a sampleRate and @a fftSize settings.
            @param lowFrequency low frequency used for the magnitude warping function,
            must be greater than zero 0 and lower than half the sampleRate.
            @param highFrequency high frequency used for the magnitude warping function,
            must be greater than zero 0 and lower than half the sampleRate.
            @param sampleRate sampleRate frequency of the input signal.
            @param fftSize size of the FFT.
        */
        SpectralNoiseSuppression( int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 6000.0, Real sampleRate = 44100.0 );

        /**
            Destroys the algorithm and frees its resources.
        */
        ~SpectralNoiseSuppression();

        void setup();
        void reset();

        /**
            Performs the estimation and suppression of the noise on each of the rows of @a spectrums.
            Puts the resulting noise spectrums and noise-suppressed spectrums in the rows of @a whitened.
            @param spectrums matrix of Real values representing one spectrum magnitude per row.
            The number of columns of @a spectrum must be equal to the fftSize / 2 + 1 where
            fftSize is specified using setFftSize().
            @param noises pointer to a matrix of Real values representing one noise spectrum per row.
            The matrix should have the same number of rows and columns as @a spectrums.
            @param suppressed pointer to a matrix of Real values representing one noise-suppressed spectrum per row.
            The matrix should have the same number of rows and columns as @a spectrums.
            Note that if the output matrices are not of the required sizes they will be resized,
            reallocating a new memory space if necessary.
        */
        void process( const MatrixXR& spectrums, MatrixXR* noises, MatrixXR* suppressed );

        /**
            Return the low frequency of the spectral whitening.
            The default is 50.0.
            @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
        */
        Real lowFrequency() const;

        /**
            Specifies the low @a frequency of the spectral whitening.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            @sa lowFrequency, highFrequency, setHighFrequency
        */
        void setLowFrequency( Real frequency, bool callSetup = true );

        /**
            Return the high frequency of the spectral whitening.
            The default is 6000.0.
            @sa lowFrequency, setLowFrequency, setHighFrequency
        */
        Real highFrequency() const;

        /**
            Specifies the high @a frequency of the spectral whitening.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            @sa lowFrequency, highFrequency, setLowFrequency
        */
        void setHighFrequency( Real frequency, bool callSetup = true );

        /**
            Return the sampleRate frequency of the input signal.
            The default is 44100.0.
            @sa setSampleRate
        */
        Real sampleRate() const;

        /**
            Specifies the sampleRate @a frequency of the input signal.
            @sa sampleRate
        */
        void setSampleRate( Real frequency, bool callSetup = true );

        /**
            Returns the size of the FFT that has been performed for the input.
            The default is 1024.
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the FFT that has been performed for the input.
            The given @a size must be higher than 0.
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );
    };

    
    class SpectralReassignment{
    protected:
        int _frameSize;
        int _fftSize;
        Real _sampleRate;
        Window::WindowType _windowType;
        
        Window _windowAlgo;
        Window _windowIntegAlgo;
        Window _windowDerivAlgo;

        FFTComplex _fftAlgo;

        MatrixXC _window;
        MatrixXC _windowInteg;
        MatrixXC _windowDeriv;

        MatrixXR _fftAbs2;
        MatrixXC _fftInteg;
        MatrixXC _fftDeriv;

        MatrixXR _time;
        MatrixXR _freq;
        
    public: 
        SpectralReassignment(int frameSize, int fftSize, Real sampleRate, Window::WindowType windowType = Window::RECTANGULAR);
        ~SpectralReassignment();
        
        void process(const MatrixXR& frames,
                    MatrixXC* fft, MatrixXR* reassignTime, MatrixXR* reassignFreq);
        
        void setup();
        void reset();

        int frameSize() const;
        int fftSize() const;

        Window::WindowType windowType() const;
    };    

    /**
    * @class SpectralWhitening
    *
    * @brief Algorithm to whiten the magnitude of spectrums represented as vectors of Real values.
    *
    * This class represents an object to perform spectral whitening on vectors 
    * of Real values.  Which is a useful technique to make the peaks of a spectrum magnitude 
    * stand out in harmonic sounds.
    *
    * This implementation consists in calculating the Mel bands and create a linear interpolation
    * between these to be used as a wheighting parameter of the spectrum's compression. 
    *
    * The sampleRate and FFT size of the input spectrum are specified using setSampleRate() and
    * setFftSize().
    *
    * The frequency limits of the Mel scale mapping are specified using setLowFrequency() and
    * setHighFrequency().
    *
    * The number of Mel bands is specified using setBandCount().
    *
    * The compression factor of the whitening process is specified by setCompressionFactor().
    *
    * @author Ricard Marxer
    *
    * @sa MelBands, Bands, PeakDetection
    */
    class SpectralWhitening {
    protected:
        int _fftSize;
        int _halfSize;
        Real _lowFrequency;
        Real _highFrequency;

        Real _sampleRate;
        Real _compressionFactor;
        int _bandCount;

        MelBands::ScaleType _scaleType;

        MatrixXR _centers;

        MatrixXR _bandEnergy;
        MatrixXR _compressionWeights;

        MelBands _bands;

    public:
        /**
            Constructs a spectral whitening object with the specified @a lowFrequency, @a highFrequency, 
            @a bandCount, @a sampleRate, @a fftSize, @a compressionFactor and @a scaleType settings.
            
            @param lowFrequency frequency of the lowest Mel band,
            must be greater than zero 0 and lower than half the sampleRate.
            
            @param highFrequency frequency of the highest Mel band,
            must be greater than zero 0 and lower than half the sampleRate.
        
            @param bandCount number of Mel bands.
            
            @param sampleRate sampleRate frequency of the input signal.
            @param fftSize size of the FFT.
            
            @param compressionFactor factor of the compression process in the whitening.
            
            @param scaleType scale used for the frequency warping.
        */  
        SpectralWhitening(int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 6000.0, Real sampleRate = 44100.0, Real compressionFactor = 0.33, int bandCount = 40, MelBands::ScaleType scaleType = MelBands::GREENWOOD);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~SpectralWhitening();

        void setup();
        void reset();

        /**
            Performs a whitening on each of the rows of @a spectrums.
            Puts the resulting whitened spectrum in the rows of @a whitened.
            
            @param spectrums matrix of Real values representing one spectrum magnitude per row.
            The number of columns of @a spectrum must be equal to the fftSize / 2 + 1 where 
            fftSize is specified using setFftSize().
            
            @param whitened pointer to a matrix of Real values representing one whitened spectrum per row.
            The matrix should have the same number of rows and columns as @a spectrums.
            
            Note that if the output matrices are not of the required sizes they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& spectrums, MatrixXR* whitened);

        /**
            Returns the number of bands to be performed.
            The default is 40.
            
            @sa setBandCount()
        */
        int bandCount() const;

        /**
            Specifies the @a count of bands to be performed.
                
            @sa bandCount()
        */
        void setBandCount( int count, bool callSetup = true );

        /**
            Return the low frequency of the spectral whitening.
            The default is 50.0.
            @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
        */  
        Real lowFrequency() const;  

        /**
            Specifies the low @a frequency of the spectral whitening.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            
            @sa lowFrequency, highFrequency, setHighFrequency
        */
        void setLowFrequency( Real frequency, bool callSetup = true );

        /**
            Return the high frequency of the spectral whitening.
            The default is 6000.0.
            @sa lowFrequency, setLowFrequency, setHighFrequency
        */  
        Real highFrequency() const;  

        /**
            Specifies the high @a frequency of the spectral whitening.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            @sa lowFrequency, highFrequency, setLowFrequency
        */
        void setHighFrequency( Real frequency, bool callSetup = true );

        /**
            Return the sampleRate frequency of the input signal.
            The default is 44100.0.
            @sa setSampleRate
        */  
        Real sampleRate() const;  

        /**
            Specifies the sampleRate @a frequency of the input signal.
            
            @sa sampleRate
        */
        void setSampleRate( Real frequency, bool callSetup = true );

        /**
            Returns the size of the FFT that has been performed for the input.
            The default is 1024.
            
            @sa setFftSize()
        */
        int fftSize() const;

        /**
            Specifies the @a size of the FFT that has been performed for the input.
            The given @a size must be higher than 0.
            
            @sa fftSize()
        */
        void setFftSize( int size, bool callSetup = true );

        /**
            Return the type of the frequency warping scale.
            
            By default it is GREENWOOD.
            
            @sa setScaleType()
        */
        MelBands::ScaleType scaleType() const;
        
        /**
            Specify the type of the frequency warping scale.
            
            @sa scaleType()
        */
        void setScaleType( MelBands::ScaleType type, bool callSetup = true );

        /**
            Return the compression factor of the whitening.
            The default is 0.33.
            @sa setCompressionFactor
        */  
        Real compressionFactor() const;  

        /**
            Specifies the compression @a factor of the whitening.
            
            @sa compressionFactor
        */
        void setCompressionFactor( Real factor, bool callSetup = true );  
        
    };    


SpectralReassignment::SpectralReassignment(int frameSize, int fftSize, Real sampleRate, Window::WindowType windowType) : 
  _frameSize( frameSize ),
  _fftSize( fftSize ),
  _sampleRate( sampleRate ),
  _windowType( windowType ),
  _windowAlgo( frameSize, windowType ), 
  _windowIntegAlgo( frameSize, Window::CUSTOM ), 
  _windowDerivAlgo( frameSize, Window::CUSTOM ), 
  _fftAlgo( frameSize, fftSize, true )
{
  
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Constructor frameSize: " << frameSize << \
        ", fftSize: " << fftSize << \
        ", sampleRate: " << sampleRate << \
        ", windowType: " << windowType);


  setup();
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Constructed");
}

SpectralReassignment::~SpectralReassignment(){}

void SpectralReassignment::setup(){
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Setting up...");
  
  // Setup the window so it gets calculated and can be reused
  _windowAlgo.setup();
  
  // Create the time vector
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Creating time vector...");
  Real timestep = 1.0 / _sampleRate;

  // The unit of the vectors is Time Sample fractions
  // So the difference between one coeff and the next is 1
  // and the center of the window must be 0, so even sized windows
  // will have the two center coeffs to -0.5 and 0.5
  // This should be a line going from [-(window_size - 1)/2 ... (window_size - 1)/2]
  _time.resize(_frameSize, 1);
  for(int i = 0; i < _time.rows(); i++){
    _time(i, 0) = (i - Real(_time.rows() - 1)/2.0);
  }
  
  // Create the freq vector
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Creating freq vector...");
  
  // The unit of the vectors is Frequency Bin fractions
  // TODO: Must rethink how the frequency vector is initialized
  // as we did for the time vector
  _freq.resize(1, _fftSize);
  range(0, _fftSize, _fftSize, &_freq);
  
  // Calculate and set the time weighted window
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Calculate time weighted window...");
  MatrixXR windowInteg = _windowAlgo.window();
  windowInteg = windowInteg.array() * _time.transpose().array();
  _windowIntegAlgo.setWindow(windowInteg);

  // Calculate and set the time derivated window
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Calculate time derivative window...");
  MatrixXR windowDeriv = _windowAlgo.window();
  for(int i = windowDeriv.cols() - 1; i > 0; i--){
    windowDeriv(0, i) = (windowDeriv(0, i) - windowDeriv(0, i - 1)) / timestep;
  }

  // TODO: Check what is the initial condition for the window
  // Should this be 0 or just the value it was originally * dt
  //windowDeriv(0, 0) = 0.0;
  _windowDerivAlgo.setWindow(windowDeriv);

  // Create the necessary buffers for the windowing
  _window.resize(1, _frameSize);
  _windowInteg.resize(1, _frameSize);
  _windowDeriv.resize(1, _frameSize);

  // Create the necessary buffers for the FFT
  _fftAbs2.resize(1, _fftSize);
  _fftInteg.resize(1, _fftSize);
  _fftDeriv.resize(1, _fftSize);
  
  // Setup the algos
  _windowIntegAlgo.setup();
  _windowDerivAlgo.setup();
  _fftAlgo.setup();

  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Finished set up...");
}


void SpectralReassignment::process(const MatrixXR& frames,
                                   MatrixXC* fft, MatrixXR* reassignTime, MatrixXR* reassignFreq){
  
  // Process the windowing
  _windowAlgo.process(frames, &_window);
  _windowIntegAlgo.process(frames, &_windowInteg);
  _windowDerivAlgo.process(frames, &_windowDeriv);
  
  // Process the FFT
  _fftAlgo.process(_window, fft);
  _fftAlgo.process(_windowInteg, &_fftInteg);
  _fftAlgo.process(_windowDeriv, &_fftDeriv);
  
  // Create the reassignment operations
  _fftAbs2 = (*fft).array().abs2();

  // Create the reassign operator matrix
  // TODO: check if the current timestamp is enough for a good reassignment
  // we might need for it to depend on past frames (if the reassignment of time
  // goes further than one)
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Processing: creating the time reassignment operation...");    
  (*reassignTime) = -((_fftInteg.array() * (*fft).conjugate().array()) / _fftAbs2.cast<Complex>().array()).real();
    
  // TODO: Check the unity of the freq reassignment, it may need to be normalized by something
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Processing: creating the freq reassignment operation...");
  (*reassignFreq) = _freq.array() + ((_fftDeriv.array() * (*fft).conjugate().array()) / _fftAbs2.cast<Complex>().array()).imag();
  
  (*reassignTime) = ((*reassignTime).array().isnan()).matrix().select(0, (*reassignTime));
  (*reassignFreq) = ((*reassignFreq).array().isnan()).matrix().select(0, (*reassignFreq));
  
  // Reassign the spectrum values
  // TODO: put this into a function and do it right
  // will have to calculate and return all the reassigned values:
  // reassignedFrequency, reassignedTime:
  //      - are calculated using Flandrin's method using the 3 DFT
  // reassignedMagnitude, reassignedPhase: 
  //      - are calculated from reassigned freq and time and the original DFT
  //        (the magnitude and phase must then be put back 
  //         in the form of a complex in the reassigned frame)
  /*
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Processing: reassigning...");
  LOUDIA_DEBUG("SPECTRALREASSIGNMENT: Processing: reassigning _reassignFreq: " << _reassignFreq.rows() << ", " << _reassignFreq.cols());
  
  for(int j = 0; j < _reassignFreq.cols(); j++){
    
    if((int)round(_reassignFreq(i, j)) >= 0 && (int)round(_reassignFreq(i, j)) < (*reassigned).cols()) {
      
      (*reassigned)(i, (int)round(_reassignFreq(i, j))) += ((1.0 - (abs(_reassignFreq(i, j) - (int)round(_reassignFreq(i,j))))) * abs(_fft(i, (int)round(_reassignFreq(i,j)))));
      
    }
  }
  */
}

void SpectralReassignment::reset(){
  _windowAlgo.reset();
  _windowIntegAlgo.reset();
  _windowDerivAlgo.reset();
  _fftAlgo.reset();
}

int SpectralReassignment::frameSize() const{
  return _frameSize;
}

int SpectralReassignment::fftSize() const{
  return _fftSize;
}

Window::WindowType SpectralReassignment::windowType() const{
  return _windowType;
}

SpectralNoiseSuppression::SpectralNoiseSuppression(int fftSize, Real lowFrequency, Real highFrequency, Real sampleRate)
{
  LOUDIA_DEBUG("SPECTRALNOISESUPPRESSION: Construction fftSize: " << fftSize
        << " sampleRate: " << sampleRate
        << " lowFrequency: " << lowFrequency
        << " highFrequency: " << highFrequency );

  setFftSize( fftSize, false );
  setSampleRate( sampleRate, false );
  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );

  setup();
}

SpectralNoiseSuppression::~SpectralNoiseSuppression(){}

void SpectralNoiseSuppression::setup(){
  LOUDIA_DEBUG("SPECTRALNOISESUPPRESSION: Setting up...");

  _k0 = (int)(( _lowFrequency / _sampleRate ) * _fftSize);
  _k1 = (int)(( _highFrequency / _sampleRate ) * _fftSize);
  
  // Prepare the bands for the moving average
  int _halfSize = (_fftSize / 2) + 1;

  int minHalfBand = (int)(100.0 / _sampleRate * _fftSize / 2.0);

  MatrixXI starts(_halfSize, 1);
  vector<MatrixXR> weights;
  weights.reserve(_halfSize);
  for ( int i = 0; i < _halfSize; i++ ) {
    int halfBandUnder = max( minHalfBand,  (int)(2.0 / 3.0 * i));
    int halfBandOver = max( minHalfBand,  (int)(i / 2.0 * 3.0));

    int begin = max( i - halfBandUnder, 0 );
    int end = min( i + halfBandOver, _halfSize - 1 );

    starts(i, 0) = begin;

    MatrixXR weight = MatrixXR::Constant(1, end - begin, 1.0 / float(end - begin));
    weights.push_back( weight );
  }
  
  _bands.setStartsWeights( starts, weights );
  _bands.setup();

  reset();

  LOUDIA_DEBUG("SPECTRALNOISESUPPRESSION: Finished setup.");
}

void SpectralNoiseSuppression::process(const MatrixXR& spectrum, MatrixXR* noise, MatrixXR* result){
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();

  (*result).resize(rows, cols);
  
  (*result) = spectrum;

  //DEBUG("SPECTRALNOISESUPPRESSION: Calculate wrapped magnitude.");
  // Calculate the wrapped magnitude of the spectrum
  _g = (1.0 / (_k1 - _k0 + 1.0) * spectrum.block(0, _k0, rows, _k1 - _k0).array().pow(1.0/3.0).rowwise().sum()).array().cube();

  //cout << (_g) << endl;
  
  for ( int i = 0; i < cols; i++ ) {
    (*result).col(i) = (((*result).col(i).array() * _g.array().inverse()) + 1.0).log();
  }
  
  //cout << (*result) << endl;
  
  //DEBUG("SPECTRALNOISESUPPRESSION: Estimate spectral noise.");
  // Estimate spectral noise
  _bands.process((*result), noise);
  
  //DEBUG("SPECTRALNOISESUPPRESSION: Suppress spectral noise.");
  // Suppress spectral noise
  (*result) = ((*result) - (*noise)).array().clipUnder();
}

void SpectralNoiseSuppression::reset(){
  // Initial values

  _bands.reset();
}

Real SpectralNoiseSuppression::lowFrequency() const{
  return _lowFrequency;
}
  
void SpectralNoiseSuppression::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real SpectralNoiseSuppression::highFrequency() const{
  return _highFrequency;
}
  
void SpectralNoiseSuppression::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real SpectralNoiseSuppression::sampleRate() const{
  return _sampleRate;
}
  
void SpectralNoiseSuppression::setSampleRate( Real frequency, bool callSetup ){
  _sampleRate = frequency;
  if ( callSetup ) setup();
}

int SpectralNoiseSuppression::fftSize() const{
  return _fftSize;
}

void SpectralNoiseSuppression::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}


SpectralWhitening::SpectralWhitening(int fftSize, Real lowFrequency, Real highFrequency, Real sampleRate, Real compressionFactor, int bandCount, MelBands::ScaleType scaleType)
{
  LOUDIA_DEBUG("SPECTRALWHITENING: Construction fftSize: " << fftSize
        << " sampleRate: " << sampleRate
        << " compressionFactor: " << compressionFactor
        << " bandCount: " << bandCount
        << " lowFrequency: " << lowFrequency
        << " highFrequency: " << highFrequency );

  setFftSize( fftSize, false );
  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );
  setBandCount( bandCount, false );
  setSampleRate( sampleRate, false );
  setScaleType( scaleType, false );
  setCompressionFactor( compressionFactor, false );

  setup();
}

SpectralWhitening::~SpectralWhitening(){}

void SpectralWhitening::setup(){
  LOUDIA_DEBUG("SPECTRALWHITENING: Setting up...");

  _halfSize = ( _fftSize / 2 ) + 1;
  
  // Setup the bands
  _bands.setLowFrequency( _lowFrequency, false );
  _bands.setHighFrequency( _highFrequency, false );
  _bands.setBandCount(_bandCount, false );
  _bands.setSampleRate( _sampleRate, false );
  _bands.setFftSize(_fftSize, false );
  _bands.setScaleType( _scaleType, false ); 
  _bands.setup();

  _bands.centers(&_centers);
  
  reset();

  LOUDIA_DEBUG("SPECTRALWHITENING: Finished setup.");
}

void SpectralWhitening::process(const MatrixXR& spectrum, MatrixXR* result){
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();

  (*result).resize( rows, cols );

  _compressionWeights.resize(rows, _halfSize);

  // Calculate the energy per band
  _bands.process( spectrum.array().square(), &_bandEnergy );

  // Calculate compress weights of bands
  _bandEnergy = (_bandEnergy / _fftSize).array().sqrt().pow(_compressionFactor - 1.0);

  // Interpolate linearly between the center frequencies of bands
  // Interpolate the region before the first frequency center
  int col = 0;
  for (; col < _centers(0, 0); col++ ) {
    _compressionWeights.col( col ) = ((Real)col * (_bandEnergy.col(0).array() - 1.0) / _centers(0, 0)) + 1.0;
  }

  // Interpolate the region between the first and last frequency centers
  for ( int band = 1; band < _bandCount; band++ ) {
    for (; col < _centers(band, 0); col++ ) {
      _compressionWeights.col(col) = (((Real)col - _centers(band - 1, 0)) * (_bandEnergy.col(band) - _bandEnergy.col(band-1)) / (_centers(band, 0) - _centers(band - 1, 0))) + _bandEnergy.col(band - 1);
    }
  }

  // Interpolate the region after the last frequency center
  for (; col < _halfSize; col++ ) {
      _compressionWeights.col(col) = (((Real)col - _centers(_bandCount - 1, 0)) * ( (-_bandEnergy).col(_bandCount - 1).array() + 1.0) / (_halfSize - _centers(_bandCount - 1, 0))) + _bandEnergy.col(_bandCount - 1).array();
  }

  // Apply compression weihgts
  (*result) = spectrum.array() * _compressionWeights.block(0, 0, spectrum.rows(), spectrum.cols()).array();
}

void SpectralWhitening::reset(){
  // Initial values

  _bands.reset();
}

Real SpectralWhitening::lowFrequency() const{
  return _lowFrequency;
}
  
void SpectralWhitening::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real SpectralWhitening::highFrequency() const{
  return _highFrequency;
}
  
void SpectralWhitening::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real SpectralWhitening::sampleRate() const{
  return _sampleRate;
}
  
void SpectralWhitening::setSampleRate( Real frequency, bool callSetup ){
  _sampleRate = frequency;
  if ( callSetup ) setup();
}

int SpectralWhitening::bandCount() const {
  return _bandCount;
}

void SpectralWhitening::setBandCount( int count, bool callSetup ) {
  _bandCount = count;
  if ( callSetup ) setup();
}

int SpectralWhitening::fftSize() const{
  return _fftSize;
}

void SpectralWhitening::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

MelBands::ScaleType SpectralWhitening::scaleType() const{
  return _scaleType;
}

void SpectralWhitening::setScaleType( MelBands::ScaleType type, bool callSetup ) {
  _scaleType = type;
  if ( callSetup ) setup();
}

Real SpectralWhitening::compressionFactor() const{
  return _compressionFactor;
}
  
void SpectralWhitening::setCompressionFactor( Real factor, bool callSetup ){
  _compressionFactor = factor;
  if ( callSetup ) setup();
}
}