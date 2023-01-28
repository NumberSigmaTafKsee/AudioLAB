#pragma once

namespace Loudia
{
    class SpectralODFBase {
    protected:
        // Internal parameters
        int _fftSize;
        int _halfSize;
    
    // Internal variables
    public:
        void setup();
        virtual void reset() = 0;
        
        virtual void process(const MatrixXC& fft, MatrixXR* odfValue) = 0;
        
        virtual ~SpectralODFBase() { };

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


    /**
    * @class SpectralODF
    *
    * @brief Algorithm to perform onset detection functions on vectors of Complex values representing the FFT of a signal.
    *
    * This class wraps several kind of spectral Onset Detection Functions (ODF).  A spectral onset detection function is 
    * a function mapping from an FFT of an audio signal to a Real value.  The onset detection function represents a value
    * proportional to the probability of the frame being an onset.
    *
    * The algorithm takes as input N-point vectors of Complex values 
    * and returns Real values.
    *
    * 5 types of ODF methods are implemented:
    * -# Spectral flux
    * -# High frequency content
    * -# Phase deviation
    * -# Weighted phase deviation
    * -# Normalized weighted phase deviation
    * -# Modified Kullback-Liebler
    * -# Complex domain
    * -# Rectified complex domain
    * -# Peak center of gravity
    *
    * The Phase deviation, Weighted phase deviation, Normalized weighted phase deviation, Complex domain and 
    * Rectified complex domain methods require the 2 past FFT frames.
    * Therefore if any of these methods are specified, the process() method would require an input matrix of at least 3 rows 
    * and will output ODF values for all the rows but the first two.
    *
    * The Spectral flux method requires the past FFT frame.
    * Therefore if this method is used, the process() method would require an input matrix of at least 2 rows 
    * and will output ODF values for all the rows but the first.
    * 
    * The ODF method can be selected using the 
    * setOdfMethod() taking as argument an ODFMethod.
    *
    * This function is often use to perform beat estimation and event segmentation tasks.
    *
    * @author Ricard Marxer
    *
    * @sa FFT
    */
    class SpectralODF : SpectralODFBase {
    public:
        /**
            @enum ODFMethod
            @brief Specifies the method for calculating the onset detection function to be used.
            @sa odfMethod
        */
        enum ODFMethod {
            FLUX = 0                          /**< Spectral flux method  */,
            HIGH_FREQUENCY_CONTENT = 1        /**< High frequency content method  */,
            PHASE_DEVIATION = 2               /**< Phase deviation method  */,
            WEIGHTED_PHASE_DEVIATION = 3      /**< Weighted phase deviation method  */,
            NORM_WEIGHTED_PHASE_DEVIATION = 4 /**< Normalized weighted phase deviation method  */,
            MODIFIED_KULLBACK_LIEBLER = 5     /**< Modified Kullback-Liebler method  */,
            COMPLEX_DOMAIN = 6                /**< Complex domain method  */,
            RECTIFIED_COMPLEX_DOMAIN = 7      /**< Rectified complex domain method  */,
            CENTER_OF_GRAVITY = 8             /**< Peak center of gravity method  */
        };

    protected:
        // Internal parameters
        ODFMethod _odfMethod;
        
        // Internal variables
        SpectralODFBase* _odf;

    public:
        /**
            Constructs a spectral onset detection function object with the specified @a fftSize and 
            @a odfMethod settings.
            
            @param fftSize size of the input FFT frames, must be > 0.
            
            @param odfMethod the onset detection method to be used
        */
        SpectralODF(int fftSize = 1024, ODFMethod odfMethod = COMPLEX_DOMAIN);
        
        /**
            Destroys the algorithm and frees its resources.
        */
        ~SpectralODF();

        void setup();
        void reset();

        /**
            Calculates the onset detection method on each of the rows of @a ffts and
            puts the resulting onset detection function values in the rows of @a odfs.
            
            @param ffts matrix of Complex values.  The number of columns of @a ffts 
            must be equal to the fftSize / 2 + 1.  Some onset detection methods require
            a minimum of 2 (or 3) rows to calculate the one onset detection function value.
            In this cases onset detection values for the first 2 (or 3 respectively) rows
            will not be output.
            
            @param odfs pointer to a single-column matrix of Real values for the output.  The matrix should
            have the same number of rows as @a ffts (minus 1 or 2 depending on the method used) and 1 single column. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXC& ffts, MatrixXR* odfs);
        
        /**
            Returns the method for calculating the spectral onset detection function.
            
            By default it is COMPLEX_DOMAIN.
        */
        ODFMethod odfMethod() const;

        /**
            Specifies the @a method for calculating the spectral onset detection function.
            
            Note that PHASE_DEVIATION, WEIGHTED_PHASE_DEVIATION, NORM_WEIGHTED_PHASE_DEVIATION, 
            COMPLEX_DOMAIN and RECTIFIED_COMPLEX_DOMAIN methods require at least 3 FFT frames, and therefore
            the input matrix to the process() method must have at least 3 rows.  In these cases
            the output matrix will be 2 rows smaller than the input, since it only calculates ODF values
            for all the rows but the first two.
            Note that SPECTRAL_FLUX method requires at least 2 FFT frames, and therefore
            the input matrix to the process() method must have at least 2 rows.  In this case
            the output matrix will be 1 row smaller than the input, since it only calculates ODF values
            for all the rows but the first.
            
            @param method the method used for calculating the spectral onset detection function.
            
            @param callSetup a flag specifying whether the setup() method must be called after setting the parameter.
        */
        void setOdfMethod( ODFMethod method, bool callSetup = true );
    };

    class SpectralODFComplex : public SpectralODFBase {
    protected:
        // Internal parameters
        int _fftSize;
        int _halfSize;
        bool _rectified;
        
        // Internal variables
        Unwrap _unwrap;

        MatrixXC _unwrappedSpectrum;
        MatrixXR _unwrappedAngle;
        MatrixXC _spectrumPredict;
        MatrixXR _predictionError;
        
        void spectralDistanceEuclidean(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg, MatrixXR* odfValues);

    public:
        SpectralODFComplex(int fftSize, bool rectified = false);

        ~SpectralODFComplex();

        void setup();
        void reset();

        void process(const MatrixXC& fft, MatrixXR* odfValue);
    };

    class SpectralODFCOG : public SpectralODFBase {
    protected:
        // Internal parameters
        int _fftSize;
        int _peakCount;
        int _bandwidth;

        // Internal variables
        MatrixXR _spectrumAbs2;
        MatrixXR _spectrumArg;
        MatrixXR _spectrumArgDeriv;

        MatrixXR _peakStarts;
        MatrixXR _peakPos;
        MatrixXR _peakEnds;
        MatrixXR _peakMag;

        MatrixXR _cog;

        PeakDetection _peaker;
        PeakCOG _peakCoger;

    public:
        SpectralODFCOG(int fftSize, int bandwidth = 1, int peakCount = -1);

        ~SpectralODFCOG();

        void setup();

        void process(const MatrixXC& fft, MatrixXR* odfValue);

        void reset();

    };

    
    class SpectralODFFlux : public SpectralODFBase {
    protected:
        // Internal parameters
        int _fftSize;
        
        // Internal variables
        MatrixXR _spectrumAbs;

    public:
        SpectralODFFlux(int fftSize);

        ~SpectralODFFlux();

        void setup();

        void process(const MatrixXC& fft, MatrixXR* odfValue);

        void reset();
    };

    class SpectralODFHFC : public SpectralODFBase {
    protected:
        // Internal parameters  

        // Internal variables
        MatrixXR _spectrumAbs;
        MatrixXR _freqBin;

    public:
        SpectralODFHFC( int fftSize = 1024 );

        ~SpectralODFHFC();

        void setup();
        void reset();

        void process(const MatrixXC& fft, MatrixXR* odfValue);

    };    

    class SpectralODFMKL : public SpectralODFBase {
    protected:
        // Internal parameters
        int _fftSize;
        Real _minSpectrum;
        
        // Internal variables
        MatrixXR _spectrumAbs;

    public:
        SpectralODFMKL(int fftSize, Real _minSpectrum = 1e-7);

        ~SpectralODFMKL();

        void setup();

        void process(const MatrixXC& fft, MatrixXR* odfValue);

        void reset();

    };

    class SpectralODFPhase : public SpectralODFBase {
    protected:
        // Internal parameters
        int _fftSize;
        int _halfSize;
        bool _weighted;
        bool _normalize;
        
        // Internal variables
        Unwrap _unwrap;
        
        MatrixXR _unwrappedAngle;
        MatrixXR _phaseDiff;
        MatrixXR _instFreq;
        
        void phaseDeviation(const MatrixXC& spectrum, const MatrixXR& spectrumArg, MatrixXR* odfValue);

    public:
        SpectralODFPhase(int fftSize, bool weighted = false, bool normalize = false);

        ~SpectralODFPhase();

        void setup();

        void process(const MatrixXC& fft, MatrixXR* odfValue);

        void reset();

    };


void SpectralODFBase::setup() {
  _halfSize = _fftSize / 2 + 1;
}

int SpectralODFBase::fftSize() const{
  return _fftSize;
}

void SpectralODFBase::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

SpectralODF::SpectralODF(int fftSize, ODFMethod odfMethod) :
  SpectralODFBase(),
  _odf( 0 )
{
  setFftSize( fftSize, false );
  setOdfMethod( odfMethod, false );

  setup();
}

SpectralODF::~SpectralODF() {
  delete _odf;
  _odf = 0;
}

void SpectralODF::setup() {

  delete _odf;
  _odf = 0;

  switch( _odfMethod ) {

  case FLUX:
    _odf = new SpectralODFFlux(_fftSize);
    break;

  case PHASE_DEVIATION:
    _odf = new SpectralODFPhase(_fftSize);
    break;

  case WEIGHTED_PHASE_DEVIATION:
    _odf = new SpectralODFPhase(_fftSize, true);
    break;

  case NORM_WEIGHTED_PHASE_DEVIATION:
    _odf = new SpectralODFPhase(_fftSize, true, true);
    break;

  case MODIFIED_KULLBACK_LIEBLER:
    _odf = new SpectralODFMKL(_fftSize);
    break;

  case COMPLEX_DOMAIN:
    _odf = new SpectralODFComplex(_fftSize);
    break;

  case RECTIFIED_COMPLEX_DOMAIN:
    _odf = new SpectralODFComplex(_fftSize, true);
    break;

  case HIGH_FREQUENCY_CONTENT:
    _odf = new SpectralODFHFC(_fftSize);
    break;

  case CENTER_OF_GRAVITY:
    _odf = new SpectralODFCOG(_fftSize);
    break;

  }

  _odf->setup();
}

void SpectralODF::process(const MatrixXC& fft, MatrixXR* odfValue) {
  if (!_odf) return;
  
  _odf->process(fft, odfValue);
}

void SpectralODF::reset() {
  if (!_odf) return;
  
  _odf->reset();
}

SpectralODF::ODFMethod SpectralODF::odfMethod() const{
  return _odfMethod;
}

void SpectralODF::setOdfMethod( ODFMethod method, bool callSetup ) {
  _odfMethod = method;
  if ( callSetup ) setup();
}

SpectralODFCOG::SpectralODFCOG(int fftSize, int peakCount, int bandwidth) :
  SpectralODFBase(),
  _fftSize( fftSize ),
  _peakCount( peakCount ),
  _bandwidth( bandwidth ),
  _peaker( peakCount, PeakDetection::BYMAGNITUDE, bandwidth ),
  _peakCoger( fftSize, bandwidth )
{
  
  LOUDIA_DEBUG("SPECTRALODFCOG: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFCOG::~SpectralODFCOG() {}


void SpectralODFCOG::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFCOG: Setting up...");

  _peaker.setup();
  _peakCoger.setup();
  
  reset();

  LOUDIA_DEBUG("SPECTRALODFCOG: Finished set up...");
}


void SpectralODFCOG::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFCOG: Processing windowed");
  const int rows = fft.rows();
  
  (*odfValue).resize(rows, 1);

  LOUDIA_DEBUG("SPECTRALODFCOG: Processing the peaks");

  _peaker.process(fft.array().abs(), &_peakStarts, &_peakPos, &_peakEnds, &_peakMag);

  _peakCoger.process(fft, _peakPos, &_cog);

  (*odfValue) = _cog.array().clipUnder().rowwise().sum();
  
  LOUDIA_DEBUG("SPECTRALODFCOG: Finished Processing");
}

void SpectralODFCOG::reset() {
  // Initial values
  _peaker.reset();
  _peakCoger.reset();

}

SpectralODFComplex::SpectralODFComplex(int fftSize, bool rectified) :
  SpectralODFBase(),
  _fftSize( fftSize ),
  _halfSize( _fftSize / 2 + 1 ),
  _rectified( rectified )
{
  
  LOUDIA_DEBUG("SPECTRALODFCOMPLEX: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFComplex::~SpectralODFComplex() {}


void SpectralODFComplex::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFCOMPLEX: Setting up...");

  _unwrap.setup();

  reset();

  LOUDIA_DEBUG("SPECTRALODFCOMPLEX: Finished set up...");
}


void SpectralODFComplex::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFCOMPLEX: Processing windowed");
  const int rows = fft.rows();
  
  if ( rows < 3 ) {
    // Throw ValueError, it must have a minimum of 3 rows
  }

  (*odfValue).resize(rows - 2, 1);

  _unwrap.process(fft.array().angle(), &_unwrappedAngle);

  LOUDIA_DEBUG("SPECTRALODFCOMPLEX: Processing unwrapped");
  
  spectralDistanceEuclidean(fft, fft.array().abs(), _unwrappedAngle, odfValue);
  
  LOUDIA_DEBUG("SPECTRALODFCOMPLEX: Finished Processing");
}

void SpectralODFComplex::spectralDistanceEuclidean(const MatrixXC& spectrum, const MatrixXR& spectrumAbs, const MatrixXR& spectrumArg, MatrixXR* odfValue) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  _spectrumPredict.resize(rows - 2, cols);
  _predictionError.resize(rows - 2, cols);

  polar(spectrumAbs.block(1, 0, rows - 2, cols), 2.0 * spectrumArg.block(1, 0, rows - 2, cols) - spectrumArg.block(0, 0, rows - 2, cols), &_spectrumPredict);
  
  _predictionError = (_spectrumPredict - spectrum.block(0, 0, rows - 2, cols)).array().abs();

  if (_rectified) {
    _predictionError = (_spectrumPredict.array().abs() <= spectrum.block(0, 0, rows - 2, cols).array().abs()).select(_predictionError, 0.0);
  }
  
  //_predictionError.col(0) = 0.0;
  
  (*odfValue) = _predictionError.rowwise().sum() / cols;
  return;
}

void SpectralODFComplex::reset() {
  // Initial values
  _unwrap.reset();
}

SpectralODFFlux::SpectralODFFlux(int fftSize) :
  SpectralODFBase(),
  _fftSize( fftSize )
{
  
  LOUDIA_DEBUG("SPECTRALODFFLUX: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFFlux::~SpectralODFFlux() {}


void SpectralODFFlux::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFFLUX: Setting up...");

  reset();

  LOUDIA_DEBUG("SPECTRALODFFLUX: Finished set up...");
}


void SpectralODFFlux::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFFLUX: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();

  if ( rows < 2 ) {
    // Throw ValueError, it must have a minimum of 2 rows
  }

  (*odfValue).resize(rows - 1, 1);
  (*odfValue).row(0).setZero();
  
  LOUDIA_DEBUG("SPECTRALODFFLUX: Spectrum resized rows: " << rows );
    
  _spectrumAbs = fft.array().abs();
  
  (*odfValue) = (_spectrumAbs.block(1, 0, rows - 1, cols) \
                 - _spectrumAbs.block(0, 0, rows - 1, cols) \
                 ).array().clipUnder().rowwise().sum() / cols;
  
  LOUDIA_DEBUG("SPECTRALODFFLUX: Finished Processing");
}

void SpectralODFFlux::reset() {
  // Initial values
}

SpectralODFHFC::SpectralODFHFC(int fftSize) :
  SpectralODFBase()
{
  
  LOUDIA_DEBUG("SPECTRALODFHFC: Constructor fftSize: " << _fftSize);
  setFftSize( fftSize, false );
  setup();
}

SpectralODFHFC::~SpectralODFHFC() {}


void SpectralODFHFC::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFHFC: Setting up...");

  SpectralODFBase::setup();

  // Create the vector with the weights (weights are the frequency bin indices)
  range(0, _halfSize, _halfSize, &_freqBin);
  
  reset();

  LOUDIA_DEBUG("SPECTRALODFHFC: Finished set up...");
}


void SpectralODFHFC::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFHFC: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();

  (*odfValue).resize(rows, 1);
  
  _spectrumAbs = fft.array().abs();

  for (int row = 0; row < rows; row ++) {
    (*odfValue).row(row) = _spectrumAbs.row(row) * _freqBin.block(0, 0, 1, cols).transpose() / cols;
  }
  
  LOUDIA_DEBUG("SPECTRALODFHFC: Finished Processing");
}

void SpectralODFHFC::reset() {
  // Initial values
}

SpectralODFMKL::SpectralODFMKL(int fftSize, Real minSpectrum) :
  SpectralODFBase(),
  _fftSize(fftSize),
  _minSpectrum(minSpectrum)
{
  
  LOUDIA_DEBUG("SPECTRALODFMKL: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFMKL::~SpectralODFMKL() {}


void SpectralODFMKL::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFMKL: Setting up...");

  reset();

  LOUDIA_DEBUG("SPECTRALODFMKL: Finished set up...");
}


void SpectralODFMKL::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFMKL: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  
  if ( rows < 2 ) {
    // Throw ValueError, it must have a minimum of 2 rows
  }

  (*odfValue).resize(rows - 1, 1);

  LOUDIA_DEBUG("SPECTRALODFMKL: Spectrum resized rows: " << rows);
  
  _spectrumAbs = fft.array().abs();

  (*odfValue) = (_spectrumAbs.block(1, 0, rows-1, cols).array() \
                 * (_spectrumAbs.block(1, 0, rows-1, cols).array() \
                    / (_spectrumAbs.block(0, 0, rows-1, cols).array().clipUnder(_minSpectrum))).array().clipUnder(_minSpectrum).logN(2.0)).rowwise().sum() / cols;
  
  LOUDIA_DEBUG("SPECTRALODFMKL: Finished Processing");
}

void SpectralODFMKL::reset() {
  // Initial values
}

SpectralODFPhase::SpectralODFPhase(int fftSize, bool weighted, bool normalize) :
  SpectralODFBase(),
  _fftSize( fftSize ),
  _halfSize( _fftSize / 2 + 1 ),
  _weighted( weighted ),
  _normalize( normalize )
{
  
  LOUDIA_DEBUG("SPECTRALODFPHASE: Constructor fftSize: " << _fftSize);
  
  setup();
}

SpectralODFPhase::~SpectralODFPhase() {}


void SpectralODFPhase::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("SPECTRALODFPHASE: Setting up...");

  _unwrap.setup();

  reset();

  LOUDIA_DEBUG("SPECTRALODFPHASE: Finished set up...");
}


void SpectralODFPhase::process(const MatrixXC& fft, MatrixXR* odfValue) {
  LOUDIA_DEBUG("SPECTRALODFPHASE: Processing windowed");
  const int rows = fft.rows();
  
  if ( rows < 3 ) {
    // Throw ValueError, it must have a minimum of 3 rows
  }

  (*odfValue).resize(rows - 2, 1);
  
  _unwrap.process(fft.array().angle(), &_unwrappedAngle);

  LOUDIA_DEBUG("SPECTRALODFPHASE: Processing unwrapped");
  
  phaseDeviation(fft, _unwrappedAngle, odfValue);
  
  LOUDIA_DEBUG("SPECTRALODFPHASE: Finished Processing");
}

void SpectralODFPhase::phaseDeviation(const MatrixXC& spectrum, const MatrixXR& spectrumArg, MatrixXR* odfValue) {
  const int rows = spectrum.rows();
  const int cols = spectrum.cols();
  
  _phaseDiff = spectrumArg.block(1, 0, rows - 1, cols) - spectrumArg.block(0, 0, rows - 1, cols);
  _instFreq = _phaseDiff.block(1, 0, rows - 2, cols) - _phaseDiff.block(0, 0, rows - 2, cols);

  if (_weighted)
    _instFreq.array() *= spectrum.block(2, 0, rows - 2, cols).array().abs();

  if (_normalize) {
    (*odfValue) = _instFreq.rowwise().sum().array() / (cols * spectrum.block(2, 0, rows - 2, cols).array().abs().rowwise().sum());
    return;
  }
  
  (*odfValue) = _instFreq.rowwise().sum() / cols;
  return;
}


void SpectralODFPhase::reset() {
  // Initial values
  _unwrap.reset();
}

}