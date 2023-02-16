#pragma once

namespace Loudia
{
    /**
    * @class PitchACF
    *
    * @brief Algorithm to estimate the most prominent pitch of a vector of
    * Real values reprsenting the FFT of an audio frame using the 
    * Autocorrelation function.
    *
    * This class represents an object to estimate the most prominent pitch 
    * of a vector of Real values reprsenting the FFT of an audio frame using 
    * the Autocorrelation function.
    *
    * The algorithm performs an autocorrelation on the input spectrums
    * and finds the first peak of a set of highest magnitude candidates.
    * The index of the peak and the values of the peak will determine the
    * pitch frequency and saliency.
    * 
    * The minimum peak width can be specified using setMinimumPeakWidth().
    * 
    * The number of candidates at the peak detection stage be
    * specified using setPeakCandidateCount().
    * 
    *
    * @author Ricard Marxer
    *
    * @sa PitchACF, PitchSaliency, PitchInvp
    */
    class PitchACF {
    protected:
        int _fftSize;
        int _halfSize;
        int _minimumPeakWidth;
        int _peakCandidateCount;

        Real _sampleRate;

        PeakDetection _peak;
        PeakInterpolation _peakInterp;
        Autocorrelation _acorr;

        MatrixXR _acorred;

        MatrixXR _starts;
        MatrixXR _ends;

    public:
        /**
            Constructs an autocorrelation based pitch estimation function 
            object with the specified @a fftSize, @a sampleRate, @a minPeakWidth
            and @a peakCandidateCount settings.
            
            @param fftSize size of the input FFT frames, must be > 0.
            
            @param sampleRate the sampleRate of the input signal.  By default
            it is 1.0, so the pitches will be returned in normalized frequencies.
            @param minimumPeakWidth the minimum width of a peak in the autocorrelation
            function for it to be detected.
            @param peakCandidateCount the number of highest magnitude candidates 
            to be considered during the peak detection process of the 
            autocorrelation function.
        */
        PitchACF(int fftSize = 1024, Real sampleRate = 1.0, int minimumPeakWidth = 6, int peakCandidateCount = 10);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~PitchACF();

        void setup();
        void reset();

        /**
            Performs a pitch estimation on each of the rows of @a spectrums.
            Puts the resulting estimated pitches in the rows of @a pitches
            and the saliencies of each pitch in the rows of @a saliencies.
            
            @param spectrums matrix of Real values representing one 
            spectrum magnitude per row.
            The number of columns of @a spectrum must be equal 
            to the fftSize / 2 + 1 where 
            fftSize is specified using setFftSize().
            
            @param pitches pointer to a matrix of Real values representing 
            the frequencies of the estimated pitches as rows.
            The matrix should have the same number of rows as @a spectrums and as many
            columns as the count of estimated pitches.  Note that this algorithm is
            only capable of detecting a single pitch at each frame, therefore @a pitches
            will be a single column matrix.
            @param saliencies pointer to a matrix of Real values representing
            the saliencies of the estimated pitches as rows.
            The matrix should have the same number of rows as @a spectrums and as many
            columns as the count of estimated pitches.  Note that this algorithm is
            only capable of detecting a single pitch at each frame, therefore @a pitches
            will be a single column matrix.
            
            Note that if the output matrices are not of the required sizes they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& spectrums, MatrixXR* pitches, MatrixXR* saliencies);
        
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
            Returns the minimum width for a peak to be detected in the
            autocorrelation function.
            The default is 6.
            
            @sa setMinimumPeakWidth()
        */
        int minimumPeakWidth() const;
        
        /**
            Specifies the minimum @a width for a peak to be detected in the
            autocorrelation function.
            
            @sa minimumPeakWidth()
        */
        void setMinimumPeakWidth( int width, bool callSetup = true );

        /**
            Returns the number of highest magnitude candidates to be considered 
            during the peak detection process of the autocorrelation function.
            Note that if the value is <= 0, then no preselection is performed
            and all detected peaks are considered as candidates.
            By default it is 6.
        */
        int peakCandidateCount() const;

        /**
            Specifies the number of highest magnitude candidates to be considered 
            during the peak detection process of the autocorrelation function.
            Note that if the value is <= 0, then no preselection is performed
            and all detected peaks are considered as candidates.
        */
        void setPeakCandidateCount( int count, bool callSetup = true );

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

    };


    class PitchInverseProblem {
    protected:
        int _fftSize;
        int _halfSize;
        Real _lowFrequency;
        Real _highFrequency;
        int _pitchCount;
        int _harmonicCount;
        int _frequencyCandidateCount;
        Real _peakWidth;

        Real _sampleRate;

        // Crop points of the frequency spectrum
        // In order to reduce the size of the inverse
        // projection matrix
        Real _lowCutFrequency;
        Real _highCutFrequency;
        int _lowBin;
        int _highBin;
        int _range;

        Real _tMin;
        Real _tMax;
        Real _alpha;
        Real _beta;
        Real _inharmonicity;
        Real _regularisation;

        MatrixXR _projectionMatrix;
        MatrixXR _inverseProjectionMatrix;

        MatrixXR _starts;
        MatrixXR _ends;

        PeakDetection _peak;
        PeakInterpolation _peakInterp;

        void harmonicWeight(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result);
        void harmonicSpread(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result);
        void harmonicPosition(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result);
        Real harmonicWeight(Real f, Real fMin, Real fMax, int harmonicIndex);
        Real harmonicSpread(Real f, Real fMin, Real fMax, int harmonicIndex);
        Real harmonicPosition(Real f, Real fMin, Real fMax, int harmonicIndex);

    public:
        PitchInverseProblem(int fftSize = 1024, Real lowFrequency = 50.0, Real highFrequency = 2100.0, Real sampleRate = 44100.0, int pitchCount = 5, int harmonicCount = 10, int frequencyCandidateCount = -1, Real peakWidth = 8);

        ~PitchInverseProblem();

        void reset();
        void setup();

        void process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies, MatrixXR* frequencies);

        void projectionMatrix(MatrixXR* matrix) const;

        /**
            Return the lowest frequency candidate.
            The default is 50.0.
            @sa lowFrequency, highFrequency, setLowFrequency, setHighFrequency
        */  
        Real lowFrequency() const;  

        /**
            Specifies the lowest @a frequency candidate.
            The given @a frequency must be in the range of 0 to the sampleRate / 2.
            
            @sa lowFrequency, highFrequency, setHighFrequency
        */
        void setLowFrequency( Real frequency, bool callSetup = true );

        /**
            Return the highest frequency candidate.
            The default is 2100.0.
            @sa lowFrequency, setLowFrequency, setHighFrequency
        */  
        Real highFrequency() const;  

        /**
            Specifies the highest @a frequency candidate.
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
            Returns the count of candidate frequencies in with which to discretize
            the frequency space.
            Note that if the value is <= 0, then fftSize / 2 + 1 is used.
            By default it is 6.
            @sa setFrequencyCandidateCount()
        */
        int frequencyCandidateCount() const;

        /**
            Specifies the number of highest magnitude candidates to be considered 
            during the frequency detection process of the autocorrelation function.
            Note that if the value is <= 0, then no preselection is performed
            and all detected frequencys are considered as candidates.
            @sa frequencyCandidateCount()
        */
        void setFrequencyCandidateCount( int count, bool callSetup = true );

        /**
            Returns the width of the harmonic peaks.
            The default is 8.
            
            @sa setPeakWidth()
        */
        int peakWidth() const;
        
        /**
            Specifies the @a width of the harmonic peaks.
            
            @sa peakWidth()
        */
        void setPeakWidth( int width, bool callSetup = true );

        /**
            Returns the maximum count of pitches to be estimated.
            By default it is 5.
            @sa setPitchCount()
        */
        int pitchCount() const;

        /**
            Specifies the maximum @a count of pitches to be estimated.
            @sa pitchCount()
        */
        void setPitchCount( int count, bool callSetup = true );

        /**
            Returns the maximum count of harmonics to be rendered in the projection matrix.
            By default it is 10.
            @sa setHarmonicCount()
        */
        int harmonicCount() const;

        /**
            Specifies the @a count of harmonics to be rendered in the projection matrix.
            @sa harmonicCount()
        */
        void setHarmonicCount( int count, bool callSetup = true );
    };

    class PitchSaliency {
    protected:
        int _fftSize;
        int _halfSize;
        Real _f0;
        Real _f1;
        Real _fPrec;
        int _numHarmonics;

        Real _tMin;
        Real _tMax;
        Real _tPrec;
        Real _alpha;
        Real _beta;

        Real _sampleRate;

        Real harmonicWeight(Real period, Real tLow, Real tUp, int harmonicIndex);

        Real saliency(Real period, Real deltaPeriod, Real tLow, Real tUp, const MatrixXR& spectrum);


    public:
        PitchSaliency(int fftSize, Real f0, Real f1, Real sampleRate = 1.0, Real fPrec = 0.01, int numHarmonics = 5);

        ~PitchSaliency();

        void setup();

        void process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies);

        void reset();
    };

PitchACF::PitchACF(int fftSize, Real sampleRate, int minimumPeakWidth, int peakCandidateCount)
{
  LOUDIA_DEBUG("PITCHACF: Construction fftSize: " << _fftSize
        << " sampleRate: " << _sampleRate );

  setFftSize( fftSize, false );
  setSampleRate( sampleRate, false );
  setMinimumPeakWidth( minimumPeakWidth, false );
  setPeakCandidateCount( peakCandidateCount, false );
  setup();
}

PitchACF::~PitchACF(){}

void PitchACF::setup(){
  LOUDIA_DEBUG("PITCHACF: Setting up...");

  _halfSize = ( _fftSize / 2 ) + 1;

  _peak.setPeakCount( 1, false );
  _peak.setSortMethod( PeakDetection::BYMAGNITUDE, false );
  _peak.setMinimumPeakWidth( _minimumPeakWidth, false );
  _peak.setCandidateCount( _peakCandidateCount, false );
  _peak.setup();

  _acorr.setInputSize( _halfSize, false );
  _acorr.setMaxLag( _halfSize, false );
  _acorr.setup();

  reset();

  LOUDIA_DEBUG("PITCHACF: Finished setup.");
}

void PitchACF::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies){
  _acorr.process(spectrum, &_acorred);
  
  _peak.process(_acorred,
                &_starts, pitches, &_ends, saliencies);
  
  _peakInterp.process(_acorred, (*pitches), (*saliencies),
                      pitches, saliencies);

  (*pitches) *= 2.0 * _sampleRate / _fftSize;
  
  (*saliencies).array() /= _acorred.col(0).array();
}

void PitchACF::reset(){
  // Initial values

}

int PitchACF::fftSize() const{
  return _fftSize;
}

void PitchACF::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

int PitchACF::minimumPeakWidth() const{
  return _minimumPeakWidth;
}

void PitchACF::setMinimumPeakWidth( int width, bool callSetup ) {
  _minimumPeakWidth = width;
  if ( callSetup ) setup();
}

int PitchACF::peakCandidateCount() const{
  return _peakCandidateCount;
}

void PitchACF::setPeakCandidateCount( int count, bool callSetup ) {
  _peakCandidateCount = count;
  if ( callSetup ) setup();
}

Real PitchACF::sampleRate() const{
  return _sampleRate;
}
  
void PitchACF::setSampleRate( Real frequency, bool callSetup ){
  _sampleRate = frequency;
  if ( callSetup ) setup();
}

PitchInverseProblem::PitchInverseProblem(int fftSize, Real lowFrequency, Real highFrequency, Real sampleRate, int pitchCount, int harmonicCount, int frequencyCandidateCount, Real peakWidth)
{
  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Construction fftSize: " << fftSize
        << " sampleRate: " << sampleRate
        << " lowFrequency: " << lowFrequency
        << " highFrequency: " << highFrequency
        << " pitchCount: " << pitchCount
        << " harmonicCount: " << harmonicCount
        << " frequencyCandidateCount: " << frequencyCandidateCount
        << " peakWidth: " << peakWidth);


  setFftSize( fftSize, false );
  setLowFrequency( lowFrequency, false );
  setHighFrequency( highFrequency, false );
  setPitchCount( pitchCount, false );
  setHarmonicCount( harmonicCount, false );
  setFrequencyCandidateCount( frequencyCandidateCount, false );
  setPeakWidth( peakWidth, false );
  setSampleRate( sampleRate, false );
  setup();
}

PitchInverseProblem::~PitchInverseProblem(){
}

void PitchInverseProblem::setup(){
  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting up...");

  _halfSize = ( _fftSize / 2 ) + 1;

  _peak.setPeakCount( _pitchCount, false );
  _peak.setSortMethod( PeakDetection::BYMAGNITUDE, false );
  _peak.setup();

  _peakInterp.setup();

  // Define the range that will be used
  _lowCutFrequency = 0;
  _highCutFrequency = 3000;
  _lowBin = (int)(_lowCutFrequency / _sampleRate * _fftSize);
  _highBin = std::min((int)(_highCutFrequency / _sampleRate * _fftSize), _halfSize);
  _range = _highBin - _lowBin;

  int frequencyCount = _frequencyCandidateCount;
  if (frequencyCount <= 0) {
    frequencyCount = _range;
  }

  _regularisation = 3.0;

  // Params taken from Klapuri ISMIR 2006
  _alpha = 27; // 27 Hz
  _beta = 320; // 320 Hz
  _inharmonicity = 0.0;

  MatrixXR freqs;
  range(_lowFrequency, _highFrequency, frequencyCount, &freqs);

  _projectionMatrix.resize(_range, freqs.cols());  // We add one that will be the noise component
  _projectionMatrix.setZero();

  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting up the projection matrix...");
  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int col = 0; col < _projectionMatrix.cols(); col++ ) {
      Real f = freqs(0, col);

      for ( int harmonicIndex = 1; harmonicIndex < _harmonicCount + 1; harmonicIndex++ ) {
        Real mu = harmonicPosition(f, _lowFrequency, _highFrequency, harmonicIndex);
        Real a = harmonicWeight(f, _lowFrequency, _highFrequency, harmonicIndex);
        Real fi = harmonicSpread(f, _lowFrequency, _highFrequency, harmonicIndex);

        _projectionMatrix(row, col) += a * gaussian(row+0.5+_lowBin, mu, fi);
      }
    }
  }
  
  /*
  MatrixXR gauss;
  MatrixXR mu;
  MatrixXR a;
  MatrixXR fi;
  //MatrixXR x;
  //range(0, _frequencyCandidateCount, frequencyCount, _projectionMatrix.rows(), &x);
  
  for ( int row = 0; row < _projectionMatrix.rows(); row++ ) {
    for ( int harmonicIndex = 1; harmonicIndex < _harmonicCount + 1; harmonicIndex++ ) {
      harmonicPosition(freqs, _lowFrequency, _highFrequency, harmonicIndex, &mu);
      harmonicWeight(freqs, _lowFrequency, _highFrequency, harmonicIndex, &a);
      harmonicSpread(freqs, _lowFrequency, _highFrequency, harmonicIndex, &fi);
      gaussian(row, mu, fi, &gauss);
      
      _projectionMatrix.row(row) += a.cwise() * gauss;
    }
  }
  */


  /*                                                                                                   
  // Use the LU inverse                                                                                
   MatrixXR sourceWeight = MatrixXR::Identity( frequencyCount, frequencyCount );                        
   MatrixXR targetWeight = MatrixXR::Identity( _range, _range );                                        
                                                                                                        
   // Since the source weights is a diagonal matrix, the inverse is the inverse of the diagonal         
   MatrixXR invSourceWeight = sourceWeight;                                                             
   invSourceWeight.diagonal() = invSourceWeight.diagonal().cwise().inverse();                           
   LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting up the inversion...");                                    
   // A = W^{-1} K^t [ K W^{-1} K^t + \lambda * I_N ]^{+}                                               
   _inverseProjectionMatrix = invSourceWeight * _projectionMatrix.transpose() * ( (_projectionMatrix * invSourceWeight * _pro
   */                                                                                                                        
                                                                                                                             
   /*                                                                                                                        
   // Use the pseudioInverse                                                                                                 
   MatrixXR sourceWeight = MatrixXR::Identity( frequencyCount, frequencyCount );                                             
   MatrixXR targetWeight = MatrixXR::Identity( _halfSize, _halfSize );                                                       
                                                                                                                             
   // Since the source weights is a diagonal matrix, the inverse is the inverse of the diagonal                              
   MatrixXR invSourceWeight = sourceWeight;                                                                                  
   invSourceWeight.diagonal() = invSourceWeight.diagonal().cwise().inverse();                                                
                                                                                                                             
   LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting up the inversion...");                                                         
   // A = W^{-1} K^t [ K W^{-1} K^t + \lambda * I_N ]^{+}                                                                    
   MatrixXR temp = (_projectionMatrix * invSourceWeight * _projectionMatrix.transpose());                                    
   temp.diagonal().cwise() += _regularisation;                                                                               
   MatrixXR pseudioInv;                                                                                                      
   pseudoInverse( temp, &pseudioInv );
   _inverseProjectionMatrix = invSourceWeight * _projectionMatrix.transpose() * pseudioInv;
   */

   // Don't use sourceWeights
   LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting up the inversion...");
   // A = K^t [ K K^t + \lambda * I_N ]^{+}
   MatrixXR temp = (_projectionMatrix * _projectionMatrix.transpose());
   temp.diagonal().array() += _regularisation;
   MatrixXR pseudioInv;
   pseudoInverse( temp, &pseudioInv );
   _inverseProjectionMatrix = _projectionMatrix.transpose() * pseudioInv;

  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting up the peak detector and interpolator...");
  _peak.setup();
  _peakInterp.setup();

  reset();

  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Finished setup.");
}


void PitchInverseProblem::harmonicWeight(MatrixXR f, Real fMin, Real fMax, int harmonicIndex, MatrixXR* result){
  (*result) = MatrixXR::Constant(f.rows(), f.cols(), (fMax + _alpha) / ((harmonicIndex * fMin) + _beta));
}

void PitchInverseProblem::harmonicPosition(MatrixXR f, Real /*fMin*/, Real /*fMax*/, int harmonicIndex, MatrixXR* result){
  (*result) = (harmonicIndex * f * sqrt(1.0 + (pow(harmonicIndex, 2.0) - 1.0) * _inharmonicity)) * (Real)_fftSize / ((Real)_sampleRate);
}

void PitchInverseProblem::harmonicSpread(MatrixXR f, Real /*fMin*/, Real /*fMax*/, int /*harmonicIndex*/, MatrixXR* result){
  (*result) = MatrixXR::Constant(f.rows(), f.cols(), _peakWidth);
}

Real PitchInverseProblem::harmonicWeight(Real f, Real /*fMin*/, Real /*fMax*/, int harmonicIndex){
  //return ((_sampleRate * fMin) + _alpha) / ((harmonicIndex * _sampleRate * fMax) + _beta);
  return ((harmonicIndex * f) + _beta) / ((harmonicIndex * f) + _alpha);
  //return _sampleRate / f / harmonicIndex;
  //return 1.0;
}

Real PitchInverseProblem::harmonicPosition(Real f, Real /*fMin*/, Real /*fMax*/, int harmonicIndex){
  return (harmonicIndex * f * sqrt(1.0 + (pow(harmonicIndex, 2.0) - 1.0) * _inharmonicity)) * (Real)_fftSize / ((Real)_sampleRate);
}

Real PitchInverseProblem::harmonicSpread(Real /*f*/, Real /*fMin*/, Real /*fMax*/, int /*harmonicIndex*/){
  // TODO: change this by a spread function which might or might not change with the position
  //       or other things such as the chirp rate or inharmonicity error
  return _peakWidth;
}


void PitchInverseProblem::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies, MatrixXR* freqs){
  const int rows = spectrum.rows();

  (*pitches).resize( rows, _pitchCount );
  (*saliencies).resize( rows, _pitchCount );
  (*freqs).resize( rows, _projectionMatrix.cols() );

  (*pitches).setZero();
  (*saliencies).setZero();

  for ( int row = 0; row < rows; row++ ) {
    LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Matrix multiplication");
    LOUDIA_DEBUG("PITCHINVERSEPROBLEM: _inverseProjectionMatrix: " << _inverseProjectionMatrix.rows() << ", " << _inverseProjectionMatrix.cols() );
    LOUDIA_DEBUG("PITCHINVERSEPROBLEM: spectrum: " << spectrum.rows() << ", " << spectrum.cols() );
    (*freqs).row( row ) = _inverseProjectionMatrix * spectrum.row( row ).segment(_lowBin, _range).transpose();
  }

  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Find peaks");
  
  _peak.process((*freqs),
                &_starts, pitches, &_ends, saliencies);

  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Interpolate peaks");
  
  MatrixXR salienciesInterp;
  _peakInterp.process((*freqs), (*pitches), (*saliencies),
                      pitches, &salienciesInterp);

  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Setting the pitches");
  
  (*pitches) = (((_highFrequency - _lowFrequency) / (_frequencyCandidateCount-1)) * (*pitches)).array() + _lowFrequency;
  
}

void PitchInverseProblem::reset(){
  // Initial values
  LOUDIA_DEBUG("PITCHINVERSEPROBLEM: Resetting...");
  _peak.reset();
  _peakInterp.reset();
  
}

void PitchInverseProblem::projectionMatrix(MatrixXR* matrix) const {
  (*matrix) = _projectionMatrix;
}

Real PitchInverseProblem::lowFrequency() const{
  return _lowFrequency;
}
  
void PitchInverseProblem::setLowFrequency( Real frequency, bool callSetup ){
  _lowFrequency = frequency;
  if ( callSetup ) setup();
}

Real PitchInverseProblem::highFrequency() const{
  return _highFrequency;
}
  
void PitchInverseProblem::setHighFrequency( Real frequency, bool callSetup ){
  _highFrequency = frequency;
  if ( callSetup ) setup();
}

Real PitchInverseProblem::sampleRate() const{
  return _sampleRate;
}
  
void PitchInverseProblem::setSampleRate( Real frequency, bool callSetup ){
  _sampleRate = frequency;
  if ( callSetup ) setup();
}

int PitchInverseProblem::fftSize() const{
  return _fftSize;
}

void PitchInverseProblem::setFftSize( int size, bool callSetup ) {
  _fftSize = size;
  if ( callSetup ) setup();
}

int PitchInverseProblem::peakWidth() const{
  return _peakWidth;
}

void PitchInverseProblem::setPeakWidth( int width, bool callSetup ) {
  _peakWidth = width;
  if ( callSetup ) setup();
}

int PitchInverseProblem::frequencyCandidateCount() const{
  return _frequencyCandidateCount;
}

void PitchInverseProblem::setFrequencyCandidateCount( int count, bool callSetup ) {
  _frequencyCandidateCount = count;
  if ( callSetup ) setup();
}

int PitchInverseProblem::pitchCount() const{
  return _pitchCount;
}

void PitchInverseProblem::setPitchCount( int count, bool callSetup ) {
  _pitchCount = count;
  if ( callSetup ) setup();
}

int PitchInverseProblem::harmonicCount() const{
  return _harmonicCount;
}

void PitchInverseProblem::setHarmonicCount( int count, bool callSetup ) {
  _harmonicCount = count;
  if ( callSetup ) setu
}

PitchSaliency::PitchSaliency(int fftSize, Real f0, Real f1, Real sampleRate, Real fPrec, int numHarmonics) :
  _fftSize( fftSize ),
  _halfSize( ( _fftSize / 2 ) + 1 ),
  _f0( f0 ),
  _f1( f1 ),
  _fPrec( fPrec ),
  _numHarmonics( numHarmonics ),
  _sampleRate( sampleRate )
{
  LOUDIA_DEBUG("PITCHSALIENCY: Construction fftSize: " << _fftSize
        << " sampleRate: " << _sampleRate
        << " f0: " << _f0
        << " f1: " << _f1
        << " fPrec: " << _fPrec
        << " numHarmonics: " << _numHarmonics );

  setup();
}

PitchSaliency::~PitchSaliency(){}

void PitchSaliency::setup(){
  LOUDIA_DEBUG("PITCHSALIENCY: Setting up...");
  
  _halfSize = ( _fftSize / 2 ) + 1;

  _tMax = _sampleRate / _f0;
  _tMin = _sampleRate / _f1;
  
  _tPrec = _fPrec;


  // Params taken from Klapuri ISMIR 2006
  _alpha = 27; // 27 Hz
  _beta = 320; // 320 Hz
  
  reset();

  LOUDIA_DEBUG("PITCHSALIENCY: Finished setup.");
}

Real PitchSaliency::harmonicWeight(Real /*period*/, Real tLow, Real tUp, int harmonicIndex){
  return ((_sampleRate / tLow) + _alpha) / ((harmonicIndex * _sampleRate / tUp) + _beta);
}

Real PitchSaliency::saliency(Real period, Real deltaPeriod, Real tLow, Real tUp, const MatrixXR& spectrum){
  const int cols = spectrum.cols();
  Real sum = 0.0;
  
  for ( int m = 1; m < _numHarmonics; m++ ) {
    
    int begin = (int)round(m * _fftSize / (period + (deltaPeriod / 2.0)));
    int end = min((int)round(m * _fftSize / (period - (deltaPeriod / 2.0))), cols - 1);

    if (begin < end) sum += harmonicWeight(period, tLow, tUp, m) * spectrum.block(0, begin, 1, end - begin).maxCoeff();
  }

  return sum;
}

void PitchSaliency::process(const MatrixXR& spectrum, MatrixXR* pitches, MatrixXR* saliencies){
  const int rows = spectrum.rows();

  (*pitches).resize( rows, 1 );
  (*saliencies).resize( rows, 1 );
  
  for ( int row = 0; row < rows; row++ ) {

    Real tLow = _tMin;
    Real tUp = _tMax;
    Real sal;

    Real tLowBest = tLow;
    Real tUpBest = tUp;
    Real salBest;
  
    Real period;
    Real deltaPeriod;

    while ( ( tUp - tLow ) > _tPrec ) {
      // Split the best block and compute new limits
      tLow = (tLowBest + tUpBest) / 2.0;
      tUp = tUpBest;
      tUpBest = tLow;
      
      // Compute new saliences for the new blocks
      period = (tLowBest + tUpBest) / 2.0;
      deltaPeriod = tUpBest - tLowBest;
      salBest = saliency(period, deltaPeriod, tLowBest, tUpBest, spectrum.row( row ));

      period = (tLow + tUp) / 2.0;
      deltaPeriod = tUp - tLow;
      sal = saliency(period, deltaPeriod, tLow, tUp, spectrum.row( row ));

      if (sal > salBest) {
        tLowBest = tLow;
        tUpBest = tUp;
        salBest = sal;
      }
    }

    period = (tLowBest + tUpBest) / 2.0;
    deltaPeriod = tUpBest - tLowBest;

    (*pitches)(row, 0) = _sampleRate / period;
    (*saliencies)(row, 0) = saliency(period, deltaPeriod, tLowBest, tUpBest, spectrum.row( row ));
  }
}

void PitchSaliency::reset(){
  // Initial values

}

}
