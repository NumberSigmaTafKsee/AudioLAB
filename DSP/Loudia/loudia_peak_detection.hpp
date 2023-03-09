#pragma once

namespace Loudia
{
    class PeakCOG {
    protected:
        // Internal parameters
        int _fftLength;
        int _bandwidth;

        // Internal variables
        MatrixXR _spectrumAbs2;
        MatrixXR _spectrumArg;
        MatrixXR _spectrumArgDeriv;
        MatrixXR _spectrumAbs2Deriv;

    public:
        PeakCOG(int fftLength, int bandwidth = 6);

        ~PeakCOG();

        void setup();

        void process(const MatrixXC& fft, const MatrixXR& peakPos, MatrixXR* peakCog);

        void reset();

    };


    /**
    * @class PeakDetection
    *
    * @brief Algorithm to find peaks in a vector of Real values.
    *
    * This class represents an object to find peaks in a vector of Real values.
    * The algorithm finds a maximum number of peaks and returns 
    * the indices of the peaks and the values of the peaks in
    * separate matrices.
    * 
    * The maximum number of peaks can be specified using setPeakCount().
    * 
    * The resulting peak arrays may be sorted by position or by magnitude. This can be
    * specified using setSortMethod().
    *
    * When sorting by position it may be interesting to specify a number of candidates, in order
    * to perform a preselection of the highest valued peaks before sorting.  This can be specified
    * using setCandidateCount().
    *
    * The implementation consists in running a sliding windows along the vector in search of 
    * indices which whose value is the maximum of the window.  The size of the window
    * defines the minimum width of the peak.
    * 
    *
    * @author Ricard Marxer
    *
    * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
    */
    class PeakDetection {
    public:
    /**
        @enum SortMethod
        @brief Specifies the way to sort the peak candidates before returning them.
        @sa sortMethod
    */
    enum SortMethod {
        NONE              = 0 /**< No sorting is performed */,
        BYMAGNITUDE       = 1 /**< Sorts the peak candidates by decreasing order of magnitude */,
        BYPOSITION        = 2 /**< Sorts the peak candidates by increasing order of position */
    };

    protected:
        // Internal parameters
        int _peakCount;
        int _minimumPeakWidth;
        int _candidateCount;
        Real _minimumPeakContrast;
            
        SortMethod _sortMethod;

        // Internal variables
        MatrixXR _magnitudes;

        struct peak{
            Real start;
            Real mid;
            Real end;
            Real mag;
            peak(const peak& other)
                :start(other.start)
                ,mid(other.mid)
                ,end(other.end)
                ,mag(other.mag)
            { }

            peak& operator=(const peak& other) {
                start = other.start;
                mid = other.mid;
                end = other.end;
                mag = other.mag;

                return *this;
            }

            peak(Real start, Real mid, Real end, Real mag)
                :start(start)
                ,mid(mid)
                ,end(end)
                ,mag(mag)
            { }

            // A peak is smaller (first in the list)
            // if it's magnitude is larger
            bool operator <(peak const& other) const {
                return mag > other.mag;
            }
        };

        struct byMagnitude{
            bool operator() (const peak& i, const peak& j) const { return ( i.mag > j.mag ); }
        } byMagnitude;

        struct byPosition{
            bool operator() (const peak& i, const peak& j) const { return ( i.mid < j.mid ); }
        } byPosition;

        void findPeaksOld(const VectorXR& a, std::vector<peak>& peaks) const;
        void findPeaks(const VectorXR& a, std::vector<peak>& peaks) const;

    public:
        /**
            Constructs a peak detection object with the given @a peakCount, @a sort method, @a minimumPeakWidth, @a candidateCount and @a minimumPeakContrast parameters given.
        */
        PeakDetection(int peakCount = -1, SortMethod sort = BYMAGNITUDE, int minimumPeakWidth = 3, int candidateCount = -1, Real minimumPeakContrast = 0);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~PeakDetection();

        void reset();
        void setup();
        
        /**
            Detects peaks on each of the rows of @a frames and
            puts the resulting peak indices and magnitudes in the rows of @a peakPositions and 
            @a peakMagnitudes respectively.
            
            @param frames matrix of Real values.
            
            @param peakPositions pointer to a matrix of Real values (but always Integers) for the peak indices.
            The matrix should have the same number of rows as @a frames and peakCount columns. 
            @param peakMagnitudes pointer to a matrix of Real values for the peak magnitudes.
            The matrix should have the same number of rows as @a frames and peakCount columns. 
            Note that if the count of peaks detect is lower than peakCount some values
            of the resulting arrays will be set to -1.0 in order to indicate that it is not
            a peak.
            Note that if the output matrices are not of the required size they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& frames,
                    MatrixXR* peakStarts, MatrixXR* peakPositions, MatrixXR* peakEnds,
                    MatrixXR* peakMagnitudes);
        /**
            Returns the maximum number of peaks to be detected by the algorithm.
            
            By default it is 1024 / 3.
        */
        int peakCount() const;

        /**
            Specifies the maximum peak @a count to be detected by the algorithm.
            If <= 0, then all possible peaks are detected.
        */
        void setPeakCount( int count, bool callSetup = true );

        /**
            Returns the minimum width of a peak for it to be detected.
            
            By default it is 3.
        */
        int minimumPeakWidth() const;

        /**
            Specifies the minimum @a width of a peak for it to be detected.
        */
        void setMinimumPeakWidth( int width, bool callSetup = true );

        /**
            Returns the number of highest value candidates to be considered before sorting.
            Note that if the value is <= 0, then no preselection is performed
            and all detected peaks are considered as candidates.
            By default it is -1.
        */
        int candidateCount() const;

        /**
            Specifies the number of highest value candidates to be considered before sorting.
            Note that if the value is <= 0, then no preselection is performed
            and all detected peaks are considered as candidates.
        */
        void setCandidateCount( int count, bool callSetup = true );

        /**
            Returns the minimum contrast of a peak for it to be detected.
            
            The contrast is considered of a peak is the maximum value minus the minimum value
            of all the points in the peak detection running window.
            
            By default it is 0.0.
        */
        Real minimumPeakContrast() const;

        /**
            Specifies the minimum contrast of a peak for it to be detected.
            
            The contrast is considered of a peak is the maximum value minus the minimum value
            of all the points in the peak detection running window.
        */
        void setMinimumPeakContrast( Real contrast, bool callSetup = true );

        /**
            Returns the method for sorting the peaks.
            
            By default it is BYMAGNITUDE.
        */
        SortMethod sortMethod() const;

        /**
            Specifies the method for sorting the peaks.
        */
        void setSortMethod( SortMethod method, bool callSetup = true );
    };


    /**
    * @class PeakDetectionComplex
    *
    * @brief Algorithm to find peaks in a vector of Complex values.
    *
    * This class represents an object to find peaks in a vector of Complex values.
    * The algorithm finds a maximum number of peaks and returns 
    * the indices of the peaks and the values of the peaks in
    * separate matrices.
    * 
    * The maximum number of peaks can be specified using setPeakCount().
    * 
    * The resulting peak arrays may be sorted by position or by magnitude. This can be
    * specified using setSortMethod().
    *
    * When sorting by position it may be interesting to specify a number of candidates, in order
    * to perform a preselection of the highest valued peaks before sorting.  This can be specified
    * using setCandidateCount
    *
    * The implementation consists in running a sliding windows along the vector in search of 
    * indices which whose value is the maximum of the window.  The size of the window
    * defines the minimum width of the peak.
    * 
    *
    * @author Ricard Marxer
    *
    * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
    */
    class PeakDetectionComplex {
    public:
    /**
        @enum SortMethod
        @brief Specifies the way to sort the peak candidates before returning them.
        @sa sortMethod
    */
    enum SortMethod {
        NONE              = 0 /**< No sorting is performed */,
        BYMAGNITUDE       = 1 /**< Sorts the peak candidates by decreasing order of magnitude */,
        BYPOSITION        = 2 /**< Sorts the peak candidates by increasing order of position */
    };

    protected:
        // Internal parameters
        int _peakCount;
        int _minimumPeakWidth;
        int _candidateCount;
        Real _minimumPeakContrast;
            
        SortMethod _sortMethod;

        // Internal variables
        MatrixXR _magnitudes;
        MatrixXR _phases;

    public:
        /**
            Constructs a peak detection object with the given @a peakCount, @a sort method, @a minimumPeakWidth, @a candidateCount and @a minimumPeakContrast parameters given.
        */
        PeakDetectionComplex(int peakCount = 1024 / 3, SortMethod sort = BYMAGNITUDE, int minimumPeakWidth = 3, int candidateCount = -1, Real minimumPeakContrast = 0);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~PeakDetectionComplex();

        void reset();
        void setup();
        
        /**
            Detects peaks on each of the rows of @a frames and
            puts the resulting peak indices and magnitudes in the rows of @a peakPositions and 
            @a peakMagnitudes respectively.
            
            @param frames matrix of Complex values.
            
            @param peakPositions pointer to a matrix of Real values (but always Integers) for the peak indices.
            The matrix should have the same number of rows as @a frames and peakCount columns. 
            @param peakMagnitudes pointer to a matrix of Real values for the peak magnitudes.
            The matrix should have the same number of rows as @a frames and peakCount columns. 
            @param peakPhases pointer to a matrix of Real values for the peak phases.
            The matrix should have the same number of rows as @a frames and peakCount columns. 
            Note that if the count of peaks detect is lower than peakCount some values
            of the resulting arrays will be set to -1.0 in order to indicate that it is not
            a peak.
            Note that if the output matrices are not of the required size they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXC& frames,
                    MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases);

        /**
            Returns the maximum number of peaks to be detected by the algorithm.
            
            By default it is 1024 / 3.
        */
        int peakCount() const;

        /**
            Specifies the maximum peak @a count to be detected by the algorithm.
            If <= 0, then all possible peaks are detected.
        */
        void setPeakCount( int count, bool callSetup = true );

        /**
            Returns the minimum width of a peak for it to be detected.
            
            By default it is 3.
        */
        int minimumPeakWidth() const;

        /**
            Specifies the minimum @a width of a peak for it to be detected.
        */
        void setMinimumPeakWidth( int width, bool callSetup = true );

        /**
            Returns the number of highest value candidates to be considered before sorting.
            Note that if the value is <= 0, then no preselection is performed
            and all detected peaks are considered as candidates.
            By default it is -1.
        */
        int candidateCount() const;

        /**
            Specifies the number of highest value candidates to be considered before sorting.
            Note that if the value is <= 0, then no preselection is performed
            and all detected peaks are considered as candidates.
        */
        void setCandidateCount( int count, bool callSetup = true );

        /**
            Returns the minimum contrast of a peak for it to be detected.
            
            The contrast is considered of a peak is the maximum value minus the minimum value
            of all the points in the peak detection running window.
            
            By default it is 0.0.
        */
        int minimumPeakContrast() const;

        /**
            Specifies the minimum contrast of a peak for it to be detected.
            
            The contrast is considered of a peak is the maximum value minus the minimum value
            of all the points in the peak detection running window.
        */
        void setMinimumPeakContrast( Real contrast, bool callSetup = true );

        /**
            Returns the method for sorting the peaks.
            
            By default it is BYMAGNITUDE.
        */
        SortMethod sortMethod() const;

        /**
            Specifies the method for sorting the peaks.
        */
        void setSortMethod( SortMethod method, bool callSetup = true );
    };

PeakDetection::PeakDetection(int peakCount, SortMethod sortMethod, int minimumPeakWidth, int candidateCount, Real minimumPeakContrast)
{
    LOUDIA_DEBUG("PEAKDETECTION: Constructor peakCount: " << peakCount
                 << ", minimumPeakWidth: " << minimumPeakWidth
                 << ", candidateCount: " << candidateCount);

    setPeakCount( peakCount, false );
    setMinimumPeakWidth( minimumPeakWidth, false );
    setCandidateCount( candidateCount, false );
    setMinimumPeakContrast( minimumPeakContrast, false );
    setSortMethod( sortMethod, false );

    setup();

    LOUDIA_DEBUG("PEAKDETECTION: Constructed");
}

PeakDetection::~PeakDetection() {
    // TODO: Here we should free the buffers
    // but I don't know how to do that with MatrixXR and MatrixXR
    // I'm sure Nico will...
}


void PeakDetection::setup(){
    // Prepare the buffers
    LOUDIA_DEBUG("PEAKDETECTION: Setting up...");
    
    reset();

    LOUDIA_DEBUG("PEAKDETECTION: Finished set up...");
}

void PeakDetection::findPeaksOld(const VectorXR& a, vector<peak>& peaks) const {
    int maxRow;
    int maxCol;

    Real maxVal;
    Real minVal;

    const int cols = a.size();

    const int halfPeakWidth = _minimumPeakWidth / 2;
    for ( int j = halfPeakWidth; j < cols - halfPeakWidth; j++) {
        // If we don't need sorting then only the first peakCount peaks are needed
        if( ( _sortMethod == NONE ) && ( (int)peaks.size() > _peakCount ) ) break;

        int inf = j - halfPeakWidth;

        // Get the maximum value and position of a region (corresponding to the min bandwidth of the peak)
        // of the spectrum
        maxVal = a.segment(inf, _minimumPeakWidth).maxCoeff( &maxRow, &maxCol );

        // If the position of the maximum value is the center, then consider it as a peak candidate
        if ( maxCol == halfPeakWidth ) {

            // Get the mininum value of the region
            minVal = a.segment(inf, _minimumPeakWidth).minCoeff();

            // If the contrast is bigger than what minimumPeakContrast says, then select as peak
            if ( maxVal - minVal >= _minimumPeakContrast ) {

                peaks.push_back( peak(-1, j, -1, a(j)) );

            }
        }
    }
}

void PeakDetection::findPeaks(const VectorXR& a, vector<peak>& peaks) const {
    Real minIn = 0;
    Real minVal = a[0];

    Real maxIn = 0;
    Real maxVal = a[0];

    Real start = 0;

    Real tol = 0;

    const int INCREASING = 0;
    const int DECREASING = 1;

    int state = INCREASING;

    int i = 1;
    while (i < a.size()) {
        if (state == INCREASING) {
            Real diffMax = a[i] - maxVal;
            if (diffMax > 0) {
                maxIn = i;
                maxVal = a[i];
            } else if (-diffMax > tol) {
                state = DECREASING;
                minIn = i;
                minVal = a[i];
            }

        } else if (state == DECREASING) {
            Real diffMin = a[i] - minVal;

            if (diffMin < 0) {
                minIn = i;
                minVal = a[i];

            } else if (diffMin > tol) {
                state = INCREASING;
                peaks.push_back(peak(start, maxIn, i-1, maxVal));

                // If we don't need sorting then only the first peakCount peaks are needed
                if( ( _sortMethod == NONE ) && ( (int)peaks.size() >= _peakCount ) ) return;

                start = i-1;
                maxIn = i;
                maxVal = a[i];
            }
        }
        i += 1;
    }

    peaks.push_back(peak(start, maxIn, i-1, maxVal));
}

void PeakDetection::process(const MatrixXR& frames,
                            MatrixXR* peakStarts, MatrixXR* peakPositions, MatrixXR* peakEnds,
                            MatrixXR* peakMagnitudes) {
    LOUDIA_DEBUG("PEAKDETECTION: Processing");

    const int rows = frames.rows();
    const int cols = frames.cols();

    if (_peakCount < 0) {
      _peakCount = cols / 2; 
    }
    
    LOUDIA_DEBUG("PEAKDETECTION: Processing, frames.shape: (" << rows << ", " << cols << ")");

    (*peakStarts).resize(rows, _peakCount);
    (*peakStarts).setConstant(-1);
    (*peakPositions).resize(rows, _peakCount);
    (*peakPositions).setConstant(-1);
    (*peakEnds).resize(rows, _peakCount);
    (*peakEnds).setConstant(-1);

    (*peakMagnitudes).resize(rows, _peakCount);
    (*peakMagnitudes).setConstant(-1);

    _magnitudes = frames.array().abs();

    LOUDIA_DEBUG("PEAKDETECTION: Processing, _magnitudes.shape: (" << rows << ", " << cols << ")");

    vector<peak> peaks;
    peaks.reserve( cols );

    for ( int i = 0 ; i < rows; i++){

        peaks.clear();
        findPeaks(_magnitudes.row(i), peaks);

        // Get the largest candidates
        int candidateCount = (int)peaks.size();
        if( _candidateCount > 0 ) {
            candidateCount = min(candidateCount, _candidateCount);
            std::sort(peaks.begin(), peaks.end(), byMagnitude);
        }

        // Sort the candidates using position or magnitude
        switch ( _sortMethod ) {
        case BYPOSITION:
            std::sort(peaks.begin(), peaks.begin() + candidateCount, byPosition);
            break;

        case BYMAGNITUDE:
            // We have not done a candidate preselection, we must do the sorting
            if (_candidateCount <= 0)
                std::sort(peaks.begin(), peaks.begin() + candidateCount, byMagnitude);
            break;

        case NONE:
        default:
            break;
        }

        // Take the first peakCount
        int peakCount = min(_peakCount, candidateCount);
        // Put the peaks in the matrices
        for( int j = 0; j < peakCount; j++ ){
            (*peakStarts)(i, j) = peaks[j].start;
            (*peakEnds)(i, j) = peaks[j].end;
            (*peakMagnitudes)(i, j) = peaks[j].mag;
            (*peakPositions)(i, j) = peaks[j].mid;
        }
    }

    LOUDIA_DEBUG("PEAKDETECTION: Finished Processing");
}

void PeakDetection::reset(){
    // Initial values
}

int PeakDetection::peakCount() const {
    return _peakCount;
}

void PeakDetection::setPeakCount( int count, bool callSetup ) {
    _peakCount = count;
    if ( callSetup ) setup();
}

int PeakDetection::candidateCount() const {
    return _candidateCount;
}

void PeakDetection::setCandidateCount( int count, bool callSetup ) {
    _candidateCount = count;
    if ( callSetup ) setup();
}

int PeakDetection::minimumPeakWidth() const {
    return _minimumPeakWidth;
}

void PeakDetection::setMinimumPeakWidth( int width, bool callSetup ) {
    _minimumPeakWidth = width;
    if ( callSetup ) setup();
}

Real PeakDetection::minimumPeakContrast() const {
    return _minimumPeakContrast;
}

void PeakDetection::setMinimumPeakContrast( Real contrast, bool callSetup ) {
    _minimumPeakContrast = contrast;
    if ( callSetup ) setup();
}

PeakDetection::SortMethod PeakDetection::sortMethod() const {
    return _sortMethod;
}

void PeakDetection::setSortMethod( SortMethod method, bool callSetup ) {
    _sortMethod = method;
    if ( callSetup ) setup();
}

PeakCOG::PeakCOG(int fftLength, int bandwidth) :
  _fftLength(fftLength),
  _bandwidth(bandwidth)
{

  LOUDIA_DEBUG("PEAKCOG: Constructor fftLength: " << _fftLength);

  setup();
}

PeakCOG::~PeakCOG() {}


void PeakCOG::setup() {
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKCOG: Setting up...");

  reset();

  LOUDIA_DEBUG("PEAKCOG: Finished set up...");
}


void PeakCOG::process(const MatrixXC& fft, const MatrixXR& peakPos, MatrixXR* peakCog) {
  LOUDIA_DEBUG("PEAKCOG: Processing windowed");
  const int rows = fft.rows();
  const int cols = fft.cols();
  const int halfCols = min((int)ceil(_fftLength / 2.0), cols);
  const int peakCount = peakPos.cols();

  LOUDIA_DEBUG("PEAKCOG: fft.shape " << fft.rows() << "," << fft.cols());
  _spectrumAbs2 = fft.block(0, 0, rows, halfCols).array().abs2();
  LOUDIA_DEBUG("PEAKCOG: Spectrum resized rows: " << rows << " halfCols: " << halfCols);

  unwrap(fft.block(0, 0, rows, halfCols).array().angle(), &_spectrumArg);
  derivate(_spectrumArg, &_spectrumArgDeriv);

  derivate(_spectrumAbs2, &_spectrumAbs2Deriv);

  (*peakCog).resize(rows, peakCount);
  (*peakCog).setZero();

  for (int row = 0; row < rows; row++) {
    for (int i = 0; i < peakCount; i++) {
      if (peakPos(row, i) == -1) {
        continue;
      }

      // Find the start and end of

      int start = peakPos(row, i);
      for (; start > 0; start-- ) {
        if (_spectrumAbs2Deriv(row, start) * _spectrumAbs2Deriv(row, start-1) < 0) {
          break;
        }
      }

      int end = peakPos(row, i);
      for (; end < _spectrumAbs2Deriv.cols()-1; end++ ) {
        if (_spectrumAbs2Deriv(row, end) * _spectrumAbs2Deriv(row, end+1) < 0) {
          break;
        }
      }

      //LOUDIA_DEBUG("peakWidth:" << end-start);

      // Calculate the actual center of gravity of the peak
      if ( (end - start) >= 3) {
        (*peakCog)(row, i) = ((-_spectrumArgDeriv).block(row, start, 1, end-start).array() * _spectrumAbs2.block(row, start, 1, end-start).array()).sum() / _spectrumAbs2.block(row, start, 1, end-start).sum();
      }

    }
  }

  LOUDIA_DEBUG("PEAKCOG: Finished Processing");
}

void PeakCOG::reset() {
  // Initial values
}

struct peakComplex{
  Real pos;
  Real mag;
  Real phase;
  peakComplex(const peakComplex& other)
    :pos(other.pos), mag(other.mag), phase(other.phase) { }

  peakComplex& operator=(const peakComplex& other) {
    pos = other.pos;
    mag = other.mag;
    phase = other.phase;

    return *this;
  }

  peakComplex(Real pos, Real mag, Real phase)
    :pos(pos), mag(mag), phase(phase) { }
  
  // A peak is smaller (first in the list)
  // if it's magnitude is larger
  bool operator <(peakComplex const& other) const {
    return mag > other.mag;
  }
};

struct byMagnitudeComp{
  bool operator() (const peakComplex& i, const peakComplex& j) const { return ( i.mag > j.mag ); }
} byMagnitudeComplex;

struct byPositionComp{
  bool operator() (const peakComplex& i, const peakComplex& j) const { return ( i.pos < j.pos ); }
} byPositionComplex;

PeakDetectionComplex::PeakDetectionComplex(int peakCount, SortMethod sortMethod, int minimumPeakWidth, int candidateCount, Real minimumPeakContrast)
{
  LOUDIA_DEBUG("PEAKDETECTION: Constructor peakCount: " << peakCount 
        << ", minimumPeakWidth: " << minimumPeakWidth
        << ", candidateCount: " << candidateCount);
  
  setPeakCount( peakCount, false );
  setMinimumPeakWidth( minimumPeakWidth, false );
  setCandidateCount( candidateCount, false );
  setMinimumPeakContrast( minimumPeakContrast, false );
  setSortMethod( sortMethod, false );
  
  setup();

  LOUDIA_DEBUG("PEAKDETECTIONCOMPLEX: Constructed");
}

PeakDetectionComplex::~PeakDetectionComplex() {
}


void PeakDetectionComplex::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKDETECTIONCOMPLEX: Setting up...");

  reset();

  LOUDIA_DEBUG("PEAKDETECTIONCOMPLEX: Finished set up...");
}


void PeakDetectionComplex::process(const MatrixXC& frames, 
                                   MatrixXR* peakPositions, MatrixXR* peakMagnitudes, MatrixXR* peakPhases){
  LOUDIA_DEBUG("PEAKDETECTIONCOMPLEX: Processing");
  
  const int rows = frames.rows();
  const int cols = frames.cols();
  
  LOUDIA_DEBUG("PEAKDETECTIONCOMPLEX: Processing, frames.shape: (" << rows << ", " << cols << ")");

  (*peakPositions).resize(rows, _peakCount);
  (*peakPositions).setConstant(-1);

  (*peakMagnitudes).resize(rows, _peakCount);
  (*peakMagnitudes).setConstant(-1);

  (*peakPhases).resize(rows, _peakCount);
  (*peakPhases).setConstant(-1);

  _magnitudes = frames.array().abs();
  _phases = frames.array().angle();

  int maxRow;
  int maxCol;
  
  Real maxVal;
  Real minVal;

  const int _halfPeakWidth = _minimumPeakWidth / 2;
  
  vector<peakComplex> peakVector;
  peakVector.reserve( cols );
  int detectedCount;
  
  for ( int i = 0 ; i < rows; i++){

    peakVector.clear();
    detectedCount = 0;
    
    for ( int j = _halfPeakWidth; j < cols - _halfPeakWidth; j++) {
      // If we don't need sorting then only the first peakCount peakVector are needed
      if ( ( _sortMethod == NONE ) && ( detectedCount >= _peakCount ) ) break;

      int inf = j - _halfPeakWidth;
      
      // Get the maximum value and position of a region (corresponding to the min bandwidth of the peak)
      // of the spectrum
      maxVal = _magnitudes.row(i).segment(inf, _minimumPeakWidth).maxCoeff( &maxRow, &maxCol );
      
      // If the position of the maximum value is the center, then consider it as a peak candidate
      if ( maxCol == _halfPeakWidth ) {

        // Get the mininum value of the region
        minVal = _magnitudes.row(i).segment(inf, _minimumPeakWidth).minCoeff();

        // If the contrast is bigger than what minimumPeakContrast says, then select as peak
        if ( (maxVal - minVal) >= _minimumPeakContrast ) {

          peakVector.push_back( peakComplex(j, _magnitudes(i, j), _phases(i, j)) );
          detectedCount ++;

        }
      }
    }
    
    // Get the largest candidates
    int candidateCount = detectedCount;
    if( _candidateCount > 0 ) {
      candidateCount = min( candidateCount, _candidateCount );
      partial_sort( peakVector.begin(), peakVector.begin() + candidateCount, peakVector.end() , byMagnitudeComplex );
    }

    // Sort and take the first peakCount peakVector
    int peakCount = min( _peakCount, candidateCount );

    // Sort the candidates using position or magnitude
    switch ( _sortMethod ) {
    case BYPOSITION:      
      partial_sort( peakVector.begin(), peakVector.begin() + peakCount, peakVector.begin() + candidateCount, byPositionComplex );
      break;
      
    case BYMAGNITUDE:
      // We have not done a candidate preselection, we must do the sorting
      if ( _candidateCount <= 0 )
        partial_sort( peakVector.begin(), peakVector.begin() + peakCount, peakVector.begin() + candidateCount, byMagnitudeComplex );
      break;
      
    case NONE:
    default:
      break;
    }
    
    // Put the peaks in the matrices
    for( int j = 0; j < peakCount; j++ ){
      (*peakMagnitudes)(i, j) = peakVector[j].mag;
      (*peakPositions)(i, j) = peakVector[j].pos;
      (*peakPhases)(i, j) = peakVector[j].phase;
    }
  }

  LOUDIA_DEBUG("PEAKDETECTIONCOMPLEX: Finished Processing");
}

void PeakDetectionComplex::reset(){
  // Initial values
}

int PeakDetectionComplex::peakCount() const {
  return _peakCount;
}

void PeakDetectionComplex::setPeakCount( int count, bool callSetup ) {
  _peakCount = count;
  if ( callSetup ) setup();  
}

int PeakDetectionComplex::candidateCount() const {
  return _candidateCount;
}

void PeakDetectionComplex::setCandidateCount( int count, bool callSetup ) {
  _candidateCount = count;
  if ( callSetup ) setup();  
}

int PeakDetectionComplex::minimumPeakWidth() const {
  return _minimumPeakWidth;
}

void PeakDetectionComplex::setMinimumPeakWidth( int width, bool callSetup ) {
  _minimumPeakWidth = width;
  if ( callSetup ) setup();
}

int PeakDetectionComplex::minimumPeakContrast() const {
  return _minimumPeakContrast;
}

void PeakDetectionComplex::setMinimumPeakContrast( Real contrast, bool callSetup ) {
  _minimumPeakContrast = contrast;
  if ( callSetup ) setup();
}

PeakDetectionComplex::SortMethod PeakDetectionComplex::sortMethod() const {
  return _sortMethod;
}

void PeakDetectionComplex::setSortMethod( SortMethod method, bool callSetup ) {
  _sortMethod = method;
  if ( callSetup ) setup();  
}
}