#pragma once

namespace Loudia
{
    class PeakSynthesize {
    protected:
        // Internal parameters
        int _windowSize;
        Window::WindowType _windowType;

        int _fftSize;

        // Internal variables

    public:
        PeakSynthesize(int windowSize, int fftSize, Window::WindowType windowType = Window::RECTANGULAR);

        ~PeakSynthesize();

        void setup();

        void process(const MatrixXR& trajPositions, const MatrixXR& trajMagnitudes,
                    MatrixXR* spectrum);

        void reset();

    };

    /**
    * @class PeakTracking
    *
    * @brief Algorithm to find peak trajectories in vectors of Complex values representing FFT, given peak
    * positions and peak magnitudes.
    *
    * The algorithm finds a maximum number of peak trajectories and returns 
    * the positions of the trajectory positions (in fractional index units) and the trajectory magnitudes (in decibel units) in
    * separate matrices.
    * 
    * The maximum number of trajectories can be specified using setTrajectoryCount().
    * 
    * The algorithm operates by matching the peaks in the current frames to the existing trajectories.
    * During the matching process a maximum frequency change of a peak can be specified using setMaximumFrequencyChange().
    *
    * The matching process also requires a trajectory to stay unmatched during a given number of frames for the trajectory to
    * disappear and leave a slot for another trajectory to be found.  The number fo silent frames can be specified using
    * silentFrameCount().
    *
    * @author Ricard Marxer
    *
    * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
    */
    class PeakTracking {
    protected:
        // Internal parameters
        int _trajectoryCount;
        Real _maximumFrequencyChange;
        int _silentFrameCount;

        // Internal variables
        MatrixXR _trajPositions, _trajMagnitudes;
        MatrixXR _pastTrajPositions, _pastTrajMagnitudes;
        
        bool createTrajectory(Real peakPos, Real peakMag,
                                MatrixXR* pastTrajPositions, MatrixXR* pastTrajMagnitudes,
                                MatrixXR* trajPositions, MatrixXR* trajMagnitudes,
                                int row);
        

    public:
        /**
            Constructs a peak tracking object with the given @a trajectoryCount, @a maximumFrequencyChange and @a silentFrameCount settings.
        */
        PeakTracking(int trajectoryCount = 20, Real maximumFrequencyChange = 3.0, int silentFrameCount = 3);

        /**
            Destroys the algorithm and frees its resources.
        */
        ~PeakTracking();

        void setup();
        void reset();

        /**
            Tracks peaks on each of the rows of @a ffts, @a peakPositions, @a peakMagnitudes and
            puts the resulting trajectory indices and magnitudes in the rows of @a trajectoryPositions and 
            @a trajectoryMagnitudes respectively.
            
            @param ffts matrix of Complex values representing the FFT frames.
            
            @param peakPositions matrix of Real values for the peak positions (in fractional index units).
            The matrix should have the same number of rows as @a ffts. 
            @param peakMagnitudes matrix of Real values for the peak magnitudes (in decibel units).
            The matrix should have the same number of rows as @a ffts. 
            @param trajectoryPositions pointer to a matrix of Real values for the trajectory positions (in fractional index units).
            The matrix should have the same number of rows as @a ffts and trajectoryCount columns. 
            @param trajectoryMagnitudes pointer to a matrix of Real values for the trajectory magnitudes (in decibel units).
            The matrix should have the same number of rows as @a ffts and trajectoryCount columns. 
            Note that if the count of trajectories detected is lower than trajectoryCount some values
            of the resulting arrays will be set to -1.0 in order to indicate that it is not
            a trajectory.
            Note that if the output matrices are not of the required size they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXC& ffts, 
                    const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
                    MatrixXR* trajectoryPositions, MatrixXR* trajectoryMagnitudes);
        
        /**
            Returns the maximum number of trajectories to be detected by the algorithm.
            
            By default it is 20.
        */
        int trajectoryCount() const;

        /**
            Specifies the maximum trajectory @a count to be detected by the algorithm.
            If <= 0, then all possible trajectories are detected.
            @param count the maximum number of trajectories to be tracked
            @param callSetup a flag specifying whether the setup() method must be called after setting the parameter.
        */
        void setTrajectoryCount( int count, bool callSetup = true );
        
        /**
            Returns the maximum frequency change of a peak for it to be matched to an existing trajectory.
            
            The change is specified in fractional index units.
            
            By default it is 3.0.
        */
        Real maximumFrequencyChange() const;

        /**
            Specifies the maximum frequency change of a peak for it to be matched to an existing trajectory.
            
            The change is specified in fractional index units.
            @param change the maximum changed allowed between a peak and an existing trajectory
            @param callSetup a flag specifying whether the setup() method must be called after setting the parameter.
        */
        void setMaximumFrequencyChange( Real change, bool callSetup = true );

        /**
            Returns the count of frames a trajectory must stay unmatched for it
            to disappear and leave the slot for another possible trajectory.
            
            By default it is 3.
        */
        int silentFrameCount() const;

        /**
            Specifies the @a count of frames a trajectory must stay unmatched for it
            to disappear and leave the slot for another possible trajectory.
            
            @param count the number of silent frames
            @param callSetup a flag specifying whether the setup() method must be called after setting the parameter.
        */
        void setSilentFrameCount( int count, bool callSetup = true );

    };


    
PeakSynthesize::PeakSynthesize(int windowSize, int fftSize, Window::WindowType windowType) :
  _windowSize( windowSize ),
  _windowType( windowType ),
  _fftSize( fftSize )

{
  LOUDIA_DEBUG("PEAKSYNTHESIZE: Constructor");
  
  setup();
  
  LOUDIA_DEBUG("PEAKSYNTHESIZE: Constructed");
}

PeakSynthesize::~PeakSynthesize() {}


void PeakSynthesize::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKSYNTHESIZE: Setting up...");

  reset();

  LOUDIA_DEBUG("PEAKSYNTHESIZE: Finished set up...");
}


void PeakSynthesize::process(const MatrixXR& trajPositions, const MatrixXR& trajMagnitudes,
                             MatrixXR* spectrum){
  
  LOUDIA_DEBUG("PEAKSYNTHESIZE: Processing");
  
  spectrum->resize(trajPositions.rows(), (int)ceil(_fftSize/2.0));
  spectrum->setZero();
  
  MatrixXR trajMags;
  dbToMag(trajMagnitudes, &trajMags);
  
  for ( int row = 0 ; row < spectrum->rows(); row++ ) {
  
    for ( int i = 0; i < trajPositions.cols(); i++ ) {
      
      // If the position is -1 do nothing since it means it is nothing
      if( trajPositions(row, i) != -1 ){
        MatrixXR windowTransform;
        int begin, end;

        switch(_windowType){

        case Window::RECTANGULAR:
          // TODO: Implement this window transform

        case Window::HANN:
        case Window::HANNING:
          hannTransform(trajPositions(row, i), trajMags(row, i), 
                        _windowSize, _fftSize, 
                        &windowTransform, &begin, &end);
          break;
          
        case Window::HAMMING:
          hammingTransform(trajPositions(row, i), trajMags(row, i), 
                           _windowSize, _fftSize, 
                           &windowTransform, &begin, &end);
          break;

        default:
          LOUDIA_DEBUG("ERROR: Unknown type of window");
          // Throw ValueError unknown window type
          break;

        }
        
        spectrum->block(0, begin, 1, windowTransform.cols()) += windowTransform.row(0);        
      }
    }
  }
  
  LOUDIA_DEBUG("PEAKSYNTHESIZE: Finished Processing");
}

void PeakSynthesize::reset(){
  // Initial values
}

PeakTracking::PeakTracking(int trajectoryCount, Real maximumFrequencyChange, int silentFrameCount)
{
  LOUDIA_DEBUG("PEAKTRACKING: Constructor");
  
  setTrajectoryCount( trajectoryCount, false );
  setMaximumFrequencyChange( maximumFrequencyChange, false );
  setSilentFrameCount( silentFrameCount, false );

  setup();
  
  LOUDIA_DEBUG("PEAKTRACKING: Constructed");
}

PeakTracking::~PeakTracking() {
}


void PeakTracking::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKTRACKING: Setting up...");
  
  reset();
  
  LOUDIA_DEBUG("PEAKTRACKING: Finished set up...");
}


void PeakTracking::process(const MatrixXC& fft,
                           const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
                           MatrixXR* trajPositions, MatrixXR* trajMagnitudes){
  
  LOUDIA_DEBUG("PEAKTRACKING: Processing");  
  
  (*trajPositions).resize(fft.rows(), _trajectoryCount);
  (*trajMagnitudes).resize(fft.rows(), _trajectoryCount);
  
  (*trajPositions) = MatrixXR::Constant(fft.rows(), _trajectoryCount, -1.0);
  (*trajMagnitudes) = MatrixXR::Constant(fft.rows(), _trajectoryCount, -120.0);

  MatrixXR currPeakPositions = peakPositions;
  MatrixXR currPeakMagnitudes = peakMagnitudes;
  
  for ( int row = 0 ; row < fft.rows(); row++ ) {

    // Find the closest peak to each of the trajectories
    for ( int i = 0 ; i < _pastTrajPositions.cols(); i++  ) {
      
      if( ! isinf( _pastTrajPositions(row, i) ) ) {
        
        int posRow, posCol;
        Real minFreqBinChange = (currPeakPositions.row(row).array() - _pastTrajPositions(row, i)).abs().minCoeff(&posRow, &posCol);
        
        if ( minFreqBinChange <= _maximumFrequencyChange ) {
          // A matching peak has been found
          LOUDIA_DEBUG("PEAKTRACKING: Processing 'Matching peak: " << posCol << "' minFreqBinChange: " << minFreqBinChange);
          
          (*trajPositions)(row, i) = currPeakPositions(row, posCol);
          (*trajMagnitudes)(row, i) = currPeakMagnitudes(row, posCol);
          
          _pastTrajPositions(0, i) = (*trajPositions)(row, i);
          _pastTrajMagnitudes(0, i) = (*trajMagnitudes)(row, i);

          currPeakPositions(row, posCol) = numeric_limits<Real>::infinity();
          currPeakMagnitudes(row, posCol) = numeric_limits<Real>::infinity();
          
        } else {
          // No matching peak has been found
          LOUDIA_DEBUG("PEAKTRACKING: Processing 'No matching peaks' minFreqBinChange: " << minFreqBinChange);
          
          if ( _pastTrajMagnitudes(0, i) <= (-120.0 - _silentFrameCount) ) {

            // The trajectory has been silent too long (resetting it)

            _pastTrajMagnitudes(0, i) = numeric_limits<Real>::infinity();
            _pastTrajPositions(0, i) = numeric_limits<Real>::infinity();

          } else if ( _pastTrajMagnitudes(0, i) <= -120.0 ) {

            // The trajectory has been silent for one more frame

            _pastTrajMagnitudes(0, i) -= 1;

          } else {

            // The first frame the trajectory is silent

            _pastTrajMagnitudes(0, i) = -120.0;
            
          }
          
          (*trajPositions)(row, i) = isinf(_pastTrajPositions(0, i)) ? -1 : _pastTrajPositions(0, i);
          (*trajMagnitudes)(row, i) = -120.0;
          
        }
      }
    }
      
    // Find those peaks that haven't been assigned and create new trajectories
    for ( int i = 0; i < currPeakPositions.cols(); i++ ) {    
      Real pos = currPeakPositions(row, i);
      Real mag = currPeakMagnitudes(row, i);
        
      if( ! isinf( pos ) ){
        bool created = createTrajectory(pos, mag, 
                                        &_pastTrajPositions, &_pastTrajMagnitudes,
                                        trajPositions, trajMagnitudes,
                                        row);
        
        if (! created ){
          LOUDIA_DEBUG("PEAKTRACKING: Processing the trajectory could not be created");
        }
      }  
    }
  }
  
  LOUDIA_DEBUG("PEAKTRACKING: Finished Processing");
}

bool PeakTracking::createTrajectory(Real peakPos, Real peakMag,
                                    MatrixXR* pastTrajPositions, MatrixXR* pastTrajMagnitudes,
                                    MatrixXR* trajPositions, MatrixXR* trajMagnitudes,
                                    int row) {

  int maxRow, maxCol;
  Real maxPos = (*pastTrajPositions).row(row).maxCoeff(&maxRow, &maxCol);

  if ( isinf( maxPos ) ) {
    //DEBUG("MAXCOL: " << maxCol);
    //DEBUG("ROW: " << row);
    
    (*pastTrajPositions)(0, maxCol) = peakPos;
    (*pastTrajMagnitudes)(0, maxCol) = peakMag;

    //DEBUG("Past: ");

    (*trajPositions)(row, maxCol) = peakPos;
    (*trajMagnitudes)(row, maxCol) = peakMag;

    return true;
  }

  return false;
}

void PeakTracking::reset(){
  // Initial values
  if ( !numeric_limits<Real>::has_infinity ) {
    // Throw PlatformError infinity not supported
  }

  Real inf = numeric_limits<Real>::infinity();
  _pastTrajPositions = MatrixXR::Constant(1, _trajectoryCount, inf);
  _pastTrajMagnitudes = MatrixXR::Constant(1, _trajectoryCount, inf);
}

int PeakTracking::trajectoryCount() const {
  return _trajectoryCount;
}

void PeakTracking::setTrajectoryCount( int count, bool callSetup ) {
  _trajectoryCount = count;
  if ( callSetup ) setup();  
}

Real PeakTracking::maximumFrequencyChange() const {
  return _maximumFrequencyChange;
}

void PeakTracking::setMaximumFrequencyChange( Real change, bool callSetup ) {
  _maximumFrequencyChange = change;
  if ( callSetup ) setup();
}

int PeakTracking::silentFrameCount() const {
  return _silentFrameCount;
}

void PeakTracking::setSilentFrameCount( int count, bool callSetup ) {
  _silentFrameCount = count;
  if ( callSetup ) setup();  
}