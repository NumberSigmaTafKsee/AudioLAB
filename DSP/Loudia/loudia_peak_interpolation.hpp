#pragma once

namespace Loudia
{
    /**
  * @class PeakInterpolation
  *
  * @brief Algorithm to interpolate peaks in a vector of Real values.
  *
  * This class represents an object to interpolate peaks in a vector of Real values.
  * The algorithm interpolates the positions and magnitudes of a set of peaks, given
  * the original frame, peak positions and peak magnidutes.
  *
  * The interpolation consists in fitting a parabola (quadratic interpolation) on the 
  * point of the peak and the two points surrounding it. 
  *
  * Note that the interpolation is performed in the decibel domain, and in order to 
  * avoid innecessary transformations the resulting interpolated peak magnitudes
  * are returned in decibels.
  *
  * @author Ricard Marxer
  *
  * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
  */
class PeakInterpolation {
protected:
  // Internal parameters
    
  // Internal variables
  MatrixXR _magnitudes;

public:
  /**
     Constructs a peak interpolation object.
  */
  PeakInterpolation();

  /**
     Destroys the algorithm and frees its resources.
  */
  ~PeakInterpolation();

  void setup();
  void reset();

  /**
     Interpolates the peaks on each of the rows of @a frames, @a peakPositions
     and @a peakMagnitudes to put the resulting peak interpolated positions and 
     magnitudes in the rows of @a peakPositions and @a peakMagnitudes respectively.
     @param frames matrix of Real values.
     
     @param peakPositions matrix of Real values (but always Integers) for the peak indices.
     The matrix must have the same number of rows as @a frames and the same number of columns
     as @a peakMagnitudes.
     
     @param peakMagnitudes pointer to a matrix of Real values (but always Integers) for the peak indices.
     The matrix must have the same number of rows as @a frames and the same number of columns
     as @a peakPositions.
     @param peakPositionsInterpolated pointer to a matrix of Real values for the peak magnitudes.
     The matrix should have the same number of rows and columns as @a peakPositions
     and @a peakMagnitudes. 
     
     @param peakMagnitudesInterpolated pointer to a matrix of Real values for the peak magnitudes.
     The matrix should have the same number of rows and columns as @a peakPositions
     and @a peakMagnitudes.
     Note that the units of this matrix are decibels.
     Note that peaks with positions values smaller than 0 are not considered peaks and will not
     be interpolated or modified.
     
     Note that if the output matrices are not of the required size they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& frames, 
               const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
               MatrixXR* peakPositionsInterpolated, MatrixXR* peakMagnitudesInterpolated);

};

/**
    * @class PeakInterpolationComplex
    *
    * @brief Algorithm to interpolate peaks in a vector of Complex values.
    *
    * This class represents an object to interpolate peaks in a vector of Complex values.
    * The algorithm interpolates the positions and magnitudes of a set of peaks, given
    * the original frame, peak positions and peak magnidutes.
    *
    * The interpolation consists in fitting a parabola (quadratic interpolation) on the 
    * point of the peak and the two points surrounding it. 
    *
    * Note that the interpolation is performed in the decibel domain, and in order to 
    * avoid innecessary transformations the resulting interpolated peak magnitudes
    * are returned in decibels.
    *
    * @author Ricard Marxer
    *
    * @sa PeakDetection, PeakDetectionComplex, PeakInterpolation, PeakInterpolationComplex, PeakTracking, PeakTrackingComplex
    */
    class PeakInterpolationComplex {
    protected:
        // Internal parameters
            
        // Internal variables
        MatrixXR _magnitudes;
        MatrixXR _phases;

    public:
        /**
            Constructs a peak interpolation object.
        */
        PeakInterpolationComplex();

        /**
            Destroys the algorithm and frees its resources.
        */
        ~PeakInterpolationComplex();

        void setup();
        void reset();

        /**
            Interpolates the peaks on each of the rows of @a frames, @a peakPositions,
            @a peakMagnitudes, @a peakPhases to put the resulting peak interpolated positions, 
            magnitudes and phases in the rows of @a peakPositions, @a peakMagnitudes and 
            @a peakMagnitudes respectively.
            @param frames matrix of Complex values.
            
            @param peakPositions matrix of Real values (but always Integers) for the peak indices.
            The matrix must have the same number of rows as @a frames and the same number of columns
            as @a peakMagnitudes.
            
            @param peakMagnitudes pointer to a matrix of Real values for the peak magnitudes.
            The matrix must have the same number of rows as @a frames and the same number of columns
            as @a peakPositions and @a peakPhases.
            @param peakPhases pointer to a matrix of Real values for the peak phases.
            The matrix must have the same number of rows as @a frames and the same number of columns
            as @a peakMagnitudes and @a peakPositions.
            @param peakPositionsInterpolated pointer to a matrix of Real values for the peak positions.
            The matrix should have the same number of rows and columns as @a peakPositions, 
            @a peakMagnitudes and @a peakPhases. 
            
            @param peakMagnitudesInterpolated pointer to a matrix of Real values for the peak magnitudes.
            The matrix should have the same number of rows and columns as @a peakPositions, 
            @a peakMagnitudes and @a peakPhases. 
            Note that the units of this matrix are decibels.
            @param peakPhasesInterpolated pointer to a matrix of Real values for the peak phases.
            The matrix should have the same number of rows and columns as @a peakPositions, 
            @a peakMagnitudes and @a peakPhases. 
            
            Note that peaks with positions values smaller than 0.0 (usually -1.0) are not considered peaks and will not
            be interpolated or modified.
            
            Note that if the output matrices are not of the required size they will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXC& frames, 
                    const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes, const MatrixXR& peakPhases,
                    MatrixXR* peakPositionsInterpolated, MatrixXR* peakMagnitudesInterpolated, MatrixXR* peakPhasesInterpolated);


    };



PeakInterpolation::PeakInterpolation() {
  LOUDIA_DEBUG("PEAKINTERPOLATION: Constructor");
  
  setup();

  LOUDIA_DEBUG("PEAKINTERPOLATION: Constructed");
}

PeakInterpolation::~PeakInterpolation() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakInterpolation::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKINTERPOLATION: Setting up...");

  reset();

  LOUDIA_DEBUG("PEAKINTERPOLATION: Finished set up...");
}


void PeakInterpolation::process(const MatrixXR& input,
                              const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes,
                              MatrixXR* peakPositionsInterp, MatrixXR* peakMagnitudesInterp) {
  
  LOUDIA_DEBUG("PEAKINTERPOLATION: Processing");  
  Real leftMag;
  Real rightMag;
  Real mag, interpFactor;
  
  (*peakPositionsInterp).resize(input.rows(), peakPositions.cols());
  (*peakMagnitudesInterp).resize(input.rows(), peakPositions.cols());
  
  _magnitudes = input.array().abs();
  
  for ( int row = 0 ; row < _magnitudes.rows(); row++ ) {
  
    for ( int i = 0; i < peakPositions.cols(); i++ ) {
      
      // If the position is -1 do nothing since it means it is nothing
      if( peakPositions(row, i) == -1 ){

        (*peakMagnitudesInterp)(row, i) = peakMagnitudes(row, i);
        (*peakPositionsInterp)(row, i) = peakPositions(row, i);
        
      } else {
        
        // Take the center magnitude in dB
        mag = 20.0 * log10( peakMagnitudes(row, i) );

        // Take the left magnitude in dB
        if( peakPositions(row, i) <= 0 ){
          
          leftMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) + 1) );
          
        } else {
          
          leftMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) - 1) );
          
        }
        
        // Take the right magnitude in dB
        if( peakPositions(row, i) >= _magnitudes.row(row).cols() - 1 ){
          
          rightMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) - 1) );
          
        } else {
          
          rightMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) + 1) );
          
        }
                
        // Calculate the interpolated position
        (*peakPositionsInterp)(row, i) = peakPositions(row, i) + 0.5 * (leftMag - rightMag) / (leftMag - 2.0 * mag + rightMag);

        interpFactor = ((*peakPositionsInterp)(row, i) - peakPositions(row, i));

        // Calculate the interpolated magnitude in dB
        (*peakMagnitudesInterp)(row, i) = mag - 0.25 * (leftMag - rightMag) * interpFactor;
      }
    }
  }
  
  LOUDIA_DEBUG("PEAKINTERPOLATION: Finished Processing");
}

void PeakInterpolation::reset(){
  // Initial values
}


PeakInterpolationComplex::PeakInterpolationComplex() {
  LOUDIA_DEBUG("PEAKINTERPOLATIONCOMPLEX: Constructor");
  
  setup();

  LOUDIA_DEBUG("PEAKINTERPOLATIONCOMPLEX: Constructed");
}

PeakInterpolationComplex::~PeakInterpolationComplex() {
  // TODO: Here we should free the buffers
  // but I don't know how to do that with MatrixXR and MatrixXR
  // I'm sure Nico will...
}


void PeakInterpolationComplex::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("PEAKINTERPOLATIONCOMPLEX: Setting up...");

  reset();

  LOUDIA_DEBUG("PEAKINTERPOLATIONCOMPLEX: Finished set up...");
}


void PeakInterpolationComplex::process(const MatrixXC& input,
                              const MatrixXR& peakPositions, const MatrixXR& peakMagnitudes, const MatrixXR& peakPhases,
                              MatrixXR* peakPositionsInterp, MatrixXR* peakMagnitudesInterp, MatrixXR* peakPhasesInterp) {
  
  LOUDIA_DEBUG("PEAKINTERPOLATIONCOMPLEX: Processing");  
  Real leftMag, leftPhase;
  Real rightMag, rightPhase;
  Real mag, interpFactor;
  
  (*peakPositionsInterp).resize(input.rows(), peakPositions.cols());
  (*peakMagnitudesInterp).resize(input.rows(), peakPositions.cols());
  (*peakPhasesInterp).resize(input.rows(), peakPositions.cols());
  
  _magnitudes = input.array().abs();
  unwrap(input.array().angle(), &_phases);
  
  for ( int row = 0 ; row < _magnitudes.rows(); row++ ) {
  
    for ( int i = 0; i < peakPositions.cols(); i++ ) {
      
      // If the position is -1 do nothing since it means it is nothing
      if( peakPositions(row, i) == -1 ){

        (*peakMagnitudesInterp)(row, i) = peakMagnitudes(row, i);         
        (*peakPhasesInterp)(row, i) = peakPhases(row, i); 
        (*peakPositionsInterp)(row, i) = peakPositions(row, i);
        
      } else {
        
        // Take the center magnitude in dB
        mag = 20.0 * log10( peakMagnitudes(row, i) );

        // Take the left magnitude in dB
        if( peakPositions(row, i) <= 0 ){
          
          leftMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) + 1) );
          
        } else {
          
          leftMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) - 1) );
          
        }
        
        // Take the right magnitude in dB
        if( peakPositions(row, i) >= _magnitudes.row(row).cols() - 1 ){
          
          rightMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) - 1) );
          
        } else {
          
          rightMag = 20.0 * log10( _magnitudes(row, (int)peakPositions(row, i) + 1) );
          
        }
                
        // Calculate the interpolated position
        (*peakPositionsInterp)(row, i) = peakPositions(row, i) + 0.5 * (leftMag - rightMag) / (leftMag - 2.0 * mag + rightMag);

        interpFactor = ((*peakPositionsInterp)(row, i) - peakPositions(row, i));

        // Calculate the interpolated magnitude in dB
        (*peakMagnitudesInterp)(row, i) = mag - 0.25 * (leftMag - rightMag) * interpFactor;

        // Calculate the interpolated phase
        leftPhase = _phases(row, (int)floor((*peakPositionsInterp)(row, i)));
        rightPhase = _phases(row, (int)floor((*peakPositionsInterp)(row, i)) + 1);
        
        interpFactor = (interpFactor >= 0) ? interpFactor : interpFactor + 1;
        
        (*peakPhasesInterp)(row, i) = (leftPhase + interpFactor * (rightPhase - leftPhase));
      }
    }
  }

  // Calculate the princarg() of the phase: remap to (-pi pi]
  (*peakPhasesInterp) = ((*peakPhasesInterp).array() != -1).select(((*peakPhasesInterp).array() + M_PI).array().modN(-2.0 * M_PI) + M_PI, (*peakPhasesInterp));
  
  LOUDIA_DEBUG("PEAKINTERPOLATIONCOMPLEX: Finished Processing");
}

void PeakInterpolationComplex::reset(){
  // Initial values
}

}