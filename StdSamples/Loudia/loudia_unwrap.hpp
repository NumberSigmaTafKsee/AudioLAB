#pragma once

namespace Loudia
{
    /**
    * @class Unwrap
    *
    * @brief Algorithm to unwrap phases vectors represented as vectors of Real values.
    *
    * This class represents an object to unwrap vectors of phases.
    * The algorithm takes as input N-point vectors of Real values 
    * and returns N-point vectors of Real values.
    *
    * Unwrapping consists in removing phase jumps larger than Pi or smaller to -Pi.
    *
    * @author Ricard Marxer
    *
    * @sa FFT
    */
    class Unwrap {
    protected:
        MatrixXR _diff;
        MatrixXR _upsteps;
        MatrixXR _downsteps;
        MatrixXR _shift;

    public:
        /**
            Constructs an unwrap object with the given @a inputSize.
        */
        Unwrap();

        /**
            Destroys the algorithm and frees its resources.
        */
        ~Unwrap();

        void setup();
        void reset();

        /**
            Performs the unwrapping on each of the rows of @a phases and
            puts the resulting unwrapped phases in the rows of @a unwrapped.
            
            @param phases matrix of Real values.
            
            @param unwrapped pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows and columns as @a phases. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */
        void process(const MatrixXR& phases, MatrixXR* unwrapped);
    };


Unwrap::Unwrap()
{
  LOUDIA_DEBUG("UNWRAP: Construction");

  setup();
}

Unwrap::~Unwrap(){}

void Unwrap::setup(){
  // Prepare the buffers
  LOUDIA_DEBUG("UNWRAP: Setting up...");
  
  reset();

  LOUDIA_DEBUG("UNWRAP: Finished setup.");
}

void Unwrap::process(const MatrixXR& input, MatrixXR* unwrapped){
  const int rows = input.rows();
  const int cols = input.cols();

  (*unwrapped).resize(rows, cols);
  
  if(input.rows() <= 1){
    (*unwrapped) = input;
  }

  _diff.resize(rows, cols);
  _upsteps.resize(rows, cols);
  _downsteps.resize(rows, cols);
  _shift.resize(rows, cols);  

  _diff << MatrixXR::Zero(1, cols), input.block(0, 0, rows-1, cols) - input.block(1, 0, rows-1, cols);
  
  _upsteps = (_diff.array() > M_PI).cast<Real>();
  _downsteps = (_diff.array() < -M_PI).cast<Real>();

  rowCumsum(&_upsteps);
  rowCumsum(&_downsteps);

  _shift =  _upsteps - _downsteps;

  (*unwrapped) = input + (2.0 * M_PI * _shift);
}

void Unwrap::reset(){
  // Initial values
}

}