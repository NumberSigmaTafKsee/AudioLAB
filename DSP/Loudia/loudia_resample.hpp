#pragma once

namespace Loudia
{
    /**
    * @class Resample
    *
    * @brief Algorithm to resample vectors represented as vectors of Real values.
    *
    * This class represents an object to resample vectors of Real values.
    * The algorithm takes as input N-point vectors of Real values 
    * and returns M-point vectors of Real values.
    *
    * This algorithm changes the sampling rate of input vectors, modifying their number
    * of samples.
    *
    * The number of samples in the input and output vectors can be specified using setInputSize()
    * and setOutputSize() respectively.
    *
    * The ratio between the output and input sampling rates can be specified
    * using setResamplingRatio().  Usually the sampling ratio will be chosen
    * close to outputSize / inputSize.
    *
    * @author Ricard Marxer
    *
    * @sa FFT
    */
    class Resample{
    public:
        /**
            @enum ResamplingMethod
            @brief Specifies the resampling method to be used.
            @sa resamplingMethod
        */
        enum ResamplingMethod {    
            SINC_BEST_QUALITY       = 0 /**< Best quality cardinal sine method  */,
            SINC_MEDIUM_QUALITY     = 1 /**< Medium quality cardinal sine method */,
            SINC_FASTEST            = 2 /**< Fastest cardinal sine method */,
            ZERO_ORDER_HOLD         = 3 /**< Hold the value until the next sample */,
            LINEAR                  = 4 /**< Linear interpolation between samples */
        };
        
    protected:
        int _inputSize;
        int _outputSize;
        Real _resamplingRatio;
        ResamplingMethod _resamplingMethod;

        SRC_DATA _resampleData;
        
    public:
        /**
            Constructs a Resampler object with the specified @a inputSize, 
            @a outputSize, @a resamplingRatio and @a resamplingMethod settings.
            
            @param inputSize size of the input frames to be resampled,
            must be > 0.
            @param outputSize size of the output resampled frames,
            must be > 0.
            @param resamplingRatio the ratio between the output sampling rate 
            and the input sampling rate
            
            @param resamplingMethod the resampling method to be used
        */
        Resample(int inputSize = 1024, int outputSize = 1024, Real resamplingRatio = 1.0, ResamplingMethod resamplingMethod = SINC_BEST_QUALITY);

        /**
            Destroys the algorithm and frees its resources.
        */  
        ~Resample();

        /**
            Performs the resampling of each of the rows of @a frames and
            puts the resulting resampled frames in the rows of @a resampled.
            
            @param frames matrix of Real values.  The number of columns of @a frames 
            must be equal to the inputSize.
            
            @param resampled pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a frames and outputSize columns. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */  
        void process(const MatrixXR& frames, MatrixXR* resampled);
        
        void setup();
        void reset();

        /**
            Returns the input size of the algorithm.
            
            By default it is 1024.
        */
        int inputSize() const;
        
        /**
            Specifies the input @a size of the algorithm.
        */
        void setInputSize( int size, bool callSetup = true );

        /**
            Returns the output size of the algorithm.
            
            By default it is 1024.
        */
        int outputSize() const;
        
        /**
            Specifies the output @a size of the algorithm.
        */
        void setOutputSize( int size, bool callSetup = true );

        /**
            Returns the ratio between the output and input sampling rate.
            Note that this value is normally around outputSize / inputSize.
            By default it is 1.0.
        */
        Real resamplingRatio() const;

        /**
            Specifies the @a ratio between the output and input sampling rate.
        */
        void setResamplingRatio( Real ratio, bool callSetup = true );

        /**
            Returns the resampling method to be used.
            
            By default it is SINC_BEST_QUALITY.
        */
        ResamplingMethod resamplingMethod() const;

        /**
            Specifies the resampling @a method to be used.
        */
        void setResamplingMethod( ResamplingMethod method, bool callSetup = true );

    };


Resample::Resample(int inputSize, int outputSize, Real resamplingRatio, ResamplingMethod resamplingMethod)
{
  LOUDIA_DEBUG("RESAMPLE: Constructor inputSize: " << inputSize 
        << ", outputSize: " << outputSize 
        << ", resamplingRatio: " << resamplingRatio);

  _resampleData.data_in = NULL;
  _resampleData.data_out = NULL;
  
  setInputSize( inputSize, false );
  setOutputSize( outputSize, false );
  setResamplingRatio( resamplingRatio, false );
  setResamplingMethod( resamplingMethod, false );

  setup();
  
  LOUDIA_DEBUG("RESAMPLE: Constructed");
}

Resample::~Resample(){
  LOUDIA_DEBUG("RESAMPLE: Destroying...");

  if ( _resampleData.data_in ) delete [] _resampleData.data_in;
  if ( _resampleData.data_out ) delete [] _resampleData.data_out;
  
  LOUDIA_DEBUG("RESAMPLE: Destroyed out");
}

void Resample::setup(){
  LOUDIA_DEBUG("RESAMPLE: Setting up...");

  _resampleData.input_frames = _inputSize;
  _resampleData.output_frames = _outputSize;
  _resampleData.src_ratio = _resamplingRatio;

  if ( _resampleData.data_in ) delete [] _resampleData.data_in;
  if ( _resampleData.data_out ) delete [] _resampleData.data_out;

  _resampleData.data_in = new float[_inputSize];
  _resampleData.data_out = new float[_outputSize];

  
  LOUDIA_DEBUG("RESAMPLE: Finished set up...");
}

void Resample::process(const MatrixXR& in, MatrixXR* out){
  const int rows = in.rows();
  const int cols = in.cols();

  if ( cols != _inputSize ) {
    // Throw ValueError, incorrect input size
  }

  (*out).resize(rows, _outputSize);

  for (int i = 0; i < rows; i++){    
    // Fill the buffer
    Eigen::Map<MatrixXR>(_resampleData.data_in, 1, _inputSize) = in;
    
    // Process the data
    int error = src_simple(&_resampleData, _resamplingMethod, 1);
    if ( error ) {
      // Throw ResampleError, src_strerror( error );
    }
    
    // Take the data from _out
    (*out).row( i ) = Eigen::Map<MatrixXR>(_resampleData.data_out, 1, _outputSize);
  }
}

void Resample::reset(){
}

int Resample::inputSize() const{
  return _inputSize;
}

void Resample::setInputSize( int size, bool callSetup ) {
  _inputSize = size;
  if ( callSetup ) setup();
}

int Resample::outputSize() const{
  return _outputSize;
}

void Resample::setOutputSize( int size, bool callSetup ) {
  _outputSize = size;
  if ( callSetup ) setup();
}

Real Resample::resamplingRatio() const{
  return _resamplingRatio;
}

void Resample::setResamplingRatio( Real ratio, bool callSetup ) {
  _resamplingRatio = ratio;
  if ( callSetup ) setup();
}

Resample::ResamplingMethod Resample::resamplingMethod() const{
  return _resamplingMethod;
}

void Resample::setResamplingMethod( ResamplingMethod method, bool callSetup ) {
  _resamplingMethod = method;
  if ( callSetup ) setup();
}

}