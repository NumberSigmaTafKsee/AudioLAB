#pragma once

namespace Loudia
{
    /**
    * @class FrameCutter
    *
    * @brief Algorithm to create frames from an input stream such as the one
    * provided by the AudioLoader.
    *
    * @author Ricard Marxer
    *
    * @sa AudioLoader
    */
    class FrameCutter{
    protected:
        int _maxInputSize;
        int _frameSize;
        int _hopSize;
        Real _defaultValue;
        int _firstSamplePosition;

        int _indexWriter;
        int _availableToWrite;
        int _availableToRead;

        int _maxFrameCount;
        
        VectorXR _buffer;
        VectorXR _row;
        
        int read(VectorXR* frame, int release);
        int write(const VectorXR& stream);

    public:
        /**
            Constructs a FrameCutter object with the specified @a inputSize, 
            @a frameSize, @a hopSize settings.
            
            @param maxInputSize maximum size of the input frames,
            must be > 0.
            @param frameSize size of the output resampled frames,
            must be > 0.
            @param hopSize the hop size of the frame cutter
        */
        FrameCutter(int maxInputSize = 1024, int frameSize = 1024, int hopSize = -1, const int firstSamplePosition = 0, Real defaultValue = 0.0);

        /**
            Destroys the algorithm and frees its resources.
        */  
        ~FrameCutter();

        /**
            Performs the frame cutting.
            
            @param stream matrix of Real values.  The number of columns of @a frames 
            must be equal to the inputSize.
            
            @param resampled pointer to a matrix of Real values for the output.  The matrix should
            have the same number of rows as @a frames and outputSize columns. 
            Note that if the output matrix is not of the required size it will be resized, 
            reallocating a new memory space if necessary.
        */  
        void process(const MatrixXR& stream, MatrixXR* frames, int* produced);
        
        void setup();
        void reset();

        /**
            Returns the maximum input size of the algorithm.
            
            By default it is 1024.
        */
        int maxInputSize() const { return _maxInputSize; };
        
        /**
            Specifies the maximum input @a size of the algorithm.
        */
        void setMaxInputSize( int size, bool callSetup = true ) { _maxInputSize = size; if (callSetup) setup(); };

        /**
            Returns the frame size of the algorithm.
            
            By default it is 1024.
        */
        int frameSize() const { return _frameSize; };
        
        /**
            Specifies the frame @a size of the algorithm.
        */
        void setFrameSize( int size, bool callSetup = true ) { _frameSize = size; if (callSetup) setup(); };

        /**
            Returns the hop size of the algorithm.
            
            By default it is 1024.
        */
        int hopSize() const {
            if ( _hopSize <= 0 ) {
            return std::max(_frameSize / 2, 1);
            }
            
            return _hopSize;
        };
        
        /**
            Specifies the hop @a size of the algorithm.
        */
        void setHopSize( int size, bool callSetup = true ) { _hopSize = size; if (callSetup) setup(); };

        void setFirstSamplePosition( int position, const bool callSetup = true ) { _firstSamplePosition = position; if ( callSetup ) setup(); };
        int firstSamplePosition() const { return _firstSamplePosition; };  
        
        int maxFrameCount() const;
    };


    FrameCutter::FrameCutter( const int maxInputSize, const int frameSize, const int hopSize, const int firstSamplePosition, const Real defaultValue ) :
    _defaultValue(defaultValue)
    {
    LOUDIA_DEBUG("FRAMECUTTER: Constructing...");

    setMaxInputSize(maxInputSize, false);

    LOUDIA_DEBUG("FRAMECUTTER: Set the maximum input size...");

    setFrameSize(frameSize, false);

    LOUDIA_DEBUG("FRAMECUTTER: Set the frame size...");

    setHopSize(hopSize, false);

    LOUDIA_DEBUG("FRAMECUTTER: Set the hop size...");

    setFirstSamplePosition(firstSamplePosition, false);

    LOUDIA_DEBUG("FRAMECUTTER: Set first sample position...");

    setup();
    }

    FrameCutter::~FrameCutter(){
    
    }

    void FrameCutter::setup(){  
    const int bufferSize = 2 * _frameSize;

    // Create the stream buffer
    // must be at least twice the size than the frame size
    // in order to keep data aligned
    _buffer.resize(bufferSize);
    _buffer.setConstant(_defaultValue);
    
    if ((_firstSamplePosition < 0) || (_firstSamplePosition > _frameSize - 1)) {
        LOUDIA_ERROR("FRAMECUTTER: The first sample position must be set between 0 and frameSize-1.");
    }
    
    _indexWriter = _firstSamplePosition;
    _availableToWrite = _frameSize - _firstSamplePosition;
    _availableToRead = _firstSamplePosition;

    // TODO: check if this is the right way to know the maxFrameCount
    _maxFrameCount = maxFrameCount();

    _row = VectorXR::Zero(_frameSize);
    }

    void FrameCutter::process(const MatrixXR& stream, MatrixXR *frames, int *produced){
    if (stream.cols() != 1) {
        LOUDIA_ERROR("FRAMECUTTER: This algorithm only accepts single channel streams.");
    }

    frames->resize(_maxFrameCount, _frameSize);
    
    int currenthopSize = hopSize();
    
    int leftSize = stream.rows();
    int inputIndex = 0;
    int inputSize;
    int framesIndex = 0;
    
    while (leftSize > 0) {
        int consumed = 1;    
        while (consumed > 0) {
        inputSize = min(_frameSize, leftSize);
        if (inputSize <= 0) break;
        
        consumed = write(stream.col(0).segment(inputIndex, inputSize));
        leftSize -= consumed;
        inputIndex += consumed;
        }  
        
        while (read(&_row, currenthopSize) > 0) {
        // TODO: try to avoid the copy to frames->row(framesIndex)
        // Maybe by passing directly frames->row(framesIndex) as a const ref
        // and const casting in read()
        frames->row(framesIndex) = _row;
        framesIndex += 1;
        }
    }
    
    (*produced) = framesIndex;
    return;
    }

    int FrameCutter::read(VectorXR* frame, int release){
    if ( frame->size() > _availableToRead) return 0;
    
    const int indexReader = ((_indexWriter - _availableToRead) + _frameSize) % _frameSize;
    (*frame) = _buffer.segment(indexReader, frame->size());
    
    // Advance reader index  
    _availableToRead -= release;
    _availableToWrite += release;
    
    return frame->size();
    }

    int FrameCutter::write(const VectorXR& stream){
    int consumed = min(min(_availableToWrite, _frameSize - _indexWriter), (int)stream.size());
    if ( consumed <= 0 ) return 0;

    _buffer.segment(_indexWriter, consumed) = stream.segment(0, consumed);
    _buffer.segment(_indexWriter + _frameSize, consumed) = stream.segment(0, consumed);
    
    // Advance writer index
    _indexWriter = (_indexWriter + consumed) % _frameSize;
    _availableToWrite -= consumed;
    _availableToRead += consumed;
    
    return consumed;
    }

    void FrameCutter::reset(){
    _indexWriter = 0;
    _availableToWrite = _frameSize;
    _availableToRead = 0;
    }

    int FrameCutter::maxFrameCount() const {
    return _maxInputSize / hopSize() + 1; 
    }
}